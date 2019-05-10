#include "neuralnetwork.h"

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

void CNeuralNetwork::Construct(Eigen::MatrixXd X_train, Eigen::MatrixXd Y_train, CParameterMap Map){
  this->X = X_train;
  this->Y = Y_train;
  this->MAP = Map;

  int Params = this->X.rows();
  int Obs = this->Y.rows();
  int Samples = this->X.cols();

  this->learning_rate = this->MAP.getD("LEARNING_RATE",0.1);
  this->mini_batch = this->MAP.getD("MINI_BATCH_SIZE",10);
  this->layers = this->MAP.getD("LAYERS",1);
  this->Layers.push_back(this->MAP.getI("LAYER0_SIZE",100));
  for(int i=1;i<layers;i++){
    std::string layer_name = "LAYER"+std::to_string(i)+"_SIZE";
    this->Layers.push_back(this->MAP.getI(layer_name,0));
  }
  this->Layers.insert(this->Layers.begin(),Params);
  this->Layers.push_back(Obs);
  this->layers = this->Layers.size();

  Eigen::MatrixXd weight;
  Eigen::MatrixXd bias;
  Eigen::MatrixXd del;
  int layer, next_layer;

  this->activations.push_back(Eigen::MatrixXd::Zero(Params,mini_batch));
  this->zs.push_back(Eigen::MatrixXd::Zero(Params,mini_batch));

  for(int i=0;i<this->layers-1;i++){
    layer=this->Layers[i];
    next_layer=this->Layers[i+1];

    weight = Eigen::MatrixXd::Zero(layer,next_layer);
    this->Weight.push_back(weight);
    this->delta_w.push_back(weight);
    this->Weight[i].setRandom();

    bias = Eigen::MatrixXd::Zero(next_layer,1);
    this->Bias.push_back(bias);
    this->delta_b.push_back(bias);
    this->Bias[i].setRandom();

    del = Eigen::MatrixXd::Zero(next_layer,mini_batch);
    this->delta.push_back(del);

    this->activations.push_back(del);
    this->zs.push_back(del);
  }
}
Eigen::MatrixXd CNeuralNetwork::Emulate(Eigen::MatrixXd X){
  Eigen::MatrixXd 
    output=X,
    temp;
  int samples=X.cols();
  for(int ilayer=0;ilayer<this->layers-1;ilayer++){
    temp = output;
    output = this->Weight[ilayer].transpose()*temp;
    temp = output;
    for(int i=0;i<samples;i++){
      temp.col(i) += this->Bias[ilayer];
    }
    output = this->Activation(temp);
  }
  //output = this->SoftMax(temp); //Softmax as final function
  return output;
}
void CNeuralNetwork::FeedForward(Eigen::MatrixXd x){
  int samples=x.cols();
  this->zs[0] = x;
  this->activations[0] = x;
  for(int ilayer=0;ilayer<this->layers-1;ilayer++){
    this->zs[ilayer+1] = this->Weight[ilayer].transpose()*this->activations[ilayer];
    for(int i=0;i<samples;i++){
      this->zs[ilayer+1].col(i) += this->Bias[ilayer];
    }
    this->activations[ilayer+1] = this->Activation(this->zs[ilayer+1]);
  }
}
void CNeuralNetwork::BackPropagation(Eigen::MatrixXd x, Eigen::MatrixXd y){
  int layer, n_layer,final_layer;
  int samples = x.cols();
  final_layer = this->layers-1;
  layer = this->Layers[final_layer-1];
  n_layer = this->Layers[final_layer];

  Eigen::MatrixXd yrun = this->Emulate(x);
  Eigen::MatrixXd cost_der = this->LossDerivative(yrun,y);
  Eigen::MatrixXd act_der = this->ActivationDerivative(zs[final_layer]);

  this->delta[final_layer-1] = Eigen::MatrixXd::Zero(n_layer,samples);

  for(int i=0;i<n_layer;i++){
    this->delta_b[final_layer-1](i,0) = 0.0; // rezero d_b
    for(int samp=0;samp<samples;samp++){
      this->delta[final_layer-1](i,samp) = cost_der(i,samp)*act_der(i,samp); //set delta
      this->delta_b[final_layer-1](i,0) += this->delta[final_layer-1](i,samp); //sum d_b across all samples
    }
    this->delta_b[final_layer-1](i,0) /= (double) samples;
    for(int j=0;j<layer;j++){
      this->delta_w[final_layer-1](j,i) = 0.0; // rezero d_w
      for(int samp=0;samp<samples;samp++){
	this->delta_w[final_layer-1](j,i) += this->delta[final_layer-1](i,samp)*activations[final_layer-1](j,samp); // sum d_w across all samples
      }
      this->delta_w[final_layer-1](j,i) /= (double) samples;
    }
  }
  
  double del=0.0;
  for(int ilay=this->layers-2;ilay>0;ilay--){
    layer=this->Layers[ilay-1];
    n_layer=this->Layers[ilay];
    this->delta[ilay-1] = Eigen::MatrixXd::Zero(n_layer,samples);
    act_der = this->ActivationDerivative(zs[ilay]);
    //for(int samp=0;samp<samples;samp++){ // d(l+1) = d(l)*W.transpose
    this->delta[ilay-1] = this->Weight[ilay]*this->delta[ilay];
      //}
    for(int i=0;i<n_layer;i++){ //calculate delta & d_b
      this->delta_b[ilay-1](i,0) = 0.0; //rezero d_b
      for(int samp=0;samp<samples;samp++){
	del=this->delta[ilay-1](i,samp);
	this->delta[ilay-1](i,samp) = del*act_der(i,samp);//coeff wise product d(l+1) = d(l+1)*f'(a(l))
	this->delta_b[ilay-1](i,0) += this->delta[ilay-1](i,samp); //sum d_b over all samples
      }
      this->delta_b[ilay-1](i,0) /= (double) samples;

      for(int j=0;j<layer;j++){
	this->delta_w[ilay-1](j,i) = 0.0; //rezero d_w
	for(int samp=0;samp<samples;samp++){ //d_w_ji = a(j)*d_l(i)
	  this->delta_w[ilay-1](j,i) += this->delta[ilay-1](i,samp)*activations[ilay-1](j,samp);
	}
	this->delta_w[ilay-1](j,i) /= (double) samples;
      }
    }
  }
}
void CNeuralNetwork::Train(int Epochs){
  printf("Loss: %f\n",this->Loss(this->Emulate(X),Y));
  int samples=this->X.cols();
  int runs=ceil((double)samples/(double)this->mini_batch);
  int params = this->X.rows();
  int obs = this->Y.rows();
  int start, finish, batch = this->mini_batch;
  Eigen::MatrixXd x;
  Eigen::MatrixXd y;
  for(int i=0;i<Epochs;i++){
    printf("Starting Epoch%d...\n",i);
    for(int j=0;j<runs;j++){
      start = j*mini_batch;
      finish = start + batch;
      if(start > samples-1) break;
      if(finish > samples) batch = samples - start;
      x = X.block(0,start,params,batch);
      y = Y.block(0,start,obs,batch);
      this->FeedForward(x);
      this->BackPropagation(x,y);
      for(int ilay=0;ilay<this->layers-1;ilay++){ //Update coeffs
	this->Weight[ilay] -= this->learning_rate*this->delta_w[ilay];
	this->Bias[ilay] -= this->learning_rate*this->delta_b[ilay];
      }
    }
    printf("Loss: %f\n",this->Loss(this->Emulate(X),Y));
  }
}
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

double CSigmoidNet::Loss(Eigen::MatrixXd AL,Eigen::MatrixXd Y){
  int rows=AL.rows(),
    cols=AL.cols();
  double loss=0.0;
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      loss += (AL(i,j) - Y(i,j))*(AL(i,j) - Y(i,j));
    }
  }
  return loss;
}
Eigen::MatrixXd CSigmoidNet::LossDerivative(Eigen::MatrixXd AL, Eigen::MatrixXd Y){
  int rows=AL.rows(),
    cols=AL.cols();
  Eigen::MatrixXd Loss = Eigen::MatrixXd::Zero(rows,cols);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      Loss(i,j) = AL(i,j) - Y(i,j);
    }
  }
  return Loss;
}
Eigen::MatrixXd CSigmoidNet::Activation(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      temp(irow,icol) = 1.0/(1.0+exp(-Z(irow,icol)));
    }
  }
  return temp;
}
Eigen::MatrixXd CSigmoidNet::ActivationDerivative(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = this->Activation(Z);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      temp(irow,icol) = temp(irow,icol)*(1.0 - temp(irow,icol));
    }
  }
  return temp;
}
Eigen::MatrixXd CSigmoidNet::Max(Eigen::MatrixXd Z){
  int rows=Z.rows(),
    cols=Z.cols(),
    index=0;
  double max=0.0;
  Eigen::MatrixXd Max = Eigen::MatrixXd::Zero(rows,cols);
  for(int icol=0;icol<cols;icol++){
    max=Z(0,icol);
    index=0;
    for(int irow=0;irow<rows;irow++){
      if(Z(irow,icol) > max){
	max = Z(irow,icol);
	index = irow;
      }
    }
    Max(index,icol) = 1.0;
  }
  return Max;
}
Eigen::MatrixXd CSigmoidNet::SoftMax(Eigen::MatrixXd Z){
  int rows=Z.rows(),
    cols=Z.cols();
  double sum=0.0;
  Eigen::MatrixXd softMax = Eigen::MatrixXd::Zero(rows,cols);
  for(int icol=0;icol<cols;icol++){
    sum=0.0;
    for(int irow=0;irow<rows;irow++){
      softMax(irow,icol) = exp(Z(irow,icol));
      sum += softMax(irow,icol);
    }
    softMax.col(icol) /= sum;
  }
  return softMax;
}
Eigen::MatrixXd CSigmoidNet::SoftMaxDerivative(Eigen::MatrixXd Z){ //NEED TO FINISH
  return Z;
}
double CSigmoidNet::Accuracy(){
  Eigen::MatrixXd ytest = this->Max(this->Emulate(this->X));
  double acc=this->X.cols();
  for(int j=0;j<this->Y.cols();j++){
    for(int i=0;i<this->Y.rows();i++){
      acc -= abs(this->Y(i,j) - ytest(i,j))/2.0;
    }
  }
  return acc;
}
double CSigmoidNet::Accuracy(Eigen::MatrixXd x, Eigen::MatrixXd y){
  Eigen::MatrixXd ytest = this->Max(this->Emulate(x));
  double acc=x.cols();
  for(int j=0;j<y.cols();j++){
    for(int i=0;i<y.rows();i++){
      acc -= abs(y(i,j) - ytest(i,j))/2.0;
    }
  }
  return acc;
}


///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

/*
Eigen::MatrixXd CNeuralNetwork::Activation(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      if(temp(irow,icol) < 0.0) temp(irow,icol) = 0.0;
      else temp(irow,icol) = Z(irow,icol);
    }
  }
  return temp;
}
Eigen::MatrixXd CNeuralNetwork::ActivationDerivative(Eigen::MatrixXd Z){
  int rows=Z.rows(), cols=Z.cols();
  Eigen::MatrixXd temp = Eigen::MatrixXd::Ones(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      if(Z(irow,icol) < 0.0) temp(irow,icol) = 0.0;
    }
  }
  return temp;
}
*/
