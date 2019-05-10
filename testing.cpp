#include<iostream>
#include<vector>
#include "neuralnetwork.h"
#include "system.h"
#include "processing.h"
#include "parametermap.h"

int main(int argc, char* argv[]){
  CSigmoidNet nn;
  CSystem system;

  //Loading csv data files
  system.setDelimiter(",");
  Eigen::MatrixXd X_train = system.LoadFile("MNIST_X_train.csv");
  Eigen::MatrixXd X_test = system.LoadFile("MNIST_X_test.csv");
  Eigen::MatrixXd y_train = system.LoadFile("MNIST_y_train.csv");
  Eigen::MatrixXd y_test = system.LoadFile("MNIST_y_test.csv");
  //Processing Data
  Eigen::MatrixXd X_train_scaled = Scale(X_train.transpose());
  Eigen::MatrixXd X_test_scaled = Scale(X_test.transpose());
  Eigen::MatrixXd y_test_vec = Vectorize(y_test.transpose(),10);
  Eigen::MatrixXd y_train_vec = Vectorize(y_train.transpose(),10);

  CParameterMap MAP;
  MAP.ReadParsFromFile("emulator.info");
  int Epochs=MAP.getI("EPOCHS",10);

  nn.Construct(X_train_scaled, y_train_vec, MAP);
  nn.Train(Epochs);

  Eigen::MatrixXd test = nn.Max(nn.Emulate(X_test_scaled));
  Eigen::MatrixXd ytest = Eigen::MatrixXd::Zero(test.cols(),1);
  for(int i=0;i<test.cols();i++){
    for(int j=0;j<10;j++){
      if(test(j,i)==1.0) ytest(i,0)=j;
    }
  }

  printf("Accuracy of test set: %d out of %d\n",(int)nn.Accuracy(X_test_scaled,y_test_vec),(int)X_test_scaled.cols());
  system.setDelimiter(" ");
  system.WriteFile("NN_output.dat",ytest);
}
