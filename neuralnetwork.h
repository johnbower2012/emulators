#ifndef __NEURAL_NETWORK_H__
#define __NEURAL_NETWORK_H__

#include "emulator.h"
#include "parametermap.h"
#include <math.h>

class CNeuralNetwork : public CEmulator{
public:
  CParameterMap MAP;
  int layers;
  double learning_rate;
  int mini_batch;

  Eigen::MatrixXd X;
  Eigen::MatrixXd Y;

  std::vector<int> Layers;
  
  std::vector<Eigen::MatrixXd> Weight;
  std::vector<Eigen::MatrixXd> Bias;

  std::vector<Eigen::MatrixXd> activations;
  std::vector<Eigen::MatrixXd> zs;

  std::vector<Eigen::MatrixXd> delta;
  std::vector<Eigen::MatrixXd> delta_w;
  std::vector<Eigen::MatrixXd> delta_b;

  CNeuralNetwork(){};
  virtual ~CNeuralNetwork(){};
  
  void Construct(Eigen::MatrixXd X_train, Eigen::MatrixXd Y_train, CParameterMap Map); //virtual from Emulator
  Eigen::MatrixXd Emulate(Eigen::MatrixXd X); //virtual from Emulator

  void FeedForward(Eigen::MatrixXd x);
  void BackPropagation(Eigen::MatrixXd x, Eigen::MatrixXd y);
  void Train(int Epochs);

  virtual double Loss(Eigen::MatrixXd AL, Eigen::MatrixXd Y) =0;
  virtual Eigen::MatrixXd LossDerivative(Eigen::MatrixXd AL, Eigen::MatrixXd Y) =0;
  virtual Eigen::MatrixXd Activation(Eigen::MatrixXd Z) =0;
  virtual Eigen::MatrixXd ActivationDerivative(Eigen::MatrixXd Z) =0;

  virtual double Accuracy() =0;
  virtual double Accuracy(Eigen::MatrixXd x, Eigen::MatrixXd y) =0;
};

///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////

class CSigmoidNet : public CNeuralNetwork{
public:
  CSigmoidNet(){};
  ~CSigmoidNet(){};
  
  double Loss(Eigen::MatrixXd AL, Eigen::MatrixXd Y);
  Eigen::MatrixXd LossDerivative(Eigen::MatrixXd AL, Eigen::MatrixXd Y);
  Eigen::MatrixXd Activation(Eigen::MatrixXd Z);
  Eigen::MatrixXd ActivationDerivative(Eigen::MatrixXd Z);

  Eigen::MatrixXd Max(Eigen::MatrixXd Z);
  Eigen::MatrixXd SoftMax(Eigen::MatrixXd Z);
  Eigen::MatrixXd SoftMaxDerivative(Eigen::MatrixXd Z);

  double Accuracy();
  double Accuracy(Eigen::MatrixXd x, Eigen::MatrixXd y);
};


#endif
