#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include<Eigen/Dense>
#include<chrono>
#include<random>
#include<cmath>
#include<iostream>
#include "parametermap.h"

class CEmulator{
 public:
  virtual ~CEmulator(){}
  virtual void Construct(Eigen::MatrixXd X_train, Eigen::MatrixXd Y_train, CParameterMap MAP) =0;
  virtual Eigen::MatrixXd Emulate(Eigen::MatrixXd testX) =0;
};


#endif
