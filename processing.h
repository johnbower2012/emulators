#ifndef __PROCESSING_H__
#define __PROCESSING_H__

#include<Eigen/Dense>

Eigen::MatrixXd Scale(Eigen::MatrixXd matrix);
Eigen::MatrixXd Vectorize(int number, int size);
Eigen::MatrixXd Vectorize(Eigen::MatrixXd number, int size); //should be samples by 1 in size

#endif
