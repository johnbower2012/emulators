#include "processing.h"

Eigen::MatrixXd Scale(Eigen::MatrixXd matrix){
  int rows=matrix.rows(),
    cols=matrix.cols();
  double max=matrix(0,0),
    min = matrix(0,0),
    range;
  Eigen::MatrixXd Scaled = Eigen::MatrixXd::Zero(rows,cols);
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      if(max < matrix(irow,icol)) max = matrix(irow,icol);
      if(min > matrix(irow,icol)) min = matrix(irow,icol);
    }
  }
  range=max-min;
  for(int irow=0;irow<rows;irow++){
    for(int icol=0;icol<cols;icol++){
      Scaled(irow,icol) = (matrix(irow,icol) - min)/range;
    }
  }
  return Scaled;
}
Eigen::MatrixXd Vectorize(int number, int size){
  Eigen::MatrixXd vector = Eigen::MatrixXd::Zero(size,1);
  vector(number,0) = 1;
  return vector;
}
Eigen::MatrixXd Vectorize(Eigen::MatrixXd number, int size){
  int cols = number.cols();
  Eigen::MatrixXd vectors = Eigen::MatrixXd::Zero(size,cols);
  for(int icol=0;icol<cols;icol++){
    vectors(number(0,icol),icol) = 1;
  }

  return vectors;
}
