#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<Eigen/Dense>
#include<string>
#include<vector>


bool BothAreSpaces(char lhs, char rhs);
void RemoveSpaces(std::string& str);

class CSystem{
 public:
  std::string delimiter;

  CSystem();
  CSystem(std::string delimiter);
  void setDelimiter(std::string Delimiter);

  void Mkdir(std::string folder);
  void MkdirLoop(std::string folder, int start, int finish);
  void Touch(std::string file);

  void PrintFile(std::string filename);
  void PrintFormattedFile(std::string filename);

  Eigen::MatrixXd LoadFile(std::string filename);
  std::vector<Eigen::MatrixXd> LoadFiles(std::string folder, std::string filename, int start, int finish);

  void WriteFile(std::string filename, Eigen::MatrixXd &Matrix);
  void WriteFile(std::string filename, Eigen::VectorXd &Vector);
  void WriteCSVFile(std::string filename, std::vector<std::string> header, Eigen::MatrixXd &Matrix);

  /*
  void LoadParamFile(std::string filename, std::vector<std::string> &Names, Eigen::MatrixXd &Matrix,std::string delimiter);
  void LoadDataFile(std::string folder, std::string filename, std::string delimiter, int start, int finish, int column, Eigen::MatrixXd &Matrix, Eigen::MatrixXd &MatrixError,Eigen::MatrixXd &MatrixVar);
  void LoadDataFiles(std::string folder, std::vector<std::string> filenames, std::string delimiter, int start, int finish, int column, std::vector<Eigen::MatrixXd> &Matrix, std::vector<Eigen::MatrixXd> &MatrixError,std::vector<Eigen::MatrixXd> &MatrixVar);

  void LoadMEDataFiles(std::vector<std::string> modelfilenames, std::vector<std::string> expfilenames,
		     std::vector<Eigen::MatrixXd> &ModelMatrix, std::vector<Eigen::MatrixXd> &ExpMatrix,
		     std::vector<Eigen::MatrixXd> &ModelError, std::vector<Eigen::MatrixXd> &ExpError,
		     std::vector<Eigen::MatrixXd> &ModelVar, std::vector<Eigen::MatrixXd> &ExpVar,
		     Eigen::VectorXd &modeldy, Eigen::VectorXd &expdy,
		     std::string foldername, std::string delimiter,
		     int start, int finish);

  void WriteParamFile(std::string fileName, std::vector<std::string> &header, std::string delimiter, Eigen::MatrixXd &file);
  void WriteParamFileLoop(std::string filename, std::string folder, int start, std::vector<std::string> &header, std::string delimiter, Eigen::MatrixXd &matrix);
  void WritePosteriorParameterFiles(std::string foldername, std::string filename, std::vector<std::string> Name, std::string paramname, std::string delimiter, int start, int finish);
  void WriteParameterFiles(std::string rangename, std::string foldername, std::string filename, std::string delimiter, int start, int finish, int ab);
  */
/*************************************************************
   loads RangeName file to create parameter files 
      from Start(0) to Finish(1000) in Foldername/run%04d/Filename
   This function also hands a copy of the parameters
      in the Parameters matrix
   AB denotes how many different quark pairs, uu, ud, us, ss = 4
*************************************************************/
  /*
  void WriteParameterFiles(std::string rangename, std::string foldername, std::string filename, std::string delimiter, int start, int finish, int ab, Eigen::MatrixXd &Parameters);
  void WriteGABFunctions(std::string outfilename, Eigen::MatrixXd Parameters, std::string delimiter, int ab);
  void WriteGABFunctions(std::string infilename, std::string delimiter, int ab);

  void LHCSampling(Eigen::MatrixXd &hypercube, int samples, int ab, Eigen::MatrixXd range);
*/
};



#endif
