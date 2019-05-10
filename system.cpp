#include "system.h"

bool BothAreSpaces(char lhs, char rhs){ 
  return (lhs == rhs) && (lhs == ' '); 
}
void RemoveSpaces(std::string& str){
  std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
  str.erase(new_end, str.end()); 
}

///////////////////////////////////////////////////////////////////////////////////////

CSystem::CSystem(){
  this->delimiter = " ";
}
CSystem::CSystem(std::string Delimiter){
  this->delimiter = Delimiter;
}
void CSystem::setDelimiter(std::string Delimiter){
  this->delimiter = Delimiter;
}

///////////////////////////////////////////////////////////////////////////////////////

void CSystem::Mkdir(std::string folder){
  char command[150];
  sprintf(command,"if [ ! -d %s ]; then mkdir -p %s; fi",folder.c_str(),folder.c_str());
  system (command);
}
void CSystem::MkdirLoop(std::string folder,int start, int finish){
  char command[150],fold[120];
  for(int i=start;i<finish;i++){
    sprintf(fold,"%s/run%04d",folder.c_str(),i);
    sprintf(command,"if [ ! -d %s ]; then mkdir -p %s; fi",fold,fold);
    system (command);
  }
}
void CSystem::Touch(std::string file){
  char command[150];
  sprintf(command,"if [ ! -f %s ]; then touch %s; fi",file.c_str(),file.c_str());
  system (command);
}
void CSystem::PrintFile(std::string filename){
  std::string line;
  std::fstream myfile(filename,std::ios::in);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      printf("%s\n",line.c_str());
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
void CSystem::PrintFormattedFile(std::string filename){
  std::string line;
  std::fstream myfile(filename,std::ios::in);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<this->delimiter.length()){
	}else{
	  if(line.substr(line.length()-this->delimiter.length(),line.length()) != this->delimiter){
	    line.append(this->delimiter);
	  }
	}
	if(line.find(this->delimiter) == 0){
	  line.erase(0,this->delimiter.length());
	}
	printf("%s\n",line.c_str());
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
Eigen::MatrixXd CSystem::LoadFile(std::string filename){
  std::string line,
    number;
  std::size_t pos;
  std::ifstream myfile(filename);
  int rows=0,row=0,
    cols=0,col=0,
    max=0;
  Eigen::MatrixXd Matrix;

  if(myfile.is_open()){
    while(getline(myfile,line)){
      cols=0;
      if( line[0]!='#' && line!="\0" ){
	RemoveSpaces(line);
	if(line.length()<this->delimiter.length()){
	}else{
	  if(line.substr(line.length()-this->delimiter.length(),line.length()) != this->delimiter){
	    line.append(this->delimiter);
	  }
	}
	if(line.find(this->delimiter) == 0){
	  line.erase(0,this->delimiter.length());
	}
	while( (pos=line.find(this->delimiter)) != std::string::npos){
	  line.erase(0,pos+this->delimiter.length());
	  cols++;
	}
	if(cols>max){
	  max=cols;
	}
	rows++;
      }
      cols=max;
    }
    myfile.close();
  }else{
    printf("Unable to open %s.\n",filename.c_str());
    exit(1);
  }

  Matrix = Eigen::MatrixXd::Zero(rows,cols);
  myfile.open(filename);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<this->delimiter.length()){
	}else{
	  if(line.substr(line.length()-this->delimiter.length(),line.length()) != this->delimiter){
	    line.append(this->delimiter);
	  }
	}
	if(line.find(this->delimiter) == 0){
	  line.erase(0,this->delimiter.length());
	}
	col=0;
	while( (pos=line.find(this->delimiter)) != std::string::npos){
	  number=line.substr(0,pos+this->delimiter.length());
	  Matrix(row,col) = std::stof( number );
	  line.erase(0,pos+this->delimiter.length());
	  col++;
	}
	row++;
      }
    }
    myfile.close();
  }else{
    printf("Unable to open %s.\n",filename.c_str());
    exit(1);
  }
  return Matrix;
}
std::vector<Eigen::MatrixXd> CSystem::LoadFiles(std::string folder, std::string filename, int start, int finish){
  std::vector<Eigen::MatrixXd> loadMatrix = std::vector<Eigen::MatrixXd> (finish-start);
  char buffer[100];
  std::string line, file;

  for(int i=start;i<finish;i++){
    sprintf(buffer,"%s/run%04d/%s",folder.c_str(),i,filename.c_str());
    file=buffer;
    loadMatrix[i]=this->LoadFile(file);
  }

  return loadMatrix;
}

void CSystem::WriteFile(std::string filename, Eigen::MatrixXd &Matrix){
  std::string line;
  int rows = Matrix.rows(),
    cols = Matrix.cols();
  std::fstream ofile(filename,std::ios::out);
  if(ofile.is_open()){
    for(int row=0;row<rows;row++){
      for(int col=0;col<cols;col++){
	ofile << Matrix(row,col) << this->delimiter;
      }
      ofile << "\n";
    }
    ofile.close();
  }else{
    printf("Unable to open %s.\n",filename.c_str());
  }  
}
void CSystem::WriteFile(std::string filename, Eigen::VectorXd &Vector){
  std::string line;
  int rows = Vector.size();
  std::fstream ofile(filename,std::ios::out);
  if(ofile.is_open()){
    for(int row=0;row<rows;row++){
      ofile << Vector(row) << this->delimiter;
    }
    ofile << "\n";
    ofile.close();
  }else{
    printf("Unable to open %s.\n",filename.c_str());
  }  
}
void CSystem::WriteCSVFile(std::string filename, std::vector<std::string> header, Eigen::MatrixXd &Matrix){
  std::string line;
  int rows = Matrix.rows(),
    cols = Matrix.cols();
  std::fstream ofile(filename,std::ios::out);
  if(ofile.is_open()){
    for(int col=0;col<cols-1;col++){
      ofile << "'" << header[col] << "',";
    }
    ofile << "'" << header[cols-1] << "'\n";
    for(int row=0;row<rows;row++){
      for(int col=0;col<cols-1;col++){
	ofile << Matrix(row,col) << ",";
      }
      ofile << Matrix(row,cols-1) << "\n";
    }
    ofile.close();
  }else{
    printf("Unable to open %s.\n",filename.c_str());
  }  
}


/*
void CSystem::LoadParamFile(std::string filename, std::vector<std::string> &Names, Eigen::MatrixXd &Matrix,std::string delimiter){
  std::string line,
    number;
  std::size_t pos,length;
  std::ifstream myfile(filename);
  int rows=0,row=0,
    cols=0,col=0,
    max=0;

  if(myfile.is_open()){
    while(getline(myfile,line)){
      cols=0;
      if( line[0]!='#' && line!="\0" ){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	while( (pos=line.find(delimiter)) != std::string::npos){
	  line.erase(0,pos+delimiter.length());
	  cols++;
	}
	if(cols>max){
	  max=cols;
	}
	rows++;
      }
      cols=max;
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }

  Names = std::vector<std::string>(rows);
  Matrix = Eigen::MatrixXd::Zero(rows,cols-2);
  myfile.open(filename);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	col=0;
	while( (pos=line.find(delimiter)) != std::string::npos ){
	  number=line.substr(0,pos);
	  if(col==0){
	  } else if(col==1){
	    Names[row] = number;
	  } else{
	    Matrix(row,col-2) = std::stof( number );
	  }
	  line.erase(0,pos+delimiter.length());
	  col++;
	}
	row++;
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}

void CSystem::WriteParamFile(std::string fileName, std::vector<std::string> &header, std::string delimiter, Eigen::MatrixXd &file){
  int rows = file.rows();
  int cols = file.cols();
  std::fstream ofile(fileName,std::ios::out);
  if(ofile.is_open()){
    for(int i=0;i<cols;i++){
      ofile << header[i] << delimiter;
      for(int j=0;j<rows;j++){
	ofile << file(j,i) << delimiter;
      }
      ofile << '\n';
    }
    ofile.close();
  }else{
    printf("Unable to open %s.\n",fileName.c_str());
  }
}
void CSystem::WriteParamFileLoop(std::string filename, std::string folder, int start,std::vector<std::string> &header, std::string delimiter, Eigen::MatrixXd &matrix){
  int rows = matrix.rows();
  int cols = matrix.cols();
  std::fstream ofile;
  char file[100];
  for(int row=0;row<rows;row++){
    sprintf(file,"%s/run%04d/%s",folder.c_str(),row+start,filename.c_str());
    ofile.open(file,std::ios::out);
    if(ofile.is_open()){
      for(int col=0;col<cols;col++){
	ofile << header[col] << delimiter << matrix(row,col) << delimiter << '\n';
      }
      ofile.close();
    }else{
      printf("Unable to open %s\n",file);
    }
  }
}
void CSystem::WriteParameterFiles(std::string rangename, std::string foldername, std::string filename, std::string delimiter, int start, int finish, int ab){
  std::vector<std::string> Name;
  Eigen::MatrixXd range,matrix;
  LoadParamFile(rangename,Name,range,delimiter);
  LHCSampling(matrix,finish-start,ab,range);
  MkdirLoop(foldername,start,finish);
  WriteParamFileLoop(filename,foldername,start,Name,delimiter,matrix);
}
void CSystem::WritePosteriorParameterFiles(std::string foldername, std::string filename, std::vector<std::string> Name, std::string paramname, std::string delimiter, int start, int finish){
  Eigen::MatrixXd matrix;
  
  LoadFile(filename,matrix,delimiter);
  MkdirLoop(foldername,start,finish);
  WriteParamFileLoop(paramname,foldername,start,Name,delimiter,matrix);
}
void CSystem::WriteParameterFiles(std::string rangename, std::string foldername, std::string filename, std::string delimiter, int start, int finish, int ab, Eigen::MatrixXd &Parameters){
  std::vector<std::string> Name;
  Eigen::MatrixXd range;
  LoadParamFile(rangename,Name,range,delimiter);
  LHCSampling(Parameters,finish-start,ab,range);
  MkdirLoop(foldername,start,finish);
  WriteParamFileLoop(filename,foldername,start,Name,delimiter,Parameters);
}
void CSystem::LHCSampling(Eigen::MatrixXd &hypercube, int samples, int ab, Eigen::MatrixXd range){
  int parameters = range.rows();
  int nmax = parameters/ab - 2;

  Eigen::VectorXd G = Eigen::VectorXd::Zero(ab);
  std::vector<int> hyperlist;
  hypercube = Eigen::MatrixXd::Zero(samples,parameters);

  for(int i=0;i<samples;i++){
    hyperlist.push_back(i);
  }

  for(int i=0;i<parameters;i++){
    std::random_shuffle(hyperlist.begin(),hyperlist.end());
    for(int j=0;j<samples;j++){
      hypercube(j,i) = hyperlist[j];
    }
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    for(int j=0;j<samples;j++){
      hypercube(j,i) = init + dx*hypercube(j,i);
    }
  }

  if(nmax>-1){
    CCosh dist(nmax);
    for(int i=0;i<samples;i++){
      for(int j=0;j<ab;j++){
	hypercube(i,j*(nmax+2)+1) = 1.0/hypercube(i,j*(nmax+2));
	for(int k=1;k<nmax+1;k++){
	  hypercube(i,j*(nmax+2)+1) -= hypercube(i,j*(nmax+2)+k+1)*dist.Z(k);
	}
	hypercube(i,j*(nmax+2)+1) /= (dist.Z(0));
      }
    }
  }
}
*/
