#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <stdlib.h>
#include <cmath>
using namespace std;

int main(int argc,char *argv[])
{
  // define refine number 
  double H = atof(argv[1]);
  // radius of circle boundary
  double R = atof(argv[2]);
	
  ofstream meshfile("D.d");
  // total number of points

  int n_circle = int(2*M_PI*R/H);
  double theta = 2*M_PI/n_circle;
  int n_point = n_circle;
  //std::cout<<"H = "<<H<<std::endl;
  
  meshfile << n_point <<'\t' << "# number of points #" << endl;

  std::cout<<" n_point = "<< n_point<<" n_circle = "<<n_circle<<std::endl;
  meshfile <<"# Nodes which define the boundary # " <<endl;
  double x,y;
  int N = 0;
  for(int i=0;i<n_circle;i++){

    x = R*cos(theta*i);
    y = R*sin(theta*i);
    meshfile << N << ":" << '\t' << x << '\t' << y <<'\t' << H <<'\t'<< 2<<endl;
    N++;
  }
  
  
  if(N==n_point)
    cout<< "it's good!"<<endl;
 
  
  int segments = n_point; 
  meshfile << segments <<'\t' << "# number of segments #" << endl;
  int i;
  meshfile <<"# Boundary segments #"<<endl;
  for(i = 0;i < n_circle-1;i++){
    meshfile << i << ":" << '\t'<< i <<'\t'<<i+1<<'\t'<<2<<endl;
  }
 
  meshfile << i << ":" << '\t'<< i <<'\t'<<0<<'\t'<<2<<endl;
  

  return 0;
}

