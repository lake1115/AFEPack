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
  // radius of circle obstacle
  double n_o_circle = atof(argv[2]);
  // radius of circle boundary
  double n_circle = atof(argv[3]);
	
  ofstream meshfile("D.d");
  // total number of points

  double R = 2;
  double r = 1;
  //int n_circle = int(2*M_PI*R/H);
  //int n_o_circle = int(2*M_PI*r/H);
  double theta = 2*M_PI/n_circle;
  double theta_o = 2*M_PI/n_o_circle;
  int n_point = n_circle+n_o_circle;
  //std::cout<<"H = "<<H<<std::endl;
  
  meshfile << n_point <<'\t' << "# number of points #" << endl;

  std::cout<<" n_point = "<< n_point<<" n_circle = "<<n_circle<<" n_o_circle = "<<n_o_circle<<std::endl;
  meshfile <<"# Nodes which define the boundary # " <<endl;
  double x,y;
  int N = 0;
  for(int i=0;i<n_circle;i++){

    x = R*cos(theta*i);
    y = R*sin(theta*i);
    meshfile << N << ":" << '\t' << x << '\t' << y <<'\t' << H <<'\t'<< 2<<endl;
    N++;
  }
  
  meshfile <<"# Nodes which define the hole #" <<endl;
  for(int i=0;i<n_o_circle;i++){

    x = r*cos(theta_o*i);
    y = r*sin(theta_o*i);
    meshfile << N << ":" << '\t' << x << '\t' << y <<'\t' << H <<'\t'<< 1<<endl;
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
  i++;
  meshfile <<"# Hole segments #"<<endl;
  int k = i;
  for(int j = n_point -1; j> k;j--,i++){
    meshfile << i << ":" << '\t'<< j <<'\t'<<j-1<<'\t'<<1<<endl;
  }
 
  meshfile << n_point-1 << ":" << '\t'<< k <<'\t'<<n_point-1<<'\t'<<1<<endl;

  return 0;
}

