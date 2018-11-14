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
  // the number of the points in the 1/4 circle  
  int n_quarter_circle = atoi(argv[1]);
  // here is rectangle obstacle
  double X = atof(argv[2]);
  double Y = atof(argv[3]);
  // the distance between fictitious surface and obstacle 
  double width = atof(argv[4]);
  // define refine number, usually is dependent by n_quarter_circle 
  double h;
  //double h = atof(argv[5]);
	
  ofstream meshfile("D.d");
  // total number of points
  int n_point = n_quarter_circle*4+4;
  double theta = 0.5*M_PI/(n_quarter_circle-1);

  h = width*sin(theta);
  std::cout<<"h = "<<h<<std::endl;
  meshfile << n_point <<'\t' << "# number of points #" << endl;

  meshfile <<"# Nodes which define the boundary # " <<endl;
  double dx,dy,x,y;
  int N = 0;
  for(int i=0;i<n_quarter_circle;i++){

    dx = -X/2;
    dy = -Y/2;
    x = -width*cos(theta*i);
    y = -width*sin(theta*i);
    meshfile << N << ":" << '\t' << x+dx << '\t' << y+dy <<'\t' << h <<'\t'<< 2<<endl;
    N++;
  }
  for(int i =0;i<n_quarter_circle;i++){
    dx = X/2;
    dy = -Y/2;
    x = width*sin(theta*i);
    y = -width*cos(theta*i); 
    meshfile << N << ":" << '\t' << x+dx << '\t' << y+dy <<'\t' << h <<'\t'<< 2<<endl;
    N++;
      }
  for(int i =0;i<n_quarter_circle;i++){
    dx = X/2;
    dy = Y/2;
    x = width*cos(theta*i);
    y = width*sin(theta*i);
    meshfile << N << ":" << '\t' << x+dx << '\t' << y+dy <<'\t' << h <<'\t'<< 2<<endl;
    N++;
  }
  for(int i =0;i<n_quarter_circle;i++){
    dx = -X/2;
    dy = Y/2;
    x = -width*sin(theta*i);
    y = width*cos(theta*i);
    meshfile << N << ":" << '\t' << x+dx << '\t' << y+dy <<'\t' << h <<'\t'<< 2<<endl;
    N++;
  }

  meshfile <<"# Nodes which define the hole #" <<endl;
  
  meshfile << N << ":" << '\t' << -X/2 << '\t' << -Y/2 <<'\t' << h <<'\t'<< 1<<endl;
  N++;
  meshfile << N << ":" << '\t' << -X/2 << '\t' << Y/2 <<'\t' << h <<'\t'<< 1<<endl;
  N++;
  meshfile << N << ":" << '\t' << X/2 << '\t' << Y/2 <<'\t' << h <<'\t'<< 1<<endl;
  N++;
  meshfile << N << ":" << '\t' << X/2 << '\t' << -Y/2 <<'\t' << h <<'\t'<< 1<<endl;
  N++;
  if(N==n_point)
    cout<< "it's good!"<<endl;
 
  
  int segments = (n_quarter_circle-1)*4+8; 
  meshfile << segments <<'\t' << "# number of segments #" << endl;
  int i;
  meshfile <<"# Boundary segments #"<<endl;
  for(i = 0;i< n_quarter_circle*4-1;i++){
    meshfile << i << ":" << '\t'<< i <<'\t'<<i+1<<'\t'<<2<<endl;
  }
 
  meshfile << i << ":" << '\t'<< i <<'\t'<<0<<'\t'<<2<<endl;
  i++;
  meshfile <<"# Hole segments #"<<endl;
 
  meshfile << i << ":" << '\t'<< i <<'\t'<<i+1<<'\t'<<1<<endl;
  meshfile << i+1 << ":" << '\t'<< i+1 <<'\t'<<i+2<<'\t'<<1<<endl;
  meshfile << i+2 << ":" << '\t'<< i+2 <<'\t'<<i+3<<'\t'<<1<<endl;
  meshfile << i+3 << ":" << '\t'<< i+3 <<'\t'<<i<<'\t'<<1<<endl;
  return 0;
}

