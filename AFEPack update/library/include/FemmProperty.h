
// File    :   FemmProperty.h
// Time    :   2024/03/12 14:54:19
// Author  :   Hu Bin
// Version :   1.0
#ifndef FEMMPROPERTY_H
#define FEMMPROPERTY_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <cstring>
#include <complex>
#include <Eigen/Dense>

typedef std::complex<double> Complex;

#define PI (4.0*atan(1.0))
#define mu0 (4e-7 * PI)  // PI * 4e-7  Vacuum permeability
#define I Complex(0,1)
#define DEG PI/180.

class LinearSystem{
 public:
  Eigen::MatrixXcd A;
  Eigen::VectorXcd b;
    
  LinearSystem(int d);
  Eigen::VectorXcd solve();
};


// Boundary Property -- structure that holds information about boundary.
class BoundaryProp {
 public:
  std::string name;
  int type = 0;  /// 0: Dirichlet, 1: Neumann, 2: Periodic, 3: Antiperiodic
  double value = 0;   /// A_0 
};

// Material Property -- structure that holds information about material.
class MaterialProp {
 public:
  std::string name;
  double mu = 1.;    /// permeabilities, relative
  double H_c = 0;    /// magnetization, A/m
  double J = 0;      /// applied current density, MA/m^2
  double sigma = 0;  /// conductivity of the material, MS/m
  double phi_h = 0;  /// hysteresis angle, degrees
  double lamfill = 1.;  /// lamination fill factor
  std::vector<double> Bdata;  /// for B-H curve
  std::vector<Complex> Hdata;
  std::vector<Complex> slope;
  
  void GetSlopes(double omega);
  void GetBHProps(double B, Complex &v, Complex &dv);
  
};

// Circuit Property -- structure that holds information about circuit.
class CircuitProp {
 public:
  std::string name;
  int type;  /// type = 1 means for block
  double total_amps = 0;
};

// Block -- structure that holds information about block.
class BlockLabel {
 public:
  double x = 0.;       /// block center x
  double y = 0.;       /// block center y
  int m_idx = 0;      /// material idx
  int c_idx = -1;     /// circuit idx, -1 means without circuit
  int turns = 0;      /// circuit turns
  double volume = 0.;  /// block volume, need get by femspace
  double magdir = 0.;  /// Magnetic direction [-180, 180]
};

struct PBCPair
{
  int p1;   /// point 1 idx
  int p2;   /// point 2 idx
  int type; /// 0: Periodic, 1: Antiperiodic
};


class FemmProp {
 public:
  double length_unit = 1.;  // centimeters as deault
  std::vector<BoundaryProp> bprop;
  std::vector<MaterialProp> mprop;
  std::vector<CircuitProp> cprop;
  std::vector<BlockLabel> block;
  double precision = 1e-5;
  bool nonlinear = false;  // for nonlinear problem
  std::vector<PBCPair> pbclist; // for periodic and antiperiodic boundary in .pdc file 

 public:
  FemmProp(){};
  FemmProp(const std::string&);
  virtual ~FemmProp();
  void reinit(const std::string&); /**< Reinitialization. */

};


#endif
