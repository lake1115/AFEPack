#ifndef PARAMETER_H
#define PARAMETER_H

#include "emdefs.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/hankel.hpp>

/*****************************************************************************
 *         
 *   Authors	: Bin Hu , Zhejiang University 
 *
 *   Project	: FEM for EM
 *
 *   $Revision: 1.0 $
 *   $Date: 2018/04/17  $
 *   
 *   @brief solve the PDE  -d/dx(c * du/dx) - d/dy(c * du/dy) + a * u = f 
 *   @input parameter c, a, f and boundary condition in parameter.h
 *
 *   Dirichlet boundary :  u = bnd;
 *   Neumann boundary :   n*c*grad(u) + q*u = g  
 *
 *****************************************************************************
 */
#define example 1
const double kappa = 2*PI;
const double lambda = 1;
const double R = 2;
const double r = 1;
const double omega = kappa*C;
const double Order = 1;

//const vec_type val_0(0);
#if example

cvaltype u_exact(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]);
  double z = kappa*r;
  cvaltype val = boost::math::cyl_hankel_1(0,z);
  return val;
}

cvec_type u_exact_prime(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]);
  double z = kappa*r;
  double j1 = boost::math::cyl_bessel_j(1,z);
  double y1 = boost::math::cyl_neumann(1,z);
  cvec_type val;
  val[0] = - (kappa*p[0]*j1)/r - (kappa*p[0]*y1*I)/r;
  val[1] = - (kappa*p[1]*j1)/r - (kappa*p[1]*y1*I)/r;
  return val;
}

cvaltype a(const double *p){
  cvaltype val = -kappa*kappa;
  return val;
}

cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}

cvaltype f(const double *p){
  cvaltype val = 0;
  double r = sqrt(p[0]*p[0]+p[1]*p[1]);
  double z = kappa*r;
  double j,y,j1,y1;
  j = boost::math::cyl_bessel_j(0,z);
  y = boost::math::cyl_neumann(0,z);
  j1 = boost::math::cyl_bessel_j(1,z);
  y1 = boost::math::cyl_neumann(1,z);

   cvaltype A = (2*kappa*j1)/r+ (kappa*y1*2*I)/r - (kappa*p[0]*p[0]*j1)/pow(r,3) - (kappa*p[1]*p[1]*j1)/pow(r,3) - (kappa*p[0]*p[0]*y1*I)/pow(r,3) - (kappa*p[1]*p[1]*y1*I)/pow(r,3) - (kappa*p[0]*((p[0]*j1)/(r*r) - (kappa*p[0]*j)/r))/r - (kappa*p[1]*((p[1]*j1)/(r*r) - (kappa*p[1]*j)/r))/r - (kappa*p[0]*((p[0]*y1)/(r*r) - (kappa*p[0]*y)/r)*I)/r - (kappa*p[1]*((p[1]*y1)/(r*r) - (kappa*p[1]*y)/r)*I)/r;

  
  val = A-kappa*kappa*u_exact(p);
  return val;

}

cvaltype q(const double *p){
  cvaltype val = 0;
}

cvaltype g(const double *p){
  double z = kappa*r;
  double j1,y1;
  j1 = boost::math::cyl_bessel_j(1,z);
  y1 = boost::math::cyl_neumann(1,z);

  cvaltype X1 = - (kappa*p[0]*j1) - (kappa*p[0]*y1*I);
  cvaltype Y1 = - (kappa*p[1]*j1) - (kappa*p[1]*y1*I);
  cvaltype val= -p[0]*X1-p[1]*Y1;
  return val;
}

cvaltype g1(const double *p){
  double z = kappa*R;
  double j1,y1;
  j1 = boost::math::cyl_bessel_j(1,z);
  y1 = boost::math::cyl_neumann(1,z);
  cvaltype X1 = - (kappa*p[0]*j1) - (kappa*p[0]*y1*I);
  cvaltype Y1 = - (kappa*p[1]*j1) - (kappa*p[1]*y1*I);
  cvaltype val= p[0]/R*X1+p[1]/R*Y1;
  return val;
}
#else 
cvaltype u_exact(const double *p){
  cvaltype val = sin(PI*p[0])*sin(PI*p[1]);
  return val;
}
cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}
cvaltype a(const double *p){
  cvaltype val = -2*PI*PI;
  return val;
}
cvaltype f(const double *p){
  cvaltype val=0;
  return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}
// for D2
cvaltype g2(const double *p){
  cvaltype val=PI/R*(p[0]*cos(PI*p[0])*sin(PI*p[1])+p[1]*sin(PI*p[0])*cos(PI*p[1]));
  return val;
}
// R = 2, bmark = 2; r = 1, bmark = 1
cvaltype g(const double *p){
  cvaltype val=-PI/r*(p[0]*cos(PI*p[0])*sin(PI*p[1])+p[1]*sin(PI*p[0])*cos(PI*p[1]));
  return val;
}

#endif

cvaltype bnd(const double *p){
  return u_exact(p);
}

#endif
