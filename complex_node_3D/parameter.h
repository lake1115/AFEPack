#ifndef PARAMETER_H
#define PARAMETER_H

#include "emdefs.h"
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
#define example 0
const double kappa = 1;
const double lambda = 1;
const double R = 2;
const double r = 1;
const double omega = kappa*C;
const double Order = 5;

//const vec_type val_0(0);
#if example

cvaltype u_exact(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype val = exp(kappa*I*r)/r;
  return val;
}

cvec_type u_exact_prime(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvec_type val;
  val[0] = - (p[0]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[0]*exp(kappa*I*r))/(r*r);
  val[1] = - (p[1]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[1]*exp(kappa*I*r))/(r*r);
  val[2] = - (p[2]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[2]*exp(kappa*I*r))/(r*r);
 
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
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype ddp0 = - exp(kappa*I*r)/pow(r,3) + (kappa*I*exp(kappa*I*r))/(r*r) + (3.*p[0]*p[0]*exp(kappa*I*r))/pow(r,5) - (kappa*p[0]*p[0]*exp(kappa*I*r)*3.*I)/pow(r,4) - (kappa*kappa*p[0]*p[0]*exp(kappa*I*r))/pow(r,3);
  cvaltype ddp1 = - exp(kappa*I*r)/pow(r,3) + (kappa*I*exp(kappa*I*r))/(r*r) + (3.*p[1]*p[1]*exp(kappa*I*r))/pow(r,5) - (kappa*p[1]*p[1]*exp(kappa*I*r)*3.*I)/pow(r,4) - (kappa*kappa*p[1]*p[1]*exp(kappa*I*r))/pow(r,3);
  cvaltype ddp2 = - exp(kappa*I*r)/pow(r,3) + (kappa*I*exp(kappa*I*r))/(r*r) + (3.*p[2]*p[2]*exp(kappa*I*r))/pow(r,5) - (kappa*p[2]*p[2]*exp(kappa*I*r)*3.*I)/pow(r,4) - (kappa*kappa*p[2]*p[2]*exp(kappa*I*r))/pow(r,3);
  
  val = -(ddp0+ddp1+ddp2)-kappa*kappa*u_exact(p);
  return val;

}

cvaltype q(const double *p){
  cvaltype val = 0;
}

cvaltype g(const double *p){
  cvaltype dx = - (p[0]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[0]*exp(kappa*I*r))/(r*r);
  cvaltype dy = - (p[1]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[1]*exp(kappa*I*r))/(r*r);
  cvaltype dz = - (p[2]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[2]*exp(kappa*I*r))/(r*r);

  cvaltype val= -p[0]*dx-p[1]*dy-p[2]*dz;
  return val;
}

cvaltype g1(const double *p){
  double r = 2.0;
  cvaltype d0 = - (p[0]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[0]*exp(kappa*I*r))/(r*r);
  cvaltype d1 = - (p[1]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[1]*exp(kappa*I*r))/(r*r);
  cvaltype d2 = - (p[2]*exp(kappa*I*r))/pow(r,3) + (kappa*I*p[2]*exp(kappa*I*r))/(r*r);
  cvaltype val= p[0]/r*d0+p[1]/r*d1+p[2]/r*d2;
  return val;
}
#else
// cube case for test

cvaltype u_exact(const double *p){
  cvaltype val = sin(PI*p[0])*sin(PI*p[1])*sin(PI*p[2]);
  return val;
}
/*
cvec_type u_exact_prime(const double *p){
  cvec_type val;
  val[0] = PI*cos(PI*p[0])*sin(PI*p[1])*sin(PI*p[2]);
  val[1] = PI*cos(PI*p[1])*sin(PI*p[2])*sin(PI*p[0]);
  val[2] = PI*cos(PI*p[2])*sin(PI*p[0])*sin(PI*p[1]);
  return val;
}
*/
cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}
cvaltype a(const double *p){
  cvaltype val = 0;
  return val;
}
cvaltype f(const double *p){
  cvaltype val=3*PI*PI*u_exact(p);
  return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}
// R = 2, bmark = 2; r = 1, bmark = 1


 //for cube n=(0,1,0)
cvaltype g2(const double *p){
  cvaltype val= PI*cos(PI*p[1])*sin(PI*p[0])*sin(PI*p[2]);
  return val;
}
/*
cvaltype g(const double *p){
  // double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  
  cvaltype dx = PI*cos(PI*p[0])*sin(PI*p[1])*sin(PI*p[2]);
  cvaltype dy = PI*cos(PI*p[1])*sin(PI*p[2])*sin(PI*p[0]);
  cvaltype dz = PI*cos(PI*p[2])*sin(PI*p[0])*sin(PI*p[1]);
  cvaltype val= -p[0]/r*dx - p[1]/r*dy - p[2]/r*dz;
  return val;
}

cvaltype g2(const double *p){
  //double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype dx = PI*cos(PI*p[0])*sin(PI*p[1])*sin(PI*p[2]);
  cvaltype dy = PI*cos(PI*p[1])*sin(PI*p[2])*sin(PI*p[0]);
  cvaltype dz = PI*cos(PI*p[2])*sin(PI*p[0])*sin(PI*p[1]);
  cvaltype val= p[0]/r*dx +p[1]/r*dy +p[2]/r*dz;
  return val;
}
*/
/*
cvaltype u_exact(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype val = 1/r;
  return val;
}

cvec_type u_exact_prime(const double *p){
  cvec_type val;
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  val[0] = -p[0]/pow(r,3);
  val[1] = -p[1]/pow(r,3);
  val[2] = -p[2]/pow(r,3);
  return val;
}

cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}
cvaltype a(const double *p){
  cvaltype val = 0;
  return val;
}
cvaltype f(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype ddx = -1/pow(r,3)+3*p[0]*p[0]/pow(r,5);
  cvaltype ddy = -1/pow(r,3)+3*p[1]*p[1]/pow(r,5);
  cvaltype ddz = -1/pow(r,3)+3*p[2]*p[2]/pow(r,5);
  
  cvaltype val= -(ddx+ddy+ddz);
  return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype g(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  
  cvaltype dx = -p[0]/pow(r,3);
  cvaltype dy = -p[1]/pow(r,3);
  cvaltype dz = -p[2]/pow(r,3);
  cvaltype val= -p[0]/r*dx - p[1]/r*dy - p[2]/r*dz;
  return val;
}

cvaltype g2(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype dx = -p[0]/pow(r,3);
  cvaltype dy = -p[1]/pow(r,3);
  cvaltype dz = -p[2]/pow(r,3);
  cvaltype val= p[0]/r*dx + p[1]/r*dy + p[2]/r*dz;
  return val;
}

/*
cvaltype u_exact(const double *p){
  cvaltype val =exp(p[0]);
  return val;
}

cvec_type u_exact_prime(const double *p){
  cvec_type val;
  val[0] = exp(p[0]);
  val[1] = 0;
  val[2] = 0;
  return val;
}

cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}
cvaltype a(const double *p){
  cvaltype val = exp(p[0]);
  return val;
}
cvaltype f(const double *p){
  cvaltype val= 0;
  return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype g(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  
  cvaltype dx = exp(p[0]);
  cvaltype dy = 0;
  cvaltype dz = 0;
  cvaltype val= -p[0]/r*dx - p[1]/r*dy - p[2]/r*dz;
  return val;
}

cvaltype g2(const double *p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  cvaltype dx = -p[0]/pow(r,3);
  cvaltype dy = -p[1]/pow(r,3);
  cvaltype dz = -p[2]/pow(r,3);
  cvaltype val= p[0]/r*dx + p[1]/r*dy + p[2]/r*dz;
  return val;
}
*/
#endif

cvaltype bnd(const double *p){
  return u_exact(p);
}

#endif
