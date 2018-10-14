#ifndef PARAMETER_H
#define PARAMETER_H

#include "emdefs.h"
/*****************************************************************************
 *         
 *   @file       parameter.h      
 *   @author	 Bin Hu , Zhejiang University <binh@zju.edu.cn> 
 *   @version    1.0 
 *   @date       2018/04/17  
 *   
 *   @brief solve the PDE  -d/dx(c * du/dx) - d/dy(c * du/dy) + a * u = f
 * 
 *   input parameter c, a, f and boundary condition in parameter.h
 *
 *   Dirichlet boundary :  u = bnd;
 *   Neumann boundary :   n*c*grad(u) + q*u = g  
 *        3
 *     -------
 *    |      |
 *  4 |      | 2      Dirichlet BC mark from 1
 *    |      |        Neumann BC mark from 101
 *    -------
 *       1
 *****************************************************************************
 */
#define example 1
double lambda = 10.;
double freq = C/lambda;
double k0 = 2*PI*lambda;

//const vec_type val_0(0);
#if example
cvaltype u_exact(const double *p){
  cvaltype val = sin(PI*p[0])*sin(PI*p[1]);
  return val;
}

cvaltype a(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}

cvaltype f(const double *p){
  cvaltype val=2.0*PI*PI*u_exact(p);
  return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype g1(const double *p){
  cvaltype val=PI*sin(PI*p[0]);
  return val;
}
cvaltype g2(const double *p){
  cvaltype val=-PI*sin(PI*p[1]);
  return val;
}
cvaltype g3(const double *p){
  cvaltype val=-PI*sin(PI*p[0]);
  return val;
}

cvaltype g4(const double *p){
  cvaltype val=PI*sin(PI*p[1]);
  return val;
}

cvaltype bnd(const double *p){
  return u_exact(p);
}

#else 
cvaltype u_exact(const double *p){
  cvaltype val = p[0]*p[0]+p[1]*p[1]-p[0]*p[1]*(p[0]+p[1]);
  return val;
}

cvaltype a(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype c(const double *p){
  cvaltype val(-1,0);
  return val;
}

cvaltype f(const double *p){
  cvaltype val=4.0-2.0*(p[0]+p[1]);
  return val;
}

cvaltype q(const double *p){
  cvaltype val(3,0);
  return val;
}

cvaltype g1(const double *p){
  cvaltype val=2*p[1]-2*p[0]*p[1]-p[0]*p[0]+q(p)*u_exact(p);
  return val;
}
cvaltype g2(const double *p){
  cvaltype val=-2*p[0]+2*p[0]*p[1]+p[1]*p[1]+q(p)*u_exact(p);
  return val;
}

cvaltype bnd(const double *p){
  return u_exact(p);
}

#endif

#endif
