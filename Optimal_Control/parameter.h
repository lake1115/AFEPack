#ifndef PARAMETER_H
#define PARAMETER_H

#include "typedefs.h"
#include <AFEPack/FEMSpace.h>
/*****************************************************************************
 *         
 *   @file       parameter.h      
 *   @author	 Hu Bin, Zhejiang University <binh@zju.edu.cn> 
 *   @version    1.0 
 *   @date       2018/10/25  
 *   
 *   @brief solve the optimal control PDE:  
 *        -div(A nabla y) + phi(y) = f + Bu
 *        min(g(y)+j(y))
 *
 *   Vartiational Formulation:
 *   (A nabla y, nabla w) + (phi(y), w) = (f + Bu, w)
 *   (A nabla p, nabla q) + (grad(phi(y)), q) = (grad(g(y)), q)
 *   (grad(j(y)) + Bp, v-u) >= 0
 *
 *   @input parameter matrix A, B and function phi(y), g(y), j(y)
 *
 *   Dirichlet boundary :  u = bnd;
 *   Neumann boundary :   n*A*grad(y) + py = g  ?  
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

#if example


cvaltype u_exact(const double *p){
  cvaltype val = sin(PI*p[0])*sin(PI*p[1]);
  return val;
}


void setA(Eigen::Matrix2d& A){
  A(0,0) = 1;
  A(0,1) = 0;
  A(1,0) = 0;
  A(1,1) = 1;
}

void setB(Eigen::Matrix2d& B){
  B(0,0) = 0;
  B(0,1) = 0;
  B(1,0) = 0;
  B(1,1) = 0;
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
  cvaltype val=1;
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
