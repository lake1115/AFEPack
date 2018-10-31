#ifndef PARAMETER_H
#define PARAMETER_H

#include "typedefs.h"
/*****************************************************************************
 *         
 *   @file       parameter.h      
 *   @author	 Hu Bin, Zhejiang University <binh@zju.edu.cn> 
 *   @version    1.0 
 *   @date       2018/10/25  
 *   
 *   @brief solve the optimal control PDE:  
 *        -div(A nabla y) + phi(y) = f + B(u)
 *        min(g(y)+j(y))
 *
 *   Vartiational Formulation:
 *   (A nabla y, nabla w) + (phi(y), w) = (f + B(u), w)
 *   (A nabla p, nabla q) + (grad(phi(y)), q) = (grad(g(y)), q)
 *   (grad(j(y)) + B*(p), v-u) >= 0
 *
 *   @input parameter matrix A and function phi(y), g(y), j(y) ,B(u),
 *  B* (the adjoint operator of B) 
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


double U_exact(const double *p){
  double Z = sin(PI*p[0])*sin(PI*p[1]);
  double val = std::max(-Z,0.0);
  return val;
}

double u_0(const double *p){
  return 0;
}

void setA(Eigen::Matrix2d& A,const double *p){
  A(0,0) = 1;
  A(0,1) = 0;
  A(1,0) = 0;
  A(1,1) = 1;
}


cvaltype f(const double *p){
  cvaltype val=1;
  return val;
}

double y_0(const double *p){
  return 0.0;
}

cvaltype bnd(const double *p){
  return 0.0;
}
#else
void setA(Eigen::Matrix2d& A,const double *p){
  A(0,0) = 1;
  A(0,1) = 0;
  A(1,0) = 0;
  A(1,1) = 1;
}


double z(const double *p){
  if(p[0]+p[1] >1.0)
    return 0.5;
  else
    return 0.0;
}
double u_0(const double *p){
  return 1.0 - sin(PI/2*p[0])-sin(PI/2*p[1])+z(p);
}

cvaltype u_exact(const double *p){
  double Z = sin(PI*p[0])*sin(PI*p[1]);
  cvaltype val = std::max(u_0(p)-Z,0.0);
  return val;
}

double y_0(const double *p){
  return 100*sqrt((p[0]-1)*(p[0]-1)+(p[1]-1)*(p[1]-1));
}

cvaltype f(const double *p){
  double Z = sin(PI*p[0])*sin(PI*p[1]);
  return 4*pow(PI,4)*Z-u_exact(p)-10000/y_0(p);
}

double Y_exact(const double *p){
  double Z = sin(PI*p[0])*sin(PI*p[1]);
  return 2*PI*PI*Z + y_0(p);
}

cvaltype bnd(const double *p){
  return 0;
}
#endif





#endif
