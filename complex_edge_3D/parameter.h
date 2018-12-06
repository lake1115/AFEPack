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
 *   @brief solve the PDE  curl (c * curl u) + a * u = f 
 *
 *   input parameter c, a, f and boundary condition in parameter.h
 *
 *   Dirichlet boundary : n times u = bnd;
 *   Neumann boundary :   n times (c * curl u) + q*n times (n times u) = g  
 *        3
 *     -------
 *    |      |
 *  4 |      | 2
 *    |      |
 *    -------
 *       1
 *****************************************************************************
 */
double lambda = 10.;
double freq = C/lambda;
double k0 = 2*PI*lambda;

const vec_type val_0(0);
#if 1
cvec_type u_exact(const double *p){
  cvec_type val;
  val[0] = sin(PI*p[1])*sin(PI*p[2]);
  val[1] = sin(PI*p[0])*sin(PI*p[2]);
  val[2] = sin(PI*p[0])*sin(PI*p[1]);
  return val;
}

cvaltype a(const double *p){
  cvaltype val(1,0);
  return val;
}

cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}

cvec_type f(const double *p){
  cvec_type val;
  val[0] = (2.*PI*PI+1.)*u_exact(p)[0];
  val[1] = (2.*PI*PI+1.)*u_exact(p)[1];
  val[2] = (2.*PI*PI+1.)*u_exact(p)[2];
 return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

// for cube n = (0,1,0) 
cvec_type g(const double *p){
  cvec_type val;
  val[0] = PI*sin(PI*p[2])*(cos(PI*p[0])-cos(PI*p[1]));
  val[1] = 0;
  val[2] = -PI*sin(PI*p[0])*(cos(PI*p[1])-cos(PI*p[2]));
  return val;
}

cvaltype bnd(const double *p){
  return 0;
}

#else 
cvec_type u_exact(const double *p){
  cvec_type val;
  val[0] = cos(PI*p[0])*sin(PI*p[1]);
  val[1] = -sin(PI*p[0])*cos(PI*p[1]);
  val[2] = 0;
  return val;
}

cvaltype a(const double *p){
  cvaltype val(1,0);
  return val;
}

cvaltype c(const double *p){
  cvaltype val(1,0);
  return val;
}

cvec_type f(const double *p){
  cvec_type val;
  val[0] = (2.*PI*PI+1.)*u_exact(p)[0];
  val[1] = (2.*PI*PI+1.)*u_exact(p)[1];
  val[2] = 0;
 return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}


cvaltype bnd(const double *p){
    return 0;
}


#endif
#endif
