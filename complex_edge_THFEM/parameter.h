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
//#if 1
cvec_type u_exact(const double *p){
  cvec_type val;
  val[0] = cos(PI*p[0])*sin(PI*p[1]);
  val[1] = -sin(PI*p[0])*cos(PI*p[1]);
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
 return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

cvec_type g1(const double *p){
  cvec_type val;
  val[0] = -2.*PI*cos(PI*p[0]);
  val[1] = 0;
  return val;
}

cvec_type g2(const double *p){
  cvec_type val;
  val[0] = 0;
  val[1] = -2.*PI*cos(PI*p[1]);
  return val;
}

cvec_type g3(const double *p){
  cvec_type val;
  val[0] = 2.*PI*cos(PI*p[0]);
  val[1] = 0;
  return val;
}

cvec_type g4(const double *p){
  cvec_type val;
  val[0] = 0;
  val[1] = 2.*PI*cos(PI*p[1]);
  return val;
}

cvaltype bnd(const double *p){
    return 0;
}

/*
#else 0
cvaltype u_exact(const double *p){
  cvaltype val = p[0]*p[0]+p[1]*p[1]-p[0]*p[1]*(p[0]+p[1]);
  return val;
}

cvaltype a(const double *p){
  cvaltype val(0,0);
  return val;
}

double c(const double *p){
  return 1;
}

cvaltype f(const double *p){
  cvaltype val=4.0-2.0*(p[0]+p[1]);
 return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype g(const double *p){
  cvaltype val=PI*sin(PI*p[1]);
  return val;
}

cvaltype bnd1(const double *p){
  return 2.0*p[0]*p[0]-p[0]+1;
}
cvaltype bnd2(const double *p){
  return 1.0-p[1];
}
cvaltype bnd3(const double *p){
  return 1.0-p[0];
}
cvaltype bnd4(const double *p){
  return 2.0*p[1]*p[1]-p[1]+1;
}

#endif
*/
/*
#else
cvaltype u_exact(const double *p){
  cvaltype val = p[0]*p[0]+p[1]*p[1]-p[0]*p[1]*(p[0]+p[1]);
  return val;
}

cvaltype a(const double *p){
  cvaltype val(0,0);
  return val;
}

double c(const double *p){
  return 1;
}

cvaltype f(const double *p){
  cvaltype val=4.0-2.0*(p[0]+p[1]);
 return val;
}

cvaltype q(const double *p){
  cvaltype val(0,0);
  return val;
}

cvaltype g1(const double *p){
  cvaltype val=-(2.0-2*p[1]-p[1]*p[1]);
  return val;
}

cvaltype g2(const double *p){
  cvaltype val=-(2.0-2*p[0]-p[0]*p[0]);
  return val;
}

cvaltype bnd1(const double *p){
  return p[0]*p[0];
}

cvaltype bnd2(const double *p){
  return p[1]*p[1];
}
#endif
*/
#endif
