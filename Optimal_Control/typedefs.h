/*****************************************************************************
 *
 *   @file       parameter.h      
 *   @author	 Hu Bin, Zhejiang University <binh@zju.edu.cn> 
 *   @version    1.0 
 *   @date       2018/10/25        
 *  
 *   @file typedefs.h
 *   @brief Constants, definitions 
 *
 *****************************************************************************
 */
#ifndef EMDEFS_H
#define EMDEFS_H

#include <cmath>
#include <complex>
#include <Eigen/Dense>
//#define vec_type nVector<vector_length,double>	/// vector_value type

typedef double valtype;
typedef std::complex<valtype> cvaltype;

const int DIM = 2;
const int DOW = DIM;
const int TDIM = DIM;

const cvaltype I(0,1);

//typedef vec_type (*vFunc)(const double*);
typedef double (*Func)(const double*);
typedef cvaltype (*CFunc)(const double*);

 

#define PI   (4.0*atan(1.0))

 

#endif
