/*****************************************************************************         
 *  
 *   @file emdefs.h
 *   @author Hu Bin <binh@zju.edu.cn>
 *   @data 2018/10/23
 *
 *   @brief Constants, definitions 
 *
 *****************************************************************************
 */
#ifndef EMDEFS_H
#define EMDEFS_H

#include <cmath>
#include <complex>


typedef double valtype;
typedef std::complex<valtype> cvaltype;

//typedef nVector<vector_length,cvaltype> cvec_type;

#ifndef COMPLEX
typedef valtype emtype;
#else
typedef cvaltype emtype;
#endif

const cvaltype I(0,1);


//typedef vec_type (*vFunc)(const double*);
typedef double (*Func)(const double*);
typedef cvaltype (*CFunc)(const double*);
#ifdef EMUNITS
// universal units :um
const double C=1.;		/// speed of light
#else
// time units are femtoseconds, space units are nanometers
const double C=300;		/// speed of light ,nm/fs
#endif

const double eps0=1./C;	/// vacuum electrical permittivity
const double mue0=1./C;	/// vacuum magnetic permeability 

#define PI   (4.0*atan(1.0))

 

#endif