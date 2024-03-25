/**
 * @file   twin_triangle.RT.1.bas_fun.cpp
 * @author Guanghui Hu
 * @date   Fri Oct 15 10:35:10 2004
 * 
 * @brief  basis functions on twin-triangle for R-T 0 element
 * 
 * 
 */

#include <cmath>
#include <Miscellaneous.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * 三角形上 1 阶的棱单元的基函数及其导数的定义。其中插值点在第 \f$ i
   * \f$ 条边上的基函数的表达式为
   *
   * \f[
   *    \vec{\psi}_i = |p_j - p_k| \left(\lambda_j (x) \nabla
   *    \lambda_k (x) - \lambda_k (x) \nabla \lambda_j (x))
   * \f]
   *
   * 其中 \f$ \lambda_i \f$ 是第 \f$ i \f$ 个面积坐标函数，\f$ j = (i
   * + 1)%3 \f$, \f$ k = (i + 2)%3 \f$。
   * 
   * 因为这个文件中用到了类 nVector<2,double> ，必须使用 C++编译器
   * 进行编译。
   *
   */

#define vector_length 2
#define vt nVector<vector_length,double>

#define ZERO (1.0e-6)
#define DISTANCE(p0, p1) sqrt(((p0)[0] - (p1)[0])*((p0)[0] - (p1)[0]) + \
                              ((p0)[1] - (p1)[1])*((p0)[1] - (p1)[1]))
#define AREA(p0, p1, p2) (((p1)[0] - (p0)[0])*((p2)[1] - (p0)[1]) - \
                          ((p1)[1] - (p0)[1])*((p2)[0] - (p0)[0]))
#define GET_L(i, j)                                                            \
  double l0 = (v[i][0] - v[j][0]);                                             \
  double l1 = (v[i][1] - v[j][1]);                                             \
  double l = sqrt(l0*l0 + l1*l1);                                              \
  if (fabs(l0) > fabs(l1)) {                                                   \
    if (l0 < 0) l = -l;                                                        \
  } else {                                                                     \
    if (l1 < 0) l = -l;                                                        \
  }
#define dlambda0_dx ((v[0][1] - v[1][1])/area)
#define dlambda0_dy ((v[1][0] - v[0][0])/area)
#define dlambda1_dx ((v[1][1] - v[2][1])/area)
#define dlambda1_dy ((v[2][0] - v[1][0])/area)
#define dlambda2_dx ((v[2][1] - v[3][1])/area)
#define dlambda2_dy ((v[3][0] - v[2][0])/area)
#define dlambda3_dx ((v[3][1] - v[0][1])/area)
#define dlambda3_dy ((v[0][0] - v[3][0])/area)
#define dlambda4_dx ((v[2][1] - v[0][1])/area) // maybe a mistake here, 
#define dlambda4_dy ((v[0][0] - v[2][0])/area) // suppose p is in triangle 0, 1, 2


void get_lambda(const double * p, const double ** v, double * lambda, double * area)
{
  area[0] = AREA(v[0], v[1], v[2]);
  lambda[0] = AREA(p, v[0], v[1])/area[0];
  lambda[1] = AREA(p, v[1], v[2])/area[0];
  lambda[2] = AREA(p, v[2], v[3])/area[0];
  lambda[3] = AREA(p, v[3], v[0])/area[0];
  lambda[4] = AREA(p, v[2], v[0])/area[0]; /// suppose p is in triangle 0, 1, 2
}

void psi_1(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  double lambda[5];
  if (p_area < -area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    get_lambda(p, v, lambda, &area);
    GET_L(0, 1);
    val[0] = l*(lambda[1]*dlambda4_dx - lambda[4]*dlambda1_dx);
    val[1] = l*(lambda[1]*dlambda4_dy - lambda[4]*dlambda1_dy);
  }
  if (p_area < area) {
    val[0] /= 2;
    val[1] /= 2;
  }
}

void psi_2(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  double lambda[5];
  if (p_area < -area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    get_lambda(p, v, lambda, &area);
    GET_L(1, 2);
    val[0] = l*(lambda[4]*dlambda0_dx - lambda[0]*dlambda4_dx);
    val[1] = l*(lambda[4]*dlambda0_dy - lambda[0]*dlambda4_dy);
  }
  if (p_area < area) {
    val[0] /= 2;
    val[1] /= 2;
  }
}

void psi_3(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  double lambda[5];
  if (p_area >  area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    get_lambda(p, v, lambda, &area);
    GET_L(2, 3);
    val[0] = l*(- lambda[3]*dlambda4_dx + lambda[4]*dlambda3_dx);
    val[1] = l*(- lambda[3]*dlambda4_dy + lambda[4]*dlambda3_dy);
  }
  if (p_area > -area) {
    val[0] /= 2;
    val[1] /= 2;
  }
}

void psi_4(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  double lambda[5];
  if (p_area >  area) {
    val[0] = 0.0;
    val[1] = 0.0;
  }
  else {
    get_lambda(p, v, lambda, &area);
    GET_L(3, 0);
    val[0] = l*(- lambda[4]*dlambda2_dx + lambda[2]*dlambda4_dx);
    val[1] = l*(- lambda[4]*dlambda2_dy + lambda[2]*dlambda4_dy);
  }
  if (p_area > -area) {
    val[0] /= 2;
    val[1] /= 2;
  }
}

void psi_5(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  double lambda[5];
  get_lambda(p, v, lambda, &area);
  if (p_area >  area) {
    GET_L(2, 0);
    val[0] = l*(lambda[0]*dlambda1_dx - lambda[1]*dlambda0_dx);
    val[1] = l*(lambda[0]*dlambda1_dy - lambda[1]*dlambda0_dy);

  }
  else if (p_area < -area){
    GET_L(0, 2);
    val[0] = l*(lambda[2]*dlambda3_dx - lambda[3]*dlambda2_dx);
    val[1] = l*(lambda[2]*dlambda3_dy - lambda[3]*dlambda2_dy);
  }
  else{
    {
      GET_L(2, 0);
      val[0] = l*(lambda[0]*dlambda1_dx - lambda[1]*dlambda0_dx);
      val[1] = l*(lambda[0]*dlambda1_dy - lambda[1]*dlambda0_dy);
    }
    {
      GET_L(0, 2);
      val[0] += l*(lambda[2]*dlambda3_dx - lambda[3]*dlambda2_dx);
      val[1] += l*(lambda[2]*dlambda3_dy - lambda[3]*dlambda2_dy);
    }
    
    val[0] /= 2.;    
    val[1] /= 2.;
  }
}


void gradient_psi_1(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  if (p_area < -area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(0, 1);
    val[0][0] = 0.0; val[0][1] = l*(dlambda1_dy*dlambda4_dx - dlambda4_dy*dlambda1_dx);
    val[1][0] = -val[0][1]; val[1][1] = 0.0;
  }
  if (p_area < area) {
    val[0][1] /= 2;
    val[1][0] /= 2;
  }
}

void gradient_psi_2(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  if (p_area < -area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(1, 2);
    val[0][0] = 0.0; val[0][1] = l*(dlambda4_dy*dlambda0_dx - dlambda0_dy*dlambda4_dx);
    val[1][0] = -val[0][1]; val[1][1] = 0.0;
  }
  if (p_area < area) {
    val[0][1] /= 2;
    val[1][0] /= 2;
  }
}

void gradient_psi_3(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  if (p_area >  area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(2, 3);
    val[0][0] = 0.0; val[0][1] = l*(-dlambda3_dy*dlambda4_dx + dlambda4_dy*dlambda3_dx);
    val[1][0] = -val[0][1]; val[1][1] = 0.0;
  }
  if (p_area > -area) {
    val[0][1] /= 2;
    val[1][0] /= 2;
  }
}

void gradient_psi_4(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  if (p_area >  area) {
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;
  }
  else {
    GET_L(3, 0);
    val[0][0] = 0.0; val[0][1] = l*(- dlambda4_dy*dlambda2_dx + dlambda2_dy*dlambda4_dx);
    val[1][0] = -val[0][1]; val[1][1] = 0.0;
  }
  if (p_area > -area) {
    val[0][1] /= 2;
    val[1][0] /= 2;
  }

}

void gradient_psi_5(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area = AREA(v[0], v[1], v[2]);
  double p_area = AREA(v[0], p, v[2])/ZERO;
  if (p_area >  area) {
    GET_L(2, 0);
    val[0][0] = 0.0; val[0][1] = l*(dlambda0_dy*dlambda1_dx - dlambda1_dy*dlambda0_dx);
    val[1][0] = -val[0][1]; val[1][1] = 0.0;

  }
  else if (p_area < -area){
    GET_L(0, 2);
    val[0][0] = 0.0; val[0][1] = l*(dlambda2_dy*dlambda3_dx - dlambda3_dy*dlambda2_dx);
    val[1][0] = -val[0][1]; val[1][1] = 0.0;
  }
  else{
    val[0][0] = 0.0; val[0][1] = 0.0;
    val[1][0] = 0.0; val[1][1] = 0.0;

  }

}
#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
