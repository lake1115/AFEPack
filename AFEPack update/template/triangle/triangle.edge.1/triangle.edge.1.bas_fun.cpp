/**
 * @file   triangle.edge.1.bas_fun.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Dec 19 13:39:28 2008
 * 
 * @brief  �������ϵ������㵥Ԫ
 * 
 * 
 */

#include <cmath>
#include <Miscellaneous.h>

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * �������� 1 �׵��ⵥԪ�Ļ��������䵼���Ķ��塣���в�ֵ���ڵ� \f$ i
   * \f$ �����ϵĻ������ı���ʽΪ
   *
   * \f[
   *    \vec{\psi}_i = |p_j - p_k| \left(\lambda_j (x) \nabla
   *    \lambda_k (x) - \lambda_k (x) \nabla \lambda_j (x))
   * \f]
   *
   * ���� \f$ \lambda_i \f$ �ǵ� \f$ i \f$ ��������꺯����\f$ j = (i
   * + 1)%3 \f$, \f$ k = (i + 2)%3 \f$��
   * 
   * ��Ϊ����ļ����õ����� nVector<2,double> ������ʹ�� C++������
   * ���б��롣
   *
   */

#define vector_length 2
#define vt nVector<vector_length,double>

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
#define dlambda0_dx ((v[1][1] - v[2][1])/area)
#define dlambda0_dy ((v[2][0] - v[1][0])/area)
#define dlambda1_dx ((v[2][1] - v[0][1])/area)
#define dlambda1_dy ((v[0][0] - v[2][0])/area)
#define dlambda2_dx ((v[0][1] - v[1][1])/area)
#define dlambda2_dy ((v[1][0] - v[0][0])/area)

void get_lambda(const double * p, const double ** v, double * lambda, double * area)
{
  area[0] = AREA(v[0], v[1], v[2]);
  lambda[0] = AREA(p, v[1], v[2])/area[0];
  lambda[1] = AREA(p, v[2], v[0])/area[0];
  lambda[2] = AREA(p, v[0], v[1])/area[0];
}

void psi_1(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3], area, l12;
  get_lambda(p, v, lambda, &area);
  GET_L(1, 2);
  val[0] = l*(lambda[1]*dlambda2_dx - lambda[2]*dlambda1_dx);
  val[1] = l*(lambda[1]*dlambda2_dy - lambda[2]*dlambda1_dy);
}

void psi_2(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3], area, l20;
  get_lambda(p, v, lambda, &area);
  GET_L(2, 0);
  val[0] = l*(lambda[2]*dlambda0_dx - lambda[0]*dlambda2_dx);
  val[1] = l*(lambda[2]*dlambda0_dy - lambda[0]*dlambda2_dy);
}

void psi_3(const double * p, const double ** v, void * value)
{
  vt& val = *((vt *)value);
  double lambda[3], area, l01;
  get_lambda(p, v, lambda, &area);
  GET_L(0, 1);
  val[0] = l*(lambda[0]*dlambda1_dx - lambda[1]*dlambda0_dx);
  val[1] = l*(lambda[0]*dlambda1_dy - lambda[1]*dlambda0_dy);
}

void gradient_psi_1(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area, l12;
  area = AREA(v[0], v[1], v[2]);
  GET_L(1, 2);
  val[0][0] = 0.0; val[0][1] = l*(dlambda1_dy*dlambda2_dx - dlambda2_dy*dlambda1_dx);
  val[1][0] = -val[0][1]; val[1][1] = 0.0;
}

void gradient_psi_2(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area, l20;
  area = AREA(v[0], v[1], v[2]);
  GET_L(2, 0);
  val[0][0] = 0.0; val[0][1] = l*(dlambda2_dy*dlambda0_dx - dlambda0_dy*dlambda2_dx);
  val[1][0] = -val[0][1]; val[1][1] = 0.0;
}

void gradient_psi_3(const double * p, const double ** v, void * value)
{
  vt * val = (vt *)value;
  double area, l01;
  area = AREA(v[0], v[1], v[2]);
  GET_L(0, 1);
  val[0][0] = 0.0; val[0][1] = l*(dlambda0_dy*dlambda1_dx - dlambda1_dy*dlambda0_dx);
  val[1][0] = -val[0][1]; val[1][1] = 0.0;
}
// by Hu Bin    
void curl_psi_1(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0],v[1],v[2]);
  GET_L(1, 2);
  val[0] = l/area;
}
  void curl_psi_2(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0],v[1],v[2]);
  GET_L(2, 0);
  val[0] = l/area;
}
  void curl_psi_3(const double * p, const double ** v, void * value)
{
  double * val = (double *)value;
  double area = AREA(v[0],v[1],v[2]);
  GET_L(0, 1);
  val[0] = l/area;
}

#ifdef __cplusplus
}
#endif

/*
 *  end of file
 **************************************************************************/
