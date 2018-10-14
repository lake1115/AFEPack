#ifndef _DATACACHE_H_
#define _DATACACHE_H_

#include <lac/full_matrix.h>

template <class value_type>
struct SolutionCache
{
  std::vector<value_type> val;
  //std::vector<std::vector<value_type> > grad;

  
};
/*
struct BoundaryCache
{
  
};
*/
template <int DIM>
struct GeometryCache
{
  Point<DIM> bc; 
  int n_quad_pnt;
  std::vector<Point<DIM> > q_pnt;
  std::vector<double> Jxw;
};

template <typename value_type, int DIM> 
struct EdgeCache : public GeometryCache<DIM>
{
  std::vector<std::vector<value_type> > basis_value;
  std::vector<std::vector<std::vector<value_type> > > basis_gradient;
  std::vector<std::vector<value_type> > un;
};

template <typename value_type, int DIM>
struct ElementCache : public GeometryCache<DIM>
{
  double ind;/// index of the element
  std::vector<std::vector<value_type> > basis_value;
  std::vector<std::vector<std::vector<value_type> > > basis_gradient;

  //ElementCache();

  //SolutionCache the_val;/// cell average
};

#endif
