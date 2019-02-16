#ifndef _DATACACHE_H_
#define _DATACACHE_H_

#include <lac/full_matrix.h>

template <class value_type>
struct SolutionCache
{
  std::vector<value_type> val;
  std::vector<std::vector<value_type> > grad;

  
};

template <int DIM>
struct GeometryCache
{
  double volume;
  Point<DIM> bc;  		/* the center of geometry */
  int n_quad_pnt;
  std::vector<Point<DIM> > q_pnt;
  std::vector<double> Jxw;
};

template <typename value_type, int DIM> 
struct EdgeCache : public GeometryCache<DIM>
{
  u_int idx;
  Element<value_type,DIM>* p_neigh;
  Element<value_type,DIM>* p_neigh2;
  std::vector<std::vector<value_type> > basis_value;
  std::vector<std::vector<std::vector<value_type> > > basis_gradient;
  std::vector<std::vector<value_type> > un;
};

template <typename value_type, int DIM>
struct ElementCache : public GeometryCache<DIM>
{
  u_int idx;/// index of the element
  std::vector<Point<DIM> > bc_list;/// bary center list of all elements in the patch
  std::vector<Point<DIM> > bnd_mid_pnt_list;
  std::vector<value_type> bnd_length_list;
  std::vector<std::vector<value_type> > bnd_un_list;
  double diameter;
  std::vector<std::vector<value_type> > basis_value;
  std::vector<std::vector<std::vector<value_type> > > basis_gradient;

  //ElementCache();

  double residual;/// residual
};



inline
void barycenter(Mesh<DIM,DIM>& mesh,
		GeometryBM& geo,
		Point<DIM>& bc)
{
  const u_int n_vtx = geo.n_vertex();
  if(n_vtx == 2){
    AFEPack::Point<DIM> p0 = mesh.point(geo.vertex(0));
    AFEPack::Point<DIM> p1 = mesh.point(geo.vertex(1));
    bc[0] = (p0[0] + p1[0])/2.;
    bc[1] = (p0[1] + p1[1])/2.;
  }
  else if(n_vtx ==3){
    AFEPack::Point<DIM> p0 = mesh.point(geo.vertex(0));
    AFEPack::Point<DIM> p1 = mesh.point(geo.vertex(1));
    AFEPack::Point<DIM> p2 = mesh.point(geo.vertex(2));
    bc[0] = (p0[0] + p1[0] + p2[0])/3.;
    bc[1] = (p0[1] + p1[1] + p2[1])/3.;

  }
};

#endif
