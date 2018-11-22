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
  Point<DIM> bc;
  double es;
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
  double ind;/// index of the element
  std::vector<Point<DIM> > bc_list;/// bary center list of all elements in the patch
  std::vector<double> es_list; /// element size list of all elements in the patch
  
  std::vector<std::vector<value_type> > basis_value;
  std::vector<std::vector<std::vector<value_type> > > basis_gradient;

  //ElementCache();
  double residual;
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
    bc[2] = (p0[2] + p1[2])/2.;
  }
  else if(n_vtx ==3){
    AFEPack::Point<DIM> p0 = mesh.point(geo.vertex(0));
    AFEPack::Point<DIM> p1 = mesh.point(geo.vertex(1));
    AFEPack::Point<DIM> p2 = mesh.point(geo.vertex(2));
    bc[0] = (p0[0] + p1[0] + p2[0])/3.;
    bc[1] = (p0[1] + p1[1] + p2[1])/3.;
    bc[2] = (p0[2] + p1[2] + p2[2])/3.;
  }
  else if(n_vtx == 4){
    AFEPack::Point<DIM> p0 = mesh.point(geo.vertex(0));
    AFEPack::Point<DIM> p1 = mesh.point(geo.vertex(1));
    AFEPack::Point<DIM> p2 = mesh.point(geo.vertex(2));
    AFEPack::Point<DIM> p3 = mesh.point(geo.vertex(3));
    bc[0] = (p0[0] + p1[0] + p2[0] + p3[0])/4.;
    bc[1] = (p0[1] + p1[1] + p2[1] + p3[1])/4.;
    bc[2] = (p0[2] + p1[2] + p2[2] + p3[2])/4.;   
  }
};

inline
void elesize(Mesh<DIM,DIM>& mesh,
		GeometryBM& geo,
		double& es)
{
  const u_int n_vtx = geo.n_vertex();
  if(n_vtx == 3){  
    Point<DIM>& p0 = mesh.point(geo.vertex(0));
    Point<DIM>& p1 = mesh.point(geo.vertex(1));     
    Point<DIM>& p2 = mesh.point(geo.vertex(2));     
    double es0 = sqrt((p1[1]-p0[1])*(p1[1]-p0[1])+(p1[0]-p0[0])*(p1[0]-p0[0]));
    double es1 = sqrt((p2[1]-p1[1])*(p2[1]-p1[1])+(p2[0]-p1[0])*(p2[0]-p1[0]));
    double es2 = sqrt((p0[1]-p2[1])*(p0[1]-p2[1])+(p0[0]-p2[0])*(p0[0]-p2[0])); 
    es = 2*(es0+es1+es2)/3.0*sqrt(3)/3.0;
  }
  else if(n_vtx == 4){
    AFEPack::Point<DIM> p0 = mesh.point(geo.vertex(0));
    AFEPack::Point<DIM> p1 = mesh.point(geo.vertex(1));
    AFEPack::Point<DIM> p2 = mesh.point(geo.vertex(2));
    AFEPack::Point<DIM> p3 = mesh.point(geo.vertex(3));
    double es0 = sqrt((p1[1]-p0[1])*(p1[1]-p0[1])+(p1[0]-p0[0])*(p1[0]-p0[0]));
    double es1 = sqrt((p2[1]-p0[1])*(p2[1]-p0[1])+(p2[0]-p0[0])*(p2[0]-p0[0]));
    double es2 = sqrt((p3[1]-p0[1])*(p3[1]-p0[1])+(p3[0]-p0[0])*(p3[0]-p0[0]));
    double es3 = sqrt((p2[1]-p1[1])*(p2[1]-p1[1])+(p2[0]-p1[0])*(p2[0]-p1[0]));
    double es4 = sqrt((p3[1]-p1[1])*(p3[1]-p1[1])+(p3[0]-p1[0])*(p3[0]-p1[0]));
    double es5 = sqrt((p3[1]-p2[1])*(p3[1]-p2[1])+(p3[0]-p2[0])*(p3[0]-p2[0]));
    es = 2*(es0+es1+es2+es3+es4+es5)/6.0*sqrt(6)/4.0;
  }
};

#endif
