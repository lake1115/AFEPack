
#include "uiexp.h"
#include "datacache.h"

template <class value_type, int DIM, int DOW=DIM,int TDIM=DIM>
void updateElementGeometryInfo(Element<value_type,DIM>& ele, 
			       u_int alg_acc, 
			       ElementCache<value_type,DIM>& ec 
			       )
{
  
  double volume = ele.templateElement().volume(); 
  const QuadratureInfo<DIM>& quad_info = ele.findQuadratureInfo(alg_acc); 
  int& n_quad_pnt = ec.n_quad_pnt; 
  n_quad_pnt = quad_info.n_quadraturePoint();
  const std::vector<Point<DIM> >& quad_pnt = quad_info.quadraturePoint(); 
  std::vector<double> jacobian = ele.local_to_global_jacobian(quad_pnt);
  std::vector<Point<DIM> >& q_pnt = ec.q_pnt; 
  q_pnt = ele.local_to_global(quad_pnt);
  
  //ec.val.resize(ec.n_quad_pnt);
  ec.Jxw.resize(n_quad_pnt);
  for (int l = 0;l < n_quad_pnt;l ++) {
    double& Jxw = ec.Jxw[l];
    Jxw = volume*jacobian[l]*quad_info.weight(l);
  }
}
/*
template <class value_type, int DIM, int DOW=DIM,int TDIM=DIM>
void updateEdgeGeometryInfo(DGElement<value_type,DIM>& edge, 
			    u_int alg_acc, 
			    EdgeCache<value_type,DIM>& ec 
			    )
{
  Element<value_type,DIM>& ele = edge.neighbourElement(0);
  
  //  FEMSpace<value_type,DIM>& sp = ele.femSpace();
  //Mesh<DIM>& mesh = sp.mesh();
  //GeometryBM& geo = edge.geometry();
  
  double volume = edge.templateElement().volume(); 
  const QuadratureInfo<DIM-1>& quad_info = edge.findQuadratureInfo(alg_acc); 
  int& n_quad_pnt = ec.n_quad_pnt; 
  n_quad_pnt = quad_info.n_quadraturePoint();
  const std::vector<Point<DIM - 1> >& quad_pnt = quad_info.quadraturePoint();
  std::vector<double> jacobian = edge.local_to_global_jacobian(quad_pnt);
  std::vector<Point<DIM> >& q_pnt = ec.q_pnt; 

  q_pnt = edge.local_to_global(quad_pnt);
  //ec.un = unitOutNormal(q_pnt, ele, edge);
  
  ec.Jxw.resize(n_quad_pnt);
  // ec.val[0].resize(n_quad_pnt);
  //ec.val[1].resize(n_quad_pnt);
  
  for (int l = 0;l < n_quad_pnt;l ++) {
    double& Jxw = ec.Jxw[l];
    Jxw = volume*jacobian[l]*quad_info.weight(l); 
  }
}
*/
void uiExperiment::updateGeometryCache(u_int alg_acc)
{
  //u_int alg_acc = 3;
  
  u_int n_ele = fem_space.n_element();
  element_cache.clear();
  element_cache.resize(n_ele);
  
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for(;the_ele != end_ele;++ the_ele){
    Element<vec_type,DIM>& ele = *the_ele;
    const u_int& ele_idx = ele.index();
    ElementCache<vec_type,DIM>& ec = element_cache[ele_idx];    
    updateElementGeometryInfo(ele, alg_acc, ec);
    ec.basis_value = ele.basis_function_value(ec.q_pnt);
    ec.basis_gradient = ele.basis_function_gradient(ec.q_pnt);
  }
}
/*
void uiExperiment::BoundaryCache(int bmark)
{
  Mesh<DIM,DIM>& mesh = fem_space.mesh();
  // calculate Neumann boundary number
  int n_side = mesh.n_geometry(DIM-1);
  int n_dg_ele = 0;
  for(int i=0;i<n_side;i++){
    // here only for mark=2 
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark)
      n_dg_ele +=1;
  }

}
*/
/*
void uiExperiment::updateDGGeometryCache(u_int alg_acc)
{
  u_int n_side = fem_space.n_DGElement();
  edge_cache.clear();
  edge_cache.resize(n_side);
  DGFEMSpace<vec_type,DIM>::DGElementIterator the_dgele = fem_space.beginDGElement();
  DGFEMSpace<vec_type,DIM>::DGElementIterator end_dgele = fem_space.endDGElement();
  for(;the_dgele != end_dgele;++ the_dgele){
    DGElement<vec_type,DIM>& edge = *the_dgele;
    const u_int& edge_idx = edge.index();
    EdgeCache<vec_type,DIM>& ec = edge_cache[edge_idx];
    updateEdgeGeometryInfo(edge, alg_acc, ec); 
    Element<vec_type,DIM>* p_neigh = edge.p_neighbourElement(0);
    ec.basis_value = p_neigh->basis_function_value(ec.q_pnt);
    ec.basis_gradient = p_neigh->basis_function_gradient(ec.q_pnt);
    p_neigh = edge.p_neighbourElement(1); 
    if (p_neigh != NULL) {
      std::cerr<<"Error: edge elements only have one neighbourelement"<<std::endl;
    }
  }
}

*/
