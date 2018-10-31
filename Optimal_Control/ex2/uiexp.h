#ifndef UIEXP_H
#define UIEXP_H


//#include <AFEPack/AMGSolver.h>       
#include <AFEPack/HGeometry.h>		 	
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/Operator.h>         
#include <AFEPack/Functional.h>       
#include <AFEPack/EasyMesh.h>
#include <AFEPack/MovingMesh2D.h>
//#include <AFEPack/SparseMatrixTool.h>

#include <Eigen/Sparse>

//#include <lac/vector.h>
//#include <lac/vector_memory.h>
//#include <lac/full_matrix.h>
//#include <lac/sparse_matrix.h>
//#include <lac/sparsity_pattern.h>
//#include <lac/compressed_sparsity_pattern.h>
//#include <lac/solver_cg.h>
//#include <lac/precondition.h>
//#include <lac/solver_bicgstab.h>
//#include <lac/sparse_ilu.h>



#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "typedefs.h"
#include "datacache.h"

typedef Eigen::Triplet<cvaltype> T;

class uiExperiment 
{
public:

 private:
 HGeometryTree<DIM> h_tree;
 IrregularMesh<DIM> ir_mesh_y;   // mesh for y and p
 IrregularMesh<DIM> ir_mesh_u; // mesh for u
 
 std::string mesh_file;		///mesh file
 
 TemplateGeometry<DIM> triangle_template_geometry;
 CoordTransform<DIM,DIM> triangle_coord_transform;
 
 TemplateDOF<DIM> triangle_template_dof_y;
 BasisFunctionAdmin<double,DIM,DIM> triangle_basis_function_y;
 UnitOutNormal<DIM> triangle_unit_out_normal_y;
 TemplateDOF<DIM> triangle_template_dof_u;
 BasisFunctionAdmin<double,DIM,DIM> triangle_basis_function_u;
 UnitOutNormal<DIM> triangle_unit_out_normal_u;

  
 TemplateGeometry<DIM> twin_triangle_template_geometry;
 CoordTransform<DIM,DIM> twin_triangle_coord_transform;
 
 TemplateDOF<DIM> twin_triangle_template_dof_y;
 BasisFunctionAdmin<double,DIM,DIM> twin_triangle_basis_function_y;
 UnitOutNormal<DIM> twin_triangle_unit_out_normal_y;
 TemplateDOF<DIM> twin_triangle_template_dof_u;
 BasisFunctionAdmin<double,DIM,DIM> twin_triangle_basis_function_u;
 UnitOutNormal<DIM> twin_triangle_unit_out_normal_u;

 TemplateGeometry<DIM-1> interval_template_geometry;
 CoordTransform<DIM-1,DIM> interval_to2d_coord_transform;
 
 std::vector<TemplateElement<double,DIM,DIM>> template_element;
 std::vector<TemplateDGElement<DIM-1,DIM>> dg_template_element;

 DGFEMSpace<double,DIM> fem_space_y;   /// FE space for y and p
 DGFEMSpace<double,DIM> fem_space_u; /// FE space for u 

 std::vector<ElementCache<double,DIM> > element_cache_y;
 std::vector<ElementCache<double,DIM> > element_cache_u;
 std::vector<EdgeCache<double,DIM> >* edge_cache;

 FEMFunction<double,DIM> y_h;    
 FEMFunction<double,DIM> p_h;    
 FEMFunction<double,DIM> u_h;    

 FEMFunction<double,DIM> y_exact;
 
 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix_y;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution_y;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs_y;
 
 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix_p;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution_p;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs_p;
 
 std::vector<T> triplets_y;
 std::vector<T> triplets_p;
 
 Indicator<DIM> indicator_y;
 Indicator<DIM> indicator_u;
 
 MeshAdaptor<DIM> mesh_adaptor;
 
 double sys_t0;		///system start time
 double dt;			/// time step

 Eigen::Matrix2d A;
 
 u_int n_bmark;
 u_int bmark_count;
 std::vector<int> bmark_list;
 
 public:
 uiExperiment(const std::string& file);

 virtual ~uiExperiment();

 // initial FEMspace
 void init();
 void buildFEMSpace();
 void buildDGFEMSpace(int bmark=11);
 void updateGeometryCache(DGFEMSpace<double,DIM>& fem_space,std::vector<ElementCache<double,DIM> >& element_cache,u_int alg_acc=3);
 void updateDGGeometryCache(DGFEMSpace<double,DIM>& fem_space,std::vector<EdgeCache<double,DIM> >& edge_cache, u_int alg_acc=3);
 void DirichletBC(Eigen::SparseMatrix<cvaltype,Eigen::RowMajor>& stiff_matrix,Eigen::Matrix<cvaltype,Eigen::Dynamic,1>& rhs,CFunc bnd,int bmark=1);
 void NeummanBC(CFunc g,int bmark=11);
 void getRhs_y();
 void getRhs_p();
 void getMat_y();
 void getMat_p();
 void getRhs_exact_y();
 void getMat_exact_y();
 void get_exact_y();
 double project_u(Point<DIM> p);
 double project_y(Point<DIM> p);
 void solve();
 void getError();
 void adaptMesh();
 void getIndicator();
 virtual void saveData();
 void run();
};

inline
bool onElement(Point<DIM> pntA,Point<DIM> pntB,Point<DIM> pntC,Point<DIM> pntP){
  std::vector<double> v0(DIM),v1(DIM),v2(DIM);
  v0[0] = pntC[0]-pntA[0];
  v0[1] = pntC[1]-pntA[1];
  v1[0] = pntB[0]-pntA[0];
  v1[1] = pntB[1]-pntA[1];
  v2[0] = pntP[0]-pntA[0];
  v2[1] = pntP[1]-pntA[1];

  double dot00 = v0[0]*v0[0]+v0[1]*v0[1];
  double dot01 = v0[0]*v1[0]+v0[1]*v1[1];
  double dot02 = v0[0]*v2[0]+v0[1]*v2[1];
  double dot11 = v1[0]*v1[0]+v1[1]*v1[1];
  double dot12 = v1[0]*v2[0]+v1[1]*v2[1];

  double u = (dot11*dot02-dot01*dot12)/(dot00*dot11-dot01*dot01);
  if(u < 0 || u > 1)
    return false;
  double v = (dot00*dot12-dot01*dot02)/(dot00*dot11-dot01*dot01);
  if(v < 0 || v > 1)
    return false;
  return u + v <= 1;
};

inline
double getDiameter(Point<DIM> p0,Point<DIM> p1,Point<DIM> p2){
  
  double l1 = sqrt((p1[1]-p0[1])*(p1[1]-p0[1])+(p1[0]-p0[0])*(p1[0]-p0[0]));
  double l2 = sqrt((p2[1]-p1[1])*(p2[1]-p1[1])+(p2[0]-p1[0])*(p2[0]-p1[0]));
  double l3 = sqrt((p0[1]-p2[1])*(p0[1]-p2[1])+(p0[0]-p2[0])*(p0[0]-p2[0]));
  double p = (l1+l2+l3)/2;
  return l1*l2*l3/(2*p*(p-l1)*(p-l2)*(p-l3));
}

#endif 
