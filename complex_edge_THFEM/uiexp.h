#ifndef UIEXP_H
#define UIEXP_H


//#include "emdefs.h"

const int DIM = 2;
const int DOW = DIM;
const int TDIM = DIM;
const int vector_length = DIM; 

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

#include "emdefs.h"
#include "datacache.h"

typedef Eigen::Triplet<cvaltype> T;

class uiExperiment 
{
public:

 private:
 HGeometryTree<DIM> h_tree;
 IrregularMesh<DIM> ir_mesh;
 
 std::string mesh_file;		///mesh file

 TemplateGeometry<DIM> triangle_template_geometry;
 CoordTransform<DIM,DIM> triangle_coord_transform;
 TemplateDOF<DIM> triangle_template_dof;
 BasisFunctionAdmin<vec_type,DIM,DIM> triangle_basis_function;
 UnitOutNormal<DIM> triangle_unit_out_normal;

 TemplateGeometry<DIM> twin_triangle_template_geometry;
 CoordTransform<DIM,DIM> twin_triangle_coord_transform;
 TemplateDOF<DIM> twin_triangle_template_dof;
 BasisFunctionAdmin<vec_type,DIM,DIM> twin_triangle_basis_function;
 UnitOutNormal<DIM> twin_triangle_unit_out_normal;

 TemplateGeometry<DIM-1> interval_template_geometry;
 CoordTransform<DIM-1,DIM> interval_to2d_coord_transform;
 
 std::vector<TemplateElement<vec_type,DIM,DIM>> template_element;
 std::vector<TemplateDGElement<DIM-1,DIM> > dg_template_element;

 //FEMSpace<vec_type,DIM> fem_space;
 DGFEMSpace<vec_type,DIM> fem_space; /// finite element space

 std::vector<ElementCache<vec_type,DIM> > element_cache;
 std::vector<EdgeCache<double,DIM> > edge_cache;
 
 FEMFunction<vec_type,DIM> u_re;    
 FEMFunction<vec_type,DIM> u_im;    

 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs;

 std::vector<T> triplets;
 
 Indicator<DIM> indicator;
 MeshAdaptor<DIM> mesh_adaptor;
 
 double sys_t0 = 0;		///system start time
 double dt = 1;			/// time step
  
 public:
 uiExperiment(const std::string& file);

 virtual ~uiExperiment();

 // initial FEMspace
 void init();
 void buildFEMSpace();
 void buildDGFEMSpace(int bmark=2);
 void updateGeometryCache(u_int alg_acc=3);
 void updateDGGeometryCache(u_int alg_acc=3);
 void DirichletBC(CFunc bnd,int bmark=1);
 void NeummanBC(CvFunc g,int bmark=2);
 void getRhs();
 void getMat();
 void solve();
 void getError();
 void adaptMesh();
 void getIndicator();
 virtual void saveData();
 void run();
};


#endif 
