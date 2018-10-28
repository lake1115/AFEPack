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

 std::vector<ElementCache<double,DIM> > element_cache;
 std::vector<EdgeCache<double,DIM> >* edge_cache;

 FEMFunction<double,DIM> y_h;    
 FEMFunction<double,DIM> u_h;    

 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix_y;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution_y;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs_y;
 
 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix_u;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution_u;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs_u;
 
 std::vector<T> triplets_y;
 std::vector<T> triplets_u;
 
 Indicator<DIM> indicator_y;
 Indicator<DIM> indicator_u;
 
 MeshAdaptor<DIM> mesh_adaptor;
 
 double sys_t0 = 0;		///system start time
 double dt = 1;			/// time step

 u_int n_bmark = 4;
 u_int bmark_count = 0;
 std::vector<int> bmark_list;
 //SparsityPattern sp;
 //SparseMatrix<double> Ra;
 //SparseMatrix<double> Rb;
 //SparseMatrix<double> Ia;
 //SparseMatrix<double> Ib;
 

 //SparsityPattern F_sp;
 //SparseMatrix<double> F_Mat;
 //Vector<double> F_Rhs;
 
 
 public:
 uiExperiment(const std::string& file);

 virtual ~uiExperiment();

 // initial FEMspace
 void init();
 void buildFEMSpace();
 void buildDGFEMSpace(int bmark=11);
 void updateGeometryCache(DGFEMSpace<double,DIM> fem_space,std::vector<ElementCache<double,DIM> >& element_cache, u_int alg_acc=3);
 void updateDGGeometryCache(std::vector<EdgeCache<double,DIM> >& edge_cache, u_int alg_acc=3);
 void DirichletBC(CFunc bnd,int bmark=1);
 void NeummanBC(CFunc g,int bmark=11);
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
