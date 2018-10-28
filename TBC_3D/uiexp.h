#ifndef UIEXP_H
#define UIEXP_H


//#include "emdefs.h"

const int DIM = 3;
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
#include <AFEPack/GmshMesh.h>
#include <AFEPack/MovingMesh3D.h>
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

 TemplateGeometry<DIM> tetrahedron_template_geometry;
 CoordTransform<DIM,DIM> tetrahedron_coord_transform;
 TemplateDOF<DIM> tetrahedron_template_dof;
 BasisFunctionAdmin<double,DIM,DIM> tetrahedron_basis_function;
 UnitOutNormal<DIM> tetrahedron_unit_out_normal;

 TemplateGeometry<DIM> twin_tetrahedron_template_geometry;
 CoordTransform<DIM,DIM> twin_tetrahedron_coord_transform;
 TemplateDOF<DIM> twin_tetrahedron_template_dof;
 BasisFunctionAdmin<double,DIM,DIM> twin_tetrahedron_basis_function;
 UnitOutNormal<DIM> twin_tetrahedron_unit_out_normal;
 
 TemplateGeometry<DIM> four_tetrahedron_template_geometry;
 CoordTransform<DIM,DIM> four_tetrahedron_coord_transform;
 TemplateDOF<DIM> four_tetrahedron_template_dof;
 BasisFunctionAdmin<double,DIM,DIM> four_tetrahedron_basis_function;
 UnitOutNormal<DIM> four_tetrahedron_unit_out_normal;

 TemplateGeometry<DIM-1> triangle_template_geometry;
 CoordTransform<DIM-1,DIM> triangle_to3d_coord_transform;
 
 std::vector<TemplateElement<double,DIM,DIM>> template_element;
 std::vector<TemplateDGElement<DIM-1,DIM>> dg_template_element;
 
 DGFEMSpace<double,DIM> fem_space; /// finite element space

 std::vector<ElementCache<double,DIM> > element_cache;
 std::vector<EdgeCache<double,DIM> >* edge_cache;
 
 FEMFunction<double,DIM> u_re;    
 FEMFunction<double,DIM> u_im;    

 FEMFunction<double,DIM> u_exact_re;
 FEMFunction<double,DIM> u_exact_im;

 FEMFunction<double,DIM> Error;
 
 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs;

 std::vector<T> triplets;
 
 Indicator<DIM> indicator;
 MeshAdaptor<DIM> mesh_adaptor;
 
 double sys_t0 = 0;		///system start time
 double dt = 1;			/// time step
 int order;

 int n_bmark;
 int bmark_count;
 std::vector<int> bmark_list;
 
 cvaltype u_hat;
 
 public:
 uiExperiment(const std::string& file);

 virtual ~uiExperiment();

 // initial FEMspace
 void init();
 void buildFEMSpace();
 void buildDGFEMSpace(int bmark);
 void updateGeometryCache(u_int alg_acc=3);
 void updateDGGeometryCache(std::vector<EdgeCache<double,DIM> >& edge_cache, u_int alg_acc=3);
 void DirichletBC(CFunc bnd,int bmark=1);
 void NeummanBC(CFunc g,int bmark);
 void getRhs();
 void getMat();
 cvaltype get_u_hat(unsigned n, int m,cvaltype cout);
 void TransparentBC(int bmark);
 void solve();
 void getError();
 void getL1Error();
 void getError_h();
 void getExactVal();
 void adaptMesh();
 void getIndicator();
 virtual void saveData();
 void run();
};


#endif 
