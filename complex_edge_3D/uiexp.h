#ifndef UIEXP_H
#define UIEXP_H

const int DIM = 3;
const int DOW = DIM;
const int TDIM = DIM;
const int vector_length = DIM; 

#include <AFEPack/HGeometry.h>		 	
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/Operator.h>         
#include <AFEPack/Functional.h>       
#include <AFEPack/GmshMesh.h>
#include <AFEPack/MovingMesh3D.h>

#include <Eigen/Sparse>

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
 IrregularMesh<DIM>* ir_mesh;
 IrregularMesh<DIM>* old_ir_mesh;
 
 std::string mesh_file;		///mesh file

 TemplateGeometry<DIM> tetrahedron_template_geometry;
 CoordTransform<DIM,DIM> tetrahedron_coord_transform;
 TemplateDOF<DIM> tetrahedron_template_dof;
 BasisFunctionAdmin<vec_type,DIM,DIM> tetrahedron_basis_function;
 UnitOutNormal<DIM> tetrahedron_unit_out_normal;
 
 TemplateGeometry<DIM-1> triangle_template_geometry;
 CoordTransform<DIM-1,DIM> triangle_to3d_coord_transform;

 std::vector<TemplateElement<vec_type,DIM,DIM>> template_element;
 std::vector<TemplateDGElement<DIM-1,DIM> > dg_template_element;

 DGFEMSpace<vec_type,DIM>* fem_space; /// finite element space
 DGFEMSpace<vec_type,DIM>* old_fem_space; 

 std::vector<ElementCache<vec_type,DIM> > element_cache;
 std::vector<EdgeCache<vec_type,DIM> >* edge_cache;
 std::vector<EdgeCache<vec_type,DIM> >* old_edge_cache;
 
 FEMFunction<vec_type,DIM> u_re;    
 FEMFunction<vec_type,DIM> u_im;    

 Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> stiff_matrix;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> solution;
 Eigen::Matrix<cvaltype,Eigen::Dynamic,1> rhs;

 std::vector<T> triplets;
 
 Indicator<DIM> indicator;
 MeshAdaptor<DIM> mesh_adaptor;
 
 double sys_t0;		///system start time
 double dt;			/// time step

 u_int n_bmark;
 u_int bmark_count;
 std::vector<int> bmark_list;
 public:
 uiExperiment(const std::string& file);

 virtual ~uiExperiment();

 // initial FEMspace
 void init();
 void buildFEMSpace();
 void buildDGFEMSpace(int bmark=2);
 void updateGeometryCache(u_int alg_acc=3);
 void updateDGGeometryCache(std::vector<EdgeCache<vec_type,DIM> >& edge_cache,u_int alg_acc=3);
 void DirichletBC(CFunc bnd,int bmark=1);
 void NeummanBC(CvFunc g,int bmark=2);
 void getRhs();
 void getMat();
 void solve();
 void getError();
 void refineMesh();
 virtual void saveData();
 void run();
};


#endif 
