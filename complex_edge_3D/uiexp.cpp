#include "uiexp.h"
#include "parameter.h"


void uiExperiment::getMat()
{
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<vec_type,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<vec_type> >& bas_val = ec.basis_value;
    // here basis_gradient means basis_rotation
    std::vector<std::vector<std::vector<vec_type> > >& bas_rot = ec.basis_gradient; 
     for(int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      cvaltype a_val = a(q_pnt[l]);
      cvaltype c_val = c(q_pnt[l]);
      for(int j = 0; j < n_ele_dof; j++){
	for(int k = 0; k < n_ele_dof; k++){
	  cvaltype cont = Jxw * (c_val*innerProduct(bas_rot[j][l][0],bas_rot[k][l][0])+a_val*innerProduct(bas_val[j][l],bas_val[k][l]));
	  triplets.push_back(T(ele_dof[j],ele_dof[k],cont));
	}
      }
      
    }
  }
}

void uiExperiment::getRhs()
{
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<vec_type,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<vec_type> >& bas_val = ec.basis_value;
    for (int l = 0;l < n_q_pnt;l ++) {
      double Jxw = ec.Jxw[l];
      cvec_type f_val = f(q_pnt[l]);
      for (int j = 0;j < n_ele_dof;j ++) {
	cvaltype cout = Jxw * (bas_val[j][l][0] * f_val[0]+bas_val[j][l][1]*f_val[1]+bas_val[j][l][2]*f_val[2]); 
	rhs(ele_dof[j]) += cout;
      }
    }
  }
}
void uiExperiment::DirichletBC(CFunc bnd,int bmark)
{
  const u_int& n_dof = fem_space->n_dof();
  for(u_int i = 0;i < n_dof;i++){
    int bm = fem_space->dofInfo(i).boundary_mark;
    if(bm != bmark)
      continue;
    const Point<DIM> point = fem_space->dofInfo(i).interp_point;
    cvaltype bndVal = bnd(point);
    rhs(i) = stiff_matrix.diagonal()[i] * bndVal;
    for(Eigen::SparseMatrix<cvaltype,Eigen::RowMajor>::InnerIterator it(stiff_matrix,i);it;++it){
      if(it.col() != it.row()){
	it.valueRef() = 0;
	u_int k = it.col();
	u_int pb = stiff_matrix.outerIndexPtr()[k];
	u_int pe = stiff_matrix.outerIndexPtr()[k+1];
	const int *p = std::find(&stiff_matrix.data().index(pb),
				 &stiff_matrix.data().index(pe),it.row());
	if(p != &stiff_matrix.data().index(pe)){
	  int l = p - &stiff_matrix.data().index(stiff_matrix.outerIndexPtr()[0]);
	  rhs(k) -= stiff_matrix.data().value(l) * rhs(i) / stiff_matrix.diagonal()[i];
	  stiff_matrix.data().value(l) = 0;
	}
      }
    }
  }
  std::cout << "Dirichlet boundary condition ... OK!"<<std::endl;
}
void uiExperiment::NeummanBC(CvFunc g,int bmark)
{
  buildDGFEMSpace(bmark);
  if(fem_space->n_DGElement() == 0)
    return;
  DGFEMSpace<vec_type,DIM>::DGElementIterator
    the_dgele = fem_space->beginDGElement(),
    end_dgele = fem_space->endDGElement();
  for (u_int i=0;the_dgele != end_dgele;++ the_dgele,++i) 
    {
      EdgeCache<vec_type,DIM>& edgec = edge_cache[bmark_count][i];
      std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
      const int& n_q_pnt = edgec.n_quad_pnt;
      std::vector<std::vector<vec_type> > bas_val = edgec.basis_value;
      std::vector<std::vector<double> > un = edgec.un;

      const std::vector<int>& dgele_dof = edgec.p_neigh->dof();
      const int& n_dgele_dof = dgele_dof.size();

      for(u_int l=0;l<n_q_pnt;l++){
	double Jxw = edgec.Jxw[l];
	cvec_type g_val = g(q_pnt[l]);
	cvaltype q_val = q(q_pnt[l]);
	for(u_int j=0;j<n_dgele_dof;j++){
	  rhs(dgele_dof[j]) += -Jxw *( bas_val[j][l][0]*g_val[0]+bas_val[j][l][1]*g_val[1]+bas_val[j][l][2]*g_val[2]);
	  for(u_int k=0;k<n_dgele_dof;k++){
	    cvaltype cout = Jxw*q_val*((un[l][1]*bas_val[j][l][2]-un[l][2]*bas_val[j][l][1])*(un[l][1]*bas_val[k][l][2]-un[l][2]*bas_val[k][l][1])
				       +(un[l][2]*bas_val[j][l][0]-un[l][0]*bas_val[j][l][2])*(un[l][2]*bas_val[k][l][0]-un[l][0]*bas_val[k][l][2])
				       +(un[l][0]*bas_val[j][l][1]-un[l][1]*bas_val[j][l][0])*(un[l][0]*bas_val[k][l][1]-un[l][1]*bas_val[k][l][0]));
	    triplets.push_back(T(dgele_dof[j],dgele_dof[k],cout));
	  }
	}
      }
    }
  std::cout << "Neumman boundary condition ... OK!"<<std::endl;
};
void uiExperiment::getError()
{  
  double L2error = 0.;

  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& ele_idx = the_ele->index();
    ElementCache<vec_type,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<vec_type> u_re_val = u_re.value(q_pnt,*the_ele);
    std::vector<vec_type> u_im_val = u_im.value(q_pnt,*the_ele);
    for (int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      cvaltype u_h_val0(u_re_val[l][0],u_im_val[l][0]);
      cvaltype u_h_val1(u_re_val[l][1],u_im_val[l][1]);
      cvaltype u_h_val2(u_re_val[l][2],u_im_val[l][2]);
      cvec_type u_val = u_exact(q_pnt[l]);
      cvaltype u_val0 = u_val[0];
      cvaltype u_val1 = u_val[1];
      cvaltype u_val2 = u_val[2];
      double df_value0 = std::norm(u_val0-u_h_val0);
      double df_value1 = std::norm(u_val1-u_h_val1);
      double df_value2 = std::norm(u_val2-u_h_val2);
      L2error += Jxw*(df_value0+df_value1+df_value2);
    }
  }
  L2error = sqrt(fabs(L2error));
  std::cerr << "\nL2 error = " << L2error << std::endl;
}
void uiExperiment::refineMesh()
{
  std::cout<<"********** Begin Refine Mesh **********"<<std::endl;
  old_ir_mesh = ir_mesh;
  ir_mesh = new IrregularMesh<DIM>(*old_ir_mesh);
  
  ir_mesh->globalRefine(1);

  ir_mesh->semiregularize();
  ir_mesh->regularize(false);

  /// build the new space
  old_fem_space = fem_space;
  old_edge_cache = edge_cache;
  buildFEMSpace();
  /// release memory
  delete old_ir_mesh;
  delete old_fem_space;
  delete[] old_edge_cache;

}
  
// read mesh file
uiExperiment::uiExperiment(const std::string& file)
{
  mesh_file = file;
  n_bmark = 2;
 
};

uiExperiment::~uiExperiment()
{};

void uiExperiment::init()
{
  //change mesh data to local mesh data
  h_tree.readMesh(mesh_file);
  ir_mesh = new IrregularMesh<DIM>(h_tree);
  //ir_mesh->globalRefine(0);
  ir_mesh->semiregularize();
  ir_mesh->regularize(false);
    
  // read the tetrahedron template
  tetrahedron_template_geometry.readData("tetrahedron.tmp_geo");
  tetrahedron_coord_transform.readData("tetrahedron.crd_trs");
  tetrahedron_template_dof.reinit(tetrahedron_template_geometry);
  tetrahedron_template_dof.readData("tetrahedron.edge.1.tmp_dof");
  tetrahedron_basis_function.reinit(tetrahedron_template_dof);
  tetrahedron_basis_function.readData("tetrahedron.edge.1.bas_fun");
  tetrahedron_unit_out_normal.readData("tetrahedron.out_nrm");
  
  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_to3d_coord_transform.readData("triangle.to3d.crd_trs");
  
  // set template element
  template_element.resize(1);
  template_element[0].reinit(tetrahedron_template_geometry,
			     tetrahedron_template_dof,
			     tetrahedron_coord_transform,
			     tetrahedron_basis_function,
			     tetrahedron_unit_out_normal); 
 
  dg_template_element.resize(1);
  dg_template_element[0].reinit(triangle_template_geometry,triangle_to3d_coord_transform);

  buildFEMSpace();
  
  std::cout << "********** Initialize Complete **********" << std::endl;
}


// build FEMSpace
void uiExperiment::buildFEMSpace()
{
  RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
  // renew FEMspace
  fem_space = new DGFEMSpace<vec_type, DIM>(mesh,template_element,dg_template_element);

  const u_int& n_element = mesh.n_geometry(DIM);
  fem_space->element().resize(n_element);
  for(int i=0; i<n_element; i++){
    fem_space->element(i).reinit(*fem_space,i,0);
  }
  fem_space->buildElement();
  fem_space->buildDof();
  fem_space->buildDofBoundaryMark();
  
  u_re.reinit(*fem_space);
  u_im.reinit(*fem_space);

  const int& n_dof = fem_space->n_dof();
  stiff_matrix.resize(n_dof,n_dof);
  rhs.setZero(n_dof);
  solution.setZero(n_dof);
  triplets.clear();
  
  std::cout << "Update Data Cache ...";
  updateGeometryCache();
  std::cout << "OK!" << std::endl;
  // for record Neumman bc bmark
  edge_cache = new std::vector<EdgeCache<vec_type,DIM> >[n_bmark];
}
void uiExperiment::buildDGFEMSpace(int bmark)
{
    //find bmark
  std::vector<int>::iterator ret;
  ret = std::find(bmark_list.begin(),bmark_list.end(),bmark);
  if(ret == bmark_list.end()){
    bmark_count = bmark_list.size();
    bmark_list.push_back(bmark);
  }
  else{
    int pos = ret - bmark_list.begin();
    bmark_count = pos;
  }

  Mesh<DIM,DIM>& mesh = fem_space->mesh();
  // calculate Neumann boundary number
  const int& n_edge = mesh.n_geometry(DIM-1);
  int n_dg_ele = 0;
  for(int i=0;i<n_edge;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark)
      n_dg_ele ++;
  }

  // build DGElement
  fem_space->dgElement().resize(n_dg_ele);
  for(int i=0,j=0;i<n_edge;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark){
      fem_space->dgElement(j).reinit(*fem_space,i,0);
      j ++;
    }
  }
  fem_space->buildDGElement();
  
  std::cout << "Update DG Data Cache ...";
  std::cout << " bmark: "<<bmark<<" ...";
  if(bmark_count == n_bmark)
    std::cerr<<"Error: not enough bmark number"<<std::endl;
  updateDGGeometryCache(edge_cache[bmark_count]);
  std::cout << "OK!" << std::endl;

}

void uiExperiment::solve()
{
  
  getMat();
  getRhs();

  NeummanBC(g,2);
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
  DirichletBC(bnd,1);
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> > Eigen_solver;

  Eigen_solver.compute(stiff_matrix);
  solution = Eigen_solver.solve(rhs);

  for(int i=0;i<fem_space->n_dof();i++){
    cvaltype value = solution(i);
    u_re(i) = value.real();
    u_im(i) = value.imag();
  }

  

};
/*
void writeMatlabData(const std::string& filename, FEMFunction<double,DIM> u_h)
{
  std::fstream file;
  file.open(filename, std::ios::in|std::ios::out|std::ios::trunc);
  const int& n_dof = u_h.size();
  const FEMSpace<double,DIM>& fem_space = u_h.femSpace();
  for(int i = 0;i < n_dof;++ i){
    const Point<DIM>& interp_point = fem_space.dofInfo(i).interp_point;
    file << interp_point[0] << " " << interp_point[1] << " " << u_h(i) << std::endl;
  }
  file.close();

}
*/
void uiExperiment::saveData()
{
  u_re.writeOpenDXData("u_re.dx");
  u_im.writeOpenDXData("u_im.dx");
  //writeMatlabData("u_r.dat",u_re);
  //writeMatlabData("u_i.dat",u_im);
};


void uiExperiment::run()
{
  init();
  do{
    solve();
    saveData();
    getError();
    
    std::cout<< " Refine? Enter or ctrl+C "<< std::endl;
    getchar();
    
    refineMesh();
  }while(1);
};

/**
 *  * end of file
 *
 */

