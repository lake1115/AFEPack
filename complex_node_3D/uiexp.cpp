#include "uiexp.h"
#include "parameter.h"


void uiExperiment::getMat()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    std::vector<std::vector<std::vector<double> > >& bas_grad = ec.basis_gradient; 
    for(int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      cvaltype a_val = a(q_pnt[l]);
      cvaltype c_val = c(q_pnt[l]);
      for(int j = 0; j < n_ele_dof; j++){
	for(int k = 0; k < n_ele_dof; k++){
	  cvaltype cont = Jxw *(c_val*innerProduct(bas_grad[j][l],bas_grad[k][l])+(a_val*bas_val[j][l]*bas_val[k][l]));
	  triplets.push_back(T(ele_dof[j],ele_dof[k],cont));
	}
      }
      
    }
  }
}

void uiExperiment::getRhs()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    for (int l = 0;l < n_q_pnt;l ++) {
      double Jxw = ec.Jxw[l];
      cvaltype f_val = f(q_pnt[l]);
      for (int j = 0;j < n_ele_dof;j ++) {
	cvaltype cout = Jxw * bas_val[j][l] * f_val; 
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
    //solution(i) = bnd(point);
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

void uiExperiment::NeummanBC(CFunc g,int bmark)
{
  buildDGFEMSpace(bmark);
  if(fem_space->n_DGElement() == 0)
    return;
  DGFEMSpace<double,DIM>::DGElementIterator
    the_dgele = fem_space->beginDGElement(),
    end_dgele = fem_space->endDGElement();
  for (u_int i =0;the_dgele != end_dgele;++ the_dgele,++ i) 
    {
      //const u_int& edge_idx = the_dgele->index();
      EdgeCache<double,DIM>& edgec = edge_cache[bmark_count][i];
      std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
      const int& n_q_pnt = edgec.n_quad_pnt;
      std::vector<std::vector<double> > bas_val = edgec.basis_value;
      const std::vector<int>& dgele_dof = edgec.p_neigh->dof();
      const int& n_dgele_dof = dgele_dof.size();
      //std::vector<std::vector<double> > un = edgec.un;
      for(u_int l=0;l<n_q_pnt;l++){
	double Jxw = edgec.Jxw[l];
	cvaltype g_val = g(q_pnt[l]);
	cvaltype q_val = q(q_pnt[l]);
	//std::vector<double>& un_val = un[l];
	for(u_int j=0;j<n_dgele_dof;j++){
	  rhs(dgele_dof[j]) += -Jxw * bas_val[j][l]*g_val;
	   for(u_int k=0;k<n_dgele_dof;k++){
	    cvaltype cout = Jxw*q_val*bas_val[j][l]*bas_val[k][l];
	    triplets.push_back(T(dgele_dof[j],dgele_dof[k],cout));
	   }
	}
      }
    }
  std::cout << "Neumman boundary condition ... OK!"<<std::endl;

};


void uiExperiment::getError()
{  
  double L2error = 0;

  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<double> u_re_val = u_re.value(q_pnt,*the_ele);
    std::vector<double> u_im_val = u_im.value(q_pnt,*the_ele);
    for (int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      cvaltype u_h_val(u_re_val[l],u_im_val[l]);
      cvaltype u_exact_val = u_exact(q_pnt[l]);
      double df_value = std::norm(u_exact_val-u_h_val);
      L2error += Jxw*df_value;
    }
  }
  L2error = sqrt(fabs(L2error));
  
  std::cerr << "\nL2 error = " << L2error << std::endl;
}

void uiExperiment::adaptMesh()
{
  std::cout<<"********** Begin Adapt Mesh **********"<<std::endl;
  old_ir_mesh = ir_mesh;
  ir_mesh = new IrregularMesh<DIM>(*old_ir_mesh);
  
  mesh_adaptor.reinit(*old_ir_mesh,*ir_mesh);
  mesh_adaptor.setIndicator(indicator);
  mesh_adaptor.tolerence() = 1.0e-4;
  mesh_adaptor.convergenceOrder() = 2;
  mesh_adaptor.refineStep() = 0;

  mesh_adaptor.adapt();

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
void uiExperiment::getIndicator()
{
  
  RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
  indicator.reinit(mesh);
#if 1
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for(;the_ele!=end_ele;++the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    for(int i = 0;i<n_dof;++i){
      Point<DIM> p = fem_space->dofInfo(i).interp_point;
      double norm = sqrt(u_re(ele_dof[i])*u_re(ele_dof[i])+u_im(ele_dof[i])*u_im(ele_dof[i]));
      indicator[ele_idx] += 1./sqrt(1.+pow(norm,2.0));
    }
  }
#else
  const int& n_face = mesh.n_geometry(DIM - 1);
  std::vector<bool> flag(n_face,false);
  std::vector<double> jump(n_face);
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for(;the_ele!=end_ele; ++ the_ele)
    {
      GeometryBM& geo = the_ele->geometry();
      for(int j=0; j<geo.n_boundary(); ++j)
	{
	  GeometryBM& bnd = mesh.geometry(DIM-1, geo.boundary(j));
	  AFEPack::Point<DIM>& p0 = mesh.point(bnd.vertex(0));
	  AFEPack::Point<DIM>& p1 = mesh.point(bnd.vertex(1));
	  AFEPack::Point<DIM>& p2 = mesh.point(bnd.vertex(2));
	  std::vector<double> u_h_grad = u_re.gradient(midpoint(p0,p1), *the_ele);
	  double a = (u_h_grad[0]*(p0[1]-p1[1])+u_h_grad[1]*(p1[0]-p0[0]));
	  if(flag[bnd.index()] == false)
	    {
	      jump[bnd.index()] = a;
	      flag[bnd.index()] = true;
	    }
	  else
	    {
	      jump[bnd.index()] -= a;
	      flag[bnd.index()] = false;
	    }
	}
    }
  the_ele = fem_space->beginElement();
  for(int i=0; the_ele!= end_ele; ++ the_ele, ++i)
    {
      GeometryBM& geo = the_ele->geometry();
      indicator[i] = 0.0;
      for(int j=0; j<geo.n_boundary(); ++j)
	{
	  GeometryBM& bnd = mesh.geometry(DIM-1,geo.boundary(j));
	  if(flag[bnd.index()])
	    continue;
	  indicator[i] += jump[bnd.index()]*jump[bnd.index()];
	}
    }
#endif

}
// read mesh file
uiExperiment::uiExperiment(const std::string& file)
{
  mesh_file = file;
  n_bmark = 3;
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
  tetrahedron_template_dof.readData("tetrahedron.1.tmp_dof");
  tetrahedron_basis_function.reinit(tetrahedron_template_dof);
  tetrahedron_basis_function.readData("tetrahedron.1.bas_fun");
  tetrahedron_unit_out_normal.readData("tetrahedron.out_nrm");
  // read the twin tetrahedron template
  twin_tetrahedron_template_geometry.readData("twin_tetrahedron.tmp_geo");
  twin_tetrahedron_coord_transform.readData("twin_tetrahedron.crd_trs");
  twin_tetrahedron_template_dof.reinit(twin_tetrahedron_template_geometry);
  twin_tetrahedron_template_dof.readData("twin_tetrahedron.1.tmp_dof");
  twin_tetrahedron_basis_function.reinit(twin_tetrahedron_template_dof);
  twin_tetrahedron_basis_function.readData("twin_tetrahedron.1.bas_fun");
  twin_tetrahedron_unit_out_normal.readData("twin_tetrahedron.out_nrm");  
  // read the four tetrahedron template
  four_tetrahedron_template_geometry.readData("four_tetrahedron.tmp_geo");
  four_tetrahedron_coord_transform.readData("four_tetrahedron.crd_trs");
  four_tetrahedron_template_dof.reinit(four_tetrahedron_template_geometry);
  four_tetrahedron_template_dof.readData("four_tetrahedron.1.tmp_dof");
  four_tetrahedron_basis_function.reinit(four_tetrahedron_template_dof);
  four_tetrahedron_basis_function.readData("four_tetrahedron.1.bas_fun");
  four_tetrahedron_unit_out_normal.readData("four_tetrahedron.out_nrm");

  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_to3d_coord_transform.readData("triangle.to3d.crd_trs");

  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  twin_triangle_to3d_coord_transform.readData("twin_triangle.to3d.crd_trs");
  
  // set template element
  template_element.resize(3);
  template_element[0].reinit(tetrahedron_template_geometry,
			     tetrahedron_template_dof,
			     tetrahedron_coord_transform,
			     tetrahedron_basis_function,
			     tetrahedron_unit_out_normal);
  
  template_element[1].reinit(twin_tetrahedron_template_geometry,
			     twin_tetrahedron_template_dof,
			     twin_tetrahedron_coord_transform,
			     twin_tetrahedron_basis_function,
			     twin_tetrahedron_unit_out_normal);
  
  template_element[2].reinit(four_tetrahedron_template_geometry,
			     four_tetrahedron_template_dof,
			     four_tetrahedron_coord_transform,
			     four_tetrahedron_basis_function,
			     four_tetrahedron_unit_out_normal);

  dg_template_element.resize(2);
  dg_template_element[0].reinit(triangle_template_geometry,triangle_to3d_coord_transform);
  dg_template_element[1].reinit(twin_triangle_template_geometry,twin_triangle_to3d_coord_transform);
  
  buildFEMSpace();

  std::cout << "********** Initialize Complete **********" << std::endl;
}


// build FEMSpace
void uiExperiment::buildFEMSpace()
{
  RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
  // renew FEMspace
  fem_space = new DGFEMSpace<double, DIM>(mesh,template_element,dg_template_element);

  const u_int& n_element = mesh.n_geometry(DIM);
  fem_space->element().resize(n_element);
  
#pragma omp parallel for
  for(u_int i = 0; i < n_element; i++){
    // if vertex is 4, use tetrahedron, 5 use twin tetrahedron and 7 use four tetrahedron
    const int& n_vtx = mesh.geometry(DIM,i).n_vertex();
    if(n_vtx == 4)
      fem_space->element(i).reinit(*fem_space,i,0);
    else if(n_vtx == 5)
      fem_space->element(i).reinit(*fem_space,i,1);
    else if(n_vtx == 7)
      fem_space->element(i).reinit(*fem_space,i,2);
    else
      {
	std::cerr<<"Error: the number of vertex is wrong with the 3D tetrahedron element"<<std::endl;
      }
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
  edge_cache = new std::vector<EdgeCache<double,DIM> >[n_bmark];
}
void uiExperiment::buildDGFEMSpace(int bmark)
{
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
  const u_int& n_edge = mesh.n_geometry(DIM-1);
  /// calculate dg_element number
  int n_dg_ele = 0;
  for(int i=0;i<n_edge;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark)
      n_dg_ele +=1;
  }

  // build DGElement
  fem_space->dgElement().resize(n_dg_ele);
  for(u_int i=0,j=0;i<n_edge;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark){
      const GeometryBM& edge_geo = mesh.geometry(DIM-1,i);
      const int& n_vtx = edge_geo.n_vertex();
      if(n_vtx == 3){
	fem_space->dgElement(j).reinit(*fem_space,i,0);
	j++;
      }
      else if (n_vtx == 4){
	fem_space->dgElement(j).reinit(*fem_space,i,1);
	j++;
      }
      else{
	std::cerr<<"Error: the number of vertex is wrong with the 2D triangle element"<<std::endl;
      }
    }
  }
  fem_space->buildDGElement();
  
  std::cout << "Update DG Data Cache ...";
  std::cout << " bmark: "<<bmark<<" ..."<< " boundary surface number: "<< n_dg_ele<<" ...";
  if(bmark_count == n_bmark)
    std::cerr<<"Error: not enough bmark number"<<std::endl;
  updateDGGeometryCache(edge_cache[bmark_count],3);
  std::cout << "OK!" << std::endl;

}

void uiExperiment::solve()
{
  
  getMat();
  getRhs();
  
  //NeummanBC(g,1); 
  //DirichletBC(bnd,2);
  
  //cube test
  NeummanBC(g2,2);
  
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

void writeMatlabData(const std::string& filename, FEMFunction<double,DIM> u_h)
{
  std::fstream file;
  file.open(filename, std::ios::in|std::ios::out|std::ios::trunc);
  const int& n_dof = u_h.size();
  const FEMSpace<double,DIM>& fem_space = u_h.femSpace();
  for(int i = 0;i < n_dof;++ i){
    const Point<DIM>& interp_point = fem_space.dofInfo(i).interp_point;
    file << interp_point[0] << " " << interp_point[1] << " " << interp_point[2] << " " << u_h(i) << std::endl;
  }
  file.close();

}

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
    
    getIndicator();
    adaptMesh();
      
    getchar();
  }while(1);
};

/**
 *  * end of file
 *
 */
