#include "uiexp.h"
#include "parameter.h"


void uiExperiment::getMat()
{
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
     ElementCache<vec_type,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<vec_type> >& bas_val = ec.basis_value;
    std::vector<std::vector<std::vector<vec_type> > >& bas_grad = ec.basis_gradient; 
     for(int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      cvaltype a_val = a(q_pnt[l]);
      for(int j = 0; j < n_ele_dof; j++){
	for(int k = 0; k < n_ele_dof; k++){
	  cvaltype cont = Jxw *((bas_grad[j][l][1][0]-bas_grad[j][l][0][1])*(bas_grad[k][l][1][0]-bas_grad[k][l][0][1])+(a_val*innerProduct(bas_val[j][l],bas_val[k][l])));
	  triplets.push_back(T(ele_dof[j],ele_dof[k],cont));
	}
      }
      
    }
  }
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
}

void uiExperiment::getRhs()
{
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
	cvaltype cout = Jxw * (bas_val[j][l][0] * f_val[0]+bas_val[j][l][1]*f_val[1]); 
	rhs(ele_dof[j]) += cout;
      }
    }
  }
}
void uiExperiment::DirichletBC(CFunc bnd,int bmark)
{
  const u_int& n_dof = fem_space.n_dof();
  for(u_int i = 0;i < n_dof;i++){
    int bm = fem_space.dofInfo(i).boundary_mark;
    if(bm != bmark)
      continue;
    const Point<DIM> point = fem_space.dofInfo(i).interp_point;
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
  // don't need to use datacache, because each boundary treat independent and only use once.
  buildDGFEMSpace(bmark);
  if(fem_space.n_DGElement() == 0)
    return;
  DGFEMSpace<vec_type,DIM>::DGElementIterator
    the_dgele = fem_space.beginDGElement(),
    end_dgele = fem_space.endDGElement();
  for (;the_dgele != end_dgele;++ the_dgele) 
    {
      double vol = the_dgele->templateElement().volume();
      const QuadratureInfo<DIM-1>& qi = the_dgele->findQuadratureInfo(3);
      int n_q_pnt = qi.n_quadraturePoint();
      std::vector<double> jac = the_dgele->local_to_global_jacobian(qi.quadraturePoint());
      std::vector<AFEPack::Point<DIM> > q_pnt = the_dgele->local_to_global(qi.quadraturePoint());
		
      Element<vec_type,DIM> *p_neigh0 = the_dgele->p_neighbourElement(0);
      Element<vec_type,DIM> *p_neigh1 = the_dgele->p_neighbourElement(1);
      // only one neighbourelement in boundary
      if(p_neigh1 != NULL){
	std::cerr<<"Error: Neumman boundary element only have one neighbour"<<std::endl;
	getchar();
      }

      std::vector<std::vector<double> > un =unitOutNormal(q_pnt,*p_neigh0,*the_dgele);
      
      std::vector<std::vector<vec_type> > bas_val = p_neigh0->basis_function_value(q_pnt);

      const std::vector<int>& dgele_dof = p_neigh0->dof();
      int n_dgele_dof = dgele_dof.size();

      for(int l=0;l<n_q_pnt;l++){
	double Jxw = vol*qi.weight(l)*jac[l];
	cvec_type g_val = g(q_pnt[l]);
	cvaltype q_val = q(q_pnt[l]);
	for(int j=0;j<n_dgele_dof;j++){
	  rhs(dgele_dof[j]) += -Jxw *( bas_val[j][l][0]*g_val[0]+bas_val[j][l][1]*g_val[1]);
	  for(int k=0;k<n_dgele_dof;k++){
	    cvaltype cout = Jxw*q_val*(un[l][0]*bas_val[j][l][1]-un[l][1]*bas_val[j][l][0])*(un[l][0]*bas_val[k][l][1]-un[l][1]*bas_val[k][l][0]);
	    triplets.push_back(T(dgele_dof[j],dgele_dof[k],cout));
	  }
	}
      }
    }
  stiff_matrix.setZero();
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
  std::cout << "Neumman boundary condition ... OK!"<<std::endl;
};

void uiExperiment::getError()
{  
  double L2error = 0.;

  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
      cvec_type u_val = u_exact(q_pnt[l]);
      cvaltype u_val0 = u_val[0];
      cvaltype u_val1 = u_val[1];
      double df_value0 = std::abs(u_val0-u_h_val0);
      double df_value1 = std::abs(u_val1-u_h_val1);
      L2error += Jxw*(df_value0*df_value0+df_value1*df_value1);
    }
  }
  L2error = sqrt(fabs(L2error));
  std::cerr << "\nL2 error = " << L2error << std::endl;
}

void uiExperiment::adaptMesh()
{
  std::cout<<"********** Begin Adapt Mesh **********"<<std::endl;
  mesh_adaptor.reinit(ir_mesh);
  mesh_adaptor.convergenceOrder() = 0;
  mesh_adaptor.refineStep() = 0;
  mesh_adaptor.setIndicator(indicator);
  mesh_adaptor.tolerence() = 5.0e-1;
  mesh_adaptor.adapt();
}
void uiExperiment::getIndicator()
{
  
  RegularMesh<DIM>& mesh = ir_mesh.regularMesh();
  indicator.reinit(mesh);
#if 1
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for(;the_ele!=end_ele;++the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    for(int i = 0;i<n_dof;++i){
      Point<DIM> p = fem_space.dofInfo(i).interp_point;
      double norm = sqrt(u_re(ele_dof[i])*u_re(ele_dof[i])+u_im(ele_dof[i])*u_im(ele_dof[i]));
      indicator[ele_idx] += 1./sqrt(1.+pow(norm,2.0));
    }
  }

#else
  int n_face = mesh.n_geometry(DIM - 1);
  std::vector<bool> flag(n_face,false);
  std::vector<double> jump(n_face);
  FEMSpace<vec_type,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for(;the_ele!=end_ele; ++ the_ele)
    {
      GeometryBM& geo = the_ele->geometry();
      for(int j=0; j<geo.n_boundary(); ++j)
	{
	  GeometryBM& bnd = mesh.geometry(DIM-1, geo.boundary(j));
	  AFEPack::Point<DIM>& p0 = mesh.point(bnd.vertex(0));
	  AFEPack::Point<DIM>& p1 = mesh.point(bnd.vertex(1));
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
  the_ele = fem_space.beginElement();
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
};

uiExperiment::~uiExperiment()
{};

void uiExperiment::init()
{
  //change mesh data to local mesh data
  //h_tree.readMesh(mesh_file);
  h_tree.readEasyMesh(mesh_file);
  ir_mesh.reinit(h_tree);
  ir_mesh.globalRefine(2);
    
  // read the triangle template
  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_coord_transform.readData("triangle.crd_trs");
  triangle_template_dof.reinit(triangle_template_geometry);
  triangle_template_dof.readData("triangle.edge.1.tmp_dof");
  triangle_basis_function.reinit(triangle_template_dof);
  triangle_basis_function.readData("triangle.edge.1.bas_fun");
  triangle_unit_out_normal.readData("triangle.out_nrm");
  // read the twin triangle template
  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  twin_triangle_coord_transform.readData("twin_triangle.crd_trs");
  twin_triangle_template_dof.reinit(twin_triangle_template_geometry);
  twin_triangle_template_dof.readData("twin_triangle.edge.1.tmp_dof");
  twin_triangle_basis_function.reinit(twin_triangle_template_dof);
  twin_triangle_basis_function.readData("twin_triangle.edge.1.bas_fun");
  twin_triangle_unit_out_normal.readData("twin_triangle.out_nrm");

  interval_template_geometry.readData("interval.tmp_geo");
  interval_to2d_coord_transform.readData("interval.to2d.crd_trs");
  
  // set template element
  template_element.resize(2);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof,
			     triangle_coord_transform,
			     triangle_basis_function,
			     triangle_unit_out_normal);
  
  template_element[1].reinit(twin_triangle_template_geometry,
			     twin_triangle_template_dof,
			     twin_triangle_coord_transform,
			     twin_triangle_basis_function,
			     twin_triangle_unit_out_normal);

  dg_template_element.resize(1);
  dg_template_element[0].reinit(interval_template_geometry,interval_to2d_coord_transform);

  
  std::cout << "********** Initialize Complete **********" << std::endl;
}


// build FEMSpace
void uiExperiment::buildFEMSpace()
{
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);
  
  RegularMesh<DIM>& mesh = ir_mesh.regularMesh();
  // renew FEMspace
  //fem_space.reinit(mesh,template_element);
  fem_space.reinit(mesh,template_element,dg_template_element);

  u_int n_element = mesh.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for(int i=0; i<n_element; i++){
    // if vertex is 3, use triangle, else use twin triangle 
    int n_vtx = mesh.geometry(DIM,i).n_vertex();
    if(n_vtx == 3)
      fem_space.element(i).reinit(fem_space,i,0);
    else if(n_vtx == 4)
      fem_space.element(i).reinit(fem_space,i,1);
    else{
      std::cerr<<"Error: the number of vertex is wrong"<<std::endl;
    }
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();
  
  u_re.reinit(fem_space);
  u_im.reinit(fem_space);

  const int& n_dof = fem_space.n_dof();
  stiff_matrix.resize(n_dof,n_dof);
  rhs.setZero(n_dof);
  solution.setZero(n_dof);
  triplets.clear();
  
  std::cout << "Update Data Cache ...";
  updateGeometryCache();
  std::cout << "OK!" << std::endl;

}
void uiExperiment::buildDGFEMSpace(int bmark)
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

  // build DGElement
  fem_space.dgElement().resize(n_dg_ele);
  for(int i=0,j=0;i<n_side;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark){
      fem_space.dgElement(j).reinit(fem_space,i,0);
      j += 1;
    }
  }
  fem_space.buildDGElement();
}

void uiExperiment::solve()
{
  
  getMat();
  getRhs();
  
  NeummanBC(g1,101);
  NeummanBC(g2,102);
  NeummanBC(g3,103);
  NeummanBC(g4,104);
  DirichletBC(bnd,1);
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> > Eigen_solver;

  Eigen_solver.compute(stiff_matrix);
  solution = Eigen_solver.solve(rhs);

  for(int i=0;i<fem_space.n_dof();i++){
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
    buildFEMSpace();
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

