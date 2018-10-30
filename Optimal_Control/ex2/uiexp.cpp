#include "uiexp.h"
#include "parameter.h"


double uiExperiment::project_u(Point<DIM> p){
  Mesh<DIM>& mesh_u = fem_space_u.mesh();
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_u.beginElement(),
    end_ele = fem_space_u.endElement();
  for(;the_ele!=end_ele;++the_ele)
    {
      GeometryBM& geo =the_ele->geometry();
      std::vector<int>& vtx = geo.vertex();
      Point<DIM>& pnt0 = mesh_u.point(vtx[0]); 
      Point<DIM>& pnt1 = mesh_u.point(vtx[1]); 
      Point<DIM>& pnt2 = mesh_u.point(vtx[2]); 
      if(onElement(pnt0,pnt1,pnt2,p))
	return u_h.value(p,*the_ele);
    }
  std::cerr<< " Can't find point in mesh of u! " << std::endl;
}

double uiExperiment::project_y(Point<DIM> p){
  Mesh<DIM>& mesh_y = fem_space_y.mesh();
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for(;the_ele!=end_ele;++the_ele)
    {
      GeometryBM& geo =the_ele->geometry();
      std::vector<int>& vtx = geo.vertex();
      Point<DIM>& pnt0 = mesh_y.point(vtx[0]); 
      Point<DIM>& pnt1 = mesh_y.point(vtx[1]); 
      Point<DIM>& pnt2 = mesh_y.point(vtx[2]); 
      if(onElement(pnt0,pnt1,pnt2,p)){
	
	return p_h.value(p,*the_ele);
      }
    }
  std::cerr<< " Can't find point in mesh of y and p! " << std::endl;
}

void uiExperiment::getMat_exact_y()
{
  stiff_matrix_y.setZero();
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    std::vector<std::vector<std::vector<double> > >& bas_grad = ec.basis_gradient; 
    for(int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      setA(A,q_pnt[l]);
      for(int j = 0; j < n_ele_dof; j++){
	for(int k = 0; k < n_ele_dof; k++){
	  std::vector<double> A_grad_y(2);
	  A_grad_y[0] = A(0,0)*bas_grad[j][l][0] + A(0,1)*bas_grad[j][l][1];
	  A_grad_y[1] = A(1,0)*bas_grad[j][l][0] + A(1,1)*bas_grad[j][l][1];
	  cvaltype u_val = 0.0;
	  //cvaltype u_val = u_exact(q_pnt[l]);
	  cvaltype cont = Jxw * (innerProduct(bas_grad[j][l],bas_grad[k][l]));
	  triplets_y.push_back(T(ele_dof[j],ele_dof[k],cont));
	}
      }
      
    }
  }
  stiff_matrix_y.setFromTriplets(triplets_y.begin(),triplets_y.end());
}

void uiExperiment::getRhs_exact_y()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    for (int l = 0;l < n_q_pnt;l ++) {
      double Jxw = ec.Jxw[l];
      cvaltype f_val = f(q_pnt[l]);
      cvaltype u_val = u_exact(q_pnt[l]);
      for (int j = 0;j < n_ele_dof;j ++) {
	cvaltype cout = Jxw * bas_val[j][l] *( f_val+u_val); 
	rhs_y(ele_dof[j]) += cout;
      }
    }
  }
}


void uiExperiment::getMat_y()
{
  stiff_matrix_y.setZero();
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    std::vector<std::vector<std::vector<double> > >& bas_grad = ec.basis_gradient; 
    for(int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      setA(A,q_pnt[l]);
      for(int j = 0; j < n_ele_dof; j++){
	for(int k = 0; k < n_ele_dof; k++){
	  std::vector<double> A_grad_y(2);
	  A_grad_y[0] = A(0,0)*bas_grad[j][l][0] + A(0,1)*bas_grad[j][l][1];
	  A_grad_y[1] = A(1,0)*bas_grad[j][l][0] + A(1,1)*bas_grad[j][l][1];
	  double u_val = 0.0;
	  //double u_val = project_u(q_pnt[l]);
	  cvaltype cont = Jxw *(innerProduct(A_grad_y,bas_grad[k][l])+u_val*bas_val[j][l]*bas_val[k][l]);
	  triplets_y.push_back(T(ele_dof[j],ele_dof[k],cont));
	}
      }
      
    }
  }
  stiff_matrix_y.setFromTriplets(triplets_y.begin(),triplets_y.end());

}

void uiExperiment::getMat_p()
{
  stiff_matrix_p.setZero();
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    std::vector<std::vector<std::vector<double> > >& bas_grad = ec.basis_gradient; 
    for(int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      setA(A,q_pnt[l]);
      for(int j = 0; j < n_ele_dof; j++){
	for(int k = 0; k < n_ele_dof; k++){
	  std::vector<double> A_grad_p(2);
	  A_grad_p[0] = A(0,0)*bas_grad[j][l][0] + A(0,1)*bas_grad[j][l][1];
	  A_grad_p[1] = A(1,0)*bas_grad[j][l][0] + A(1,1)*bas_grad[j][l][1];
	  double u_val = 0.0;
	  // double u_val = project_u(q_pnt[l]);
	  cvaltype cont = Jxw *(innerProduct(A_grad_p,bas_grad[k][l])+u_val*bas_val[j][l]*bas_val[k][l]);
	  triplets_p.push_back(T(ele_dof[j],ele_dof[k],cont));
	}
      }
      
    }
  }
  stiff_matrix_p.setFromTriplets(triplets_p.begin(),triplets_p.end());
}

void uiExperiment::getRhs_y()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    for (int l = 0;l < n_q_pnt;l ++) {
      double Jxw = ec.Jxw[l];
      cvaltype f_val = f(q_pnt[l]);
      double u_val = project_u(q_pnt[l]);
      for (int j = 0;j < n_ele_dof;j ++) {
	cvaltype cout = Jxw * bas_val[j][l] *( f_val+u_val); 
        //std::cout<<" c " << cout << std::endl;
	rhs_y(ele_dof[j]) += cout;
      }
    }
  }
}


void uiExperiment::getRhs_p()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& n_ele_dof = ele_dof.size();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<std::vector<double> >& bas_val = ec.basis_value;
    for (int l = 0;l < n_q_pnt;l ++) {
      double Jxw = ec.Jxw[l];
      double y_val = y_h.value(q_pnt[l],*the_ele);
      double y0_val = y_0(q_pnt[l]);
      //double y0_val = y_exact.value(q_pnt[l],*the_ele);
      for (int j = 0;j < n_ele_dof;j ++) {
	cvaltype cout = Jxw * bas_val[j][l] * (y_val - y0_val); 
	rhs_p(ele_dof[j]) += cout;
      }
    }
  }
}

void uiExperiment::DirichletBC(Eigen::SparseMatrix<cvaltype,Eigen::RowMajor>& stiff_matrix,Eigen::Matrix<cvaltype,Eigen::Dynamic,1>& rhs,CFunc bnd,int bmark) 
{
  const u_int& n_dof = fem_space_y.n_dof();
  for(u_int i = 0;i < n_dof;i++){
    int bm = fem_space_y.dofInfo(i).boundary_mark;
    if(bm != bmark)
      continue;
    const Point<DIM>& point = fem_space_y.dofInfo(i).interp_point;
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
/*
  void uiExperiment::NeummanBC(CFunc g,int bmark)
  {
  buildDGFEMSpace(bmark);
  if(fem_space.n_DGElement() == 0)
  return;
  DGFEMSpace<double,DIM>::DGElementIterator
  the_dgele = fem_space.beginDGElement(),
  end_dgele = fem_space.endDGElement();
  for (u_int i = 0;the_dgele != end_dgele;++ the_dgele,++ i) 
  {
  //const u_int& edge_idx = the_dgele->index();
  EdgeCache<double,DIM>& edgec = edge_cache[bmark_count][i];
  std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
  const int& n_q_pnt = edgec.n_quad_pnt;
  std::vector<std::vector<double> > bas_val = edgec.basis_value;
  const std::vector<int>& dgele_dof = edgec.p_neigh->dof();
  int n_dgele_dof = dgele_dof.size();
  for(int l=0;l<n_q_pnt;l++){
  double Jxw = edgec.Jxw[l];
  cvaltype g_val = g(q_pnt[l]);
  cvaltype q_val = q(q_pnt[l]);
  for(int j=0;j<n_dgele_dof;j++){
  rhs(dgele_dof[j]) += Jxw * bas_val[j][l]*g_val;
  for(int k=0;k<n_dgele_dof;k++){
  cvaltype cont = Jxw*q_val*bas_val[j][l]*bas_val[k][l];
  triplets.push_back(T(dgele_dof[j],dgele_dof[k],cont));
  }
  }
  }
  }
  stiff_matrix.setZero();
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
  std::cout << "Neumman boundary condition ... OK!"<<std::endl;
  };
*/
void uiExperiment::getError()
{  
  double L2error = 0;

  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space_y.beginElement(),
    end_ele = fem_space_y.endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache_y[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    for (int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      double y_h_val = y_h.value(q_pnt[l],*the_ele);
      double y_exact_val = Y_exact(q_pnt[l]);
      //double y_exact_val = y_exact.value(q_pnt[l],*the_ele);
      double df_value = std::norm(y_exact_val-y_h_val);
      L2error += Jxw*df_value;
    }
  }
  L2error = sqrt(L2error);
  std::cout << "\nL2 error = " << L2error << std::endl;
}
/*
  void uiExperiment::adaptMesh()
  {
  std::cout<<"********** Begin Adapt Mesh **********"<<std::endl;
  mesh_adaptor.reinit(ir_mesh);
  mesh_adaptor.convergenceOrder() = 0;
  mesh_adaptor.refineStep() = 0;
  mesh_adaptor.setIndicator(indicator);
  mesh_adaptor.tolerence() = 0.5;
  mesh_adaptor.adapt();
  }
  void uiExperiment::getIndicator()
  {
  
  RegularMesh<DIM>& mesh = ir_mesh.regularMesh();
  indicator.reinit(mesh);
  #if 1
  FEMSpace<double,DIM>::ElementIterator
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
  FEMSpace<double,DIM>::ElementIterator
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
*/
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
  ir_mesh_y.reinit(h_tree);
  ir_mesh_y.globalRefine(2);
    
  ir_mesh_u.reinit(h_tree);
  ir_mesh_u.globalRefine(2);
  // read the triangle template
  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_coord_transform.readData("triangle.crd_trs");
  
  triangle_template_dof_y.reinit(triangle_template_geometry);
  triangle_template_dof_y.readData("triangle.1.tmp_dof");
  triangle_basis_function_y.reinit(triangle_template_dof_y);
  triangle_basis_function_y.readData("triangle.1.bas_fun");
  triangle_unit_out_normal_y.readData("triangle.out_nrm");
  
  triangle_template_dof_u.reinit(triangle_template_geometry);
  triangle_template_dof_u.readData("triangle.1.D.tmp_dof");
  triangle_basis_function_u.reinit(triangle_template_dof_u);
  triangle_basis_function_u.readData("triangle.1.D.bas_fun");
  triangle_unit_out_normal_u.readData("triangle.out_nrm");
  // read the twin triangle template
  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  twin_triangle_coord_transform.readData("twin_triangle.crd_trs");
  
  twin_triangle_template_dof_y.reinit(twin_triangle_template_geometry);
  twin_triangle_template_dof_y.readData("twin_triangle.1.tmp_dof");
  twin_triangle_basis_function_y.reinit(twin_triangle_template_dof_y);
  twin_triangle_basis_function_y.readData("twin_triangle.1.bas_fun");
  twin_triangle_unit_out_normal_y.readData("twin_triangle.out_nrm");
  
  twin_triangle_template_dof_u.reinit(twin_triangle_template_geometry);
  twin_triangle_template_dof_u.readData("twin_triangle.1.D.tmp_dof");
  twin_triangle_basis_function_u.reinit(twin_triangle_template_dof_u);
  twin_triangle_basis_function_u.readData("twin_triangle.1.D.bas_fun");
  twin_triangle_unit_out_normal_u.readData("twin_triangle.out_nrm");
  
  interval_template_geometry.readData("interval.tmp_geo");
  interval_to2d_coord_transform.readData("interval.to2d.crd_trs");
  
  // set template element
  template_element.resize(4);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof_y,
			     triangle_coord_transform,
			     triangle_basis_function_y,
			     triangle_unit_out_normal_y);
  
  template_element[1].reinit(twin_triangle_template_geometry,
			     twin_triangle_template_dof_y,
			     twin_triangle_coord_transform,
			     twin_triangle_basis_function_y,
			     twin_triangle_unit_out_normal_y);
  
  template_element[2].reinit(triangle_template_geometry,
			     triangle_template_dof_u,
			     triangle_coord_transform,
			     triangle_basis_function_u,
			     triangle_unit_out_normal_u);
  
  template_element[3].reinit(twin_triangle_template_geometry,
			     twin_triangle_template_dof_u,
			     twin_triangle_coord_transform,
			     twin_triangle_basis_function_u,
			     twin_triangle_unit_out_normal_u);

  dg_template_element.resize(1);
  dg_template_element[0].reinit(interval_template_geometry,interval_to2d_coord_transform);
  
  
  std::cout << "********** Initialize Complete **********" << std::endl;
}


// build FEMSpace
void uiExperiment::buildFEMSpace()
{
  //for y and p FEMSpace
  ir_mesh_y.semiregularize();
  ir_mesh_y.regularize(false);
  
  RegularMesh<DIM>& mesh_y = ir_mesh_y.regularMesh();
  // renew FEMspace
  fem_space_y.reinit(mesh_y,template_element,dg_template_element);

  u_int n_element_y = mesh_y.n_geometry(DIM);
  fem_space_y.element().resize(n_element_y);
  for(int i=0; i<n_element_y; i++){
    // if vertex is 3, use triangle, else use twin triangle 
    int n_vtx = mesh_y.geometry(DIM,i).n_vertex();
    if(n_vtx == 3)
      fem_space_y.element(i).reinit(fem_space_y,i,0);
    else if(n_vtx == 4)
      fem_space_y.element(i).reinit(fem_space_y,i,1);
    else{
      std::cerr<<"Error: the number of vertex is wrong"<<std::endl;
    }
  }
  fem_space_y.buildElement();
  fem_space_y.buildDof();
  fem_space_y.buildDofBoundaryMark();
  //for u FEMSpace
  ir_mesh_u.semiregularize();
  ir_mesh_u.regularize(false);
  RegularMesh<DIM>& mesh_u = ir_mesh_u.regularMesh();
  // renew FEMspace
  fem_space_u.reinit(mesh_u,template_element,dg_template_element);

  u_int n_element_u = mesh_u.n_geometry(DIM);
  fem_space_u.element().resize(n_element_u);
  for(int i=0; i<n_element_u; i++){
    // if vertex is 3, use triangle, else use twin triangle 
    int n_vtx = mesh_u.geometry(DIM,i).n_vertex();
    if(n_vtx == 3)
      fem_space_u.element(i).reinit(fem_space_u,i,2);
    else if(n_vtx == 4)
      fem_space_u.element(i).reinit(fem_space_u,i,3);
    else{
      std::cerr<<"Error: the number of vertex is wrong"<<std::endl;
    }
  }
  fem_space_u.buildElement();
  fem_space_u.buildDof();
  fem_space_u.buildDofBoundaryMark();

  
  y_h.reinit(fem_space_y);
  p_h.reinit(fem_space_y);
  u_h.reinit(fem_space_u);

  y_exact.reinit(fem_space_y);
  
  const int& n_dof = fem_space_y.n_dof();
  stiff_matrix_y.resize(n_dof,n_dof);
  stiff_matrix_p.resize(n_dof,n_dof);
  rhs_y.setZero(n_dof);
  rhs_p.setZero(n_dof);
  solution_y.setZero(n_dof);
  solution_p.setZero(n_dof);
  triplets_y.clear();
  triplets_p.clear();
  
  std::cout << "Update Data Cache ...";
  updateGeometryCache(fem_space_y,element_cache_y);
  //updateGeometryCache(fem_space_u);
  std::cout << "OK!" << std::endl;

  //initial u value
  for(int i=0;i<fem_space_u.n_dof();i++){
    const Point<DIM>& point = fem_space_u.dofInfo(i).interp_point;
    u_h(i) = u_0(point);
  }
  for(int i=0;i<fem_space_y.n_dof();i++){
    const Point<DIM>& point = fem_space_y.dofInfo(i).interp_point;
    y_exact(i) = Y_exact(point);
  }

  
  // for record Neumman bc bmark
  //edge_cache = new std::vector<EdgeCache<double,DIM> >[n_bmark];
}
/*
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

  Mesh<DIM,DIM>& mesh = fem_space.mesh();
  // calculate Neumann boundary number
  int n_side = mesh.n_geometry(DIM-1);
  int n_dg_ele = 0;
  for(int i=0;i<n_side;i++){
  // begin from mark 11  
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
  
  std::cout << "Update DG Data Cache ...";
  std::cout << " bmark: "<<bmark<<" ...";
  if(bmark_count == n_bmark)
  std::cerr<<"Error: not enough bmark number"<<std::endl;
  updateDGGeometryCache(edge_cache[bmark_count],3);
  std::cout << "OK!" << std::endl;
  }
*/
void uiExperiment::get_exact_y(){
  getMat_exact_y();
  getRhs_exact_y();
  DirichletBC(stiff_matrix_y,rhs_y,bnd,1);
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> > Eigen_solver;
  Eigen_solver.compute(stiff_matrix_y);
  solution_y = Eigen_solver.solve(rhs_y);
  for(int i=0;i<fem_space_y.n_dof();i++){
    cvaltype value = solution_y(i);
    y_exact(i) = value.real();
  }
  //clear
  const int& n_dof = fem_space_y.n_dof();
  stiff_matrix_y.setZero();
  rhs_y.setZero(n_dof);
  solution_y.setZero(n_dof);
  triplets_y.clear();
}

void uiExperiment::solve()
{
  // get_exact_y();
  
  // y_exact.writeOpenDXData("y_exact.dx");
  double error;
  double rho = 0.1;
  double tolerence = 1.0e-5;
  do{
    //backup old u
    FEMFunction<double,DIM> old_u_h;
    old_u_h = u_h;
    // First: solve y 
    getMat_y();
    getRhs_y();

    DirichletBC(stiff_matrix_y,rhs_y,bnd,1);
  
    Eigen::BiCGSTAB<Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> > Eigen_solver;

    Eigen_solver.compute(stiff_matrix_y);
    solution_y = Eigen_solver.solve(rhs_y);
  
    for(int i=0;i<fem_space_y.n_dof();i++){
      cvaltype value = solution_y(i);
      y_h(i) = value.real();
    }
    // Second: solve p
    getMat_p();
    getRhs_p();

    DirichletBC(stiff_matrix_p,rhs_p,bnd,1);
    Eigen_solver.compute(stiff_matrix_p); 
    solution_p = Eigen_solver.solve(rhs_p);
  
    for(int i=0;i<fem_space_y.n_dof();i++){
      cvaltype value = solution_p(i);
      p_h(i) = value.real();
    }
    // Third: update u
    for(int i=0;i<fem_space_u.n_dof();i++){
      const Point<DIM>& point = fem_space_u.dofInfo(i).interp_point;
      double du = project_y(point);
      // std::cout<<" du is "<<du<<std::endl;
      u_h(i) = old_u_h(i) - rho*(old_u_h(i)-u_0(point)+du);
      u_h(i) = std::max(u_h(i),0.0);
    }
    // Forth: calculate error
    getError();
    old_u_h.add(-1,u_h);
    error = Functional::L2Norm(old_u_h,3);
    std::cout<<" error_u = " << error << std::endl;
    //getchar();
  }while(error>tolerence);
  
  std::cout<<"********** Finish Iterator *********" <<std::endl;
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
  u_h.writeOpenDXData("u_h.dx");
  y_h.writeOpenDXData("y_h.dx");
  y_exact.writeOpenDXData("y_exact.dx");
  p_h.writeOpenDXData("p_h.dx");
  // writeMatlabData("u_r.dat",u_re);
  //writeMatlabData("u_i.dat",u_im);
};


void uiExperiment::run()
{
  init();
  do{
    buildFEMSpace();
    solve();
    saveData();
    //getError();

    //getIndicator();
    //adaptMesh();
      
    getchar();
  }while(1);
};

/**
 *  * end of file
 *
 */

