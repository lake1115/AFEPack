#include "uiexp.h"
#include "parameter.h"
//#include <boost/math/special_functions/bessel.hpp>
//#include <boost/math/special_functions/bessel_prime.hpp>
//#include <boost/math/special_functions/hankel.hpp>

cvaltype get_h(int n)
{
  cvaltype h;
  double z = kappa*R;
  cvaltype H_n = boost::math::cyl_hankel_1(n,z);
  double j_prime,y_prime;
  j_prime = boost::math::cyl_bessel_j_prime(n,z);
  y_prime = boost::math::cyl_neumann_prime(n,z);
  cvaltype H_n_prime(j_prime,y_prime);
  h = z*H_n_prime/H_n;
  return h;
}

cvaltype uiExperiment::get_u_hat(int n,cvaltype cout)
{
  // /int_0^{2PI} f(x) = 2*PI/M*sum(f(x_i)) 
  double M = 0;
  double theta;
  cvaltype u_hat = 0;
  double N = 1.0*n;
  //cvaltype H_n = boost::math::cyl_hankel_1(0,kappa*R);
  const u_int& n_dof = fem_space.n_dof();
  for(int i = 0; i<n_dof; i++){
    int bm = fem_space.dofInfo(i).boundary_mark;
    if(bm != 2)
      continue;
    const Point<DIM> point = fem_space.dofInfo(i).interp_point;
    theta = atan2(point[1],point[0]);
    if(theta < 0)
      theta += 2*PI;
    //cvaltype u_ex = u_exact(point);
    //u_hat += u_ex*exp(-1.0*I*N*theta);
    u_hat += cout*exp(-1.0*I*N*theta);
    M++;
  }
  //u_hat = u_hat/M;
  return u_hat;
}


void uiExperiment::getMat()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
	  cvaltype cout = Jxw *(c_val*innerProduct(bas_grad[j][l],bas_grad[k][l])+(a_val*bas_val[j][l]*bas_val[k][l]));
	  triplets.push_back(T(ele_dof[j],ele_dof[k],cout));
	}
      }
      
    }
  }
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
}

void uiExperiment::getRhs()
{
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
  const u_int& n_dof = fem_space.n_dof();
  for(u_int i = 0;i < n_dof;i++){
    int bm = fem_space.dofInfo(i).boundary_mark;
    if(bm != bmark)
      continue;
    const Point<DIM> point = fem_space.dofInfo(i).interp_point;
    // cvaltype bndVal = bnd(point);
    solution(i) = bnd(point);
    rhs(i) = stiff_matrix.diagonal()[i] * solution(i);
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
	    cvaltype cout = Jxw*q_val*bas_val[j][l]*bas_val[k][l];
	    triplets.push_back(T(dgele_dof[j],dgele_dof[k],cout));
	  }
	}
      }
    }
  stiff_matrix.setZero();
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
  std::cout << "Neumman boundary condition ... OK!"<<std::endl;

};

void uiExperiment::TransparentBC(int bmark)
{
  // don't need to use datacache, because each boundary treat independent and only use once.
  buildDGFEMSpace(bmark);
  if(fem_space.n_DGElement() == 0)
    return;
  cvaltype H,u_hat;
  double theta,theta2;
  for(int n=-order;n<=order;n++){
    std::cout<<" n = "<<n<<std::endl;
    H = get_h(n);
    cvaltype u_hat2 = get_u_hat(n,1);
    double N = 1.0*n;
    DGFEMSpace<double,DIM>::DGElementIterator
      the_dgele = fem_space.beginDGElement(),
      end_dgele = fem_space.endDGElement();
    DGFEMSpace<double,DIM>::DGElementIterator
      the_dgele2 = fem_space.beginDGElement();
    for (u_int i = 0; the_dgele != end_dgele;++ the_dgele,++i) 
      {
	EdgeCache<double,DIM>& edgec = edge_cache[bmark_count][i];
	std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
	const int& n_q_pnt = edgec.n_quad_pnt;
	std::vector<std::vector<double> > bas_val = edgec.basis_value;
	const std::vector<int>& dgele_dof = edgec.p_neigh->dof();
	int n_dgele_dof = dgele_dof.size();
	
	for(int l=0;l<n_q_pnt;l++){
	  double Jxw = edgec.Jxw[l];
	  theta = atan2(q_pnt[l][1],q_pnt[l][0]);
	  if(theta<0)
	    theta += 2*PI;
	  for(int j=0;j<n_dgele_dof;j++){
	    u_hat = Jxw*bas_val[j][l]*exp(-1.0*I*N*theta)/2.0/PI;
	    the_dgele2 = fem_space.beginDGElement();
	    for(u_int t = 0; the_dgele2 !=end_dgele;++the_dgele2,++t){
	      EdgeCache<double,DIM>& edgec2 = edge_cache[bmark_count][t];
	      std::vector<Point<DIM> >& q_pnt2 = edgec2.q_pnt;
	      const int& n_q_pnt2 = edgec2.n_quad_pnt;
	      std::vector<std::vector<double> > bas_val2 = edgec2.basis_value;
	      const std::vector<int>& dgele_dof2 = edgec2.p_neigh->dof();
	      int n_dgele_dof2 = dgele_dof2.size();
	      for(int l=0;l<n_q_pnt2;l++){
		double Jxw = edgec2.Jxw[l];
		theta2 = atan2(q_pnt2[l][1],q_pnt2[l][0]);
 		if(theta2<0)
		  theta2 += 2*PI;
		for(int k=0;k<n_dgele_dof2;k++){
		  cvaltype cout = Jxw*exp(I*N*theta2)*bas_val2[k][l];
		  cvaltype value = -1/R*H*u_hat*cout;
		  // I don't know why I calculate twice
		  if(dgele_dof[j]!=dgele_dof2[k])
		     value = value /2.0;
		  triplets.push_back(T(dgele_dof[j],dgele_dof2[k],value));
		}
	      }
	    }	 
	  }
	}     	 
      }
  }
  
  stiff_matrix.setZero();
  stiff_matrix.setFromTriplets(triplets.begin(),triplets.end());
  std::cout << "Transparent boundary condition ... OK!"<<std::endl;
};


void uiExperiment::getExactVal()
{
  u_exact_re.reinit(fem_space);
  u_exact_im.reinit(fem_space);
  Error.reinit(fem_space);
  for(int i=0;i<fem_space.n_dof();i++){
    const Point<DIM> point = fem_space.dofInfo(i).interp_point;
    u_exact_re(i) = u_exact(point).real();
    u_exact_im(i) = u_exact(point).imag();
    Error(i) = std::abs(u_exact_re(i)+I*u_exact_im(i)-u_re(i)-I*u_im(i));
    
  }
}

void uiExperiment::getL1Error()
{
  double L1error = 1E-9;
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
      double df_value =std::abs(Jxw*( u_exact(q_pnt[l])-u_h_val));
      if(df_value > L1error){
	L1error = df_value;
      }
    }
  }
  
  std::cerr << "\n1 error = " << L1error << std::endl;
}

void uiExperiment::getError()
{  
  double L2error = 0;
  double L2error_grad = 0;

 FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (;the_ele != end_ele; ++ the_ele){
    const std::vector<int>& ele_dof = the_ele->dof();
    const int& ele_idx = the_ele->index();
    ElementCache<double,DIM>& ec = element_cache[ele_idx];
    std::vector<Point<DIM> >& q_pnt = ec.q_pnt;
    const int& n_q_pnt = ec.n_quad_pnt;
    std::vector<double> u_re_val = u_re.value(q_pnt,*the_ele);
    std::vector<double> u_im_val = u_im.value(q_pnt,*the_ele);
    std::vector<std::vector<double> > u_re_grad = u_re.gradient(q_pnt,*the_ele);
    std::vector<std::vector<double> > u_im_grad = u_im.gradient(q_pnt,*the_ele);
    cvaltype residual = 0.0;
    for (int l = 0; l < n_q_pnt; l++){
      double Jxw = ec.Jxw[l];
      cvaltype a_val = a(q_pnt[l]);
      cvaltype c_val = c(q_pnt[l]); 
      cvaltype u_h_val(u_re_val[l],u_im_val[l]);
      std::vector<cvaltype> u_h_grad(2);
      u_h_grad[0] = u_re_grad[l][0]+I*u_im_grad[l][0];
      u_h_grad[1] = u_re_grad[l][1]+I*u_im_grad[l][1];   
      cvaltype u_exact_val = u_exact(q_pnt[l]);
      //cvec_type u_exact_grad = u_exact_prime(q_pnt[l]);
      
      double df_value = std::norm(u_exact_val-u_h_val);
      //double df_grad = std::norm(u_exact_grad[0]-u_h_grad[0])+std::norm(u_exact_grad[1]-u_h_grad[1]);
      
      L2error += Jxw*df_value;
      // L2error_grad += Jxw*df_grad;
      
      residual += Jxw*(-c_val*(u_h_grad[0]*u_h_grad[0]+u_h_grad[1]*u_h_grad[1])+a_val*u_h_val);
    }
    ec.residual = std::abs(residual);
  }
  L2error = sqrt(fabs(L2error));
  //L2error_grad = sqrt(fabs(L2error_grad));
  
  std::cout << "\nL2 error = " << L2error << std::endl;
  //std::cout << "\nL2 gradient error = " << L2error_grad << std::endl;
}

void uiExperiment::getError_h()
{
  RegularMesh<DIM>& mesh = ir_mesh.regularMesh();
  u_int n_ele = mesh.n_geometry(DIM);
  double error_h = 0.0;
  for(u_int i=0;i<n_ele;i++){
    error_h += indicator[i]*indicator[i];
  }
  error_h = sqrt(error_h);
  std::cout << "\nL2 error_h = " << error_h << std::endl;
}

void uiExperiment::adaptMesh()
{
  std::cout<<"********** Begin Adapt Mesh **********"<<std::endl;
  mesh_adaptor.reinit(ir_mesh);
  mesh_adaptor.convergenceOrder() = 0;
  mesh_adaptor.refineStep() = 0;
  mesh_adaptor.setIndicator(indicator);
  mesh_adaptor.tolerence() = 0.0008;
  mesh_adaptor.adapt();
}
void uiExperiment::getIndicator()
{
  
  RegularMesh<DIM>& mesh = ir_mesh.regularMesh();
  indicator.reinit(mesh);

  int n_face = mesh.n_geometry(DIM - 1);
  std::vector<double> jump(n_face);

  for(u_int k=0;k<edge_cache[0].size();k++){
    EdgeCache<double,DIM>& edgec = edge_cache[0][k];
    std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
    Element<double,DIM>* ele1 = edgec.p_neigh;
    Element<double,DIM>* ele2 = edgec.p_neigh2;
    const int& n_q_pnt = edgec.n_quad_pnt;
    for(u_int l=0;l<n_q_pnt;l++){
      double Jxw = edgec.Jxw[l];
      std::vector<double> u_re_grad = u_re.gradient(q_pnt[l], *ele1);
      std::vector<double> u_im_grad = u_im.gradient(q_pnt[l], *ele1);
      std::vector<double> u_re_grad2 = u_re.gradient(q_pnt[l], *ele2);
      std::vector<double> u_im_grad2 = u_im.gradient(q_pnt[l], *ele2);
      std::vector<cvaltype> u_h_grad(2),u_h_grad2(2);
      u_h_grad[0] = u_re_grad[0]+I*u_im_grad[0];
      u_h_grad[1] = u_re_grad[1]+I*u_im_grad[1];      
      u_h_grad2[0] = u_re_grad2[0]+I*u_im_grad2[0];
      u_h_grad2[1] = u_re_grad2[1]+I*u_im_grad2[1];
      
      jump[edgec.idx] += Jxw*std::norm(((u_h_grad[0]*edgec.un[l][0]+u_h_grad[1]*edgec.un[l][1])-(u_h_grad2[0]*edgec.un[l][0]+u_h_grad2[1]*edgec.un[l][1])));
    }
    //  std::cout<<" i= "<< edgec.idx<< " indicator "<< jump[edgec.idx] <<std::endl; 
  }
  //in Neumann boundary
  for(u_int k=0;k<edge_cache[1].size();++k){
    EdgeCache<double,DIM>& edgec = edge_cache[1][k];
    std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
    Element<double,DIM>* ele = edgec.p_neigh;
    const int& n_q_pnt = edgec.n_quad_pnt;
    for(int l=0;l<n_q_pnt;l++){
      double Jxw = edgec.Jxw[l];
      std::vector<double> u_re_grad = u_re.gradient(q_pnt[l], *ele);
      std::vector<double> u_im_grad = u_im.gradient(q_pnt[l], *ele);
      std::vector<cvaltype> u_h_grad(2);
      u_h_grad[0] = u_re_grad[0]+I*u_im_grad[0];
      u_h_grad[1] = u_re_grad[1]+I*u_im_grad[1];      
      cvaltype g_val = g(q_pnt[l]);
      cvaltype a = (u_h_grad[0]*edgec.un[l][0]+u_h_grad[1]*edgec.un[l][1]);
      jump[edgec.idx] += Jxw*std::norm(2.0*(a-g_val));
    }
    //std::cout<<" i= "<< edgec.idx<< " indicator "<< jump[edgec.idx] <<std::endl; 
  }
  //in Transparent boundary
  for(u_int k=0;k<edge_cache[2].size();++k){
    EdgeCache<double,DIM>& edgec = edge_cache[2][k];
    std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
    Element<double,DIM>* ele = edgec.p_neigh;
    const int& n_q_pnt = edgec.n_quad_pnt;
    for(int l=0;l<n_q_pnt;l++){
      double Jxw = edgec.Jxw[l];
      double u_re_val = u_re.value(q_pnt[l],*ele);
      double u_im_val = u_im.value(q_pnt[l],*ele);
      std::vector<double> u_re_grad = u_re.gradient(q_pnt[l], *ele);
      std::vector<double> u_im_grad = u_im.gradient(q_pnt[l], *ele);
      std::vector<cvaltype> u_h_grad(2);
      u_h_grad[0] = u_re_grad[0]+I*u_im_grad[0];
      u_h_grad[1] = u_re_grad[1]+I*u_im_grad[1];
      cvaltype a = (u_h_grad[0]*edgec.un[l][0]+u_h_grad[1]*edgec.un[l][1]);
      cvaltype u_h_val(u_re_val,u_im_val);
      cvaltype T_val=0.0;
      for(int n=-order;n<=order;n++){
	cvaltype H = get_h(n);
	double N = 1.0*n;
	cvaltype u_hat = get_u_hat(n,u_h_val);
	double theta = atan2(q_pnt[l][1],q_pnt[l][0]);
	if(theta<0)
	  theta += 2*PI;
	T_val += 1/R*H*u_hat*exp(I*N*theta);
      }
      jump[edgec.idx] += Jxw*std::norm(2.0*(a-T_val));
    }
    //std::cout<<" i= "<< edgec.idx<< " indicator "<< jump[edgec.idx] <<std::endl; 
  }
	      
  // once more cycle for calculate indicator
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for(int i=0; the_ele!= end_ele; ++ the_ele, ++i)
    {
      GeometryBM& geo = the_ele->geometry();
      const int& ele_idx = the_ele->index();
      ElementCache<double,DIM>& ec = element_cache[ele_idx];
      // the residual
      double A = ec.diameter*ec.residual;
      double B = 0.0;
      for(int j=0; j<geo.n_boundary(); ++j){
	GeometryBM& bnd = mesh.geometry(DIM-1,geo.boundary(j));
	B += 0.5*ec.bnd_length_list[j]*jump[bnd.index()];
      }
      indicator[i] = sqrt(A*A + B);
      // std::cout<<" i= "<< the_ele->index()<< " indicator "<< indicator[i] <<std::endl; 
    } 

}
// read mesh file
uiExperiment::uiExperiment(const std::string& file)
{
  mesh_file = file;
  order = Order;
  n_bmark = 3;
};

uiExperiment::~uiExperiment()
{};

void uiExperiment::init()
{
  //change mesh data to local mesh data
  //h_tree.readMesh(mesh_file);
  h_tree.readEasyMesh(mesh_file);
  ir_mesh.reinit(h_tree);
  // ir_mesh.globalRefine(2);
    
  // read the triangle template
  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_coord_transform.readData("triangle.crd_trs");
  triangle_template_dof.reinit(triangle_template_geometry);
  triangle_template_dof.readData("triangle.1.tmp_dof");
  triangle_basis_function.reinit(triangle_template_dof);
  triangle_basis_function.readData("triangle.1.bas_fun");
  triangle_unit_out_normal.readData("triangle.out_nrm");
  // read the twin triangle template
  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  twin_triangle_coord_transform.readData("twin_triangle.crd_trs");
  twin_triangle_template_dof.reinit(twin_triangle_template_geometry);
  twin_triangle_template_dof.readData("twin_triangle.1.tmp_dof");
  twin_triangle_basis_function.reinit(twin_triangle_template_dof);
  twin_triangle_basis_function.readData("twin_triangle.1.bas_fun");
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
  updateGeometryCache(3);
  std::cout << "OK!" << std::endl;
  // for record Neumman bc bmark
  edge_cache = new std::vector<EdgeCache<double,DIM> >[n_bmark];
  // dg_dof = new std::vector<int>[n_bmark];
  buildDGFEMSpace(0);
  
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
  
  Mesh<DIM,DIM>& mesh = fem_space.mesh();
  // calculate Neumann boundary number
  int n_side = mesh.n_geometry(DIM-1);
  int n_dg_ele = 0;
  for(int i=0;i<n_side;i++){
    // here only for mark=2 
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark)
      n_dg_ele +=1;
  }

  //dg_dof[bmark_count].resize(n_dg_ele);
  // build DGElement
  fem_space.dgElement().resize(n_dg_ele);
  for(int i=0,j=0;i<n_side;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark){
      fem_space.dgElement(j).reinit(fem_space,i,0);
      // dg_dof[bmark_count][j] = i;
      j += 1;
    }
  }
  fem_space.buildDGElement();
  
  std::cout << "Update DG Data Cache ...";
  std::cout << " bmark: "<<bmark<<" ..."<< " side: "<< n_dg_ele<<" ...";
  if(bmark_count == n_bmark)
    std::cerr<<"Error: not enough bmark number"<<std::endl;
  updateDGGeometryCache(edge_cache[bmark_count],3);
  std::cout << "OK!" << std::endl;
}

void uiExperiment::solve()
{
  
  getMat();
  getRhs();
  
  NeummanBC(g,1);

  TransparentBC(2);
  
  //DirichletBC(bnd,2);
  //NeummanBC(g2,2);
  //DirichletBC(bnd,1);
  Eigen::BiCGSTAB<Eigen::SparseMatrix<cvaltype,Eigen::RowMajor> > Eigen_solver;

  Eigen_solver.compute(stiff_matrix);
  solution = Eigen_solver.solve(rhs);

  for(int i=0;i<fem_space.n_dof();i++){
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
    file << interp_point[0] << " " << interp_point[1] << " " << u_h(i) << std::endl;
  }
  file.close();

}

void uiExperiment::saveData()
{
  u_re.writeOpenDXData("u_re.dx");
  u_im.writeOpenDXData("u_im.dx");
  u_exact_re.writeOpenDXData("u_exact_re.dx");
  Error.writeOpenDXData("error.dx");
  writeMatlabData("u_r.dat",u_re);
  writeMatlabData("u_i.dat",u_im);
};


void uiExperiment::run()
{
  init();
  do{
    buildFEMSpace();
    solve();
    getExactVal();
    getError();  
    saveData();
    getIndicator();
    getError_h();
    adaptMesh();
      
    getchar();
  }while(1);
};

/**
 *  * end of file
 *
 */

