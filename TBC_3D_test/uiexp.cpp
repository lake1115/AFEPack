#include "uiexp.h"
#include "parameter.h"
#include <ctime>
cvaltype get_Theta(int n)
{
  cvaltype Theta;
  double z = kappa*R;
  cvaltype H_n = boost::math::sph_hankel_1(n,z);
  double j_prime,y_prime;
  j_prime = boost::math::sph_bessel_prime(n,z);
  y_prime = boost::math::sph_neumann_prime(n,z);
  cvaltype H_n_prime(j_prime,y_prime);
  Theta = z*H_n_prime/H_n;
  return Theta;
}

cvaltype uiExperiment::get_u_hat(int n, int m)
{
  
  double M = dg_dof.size();
  double theta,phi;
  cvaltype u_hat = 0;
  double N = 1.0*n;
  std::cout << "M = "<<M <<std::endl;
  for(int i = 0; i < M; i++){
    const Point<DIM> point = fem_space.dofInfo(dg_dof[i]).interp_point;
    theta = acos(point[2]/R);
    phi = atan2(point[1],point[0]);
    if(phi < 0)
      phi += 2*PI;
    cvaltype Y = boost::math::spherical_harmonic(n,m,theta,phi);
    cvaltype cout = u_exact(point);
    //cvaltype cout = u_re(dg_dof[i]) + I*u_im(dg_dof[i]);
    u_hat += cout*std::conj(Y);  
  }
  u_hat = 4*PI*u_hat/M;
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
	  cvaltype cont = Jxw *(c_val*innerProduct(bas_grad[j][l],bas_grad[k][l])+(a_val*bas_val[j][l]*bas_val[k][l]));
	  triplets.push_back(T(ele_dof[j],ele_dof[k],cont));
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
  for (u_int i =0;the_dgele != end_dgele;++ the_dgele,++ i) 
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
  buildDGFEMSpace(bmark);
  if(fem_space.n_DGElement() == 0)
    return;
  double M = 0;
  for(int i=0;i<fem_space.n_dof(); i++){
    int bm = fem_space.dofInfo(i).boundary_mark;
    if(bm != 2)
      continue;
    M++;
  }
  std::cout<<" M ="<<M<<std::endl;
  dg_dof.resize(M);
  for(int i = 0,j = 0; i < fem_space.n_dof(); i++){
    int bm = fem_space.dofInfo(i).boundary_mark;
    if(bm != 2)
      continue;
    dg_dof[j] = i;
    j++;
  }

  cvaltype Theta,u_hat;
  double theta,phi;
  double theta2,phi2;
  cvaltype Y,Y2;
  clock_t start,end;
  start = clock();
  for(int n = 0; n <= order; n++){
    Theta = get_Theta(n);
    for(int m = -n; m <= n; m++){
      std::cout<<" n = "<<n <<" m = "<<m<< std::endl;
      // u_hat = get_u_hat(n,m);
      /*
	for(int i=0; i<M; i++){
	const Point<DIM> point = fem_space.dofInfo(dg_dof[i]).interp_point;
	theta = acos(point[2]/R);
	phi = atan2(point[1],point[0]);
	if(phi < 0)
	phi += 2*PI;
	Y = boost::math::spherical_harmonic(n,m,theta,phi); 
	u_hat = 4*PI/M*std::conj(Y);
      */

      DGFEMSpace<double,DIM>::DGElementIterator
	the_dgele = fem_space.beginDGElement(),
	end_dgele = fem_space.endDGElement();

      DGFEMSpace<double,DIM>::DGElementIterator
	the_dgele2 = fem_space.beginDGElement();
      //////////////////////////////// 
      for (u_int i = 0; the_dgele != end_dgele; ++ the_dgele,++i)
	{
	  EdgeCache<double,DIM>& edgec = edge_cache[bmark_count][i];
	  std::vector<Point<DIM> >& q_pnt = edgec.q_pnt;
	  const int& n_q_pnt = edgec.n_quad_pnt;
	  std::vector<std::vector<double> > bas_val = edgec.basis_value;
	  const std::vector<int>& dgele_dof = edgec.p_neigh->dof();
	  int n_dgele_dof = dgele_dof.size();
	  for(int l=0;l<n_q_pnt;l++){
	    double Jxw = edgec.Jxw[l];
	    theta = acos(q_pnt[l][2]/R);
	    phi = atan2(q_pnt[l][1],q_pnt[l][0]);
	    if(phi < 0)
	      phi += 2*PI;
	    cvaltype Y = boost::math::spherical_harmonic(n,m,theta,phi);
	    for(int j=0;j<n_dgele_dof;j++){
	      u_hat = Jxw*bas_val[j][l]*std::conj(Y);
	      ///////////////////////////////  
	      the_dgele2 = fem_space.beginDGElement();
	      for (u_int t = 0;the_dgele2 != end_dgele;++the_dgele2,++t) 
		{
		  EdgeCache<double,DIM>& edgec2 = edge_cache[bmark_count][t];
		  std::vector<Point<DIM> >& q_pnt2 = edgec2.q_pnt;
		  const int& n_q_pnt2 = edgec2.n_quad_pnt;
		  std::vector<std::vector<double> > bas_val2 = edgec2.basis_value;
		  const std::vector<int>& dgele_dof2 = edgec2.p_neigh->dof();
		  int n_dgele_dof2 = dgele_dof2.size();

		  for(int l2=0;l2<n_q_pnt2;l2++){
		    double Jxw2 = edgec2.Jxw[l2];
		    theta2 = acos(q_pnt2[l2][2]/R);
		    phi2 = atan2(q_pnt2[l2][1],q_pnt2[l2][0]);
		    if(phi2 < 0)
		      phi2 += 2*PI;
	  
		    Y2 = boost::math::spherical_harmonic(n,m,theta2,phi2);
		    for(int k=0;k<n_dgele_dof2;k++){
		      ///////////////////
		      //  rhs(dgele_dof2[k]) += -1/R*Theta*u_hat*Y2*Jxw2*bas_val2[k][l2];
		      //////////////////
		      cvaltype cout = Jxw2*Y2*bas_val2[k][l2];
		      cvaltype value = 1/R*Theta*cout*u_hat;
		      triplets.push_back(T(dgele_dof[j],dgele_dof2[k],value));
		      //triplets.push_back(T(dg_dof[i],dgele_dof2[k],value));
		    } 	      
		  }
		}
	    }
	  }
	}
    }
  }
  end = clock();
  std::cout<< " time: " << (double)(end-start)/CLOCKS_PER_SEC<<std::endl;
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
      std::vector<cvaltype> u_h_grad(3);
      u_h_grad[0] = u_re_grad[l][0]+I*u_im_grad[l][0];
      u_h_grad[1] = u_re_grad[l][1]+I*u_im_grad[l][1];
      u_h_grad[2] = u_re_grad[l][2]+I*u_im_grad[l][2];
      cvaltype u_exact_val = u_exact(q_pnt[l]);
      cvec_type u_exact_grad = u_exact_prime(q_pnt[l]);
      
      double df_value = std::norm(u_exact_val-u_h_val);
      double df_grad = std::norm(u_exact_grad[0]-u_h_grad[0])+std::norm(u_exact_grad[1]-u_h_grad[1])+std::norm(u_exact_grad[2]-u_h_grad[2]);
      
      L2error += Jxw*df_value;
      L2error_grad += Jxw*df_grad;

      residual += Jxw*(-c_val*(u_h_grad[0]*u_h_grad[0]+u_h_grad[1]*u_h_grad[1]+u_h_grad[2]*u_h_grad[2])+a_val*u_h_val);
    }
    ec.residual = std::abs(residual);
  }
  L2error = sqrt(fabs(L2error));
  L2error_grad = sqrt(fabs(L2error_grad));
  
  std::cerr << "\nL2 error = " << L2error << std::endl;
  std::cerr << "\nL2 gradient error = " << L2error_grad << std::endl;
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
  mesh_adaptor.tolerence() = 0.0003;
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
      //std::cout<<indicator[ele_idx]<<std::endl;
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
// read mesh file
uiExperiment::uiExperiment(const std::string& file)
{
  mesh_file = file;
  order = Order;
  n_bmark = 2;
  //u_hat.resize((order+1)*(order+1));
};

uiExperiment::~uiExperiment()
{};

void uiExperiment::init()
{
  //change mesh data to local mesh data
  h_tree.readMesh(mesh_file);
  ir_mesh.reinit(h_tree);
    
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

  dg_template_element.resize(1);
  dg_template_element[0].reinit(triangle_template_geometry,triangle_to3d_coord_transform);
  
  
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
    // if vertex is 4, use tetrahedron, 5 use twin tetrahedron and 7 use four tetrahedron
    int n_vtx = mesh.geometry(DIM,i).n_vertex();
    if(n_vtx == 4)
      fem_space.element(i).reinit(fem_space,i,0);
    else if(n_vtx == 5)
      fem_space.element(i).reinit(fem_space,i,1);
    else if(n_vtx == 7)
      fem_space.element(i).reinit(fem_space,i,2);
    else
      {
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

  Mesh<DIM,DIM>& mesh = fem_space.mesh();
  // calculate Neumann boundary number
  int n_surface = mesh.n_geometry(DIM-1);
  int n_dg_ele = 0;
  for(int i=0;i<n_surface;i++){
    // here only for mark=2 
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark)
      n_dg_ele +=1;
  }

  // build DGElement
  fem_space.dgElement().resize(n_dg_ele);
  for(int i=0,j=0;i<n_surface;i++){
    if(mesh.geometry(DIM-1,i).boundaryMark() == bmark){
      fem_space.dgElement(j).reinit(fem_space,i,0);
      j += 1;
    }
  }
  fem_space.buildDGElement();
  
  std::cout << "Update DG Data Cache ...";
  std::cout << " bmark: "<<bmark<<" ..."<< " surface: "<< n_dg_ele<<" ...";
  if(bmark_count == n_bmark)
    std::cerr<<"Error: not enough bmark number"<<std::endl;
  updateDGGeometryCache(edge_cache[bmark_count]);
  std::cout << "OK!" << std::endl;

}

void uiExperiment::solve()
{
  
  getMat();
  getRhs();
  
  NeummanBC(g,1);

  TransparentBC(2);
  
  // DirichletBC(bnd,2);
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
    file << interp_point[0] << " " << interp_point[1] << " " << interp_point[2] << " " << u_h(i) << std::endl;
  }
  file.close();

}

void uiExperiment::saveData()
{
  u_re.writeOpenDXData("u_re.dx");
  u_im.writeOpenDXData("u_im.dx");
  u_exact_re.writeOpenDXData("u_exact_re.dx");
  Error.writeOpenDXData("error.dx");
  //writeMatlabData("u_r.dat",u_re);
  //writeMatlabData("u_i.dat",u_im);
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
    adaptMesh();
    delete[] edge_cache;
    // edge_cache = NULL;
    getchar();
  }while(1);
};

/**
 *  * end of file
 *
 */

