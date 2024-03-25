// File    :   FemmProperty.cpp
// Time    :   2024/03/12 15:42:07
// Author  :   Hu Bin 
// Version :   1.0


#include "FemmProperty.h"

/////////////////////////////////////////////
#define THIS LinearSystem
THIS::LinearSystem(int d){
  A = Eigen::MatrixXcd::Zero(d, d);
  b = Eigen::VectorXcd::Zero(d);
}

Eigen::VectorXcd THIS::solve(){
  return A.llt().solve(b);
}
#undef THIS

/////////////////////////////////////////////
#define THIS MaterialProp

void THIS::GetSlopes(double omega){
  if (!slope.empty()) return;
  int n = Bdata.size(); // BH points
  int i,k;
  double l1,l2;
  bool CurveOK = false;
  std::vector<double> bn(n);
  std::vector<Complex> hn(n);
  
  LinearSystem L(n);

	// strip off some info that we can use during the first nonlinear iteration;
	mu = Bdata[1]/(mu0*abs(Hdata[1]));

	// first, we need to doctor the curve if the problem is being evaluated at a nonzero frequency.
  if(omega >0){
    // Make an effective B-H curve for harmonic problems.
		// this one convolves B(H) where H is sinusoidal with a sine at the same frequency to get the effective amplitude of B
    double mu_now, mu_max=0;

		for(i=1;i<n;i++)
		{

			for(k=1,bn[i]=0;k<=i;k++)
			{
				bn[i] += real((4.*(Hdata[k]*Bdata[k-1] -
					Hdata[k-1]*Bdata[k])*(-cos((Hdata[k-1]*PI)/(2.*Hdata[i])) +
					cos((Hdata[k]*PI)/(2.*Hdata[i]))) + (-Bdata[k-1] +
					Bdata[k])*((Hdata[k-1] - Hdata[k])*PI +
					Hdata[i]*(-sin((Hdata[k-1]*PI)/Hdata[i]) +
					sin((Hdata[k]*PI)/Hdata[i]))))/ ((Hdata[k-1] - Hdata[k])*PI));
			}
		}
    for(i=1;i<n;i++){
      Bdata[i] = bn[i];
      mu_now = real(Bdata[i]/Hdata[i]);
      if (mu_now > mu_max){
        mu_max = mu_now;
      }
    }
    for(i=1;i<n;i++)
		{
			Hdata[i] *= exp(I*Bdata[i]*phi_h*DEG/(Hdata[i]*mu_max));
		} 

  }

  while(CurveOK != true){
    // impose natural BC on the `left'
    l1 = Bdata[1] - Bdata[0];
		L.A(0, 0) = 4./l1;
		L.A(0, 1) = 2./l1;
		L.b(0) = 6.*(Hdata[1]-Hdata[0])/(l1*l1);

		// impose natural BC on the `right'
		l1 = Bdata[n-1] - Bdata[n-2];
		L.A(n-1, n-1) = 4./l1;
		L.A(n-1, n-2) = 2./l1;
		L.b(n-1) = 6.*(Hdata[n-1]-Hdata[n-2])/(l1*l1);

		for(i=1;i<n-1;i++)
		{
			l1 = Bdata[i]-Bdata[i-1];
	    l2 = Bdata[i+1]-Bdata[i];

			L.A(i, i-1) = 2./l1;
			L.A(i, i) = 4.*(l1+l2)/(l1*l2);
			L.A(i, i+1) = 2./l2;

			L.b(i) = 6.*(Hdata[i]-Hdata[i-1])/(l1*l1) +
				   6.*(Hdata[i+1]-Hdata[i])/(l2*l2);
		}    
    Eigen::VectorXcd x = L.solve();
    std::vector<Complex> vec_x(&x[0], x.data()+x.cols()*x.rows());
    slope = vec_x;
    // slope = L.solve();

		// now, test to see if there are any "bad" segments in there.
		// it is probably sufficient to do this test just on the real part of the BH curve... 
    for(i=1,CurveOK=true;i<n;i++)
		{
			double L,c0,c1,c2,d0,d1,u0,u1,X0,X1;

			// it is probably sufficient to do this test just on the 
			// real part of the BH curve.  We do the test on just the
			// real part of the curve by taking the real parts of slope and Hdata
			d0 = real(slope[i-1]);
			d1 = real(slope[i]);
			u0 = real(Hdata[i-1]);
			u1 = real(Hdata[i]);
			L = Bdata[i]-Bdata[i-1];

			c0 = d0;
			c1 = -(2.*(2.*d0*L + d1*L + 3.*u0 - 3.*u1))/(L*L);
			c2 = (3.*(d0*L + d1*L + 2.*u0 - 2.*u1))/(L*L*L);
			X0 = -1.;
			X1 = -1.;

			u0 = c1*c1-4.*c0*c2;
			// check out degenerate cases
			if(c2==0)
			{
				if(c1!=0) X0 = -c0/c1;
			}
			else if(u0>0)
			{
				u0=sqrt(u0);
				X0 = -(c1 + u0)/(2.*c2);
				X1 = (-c1 + u0)/(2.*c2);
			}

			//now, see if we've struck gold!
			if (((X0>=0.)&&(X0<=L))||((X1>=0.)&&(X1<=L)))
				CurveOK=false;
		}

    if(CurveOK!=true)  //remedial action
		{
			// Smooth out input points
			// to get rid of rapid transitions;
			// Uses a 3-point moving average
			for(i=1;i<n-1;i++)
			{
				bn[i] = (Bdata[i-1]+Bdata[i]+Bdata[i+1])/3.;
				hn[i] = (Hdata[i-1]+Hdata[i]+Hdata[i+1])/3.;
			}
	
			for(i=1;i<n-1;i++){
				Hdata[i] = hn[i];
				Bdata[i] = bn[i];
			}
		}

  }
}

void THIS::GetBHProps(double B, Complex &v, Complex &dv)
{
	double b,z,z2,l;
	Complex h,dh;
	int i;

	b=fabs(B);
	int n = Bdata.size();
	if(n==0){
		v=mu;
		dv=0;
		return;
	}

	if(b==0){
		v=slope[0];
		dv=0;
		return;
	}

	if(b>Bdata[n-1]){
		h=(Hdata[n-1] + slope[n-1]*(b-Bdata[n-1]));
		dh=slope[n-1];
		v=h/b;
		dv=0.5*(dh/(b*b) - h/(b*b*b));
		return;
	}

	for(i=0;i<n-1;i++)
		if((b>=Bdata[i]) && (b<=Bdata[i+1])){
			l=(Bdata[i+1]-Bdata[i]);
			z=(b-Bdata[i])/l;
			z2=z*z;
			h=(1.-3.*z2+2.*z2*z)*Hdata[i] +
			  z*(1.-2.*z+z2)*l*slope[i] +
			  z2*(3.-2.*z)*Hdata[i+1] +
			  z2*(z-1.)*l*slope[i+1];
			dh=6.*z*(z-1.)*Hdata[i]/l +
			  (1.-4.*z+3.*z*z)*slope[i] +
			  6.*z*(1.-z)*Hdata[i+1]/l +
			  z*(3.*z-2.)*slope[i+1];
			v=h/b;
			dv=0.5*(dh/(b*b) - h/(b*b*b));
			return; 
		}
}

#undef THIS
/////////////////////////////////////////////
#define THIS FemmProp
std::string StripKey(std::string s, bool front=0) {
  std::stringstream ss(s);
  std::string v;
  std::vector<std::string> list;
  while (ss >> v) {
    list.push_back(v);
  }  // read the last word in string
  if (list.empty())
    return "";
  if (front)
    return list.front();
  else
    return list.back();
}

std::map<std::string, double> mapLengthUnits;
std::map<std::string, int> mapBoundayType;

void register_map() {
  // define length units
  mapLengthUnits.insert(std::pair<std::string, double>("inches", 2.54));
  mapLengthUnits.insert(std::pair<std::string, double>("millimeters", 0.1));
  mapLengthUnits.insert(std::pair<std::string, double>("centimeters", 1.));
  mapLengthUnits.insert(std::pair<std::string, double>("meters", 100));
  mapLengthUnits.insert(std::pair<std::string, double>("mils", 0.00254));
  mapLengthUnits.insert(std::pair<std::string, double>("microns", 1.e-04));
  // define boundary type
  mapBoundayType.insert(std::pair<std::string, int>("0", 0));   // Dirichlet 
  mapBoundayType.insert(std::pair<std::string, int>("4", 2));   // Periodic
  mapBoundayType.insert(std::pair<std::string, int>("5", 3));   // Antiperiodic
}



THIS::~FemmProp(){};

THIS::FemmProp(const std::string& filename) {
  reinit(filename);
}

void THIS::reinit(const std::string& filename) {
  std::string s, v;
  int boundary_num, BHpoints_num, material_num, circuit_num, block_num, pbc_num;
  int idx;
  int m_idx, c_idx;
  double temp;
  double freq;
  std::ifstream is;
  /// read .fem
  is.open((filename + ".fem").c_str());

  if (!is.is_open()) {
    std::cerr << "Couldn't read from specified .fem file" << std::endl;
    return;
  } else {
    std::cout << "Reading Femm parameter data file ...";
  }

  register_map();
  while (getline(is, s)) {
    s.erase(0, s.find_first_not_of(' '));
    v = StripKey(s,1);
    if (v == "[Frequency]"){
      freq = std::stod(StripKey(s));
    }
    else if (v == "[Precision]"){
      precision = std::stod(StripKey(s));
    }
    else if (v == "[LengthUnits]"){
      length_unit = mapLengthUnits[StripKey(s)];
    }
    else if (v == "[BdryProps]"){  /// boundary property
      boundary_num = std::stoi(StripKey(s));
      if (boundary_num > 0) {
        bprop.resize(boundary_num);
        idx = -1;
      }      
    }
    else if (v == "<BeginBdry>"){
      idx++;
    }
    else if (v == "<BdryName>"){
      bprop[idx].name = StripKey(s);
    }
    else if (v == "<BdryType>"){
      bprop[idx].type = mapBoundayType[StripKey(s)];
    }
    else if (v == "<A_0>"){
      bprop[idx].value = std::stod(StripKey(s));
    }
    else if (v == "[BlockProps]"){ /// block property
      material_num = std::stoi(StripKey(s));
      if (material_num > 0) {
        mprop.resize(material_num);
        idx = -1;
      }     
    }
    else if (v == "<BeginBlock>"){
      idx++;
    }
    else if (v == "<BlockName>"){
      mprop[idx].name = StripKey(s);
    }
    else if (v == "<Mu_x>"){
      mprop[idx].mu = std::stod(StripKey(s));
    }
    else if (v == "<H_c>"){
      mprop[idx].H_c = std::stod(StripKey(s));
    }
    else if (v == "<J_re>"){
      mprop[idx].J = std::stod(StripKey(s));
    }
    else if (v == "<Sigma>"){
      mprop[idx].sigma = std::stod(StripKey(s));
    }
    else if (v == "<Phi_h>"){
      mprop[idx].phi_h = std::stod(StripKey(s));
    }
    else if (v == "<Phi_h>"){
      mprop[idx].lamfill = std::stod(StripKey(s));
    }
    else if (v == "<BHPoints>"){
      BHpoints_num = std::stoi(StripKey(s));
      if (BHpoints_num > 0) {
        mprop[idx].Bdata.resize(BHpoints_num);
        mprop[idx].Hdata.resize(BHpoints_num);
        for (int i=0; i < BHpoints_num; i++){
          is >> mprop[idx].Bdata[i] >> mprop[idx].Hdata[i];
        }
        mprop[idx].GetSlopes(freq);
        nonlinear = true;
      } 
    }
    else if (v == "[CircuitProps]"){ /// circuit property
      circuit_num = std::stoi(StripKey(s));
      if (circuit_num > 0) {
        cprop.resize(circuit_num);
        idx = -1;
      }
    }
    else if (v == "<BeginCircuit>"){
      idx++;
    }
    else if (v == "<CircuitName>"){
      cprop[idx].name = StripKey(s);
    }
    else if (v == "<TotalAmps_re>"){
      cprop[idx].total_amps = std::stod(StripKey(s));
    }
    else if (v == "<CircuitType>"){
      cprop[idx].type = std::stoi(StripKey(s));
    }
    else if (v == "[NumBlockLabels]"){ /// block
      block_num = std::stoi(StripKey(s));
      block.resize(block_num);
      for (int i = 0; i < block_num; i++) {
        is >> block[i].x >> block[i].y >> m_idx >> temp >> c_idx >>
            block[i].magdir >> temp >> block[i].turns >> temp;
        block[i].m_idx = m_idx - 1;
        c_idx = c_idx - 1;  /// -1 means without circuit
        if (c_idx >= 0) {
          if (cprop[c_idx].type == 1) {
            CircuitProp b2cprop;
            b2cprop.name = cprop[c_idx].name;
            b2cprop.type = cprop[c_idx].type;
            b2cprop.total_amps = cprop[c_idx].total_amps * block[i].turns * 0.01;  /// why here has 0,01
            cprop.push_back(b2cprop);
            block[i].c_idx = cprop.size() - 1;
          }
        }
      }
    }
  }
  is.close();
  std::cout << " OK!" << std::endl; 

  /// read .pbc if exist
  is.open((filename + ".pbc").c_str());
  if (!is.is_open()) {
    std::cerr << "Couldn't read from specified .pbc file" << std::endl;
    return;
  } else {
    std::cout << "Reading Femm PBC data file ...";
  }
  is >> pbc_num;
  pbclist.resize(pbc_num);
  for (int i = 0;i < pbc_num;i ++) {
    is >> idx >> pbclist[i].p1 >> pbclist[i].p2 >> pbclist[i].type;
  }
  is.close();
  std::cout << " OK!" << std::endl; 

}
#undef THIS