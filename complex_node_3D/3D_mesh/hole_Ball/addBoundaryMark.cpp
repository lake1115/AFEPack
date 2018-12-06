/**
 * @file   addBoundaryMark.cpp
 * @author Hu Guanghui <gary@Brin>
 * @date   Sat Mar 12 18:26:37 2016
 * 
 * @brief  add boundary mark on the mesh
 *
 * ./addBoundaryMark input_mesh output_mesh
 * 
 */
#include <iostream>
#include <AFEPack/Geometry.h>

#define DIM 3

const double r = 1.0;
const double R = 2.0;

int onBoundary(const double * p,u_int i)
{
  
  if(fabs(sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - r) < 1.0e-03){
    std::cout << "Point: "<<i<< ", The length is " << sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) << ", and the boundary mark is " << 1 << std::endl;   
    return 1;
  }
  else if(fabs(sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - R) < 1.0e-03){
    std::cout << "Point: "<<i<< ", The length is " << sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) << ", and the boundary mark is " << 2 << std::endl;
    return 2;
  } 
  else{
    std::cout << "Point: "<<i<<", The length is " << sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) << ", and the boundary mark is " << 0 << std::endl;
    return 0;
  }
}

int main(int argc, char* argv[])
{
  Mesh<DIM> mesh;
  mesh.readData(argv[1]);
  
  const int& n_pnt = mesh.n_point();

  for(int i = 0;i < n_pnt;++ i){
    GeometryBM& pnt_geo = mesh.geometry(0, i);
    Point<DIM>& pnt = mesh.point(i);

    if(onBoundary(pnt,i) == 1){
      pnt_geo.boundaryMark() = 1;
    }
    else if(onBoundary(pnt,i) == 2){
      pnt_geo.boundaryMark() = 2;
    }
    else{
      pnt_geo.boundaryMark() = 0;
    }

  }

  std::cout << "boundary points are done..." << std::endl;
  
  const int& n_line = mesh.n_geometry(1);

  for(int i = 0;i < n_line;++ i){
    GeometryBM& line_geo = mesh.geometry(1, i);
    line_geo.boundaryMark() = 0;/// assume it is not boundary
    std::vector<int>& vtx = line_geo.vertex();
    
    if(mesh.geometry(0, vtx[0]).boundaryMark() == 1 && mesh.geometry(0,vtx[1]).boundaryMark() == 1){
      line_geo.boundaryMark() = 1;
    }
    else if(mesh.geometry(0, vtx[0]).boundaryMark() == 2 && mesh.geometry(0,vtx[1]).boundaryMark() == 2){
      line_geo.boundaryMark() = 2;
    }
    
  }

  std::cout << "boundary lines are done..." << std::endl;
  
  const int& n_tri = mesh.n_geometry(2);

  for(int i = 0;i < n_tri;++ i){
    GeometryBM& tri_geo = mesh.geometry(2, i);
    tri_geo.boundaryMark() = 0;
    std::vector<int>& vtx = tri_geo.vertex();

    if(mesh.geometry(0,vtx[0]).boundaryMark() == 1 && mesh.geometry(0,vtx[1]).boundaryMark() == 1 && mesh.geometry(0,vtx[2]).boundaryMark() == 1){
      tri_geo.boundaryMark() = 1;
    }
    else if(mesh.geometry(0,vtx[0]).boundaryMark() == 2 && mesh.geometry(0,vtx[1]).boundaryMark() == 2 && mesh.geometry(0,vtx[2]).boundaryMark() == 2){
      tri_geo.boundaryMark() = 2;
    }
  }

  std::cout << "boundary triangle are done..." << std::endl;
  
  const int& n_tet = mesh.n_geometry(3);

  for(int i = 0;i < n_tet;++ i){
    GeometryBM& tet_geo = mesh.geometry(3, i);
    
    tet_geo.boundaryMark() = 0;

/*
    const int& n_bnd = tet_geo.n_boundary();
    for(int j = 0;j < n_bnd;++ j){
      if(mesh.geometry(2, tet_geo.boundary()[j]).boundaryMark() == 1){
	tet_geo.boundaryMark() = 1;
	break;
      }
      else if(mesh.geometry(2, tet_geo.boundary()[j]).boundaryMark() == 2){
	tet_geo.boundaryMark() = 2;
	break;
      }
    }
*/
  }

  std::cout << "boundary tetrahedron are done..." << std::endl;  
  
  mesh.writeData(argv[2]);
  
  return 0;
}
/**
 * end of file
 * 
 */
