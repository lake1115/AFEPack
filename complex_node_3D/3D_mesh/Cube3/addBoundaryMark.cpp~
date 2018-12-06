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

const double radius_cube = 1.;

int onBoundary(const double * p)
{

  if (p[1] == 1.0){
    return 2;
  }
  else if (p[0] == 0.0){
    return 1;
  }
  else if (p[0] == 1.0){
    return 1;
  }
  else if (p[1] == 0.0){
    return 1;
  }
  else if (p[2] == 0.0){
    return 1;
  }
  else if (p[2] == 1.0){
    return 1;
  }
  else{
    return 0;
  }

}

int main(int argc, char* argv[])
{
  Mesh<DIM> mesh;
  mesh.readData(argv[1]);
  

  ////////////////////////
  //count
  int point1 = 0;
  int point2 = 0;
  /////////////////////////

  const int& n_pnt = mesh.n_point();

  for(int i = 0;i < n_pnt;++ i){
    GeometryBM& pnt_geo = mesh.geometry(0, i);
    Point<DIM>& pnt = mesh.point(i);

    if(onBoundary(pnt) == 1){
      pnt_geo.boundaryMark() = 1;
	point1++;
    }
    else if(onBoundary(pnt) == 2){
      pnt_geo.boundaryMark() = 2;
	point2++;
    }
    else{
      pnt_geo.boundaryMark() = 0;
    }

  }

  std::cout << "boundary points are done..." << std::endl;
  std::cout << "point1: " << point1 <<" point2: " << point2  << std::endl;
  const int& n_line = mesh.n_geometry(1);

  for(int i = 0;i < n_line;++ i){
    GeometryBM& line_geo = mesh.geometry(1, i);
    line_geo.boundaryMark() = 0;/// assume it is boundary
    std::vector<int>& vtx = line_geo.vertex();
    
    Point<DIM> p0 = mesh.point(vtx[0]);
    Point<DIM> p1 = mesh.point(vtx[1]);
    
    if(p0[1] == 1.0  && p1[1] == 1.0 ){
      line_geo.boundaryMark() = 2;
    }
    else if(p0[0] == 0.0 && p1[0] == 0.0 ){
      line_geo.boundaryMark() = 1;
    }
    else if(p0[0] == 1.0 && p1[0] == 1.0 ){
      line_geo.boundaryMark() = 1;
    }
    else if(p0[1] == 0.0 && p1[1] == 0.0 ){
      line_geo.boundaryMark() = 1;
    }
    else if(p0[2] == 0.0 && p1[2] == 0.0 ){
      line_geo.boundaryMark() = 1;
    }
    else if(p0[2] == 1.0 && p1[2] == 1.0 ){
      line_geo.boundaryMark() = 1;
    }
  }

  std::cout << "boundary lines are done..." << std::endl;
  
  const int& n_tri = mesh.n_geometry(2);

  ////////////////////////
  //count
  int boundary1 = 0;
  int boundary2 = 0;
  /////////////////////////
  
  for(int i = 0;i < n_tri;++ i){
    GeometryBM& tri_geo = mesh.geometry(2, i);

    tri_geo.boundaryMark() = 0;
    std::vector<int>& vtx = tri_geo.vertex();

    Point<DIM> p0 = mesh.point(vtx[0]);
    Point<DIM> p1 = mesh.point(vtx[1]);
    Point<DIM> p2 = mesh.point(vtx[2]);
    
    if(p0[1] == 1.0  && p1[1] == 1.0 && p2[1] == 1.0 ){
      tri_geo.boundaryMark() = 2;
      boundary2++;
    }
    else if( p0[0] == 0.0  && p1[0] == 0.0 && p2[0] == 0.0 ){
      tri_geo.boundaryMark() = 1;
      boundary1++;
    }
    else if(p0[0] == 1.0  && p1[0] == 1.0 && p2[0] == 1.0 ){
      tri_geo.boundaryMark() = 1;
      boundary1++;
    }
    else if(p0[1] == 0.0  && p1[1] == 0.0 && p2[1] == 0.0  ){
      tri_geo.boundaryMark() = 1;
      boundary1++;
    }
    else if(p0[2] == 0.0 && p1[2] == 0.0 && p2[2] == 0.0 ){
      tri_geo.boundaryMark() = 1;
      boundary1++;
    }
    else if(p0[2] == 1.0 && p1[2] == 1.0 && p2[2] == 1.0){
      tri_geo.boundaryMark() = 1;
      boundary1++;
    }
    
  }

  std::cout << "boundary triangle are done..." << std::endl;
  std::cout << "boundary1: " << boundary1 <<" boundary2: " << boundary2  << std::endl;
  
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
