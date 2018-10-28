/**
 * @file   GmshMesh.cpp
 * @author Bin Hu
 * @date   Sat Oct 27 2018
 * 
 * @brief  本程序作用: 读取 Gmsh 的 网格文件 "*.msh", 并将其转化成 AFEPack 适用的网格文件 "*.mesh".
 *         注意: Gmsh 有3种文件格式. 本程序适用于 Version 2.*.
 *         两种格式的区别参见 Gmsh 的 Reference manual. 该文档可以从 http://geuz.org/gmsh 下载.
 * 
 * 
 */

#include "GmshMesh.h"

AFEPACK_OPEN_NAMESPACE

#define N_POINT_NODE		1
#define N_LINE_NODE		2
#define N_TRIANGLE_NODE		3
#define N_QUADRANGLE_NODE	4
#define N_TETRAHEDRON_NODE	4
#define N_HEXAHEDRON_NODE	8
#define N_PRISM_NODE		6
#define N_PYRAMID_NODE		5


GmshMesh::GmshMesh()
{
  setTemplateGeometry(*(new std::vector<TemplateGeometry<3> >(1)));
  templateGeometry(0).readData("tetrahedron.tmp_geo");
  //   templateGeometry(1).readData("hexahedron.tmp_geo");
  //   templateGeometry(2).readData("prism.tmp_geo");
  //   templateGeometry(3).readData("pyramid.tmp_geo");
};

GmshMesh::~GmshMesh()
{
  delete &(templateGeometry());
};

void GmshMesh::readData(const std::string& filename)
{
  int i,j,k;
  double k;
  int n_point, n_line, n_surface, n_element;
  std::string dummy;
  std::ifstream is(filename.c_str());

  std::cout << "Reading gmsh date file ..." << std::endl;
  is >> dummy;/// $MeshFormat
  is >> dummy;/// integer : for example, 2.2 means the version of Gmsh is 2.2
  is >> dummy;/// integer : for example, 0 means the type of data file is ASCII
  is >> dummy;/// integer : the size of the floating point numbers used in the file.
  is >> dummy;/// $EndMeshFormat$
  is >> dummy;/// $Nodes

  std::cout << "\t reading node date ..." << std::endl;
  is >> n_point;
  this->Mesh<3,DOW>::point().resize(n_node);
  this->Mesh<3,DOW>::geometry(0)
}
