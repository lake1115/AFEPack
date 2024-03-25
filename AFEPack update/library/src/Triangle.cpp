// File    :   Triangle.cpp
// Time    :   2024/02/27 16:21:23
// Author  :   Hu Bin 
// Version :   1.0
// Brief   :   read triangle file .node, .edge, .ele to AFEPack mesh

#include "Triangle.h"

AFEPACK_OPEN_NAMESPACE

void Triangle::readData(const std::string& filename)
{
    int i, j;
    int DIM;
    int n_node, n_side, n_element;
    bool attribute_flag, marker_flag;
    int temp;

    std::cout << "Reading triangle data file ..." << std::endl;	
    std::cout << "\treading node data ..." << std::flush;
    std::ifstream is((filename + ".node").c_str());
    is >> n_node >> DIM >> attribute_flag >> marker_flag;
    Assert(DIM == 2, ExcMeshData("triangle only support dim = 2"));
    this->Mesh<2,2>::point().resize(n_node);
    this->Mesh<2,2>::geometry(0).resize(n_node);
    for (i = 0;i < n_node;i ++) {
        is  >> j
            >> this->Mesh<2,2>::point(i)
            >> this->Mesh<2,2>::boundaryMark(0,i);
        this->Mesh<2,2>::geometry(0,i).index() = j;
        this->Mesh<2,2>::geometry(0,i).vertex().resize(1,j);
        this->Mesh<2,2>::geometry(0,i).boundary().resize(1,j);
    }
    is.close();
    std::cout << " OK!" << std::endl;
        
    std::cout << "\treading side data ..." << std::flush;
    is.open((filename + ".edge").c_str());
    is >> n_side >> marker_flag;
    this->Mesh<2,2>::geometry(1).resize(n_side);
    for (i = 0;i < n_side;i ++) {
        Geometry& g = this->Mesh<2,2>::geometry(1,i);
        g.vertex().resize(2);
        is  >> g.index()
            >> g.vertex(0) >> g.vertex(1)
            >> this->Mesh<2,2>::boundaryMark(1,i);
        this->Mesh<2,2>::geometry(1,i).boundary() = this->Mesh<2,2>::geometry(1,i).vertex();
    }
    is.close();
    std::cout << " OK!" << std::endl;
	
    std::cout << "\treading element data ..." << std::flush;
    is.open((filename + ".ele").c_str());
    is >> n_element >> temp >> attribute_flag;
    this->Mesh<2,2>::geometry(2).resize(n_element);
    for (i = 0;i < n_element;i ++) {
        Geometry& g = this->Mesh<2,2>::geometry(2,i);
        g.vertex().resize(3);
        g.boundary().resize(3);

        is  >> g.index()
            >> g.vertex(0) >> g.vertex(1) >> g.vertex(2)
            >> this->Mesh<2,2>::boundaryMark(2,i);
        // find element sides by side data  
        temp = 0;
        for (j = 0; j < n_side; j++) {
            if (temp == 3)
                break;
            Geometry& g_side = this->Mesh<2,2>::geometry(1,j);
            if ((g_side.vertex(0) == g.vertex(0) && g_side.vertex(1) == g.vertex(1)) 
            || (g_side.vertex(1) == g.vertex(0) && g_side.vertex(0) == g.vertex(1))){
                g.boundary(2) = j;
                temp ++;
                continue;
            }             
            if ((g_side.vertex(0) == g.vertex(0) && g_side.vertex(1) == g.vertex(2))
            || (g_side.vertex(1) == g.vertex(0) && g_side.vertex(0) == g.vertex(2))){
                g.boundary(1) = j;
                temp ++;
                continue;
            }
            if ((g_side.vertex(0) == g.vertex(1) && g_side.vertex(1) == g.vertex(2))
            || (g_side.vertex(1) == g.vertex(1) && g_side.vertex(0) == g.vertex(2))){
                g.boundary(0) = j;
                temp ++;
                continue;
            }
        } 
     
        const Point<2>& p0 = this->Mesh<2,2>::point(g.vertex(0));
        const Point<2>& p1 = this->Mesh<2,2>::point(g.vertex(1));
        const Point<2>& p2 = this->Mesh<2,2>::point(g.vertex(2));
        if ((p1[0] - p0[0])*(p2[1] - p0[1]) - (p1[1] - p0[1])*(p2[0] - p0[0]) < 0) {
            std::cerr << "vertices of element " << i 
                << " is not correctly ordered." << std::endl;
            j = g.vertex(2);
            g.vertex(2) = g.vertex(1);
            g.vertex(1) = j;
            j = g.boundary(2);
            g.boundary(2) = g.boundary(1);
            g.boundary(1) = j;
        }
   
    }
    is.close();
    std::cout << " OK!" << std::endl;
}

void Triangle::writeOpenDXData(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(12);
  os.setf(std::ios::scientific, std::ios::floatfield);
  u_int n_node = n_point();
	
  os << "object 1 class array type float rank 1 shape 2 item " 
     << n_node << " data follows\n";
  for (u_int i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }
	
  u_int n_element = n_geometry(2);
  os << "\nobject 2 class array type int rank 1 shape 3 item "
     << n_element << " data follows\n";
  for (u_int i = 0;i < n_element;i ++) {
    os << geometry(2,i).vertex(0) << "\t"
       << geometry(2,i).vertex(1) << "\t"
       << geometry(2,i).vertex(2) << "\t\n";
  }
  os << "attribute \"element type\" string \"triangles\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction-2d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "end\n";
  os.close();
}


AFEPACK_CLOSE_NAMESPACE
