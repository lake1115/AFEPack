/*****************************************************************************
 *         
 *   Authors	: Bin Hu , Zhejiang University 
 *
 *   Project	: FEM for EM
 *
 *   $Revision: 1.0 $
 *   $Date: 2018/04/17  $
 *   
 *   @brief solve the PDE  -d/dx(c * du/dx) - d/dy(c * du/dy) + a * u = f 
 *   @input parameter c, a, f and boundary condition in parameter.h
 *
 *   Dirichlet boundary :  u = bnd;
 *   Neumann boundary :   n*c*grad(u) + q*u = g  
 *
 *****************************************************************************
 */
#define EMUNITS
#include "uiexp.h"
using namespace std;
int main(int argc, char *argv[])
{
  
  // require mesh file   
  if (argc != 2){
    std::cout << "Usage: " << argv[0] << " meshfile" << std::endl;
    return 1;
  }
  uiExperiment task(argv[1]);
  
  task.run();
  /* 
  task.init();
  int iterations = 2;
  do{
    task.buildFEMSpace();
    //solve linear system
    task.solve();
    //save for Opendx
    // task.saveData();
    //exact solution u, can get L2 error
    //task.getError();

    //begin adaptmesh
    task.adaptMesh();
    
     getchar();
    }while(1);
 
  */
  return 0;
}

