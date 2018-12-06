/*****************************************************************************
 *         
 *   @file main.cpp
 *   @author Hu Bin <binh@zju.edu.cn>
 *   @data 2018/10/23
 *
 *   @brief solve the PDE  -d/dx(c * du/dx) - d/dy(c * du/dy) + a * u = f 
 *   @input parameter c, a, f and boundary condition in parameter.h
 *
 *   Dirichlet boundary :  u = bnd;
 *   Neumann boundary :   n*c*grad(u) + q*u = g  
 *
 *   ./main meshfile
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
  return 0;
}

