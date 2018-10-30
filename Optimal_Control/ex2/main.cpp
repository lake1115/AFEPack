/*****************************************************************************
 *         
 *   @file main.cpp
 *   @author Hu Bin, Zhejiang University <binh@zju.edu.cn>
 *   @data 2018/10/25
 *
 *   @brief solve the optimal control PDE:  
 *        -div(A nabla y) + phi(y) = f + Bu
 *        min(g(y)+j(y))
 *
 *   @input parameter A, B and function phi(y), g(y),j(y)
 *   boundary condition in parameter.h
 *
 *   Dirichlet boundary :  y = bnd;
 *   Neumann boundary :   n*A*grad(y) + py = g  ?  
 *
 *   ./main meshfile *
 *****************************************************************************
 */
#include "uiexp.h"
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

