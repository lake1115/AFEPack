///////////////////////////////////////////////////////////////////////////
// DGFEMSpace.cpp : by Hu Bin
//

#include "DGFEMSpace.templates.h"

AFEPACK_OPEN_NAMESPACE

#define DOW 3
#define DIM 3
#define TDIM 3
#define TDIM1 2

template class TemplateDGElement<TDIM1,DOW>;

#define vector_length 1

#define value_type nVector<vector_length,double>
template class DGElement<value_type,DIM,DOW,TDIM,TDIM1>;
template class DGFEMSpace<value_type,DIM,DOW,TDIM,TDIM1>;

template std::vector<double> unitOutNormal(const Point<DIM>&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DIM> >&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template struct DGElementAdditionalData<value_type,DIM>;

#undef value_type
#undef vector_length


#define vector_length 2

#define value_type nVector<vector_length,double>
template class DGElement<value_type,DIM,DOW,TDIM,TDIM1>;
template class DGFEMSpace<value_type,DIM,DOW,TDIM,TDIM1>;

template std::vector<double> unitOutNormal(const Point<DIM>&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DIM> >&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template struct DGElementAdditionalData<value_type,DIM>;

#undef value_type
#undef vector_length


#define vector_length 3

#define value_type nVector<vector_length,double>
template class DGElement<value_type,DIM,DOW,TDIM,TDIM1>;
template class DGFEMSpace<value_type,DIM,DOW,TDIM,TDIM1>;

template std::vector<double> unitOutNormal(const Point<DIM>&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DIM> >&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template struct DGElementAdditionalData<value_type,DIM>;

#undef value_type
#undef vector_length


#undef TDIM1
#undef TDIM
#undef DIW
#undef DOW


AFEPACK_CLOSE_NAMESPACE

//
// end of file
///////////////////////////////////////////////////////////////////////////
