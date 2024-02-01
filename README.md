# [AFEPack for TBC]

==============================

## Table of contents
 1. [安装AFEPack](#安装AFEPack)
 2. [运行AFEPack](#运行AFEPack)
 3. [调试AFEPack](#调试AFEPack)

## 安装AFEPack
具体请参考https://github.com/wangheyu/AFEPack/tree/master  非常详细。目前在Ubuntu 18.04上测试过，下面对一些可能的问题做补充

首先安装如下依赖
```
sudo apt-get install cmake g++ gcc automake autoconf mpich
sudo apt-get install liblapack-dev libboost-all-dev libtbb-dev libbz2-dev
sudo apt-get install libmumps-dev trilinos-all-dev libsuitesparse-dev libarpack2-dev
```

做如下链接修改， 这里不确定这样做是否最好， 但至少可以继续下去。 （或者定义到你自己指定的一个目录下）,注意先check一下源文件是否存在，可能名称不太一样 
```
sudo ln -s /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so /usr/local/lib/libscalapack.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libptscotch-6.1.3.so /usr/lib/x86_64-linux-gnu/libptscotch.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libptscotcherr-6.1.3.so /usr/lib/x86_64-linux-gnu/libptscotcherr.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libscotch-6.1.3.so /usr/lib/x86_64-linux-gnu/libscotch.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libscotcherr-6.1.3.so /usr/lib/x86_64-linux-gnu/libscotcherr.so
```

**Boost-1.50.0**
```
wget https://phoenixnap.dl.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
tar -xjf boost_1_50_0.tar.bz2
cd boost_1_50_0 && ./bootstrap.sh --prefix=/path/to/install
./b2
sudo ./b2 install
```
- 遇到error: invalid conversion from ‘const void*’ to ‘void*’ [-fpermissive]

在boost/libs/locale/src/icu/formatter.cpp, 61行处:
```
icu_fmt_->format(value,tmp); -> icu_fmt_->format(::int64_t(value),tmp);
```
- 提示找不到pyconfig.h
  
在project-config.jam, 18行处改为:
```
Using python : 2.7 : /usr/include/python2.7 ;
```
然后执行
```
./bootstrap.sh --with-python=python2.7
./b2
sudo ./b2 install
```
**Deal.II-8.1.0**
```
tar -xzf dealii-8.1.0.tar.gz
cd deal.II
mkdir build && cd build
```
在deal.II/include/deal.II/lac/sparsity_pattern.h, 在最开始增加一行：
```
#include <algorithm>
```

1.根据 boost 安装路径修改
```
export BOOST_DIR=/path/to/boost/install //默认是在/usr/local/include下
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 -DDEAL_II_WITH_MPI=on -DDEAL_II_WITH_THREADS=off -DCMAKE_C_COMPILER=/usr/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_Fortran_COMPILER=/usr/bin/mpifort -DDEAL_II_WITH_PETSC=OFF -DHDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi ..
```
- 提示调用系统或者anaconda下cmake的高版本boost，导致找不到安装的boost-1.50

在deal.II/cmake/configure/configure_boost.cmake, line 63, 插入：
```
set(BOOST_ROOT "/usr/local/include")
set(Boost_NO_SYSTEM_PATHS ON) 
FIND_PACKAGE(Boost 1.50 COMPONENTS ${_boost_components})
```
为了验证可以在最后插入:
```
MESSAGE( STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
MESSAGE( STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
MESSAGE( STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")
```

2. make
```
make
sudo make install
```

在Ubuntu 22.04上由于GCC版本 和 OpenMPI版本不同，会出现很多错误：

 - 出现GCC 11 build errors: 'numeric_limits' is not a member of 'std'

在deal.II/include/deal.II/lac/vector.templates.h, 加入一行：
```
#include <limits>
```

 - error: ‘MPI_Type_struct’ was not declared in this scope

在deal.II/source/base/mpi.cc, line 253, 如下修改:
```
ierr = MPI_Type_struct(2, lengths, displacements, types, &type); -> ierr = MPI_Type_create_struct(2, lengths, displacements, types, &type);
```

 - error: ‘const ptree’ {aka ‘const class boost::property_tree::basic_ptree<std::__cxx11::basic_string<char>, std::__cxx11::basic_string<char> >’} has no member named ‘get_optionalstd’;
在deal.II/source/base/parameter_handler.cc, line 1278, 如下修改:
```
return (p.get_optionalstd::string("value")); -> return bool(p.get_optional<std::string>("value"));
```

3.设置一下deal.ii的库链接：
```
sudo ln -s /usr/local/dealii-8.1.0/lib/* /usr/local/lib
```
**AFEPack**
```
cd /usr/local/AFEPack
aclocal
autoconf
automake
CC="/usr/bin/mpicc" CXX="/usr/bin/mpic++" CFLAGS="-I/usr/include/trilinos" CPPFLAGS="-I/usr/include/trilinos" CXXFLAGS="-I/usr/include/trilinos" ./configure --with-dealii="/usr/local/dealii-8.1.0"
make
sudo make install
```

- deal.II依然出现'numeric_limits' is not a member of 'std'问题:
在/usr/local/dealii-8.1.0/include/deal.II/lac/solver_gmres.h, 加入一行：
```
#include <limits>
```

需要在.bashrc中增加AFEPack path
```
export AFEPACK_PATH="/usr/local/AFEPack"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
export AFEPACK_TEMPLATE_PATH="$AFEPACK_PATH/template/triangle:$AFEPACK_PATH/template/twin_triangle:$AFEPACK_PATH/template/interval:$AFEPACK_PATH/template/tetrahedron:$AFEPACK_PATH/template/twin_tetrahedron:$AFEPACK_PATH/template/four_tetrahedron"
```

## 运行AFEPack
运行按如下步骤实现
```
make
easymesh D.d
./main D
```
结果可以通过opendx来看， Ubuntu 通过sudo apt install dx dxsample 安装

终端运行dx，打开界面，然后在Import Data...中选择运行结果u.dx文件，得到可视化结果


## 调试AFEPack
主要通过vscode来调试代码。需要在vscode扩展中安装C/C++ 和C/C++ Extension Pack

在.vscode中放入vscode文件夹下的launch.json, task.json, c_cpp_properties.json 即可
