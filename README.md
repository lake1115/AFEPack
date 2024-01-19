# [AFEPack for TBC]

==============================

## Table of contents
 1. [安装AFEPack](#安装AFEPack)
 2. [运行AFEPack](#运行AFEPack)

## 安装AFEPack
具体请参考https://github.com/wangheyu/AFEPack/tree/master， 非常详细。目前在Ubuntu 18.04上测试过，下面对一些可能的问题做补充

**Boost-1.50.0**
```
wget https://phoenixnap.dl.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
tar -xjf boost_1_50_0.tar.bz2
cd boost_1_50_0 && ./bootstrap.sh --prefix=/path/to/install
./b2
sudo ./b2 install
```
- 遇到error: invalid conversion from ‘const void*’ to ‘void*’ [-fpermissive]

需要修改文件boost/libs/locale/src/icu/formatter.cpp
```
61行处： icu_fmt_->format(value,tmp); 改为
 icu_fmt_->format(::int64_t(value),tmp);
```
- 提示找不到pyconfig.h
  
手动修改 project-config.jam文件

18行处改为
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
根据 boost 安装路径修改
```
export BOOST_DIR=/path/to/boost/install
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 -DDEAL_II_WITH_MPI=on -DDEAL_II_WITH_THREADS=off -DCMAKE_C_COMPILER=/usr/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_Fortran_COMPILER=/usr/bin/mpifort -DDEAL_II_WITH_PETSC=OFF -DHDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi ..
```
这里注意不要打开anaconda环境，不然会默认调用anaconda下cmake的boost-1.73.0
```
make
sudo make install
```
设置一下deal.ii的库链接：
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
