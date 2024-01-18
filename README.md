# AFEPack for TBC

安装AFEPack 请参考https://github.com/wangheyu/AFEPack/tree/master

# Boost-1.50.0

* 遇到error: invalid conversion from ‘const void*’ to ‘void*’ [-fpermissive]

需要修改文件boost/libs/locale/src/icu/formatter.cpp

61行处： icu_fmt_->format(value,tmp); 改为
 icu_fmt_->format(::int64_t(value),tmp);

* 提示找不到pyconfig.h
  
手动修改 project-config.jam文件
18行处改为
Using python : 2.7 : /usr/include/python2.7 ;
然后执行
./bootstrap.sh --with-python=python2.7

# Deal.II-8.1.0

tar -xzf dealii-8.1.0.tar.gz
cd deal.II
mkdir build && cd build

根据 boost 安装路径修改
export BOOST_DIR=/path/to/boost/install
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 -DDEAL_II_WITH_MPI=on -DDEAL_II_WITH_THREADS=off -DCMAKE_C_COMPILER=/usr/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_Fortran_COMPILER=/usr/bin/mpifort -DDEAL_II_WITH_PETSC=OFF -DHDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi ..

这里注意不要打开anaconda环境，不然会默认调用anaconda下cmake的boost-1.73.0
make
sudo make install

设置一下deal.ii的库链接：
sudo ln -s /usr/local/dealii-8.1.0/lib/* /usr/local/lib

# AFEPack运行

make
easymesh D.d
./main D

