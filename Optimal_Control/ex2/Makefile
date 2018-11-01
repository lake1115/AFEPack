CXX = g++
CXXFLAGS = -O2 -g -std=c++0x -pthread -ftemplate-depth-128 -I/usr/local/include/AFEPack -I/home/binh/afepack

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS = -L/usr/local/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7
LIBS = -lAFEPack  -llac -lbase -ltbb  -pthread -llapack -lblas -lgfortran -ldl -lm  


all : main

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS) 

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx

.PHONY : default clean
