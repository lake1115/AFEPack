# by R.Lie, Nov 01, 2002

include /usr/local/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS +=  -lAFEPack

all : addBoundaryMark

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@  

addBoundaryMark : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS)
#	$(CXX) -o $@ $(object) $< $(CXXFLAGS)

.PHONY : default clean

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx


