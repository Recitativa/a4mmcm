CXXFLAGS= -I/usr/include/boost-1_33_1/ -I/usr/local/include  
LDFLAGS	= -L/usr/local/lib -lgsl -lgslcblas 
BOOSTFLAGS = -lboost_unit_test_framework-gcc-mt

default: run

test: theoretic_test.exe
	theoretic_test

run: simquantile.exe
	simquantile

simquantile.exe: simquantile.cpp theoretic.cpp theoretic.hpp
	g++ -Wall -mtune=prescott -march=prescott -o $@ $< theoretic.cpp $(CXXFLAGS) $(LDFLAGS)

theoretic_test.exe: theoretic_test.cpp theoretic.cpp theoretic.hpp 
	g++ -Wall -g -o $@ $< theoretic.cpp $(CXXFLAGS) $(LDFLAGS) $(BOOSTFLAGS)


clean:
	rm -f *.o 
	rm -f *.exe

astyle:
	sh -c "astyle --style=kr --pad=oper --indent=spaces=2  *.cpp *.hpp"