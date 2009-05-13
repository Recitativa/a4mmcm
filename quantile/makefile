CXXFLAGS= -I/usr/include/boost-1_33_1/ -I/usr/local/include  -L/usr/local/lib
LDFLAGS	= -lgsl -lgslcblas -lboost_unit_test_framework-gcc-mt

test: theoretic_test.exe
	theoretic_test

theoretic_test.exe: theoretic_test.cpp theoretic.cpp theoretic.hpp 
	g++ -o $@ $< $(CXXFLAGS) $(LDFLAGS)

default: test

clean:
	rm -f *.o 
	rm -f *.exe