CXXFLAGS= -I/usr/local/include  
LDFLAGS	= -L/usr/local/lib -lgsl 


simquantile.exe: simquantile.cpp theoretic.cpp theoretic.hpp
	g++ -Wall -o $@ $< theoretic.cpp $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f *.o 
	rm -f *.exe

astyle:
	sh -c "astyle --style=kr --pad=oper --indent=spaces=2  *.cpp *.hpp"

