CC = g++-4
run: biit biam simqq bfsg

biit: bitreeit.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

biam: biam.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

simqq: simquantile.cpp theoretic.cpp theoretic.hpp
	$(CC)  -o $@ $< theoretic.cpp -lgsl

bfsg: bfsg.cpp
	$(CC) -o $@ $< 
clean:
	rm -f *.o 
