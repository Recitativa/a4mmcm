CC = g++-4
run: bitree_rec bitreeit biam simqq bfsg

bitree_rec: bitree_rec.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

bitreeit: bitreeit.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

biam: biam.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

simqq: simquantile.cpp theoretic.cpp theoretic.hpp
	$(CC)  -o $@ $< theoretic.cpp -lgsl

bfsg: bfsg.cpp
	$(CC) -o $@ $< 
clean:
	rm -f *.o 

.PHONY: clean bitree_rec test run