
run: bitree_rec
	./bitree_rec

bitree_rec: bitree_rec.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

bitreeit: bitreeit.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

biam: biam.cpp
	$(CC) -fopenmp -lgomp -o  $@ $<

simqq: simquantile.cpp theoretic.cpp theoretic.hpp
	g++ -o $@ $< theoretic.cpp 

clean:
	rm -f *.o 

.PHONY: clean bitree_rec test run