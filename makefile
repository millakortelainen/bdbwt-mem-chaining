CFLAGS =  -Iinclude
BDBWT = 
SDSL = -I include -I $(BDBWT)sdsl-lite/include -I ./$(BDBWT)sdsl-lite/build/external/libdivsufsort/include -L $(BDBWT)sdsl-lite/build/lib -lsdsl -L $(BDBWT)sdsl-lite/build/external/libdivsufsort/lib/ -ldivsufsort -ldivsufsort64
EDLIB = $(BDBWT)edlib/src/edlib.cpp -o
HEADERS = $(BDBWT)io.hh $(BDBWT)util.hh $(BDBWT)rsa1d.hh $(BDBWT)mem.hh $(BDBWT)minimizer.hh $(BDBWT)driver.hh $(BDBWT)chaining.hh
RMAX = $(BDBWT)RMaxQTree.h $(BDBWT)RMaxQTree.cpp

all:
	g++  -fopenmp $(EDLIB) $(RMAX) $(HEADERS) $(BDBWT)main.cpp $(SDSL) -I ./$(BDBWT)edlib/include -o main -std=c++14
debug:
	g++  -fopenmp $(EDLIB) $(RMAX) $(BDBWT)io.hh -g $(BDBWT)util.hh -g $(BDBWT)rsa1d.hh -g -g $(BDBWT)driver.hh -g $(BDBWT)mem.hh -g $(BDBWT)minimizer.hh -g $(BDBWT)chaining.hh -g $(BDBWT)main.cpp -g $(SDSL) -I ./$(BDBWT)edlib/include -o main -std=c++14
test:
	g++  -fopenmp $(EDLIB) $(HEADERS) $(BDBWT)testSuite.cpp $(SDSL) -I ./$(BDBWT)edlib/include -o testSuite -std=c++14
prof:
	g++  -fopenmp $(EDLIB) $(BDBWT)io.hh -g $(BDBWT)util.hh $(BDBWT)rsa1d.hh -g $(BDBWT)mem.hh $(BDBWT)minimizer.hh -g $(BDBWT)main.cpp -g $(SDSL) -I ./$(BDBWT)edlib/include -o main -std=c++14 -pg
mutator:
	g++	 $(BDBWT)textMutator.cpp -o textMutator -std=c++14
