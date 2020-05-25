CFLAGS =  -Iinclude
BDBWT = 
SDSL = -I include -I $(BDBWT)sdsl-lite/include -I ./$(BDBWT)sdsl-lite/build/external/libdivsufsort/include -L $(BDBWT)sdsl-lite/build/lib -lsdsl -L $(BDBWT)sdsl-lite/build/external/libdivsufsort/lib/ -ldivsufsort -ldivsufsort64
all:
	g++ $(BDBWT)rsa1d.hh $(BDBWT)rsa2d.hh $(BDBWT)mem.cpp -g $(SDSL) -o mem -std=c++14


