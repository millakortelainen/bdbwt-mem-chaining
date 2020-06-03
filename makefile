CFLAGS =  -Iinclude
BDBWT = 
SDSL = -I include -I $(BDBWT)sdsl-lite/include -I ./$(BDBWT)sdsl-lite/build/external/libdivsufsort/include -L $(BDBWT)sdsl-lite/build/lib -lsdsl -L $(BDBWT)sdsl-lite/build/external/libdivsufsort/lib/ -ldivsufsort -ldivsufsort64
all:
	g++  $(BDBWT)util.hh $(BDBWT)sortUtil.hh $(BDBWT)rsa1d.hh $(BDBWT)rsa2d.hh -g $(BDBWT)mem.hh -g $(BDBWT)main.cpp -g $(SDSL) -o main -std=c++14


