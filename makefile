CFLAGS =  -Iinclude
BDBWT = 
SDSL = -I include -I $(BDBWT)sdsl-lite/include -I ./$(BDBWT)sdsl-lite/build/external/libdivsufsort/include -L $(BDBWT)sdsl-lite/build/lib -lsdsl -L $(BDBWT)sdsl-lite/build/external/libdivsufsort/lib/ -ldivsufsort -ldivsufsort64
EDLIB = $(BDBWT)edlib/src/edlib.cpp -o
all:
	g++  $(EDLIB) $(BDBWT)io.hh -g $(BDBWT)util.hh $(BDBWT)sortUtil.hh $(BDBWT)rsa1d.hh $(BDBWT)rsa2d.hh -g $(BDBWT)mem.hh -g $(BDBWT)main.cpp -g $(SDSL) -I ./$(BDBWT)edlib/include -o main -std=c++14


