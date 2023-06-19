
default: all

all: sbva

sbva: sbva.cc
	clang++ -I eigen-3.4.0/ -std=c++11 -O3 -o sbva sbva.cc

clean:
	rm sbva
