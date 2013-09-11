CXX = g++

all : dna2bit

dna2bit :  
	$(CXX) -Wall -g -o dna2bit dna2bit.cpp

clean :
	rm -f dna2bit
