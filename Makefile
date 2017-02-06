
GITDIR=/home/michel/git
CXX=g++
OMP=-fopenmp
#CFLAGS=-O3 -march=native -xMIC-AVX512 -qopt-report=5 -DHBW 
CFLAGS=-g -O3 
CFLAGS+= -I$(GITDIR)/libxsmm/include
#-DMEMKIND #-DHBW 
LFLAGS=-L$(GITDIR)/libxsmm/lib -lxsmm -ldl -lxsmmnoblas 

all: mxm_openmp

mxm_openmp: mxm_openmp.o
	$(CXX) $(OMP) $(CFLAGS) -o $@ $< $(LFLAGS) 

%.o: %.cpp
	$(CXX) $(OMP) $(CFLAGS) -c $<

clean:
	-rm mxm_openmp *.o

.PHONY: clean
