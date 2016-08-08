

CXX=icc
OMP=-qopenmp
#CFLAGS=-O3 -march=native -xMIC-AVX512 -qopt-report=5 -DHBW 
CFLAGS=-O3 -march=native -qopt-report=5 -restrict -DMEMKIND #-DHBW 
LFLAGS=-lmemkind

all: mxm_openmp

mxm_openmp: mxm_openmp.o
	$(CXX) $(OMP) $(CFLAGS) $(LFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) $(OMP) $(CFLAGS) -c $<

clean:
	-rm mxm_openmp *.o

.PHONY: clean
