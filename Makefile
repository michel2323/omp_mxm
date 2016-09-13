

CXX=icc
OMP=-qopenmp
#CFLAGS=-O3 -march=native -xMIC-AVX512 -qopt-report=5 -DHBW 
CFLAGS=-g -O3 -xMIC-AVX512 -qopt-report=5 -restrict #-DHBW
CFLAGS+= -I/homes/schanen/git/libxsmm/include
#-DMEMKIND #-DHBW 
LFLAGS=-L/homes/schanen/git/libxsmm/lib -lmemkind -lxsmm 

all: mxm_openmp

mxm_openmp: mxm_openmp.o
	$(CXX) $(OMP) $(CFLAGS) $(LFLAGS) -o $@ $< /homes/schanen/git/libxsmm/lib/libxsmm.a /homes/schanen/git/libxsmm/lib/libxsmmext.a 

%.o: %.cpp
	$(CXX) $(OMP) $(CFLAGS) -c $<

clean:
	-rm mxm_openmp *.o

.PHONY: clean
