#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <cstdio> 
#include <unistd.h>
#include <numa.h>
#include <hbwmalloc.h>
#include <libxsmm.h>
#define DEF_RESTRICT restrict


using namespace std;

#define GEMM(REAL) LIBXSMM_FSYMBOL(LIBXSMM_TPREFIX(REAL, gemm))

int main (int argc, char * argv[])

{
  if(argc < 5) {
    printf("Please enter size as argument\n");
    exit(1);
  }
  long int n=atoi(argv[1]);
  long int m=atoi(argv[2]);
  long int l=atoi(argv[3]);
  long int ram=atoi(argv[4]);
  void libxsmm_init(void);
  libxsmm_dmmfunction xsmm_dgemm;
  xsmm_dgemm=libxsmm_dmmdispatch((int) n, (int) m, (int) l,
        NULL, NULL, NULL, NULL, NULL,NULL, NULL);
  long int nel=(long int) ram/(((double) (m*l*8))/(1024.*1024.));
  double angle;
  long int i,j,k;
  double pi = 3.141592653589793;
  double s;
  int thread_num;

  int num_kind;
#ifdef HBW
  double * DEF_RESTRICT a=(double*) hbw_malloc(sizeof(double)*n*m); 
  double * DEF_RESTRICT b=(double*) hbw_malloc(sizeof(double)*nel*m*l); 
  double * DEF_RESTRICT c=(double*) hbw_malloc(sizeof(double)*nel*n*l); 
#else
  double * DEF_RESTRICT a=(double*) malloc(sizeof(double)*n*m); 
  double * DEF_RESTRICT b=(double*) malloc(sizeof(double)*nel*m*l); 
  double * DEF_RESTRICT c=(double*) malloc(sizeof(double)*nel*n*l); 
  //double * DEF_RESTRICT a=(double*) numa_alloc_onnode(sizeof(double)*n*m,4); 
  //double * DEF_RESTRICT b=(double*) numa_alloc_onnode(sizeof(double)*nel*m*l,4); 
  //double * DEF_RESTRICT c=(double*) numa_alloc_onnode(sizeof(double)*nel*n*l,4); 
#endif
  cout << "\n";
  cout << "  Compute matrix product C = A * B.\n";

  thread_num = omp_get_max_threads ( );

  cout << "\n";
  cout << "The number of threads used=" << thread_num <<  "\n";

  cout << "nel= " << nel << "\n";
  cout << "n= " << n << "\n";
  cout << "m= " << m << "\n";
  cout << "l= " << l << "\n";
  cout << "A uses " << (n*m*8)/(1024*1024) << "MB of RAM\n";
  cout << "B uses " << (nel*m*l*8)/(1024*1024) << "MB of RAM\n";
  cout << "C uses " << (nel*n*l*8)/(1024*1024) << "MB of RAM\n";
  s = 1.0 / sqrt ( ( double ) ( n ) );
  double wtime0 = omp_get_wtime ( );
  double wtime1;

# pragma omp parallel 
  {
# pragma omp for 
    for(i=0;i<n*m;i++) {
        a[i]=2.;
    }
# pragma omp for
    for(int i=0;i<nel*m*l;i++) { 
          b[i]=3.;
        }
  }
    wtime1=omp_get_wtime ( );
# pragma omp parallel
    {
      for(int i=0;i<10;i++) {
# pragma omp for 
        for(int iel=0;iel<nel;iel++) { 
          //xsmm_dgemm(a, &b[iel*m*l], &c[iel*l*n]);
        }
      }
    }
  double wtime2 = omp_get_wtime ( );
  cout << "Elapsed seconds w1=" << wtime1-wtime0 << " w2="<< wtime2-wtime1 << " w3=" << " total=" << wtime2-wtime0 << "\n";
  cout << "GB/s:" << (2*nel*(m*l+l*n))/((wtime2-wtime1)*1024.*1024.*1024.) << endl;
  cout << "GFLOPS/s:" << (nel*(n*l*m))/((wtime2-wtime1)*1024.*1024.*1024.) << endl;
  cout << "\n";
#ifdef HBW
  hbw_free(a);
  hbw_free(b);
  hbw_free(c);
#else
  free(a);
  free(b);
  free(c);
  //numa_free(a,sizeof(double)*n*m);
  //numa_free(b,sizeof(double)*nel*m*l);
  //numa_free(c,sizeof(double)*nel*n*l);
#endif
return 0;
}
