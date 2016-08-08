#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <cstdio> 
#include <unistd.h>
#include <memkind.h>
#include <hbwmalloc.h>
#define DEF_RESTRICT restrict

using namespace std;

int main (int argc, char * argv[])

{
  if(argc < 2) {
    printf("Please enter size as argument\n");
    exit(1);
  }
  int n=atoi(argv[1]);
  //if(hbw_set_policy(HBW_POLICY_INTERLEAVE)) {
    //cout << "Error on hbw policy setting" << endl;
  //}
  //cout << hbw_get_policy() << endl;
  double angle;
  int i;
  int j;
  int k;
  double pi = 3.141592653589793;
  double s;
  int thread_num;

  int num_kind;
  memkind_t kind;
  memkind_get_num_kind(&num_kind); 
  cout << "Kinds of memory: " << num_kind << endl;
  for(i=0;i<num_kind;i++) {
    size_t total, free;
    memkind_get_kind_by_partition(i,&kind);
    memkind_get_size(kind, &total, &free); 
    cout << "Partition=" << kind << ": \t total=" << total << ", \tfree=" << free << endl;
  }
#ifdef HBW
  double * DEF_RESTRICT a=(double*) hbw_malloc(sizeof(double)*n*n); 
  double * DEF_RESTRICT b=(double*) hbw_malloc(sizeof(double)*n*n); 
  double * DEF_RESTRICT c=(double*) hbw_malloc(sizeof(double)*n*n); 
#else
#ifdef MEMKIND
  kind=MEMKIND_HBW;
  //memkind_get_kind_by_partition(part,&kind);
  //int part=atoi(argv[2]);
  cout << "Allocating on partition " << kind << "." << endl;
  double * DEF_RESTRICT a=(double*) memkind_malloc(kind,sizeof(double)*n*n); 
  double * DEF_RESTRICT b=(double*) memkind_malloc(kind,sizeof(double)*n*n); 
  double * DEF_RESTRICT c=(double*) memkind_malloc(kind,sizeof(double)*n*n); 
#else
  double * DEF_RESTRICT a=(double*) malloc(sizeof(double)*n*n); 
  double * DEF_RESTRICT b=(double*) malloc(sizeof(double)*n*n); 
  double * DEF_RESTRICT c=(double*) malloc(sizeof(double)*n*n); 
#endif
#endif
  cout << "After allocation" << endl << "---------------------------\n";
  for(i=0;i<num_kind;i++) {
    memkind_t kind;
    size_t total, free;
    memkind_get_kind_by_partition(i,&kind);
    memkind_get_size(kind, &total, &free); 
    cout << "Partition=" << i << ": \t total=" << total << ", \tfree=" << free << endl;
  }

  cout << "\n";
  cout << "  Compute matrix product C = A * B.\n";

  thread_num = omp_get_max_threads ( );

  cout << "\n";
  cout << "The number of threads used=" << thread_num <<  "\n";

  cout << "n= " << n << "\n";
  s = 1.0 / sqrt ( ( double ) ( n ) );
  double wtime0 = omp_get_wtime ( );
  double wtime1;

# pragma omp parallel 
  {
# pragma omp for nowait 
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        angle=2.0*pi*i*j/(double)n;
        a[i*n+j]=s*(sin(angle)+cos(angle));
      }
    }
# pragma omp for nowait
    for(i=0;i<n*n;i++) c[i] = 0.0;
# pragma omp for
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        b[i*n+j]=a[i*n+j];
      }
    }
    wtime1=omp_get_wtime ( );
# pragma omp for
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
        for(k=0;k<n;k++) {
          c[i*n+j]=c[i*n+j]+a[i*n+k]*b[k*n+j];
        }
      }
    }
  }
  double wtime2 = omp_get_wtime ( );
  cout << "Elapsed seconds w1=" << wtime1-wtime0 << " w2="<< wtime2-wtime1 << " w3=" << " total=" << wtime2-wtime0 << "\n";
  cout << "C(100,100)=" << c[99*n+99] << "\n";
  cout << "\n";
#ifdef HBW
  hbw_free(a);
  hbw_free(b);
  hbw_free(c);
#else
#ifdef MEMKIND
  memkind_free(kind,a);
  memkind_free(kind,b);
  memkind_free(kind,c);
#else
  free(a);
  free(b);
  free(c);
#endif
#endif
return 0;
}
