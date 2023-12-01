#include <string.h>
#include <cuda_runtime.h>

 
extern "C" {
cudaError_t c_cudamemgetinfo( int *meg_free, int *meg_total)
{

  size_t l_free  = 0;
  size_t l_total = 0;
  
  cudaError_t error_memgetinfo;
  error_memgetinfo = cudaMemGetInfo(&l_free, &l_total);
  
  long long ll_free = (long long) l_free;
  long long ll_total = (long long) l_total;
  
  *meg_free  = (int) (ll_free  / 1048576);
  *meg_total = (int) (ll_total / 1048576);
  
  return error_memgetinfo;
}
 
} 
