#include <stdio.h>
#include "execute_plan_ffth.hip.h"

extern "C" {

void execute_plan_ffth_c_(int ISIGNp, int N, DATA_TYPE *data_in, DATA_TYPE *data_out, long *iplan)
{
    /*printf("CPU: N=%d\n",N);
    for (int i = 0; i < N; i++)
    {
        printf("CPU: input[%d]=%2.4f\n",i+1,data_in[2*i]);
    }
#pragma omp target data map(to: data_in[0:2*N])  map(from: data_out[0:2*N])
{	
#pragma omp target data use_device_ptr(data_in,data_out)
{*/
    hipfunction(ISIGNp,N,data_in,data_out,iplan);
/*}}
    for (int i = 0; i < N; i++)
    {
        printf("CPU: output[%d]=(%2.4f,%2.4f)\n",i+1,data_out[2*i],data_out[2*i+1]);
    }*/
}}