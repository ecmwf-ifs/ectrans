#ifdef TRANS_SINGLE
typedef float DATA_TYPE;
#else
typedef double DATA_TYPE;
#endif

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
void hipfunction(int ISIGNp, int N, DATA_TYPE *data_in, DATA_TYPE *data_out, long *iplan);
#endif
