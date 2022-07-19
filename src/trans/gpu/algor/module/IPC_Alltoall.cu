// AlltoAll for NVLINK-connected GPUs ithin a single server, using CUDA IPC
// All pairs of (distinct) GPUs must return 1 for canAccessPeer
// Alan Gray, NVIDIA

#include <stdio.h>
#include <assert.h>

#include <stdlib.h>
#include "mpi.h"
#include "cuda_profiler_api.h"
#include <cuda_runtime_api.h>

#include <unistd.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <linux/version.h>


// maximum number of devices supported
#define MAX_DEVICES          32

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}

// data structure required for IPC setup 

typedef struct ipcDevices_st
{
  int count;
  int ordinals[MAX_DEVICES];
} ipcDevices_t;



// we need a seperate CUDA stream for each target GPU
static cudaStream_t streams[MAX_DEVICES];

// structure to contain pointers to remote array data, and offsets into each for destinatin data.
// we maintain 2 copies (first array dimension), corresponding to MTOL and LTOM trans comms.
// This allows us to only perform setup steps the first time, and re-use all following times. 
static double* outputptrall[2][MAX_DEVICES];
static int roff_remote[2][MAX_DEVICES];


// Initialize IP communicatins
extern "C" void initIPC(double* output_d,int* roff, int mtol_or_ltom){


  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  for (int i = 0; i < MAX_DEVICES; i++) {
    cudaStreamCreate(&streams[i]);
  }


  ipcDevices_t *s_devices = (ipcDevices_t *) mmap(NULL, sizeof(*s_devices),
						  PROT_READ | PROT_WRITE, 
						  MAP_SHARED | MAP_ANONYMOUS, 0, 0);
  assert(MAP_FAILED != s_devices);


  // get memory handle for data array local to this MPI process
  cudaIpcMemHandle_t outputhandle;

  cudaIpcGetMemHandle((cudaIpcMemHandle_t *) &outputhandle, (void *) output_d);
  cudaCheckError();


  // create a "full" structure to contain handles for all processes.
  static cudaIpcMemHandle_t outputhandleall[MAX_DEVICES];


  //copy handle for this processor to full structure
  outputhandleall[rank]=outputhandle;

  int irank;
 
  //each rank broadcasts its handle to all other ranks
  for (irank=0;irank<size;irank++){
    MPI_Bcast( &outputhandleall[irank], sizeof(cudaIpcMemHandle_t), 
	       MPI_BYTE, irank, 
    	       MPI_COMM_WORLD );
    }


  //now, using memory handles, populate structure that contains actual pointers to data 
  for (irank=0;irank<size;irank++){
    
    
    if (irank==rank){
      //pointer is locel
      outputptrall[mtol_or_ltom][irank]=output_d;
    }
    else{
      //pointer is remote - use corresponding memory handle
      cudaIpcOpenMemHandle((void **) &outputptrall[mtol_or_ltom][irank], 
			   outputhandleall[irank],
			   cudaIpcMemLazyEnablePeerAccess);
      cudaCheckError();
      
    }
    
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  cudaCheckError();
  

  int targetrank;
  MPI_Status status;
  MPI_Request sendrequest[MAX_DEVICES],recvrequest[MAX_DEVICES];

  
  // each MPI task needs to get remote offset info from loacl info on remote processes
  for(targetrank=0;targetrank<size;targetrank++){


    MPI_Irecv(&(roff_remote[mtol_or_ltom][targetrank]),1,MPI_INTEGER,targetrank,
	      MPI_ANY_TAG,MPI_COMM_WORLD,&recvrequest[targetrank]);
    MPI_Isend(&(roff[targetrank]),1,MPI_INTEGER,targetrank
	      ,0,MPI_COMM_WORLD,&sendrequest[targetrank]);
    
  }

for(targetrank=0;targetrank<size;targetrank++){
    MPI_Wait(&sendrequest[targetrank],&status);
    MPI_Wait(&recvrequest[targetrank],&status);
 }




}


static bool already_initialized[2]={0,0};

static bool notFullPeerAccess=0;

// main externally visible routine for performing AlltoAll comms
extern "C" int Alltoallv_CUDAIPC(double* input, int* len, int* soff, 
				  double* output, int* roff,int mtol_or_ltom){


  int rank,size;
  int targetrank;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  if(notFullPeerAccess)
    return 1;
  
  if (!already_initialized[mtol_or_ltom]){

    //make sure we have full peer access, if not return 1
  
    for(targetrank=0;targetrank<size;targetrank++){
	int canAccessPeer=0;
	if (targetrank != rank){
	  cudaDeviceCanAccessPeer (&canAccessPeer,rank,targetrank);
	  
	  if (!canAccessPeer){
	    notFullPeerAccess=1;
	    return 1;
	  }
	}
	
    }
    

    //perform initialization
    initIPC(output,roff,mtol_or_ltom);
    already_initialized[mtol_or_ltom]=1;

  }
  else
    {

      //check initialization is consistent with data poiiinter being passed in

      if (output!=outputptrall[mtol_or_ltom][rank]){
	printf("Error: IPC_Alltoall currently only supports being called with 2 distinct output arrays (i.e. from mtol or ltom)\n");

      }

    }


  //Perform AlltoAll: each MPI process pushes data to local or remote location (as already determined
  // by setup of outputptrall array). Each data transfer done in a seperate stream to allow overlap.


  for(targetrank=0;targetrank<size;targetrank++){

    int targetrank_=(rank+targetrank)%size; //for better balance

    cudaMemcpyAsync(&(outputptrall[mtol_or_ltom][targetrank_][roff_remote[mtol_or_ltom][targetrank_]]),&input[soff[targetrank_]],len[targetrank_]*sizeof(double),cudaMemcpyDeviceToDevice,streams[targetrank_]);

  }
  
  // sync all streams
  for(targetrank=0;targetrank<size;targetrank++)
    cudaStreamSynchronize(streams[targetrank]);
  
  
  cudaCheckError();

  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0; 
  
  
}


