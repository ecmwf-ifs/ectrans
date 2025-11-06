---
title: CMake Features
---

Different build-time features of ecTrans are enabled or disabled through CMake's "features",
functionality. Here we list all of the features which can be enabled by passing
`-DENABLE_<feature name>=ON` (and disabled by specifying `OFF`).

- `MPI`: Enable Message Passing Interface functionality, required for distributed-memory
  parallelism.
- `OMP`: Enable OpenMP parallelism. Default `ON`. (Note: currently this controls both OpenMP multithreading for the
  CPU version of ecTrans but also OpenMP offloading for the GPU version. In future we may split
  these into two separate features.)
- `ACC`: Enable OpenACC parallelism. Default `ON`.
- `DOUBLE_PRECISION`: Enable double-precision build of ecTrans. Default `ON`.
- `SINGLE_PRECISION`: Enable single-precision build of ecTrans. Default `ON`.
- `MKL`: Use Intel's Math Kernel Library for the BLAS and FFTW calls in ecTrans. Default `ON`.
- `CPU`: Enable CPU build of ecTrans. Default `ON`.
- `TRANSI`: Enable transi C interface for ecTrans. Requires `DOUBLE_PRECISION` and `CPU` to be
  enabled. Default `ON`.
- `GPU`: Enable GPU build of ecTrans. Requires `ACC` or `OMP` to be enabled. Default `OFF`.
- `CUTLASS`: Use Nvidia's CUTLASS library to accelerate the Legendre transform BLAS operations.
  Requires `GPU` and `CUDA` to be enabled. Default `OFF`.
- `CUTLASS_3XTF32`: Use 3XTF32 kernel to accelerate the Legendre transform GEMM. Requires
  `SINGLE_PRECISION` and `CUTLASS` to be enabled. Default `ON`.
- `GPU_AWARE_MPI`: Use direct GPU-GPU communication for the transpose routines. Requires `GPU` to be
  enabled. Default `ON`.
- `GPU_GRAPHS_GEMM`: Enable graph optimisation of Legendre transform GEMM kernels. Requires `GPU` to
  be enabled. Default `ON`.
- `GPU_GRAPHS_FFT`: Enable graph optimisation of FFT kernels. Requires `GPU` to be enabled. Default
  `ON`.
- `GPU_STATIC`: Compile GPU library as a static library. Default `ON` when `BUILD_SHARED_LIBS`
  evalutes to true, otherwise `OFF`.
- `ETRANS`: Build limited-area version of ecTrans ("etrans"). Default `OFF`.
- `ECTRANS4PY`: Build ectrans4py Python interface to ecTrans. Requires `ETRANS` and
  `DOUBLE_PRECISION` to be enabled. Default `OFF`.