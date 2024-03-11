ecTrans
*******

Introduction
============

ecTrans is the global spherical Harmonics transforms library, extracted from the IFS.
It contains both CPU and GPU (Nvidia) code-paths.
The CPU version uses a hybrid of MPI and OpenMP parallelisation strategies, while the GPU version combines MPI and OpenACC.
A default installation builds both CPU libraries (trans_sp, trans_dp) and various flavours of GPU libraries in (trans_gpu_{sp/dp} shared library, trans_gpu_static_{sp/dp} static library, trans_gpu_static_CA_{sp/dp} static library requiring CUDA-aware MPI implementation), as well as a C interface to the double-precision version (transi_dp). A simple benchmark driver is also built against each of these libraries, allowing simple testing of the transforms.

License
=======

ecTrans is distributed under the Apache License Version 2.0.
See `LICENSE` file for details.

Installing ecTrans
==================

Supported Platforms
-------------------

- Linux
- Apple MacOS

Other UNIX-like operating systems may work too out of the box.

The GPU codepath has only been tested with NVHPC compilers on recent Nvidia GPUs.

Requirements
------------
- Fortran compiler with OpenMP support
- C compiler
- FIAT (see https://github.com/ecmwf-ifs/fiat )
- ecBuild (see https://github.com/ecmwf/ecbuild)
- CMake (see https://cmake.org)
- A BLAS library

Further optional recommended dependencies:
- FFTW (http://www.fftw.org)

For the GPU libraries :
- Fortran compiler with OpenACC support
- CUDA toolkit (compiler, and CUBLAS and CUFFT libraries)

Building ecTrans
----------------

Building and installing Trans happens via CMake, which provides automatic detection for
third-party libraries in standard locations and helps cross-plaform portability.
There are multiple ways to help CMake discover packages in non-standard locations.
One explicit way is to e.g. set environment variables for each dependency.

Environment variables 

    $ export ecbuild_ROOT=<path-to-ecbuild>
    $ export fiat_ROOT=<path-to-fiat>
    $ export CC=<path-to-C-compiler>
    $ export FC=<path-to-Fortran-compiler>

You must compile FIAT out-of-source, so create a build-directory (anywhere)

    $ mkdir build && cd build
 
Configuration of the build happens through standard CMake

    $ cmake

Extra options can be added to the `cmake` command to control the build:

 - `-DCMAKE_BUILD_TYPE=<Debug|RelWithDebInfo|Release|Bit>` default=RelWithDebInfo (typically `-O2 -g`)
 - `-DENABLE_TESTS=<ON|OFF>`            default=ON
 - `-DENABLE_SINGLE_PRECISION=<ON|OFF>` default=ON
 - `-DENABLE_DOUBLE_PRECISION=<ON|OFF>` default=ON
 - `-DENABLE_TRANSI=<ON|OFF>`           default=ON
 - `-DENABLE_MKL=<ON|OFF>`              default=ON
 - `-DENABLE_FFTW=<ON|OFF>`             default=ON
 - `-DENABLE_GPU=<ON|OFF>`              default=OFF
 - `-DCMAKE_INSTALL_PREFIX=<install-prefix>`

Specific extra options exist for GPU installation:
 - `-DENABLE_GPU_AWARE_MPI=<ON|OFF>`    default=OF
 - `-DENABLE_GPU_GRAPHS_GEMM=<ON|OFF>`  default=ON
 - `-DENABLE_CUTLASS=<ON|OFF>`          default=OFF
 - `-DENABLE_3XTF32=<ON|OFF>`           default=OFF

GPU-aware MPI allows buffers residing on GPU to be passed to MPI communication calls directly. This requires a compatible MPI installation.
Graph work-flows allow a series of GPU operations to be scheduled in an efficient manner. 
This is useful both for the batched FFTs and the batched GEMMs on which ecTrans relies, although for FFTs this is currently relied upon.
Cutlass is an Nvidia library of templates for GEMM operations. 3xTF32 is a specific acceleration for single precision operations, enabled by Cutlass.

More options to control compilation flags, only when defaults are not sufficient

 - `-DCMAKE_Fortran_FLAGS=<fortran-flags>`
 - `-DCMAKE_C_FLAGS=<c-flags>`

Once this has finished successfully, run ``make`` and ``make install``.

Optionally, tests can be run to check succesful compilation, when the feature TESTS is enabled (`-DENABLE_TESTS=ON`, default ON)

    $ ctest

The benchmark drivers are found in the bin directory.
A brief description of available command-line arguments can be obtained with e.g.
ectrans-benchmark-sp --help

Reporting Bugs
==============

TODO

Contributing
============

TODO

