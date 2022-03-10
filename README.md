ecTrans
*******

Introduction
============

ecTrans is the global spherical Harmonics transforms library, extracted from the IFS.
It contains both CPU and GPU (Nvidia) code-paths.
The CPU version uses a hybrid of MPI and OpenMP parallelisation strategies, while the GPU version combines MPI and OpenACC.
A default installation builds both CPU libraries (trans_sp, trans_dp) and various flavours of GPU libraries in (trans_gpu_{sp/dp} shared library, trans_gpu_static_{sp/dp} static library, trans_gpu_static_CA_{sp/dp} static library requiring CUDA-aware MPI implementation),
as well as a C interface to the double-precision version (transi_dp).

License
=======

Trans is distributed under the Apache License Version 2.0.
See `LICENSE` file for details.

Installing ecTrans
==================

Supported Platforms
-------------------

The GPU codepath has only been tested with NVHPC compilers on recent Nvidia GPUs.

Requirements
------------
- MPI-enabled Fortran compiler with OpenACC support
- FIAT (see https://github.com/ecmwf-ifs/fiat )
- ecBuild (see https://github.com/ecmwf/ecbuild)
- CMake (see https://cmake.org)
- CUDA toolkit (compiler, and CUBLAS and CUFFT libraries)
- ecCodes (see https://github.com/ecmwf/eccodes ) if building the drivers

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
 -  -DENABLE_GPU=ON
 -  -DENABLE_FFTW=OFF
 - `-DCMAKE_INSTALL_PREFIX=<install-prefix>`

More options to control compilation flags, only when defaults are not sufficient

 - `-DCMAKE_Fortran_FLAGS=<fortran-flags>`
 - `-DCMAKE_C_FLAGS=<c-flags>`

Once this has finished successfully, run ``make`` and ``make install``.

Optionally, tests can be run to check succesful compilation, when the feature TESTS is enabled (`-DENABLE_TESTS=ON`, default ON)

    $ ctest

Reporting Bugs
==============

TODO

Contributing
============

TODO

