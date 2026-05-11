---
title: Installation
---

We currently support three methods for installing ecTrans:

1. Through [Spack](https://github.com/spack/spack).
2. Manually through CMake.
3. Through [ecBundle](https://github.com/ecmwf/ecbundle) (develop branch only until the next release).

# 1. Installing ecTrans with Spack

A recipe for building ecTrans is provided in the main [spack-packages repository](https://github.com/spack/spack-packages/tree/develop/repos/spack_repo/builtin/packages/ectrans).
For most users a simple
```bash
spack install ectrans
```
should be sufficient. Note that, by default,
ecTrans' dependency FIAT will be built with fckit support. This is not necessary for ecTrans and
brings in several other dependencies, so it's recommend to disable it when installing ecTrans:
```bash
spack install ectrans %fiat~fckit
```

# 2. Installing ecTrans manually

ecTrans relies on CMake for building and for unit testing. It follows standard CMake procedures for
building out-of-source, but below we describe this explicitly for newcomers to CMake.

## Requirements

ecTrans has the following requirements:

- A [CMake](https://cmake.org/) with version >= 3.25
- [ecBuild](https://github.com/ecmwf/ecbuild.git) (a collection of ECMWF-specific CMake macros)
- A Fortran compiler with OpenMP support and a C compiler. Officially we support:
    - Classic Intel (i.e. ifort and icc)
    - LLVM Intel (i.e. ifx)
    - GNU
    - NVHPC
    - Cray
    - AMD ROCm AFAR
- [FIAT: the Fortran IFS and Arpege Toolkit](https://github.com/ecmwf-ifs/fiat)
- [FFTW](https://www.fftw.org/)
- A library containing standard BLAS routines such as [LAPACK](https://www.netlib.org/lapack/)

Note that you will also of course need a MPI library if you want to run ecTrans with distributed
memory parallelism, but this is not an explicit dependency of ecTrans. Instead MPI functionality is
provided through a wrapper library provided by FIAT. It is the latter that must be built with MPI
support. In any case, ecTrans can still be tested without MPI with only shared memory parallelism
support.

For all of these except for FIAT and ecBuild we cannot give general instructions as it depends
entirely on your system and software environment. Most modern high-performance computing systems
should have these installed already, so we will assume this is the case for you as well.

ecBuild can simply be cloned from GitHub like so:

```bash
git clone https://github.com/ecmwf/ecbuild.git --branch 3.12.0 --single-branch
```

It does not require a build or installation step. Simply export a variable `ecbuild_DIR` pointing to
the cloned repository. 

@note
We are always willing to add other compilers to our list, so please
[raise an issue](https://github.com/ecmwf-ifs/ectrans/issues) if you are encountering
difficulties with a particular compiler.
@endnote

## (Optional) Prepare env and toolchain files

Preparing these files is not strictly necessary, but is recommended to ease reproducibility. An env
file (called env.sh below) should contain all of the steps required to make the dependencies of FIAT
and ecTrans available on the path, and otherwise prepare the environment for their building and
execution. For example, the env file on a system that uses modules might load the modules for CMake,
the compiler suite, FFTW, and the BLAS library. The toolchain file (called toolchain.cmake below) 
contains definitions of CMake variables which will be needed to build FIAT and ecTrans correctly.
The most important are the compilers. The toolchain file should set the values of
`CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`, and `CMAKE_Fortran_COMPILER`. The toolchain file for
NVHPC will look like this, for example:

```cmake
set(CMAKE_C_COMPILER nvc)
set(CMAKE_CXX_COMPILER nvc++)
set(CMAKE_Fortran_COMPILER nvfortran)
```

To get CMake to always detect ecBuild, add this to your env file:

```bash
export ecbuild_DIR=<path to ecBuild>
```

We will assume both the env file (env.sh) and the toolchain file (toolchain.cmake) are in the
current directory. Now source the env file:

```bash
source env.sh
```

## Building FIAT

First clone the latest version of the FIAT repository:

```bash
git clone https://github.com/ecmwf-ifs/fiat.git -b 1.6.1
```

Then run the configure step for FIAT (you can leave out `-DCMAKE_TOOLCHAIN_FILE` if you don't want
to bother with toolchain files):

```bash
cmake -S fiat -B fiat/build -DCMAKE_TOOLCHAIN_FILE=$PWD/toolchain.cmake
```

Now run the build step (adjusting the level of multithreading as required with `-j`):

```bash
cmake --build fiat/build -j 4
```

## Building ecTrans

Clone the latest version of the ecTrans repository:

```bash
git clone https://github.com/ecmwf-ifs/ectrans.git -b 1.6.2
```

Then run the configure step, making sure to pass the location of the FIAT build directory to CMake:

```bash
cmake -S ectrans -B ectrans/build -DCMAKE_TOOLCHAIN_FILE=$PWD/toolchain.cmake \
  -Dfiat_ROOT=$PWD/fiat/build
```

Now run the build step:

```bash
cmake --build ectrans/build -j 4
```

You should now find in the ecTrans build directory the ecTrans libraries in the lib subdirectory
(single- and double-precision versions) and the standalone ecTrans benchmarking binary in the bin
subdirectory. We strongly recommend you run the full CTest suite to check everything's installed
correctly:

```bash
cd ectrans/build
ctest
```
If any of the tests fail, please let us know by
[raising an issue](https://github.com/ecmwf-ifs/ectrans/issues)!

# 3. Installing ecTrans through ecBundle

@note
This functionality is available in the develop branch only at this stage, until the next release of
ecTrans is made.
@endnote

[ecBundle](https://github.com/ecmwf/ecbundle) is a package developed by ECMWF for "bundling"
separate CMake projects into one, with simple commands for fetching source code and building the
unified project. The unified project is defined in a bundle YAML file, usually called `bundle.yml`.
This defines all of the separate CMake repositories, their git versions (e.g. a branch or commit
hash), URLs for fetching those repositories, and defines any command line arguments that can be used
to enable or disable certain CMake features.

ecTrans comes with a bundle in the `package` subdirectory. To build ecTrans in this
manner, simply call the "create" step to first fetch the source code required (including ecTrans'
dependencies ecBuild and FIAT):

```bash
./package/bundle/ectrans-bundle create --bundle ./package/bundle/bundle.yml --src-dir=./package/bundle/source
```

This checks out dependency packages into a directory `source` in the same directory as `bundle.yml`.

You can then build ecTrans by calling the `build` subcommand:

```bash
./package/bundle/ectrans-bundle build --build-dir=./package/bundle/build --src-dir=./package/bundle/source [--build-type=<build-type>] [--arch=<path-to-arch>] [--any --other --options]
```

The build type refers to the CMake build type (`BIT`, `RELEASE` etc.). The `arch` option allows one
to specify environment and toolchain files for the platform on which you would like to build. For
example, archs for ECMWF's current HPC system are provided in
`package/bundle/ectrans-bundle/arch/ecmwf/hpc2020`. To build on the ECMWF Atos system using Intel
compilers and the hpcx-openmpi `MPI` library, use:

```bash
--arch=package/bundle/arch/ecmwf/hpc2020/intel/2021.4.0/hpcx-openmpi/2.9.0
```

or equivalently

```bash
--arch ecmwf/hpc2020/intel/2021.4.0/hpcx-openmpi/2.9.0
```

Following these one may provide additional arguments to enable or disable certain features, such as
`--without-mpi` (disable MPI), `--without-single-precision` (only build double-precision variant of
ecTrans) etc. These are all defined in `bundle.yml`.

 Finally, additional `CMake` options can also be set during the bundle build step:

`--cmake="OPTION=<arg>"`

Once the build completes you can navigate to the build directory and run `ctest`. This will execute
a combined test suite for all packages in the bundle (in this case FIAT and ecTrans).
