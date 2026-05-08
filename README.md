ecTrans
=======

ecTrans is a library for performing efficient and scalable spectral transformations. It is used for transforming fields from a grid point space on the sphere (e.g. latitude-longitude) to a spectral space based on spherical harmonics (for global transformations) or bifourier harmonics (for limited area transformations), which constitutes a direct transform. A corresponding inverse transform can also be performed. A transform consists of a Fourier transform in the longitudinal direction and either a Legendre transform (global) or another Fourier transform (limited area) in the latitudinal direction. ecTrans can also operate on fields which are distributed across separate MPI tasks and performs the necessary communication to ensure all data needed for a particular transform are resident on a local task.

After co-development as part of the [Integrated Forecasting System (IFS)](https://www.ecmwf.int/en/forecasts/documentation-and-support/changes-ecmwf-model) atmospheric model of the [European Centre for Medium-Range Weather Forecasts](https://www.ecmwf.int/) for several decades, ecTrans became a standalone software package in 2022. It constitutes one of the most important and expensive parts of the IFS and neatly encapsulates both computational and communicational paradigms and bottlenecks exhibited by the IFS model as a whole.

ecTrans primarily targets conventional CPU platforms, requiring FFTW- and BLAS-implementing libraries. It can also operate efficiently on GPU accelerators making use of offloading directives (either OpenACC or OpenMP) and vendor library routines (cuBLAS/cuFFT or hipBLAS/hipFFT). ecTrans performs efficiently and stably on Nvidia platforms but is currently less mature on AMD platforms.

To learn more about ecTrans, please consult the [documentation](https://sites.ecmwf.int/docs/ectrans/page/index.html) (which is under construction).

License
-------

ecTrans is distributed under the Apache License Version 2.0.
See `LICENSE` file for details.

Requirements
------------

Generally, ecTrans has the following requirements:
- [CMake](https://cmake.org/) >= 3.25
- [ecBuild](https://github.com/ecmwf/ecbuild) >= 3.8
- C, C++, and Fortran compilers. Officially we support:
  - Classic Intel (i.e. ifort and icc) >= 19.0.5
  - LLVM Intel (i.e. ifx) >= 2023.2.0
  - GNU Compiler Collection >= 8.5.0
  - NVHPC >= 22.11
  - Cray Compiler Environment >= 19.0.0
  - AMD ROCm AFAR >= 22.3.0
  Earlier versions may work just fine, but without the means to test these versions, we cannot offer
  support.
- [FIAT](https://github.com/ecmwf-ifs/fiat) >= 2.0.0 (earlier versions are likely to work, but we
  only offer support for the latest version of FIAT, since it is straightforward to build.)
- A library implementing standard BLAS routines

Builds targeting CPU execution have the following additional requirement:
- An FFTW-compatible library, such as [FFTW itself](https://www.fftw.org/) or Intel MKL.

Builds targeting GPU execution have the following additional requirements:
- A compiler compatible with OpenACC offload or OpenMP offload.
- CUDA or HIP.

Installing ecTrans
------------------

Please consult the [documentation](https://sites.ecmwf.int/docs/ectrans/page/installation.html).

## bundle build, text to be moved to sites.ecmwf.int/[]/installation.html
Build using ecTrans bundle

Another way of building ecTrans is to use the bundle definition included in `package/bundle`:

    $ ./package/bundle/ectrans-bundle create  --bundle package/bundle/bundle.yml # Checks out dependency packages into "source" directory
    $ ./package/bundle/ectrans-bundle build [--build-type=<build-type>] [--arch=<path-to-arch>] [--option]

The bundle also facilitates setting environment variables and compiler flags relevant to certain architectures by specifying the corresponding arch file at the build step. For example, to build on the ECMWF Atos system using Intel compilers and the hpcx-openmpi `MPI` library:

`--arch=package/bundle/arch/ecmwf/hpc2020/intel/2021.4.0/hpcx-openmpi/2.9.0`
or equivalently
`--arch ecmwf/hpc2020/intel/2021.4.0/hpcx-openmpi/2.9.0`

A number of options can also be configured during the bundle build step, including:
 - `--without-mpi` - Disable MPI
 - `--without-single-precision` - Only build double-precision variant of ecTrans
All options can be seen at the end of the `package/bundle/bundle.yml` file.

 Finally, additional `CMake` options can also be set during the bundle build step:

`--cmake="OPTION=<arg>"`

Reporting Bugs
--------------

Please report bugs using a [GitHub issue](https://github.com/ecmwf-ifs/ectrans/issues). Support is given on a best-effort basis by package developers.

Contributing
------------

Contributions to ecTrans are welcome. In order to do so, please open a [GitHub issue](https://github.com/ecmwf-ifs/ectrans/issues) where a feature request or bug can be discussed. Then create a [pull request](https://github.com/ecmwf-ifs/ectrans/pulls) to the develop branch (not the main branch) with your contribution. All contributors to the pull request need to sign the [contributors license agreement (CLA)](https://bol-claassistant.ecmwf.int/ecmwf-ifs/ectrans).

