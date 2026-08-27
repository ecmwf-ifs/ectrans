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
- [ecBuild](https://github.com/ecmwf/ecbuild) >= 3.14
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

Please consult the [documentation](https://sites.ecmwf.int/docs/ectrans/page/installation.html)
for general instructions on required dependencies, their installation, and use of package managers.


Building ecTrans
----------------

ecTrans can be built either as a standalone CMake project, using already installed dependencies,
or through the bundle helper in [`package/bundle`](package/bundle), which can build ecTrans and
its bundled dependencies together.

### Standalone repository

For a standalone build, provided all dependencies can be found via CMake search paths:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Bit
cmake --build build
ctest --test-dir build --output-on-failure
```

Common configuration options include:

- `-DENABLE_ETRANS=ON` to build the limited-area transform implementation.
- `-DENABLE_FIELD_API=ON` to build the field-api interface.
- `-DENABLE_GPU=ON` with `-DENABLE_ACC=ON` or `-DENABLE_OMP=ON` for GPU builds.
- `-DENABLE_MKL=ON` to use Intel MKL where available.
- `-DENABLE_FFTW=ON` or package/toolchain settings to select an FFTW-compatible library.
  (e.g. when MKL is available MKL will be used as it contains an FFTW-compatibility layer)

Options that configure testing include:

- `-DENABLE_TESTS=ON|OFF` to disable all tests. These are enabled by default.
- `-DENABLE_BITID_TESTS=ON|OFF` to enable or disable bit-identicality tests. These are disabled
  by default.
- `-DENABLE_REFERENCE_TESTS=ON|OFF` to enable or disable comparisons with stored checksum
  references. These are disabled by default.
- `-DECTRANS_TEST_REFERENCE_DIR=<reference-dir>` to use a non-default stored reference directory.

### Bundle build

The repository also contains an ecbundle wrapper:

```sh
package/bundle/ectrans-bundle create [options]
package/bundle/ectrans-bundle build --build-dir=<build-dir> [options]
```

See [`package/bundle/README.md`](package/bundle/README.md) for bundle-specific options,
architecture environment files, and example commands.

Testing ecTrans
---------------

After configuring and building with tests enabled, run the test suite with CTest:

```sh
ctest --test-dir <build-dir> --output-on-failure
```

When `BITID_TESTS` feature is enabled, extra tests are present that verify
bit-identicality between test cases that are expected to produce identical checksum files. This
can only be guaranteed with BLAS libraries and runtime settings that provide reproducible results.
One good choice is the [reference BLAS implementation from Netlib](https://www.netlib.org/blas/).
Another known configuration is Intel MKL 19.0.5 with:

```sh
export MKL_CBWR=AUTO,STRICT
```

To enable these tests when the BLAS/runtime combination is bit-reproducible, configure with:

```sh
cmake -S <src-dir> -B <build-dir> -DENABLE_BITID_TESTS=ON
```

Stored reference tests can be enabled separately with the feature `REFERENCE_TESTS`:

```sh
cmake -S <src-dir> -B <build-dir> -DENABLE_REFERENCE_TESTS=ON
```

By default, ecTrans chooses a reference directory based on the configured build. To configure a
non-default directory, set the CMake variable `ECTRANS_TEST_REFERENCE_DIR`:

```sh
cmake -S <src-dir> -B <build-dir> \
  -DENABLE_REFERENCE_TESTS=ON \
  -DECTRANS_TEST_REFERENCE_DIR=<reference-dir>
```

The `ECTRANS_TEST_REFERENCE_DIR` environment variable can also be set when running tests or
updating references to select a different `<reference-dir>` than the one configured by CMake.

Reference tests compare generated checksum files with stored references and skip individual
comparisons when the corresponding reference file does not exist.

To generate and store the references into the configured `<reference-dir>`, we can use the build-target `ectrans-update-references`:

```sh
cmake --build <build-dir> --target ectrans-update-references
```

Here it is possible to override the configured `<reference-dir>` to a `<custom-reference-dir>`
via environment variable `ECTRANS_TEST_REFERENCE_DIR`:

```sh
ECTRANS_TEST_REFERENCE_DIR=<custom-reference-dir> \
cmake --build <build-dir> --target ectrans-update-references
```

This also provides a way to check bit-reproducibility across branches:
generate or update references on control branch, then configure another branch with
`-DENABLE_REFERENCE_TESTS=ON` and point the `ECTRANS_TEST_REFERENCE_DIR` variable
to the control `<reference-dir>`.

Reporting Bugs
--------------

Please report bugs using a [GitHub issue](https://github.com/ecmwf-ifs/ectrans/issues). Support is given on a best-effort basis by package developers.

Contributing
------------

Contributions to ecTrans are welcome. In order to do so, please open a [GitHub issue](https://github.com/ecmwf-ifs/ectrans/issues) where a feature request or bug can be discussed. Then create a [pull request](https://github.com/ecmwf-ifs/ectrans/pulls) to the develop branch (not the main branch) with your contribution. All contributors to the pull request need to agree to the [contributors license agreement (CLA)](contributor_license_agreement.md).

