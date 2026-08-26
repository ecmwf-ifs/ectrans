ecTrans Bundle
==============

The ecTrans bundle uses [ecbundle](https://github.com/ecmwf/ecbundle) to build ecTrans together
with the dependencies described in [`bundle.yml`](bundle.yml). The wrapper script
[`ectrans-bundle`](ectrans-bundle) downloads ecbundle locally when it is not already available and
then forwards all arguments to `ecbundle`.

Quick Start
-----------

Run the bundle commands from the ecTrans source tree:

```sh
package/bundle/ectrans-bundle create [options]
package/bundle/ectrans-bundle build [options]
```

The `create` command prepares the bundle build tree and dependency checkout layout. The `build`
command configures, builds, and installs the projects described by the bundle.

For example:

```sh
package/bundle/ectrans-bundle create --bundle package/bundle/bundle.yml --src-dir=package/bundle/source
package/bundle/ectrans-bundle build --with-etrans --arch=ecmwf/hpc2020 --src-dir=package/bundle/source --build-dir=package/bundle/build
```

Architecture Environments
-------------------------

Platform-specific environment and toolchain files are kept under [`arch`](arch). On supported
systems, source the matching environment file before creating or building the bundle. For example:

```sh
source package/bundle/arch/ecmwf/hpc2020/intel/2023.2.0/hpcx-openmpi/2.9.0/env.sh
package/bundle/ectrans-bundle create --bundle package/bundle/bundle.yml --src-dir=package/bundle/source
package/bundle/ectrans-bundle build --with-mkl --with-etrans --arch=ecmwf/hpc2020/intel/2023.2.0/hpcx-openmpi/2.9.0 --src-dir=package/bundle/source --build-dir=package/bundle/build
```

The architecture environment files may set compiler wrappers, module settings, CMake toolchain
files, and reproducibility-related runtime variables such as `MKL_CBWR=AUTO,STRICT`.

Bundle Options
--------------

The available bundle options are defined in [`bundle.yml`](bundle.yml). Common options include:

- `--with-etrans` to enable the limited-area transform implementation.
- `--with-mkl` to enable Intel MKL for BLAS and FFT calls.
- `--with-fftw` to enable FFTW calls.
- `--without-mpi` to disable MPI.
- `--without-omp` to disable OpenMP.
- `--with-bitid-tests` to enable bit-identicality tests.
- `--with-reference-tests` to enable comparisons with stored checksum references.
- `--reference-test-dir=<reference-dir>` to use a non-default stored reference directory.
- `--without-single-precision` or `--without-double-precision` to build only one precision.
- `--with-gpu-acc` or `--with-gpu-omp` to build GPU libraries with OpenACC or OpenMP offload.
- `--without-gpu-mpi` to disable GPU-aware MPI.
- `--with-cutlass` to enable CUTLASS GEMM kernels on supported CUDA/NVHPC builds.
- `--without-cutlass-3XTF32` to disable 3XTF32 usage on compatible CUDA/NVHPC builds.
- `--without-gpu-graphs` to disable GPU graph optimisations.
- `--with-static-linking` to build static libraries by default.
- `--with-rocfft` to use rocFFT on AMD GPU builds.

Run the wrapper with ecbundle help options, or inspect [`bundle.yml`](bundle.yml), for the full
set of supported options in the current source tree.

Testing
-------

After a bundle build, go to the generated build directory and source its `env.sh` file before
running tests. The exact build directory name depends on the ecbundle configuration:

```sh
cd <build-dir>
source env.sh
ctest -R ectrans --output-on-failure
```

Run CTest as described in the top-level [`README.md`](../../README.md)