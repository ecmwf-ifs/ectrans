---
title: ecTrans User Guide
ordered_subpage: design.md
ordered_subpage: installation.md
ordered_subpage: usage.md
ordered_subpage: benchmarking.md
ordered_subpage: api.md
ordered_subpage: transi.md
ordered_subpage: gpu.md
ordered_subpage: license.md
copy_subdir: img
---

ecTrans is a high-performance numerical library for transforming meteorological fields between
global grid-point space representation and a spectral representation based on spherical harmonics.
It is a fundamental part of the
[European Centre for Medium-Range Weather Forecasts'](https://www.ecmwf.int/)
[Integrated Forecasting System (IFS)](https://www.ecmwf.int/en/forecasts/documentation-and-support/
changes-ecmwf-model), a global numerical weather prediction suite. Indeed, ecTrans
was previously part of the IFS source code itself. It therefore benefits from over 30 years of
development and optimisation. In 2022, ecTrans was split out from the IFS source code and released
as the first open-source component of the IFS as its own project.

ecTrans is engineered to work efficiently running on many hundreds, or even thousands, of compute
nodes. This is achieved through a significant optimisation of the constituent compute kernels,
making use of the FFTW library for the Fourier transform in the longitudinal direction and BLAS
GEMMs for the Legendre transforms in the latitudinal direction. However, given that transformed
fields are distributed across compute tasks, great care has been taken to ensure that parallelism
can be fully exploited at all stages in the algorithm. This is achieved through data exchange steps
interleaved between the Fourier and Legendre transforms, which are implemented using the Message
Passing Interface (MPI).

The result is an algorithm which stresses a high-performance computing system both on the
node level _and_ the network level, serving as an excellent overall benchmark and a target for
optimisation of IFS execution speed.

This user guide contains the following sections:

- [Design](design.html)
- [Installation](installation.html)
- [Usage](usage.html)
- [Benchmarking](benchmarking.html)
- [API](api.html)
- [Interfacing with C](transi.html)
- [GPU offloading](gpu.html)
- [License](license.html)