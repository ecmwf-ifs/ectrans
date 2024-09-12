---
title: User Guide
ordered_subpage: introduction.md
ordered_subpage: installation.md
ordered_subpage: usage.md
ordered_subpage: benchmarking.md
ordered_subpage: api.md
ordered_subpage: transi.md
ordered_subpage: license.md
copy_subdir: img
---

@warning
Page under construction.
@endwarning

## Introduction

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

## The spectral transform

ecTrans transforms a batch of meteorological fields from a grid point space representation
\( X_k(\lambda_i, \phi_j) \), where \( \lambda_i \) is the \( i^{\text{th}} \) longitude,
\( \phi_j \) is the \( j^{\text{th}} \) latitude, and \( k \) is the index which ranges over the
batch of fields, to a spectral space representation \( X_{m,n,k} \), where \( m \) is the zonal
wavenumber, and \( n \) is the total wavenumber. This constitutes a direct spectral transform.
ecTrans can also carry out the inverse spectral transform.

Beginning with the inverse spectral transform (spectral space to grid point space), this is
accomplished in two computational steps. Firstly, an inverse Legendre transform is performed in the
latitudinal direction,

\[
X_{m,k}(\phi_j) = \sum_{n=|m|}^{N} X_{m,n,k} P_{m,n}(\sin(\phi_j)).
\]

Then, an inverse Fourier transform is performed in the longitudinal direction,

\[
X_k(\lambda_i, \phi_j) = \sum_{n=-N}^{N} X_{m,k}(\phi_j) e^{im\lambda_i}.
\]

## Parallelizing a spectral transform

## Basic usage of ecTrans
