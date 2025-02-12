---
title: Design
---

@warning
Page under construction.
@endwarning

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
