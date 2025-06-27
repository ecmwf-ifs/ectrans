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

Beginning with the direct spectral transform (grid point space to spectral space), this is
accomplished in two computational steps. Firstly, a Fourier transform is performed in the
longitudinal direction at each latitude and for each field, independently,

\[
X_{m,k}(\phi_j) = \sum_{i=0}^{N_j} X_k(\lambda_i, \phi_j) \exp{\left(-I2\pi\frac{m}{N_j}i\right)}.
\]

Here, \(N_j\) is the number of longitudinal points at latitude \(j\) and \(I\) is the imaginary
identity. The input meteorological fields, \( X_k(\lambda_i, \phi_j) \), are always entirely real.
Thus, although a full Fourier transform would be evaluated for \( (m = 0, N_j - 1) \), terms
\( (m = N_j / 2 + 1, N_j - 1) \) are simply the conjugate of terms \( (m = 1, N_j / 2) \) so we
needn't calculate them. This is referred to as an "RFFT" algorithm.

The second step is the Legendre transform which is performed in the latitudinal direction. This is
performed by Gaussian quadrature with Gaussian weights \( w_j \),

\[
    X_{m,n,k} = \sum_{j=0}^{N_\text{lat}} X_{m,n}(\phi_j) w_j P_{m,n}(\sin(\phi_j)),
\]

where \( P_{m,n} (\sin(\phi_j)) \) is the associated Legendre polynomial of zonal wavenumber \(m\)
and total wavenumber \(n\) which is evaluated at the point \(\sin(\phi_j)\). The Legendre transform
is implemented as a matrix-matrix multiplication, collapsing over the latitudinal dimension with the
fields and total wavenumber dimensions "free". This operation is performed independently at each
zonal wavenumber \(m\). Notably, the fact that \( P_{m,n} \)s is symmetric when \(n - m\) is even
and is antisymmetric when \(n - m\) is odd is exploited to halve the required number of floating-
point operations. The Legendre transform is split into two transforms over the Northern hemisphere:
one for each of the even and odd polynomials.


<!-- Beginning with the inverse spectral transform (spectral space to grid point space), this is
accomplished in two computational steps. Firstly, an inverse Legendre transform is performed in the
latitudinal direction,

\[
X_{m,k}(\phi_j) = \sum_{n=m}^{N} X_{m,n,k} P_{m,n} (\sin(\phi_j)),
\]

where \( P_{m,n} (\sin(\phi_j)) \) is the associated Legendre polynomial of zonal wavenumber \(m\)
and total wavenumber \(n\) which is evaluated at the point \(\sin(\phi_j)\). This operation is
performed independently for each zonal wavenumber, \(m\), and each field, \(k\). Zonal wavenumbers
range from \(0\) up to \(N_{\text{trunc}}\), where the latter denotes the highest zonal and total
wavenumber. For each \(m\), the inverse Legendre transform is performed efficiently as a matrix-
matrix multiplication.

Then, an inverse Fourier transform is performed in the longitudinal direction,

\[
X_k(\lambda_i, \phi_j) = \sum_{m=0}^{N} X_{m,k}(\phi_j) e^{im\lambda_i}.
\]

This operation is performed independently at each latitude using a fast Fourier transform (FFT)
algorithm. Note that, as suggested by the summation indices, a "complex-to-real" inverse FFT
algorithm is used. This means that negative zonal wavenumbers do not need to be considered. This is
because the field values in grid point space (\(X_k\)) are always all real in a meterological
context. -->

## Parallelizing a spectral transform

## Basic usage of ecTrans
