---
title: Introduction
---

@warning
Page under construction.
@endwarning

ecTrans is a high-performance numerical library for transforming meteorological fields between
global grid-point space representation and a spectral representation based on spherical harmonics.
It consists of a Fourier transform in the longitudinal (east-west) direction and a Legendre
transform in the latitudinal (north-south) direction. It also consists of transposition routines interleaved with these two transforms to ensure that the data required has the correct layout.

ecTrans is a fundamental part of the European Centre for Medium-Range Weather Forecasts' (ECMWF's)
Integrated Forecasting System (IFS), a global numerical weather prediction suite. Indeed, ecTrans
was previously part of the IFS source code itself. It therefore benefits from over 30 years of
development and optimisation. In 2022, ecTrans was split out from the IFS source code and released
as the first open-source component of the IFS as its own project.

## The spectral transform algorithm



## Parallelizing a spectral transform

## Basic usage of ecTrans
