---
project: ecTrans
src_dir: ../src
exclude_dir: ../src/transi
page_dir: content_processed
output_dir: site
css: css/ectrans.css
author: ECMWF
graph: true
website: https://www.ecmwf.int
favicon: img/favicon.ico
privacy_policy_url: https://www.ecmwf.int/en/privacy
terms_of_service_url: https://www.ecmwf.int/en/terms-use
version: 1.7.0
html_template_dir: html
---

## Efficient and _scalable_ spectral transforms.

ecTrans is a library for performing efficient and scalable spectral transformations. It is used for
transforming fields from a grid point space on the sphere (e.g. latitude-longitude) to a spectral
space based on spherical harmonics (for global transformations) or bifourier harmonics (for limited
area transformations), which constitutes a direct transform. A corresponding inverse transform can
also be performed. A transform consists of a Fourier transform in the longitudinal direction and
either a Legendre transform (global) or another Fourier transform (limited area) in the latitudinal
direction. ecTrans can also operate on fields which are distributed across separate MPI tasks and
performs the necessary communication to ensure all data needed for a particular transform are
resident on a local task.

After co-development as part of the Integrated Forecasting System (IFS) atmospheric model of the
European Centre for Medium-Range Weather Forecasts for several decades, ecTrans became a standalone
software package in 2022. It constitutes one of the most important and expensive parts of the IFS
and neatly encapsulates both computational and communicational paradigms and bottlenecks exhibited
by the IFS model as a whole.

ecTrans primarily targets conventional CPU platforms, requiring FFTW- and BLAS-implementing
libraries. It can also operate efficiently on GPU accelerators making use of offloading directives
(either OpenACC or OpenMP) and vendor library routines (cuBLAS/cuFFT or hipBLAS/hipFFT). ecTrans
performs efficiently and stably on Nvidia platforms but is currently less mature on AMD platforms.
