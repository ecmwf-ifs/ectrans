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
version: 1.4.0-prerelease
---

## Efficient and _scalable_ spectral transforms.

ecTrans is the global spherical Harmonics transforms library, extracted from the IFS. It is using a hybrid of MPI and OpenMP parallelisation strategies. The package contains both single- and double precision Fortran libraries (trans_sp, trans_dp), as well as a C interface to the double-precision version (transi_dp)

<p><a class="btn btn-primary" href="https://github.com/ecmwf-ifs/ectrans" role="button">Browse source code on GitHub</a></p>