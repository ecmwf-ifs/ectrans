---
project: ecTrans
src_dir: ../src
exclude_dir: ../src/transi
page_dir: content
output_dir: site
project_github: https://github.com/ecmwf-ifs/ectrans
css: css/ectrans.css
author: ECMWF
graph: true
---

ecTrans is the global spherical Harmonics transforms library, extracted from the IFS. It is using a hybrid of MPI and OpenMP parallelisation strategies. The package contains both single- and double precision Fortran libraries (trans_sp, trans_dp), as well as a C interface to the double-precision version (transi_dp)
