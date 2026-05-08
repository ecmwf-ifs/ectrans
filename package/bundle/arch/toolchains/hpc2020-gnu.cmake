####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "-fconvert=big-endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffree-line-length-none")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fbacktrace")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffpe-trap=invalid,overflow,zero")

set(ECBUILD_Fortran_FLAGS_BIT "-g -O2")
set(ECBUILD_C_FLAGS_BIT "-g -O2")
set(ECBUILD_CXX_FLAGS_BIT "-g -O2")

