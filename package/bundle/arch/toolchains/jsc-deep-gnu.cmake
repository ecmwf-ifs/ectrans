####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI ON )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffpe-trap=invalid,zero,overflow")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fstack-arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fconvert=big-endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fbacktrace")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fno-second-underscore")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffree-form")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffast-math")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fno-unsafe-math-optimizations")
