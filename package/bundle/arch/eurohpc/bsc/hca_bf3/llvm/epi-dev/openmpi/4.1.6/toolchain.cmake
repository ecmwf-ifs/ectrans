# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_CXX_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-fopenmp" CACHE STRING "" )
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -O3 -ffast-math -mepi -mllvm -combiner-store-merging=0 -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -mllvm -vectorizer-use-vp-strided-load-store -mcpu=avispado -mllvm -disable-loop-idiom-memcpy -mllvm -disable-loop-idiom-memset -Rpass-missed=loop-vectorize -Xflang -target-feature -Xflang +does-not-implement-vszext -Xflang -target-feature -Xflang +does-not-implement-tu -mllvm -riscv-uleb128-reloc=0 -fno-slp-vectorize")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -I/apps/x86/rave/rave-latest/interfaces") # Make RAVE available (lib is also added in cmake)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -I/apps/x86/rave/rave-frozen-09-2025/interfaces") # Make RAVE available (lib is also added in cmake)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -I/apps/riscv/ubuntu/openmpi/4.1.6_llvm1.0/lib") # Quick hack to satisfy broken libs
####################################################################
# COMMON FLAGS
####################################################################

set(CMAKE_EXE_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -Wl,--whole-archive /apps/riscv/llvm/EPI/development/lib/libprovector-vecclonevp.a -Wl,--no-whole-archive -fopenmp") # Enable the vector maths library
set(CMAKE_SHARED_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")

