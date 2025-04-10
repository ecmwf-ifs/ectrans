# (using CMAKE_CURRENT_SOURCE_DIR is necessary because sources are in a different directory than the target library (trans_${prec}))
ecbuild_list_add_pattern(
    LIST ectrans4py_src
    GLOB ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
    QUIET
)

# Add ectrans4py library
ecbuild_add_library(
      TARGET           ectrans4py_dp
      LINKER_LANGUAGE  Fortran
      SOURCES          ${ectrans4py_src}
      PUBLIC_LIBS      fiat trans_dp parkind_dp etrans_dp
)
