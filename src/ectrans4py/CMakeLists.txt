# (using CMAKE_CURRENT_SOURCE_DIR is necessary because sources are in a different directory than the target library (trans_${prec}))
ecbuild_list_add_pattern(
    LIST ectrans4py_src
    GLOB ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
    QUIET
)

# Add sources
target_sources(trans_dp PRIVATE ${ectrans4py_src})

