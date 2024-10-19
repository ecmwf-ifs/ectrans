if(HAVE_ETRANS)
    # (using CMAKE_CURRENT_SOURCE_DIR is necessary because sources are in a different directory than the target library (trans_${prec}))
    ecbuild_list_add_pattern(
        LIST ectrans4py_src
        GLOB ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
        QUIET
    )

    set(HAVE_dp ${HAVE_DOUBLE_PRECISION})
    set(prec dp)

    if(HAVE_${prec})
        # Add sources
        target_sources(trans_${prec} PRIVATE ${ectrans4py_src})
    endif()

else()
    ecbuild_critical("To activate the ectrans Python interface, you must enable the ETRANS option.")
endif()

