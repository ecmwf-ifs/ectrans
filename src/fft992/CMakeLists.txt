# fft992 library
# Legacy FFT implementations (FFT992 and Bluestein algorithms)

set(fft992_sources
  fft992.F90
  fft992_cc.F90
  set99.F90
  set99b.F90
  bluestein_mod.F90
)

foreach(prec dp sp)
  if(HAVE_${prec})
    set(target fft992_${prec})

    ecbuild_add_library(
      TARGET ${target}
      SOURCES ${fft992_sources}
      PUBLIC_LIBS parkind_${prec} fiat
      INSTALL_HEADERS LISTED
      HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/fft992
    )

    ecbuild_target_fortran_module_directory(
      TARGET ${target}
      MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${target}
      INSTALL_MODULE_DIRECTORY module/${target}
    )

    target_include_directories(${target} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/ectrans>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/fiat>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/${target}>
      $<INSTALL_INTERFACE:module/ectrans>
      $<INSTALL_INTERFACE:module/fiat>
      $<INSTALL_INTERFACE:module/${target}>
    )
  endif()
endforeach()
