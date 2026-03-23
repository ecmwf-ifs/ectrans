# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  if [ "${2:-none}" == "ECBUILD_CONFIGURE_ONLY" ]; then
    if [ -n "${ECBUILD_CONFIGURE}" ]; then
      echo "+ module load $1"
      module load $1
    else
      echo " WARNING: Module $1 not loaded (only during configuration)"
    fi
  else
    echo "+ module load $1"
    module load $1
  fi
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

module_load cmake/3.31.6
module_load ninja
module_load prgenv/nvidia
module_load nvidia/24.5
module_load hpcx-openmpi/2.19.0-cuda
module_load fftw
module_load intel-mkl

# Even for nvhpc, we use MKL to ensure bit-reproducibility for CPU builds
export MKL_CBWR=AUTO,STRICT

# Get path this env.sh file is located in:
ARCH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"
export CMAKE_TOOLCHAIN_FILE=${ARCH_DIR}/toolchain.cmake

