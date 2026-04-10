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

toload=""
module_load() {
  echo "+ module load $*"
  module load $*
}
module_unload() {
  echo "+ module unload $*"
  module unload $*
}

# Unload all modules to be certain
module --force purge

# Load modules
module_load PrgEnv-cray/8.6.0
module_load LUMI/25.03
module_load cce/19.0.0
module_load craype-accel-amd-gfx90a
module_load rocm/6.3.4
module_load cray-fftw
module_load buildtools
module_load craype
module_load cray-mpich/8.1.32

# Get path this env.sh file is located in:
ARCH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"
export CMAKE_TOOLCHAIN_FILE=${ARCH_DIR}/toolchain.cmake

