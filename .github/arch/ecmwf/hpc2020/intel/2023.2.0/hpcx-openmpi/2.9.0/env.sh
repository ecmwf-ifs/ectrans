# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Unload all modules to be certain
module purge

# Load modules
module load prgenv/intel
module load intel/2023.2.0
module load hpcx-openmpi/2.9.0
module load intel-mkl/19.0.5
module load fftw/3.3.10
module load cmake/3.25.2
module load ninja/1.10.0

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Get path this env.sh file is located in:
ARCH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"
export CMAKE_TOOLCHAIN_FILE=${ARCH_DIR}/toolchain.cmake

