# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  module load $1
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

module purge

module_load llvm/EPI-development
module_load openmpi/riscv/4.1.6_llvm1.0
module_load hdf5/ubuntu/1.14.6_llvmEPIdev
module_load cmake/3.28.1
module_load openBLAS/ubuntu/0.3.30_vlen256_llvmEPI1.0
module_load fftw/ubuntu/3.3.10_ompi4.1.6_llvm1.0

export FC=flang
export CC=clang
export CXX=clang++
set -x
ulimit -s unlimited

export PMIX_MCA_pcompress_base_silence_warning=1

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"

