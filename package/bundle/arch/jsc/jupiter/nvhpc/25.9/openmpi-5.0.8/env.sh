# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $*"
  module load $*
}
module_unload() {
  echo "+ module unload $*"
  module unload $*
}
module_purge() {
  echo "+ module purge"
  module purge
}

# Unload all modules to be certain
[[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]] && module_purge

# Load modules
module_load Stages/2026
module_load OpenSSL/3
module_load CUDA/13
module_load StdEnv/2026
module_load nvidia-compilers/25.9-CUDA-13
module_load git/2.50.1
module_load CMake/4.0.3
module_load OpenMPI/5.0.8
module_load FFTW/3.3.10
[[ ${ECTRANS_GPU_AWARE_MPI:-unset} == "1" ]] && module_load MPI-settings/CUDA

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

