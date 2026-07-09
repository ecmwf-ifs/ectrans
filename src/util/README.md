This directory contains the internal `ectrans_util` support library used by test code and benchmark-style programs in this repository.

It is not installed and is not intended to be installed as part of the public ectrans package. The build definition in [CMakeLists.txt](./CMakeLists.txt) marks both utility library variants with `NOINSTALL`, and the target is treated as repository-internal tooling rather than a downstream API.

The main consumers are the programs under `src/programs` and the test code under `tests`, where these modules provide common setup, argument parsing, logging, timing, MPI detection, memory handling, and a few convenience helpers.

## Main entry points

### Program setup

- `ectrans_program_init` and `ectrans_program_end` in `ectrans_program_mod.F90`
  Set up and tear down the common runtime environment for repository tools. This includes MPI initialization through MPL when needed, thread-count discovery, logging setup, allocator configuration, GSTATS setup, and DR_HOOK initialization.

- `nthread`, `nproc`, and `myproc` in `ectrans_program_mod.F90`
  Expose the detected OpenMP thread count and the active MPI process layout for utility programs.

### Command-line parsing and user-facing errors

- `ectrans_command_line_parser` in `ectrans_command_line_parser_mod.F90`
  Small stateful command-line parser used by benchmark and test executables. It provides sequential argument iteration, typed value extraction for integers and strings, and help/error dispatch.

- `ectrans_error_parsing_failed` in `ectrans_error_mod.F90`
  Emits concise parsing/configuration errors and aborts in a way that avoids noisy duplicated output from all MPI ranks.

### MPI and process decomposition helpers

- `ectrans_mpi_detect`, `ectrans_mpi_enabled`, `ectrans_mpi_world_size`, `ectrans_mpi_world_rank`, and `ectrans_mpi_world_comm` in `ectrans_mpi_mod.F90`
  Detect whether the current process was launched under MPI, and expose the inferred world communicator metadata without requiring callers to initialize MPI first.

- `ectrans_spectral_decomposition`, `ectrans_gridpoint_decomposition`, and `ectrans_make_spectral_distribution` in `ectrans_decomposition_mod.F90`
  Choose or derive process decompositions for spectral and gridpoint work distribution in tests and benchmark drivers.

### Device and memory helpers

- `ectrans_device_init` and `ectrans_device_is_host` in `ectrans_device_mod.F90`
  Perform minimal accelerator runtime setup where needed and report whether execution is effectively host-side, which is used for choices such as memory pinning.

- `allocator` in `ectrans_memory_mod.F90`
  Repository-local allocation wrapper with generic `allocate` and `deallocate` bindings for rank-1 to rank-4 real arrays in single and double precision. It also provides `set_pinning`, `set_logging`, and `set_logging_output_unit` to control the underlying C allocator in `ectrans_memory.c`.

- `ectrans_memory.c`
  Implements the low-level allocation backend, including optional host-memory pinning for CUDA or HIP builds and optional allocation logging.

### Logging, timing, and statistics

- `ectrans_log_init` and `ectrans_log_end` in `ectrans_log_mod.F90`
  Set the output units and suppress most non-root-rank output for utility programs.

- `ectrans_timer` and `ectrans_timings` in `ectrans_timer_mod.F90`
  Provide simple wall-clock timing plus local and MPI-global summary statistics such as min, max, average, and median.

- `ectrans_gstats_init`, `ectrans_gstats_end`, `ectrans_gstats_new_region`, `ectrans_gstats_labels`, `ectrans_gstats_enable`, and `ectrans_gstats_configuration` in `ectrans_gstats_mod.F90`
  Integrate utility programs with GSTATS region registration and configuration.

### Miscellaneous helpers

- `cubic_octahedral_gaussian_grid` and `parse_gaussian_grid` in `ectrans_grids_mod.F90`
  Convert compact Gaussian-grid identifiers such as `F128` or `O128` into the grid metadata used by the benchmark drivers.

- `ECTRANS_ENABLE_FPE` in `ectrans_fpe_trapping.F90`
  Enables floating-point exception trapping through the companion C++ helper when the platform supports it.

## Scope

This utility layer is intentionally pragmatic and local to the repository. It exists to keep test and benchmark code readable and to avoid duplicating setup logic across internal executables. It should be treated as internal support code rather than a stable external interface.