---
title: Benchmarking ecTrans
---

@warning
Page under construction.
@endwarning

# Benchmarking ecTrans

A ["benchmark driver" program](https://sites.ecmwf.int/docs/ectrans/sourcefile/ectrans-benchmark.
f90.html) is bundled with ecTrans. This program performs a loop of inverse and
direct spectral transforms over and over a specified number of times and collects timing statistics
to provide an assessment of the overall performance of ecTrans. It is designed to mimic the use of
ecTrans from within the IFS atmospheric model, in which inverse and direct spectral transforms are
carried out on every model timestep. The benchmark program also includes a simple error checking
algorithm for verifying that the transforms are performing with correct numerics. This latter
feature is in fact used for the ecTrans CTest suite.

Here we describe how to write a benchmark suite for ecTrans.

## Installing ecTrans

First follow the [instructions for installing ecTrans](installation.html) on your system. Verify
that the benchmark programs (one for single and double precision) exist in your build's bin
directory. You should see

```bash
ectrans-benchmark-cpu-sp  ectrans-benchmark-cpu-dp
```

Here we assume you have only enabled the `CPU` feature of ecTrans (which is on by default). If you
also enabled the `GPU` feature, you'll also see GPU versions of these two programs. We'll just focus
on CPUs here.

## Using the benchmark program

The benchmark program has many arguments for running ecTrans in different configurations. You can
see the full set by running one of the benchmark programs with the `--help` option:

```
NAME    ectrans-benchmark-cpu-sp

DESCRIPTION
        This program tests ecTrans by transforming fields back and forth between spectral
        space and grid-point space (single-precision version)

USAGE
        ectrans-benchmark-cpu-sp [options]

OPTIONS
    -h, --help          Print this message
    -v                  Run with verbose output
    -t, --truncation T  Run with this triangular spectral truncation (default = 79)
    -g, --grid GRID     Run with this grid. Possible values: O<N>, F<N>
                        If not specified, O<N> is used with N=truncation+1 (cubic relation)
    -n, --niter NITER   Run for this many inverse/direct transform iterations (default = 10)
    -f, --nfld NFLD     Number of scalar fields (default = 1)
    -l, --nlev NLEV     Number of vertical levels (default = 1)
    --vordiv            Also transform vorticity-divergence to wind
    --scders            Compute scalar derivatives (default off)
    --uvders            Compute uv East-West derivatives (default off). Only when also --vordiv is given
    --flt               Run with fast Legendre transforms (default off)
    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)
    --norms             Calculate and print spectral norms of transformed fields
                        The computation of spectral norms will skew overall timings
    --meminfo           Show diagnostic information from FIAT's ec_meminfo subroutine on memory usage, thread-binding etc.
    --nprtrv            Size of V set in spectral decomposition
    --nprtrw            Size of W set in spectral decomposition
    -c, --check VALUE   The multiplier of the machine epsilon used as a tolerance for correctness checking

DEBUGGING
    --dump-values       Output gridpoint fields in unformatted binary file
```

Some of these options (e.g. `-nprtrv`) require a detailed understanding of how fields are
distributed across MPI tasks, so we won't describe them in detail here. The most important arguments
are the following:
- `-t, --truncation T`: this sets the overall resolution of the benchmark. The truncation T refers  
  to the highest zonal and total wavenumber that can be kept in spectral space. By default, a  
  suitable grid point resolution (i.e. a suitable number of latitudes on the octahedral grid) will  
  be chosen for spectral space. This single number then determines the overall problem size of the  
  spectral transform.  The higher this number, the larger the problem size. As of August 2024, the  
  "HRES" (high-resolution, deterministic) forecast of ECMWF uses a spectral truncation of 1279,  
  combined with an octahedral grid of 2560 latitudes, which gives a grid point resolution of  
  approximately 8 km.
- `-n, --niter NITER`: this determines how many iterations to perform in the spectral transform.  
  The more interations you perform, the more reliable the timing statistics you gather. Note that  
  two additional iterations are always performed at the start. This is because (at least for the  
  GPU version of ecTrans) the first two iterations include some initialisation costs which  
  shouldn't be included in any timing statistics.
- `-l, --nlev NLEV`: this determines the number of vertical levels for three-dimensional fields  
  such as U and V wind (or vorticity and divergence). ecTrans can operate on a batch of vertical  
  levels with a single call and this determines the size of this batch (though by default, fields  
  are distributed across MPI tasks on the vertical dimension at some stages in the spectral  
  transform)
- `--vordiv --scders --uvders`: these options enable some auxiliary code paths when calling the  
  inverse transform. `--vordiv` calculates grid point vorticity and divergence, `--scders`  
  calculates derivatives of scalar fields in grid point space, and `--uvders` calculates gradients  
  of the U and V wind in grid point space. For testing code changes, it's good to include these  
  options so as many code paths as possible are verified.
- `--norms`: this option enables error norms, which are printed aggregated over all fields at the  
  end of the benchmark. The errors are computed in spectral space with respect to the initial  
  values of the fields. This is useful to get a good idea that the benchmark is numerically  
  correct.









When inspecting the program, you will notice that it is significantly more complex than, say, the
example program described in our [usage guide](usage.html). This additional complexity comes not
just from the instrumentation code for the timings, but notably also from the infrastructure to
permit transforms of distributed fields. As explained in the [introduction](introduction.html),
ecTrans can operate on fields distributed across MPI tasks, and the dimension across which fields
are split is different for spectral space and grid point space. As such, the benchmark program
includes infrastructure for specifying which elements of the relevant decomposed dimension belong
to which MPI task.


