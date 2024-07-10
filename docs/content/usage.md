---
title: Usage
---

# Usage

Here we describe how to construct from scratch a minimal working program that uses ecTrans. The
instructions here assume you have already installed FIAT and ecTrans according to
[our instructions](installation.html). Therefore, in the current directory you should see fiat/build
and ectrans/build.

## Getting started

We begin with the simplest Fortran program that calls an ecTrans subroutine:

```f90
PROGRAM ECTRANS_DEMONSTRATOR

IMPLICIT NONE

#include "setup_trans0.h"

CALL SETUP_TRANS0(LDMPOFF=.TRUE.)

END PROGRAM ECTRANS_DEMONSTRATOR
```

Save this in a file called ectrans_demonstrator.F90.

In order to use ecTrans within a Fortran program, you must include the header files corresponding to
the subroutines you intend to call. Only some of the ecTrans subroutines are callable from an
external program and we refer to these as the "external" subroutines. They are documented on the
[API page](api.html). Each external subroutine exists in its own source file with its own header
file. The header files contains an `INTERFACE` block so the compiler knows what the signature of the
subroutine is. The `#include` statements must be placed between the `IMPLICIT NONE` and the main
body of the source code.

All arguments to `SETUP_TRANS0` are optional, so for now let's just provide one: `LDMPOFF` which
turns off MPI. This means we don't have to build with MPI support for this simple test.

We can now try building our program just to see if it at least runs.

## Building our program with CMake

CMake might seem overkill for a simple program, but we find it is still easier than manually
crafting a Makefile. Since ecTrans is built with CMake, much of the hard work has already been done.

Here is a minimal toplevel CMakeLists.txt for this guide:

```cmake
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ectrans_demonstrator LANGUAGES Fortran )

find_package( fiat REQUIRED )
find_package( ectrans REQUIRED )

add_executable( ectrans_demonstrator ectrans_demonstrator.F90 )

target_include_directories( ectrans_demonstrator 
                            PRIVATE ${CMAKE_SOURCE_DIR}/ectrans/src/trans/include/ectrans )
target_include_directories( ectrans_demonstrator PRIVATE ${fiat_ROOT}/module/parkind_sp )
target_link_libraries( ectrans_demonstrator PRIVATE trans_sp )
```

ecTrans supports single- and double-precision arithmetic. Currently you must decide which one you
want to use at compile time, as there is a separate library for each kind. Here we link against the
single-precision library which means all `REAL` variables passed to `INV_TRANS` and `DIR_TRANS` must
be single precision. 

With this file you can follow the standard procedure for configuring and building with CMake:

```bash
mkdir build && cd build
cmake .. -Dfiat_ROOT=../fiat/build -Dectrans_ROOT=../ectrans/build
make
```

Now try running the program with `./ectrans_demonstrator`. You should see output like the following:

```
ecTrans at version: 1.3.2
commit: 023f5e320846f7a2ef538166366f68ce2046efa1
```

## Building our program without CMake

It is still possibly to build ecTrans without CMake by manually writing a Makefile. As with any
library, we have to ensure all header files are on the include path and the ecTrans library is
included at the link time of the program. That looks like this:

```Makefile
FC = <your compiler>

ectrans_demonstrator: ectrans_demonstrator.o
	$(FC) $^ -o ectrans_demonstrator

ectrans_demonstrator.o: ectrans_demonstrator.F90

%.o: %.f90
	$(FC) -c $< -o $@
```

## Setting up ecTrans

Let's make our program more interesting. Let's use ecTrans to visualise some of the spherical
harmonics! We will use a spectral truncation of 79 which corresponds to a full Gaussian grid with 
160 latitudes and 320 points per latitude. First we need to add a call to `SETUP_TRANS` so that
everything is initialised. We can rely mostly on the default parameters, but we always need to
provide `KSMAX`, the spectral truncation, and `KDGL`, the number of Gaussian latitudes pole to pole.

```f90
CALL SETUP_TRANS(KSMAX=79, KDGL=160)
```

Remember to also add "setup_trans.h" to the list of included files at the top.

Now we need to allocate arrays to store our field in spectral space and grid point space. As
mentioned earlier, all variables passed to `INV_TRANS` and `DIR_TRANS` must be single precision in
this test. We recommend importing the variable `JPRM` from the module `PARKIND1` (part of FIAT) and
using this as the `KIND` when declaring all `REAL` variables. Those variables will then be 4 bytes
wide as expected by ecTrans. While we're here, let's also import the standard integer kind `JPIM` as
well. So in other words, add this before your `IMPLICIT NONE`:

```f90
USE PARKIND1, ONLY: JPRM, JPIM
```

Now we need to create an array to store spectral coefficients. This will be an `ALLOCATABLE`, but
what size should we allocate for it? We can use the `TRANS_INQ` subroutine to inquire about how many
spectral coefficients are required to represent a field in spectral space with a truncation of 79:

```f90
CALL TRANS_INQ(KSPEC2=NSPEC2)
```

Make sure you declare `NSPEC2` as `INTEGER(KIND=JPIM)`. This variable should take a value of 6480,
which is \( 2\times(N+1)(N+2)/2 \), where \(N\) is the truncation (79), the first factor of two is
because spectral coefficients have an imaginary and real part and the division by two is because
ecTrans assumes all grid-point space variables are real. That means the negative zonal modes in
spectral space must be the complex conjugates of the positive ones, so we don't need to store them
explicitly. This reduces the number of coefficients by half.

Now we're ready to perform a spectral transform. 