---
title: Usage
---

# Usage

Here we describe how to construct from scratch a minimal working program that uses ecTrans.

## Getting started

We begin with a standard skeleton for a Fortran program:

```f90
PROGRAM ECTRANS_DEMONSTRATION

IMPLICIT NONE

END PROGRAM ECTRANS_DEMONSTRATION
```

Save this in a file called ectrans_demonstrator.F90.

In order to use ecTrans within a Fortran program, you must include the header files corresponding to
the subroutines you intend to call. Here we will call `SETUP_TRANS0`, `SETUP_TRANS`, `INV_TRANS`,
`DIR_TRANS`, and `TRANS_INQ`. These contain `INTERFACE` blocks, so we place them at the top of the
program source file between the `IMPLICIT NONE` statement and the main body of the source code:

```f90
#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
```

Now let's start calling ecTrans subroutines straight away. Add a call to `SETUP_TRANS0`:

```f90
CALL SETUP_TRANS0
```

All arguments here are optional, so let's just leave them out for now.

We can now try building our program just to see if it at least runs.

## Building ecTrans with CMake

CMake might seem overkill for a simple program, but we find it is still easier than manually
crafting a Makefile. Since ecTrans is built with CMake, much of the hard work has already been done.

## Building ecTrans without CMake

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