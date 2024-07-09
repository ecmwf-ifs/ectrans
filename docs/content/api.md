---
title: API
---

# ecTrans API

## General notes

NOTE: ecTrans is a legacy code with an accumulated 30 years of history. Over this time certain
features enabled through optional arguments will have fallen out of use. We are currently reviewing
all options to identify those that can be safely deleted, but this takes time. In the mean time, we
have tagged below all options we deem to be "potentially deprecatable".

### Variable names

ecTrans _in principle_ follows the coding standard and conventions outlined in the [IFS
Documentation - Part VI: Technical and Computational Procedures](https://www.ecmwf.int/en/elibrary/
81372-ifs-documentation-cy48r1-part-vi-technical-and-computational-procedures) section 1.5.
Following these standards, all variable names must begin with a one- or two-character prefix
denoting their scope (module level, dummy argument, local variables, loop index, or parameter) and
type. These are outlined in Table 1.2. Dummy variables have the following prefixes:

- `K` - integer
- `P` - real (single or double precision)
- `LD` - logical
- `CD` - character
- `YD` - derived type

### `KIND` parameters

As with the IFS, integer and real variables in ecTrans always have an explicit `KIND` specification.
These are defined in the [`PARKIND1` module](https://github.com/ecmwf-ifs/fiat/blob/main/src/parkind
/parkind1.F90) which is part of the Fiat library (a dependency of ecTrans). To understand the
subroutines described here, only two must be considered:

- `INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)` (i.e. 4-byte integer)
- `INTEGER, PARAMETER :: JPRD = SELECTED_REAL_KIND(13,300)` (i.e. 8-byte float)

## `SETUP_TRANS0`

### Signature

```f90
SUBROUTINE SETUP_TRANS0(KOUT, KERR, KPRINTLEV, KMAX_RESOL, KPROMATR, KPRGPNS, KPRGPEW, &
  &                     KPRTRW, LDMPOFF, LDSYNC_TRANS, KTRANS_SYNC_LEVEL, LDEQ_REGIONS, &
  &                     K_REGIONS_NS, K_REGIONS_EW, K_REGIONS, PRAD, LDALLOPERM, &
  &                     KOPT_MEMORY_TR)
```

### Purpose

This subroutine initialises all resolution-agnostic aspects of the ecTrans library. It must be
called before any other ecTrans subroutines.

### `OPTIONAL, INTENT(IN)` arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOUT`  
  Fortran unit number for listing output.  
  *Default*: `6` (i.e. STDOUT)
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KERR`  
  Unit number for error messages.  
  *Default*: `0` (i.e. STDERR)
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRINTLEV`  
  Verbosity of standard output (0: output disabled, 1: normal, 2: debug).  
  *Default*: `0`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMAX_RESOL`  
  Maximum number of different resolutions for which spectral transforms will be performed. This is
  required in advance so ecTrans knows how to correctly size the work arrays.  
  *Default*: `1`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMATR`  
  Batch size for splitting of vertical field set. This allows one to transform one batch of vertical  
  levels at a time rather than transforming all in one go. A value of 0 disables this feature.  
  *Default*: `0`  
  **Potentially deprecatable**
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRGPNS`  
  Splitting level in North-South direction in grid-point space.  
  *Default*: `1`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRGPEW`  
  Splitting level in East-West direction in grid-point space.  
  *Default*: `1`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRTRW`  
  Splitting level in wave direction in spectral space.  
  *Default*: `1`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDMPOFF`  
  Switch off message passing.  
  *Default*: `.FALSE.`  
  **Potentially deprecatable**
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDSYNC_TRANS`  
  Switch to activate barriers in M->L and L->M transposition routines.  
  *Default*: `.FALSE.`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KTRANS_SYNC_LEVEL`  
  Control parameter for synchronization and blocking in communication routines.  
  *Default*: `0`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDEQ_REGIONS`  
  Enable EQ_REGIONS partitioning mode.  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDALLOPERM`  
  Allocate certain arrays permanently.  
  *Default*: `.FALSE.`
- `REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: PRAD`  
  Radius of the planet (Earth) in metres.  
  *Default*: `6371229.0`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOPT_MEMORY_TR`  
  Memory strategy for gridpoint transpositions (0: heap, 1: stack).  
  *Default*: `0`

### `OPTIONAL, INTENT(OUT)` arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS(:)`  
  Number of regions returned by EQ_REGIONS algorithm.
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS_NS`  
  Maximum number of North-South partitions, as determined by EQ_REGIONS algorithm.
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS_EW`  
  Maximum number of East-West partitions, as determined by EQ_REGIONS algorithm.

## `SETUP_TRANS`

### Signature

```f90
SUBROUTINE SETUP_TRANS(KSMAX, KDGL, KDLON, KLOEN, LDSPLIT, PSTRET, KTMAX, KRESOL, &
  &                    PWEIGHT, LDGRIDONLY, LDUSERPNM, LDKEEPRPNM, LDUSEFLT, &
  &                    LDSPSETUPONLY, LDPNMONLY, LDUSEFFTW, LDLL, LDSHIFTLL, &
  &                    CDIO_LEGPOL, CDLEGPOLFNAME, KLEGPOLPTR, KLEGPOLPTR_LEN)
```

### Purpose

This subroutine establishes the environment for performing a spectral transform. Each call to this
subroutine creates work arrays for a certain resolution. One can do this for `NMAX_RESOL` different
resolutions (internal parameter corresponding to `KMAX_RESOL` argument of `SETUP_TRANS0`). One
must call `SETUP_TRANS0` before calling this subroutine for the first time.

### `INTENT(IN)` arguments

- `INTEGER(KIND=JPIM), INTENT(IN) :: KSMAX`  
  Spectral truncation to perform the transform up to.
- `INTEGER(KIND=JPIM), INTENT(IN) :: KDGL`  
  The number of Gaussian latitudes in grid point space.

### `OPTIONAL, INTENT(IN)` arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KDLON`  
  Maximum number of points on any latitude (usually the latitude nearest to the Equator).  
  If not provided, it will take the value of `2 * KDGL`, unless `LDLL` is `.TRUE.` in which case  
  it will take the value of `2 * KDGL + 2`, unless `KLOEN` is provided in which case it will take  
  the value of `MAXVAL(KLOEN)`.
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KLOEN(:)`  
  An array giving the number of points on each Gaussian latitude (dimension)  
  If this is not provided it will be assumed that a full Gaussian grid is used in grid point space,  
  for which every latitude will have `KDLON` points.
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDSPLIT`  
  Split latitude in grid point space.  
  *Default*: `.FALSE.`
- `REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: PSTRET`  
  Stretching factor for when the Legendre polynomials are computed on the stretched sphere.  
  *Default*: `1.0`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KTMAX`  
  Spectral truncation to be applied for tendencies.  
  *Default*: `KSMAX`
- `REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: PWEIGHT(:)`  
  The weight per grid-point (for a weighted distribution).
  If this argument is not provided, a weighted distribution will not be used.
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDGRIDONLY`  
  Only provide grid point space results.  
  *Default*: `.FALSE.`  
  **Potentially deprecatable**
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDUSEFLT`  
  Use the fast Legendre transform algorithm.  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDUSERPNM`  
  Use the Belusov algorithm to compute the Legendre polynomials.  
  *Default*: `.TRUE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDKEEPRPNM`  
  Store Legendre Polynomials instead of recomputing (only applicable when using the fast Legendre
  transform, otherwise these are always stored).  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDSPSETUPONLY`  
  Only create distributed spectral setup.  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDPNMONLY`  
  Compute the Legendre polynomials only, not the FFTs.  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDLL`  
  Setup second set of input/output latitudes.  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDSHIFTLL`  
  Shift output lon/lat data by 0.5\*dx and 0.5\*dy.  
  *Default*: `.FALSE.`
- `CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDIO_LEGPOL`  
  String config argument to determine several options relating to I/O operations relevant to the
  Legendre polynomials. If `'readf'` is passed, Legendre polynomials will be read from the file
  given to `CDLEGPOLFNAME`. If `'writef'` is passed, Legendre polynomials will be written to the
  file given to `CDLEGPOLFNAME`. If `'membuf'` is passed, Legendre polynomials will be read from
  shared memory.
- `CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDLEGPOLFNAME`  
  Filename from which to read or to which to write the Legendre polynomials. See `CDIO_LEGPOL`.
- `TYPE(C_PTR), OPTIONAL, INTENT(IN) :: KLEGPOLPTR`  
  Pointer to memory segment containing Legendre polynomials. See `CDIO_LEGPOL`.
- `INTEGER(C_SIZE_T), OPTIONAL, INTENT(IN) :: KLEGPOLPTR_LEN`  
  Length of memory segment containing Legendre polynomials. See `CDIO_LEGPOL`.

### `OPTIONAL, INTENT(OUT)` arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KRESOL`  
  A unique integer identifier to act as a handle for the work arrays corresponding to the provided
  transform resolution.

## `DIR_TRANS`

### Signature

```f90
SUBROUTINE DIR_TRANS(PSPVOR, PSPDIV, PSPSCALAR, PSPSC3A, PSPSC3B, PSPSC2, LDLATLON, &
  &                  KPROMA, KVSETUV, KVSETSC, KRESOL, KVSETSC3A, KVSETSC3B, KVSETSC2, &
  &                  PGP, PGPUV, PGP3A, PGP3B, PGP2)
```

### Purpose

### `OPTIONAL, INTENT(IN)` arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMA`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETUV(:)`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC(:)`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC3A(:)`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC3B(:)`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC2(:)`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDLATLON`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP(:,:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGPUV(:,:,:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP3A(:,:,:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP3B(:,:,:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP2(:,:,:)`

### `OPTIONAL, INTENT(OUT)` arguments

- `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)`

## `INV_TRANS`

### Signature

### Purpose

## `TRANS_RELEASE`

### Signature

### Purpose

## `TRANS_END`

### Signature

### Purpose

## `TRANS_INQ`

### Signature

### Purpose

## `SPECNORM`

### Signature

### Purpose

## `DIR_TRANSAD`

### Signature

### Purpose

## `INV_TRANSAD`

### Signature

### Purpose

## `DIST_GRID`

### Signature

### Purpose

## `DIST_GRID32`

### Signature

### Purpose

## `GATH_GRID`

### Signature

### Purpose

## `GATH_GRID_32`

### Signature

### Purpose

## `GATH_SPEC`

### Signature

### Purpose

## `GET_CURRENT`

### Signature

### Purpose

## `GPNORM_TRANS`

### Signature

### Purpose

## `INI_SPEC_DIST`

### Signature

### Purpose

## `SUGAWC`

### Signature

### Purpose

## `TRANS_PNM`

### Signature

### Purpose

## `VORDIV_TO_UV`

### Signature

### Purpose
