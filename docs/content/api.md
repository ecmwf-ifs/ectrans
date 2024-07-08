# ecTrans API

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

NOTE: ecTrans is a legacy code with an accumulated 30 years of history. Over this time certain
features enabled through optional arguments will have fallen out of use. We are currently reviewing
all options to identify those that can be safely deleted, but this takes time. In the mean time, we
have tagged below all options we deem to be "potentially deprecatable".

## `SETUP_TRANS0`

### Signature

```
SUBROUTINE SETUP_TRANS0(KOUT, KERR, KPRINTLEV, KMAX_RESOL, KPROMATR, &
  &                     KPRGPNS, KPRGPEW, KPRTRW, &
  &                     LDMPOFF, LDSYNC_TRANS, KTRANS_SYNC_LEVEL, &
  &                     LDEQ_REGIONS, K_REGIONS_NS, K_REGIONS_EW, K_REGIONS, &
  &                     PRAD, LDALLOPERM, KOPT_MEMORY_TR)
```

### Purpose

This subroutine initialises all resolution-agnostic aspects of the ecTrans library. It must be
called before any other ecTrans subroutines.

### `INTENT(IN)` arguments

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
  Enable EQ_REGIONS partioning mode.  
  *Default*: `.FALSE.`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDALLOPERM`  
  Allocate certain arrays permanently.  
  *Default*: `.FALSE.`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PRAD`  
  Radius of the planet (Earth) in metres.  
  *Default*: `6371229.0`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOPT_MEMORY_TR`  
  Memory strategy for gridpoint transpositions (0: heap, 1: stack).  
  *Default*: `0`

### `INTENT(OUT)` arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS(:)`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS_NS`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS_EW`

## `SETUP_TRANS`

## `DIR_TRANS`

## `INV_TRANS`

## `TRANS_RELEASE`

## `TRANS_END`

## `TRANS_INQ`

## `SPECNORM`

## `DIR_TRANSAD`

## `INV_TRANSAD`

## `DIST_GRID`

## `DIST_GRID32`

## `GATH_GRID`

## `GATH_GRID_32`

## `GATH_SPEC`

## `GET_CURRENT`

## `GPNORM_TRANS`

## `INI_SPEC_DIST`

## `SUGAWC`

## `TRANS_PNM`

## `VORDIV_TO_UV`