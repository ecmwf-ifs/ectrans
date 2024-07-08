# ecTrans API

## `SETUP_TRANS0`

### Signature

```
SUBROUTINE SETUP_TRANS0(KOUT, KERR, KPRINTLEV, KMAX_RESOL, KPROMATR, &
  &                     KPRGPNS, KPRGPEW, KPRTRW, KCOMBFLEN, &
  &                     LDMPOFF, LDSYNC_TRANS, KTRANS_SYNC_LEVEL, &
  &                     LDEQ_REGIONS, K_REGIONS_NS, K_REGIONS_EW, K_REGIONS, &
  &                     PRAD, LDALLOPERM, KOPT_MEMORY_TR)
```

### Purpose

This subroutine initialises all resolution-agnostic aspects of the ecTrans library. It must be
called before any other ecTrans subroutines.

### Arguments

- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOUT`  
  Unit number for listing output.  
  *Default*: 6
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KERR`  
  Unit number for error messages.  
  *Default*: 0
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRINTLEV`  
  Verbosity of standard output (0: output disabled, 1: normal, 2: debug).  
  *Default*: 0
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMAX_RESOL`  
  Maximum number of different resolutions for which spectral transforms will be performed.
  *Default*: 1
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMATR`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRGPNS`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRGPEW`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRTRW`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KCOMBFLEN`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDMPOFF`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDSYNC_TRANS`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KTRANS_SYNC_LEVEL`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDEQ_REGIONS`
- `LOGICAL, OPTIONAL, INTENT(IN) :: LDALLOPERM`
- `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PRAD`
- `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOPT_MEMORY_TR`
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