! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KDLON,KLOEN,LDSPLIT,PSTRET,&
&KTMAX,KRESOL,PWEIGHT,LDGRIDONLY,LDUSERPNM,LDKEEPRPNM,LDUSEFLT,&
&LDSPSETUPONLY,LDPNMONLY,LDUSEFFTW,&
&LDLL,LDSHIFTLL,CDIO_LEGPOL,CDLEGPOLFNAME,KLEGPOLPTR,KLEGPOLPTR_LEN)

! begin_doc_block
! ## `SETUP_TRANS`
!
! ### Signature
!
! ```f90
! SUBROUTINE SETUP_TRANS(KSMAX, KDGL, KDLON, KLOEN, LDSPLIT, PSTRET, KTMAX, KRESOL, &
!   &                    PWEIGHT, LDGRIDONLY, LDUSERPNM, LDKEEPRPNM, LDUSEFLT, &
!   &                    LDSPSETUPONLY, LDPNMONLY, LDUSEFFTW, LDLL, LDSHIFTLL, &
!   &                    CDIO_LEGPOL, CDLEGPOLFNAME, KLEGPOLPTR, KLEGPOLPTR_LEN)
! ```
!
! ### Purpose
!
! This subroutine establishes the environment for performing a spectral transform. Each call to this
! subroutine creates work arrays for a certain resolution. One can do this for `NMAX_RESOL` different
! resolutions (internal parameter corresponding to `KMAX_RESOL` argument of `SETUP_TRANS0`). One
! must call `SETUP_TRANS0` before calling this subroutine for the first time.
!
! ### `INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KSMAX`  
!   Spectral truncation to perform the transform up to.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KDGL`  
!   The number of Gaussian latitudes in grid point space.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KDLON`  
!   Maximum number of points on any latitude (usually the latitude nearest to the Equator).  
!   If not provided, it will take the value of `2 * KDGL`, unless `LDLL` is `.TRUE.` in which case  
!   it will take the value of `2 * KDGL + 2`, unless `KLOEN` is provided in which case it will take  
!   the value of `MAXVAL(KLOEN)`.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KLOEN(:)`  
!   An array giving the number of points on each Gaussian latitude (dimension)  
!   If this is not provided it will be assumed that a full Gaussian grid is used in grid point space,  
!   for which every latitude will have `KDLON` points.
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDSPLIT`  
!   Split latitude in grid point space.  
!   *Default*: `.FALSE.`
! - `REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: PSTRET`  
!   Stretching factor for when the Legendre polynomials are computed on the stretched sphere.  
!   *Default*: `1.0`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KTMAX`  
!   Spectral truncation to be applied for tendencies.  
!   *Default*: `KSMAX`
! - `REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: PWEIGHT(:)`  
!   The weight per grid-point (for a weighted distribution).
!   If this argument is not provided, a weighted distribution will not be used.
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDGRIDONLY`  
!   Only provide grid point space results.  
!   *Default*: `.FALSE.`  
!   **Potentially deprecatable.**
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDUSEFLT`  
!   Use the fast Legendre transform algorithm.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDUSERPNM`  
!   Use the Belusov algorithm to compute the Legendre polynomials.  
!   *Default*: `.TRUE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDKEEPRPNM`  
!   Store Legendre Polynomials instead of recomputing (only applicable when using the fast Legendre
!   transform, otherwise these are always stored).  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDSPSETUPONLY`  
!   Only create distributed spectral setup.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDPNMONLY`  
!   Compute the Legendre polynomials only, not the FFTs.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDLL`  
!   Setup second set of input/output latitudes.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDSHIFTLL`  
!   Shift output lon/lat data by 0.5\*dx and 0.5\*dy.  
!   *Default*: `.FALSE.`
! - `CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDIO_LEGPOL`  
!   String config argument to determine several options relating to I/O operations relevant to the
!   Legendre polynomials. If `'readf'` is passed, Legendre polynomials will be read from the file
!   given to `CDLEGPOLFNAME`. If `'writef'` is passed, Legendre polynomials will be written to the
!   file given to `CDLEGPOLFNAME`. If `'membuf'` is passed, Legendre polynomials will be read from
!   shared memory.
! - `CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDLEGPOLFNAME`  
!   Filename from which to read or to which to write the Legendre polynomials. See `CDIO_LEGPOL`.
! - `TYPE(C_PTR), OPTIONAL, INTENT(IN) :: KLEGPOLPTR`  
!   Pointer to memory segment containing Legendre polynomials. See `CDIO_LEGPOL`.
! - `INTEGER(C_SIZE_T), OPTIONAL, INTENT(IN) :: KLEGPOLPTR_LEN`  
!   Length of memory segment containing Legendre polynomials. See `CDIO_LEGPOL`.
!
! ### `OPTIONAL, INTENT(OUT)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KRESOL`  
!   A unique integer identifier to act as a handle for the work arrays corresponding to the provided
!   transform resolution.
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRD
    USE, INTRINSIC :: ISO_C_BINDING, ONLY:  C_PTR, C_INT,C_ASSOCIATED,C_SIZE_T


IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) ,INTENT(IN) :: KSMAX,KDGL
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KDLON
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KTMAX
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT):: KRESOL
REAL(KIND=JPRD)    ,OPTIONAL,INTENT(IN) :: PWEIGHT(:)
REAL(KIND=JPRD)    ,OPTIONAL,INTENT(IN) :: PSTRET
LOGICAL   ,OPTIONAL,INTENT(IN):: LDGRIDONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFLT
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSERPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDKEEPRPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDPNMONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSPSETUPONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFTW
LOGICAL   ,OPTIONAL,INTENT(IN):: LDLL
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSHIFTLL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDIO_LEGPOL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDLEGPOLFNAME
TYPE(C_PTR) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR
INTEGER(C_SIZE_T) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR_LEN


END SUBROUTINE SETUP_TRANS


END INTERFACE
