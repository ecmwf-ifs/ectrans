! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FLT

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT, JPRD
USE BUTTERFLY_ALG_MOD,ONLY : BUTTERFLY_STRUCT
USE SEEFMM_MIX
IMPLICIT NONE

SAVE


TYPE FLT_TYPE
INTEGER(KIND=JPIM) :: NSPOLEGL
INTEGER(KIND=JPIM) :: NDGNH
INTEGER(KIND=JPIM) :: INS2
INTEGER(KIND=JPIM) :: INA2
REAL(KIND=JPRBT) ,POINTER :: RPNMS(:,:) ! Legendre polynomials
REAL(KIND=JPRBT) ,POINTER :: RPNMA(:,:) ! Legendre polynomials
REAL(KIND=JPRD) ,POINTER :: RPNMDS(:,:) ! Legendre polynomials
REAL(KIND=JPRD) ,POINTER :: RPNMDA(:,:) ! Legendre polynomials
REAL(KIND=JPRBT) :: RCS
REAL(KIND=JPRBT) :: RCA
!REAL(KIND=JPRBT) ,POINTER :: RPNMCDO(:,:) ! Legendre polynomials for C-D formula at orig roots
!REAL(KIND=JPRBT) ,POINTER :: RPNMCDD(:,:) ! Legendre polynomials for C-D formula at dual roots
REAL(KIND=JPRBT) ,POINTER :: RPNMWI(:,:) ! special weights
REAL(KIND=JPRBT) ,POINTER :: RPNMWO(:,:) ! special weights
INTEGER(KIND=JPIM) :: ISLD ! starting latitude dual

! Butterfly

INTEGER(KIND=JPIM) :: MAXCOLS
TYPE(BUTTERFLY_STRUCT) :: YBUT_STRUCT_S,YBUT_STRUCT_A

END TYPE FLT_TYPE

TYPE FLT_TYPE_WRAP
TYPE(FLT_TYPE),ALLOCATABLE :: FA(:)
LOGICAL :: LDLL
LOGICAL :: LSHIFTLL
LOGICAL :: LUSEFLT
LOGICAL :: LUSE_BELUSOV
LOGICAL :: LKEEPRPNM
LOGICAL :: LSOUTHPNM ! .TRUE. to compute Legendre polynomials on southern hemisphere
INTEGER(KIND=JPIM) :: IMLOC
INTEGER(KIND=JPIM) :: ITHRESHOLD
INTEGER(KIND=JPIM) :: NDGNHD ! dual set dimension
INTEGER(KIND=JPIM) :: NDLON  ! dual number of longitudes
INTEGER(KIND=JPIM) :: NDGL   ! dual number of latitudes
LOGICAL :: LSYM
TYPE(FMM_TYPE),POINTER :: FMM_INTI ! FMM interpolation

END TYPE FLT_TYPE_WRAP

TYPE(FLT_TYPE_WRAP),ALLOCATABLE,TARGET :: FLT_RESOL(:)
TYPE(FLT_TYPE_WRAP),POINTER     :: S


END MODULE TPM_FLT







