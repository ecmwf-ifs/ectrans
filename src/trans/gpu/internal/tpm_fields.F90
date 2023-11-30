! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

#include "renames.inc"
MODULE TPM_FIELDS

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE ISO_C_BINDING

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRD)   ,ALLOCATABLE :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRD)   ,ALLOCATABLE :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRB)   ,ALLOCATABLE :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRB)   ,ALLOCATABLE :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRB)   ,ALLOCATABLE :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRB)   ,ALLOCATABLE :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRB)   ,ALLOCATABLE :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRB)   ,ALLOCATABLE :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTN(:)   ! R%NTMAX+2-JN

REAL(KIND=JPRB)   ,ALLOCATABLE :: RMU2(:)   ! sin(theta) for dual input/output latitudes
REAL(KIND=JPRB)   ,ALLOCATABLE :: RACTHE2(:)! 1./SQRT(R1MU2), 1/(cos(theta)) dual input/output latitudes
END TYPE FIELDS_TYPE

!flat copies of the above
REAL(KIND=JPRB) ,ALLOCATABLE :: F_RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRB) ,ALLOCATABLE :: F_RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRB) ,ALLOCATABLE :: F_RLAPIN(:) ! eigen-values of the inverse Laplace operator
REAL(KIND=JPRB) ,ALLOCATABLE :: F_RACTHE(:) ! eigen-values of the inverse Laplace operator

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

! scratch arrays for ltinv and ltdir and associated dimension variables

REAL(KIND=JPRB),ALLOCATABLE :: ZAA(:,:,:)  !! JPRL for 1/2
REAL(KIND=JPRB),ALLOCATABLE :: ZAS(:,:,:)  !! JPRL for 1/2

REAL(KIND=JPRB), POINTER :: IZBA(:,:,:)    !! JPRL for 1/2
!!origSam REAL(KIND=JPRB),ALLOCATABLE,TARGET  :: IZBS(:,:,:) !! JPRL for 1/2
REAL(KIND=JPRB),ALLOCATABLE :: IZBS(:)  !! from working RAPS
REAL(KIND=JPRB),ALLOCATABLE :: IZCA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: IZCS(:,:,:)
!REAL(KIND=JPRB),ALLOCATABLE :: IZCAT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: IZCST(:)

REAL(KIND=JPRB),ALLOCATABLE :: DZBA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZBS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZBAT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZBST(:) !! JPRL for 1/2
REAL(KIND=JPRB),ALLOCATABLE :: DZCA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZCS(:,:,:)
!REAL(KIND=JPRB),POINTER :: DZCAT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZCAT(:)
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: DZCST(:)

! Arrays used for rescaling to allow half-precision Legende transforms
REAL(KIND=JPRB), ALLOCATABLE :: ZAMAX(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSMAX(:,:)

! for m=0 in ledir_mod:
REAL(KIND=JPRD),ALLOCATABLE :: ZAA0(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAS0(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: DZBST0(:)
REAL(KIND=JPRD),ALLOCATABLE :: DZCAT0(:)
REAL(KIND=JPRD),ALLOCATABLE :: DZCST0(:)
INTEGER(KIND=JPIM) :: KMLOC0

INTEGER(KIND=JPIM) :: LDZAA
INTEGER(KIND=JPIM) :: LDZAS
INTEGER(KIND=JPIM) :: TDZAA
INTEGER(KIND=JPIM) :: TDZAS

INTEGER(KIND=JPIM) :: ILDZBA
INTEGER(KIND=JPIM) :: ILDZBS
INTEGER(KIND=JPIM) :: ILDZCA
INTEGER(KIND=JPIM) :: ILDZCS



INTEGER(KIND=JPIM) :: DLDZBA
INTEGER(KIND=JPIM) :: DLDZBS
INTEGER(KIND=JPIM) :: DLDZCA
INTEGER(KIND=JPIM) :: DLDZCS

! enable calling setup_trans with a different set of fields than inv_trans and dir_trans:
! IF_FS_INV0: size used for the allocation in setup_trans
! IF_FS_INV: size used in inv_trans and dir_Trans, needs to be <= IF_FS_INV0 
INTEGER(KIND=JPIM) :: IF_FS_INV, IF_FS_INV0
INTEGER(KIND=JPIM) :: IF_FS_DIR, IF_FS_DIR0
INTEGER(KIND=JPIM) :: NFLEV, NFLEV0
INTEGER(KIND=JPIM) :: ITDZBA, ITDZBA0
INTEGER(KIND=JPIM) :: ITDZBS, ITDZBS0
INTEGER(KIND=JPIM) :: DTDZBA, DTDZBA0
INTEGER(KIND=JPIM) :: DTDZBS, DTDZBS0
INTEGER(KIND=JPIM) :: DTDZCA, DTDZCA0
INTEGER(KIND=JPIM) :: DTDZCS, DTDZCS0
INTEGER(KIND=JPIM) :: ITDZCA, ITDZCA0
INTEGER(KIND=JPIM) :: ITDZCS, ITDZCS0

REAL(KIND=JPRB),ALLOCATABLE, TARGET :: ZIA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZEPSNM(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSOA1(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZAOA1(:,:,:)
INTEGER(KIND=JPIM),ALLOCATABLE :: ISTAN(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE :: ISTAS(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSIA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZAIA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE, TARGET :: ZOA1(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE, TARGET :: ZOA2(:,:,:)

END MODULE TPM_FIELDS
