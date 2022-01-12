! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FIELDS

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD, JPRL

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRD) ,ALLOCATABLE :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRD) ,ALLOCATABLE :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRB) ,ALLOCATABLE :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRB) ,ALLOCATABLE :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRB) ,ALLOCATABLE :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRB) ,ALLOCATABLE :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRB) ,ALLOCATABLE :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRB) ,ALLOCATABLE :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTN(:) ! R%NTMAX+2-JN

REAL(KIND=JPRB) ,ALLOCATABLE :: RMU2(:)    ! sin(theta) for dual input/output latitudes
REAL(KIND=JPRB) ,ALLOCATABLE :: RACTHE2(:) ! 1./SQRT(R1MU2), 1/(cos(theta)) dual input/output latitudes
END TYPE FIELDS_TYPE

!flat copies of the above
REAL(KIND=JPRB) ,ALLOCATABLE :: F_RW(:)     ! Weights of the Gaussian quadrature

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

! scratch arrays for ltinv and ltdir and associated dimension variables

REAL(KIND=JPRL),ALLOCATABLE :: ZAA(:,:,:)
REAL(KIND=JPRL),ALLOCATABLE :: ZAS(:,:,:)

REAL(KIND=JPRL), POINTER :: IZBA(:,:,:)
REAL(KIND=JPRL),ALLOCATABLE,TARGET  :: IZBS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: IZCA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: IZCS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: IZCAT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: IZCST(:,:,:)

REAL(KIND=JPRB),ALLOCATABLE :: DZBA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZBS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZBAT(:,:,:)
REAL(KIND=JPRL),ALLOCATABLE :: DZBST(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZCA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZCS(:,:,:)
!REAL(KIND=JPRB),POINTER :: DZCAT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: DZCAT(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: DZCST(:,:,:)

! Arrays used for rescaling to allow half-precision Legende transforms
REAL(KIND=JPRB), ALLOCATABLE :: ZAMAX(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSMAX(:,:)

INTEGER :: LDZAA
INTEGER :: LDZAS
INTEGER :: TDZAA
INTEGER :: TDZAS

INTEGER :: ILDZBA
INTEGER :: ILDZBS
INTEGER :: ITDZBA
INTEGER :: ITDZBS
INTEGER :: ILDZCA
INTEGER :: ILDZCS
INTEGER :: ITDZCA
INTEGER :: ITDZCS



INTEGER :: DLDZBA
INTEGER :: DLDZBS
INTEGER :: DTDZBA
INTEGER :: DTDZBS
INTEGER :: DLDZCA
INTEGER :: DLDZCS
INTEGER :: DTDZCA
INTEGER :: DTDZCS

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
