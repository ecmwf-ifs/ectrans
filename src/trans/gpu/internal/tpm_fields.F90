! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FIELDS

USE PARKIND_ECTRANS  ,ONLY : JPIM, JPIB, JPRB, JPRBT, JPRD
USE ISO_C_BINDING

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRD)    ,ALLOCATABLE :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRD)    ,ALLOCATABLE :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRBT)   ,ALLOCATABLE :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRBT)   ,ALLOCATABLE :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRBT)   ,ALLOCATABLE :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRBT)   ,ALLOCATABLE :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRBT)   ,ALLOCATABLE :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRBT)   ,ALLOCATABLE :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTN(:)   ! R%NTMAX+2-JN

REAL(KIND=JPRBT)   ,ALLOCATABLE :: RMU2(:)   ! sin(theta) for dual input/output latitudes
REAL(KIND=JPRBT)   ,ALLOCATABLE :: RACTHE2(:)! 1./SQRT(R1MU2), 1/(cos(theta)) dual input/output latitudes
END TYPE FIELDS_TYPE

!flat copies of the above
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RLAPIN(:) ! eigen-values of the inverse Laplace operator
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RACTHE(:) ! eigen-values of the inverse Laplace operator

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

! scratch arrays for ltinv and ltdir and associated dimension variables

REAL(KIND=JPRBT),ALLOCATABLE :: ZAA(:,:,:)  !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: ZAS(:,:,:)  !! JPRL for 1/2

! for m=0 in ledir_mod:
REAL(KIND=JPRD),ALLOCATABLE :: ZAA0(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAS0(:,:)
INTEGER(KIND=JPIM) :: KMLOC0

REAL(KIND=JPRBT),ALLOCATABLE :: ZEPSNM(:,:)

END MODULE TPM_FIELDS
