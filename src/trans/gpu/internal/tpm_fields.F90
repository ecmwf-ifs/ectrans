! (C) Copyright 2000- ECMWF.
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

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRD) ,ALLOCATABLE :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRD) ,ALLOCATABLE :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRBT) ,ALLOCATABLE :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRBT) ,ALLOCATABLE :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRBT) ,ALLOCATABLE :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRBT) ,ALLOCATABLE :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRBT) ,ALLOCATABLE :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRBT) ,ALLOCATABLE :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTN(:) ! R%NTMAX+2-JN

REAL(KIND=JPRBT) ,ALLOCATABLE :: RMU2(:)    ! sin(theta) for dual input/output latitudes
REAL(KIND=JPRBT) ,ALLOCATABLE :: RACTHE2(:) ! 1./SQRT(R1MU2), 1/(cos(theta)) dual input/output latitudes
END TYPE FIELDS_TYPE

!flat copies of the above
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RW(:)     ! Weights of the Gaussian quadrature

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

! scratch arrays for ltinv and ltdir and associated dimension variables

REAL(KIND=JPRBT),ALLOCATABLE :: ZAA(:,:,:)  !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: ZAS(:,:,:)  !! JPRL for 1/2
! for m=0 in ledir_mod:
REAL(KIND=JPRD),ALLOCATABLE :: ZAA0(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAS0(:,:)
INTEGER(KIND=JPIM) :: KMLOC0

INTEGER(KIND=JPIM) :: TDZAA
INTEGER(KIND=JPIM) :: TDZAS

! enable calling setup_trans with a different set of fields than inv_trans and dir_trans:
! IF_FS_INV0: size used for the allocation in setup_trans
! IF_FS_INV: size used in inv_trans and dir_Trans, needs to be <= IF_FS_INV0 
INTEGER(KIND=JPIM) :: IF_FS_INV, IF_FS_INV0
INTEGER(KIND=JPIM) :: IF_FS_DIR, IF_FS_DIR0
INTEGER(KIND=JPIM) :: NFLEV, NFLEV0

REAL(KIND=JPRB),ALLOCATABLE, TARGET :: ZIA(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZEPSNM(:,:)
REAL(KIND=JPRBT),ALLOCATABLE, TARGET :: ZOA1(:,:,:)
REAL(KIND=JPRBT),ALLOCATABLE, TARGET :: ZOA2(:,:,:)

! TODO find a better place for this
INTEGER(JPIM) :: ZGTF_START(8)
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_VOR = 1
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_DIV = 2
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_UV = 3
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_SCALAR = 4
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_NSDERS = 5
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_UVDERS = 6
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_EWDERS = 7
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_END = 8

END MODULE TPM_FIELDS
