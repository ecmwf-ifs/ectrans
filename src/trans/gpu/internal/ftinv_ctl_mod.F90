! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_CTL_MOD
CONTAINS
SUBROUTINE FTINV_CTL(KF_UV_G,KF_SCALARS_G,&
 & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,KVSETUV,KVSETSC,KPTRGP, &
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *FTINV_CTL - Inverse Fourier transform control

!     Purpose. Control routine for Fourier to Gridpoint transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV_CTL(..)

!        Explicit arguments :
!        --------------------
!        PGP     -  gridpoint array
!        KF_UV_G      - global number of spectral u-v fields
!        KF_SCALARS_G - global number of scalar spectral fields
!        KF_UV        - local number of spectral u-v fields
!        KF_SCALARS   - local number of scalar spectral fields
!        KF_SCDERS    - local number of derivatives of scalar spectral fields
!        KF_GP        - total number of output gridpoint fields
!        KF_FS        - total number of fields in fourier space
!        KF_OUT_LT    - total number of fields coming out from inverse LT
!        KVSETUV - "B"  set in spectral/fourier space for
!                   u and v variables
!        KVSETSC - "B" set in spectral/fourier space for
!                  scalar variables
!        KPTRGP - pointer array to fi3elds in gridpoint space

!     Method.
!     -------

!     Externals.  TRLTOG      - transposition routine
!     ----------  FOURIER_IN  - copy fourier data from Fourier buffer
!                 FTINV       - fourier transform
!                 FSC         - Fourier space computations

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,   JPRBT

USE TPM_GEN         ,ONLY : NERR, nout
!USE TPM_DIM
!USE TPM_GEOMETRY
USE TPM_TRANS       ,ONLY : FOUBUF, LDIVGP, LSCDERS, LUVDER, LVORGP,LATLON, ZGTF
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC
USE TPM_FLT         ,ONLY : S
USE FOURIER_IN_MOD  ,ONLY : FOURIER_IN
USE FSC_MOD         ,ONLY : FSC
USE FTINV_MOD       ,ONLY : FTINV
USE TRLTOG_MOD      ,ONLY : TRLTOG, TRLTOG_CUDAAWARE
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TPM_FIELDS      ,ONLY : ZGTF_START, ZGTF_START_INDEX_UV, ZGTF_START_INDEX_NSDERS, &
    & ZGTF_START_INDEX_UVDERS, ZGTF_START_INDEX_EWDERS, ZGTF_START_INDEX_SCALAR, &
    & ZGTF_START_INDEX_VOR, ZGTF_START_INDEX_DIV, ZGTF_START_INDEX_END
use ieee_arithmetic
!

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCDERS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP2(:,:,:)

INTEGER(KIND=JPIM) :: IST
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: J3,JGL,IGL,IOFF,IFGP2,IFGP3A,IFGP3B,IGP3APAR,IGP3BPAR

INTEGER(KIND=JPIM) :: JF_FS


!     ------------------------------------------------------------------

!    1.  Copy Fourier data to local array

CALL GSTATS(107,0)

! Figure out where we want to store data in ZGTF
IST = 1
ZGTF_START(ZGTF_START_INDEX_VOR) = IST
IF (LVORGP) THEN
  IST = IST+KF_UV
ENDIF
ZGTF_START(ZGTF_START_INDEX_DIV) = IST
IF (LDIVGP) THEN
  IST = IST+KF_UV
ENDIF
ZGTF_START(ZGTF_START_INDEX_UV) = IST
IST = IST+2*KF_UV
ZGTF_START(ZGTF_START_INDEX_SCALAR) = IST
IST = IST+KF_SCALARS
ZGTF_START(ZGTF_START_INDEX_NSDERS) = IST
IST = IST+KF_SCDERS
ZGTF_START(ZGTF_START_INDEX_UVDERS) = IST
IF (LUVDER) THEN
  IST = IST+2*KF_UV
ENDIF
ZGTF_START(ZGTF_START_INDEX_EWDERS) = IST
IF (KF_SCDERS > 0) THEN
  IST = IST+KF_SCDERS
ENDIF
ZGTF_START(ZGTF_START_INDEX_END) = IST


CALL GSTATS(1639,0)
! from FOUBUF to ZGTF
CALL FOURIER_IN(ZGTF,KF_OUT_LT)

!    2.  Fourier space computations

IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
    CALL FSC
ENDIF

!   3.  Fourier transform
IF(KF_FS > 0) THEN
  ! from ZGTF to ZGTF
  CALL FTINV(ZGTF,size(zgtf,1))
ENDIF

CALL GSTATS(1639,1)

CALL GSTATS(107,1)

!   4.  Transposition

IF (PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IVSETSC(:)=-1
IF (PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF (PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF (PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    IGP3APAR=UBOUND(PGP3A,3)
    IF (LSCDERS) IGP3APAR=IGP3APAR/3
    DO J3=1,IGP3APAR
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF (PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    IGP3BPAR=UBOUND(PGP3B,3)
    IF (LSCDERS) IGP3BPAR=IGP3BPAR/3
    DO J3=1,IGP3BPAR
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
  IF (IOFF > 0 .AND. IOFF /= KF_SCALARS_G ) THEN
    WRITE(NERR,*)'FTINV:IOFF,KF_SCALARS_G ',IOFF,KF_SCALARS_G
    CALL ABORT_TRANS('FTINV_CTL_MOD:IOFF /= KF_SCALARS_G')
  ENDIF
ENDIF

IST = 1
IF (KF_UV_G > 0) THEN
  IF (LVORGP) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF ( LDIVGP) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF (KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
  IF (LSCDERS) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF
ENDIF
IF (KF_UV_G > 0 .AND. LUVDER) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF (KF_SCALARS_G > 0) THEN
  IF (LSCDERS) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF
ENDIF

CALL GSTATS(157,0)
JF_FS=KF_FS-D%IADJUST_I
#ifdef USE_CUDA_AWARE_MPI_FT
WRITE(NOUT,*) 'ftinv_ctl:TRLTOG_CUDAAWARE'
CALL TRLTOG_CUDAAWARE(ZGTF,JF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
#else
!WRITE(NOUT,*) 'ftinv_ctl:TRLTOG'
!$ACC UPDATE HOST(ZGTF)
CALL TRLTOG(ZGTF,JF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
#endif
CALL GSTATS(157,1)
!     ------------------------------------------------------------------

!DEALLOCATE(ZGTF)

END SUBROUTINE FTINV_CTL
END MODULE FTINV_CTL_MOD
