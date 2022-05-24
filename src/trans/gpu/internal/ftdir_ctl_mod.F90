! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_CTL_MOD
CONTAINS
SUBROUTINE FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS, &
 & KVSETUV,KVSETSC,KPTRGP,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *FTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTL(..)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     PGP     -  gridpoint array
!     KVSETUV - "B" set in spectral/fourier space for
!                u and v variables
!     KVSETSC - "B" set in spectral/fourier space for
!                scalar variables
!     KPTRGP  -  pointer array to fields in gridpoint space

!     Method.
!     -------

!     Externals.  TRGTOL      - transposition routine
!     ----------  FOURIER_OUT - copy fourier data to Fourier buffer
!                 FTDIR       - fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT

USE TPM_GEN, only: nout
!USE TPM_DIM
!USE TPM_GEOMETRY
USE TPM_TRANS       ,ONLY : ZGTF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC

USE TRGTOL_MOD      ,ONLY : TRGTOL, TRGTOL_CUDAAWARE
USE FTDIR_MOD       ,ONLY : FTDIR
use ieee_arithmetic
!

IMPLICIT NONE


INTERFACE
  SUBROUTINE cudaProfilerStart() BIND(C,name='cudaProfilerStart')
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
    IMPLICIT NONE
  END SUBROUTINE cudaProfilerStart
END INTERFACE

INTERFACE
  SUBROUTINE cudaProfilerStop() BIND(C,name='cudaProfilerStop')
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
    IMPLICIT NONE
  END SUBROUTINE cudaProfilerStop
END INTERFACE

! Dummy arguments

INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)

! Local variables
!REAL(KIND=JPRBT),ALLOCATABLE :: ZGTF(:,:)

INTEGER(KIND=JPIM) :: IST,JGL,IGL,JF_FS
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
INTEGER(KIND=JPIM) :: ISIZE,IFIELDS,ICHUNK,ICHUNKS,JK

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space

!call cudaProfilerStart()

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IVSETSC(:) = -1
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IOFF=0
  IF(PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF
  IF(PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    DO J3=1,UBOUND(PGP3A,3)
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF
  IF(PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    DO J3=1,UBOUND(PGP3B,3)
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF
ENDIF

IST = 1
IF(KF_UV_G > 0) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
ENDIF

! Transposition

CALL GSTATS(158,0)

! needed ??? JF_FS=KF_FS-D%IADJUST_D
#ifdef USE_CUDA_AWARE_MPI_FT
CALL GSTATS(430,0)
!$ACC DATA IF(PRESENT(PGP))   COPYIN(PGP)
!$ACC DATA IF(PRESENT(PGPUV)) COPYIN(PGPUV)
!$ACC DATA IF(PRESENT(PGP2))  COPYIN(PGP2)
!$ACC DATA IF(PRESENT(PGP3A)) COPYIN(PGP3A)
!$ACC DATA IF(PRESENT(PGP3B)) COPYIN(PGP3B)
CALL GSTATS(430,1)
CALL TRGTOL_CUDAAWARE(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
#else
CALL TRGTOL(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)
#endif

CALL GSTATS(158,1)
CALL GSTATS(106,0)

! Fourier transform

!write(301,*) 'sizey: ', myproc, size(zgtf,1), KF_FS

CALL GSTATS(1640,0)
IF (KF_FS > 0) THEN
  CALL FTDIR(SIZE(ZGTF,1),KF_FS)
ENDIF

CALL GSTATS(1640,1)
!DEALLOCATE(ZGTF)
CALL GSTATS(106,1)
!     ------------------------------------------------------------------
!call cudaProfilerStop()
END SUBROUTINE FTDIR_CTL
END MODULE FTDIR_CTL_MOD

