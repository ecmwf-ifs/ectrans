! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_VIEW_CTL_MOD
CONTAINS
SUBROUTINE FTDIR_VIEW_CTL(KF_GP,KF_FS,YDGP)


!**** *FTDIR_VIEW_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_VIEW_CTL(..)

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
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!      R. El Khatib 01-Jun-2022 contiguous pointer
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN   ,ONLY : NSTACK_MEMORY_TR
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC

USE TRGTOL_VIEW_MOD      ,ONLY : TRGTOL_VIEW
USE FOURIER_OUT_MOD ,ONLY : FOURIER_OUT
USE FTDIR_MOD       ,ONLY : FTDIR
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : GRID_VIEW

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM),INTENT(IN) ::KF_GP,KF_FS
TYPE(GRID_VIEW),INTENT(INOUT) :: YDGP(:)
! Local variables
REAL(KIND=JPRB),TARGET  :: ZGTF_STACK(KF_FS*MIN(1,MAX(0,NSTACK_MEMORY_TR)),D%NLENGTF)
REAL(KIND=JPRB),TARGET, ALLOCATABLE :: ZGTF_HEAP(:,:)
REAL(KIND=JPRB),POINTER, CONTIGUOUS :: ZGTF(:,:)

INTEGER(KIND=JPIM) :: JGL,IBLEN
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: I, KI

!     ------------------------------------------------------------------

IF (NSTACK_MEMORY_TR == 1) THEN
  ZGTF => ZGTF_STACK(:,:)
ELSE
  ALLOCATE(ZGTF_HEAP(KF_FS,D%NLENGTF))
! Now, force the OS to allocate this shared array right now, not when it starts
! to be used which is an OPEN-MP loop, that would cause a threads
! synchronization lock :
  IF (KF_FS > 0 .AND. D%NLENGTF > 0) THEN
    ZGTF_HEAP(1,1)=HUGE(1._JPRB)
  ENDIF
  ZGTF => ZGTF_HEAP(:,:)
ENDIF

DO I=1,SIZE(YDGP)
  IVSET(I) = YDGP(I)%IVSET  
ENDDO

! Transposition

CALL GSTATS(158,0)
CALL TRGTOL_VIEW(ZGTF,KF_FS,KF_GP,IVSET,YDGP)
CALL GSTATS(158,1)
CALL GSTATS(106,0)

! Fourier transform

IBLEN=D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
ENDIF

CALL GSTATS(1640, 0)
! If this rank has any Fourier fields, Fourier transform them
IF (KF_FS > 0) THEN
  ! Loop over latitudes
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL)
  DO JGL = 1, D%NDGL_FS
    ! Fourier transform
    CALL FTDIR(ZGTF, KF_FS, JGL)

    ! Save Fourier data in FOUBUF_IN
    CALL FOURIER_OUT(ZGTF, KF_FS, JGL)
  ENDDO
  !$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1640, 1)

CALL GSTATS(106,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_VIEW_CTL
END MODULE FTDIR_VIEW_CTL_MOD



