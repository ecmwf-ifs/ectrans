! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_VIEW_CTL_MOD
CONTAINS

SUBROUTINE FTINV_VIEW_CTL(KF_GP,KF_FS,KF_OUT_LT,&
 & YDGVU,YDGVV,YDGVSCALAR,&
 & YDGVVOR,YDGVDIV,&
 & YDGVU_EW,YDGVV_EW,YDGVSCALAR_NS,YDGVSCALAR_EW)


!**** *FTINV_VIEW_CTL - Inverse Fourier transform control

!     Purpose. Control routine for Fourier to Gridpoint transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV_VIEW_CTL(..)

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
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NERR   ,NSTACK_MEMORY_TR, LALLOPERM
USE TPM_TRANS       ,ONLY : FOUBUF, LDIVGP, LSCDERS, LUVDER, LVORGP,LATLON
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC, MYSETV
USE TPM_FLT         ,ONLY : S
USE FOURIER_IN_MOD  ,ONLY : FOURIER_IN
USE FSC_VIEW_MOD    ,ONLY : FSC_VIEW
USE FTINV_MOD       ,ONLY : FTINV
USE TRLTOG_VIEW_MOD      ,ONLY : TRLTOG_VIEW
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD ,ONLY : GRID_VIEW, SPEC_VIEW

IMPLICIT NONE
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_OUT_LT

TYPE(GRID_VIEW),INTENT(INOUT) :: YDGVU(:),YDGVV(:)
TYPE(GRID_VIEW),INTENT(INOUT) :: YDGVSCALAR(:)
TYPE(GRID_VIEW),INTENT(INOUT) :: YDGVVOR(:),YDGVDIV(:)

TYPE(GRID_VIEW),INTENT(IN) :: YDGVU_EW(:),YDGVV_EW(:)
TYPE(GRID_VIEW),INTENT(IN) :: YDGVSCALAR_NS(:), YDGVSCALAR_EW(:)

! Local variables

TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPU(:),YLSPV(:),YLSPVSCALAR(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVOR(:),YLSPDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPU_EW(:),YLSPV_EW(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVSCALAR_NS(:), YLSPVSCALAR_EW(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPUV(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPUV_EW(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGV(:)

INTEGER(KIND=JPIM) :: I, KI
INTEGER(KIND=JPIM) :: JGL,IGL,IOFF
INTEGER(KIND=JPIM) :: IVSET(KF_GP)

REAL(KIND=JPRB),TARGET  :: ZGTF_STACK(KF_FS*MIN(1,MAX(0,NSTACK_MEMORY_TR)),D%NLENGTF)
REAL(KIND=JPRB),TARGET, ALLOCATABLE :: ZGTF_HEAP(:,:)
REAL(KIND=JPRB),POINTER  :: ZGTF(:,:)
INTEGER(KIND=JPIM) :: ICOUNT

!     ------------------------------------------------------------------

!    1.  Copy Fourier data to local array

CALL GSTATS(107,0)

IF (NSTACK_MEMORY_TR == 1) THEN
  ZGTF => ZGTF_STACK(:,:)
ELSE
  ALLOCATE(ZGTF_HEAP(KF_FS,D%NLENGTF))
! Now, force the OS to allocate this shared array right now, not when it starts
! to be used which is an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (KF_FS > 0 .AND. D%NLENGTF > 0) THEN
    ZGTF_HEAP(1,1)=HUGE(1._JPRB)
  ENDIF
  ZGTF => ZGTF_HEAP(:,:)
ENDIF

CALL GSTATS(1639,0)

! Keep the order: ordering in ZGTF has to follow the classical PGP array ordering
ICOUNT = 1
YLSPVOR = WRAP_G2S(YDGVVOR,ZGTF,ICOUNT)
YLSPDIV = WRAP_G2S(YDGVDIV,ZGTF,ICOUNT)
YLSPU = WRAP_G2S(YDGVU,ZGTF,ICOUNT)
YLSPV = WRAP_G2S(YDGVV,ZGTF,ICOUNT)
YLSPVSCALAR = WRAP_G2S(YDGVSCALAR,ZGTF,ICOUNT)
YLSPVSCALAR_NS = WRAP_G2S(YDGVSCALAR_NS,ZGTF,ICOUNT)
YLSPU_EW = WRAP_G2S(YDGVU_EW,ZGTF,ICOUNT)
YLSPV_EW = WRAP_G2S(YDGVV_EW,ZGTF,ICOUNT)
YLSPVSCALAR_EW = WRAP_G2S(YDGVSCALAR_EW,ZGTF,ICOUNT)

IF (ICOUNT - 1 /= KF_FS) THEN
  WRITE(NERR, *) 'FTINV_VIEW_CTL: mismatched local field count: ICOUNT-1=', ICOUNT-1, ' KF_FS=', KF_FS
  CALL ABORT_TRANS('FTINV_VIEW_CTL: local KF_FS mismatch after WRAP_G2S filtering')
ENDIF

YLGV = [YDGVVOR,YDGVDIV,YDGVU,YDGVV,YDGVSCALAR,YDGVSCALAR_NS,YDGVU_EW,YDGVV_EW,YDGVSCALAR_EW]

DO I=1,SIZE(YLGV)
  IVSET(I) = YLGV(I)%IVSET
ENDDO

YLSPUV = [YLSPU, YLSPV]
YLSPUV_EW = [YLSPU_EW, YLSPV_EW]

! Loop over latitudes
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,IGL)
DO JGL = 1, D%NDGL_FS
  IGL = JGL
  CALL FOURIER_IN(ZGTF,KF_OUT_LT,IGL)

!    2.  Fourier space computations

  IF (SIZE (YDGVU) > 0 .OR. SIZE (YDGVSCALAR_NS) > 0 .OR. (LATLON.AND.S%LDLL) ) THEN
    CALL FSC_VIEW(IGL,&
                 & YLSPUV,YLSPVSCALAR,&
                 & YLSPUV_EW,YLSPVSCALAR_NS,YLSPVSCALAR_EW)
  ENDIF

!    3.  Fourier transform

  IF (KF_FS > 0) THEN
    CALL FTINV(ZGTF,KF_FS,IGL) ! Watch out failures here (Cray CCE 8.6.2 ? Intel 18.0.1 ?)
  ENDIF
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1639,1)

CALL GSTATS(107,1)

CALL GSTATS(157,0)
CALL TRLTOG_VIEW(ZGTF,KF_FS,KF_GP,IVSET,YLGV)
CALL GSTATS(157,1)

IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)

!     ------------------------------------------------------------------

!DEALLOCATE(ZGTF)
CONTAINS

FUNCTION WRAP_G2S(YGP,PZGTF,KII) RESULT(YSP)
  TYPE(GRID_VIEW), INTENT(IN) :: YGP(:)
  TYPE(SPEC_VIEW), ALLOCATABLE :: YSP(:)
  REAL(KIND=JPRB),TARGET :: PZGTF(:,:)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KII
  INTEGER(KIND=JPIM) :: II, ILOC, INLOC
  INLOC = 0
  DO II=1,SIZE(YGP)
    IF (YGP(II)%IVSET == MYSETV .OR. YGP(II)%IVSET == -1) INLOC = INLOC + 1
  ENDDO
  ALLOCATE(YSP(INLOC))
  ILOC = 0
  DO II=1,SIZE(YGP)
    IF (YGP(II)%IVSET == MYSETV .OR. YGP(II)%IVSET == -1) THEN
      ILOC = ILOC + 1
      YSP(ILOC)%P => PZGTF(KII,:)
      KII = KII + 1
    ENDIF
  ENDDO
END FUNCTION WRAP_G2S


END SUBROUTINE FTINV_VIEW_CTL
END MODULE FTINV_VIEW_CTL_MOD
