MODULE FTINV_CTL_MOD
CONTAINS
SUBROUTINE FTINV_CTL(PGP,KF_UV_G,KF_SCALARS_G,&
 & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,KVSETUV,KVSETSC,KPTRGP)


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

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_TRANS
USE TPM_DISTR

USE FOURIER_IN_MOD
USE FSC_MOD
USE FTINV_MOD
USE TRLTOG_MOD

IMPLICIT NONE

REAL_B , INTENT(OUT) :: PGP(:,:,:)
INTEGER_M ,INTENT(IN) :: KF_UV_G,KF_SCALARS_G
INTEGER_M ,INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPTRGP(:)

REAL_B,ALLOCATABLE,TARGET  :: ZGTF(:,:)
REAL_B,TARGET  :: ZDUM(1,D%NLENGTF)
REAL_B,POINTER :: ZUV(:,:)
REAL_B,POINTER :: ZSCALAR(:,:)
REAL_B,POINTER :: ZNSDERS(:,:)
REAL_B,POINTER :: ZEWDERS(:,:)
REAL_B,POINTER :: ZUVDERS(:,:)

INTEGER_M :: IST
INTEGER_M :: IVSETUV(KF_UV_G)
INTEGER_M :: IVSETSC(KF_SCALARS_G)
INTEGER_M :: IVSET(KF_GP)
INTEGER_M :: J1,J2,JGL,IGL

!     ------------------------------------------------------------------

!    1.  Copy Fourier data to local array

CALL GSTATS(107,0)
IF(LALLOPERM) THEN
  IF(ALLOCATED(ZGTF)) THEN
    IF( .NOT. (SIZE(ZGTF,1) == KF_FS .AND. SIZE(ZGTF,2) ==  D%NLENGTF) ) THEN
      DEALLOCATE(ZGTF)
    ENDIF
  ENDIF
ENDIF
IF(.NOT.ALLOCATED(ZGTF)) THEN
  ALLOCATE(ZGTF(KF_FS,D%NLENGTF))
ENDIF

IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  IST = 1
  IF(LVORGP) THEN
    IST = IST+KF_UV
  ENDIF
  IF(LDIVGP) THEN
    IST = IST+KF_UV
  ENDIF
  ZUV => ZGTF(IST:IST+2*KF_UV-1,:)
  IST = IST+2*KF_UV
  ZSCALAR => ZGTF(IST:IST+KF_SCALARS-1,:)
  IST = IST+KF_SCALARS
  ZNSDERS => ZGTF(IST:IST+KF_SCDERS-1,:)
  IST = IST+KF_SCDERS
  IF(LUVDER) THEN
    ZUVDERS => ZGTF(IST:IST+2*KF_UV-1,:)
    IST = IST+2*KF_UV
  ELSE
    ZUVDERS => ZDUM(1:1,:)  
  ENDIF
  IF(KF_SCDERS > 0) THEN
    ZEWDERS => ZGTF(IST:IST+KF_SCDERS-1,:)
  ELSE
    ZEWDERS => ZDUM(1:1,:)
  ENDIF
ENDIF

CALL GSTATS(1039,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JGL,IGL)
DO JGL=1,D%NDGL_FS
  IGL = JGL
  CALL FOURIER_IN(ZGTF,KF_OUT_LT,IGL)

!    2.  Fourier space computations


  IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
    CALL FSC(IGL,KF_UV,KF_SCALARS,KF_SCDERS,&
     & ZUV,ZSCALAR,ZNSDERS,ZEWDERS,ZUVDERS)
  ENDIF

!   3.  Fourier transform
  IF(KF_FS > 0) THEN
    CALL FTINV(ZGTF,KF_FS,IGL)
  ENDIF
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1039,1)

IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  NULLIFY(ZUV)
  NULLIFY(ZSCALAR)
  NULLIFY(ZNSDERS)
  NULLIFY(ZUVDERS)
  NULLIFY(ZEWDERS)
ENDIF
CALL GSTATS(107,1)

!   4.  Transposition

IF(PRESENT(KVSETUV)) THEN
  IVSETUV(:) = KVSETUV(:)
ELSE
  IVSETUV(:) = -1
ENDIF
IF(PRESENT(KVSETSC)) THEN
  IVSETSC(:) = KVSETSC(:)
ELSE
  IVSETSC(:) = -1
ENDIF

IST = 1
IF(KF_UV_G > 0) THEN
  IF( LVORGP) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF( LDIVGP) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
  IF(LSCDERS) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF
ENDIF
IF(KF_UV_G > 0 .AND. LUVDER) THEN
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IF(LSCDERS) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF
ENDIF

CALL GSTATS(157,0)
CALL TRLTOG(ZGTF,PGP,KF_FS,KF_GP,IVSET,KPTRGP)
CALL GSTATS(157,1)

IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(ZGTF)
  DEALLOCATE(FOUBUF)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTINV_CTL
END MODULE FTINV_CTL_MOD




