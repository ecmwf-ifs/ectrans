MODULE FTDIR_CTLAD_MOD
CONTAINS
SUBROUTINE FTDIR_CTLAD(PGP,KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS, &
 & KVSETUV,KVSETSC,KPTRGP)


!**** *FTDIR_CTLAD - Direct Fourier transform control - adjoint

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTLAD(..)

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

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_TRANS
USE TPM_DISTR

USE TRLTOG_MOD
USE FOURIER_OUTAD_MOD
USE FTDIRAD_MOD

IMPLICIT NONE

! Dummy arguments
REAL_B , INTENT(OUT) :: PGP(:,:,:)
INTEGER_M,INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPTRGP(:)

! Local variables
REAL_B,ALLOCATABLE :: ZGTF(:,:)


INTEGER_M :: IST
INTEGER_M :: IVSETUV(KF_UV_G)
INTEGER_M :: IVSETSC(KF_SCALARS_G)
INTEGER_M :: IVSET(KF_GP)
INTEGER_M :: JGL,IGL

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space

CALL GSTATS(133,0)
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

! Save Fourier data in FOUBUF_IN

CALL GSTATS(1642,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,IGL)
DO JGL=1,D%NDGL_FS
  IGL = JGL
  CALL FOURIER_OUTAD(ZGTF,KF_FS,IGL)

! Fourier transform

  IF(KF_FS>0) THEN
    CALL FTDIRAD(ZGTF,KF_FS,IGL)
  ENDIF
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1642,1)
CALL GSTATS(133,1)

IF(.NOT.LALLOPERM) THEN
  DEALLOCATE(FOUBUF_IN)
ENDIF

! Transposition

CALL GSTATS(183,0)
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
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
  IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
  IST = IST+KF_UV_G
ENDIF
IF(KF_SCALARS_G > 0) THEN
  IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+KF_SCALARS_G
ENDIF
CALL TRLTOG(ZGTF,PGP,KF_FS,KF_GP,IVSET,KPTRGP)

IF(.NOT.LALLOPERM)THEN
  DEALLOCATE(ZGTF)
ENDIF
CALL GSTATS(183,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTLAD
END MODULE FTDIR_CTLAD_MOD




