MODULE FTDIR_CTL_MOD
CONTAINS
SUBROUTINE FTDIR_CTL(PGP,KVSETUV,KVSETSC)


!**** *FTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR_CTL(..)

!        Explicit arguments :  PGP     -  gridpoint array
!        --------------------  KVSETUV - "B" set in spectral/fourier space for
!                                         u and v variables
!                              KVSETSC - "B" set in spectral/fourier space for
!                                         scalar variables

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

USE TRGTOL_MOD
USE FOURIER_OUT_MOD
USE FTDIR_MOD

IMPLICIT NONE

! Dummy arguments
REAL_B , INTENT(IN) :: PGP(:,:,:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)

! Local variables
REAL_B,ALLOCATABLE :: ZGTF(:,:)

INTEGER_M :: IST
INTEGER_M :: IVSETUV(NF_UV_G)
INTEGER_M :: IVSETSC(NF_SCALARS_G)
INTEGER_M :: IVSET(NF_GP)

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space

IF(LALLOPERM) THEN
  IF(ALLOCATED(ZGTF)) THEN
    IF( .NOT. (SIZE(ZGTF,1) == NF_FS .AND. SIZE(ZGTF,2) ==  D%NLENGTF) ) THEN
      DEALLOCATE(ZGTF)
    ENDIF
  ENDIF
ENDIF
IF(.NOT.ALLOCATED(ZGTF)) THEN
  ALLOCATE(ZGTF(NF_FS,D%NLENGTF))
ENDIF
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
IF(NF_UV_G > 0) THEN
  IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
  IST = IST+NF_UV_G
  IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
  IST = IST+NF_UV_G
ENDIF
IF(NF_SCALARS_G > 0) THEN
  IVSET(IST:IST+NF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+NF_SCALARS_G
ENDIF

! Transposition

CALL GSTATS(158,0)
CALL TRGTOL(ZGTF,PGP,IVSET)
CALL GSTATS(158,1)

! Fourier transform

CALL GSTATS(106,0)
CALL FTDIR(ZGTF,NF_FS)

! Save Fourier data in FOUBUF_IN

IF(LALLOPERM) THEN
  IF(ALLOCATED(FOUBUF_IN)) THEN
    IF(SIZE(FOUBUF_IN) < D%NLENGT0B*2*NF_FS) THEN
      DEALLOCATE(FOUBUF_IN)
    ENDIF
  ENDIF
ENDIF
IF(.NOT.ALLOCATED(FOUBUF_IN)) ALLOCATE(FOUBUF_IN(D%NLENGT0B*2*NF_FS))
CALL FOURIER_OUT(ZGTF,NF_FS)
CALL GSTATS(106,1)
IF(.NOT.LALLOPERM) THEN
  DEALLOCATE(ZGTF)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTL
END MODULE FTDIR_CTL_MOD




