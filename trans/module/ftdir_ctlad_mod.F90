MODULE FTDIR_CTLAD_MOD
CONTAINS
SUBROUTINE FTDIR_CTLAD(PGP,KVSETUV,KVSETSC)


!**** *FTDIR_CTLAD - Direct Fourier transform control - adjoint

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR_CTLAD(..)

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

USE TRLTOG_MOD
USE FOURIER_OUTAD_MOD
USE FTDIRAD_MOD

IMPLICIT NONE

! Dummy arguments
REAL_B , INTENT(OUT) :: PGP(:,:,:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)

! Local variables
REAL_B    :: ZGTF(NF_FS,D%NLENGTF)

INTEGER_M :: IST
INTEGER_M :: IVSETUV(NF_UV_G)
INTEGER_M :: IVSETSC(NF_SCALARS_G)
INTEGER_M :: IVSET(NF_GP)

!     ------------------------------------------------------------------

! Field distribution in Spectral/Fourier space

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


! Save Fourier data in FOUBUF_IN

CALL GSTATS(133,0)
CALL FOURIER_OUTAD(ZGTF,NF_FS)

IF(.NOT.LALLOPERM) THEN
  DEALLOCATE(FOUBUF_IN)
ENDIF

! Fourier transform

CALL FTDIRAD(ZGTF,NF_FS)
CALL GSTATS(133,1)


! Transposition

CALL GSTATS(183,0)
CALL TRLTOG(ZGTF,PGP,IVSET)
CALL GSTATS(183,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTLAD
END MODULE FTDIR_CTLAD_MOD




