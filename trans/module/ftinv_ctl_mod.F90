MODULE FTINV_CTL_MOD
CONTAINS
SUBROUTINE FTINV_CTL(PGP,KVSETUV,KVSETSC)


!**** *FTINV_CTL - Inverse Fourier transform control

!     Purpose. Control routine for Fourier to Gridpoint transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV_CTL(..)

!        Explicit arguments :  PGP     -  gridpoint array
!        --------------------  KVSETUV - "B" set in spectral/fourier space for
!                                         u and v variables
!                              KVSETSC - "B" set in spectral/fourier space for
!                                         scalar variables

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
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)

REAL_B,ALLOCATABLE,TARGET  :: ZGTF(:,:)
REAL_B,TARGET  :: ZDUM(1,D%NLENGTF)
REAL_B,POINTER :: ZUV(:,:)
REAL_B,POINTER :: ZSCALAR(:,:)
REAL_B,POINTER :: ZNSDERS(:,:)
REAL_B,POINTER :: ZEWDERS(:,:)
REAL_B,POINTER :: ZUVDERS(:,:)

INTEGER_M :: IST
INTEGER_M :: IVSETUV(NF_UV_G)
INTEGER_M :: IVSETSC(NF_SCALARS_G)
INTEGER_M :: IVSET(NF_GP)
INTEGER_M :: J1,J2

!     ------------------------------------------------------------------

!    1.  Copy Fourier data to local array

CALL GSTATS(107,0)
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
!$OMP PARALLEL DO PRIVATE(J,I)
DO J1=1,D%NLENGTF
  DO J2=1,NF_FS
    ZGTF(J2,J1) = _ZERO_
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL FOURIER_IN(ZGTF,NF_OUT_LT)

!    2.  Fourier space computations

IF(NF_UV > 0 .OR. NF_SCDERS > 0) THEN
  IST=1
  IF(LVORGP) THEN
    IST = IST+NF_UV
  ENDIF
  IF(LDIVGP) THEN
    IST = IST+NF_UV
  ENDIF
  ZUV => ZGTF(IST:IST+2*NF_UV-1,:)
  IST = IST+2*NF_UV
  ZSCALAR => ZGTF(IST:IST+NF_SCALARS-1,:)
  IST = IST+NF_SCALARS
  ZNSDERS => ZGTF(IST:IST+NF_SCDERS-1,:)
  IST = IST+NF_SCDERS
  IF(LUVDER) THEN
    ZUVDERS => ZGTF(IST:IST+2*NF_UV-1,:)
    IST = IST+2*NF_UV
  ELSE
    ZUVDERS => ZDUM(1:1,:)  
  ENDIF
  IF(NF_SCDERS > 0) THEN
    ZEWDERS => ZGTF(IST:IST+NF_SCDERS-1,:)
  ELSE
    ZEWDERS => ZDUM(1:1,:)
  ENDIF

  CALL FSC(ZUV,ZSCALAR,ZNSDERS,ZEWDERS,ZUVDERS)

  NULLIFY(ZUV)
  NULLIFY(ZSCALAR)
  NULLIFY(ZNSDERS)
  NULLIFY(ZUVDERS)
  NULLIFY(ZEWDERS)
ENDIF

!   3.  Fourier transform
CALL FTINV(ZGTF,NF_FS)
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
IF(NF_UV_G > 0) THEN
  IF( LVORGP) THEN
    IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
    IST = IST+NF_UV_G
  ENDIF
  IF( LDIVGP) THEN
    IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
    IST = IST+NF_UV_G
  ENDIF
  IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
  IST = IST+NF_UV_G
  IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
  IST = IST+NF_UV_G
ENDIF
IF(NF_SCALARS_G > 0) THEN
  IVSET(IST:IST+NF_SCALARS_G-1) = IVSETSC(:)
  IST = IST+NF_SCALARS_G
  IF(NF_SCDERS > 0) THEN
    IVSET(IST:IST+NF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+NF_SCALARS_G
  ENDIF
ENDIF
IF(NF_UV_G > 0 .AND. LUVDER) THEN
  IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
  IST = IST+NF_UV_G
  IVSET(IST:IST+NF_UV_G-1) = IVSETUV(:)
  IST = IST+NF_UV_G
ENDIF
IF(NF_SCALARS_G > 0) THEN
  IF(NF_SCDERS > 0) THEN
    IVSET(IST:IST+NF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+NF_SCALARS_G
  ENDIF
ENDIF

CALL GSTATS(157,0)
CALL TRLTOG(ZGTF,PGP,IVSET)
CALL GSTATS(157,1)

IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(ZGTF)
  DEALLOCATE(FOUBUF)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTINV_CTL
END MODULE FTINV_CTL_MOD




