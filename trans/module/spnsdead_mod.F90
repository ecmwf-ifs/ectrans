MODULE SPNSDEAD_MOD
CONTAINS
SUBROUTINE SPNSDEAD(KM,KF_SCALARS,PEPSNM,PF,PNSD)

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_FIELDS
USE TPM_TRANS

#ifdef DOC

!**** *SPNSDEAD* - Compute North-South derivative in spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the North-south derivative

!**   Interface.
!     ----------
!        CALL SPNSDEAD(...)

!        Explicit arguments :
!        -------------------- 
!        KM -zonal wavenumber (input-c)
!        PEPSNM - REPSNM for wavenumber KM (input-c)
!        PF  (NLEI1,2*KF_SCALARS) - input field (input)
!        PNSD(NLEI1,2*KF_SCALARS) - N-S derivative (output)

!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  YOMLAP
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From SPNSDEAD in IFS CY22R1

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER_M, INTENT(IN)    :: KM,KF_SCALARS
REAL_B,    INTENT(IN)    :: PEPSNM(0:R%NTMAX+2)
REAL_B,    INTENT(INOUT) :: PF(:,:)
REAL_B,    INTENT(IN)    :: PNSD(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IJ, ISKIP, J, JN,II


!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------


!*       1.1      COMPUTE


IF(KM == 0) THEN
  ISKIP = 2
ELSE
  ISKIP = 1
ENDIF

DO JN=KM,R%NSMAX+1
  IJ = R%NSMAX+3-JN
  DO J=1,2*KF_SCALARS,ISKIP
    PF(IJ+1,J) = PF(IJ+1,J)-F%RN(JN-1)*PEPSNM(JN  )*PNSD(IJ,J)
    PF(IJ-1,J) = PF(IJ-1,J)+F%RN(JN+2)*PEPSNM(JN+1)*PNSD(IJ,J)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SPNSDEAD
END MODULE SPNSDEAD_MOD
