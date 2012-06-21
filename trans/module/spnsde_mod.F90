MODULE SPNSDE_MOD
CONTAINS
SUBROUTINE SPNSDE(KM,KF_SCALARS,PEPSNM,PF,PNSD)

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_FIELDS
USE TPM_TRANS

#ifdef DOC

!**** *SPNSDE* - Compute North-South derivative in spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the North-south derivative

!**   Interface.
!     ----------
!        CALL SPNSDE(...)

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
!        Original : 00-02-01 From SPNSDE in IFS CY22R1

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER_M, INTENT(IN)  :: KM,KF_SCALARS
REAL_B,    INTENT(IN)  :: PEPSNM(0:R%NTMAX+2)
REAL_B,    INTENT(IN)  :: PF(:,:)
REAL_B,    INTENT(OUT) :: PNSD(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IJ, ISKIP, J, JN,II,JI,ISMAX
REAL_B :: ZEPSNM(-1:R%NSMAX+4)
REAL_B :: ZN(-1:R%NTMAX+4)


!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------


!*       1.1      COMPUTE

ISMAX = R%NSMAX
DO JN=KM-1,ISMAX+2
  IJ = ISMAX+3-JN
  ZN(IJ) = F%RN(JN)
  ZEPSNM(IJ) = PEPSNM(JN)
ENDDO
ZN(0) = F%RN(ISMAX+3)

IF(KM == 0) THEN
  ISKIP = 2
ELSE
  ISKIP = 1
ENDIF

DO J=1,2*KF_SCALARS,ISKIP
  DO JI=2,ISMAX+3-KM
    PNSD(JI,J) = -ZN(JI+1)*ZEPSNM(JI)*PF(JI+1,J)+&
     &ZN(JI-2)*ZEPSNM(JI-1)*PF(JI-1,J)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SPNSDE
END MODULE SPNSDE_MOD
