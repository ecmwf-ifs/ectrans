MODULE UVTVD_MOD
CONTAINS
SUBROUTINE UVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!**** *UVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX.

!**   Interface.
!     ----------
!        CALL UVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KFIELD - number of fields (levels)
!                              PEPSNM - REPSNM for wavenumber KM
!                              PU - u wind component for zonal
!                                   wavenumber KM
!                              PV - v wind component for zonal
!                                   wavenumber KM
!                              PVOR - vorticity for zonal
!                                     wavenumber KM
!                              PDIV - divergence for zonal
!                                     wavenumber KM


!     Method.  See ref.
!     -------

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 91-07-01
!        D. Giard : NTMAX instead of NSMAX
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_FIELDS
USE TPM_DISTR

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER_M, INTENT(IN)  :: KFIELD
INTEGER_M, INTENT(IN)  :: KM

REAL_B, INTENT(IN)     :: PEPSNM(0:R%NTMAX+2)
REAL_B, INTENT(OUT)    :: PVOR(:,:),PDIV(:,:)
REAL_B, INTENT(INOUT)  :: PU  (:,:),PV  (:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: II, IN, IR, J, JN

!     LOCAL REAL SCALARS
REAL_B :: ZKM


!     ------------------------------------------------------------------


!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ZKM = KM

!*       1.1      SET N=KM-1 COMPONENT TO 0 FOR U AND V

IN = F%NLTN(KM-1)
DO J=1,2*KFIELD
  PU(IN,J) = _ZERO_
  PV(IN,J) = _ZERO_
ENDDO

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

IF(KM /= 0) THEN
  DO JN=KM,R%NTMAX
    IN=F%NLTN(JN)
!DIR$ IVDEP
!OCL NOVREC
    DO J=1,KFIELD
      IR=2*J-1
      II=IR+1
      PVOR(IN,IR)=-ZKM*PV(IN,II)-&
       &F%RN(JN)*PEPSNM(JN+1)*PU(IN-1,IR)+&
       &F%RN(JN+1)*PEPSNM(JN)*PU(IN+1,IR)
      PVOR(IN,II)=+ZKM*PV(IN,IR)-&
       &F%RN(JN)*PEPSNM(JN+1)*PU(IN-1,II)+&
       &F%RN(JN+1)*PEPSNM(JN)*PU(IN+1,II)
      PDIV(IN,IR)=-ZKM*PU(IN,II)+&
       &F%RN(JN)*PEPSNM(JN+1)*PV(IN-1,IR)-&
       &F%RN(JN+1)*PEPSNM(JN)*PV(IN+1,IR)
      PDIV(IN,II)=+ZKM*PU(IN,IR)+&
       &F%RN(JN)*PEPSNM(JN+1)*PV(IN-1,II)-&
       &F%RN(JN+1)*PEPSNM(JN)*PV(IN+1,II)
    ENDDO
  ENDDO
ELSE
  DO JN=KM,R%NTMAX
    IN=F%NLTN(JN)
    DO J=1,KFIELD
      IR=2*J-1
      PVOR(IN,IR)=-&
       &F%RN(JN)*PEPSNM(JN+1)*PU(IN-1,IR)+&
       &F%RN(JN+1)*PEPSNM(JN)*PU(IN+1,IR)
      PDIV(IN,IR)=&
       &F%RN(JN)*PEPSNM(JN+1)*PV(IN-1,IR)-&
       &F%RN(JN+1)*PEPSNM(JN)*PV(IN+1,IR)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UVTVD
END MODULE UVTVD_MOD
