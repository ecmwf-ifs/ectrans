MODULE VDTUV_MOD
CONTAINS
SUBROUTINE VDTUV(KM,KFIELD,PEPSNM,PVOR,PDIV,PU,PV)

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_FIELDS

#ifdef DOC

!**** *VDTUV* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL VDTUV(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!                              PVOR(NLEI1,2*KFIELD) - vorticity (input)
!                              PDIV(NLEI1,2*KFIELD) - divergence (input)
!                              PU(NLEI1,2*KFIELD)   - u wind (output)
!                              PV(NLEI1,2*KFIELD)   - v wind (output)
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

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

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
!        Original : 00-02-01 From VDTUV in IFS CY22R1

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER_M, INTENT(IN) :: KM,KFIELD
REAL_B, INTENT(IN)    :: PEPSNM(0:R%NTMAX+2)
REAL_B, INTENT(IN)    :: PVOR(:,:),PDIV(:,:)
REAL_B, INTENT(OUT)   :: PU  (:,:),PV  (:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: II, IJ, IR, J, JN, ISMAX

!     LOCAL REAL SCALARS
REAL_B :: ZKM
REAL_B :: ZN(-1:R%NTMAX+3)
REAL_B :: ZLAPIN(-1:R%NSMAX+2)



!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ZKM = KM
ISMAX = R%NSMAX
DO JN=KM-1,ISMAX+2
  ZN(JN) = F%RN(JN)
  ZLAPIN(JN) = F%RLAPIN(JN)
ENDDO
ZN(ISMAX+3) = F%RN(ISMAX+3)

!*       1.1      U AND V (KM=0) .

IF(KM == 0) THEN
  DO JN=KM,ISMAX+1
    IJ = ISMAX+3-JN
    DO J=1,KFIELD
      IR = 2*J-1
      PU(IJ,IR) = +&
       &ZN(JN-1)*PEPSNM(JN)*ZLAPIN(JN-1)*PVOR(IJ+1,IR)-&
       &ZN(JN+2)*PEPSNM(JN+1)*ZLAPIN(JN+1)*PVOR(IJ-1,IR)
      PV(IJ,IR) = -&
       &ZN(JN-1)*PEPSNM(JN)*ZLAPIN(JN-1)*PDIV(IJ+1,IR)+&
       &ZN(JN+2)*PEPSNM(JN+1)*ZLAPIN(JN+1)*PDIV(IJ-1,IR)
    ENDDO
  ENDDO

!*       1.2      U AND V (KM!=0) .

ELSE
  DO JN=KM,ISMAX+1
    IJ = ISMAX+3-JN
!DIR$ IVDEP
!OCL NOVREC
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1
      PU(IJ,IR) = -ZKM*ZLAPIN(JN)*PDIV(IJ,II)+&
       &ZN(JN-1)*PEPSNM(JN)*ZLAPIN(JN-1)*PVOR(IJ+1,IR)-&
       &ZN(JN+2)*PEPSNM(JN+1)*ZLAPIN(JN+1)*PVOR(IJ-1,IR)
      PU(IJ,II) = +ZKM*ZLAPIN(JN)*PDIV(IJ,IR)+&
       &ZN(JN-1)*PEPSNM(JN)*ZLAPIN(JN-1)*PVOR(IJ+1,II)-&
       &ZN(JN+2)*PEPSNM(JN+1)*ZLAPIN(JN+1)*PVOR(IJ-1,II)
      PV(IJ,IR) = -ZKM*ZLAPIN(JN)*PVOR(IJ,II)-&
       &ZN(JN-1)*PEPSNM(JN)*ZLAPIN(JN-1)*PDIV(IJ+1,IR)+&
       &ZN(JN+2)*PEPSNM(JN+1)*ZLAPIN(JN+1)*PDIV(IJ-1,IR)
      PV(IJ,II) = +ZKM*ZLAPIN(JN)*PVOR(IJ,IR)-&
       &ZN(JN-1)*PEPSNM(JN)*ZLAPIN(JN-1)*PDIV(IJ+1,II)+&
       &ZN(JN+2)*PEPSNM(JN+1)*ZLAPIN(JN+1)*PDIV(IJ-1,II)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE VDTUV
END MODULE VDTUV_MOD
