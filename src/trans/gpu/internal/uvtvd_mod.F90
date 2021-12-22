MODULE UVTVD_MOD
CONTAINS
SUBROUTINE UVTVD(KFIELD)
!SUBROUTINE UVTVD(KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

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

USE PARKIND1  ,ONLY : JPIM     ,JPRBT

USE TPM_DIM         ,ONLY : R, R_NTMAX
USE TPM_FIELDS      ,ONLY : F
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_FIELDS      ,ONLY : ZOA1,ZOA2,ZEPSNM
!

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELD
INTEGER(KIND=JPIM)  :: KM, KMLOC

!REAL(KIND=JPRBT), INTENT(IN)     :: PEPSNM(1:d%nump,0:R%NTMAX+2)
!REAL(KIND=JPRBT), INTENT(OUT)    :: PVOR(:,:,:),PDIV(:,:,:)
!REAL(KIND=JPRBT), INTENT(INOUT)  :: PU  (:,:,:),PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IN, IR, J, JN, ITMAX
INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM
REAL(KIND=JPRBT) :: ZN(-1:R%NTMAX+3)
REAL(KIND=JPRBT), POINTER :: PU(:,:,:),PV(:,:,:),PVOR(:,:,:),PDIV(:,:,:)

IUS = 1
IUE = 2*KFIELD
IVS = 2*KFIELD+1
IVE = 4*KFIELD
IVORS = 1
IVORE = 2*KFIELD
IDIVS = 2*KFIELD+1
IDIVE = 4*KFIELD

!$acc data&
!$acc& create(ZN) &
!$acc& copy(D_MYMS,D_NUMP,R_NTMAX) &
!$acc& copy(F,F%RN,F%NLTN) &
!$acc& present(ZEPSNM,ZOA1,ZOA2)

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

PU => ZOA1(IUS:IUE,:,:)
PV => ZOA1(IVS:IVE,:,:)
PVOR => ZOA2(IVORS:IVORE,:,:)
PDIV => ZOA2(IDIVS:IDIVE,:,:)

!$acc parallel loop
DO J=-1,R_NTMAX+3
  ZN(j) = F%RN(j)
ENDDO
!*       1.1      SET N=KM-1 COMPONENT TO 0 FOR U AND V

!$ACC parallel loop collapse(2) private(KM,IN)
DO KMLOC=1,D_NUMP
  DO J=1,2*KFIELD
    KM = D_MYMS(KMLOC)
    IN = F%NLTN(KM-1)
!    IN=R_NTMAX+3-KM
    PU(J,IN,KMLOC) = 0.0_JPRBT
    PV(J,IN,KMLOC) = 0.0_JPRBT
  ENDDO
ENDDO

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

!$ACC parallel loop collapse(3) private(IR,II,IN,KM,ZKM)
DO KMLOC=1,D_NUMP
  DO JN=0,R_NTMAX
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1
      KM = D_MYMS(KMLOC)
      ZKM = REAL(KM,JPRBT)
      IN = R_NTMAX+2-JN

      IF(KM /= 0 .and. JN.GE.KM) THEN
      PVOR(IR,IN,kmloc) = -ZKM*PV(II,IN,kmloc)-&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PU(IR,IN-1,kmloc)+&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PU(IR,IN+1,kmloc)
      PVOR(II,IN,kmloc) = +ZKM*PV(IR,IN,kmloc)-&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PU(II,IN-1,kmloc)+&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PU(II,IN+1,kmloc)
      PDIV(IR,IN,kmloc) = -ZKM*PU(II,IN,kmloc)+&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PV(IR,IN-1,kmloc)-&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PV(IR,IN+1,kmloc)
      PDIV(II,IN,kmloc) = +ZKM*PU(IR,IN,kmloc)+&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PV(II,IN-1,kmloc)-&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PV(II,IN+1,kmloc)
      ELSE
        IF(KM == 0) THEN
         PVOR(IR,IN,kmloc) = -&
         &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PU(IR,IN-1,kmloc)+&
         &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PU(IR,IN+1,kmloc)
         PDIV(IR,IN,kmloc) = &
         &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PV(IR,IN-1,kmloc)-&
         &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PV(IR,IN+1,kmloc)
        ENDIF
      ENDIF
   ENDDO
  ENDDO
ENDDO
!$acc end data
!     ------------------------------------------------------------------

END SUBROUTINE UVTVD
END MODULE UVTVD_MOD
