! (C) Copyright 1991- ECMWF.
! (C) Copyright 1991- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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

USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRBT

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

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

PU => ZOA1(IUS:IUE,:,:)
PV => ZOA1(IVS:IVE,:,:)
PVOR => ZOA2(IVORS:IVORE,:,:)
PDIV => ZOA2(IDIVS:IDIVE,:,:)

#ifdef ACCGPU
!$ACC DATA&
!$ACC& CREATE(ZN) &
!$ACC& COPY(D_MYMS,D_NUMP,R_NTMAX) &
!$ACC& COPY(F,F%RN,F%NLTN) &
!$ACC& PRESENT(ZEPSNM,PU,PV,PVOR,PDIV)
#endif
#ifdef OMPGPU
!WARNING: following line should be PRESENT,ALLOC but causes issues with AMD compiler!
!$OMP TARGET DATA&
!$OMP& MAP(ALLOC:ZN) &
!$OMP& MAP(TO:D_MYMS,D_NUMP,R_NTMAX) &
!$OMP& MAP(TO:F,F%RN,F%NLTN) &
!$OMP& MAP(ALLOC:ZEPSNM,PU,PV,PVOR,PDIV)
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!! DEFAULT(NONE) SHARED(R_NTMAX,ZN,F)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP DEFAULT(NONE) PRESENT(R_NTMAX,ZN,F)
#endif
DO J=-1,R_NTMAX+3
  ZN(j) = F%RN(j)
ENDDO
!*       1.1      SET N=KM-1 COMPONENT TO 0 FOR U AND V

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
!! PRIVATE(KM,IN) DEFAULT(NONE) &
!!$OMP& SHARED(D_NUMP,KFIELD,D_MYMS,F,PU,PV)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IN) DEFAULT(NONE) &
!$ACC& COPYIN(KFIELD) &
!$ACC& PRESENT(D_NUMP,D_MYMS,F,F%NLTN,PU,PV)
#endif
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

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(IR,II,IN,KM,ZKM) DEFAULT(NONE) &
!$OMP& SHARED(D_NUMP,R_NTMAX,KFIELD,D_MYMS,PVOR,PV,PU,ZN,PDIV,ZEPSNM)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IR,II,IN,KM,ZKM) DEFAULT(NONE) &
!$ACC& COPYIN(KFIELD) &
!$ACC& PRESENT(D_NUMP,R_NTMAX,D_MYMS,PVOR,PV,PU,ZN,PDIV,ZEPSNM)
#endif
DO KMLOC=1,D_NUMP
  DO JN=0,R_NTMAX
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1
      KM = D_MYMS(KMLOC)
      ZKM = REAL(KM,JPRBT)
      IN = R_NTMAX+2-JN

      IF(KM /= 0 .AND. JN.GE.KM) THEN
      PVOR(IR,IN,KMLOC) = -ZKM*PV(II,IN,KMLOC)-&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PU(IR,IN-1,KMLOC)+&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PU(IR,IN+1,KMLOC)
      PVOR(II,IN,KMLOC) = +ZKM*PV(IR,IN,KMLOC)-&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PU(II,IN-1,KMLOC)+&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PU(II,IN+1,KMLOC)
      PDIV(IR,IN,KMLOC) = -ZKM*PU(II,IN,KMLOC)+&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PV(IR,IN-1,KMLOC)-&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PV(IR,IN+1,KMLOC)
      PDIV(II,IN,KMLOC) = +ZKM*PU(IR,IN,KMLOC)+&
       &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PV(II,IN-1,KMLOC)-&
       &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PV(II,IN+1,KMLOC)
      ELSE
        IF(KM == 0) THEN
         PVOR(IR,IN,KMLOC) = -&
         &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PU(IR,IN-1,KMLOC)+&
         &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PU(IR,IN+1,KMLOC)
         PDIV(IR,IN,KMLOC) = &
         &ZN(JN)*ZEPSNM(KMLOC,JN+1)*PV(IR,IN-1,KMLOC)-&
         &ZN(JN+1)*ZEPSNM(KMLOC,JN)*PV(IR,IN+1,KMLOC)
        ENDIF
      ENDIF
   ENDDO
  ENDDO
ENDDO
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif
!     ------------------------------------------------------------------

END SUBROUTINE UVTVD
END MODULE UVTVD_MOD
