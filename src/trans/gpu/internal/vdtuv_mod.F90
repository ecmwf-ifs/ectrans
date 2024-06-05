! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VDTUV_MOD
CONTAINS
SUBROUTINE VDTUV(KFIELD)

USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRB,  JPRBT

USE TPM_DIM         ,ONLY : R, R_NSMAX
USE TPM_FIELDS      ,ONLY : F_RN, F_RLAPIN, ZEPSNM, ZIA
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_GEN         ,ONLY : NOUT


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

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KM, kmloc
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
!REAL(KIND=JPRBT), INTENT(IN)    :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
!REAL(KIND=JPRB), INTENT(IN)    :: PVOR(:,:,:),PDIV(:,:,:)
!REAL(KIND=JPRB), INTENT(OUT)   :: PU  (:,:,:),PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, ISMAX,JI

INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU
INTEGER(KIND=JPIM) :: ZKM

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZN(-1:R%NTMAX+4)
REAL(KIND=JPRBT) :: ZLAPIN(-1:R%NSMAX+4)
REAL(KIND=JPRBT) :: ZZEPSNM(-1:R%NSMAX+4)
REAL(KIND=JPRBT), POINTER :: PU(:,:,:),PV(:,:,:),PVOR(:,:,:),PDIV(:,:,:)

IVORL = 1
IVORU = 2*KFIELD
IDIVL = 2*KFIELD+1
IDIVU = 4*KFIELD
IUL   = 4*KFIELD+1
IUU   = 6*KFIELD
IVL   = 6*KFIELD+1
IVU   = 8*KFIELD

PU => ZIA(IUL:IUU,:,:)
PV => ZIA(IVL:IVU,:,:)
PVOR => ZIA(IVORL:IVORU,:,:)
PDIV => ZIA(IDIVL:IDIVU,:,:)

#ifdef ACCGPU
!$ACC DATA                                     &
!$ACC&      CREATE (ZZEPSNM, ZN, ZLAPIN)         &
!$ACC&      COPYIN (R_NSMAX, D_MYMS,D_NUMP,F_RLAPIN,F_RN)   &
!$ACC&      PRESENT(ZEPSNM, PVOR, PDIV)         &
!$ACC&      PRESENT(PU, PV)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA                              &
!$OMP&      MAP (PRESENT,ALLOC:ZEPSNM, ZN, ZLAPIN)     &
!$OMP&      MAP (TO:R_NSMAX, D_MYMS,D_NUMP,F_RLAPIN,F_RN)  &
!$OMP&      MAP(PRESENT,ALLOC:ZEPSNM, PVOR, PDIV)      &
!$OMP&      MAP(PRESENT,ALLOC:PU, PV)
#endif

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ISMAX = R_NSMAX
DO KMLOC=1,D_NUMP
  ZKM = D_MYMS(KMLOC)
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) PRIVATE(IJ) &
  !$OMP&    SHARED(ZKM,ISMAX,ZN,ZLAPIN,ZZEPSNM,ZEPSNM,KMLOC,F_RN,F_RLAPIN)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IJ) &
  !$ACC&    PRESENT(F_RN,F_RLAPIN,ZN,ZLAPIN,ZZEPSNM,ZEPSNM) &
  !$ACC&    COPYIN(ISMAX,ZKM,KMLOC)
#endif
  DO JN=ZKM-1,ISMAX+2
    IJ = ISMAX+3-JN
    ZN(IJ) = F_RN(JN)
    ZLAPIN(IJ) = F_RLAPIN(JN)
    ZN(0) = F_RN(ISMAX+3)
    IF( JN >= 0 ) THEN
        ZZEPSNM(IJ) = ZEPSNM(KMLOC,JN)
    ELSE
        ZZEPSNM(IJ) = 0._JPRBT
    ENDIF
  ENDDO
  !!! alternative !$ACC UPDATE DEVICE(ZN, ZLAPIN, ZZEPSNM)

!*       1.1      U AND V (KM=0) .

IF(ZKM == 0) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) PRIVATE(IR) &
  !$OMP&   SHARED(KFIELD,ISMAX,KMLOC,PU,ZN,ZLAPIN,PVOR,PV,PDIV,ZZEPSNM)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IR) &
  !$ACC&   COPYIN(ZN,ZLAPIN,ZZEPSNM) &
  !$ACC&   PRESENT(PU,PVOR,PV,PDIV) &
  !$ACC&   COPYIN(KFIELD,ISMAX,KMLOC)
#endif
  DO J=1,KFIELD
    IR = 2*J-1
    DO JI=2,ISMAX+3
      PU(IR,JI,KMLOC) = +&
      &ZN(JI+1)*ZZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(IR,JI+1,KMLOC)-&
      &ZN(JI-2)*ZZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(IR,JI-1,KMLOC)
      PV(IR,JI,KMLOC) = -&
      &ZN(JI+1)*ZZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(IR,JI+1,KMLOC)+&
      &ZN(JI-2)*ZZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(IR,JI-1,KMLOC)
    ENDDO
  ENDDO
ELSE
!*       1.2      U AND V (KM!=0) .

#ifdef OMPGPU
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) PRIVATE(IR,II) &
    !$OMP&   SHARED(KFIELD,ZKM,ISMAX,KMLOC,PU,ZN,ZLAPIN,PVOR,PV,PDIV,ZZEPSNM)
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(IR,II) &
    !$ACC&   PRESENT(PU,PVOR,PV,PDIV) &
    !$ACC&   COPYIN(ZN,ZLAPIN,ZZEPSNM) &
    !$ACC&   COPYIN(KFIELD,ZKM,ISMAX,KMLOC)
#endif
    DO J=1,KFIELD
      DO JI=2,ISMAX+3-ZKM
        !ZKM = D_MYMS(KMLOC)
        IR = 2*J-1
        II = IR+1
        !IF (ZKM>0 .AND. JI<=ISMAX+3-zKM) THEN
          PU(IR,JI,KMLOC) = -ZKM*ZLAPIN(JI)*PDIV(II,JI,KMLOC)+&
          &ZN(JI+1)*ZZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(IR,JI+1,KMLOC)-&
          &ZN(JI-2)*ZZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(IR,JI-1,KMLOC)
          PU(II,JI,KMLOC) = +ZKM*ZLAPIN(JI)*PDIV(IR,JI,KMLOC)+&
          &ZN(JI+1)*ZZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(II,JI+1,KMLOC)-&
          &ZN(JI-2)*ZZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(II,JI-1,KMLOC)
          PV(IR,JI,KMLOC) = -ZKM*ZLAPIN(JI)*PVOR(II,JI,KMLOC)-&
          &ZN(JI+1)*ZZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(IR,JI+1,KMLOC)+&
          &ZN(JI-2)*ZZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(IR,JI-1,KMLOC)
          PV(II,JI,KMLOC) = +ZKM*ZLAPIN(JI)*PVOR(IR,JI,KMLOC)-&
          &ZN(JI+1)*ZZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(II,JI+1,KMLOC)+&
          &ZN(JI-2)*ZZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(II,JI-1,KMLOC)
        !ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDDO

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif
!     ------------------------------------------------------------------

END SUBROUTINE VDTUV
END MODULE VDTUV_MOD

