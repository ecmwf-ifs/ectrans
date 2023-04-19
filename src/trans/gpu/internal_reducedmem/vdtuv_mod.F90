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
SUBROUTINE VDTUV(KFIELD,PEPSNM,PVOR,PDIV,PU,PV)

USE PARKIND_ECTRANS ,ONLY : JPIM, JPRB, JPRBT

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
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
INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT), INTENT(IN)    :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
REAL(KIND=JPRB),  INTENT(IN)    :: PVOR(:,:,:)
REAL(KIND=JPRB),  INTENT(IN)    :: PDIV(:,:,:)
REAL(KIND=JPRB),  INTENT(OUT)   :: PU  (:,:,:)
REAL(KIND=JPRB),  INTENT(OUT)   :: PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, ISMAX,JI

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM
REAL(KIND=JPRBT) :: ZZN(-1:R%NTMAX+4)
REAL(KIND=JPRBT) :: ZZLAPIN(-1:R%NSMAX+4)
REAL(KIND=JPRBT) :: ZZEPSNM(-1:R%NSMAX+4)

#ifdef ACCGPU
!$ACC DATA                                     &
!$ACC      CREATE (ZZEPSNM, ZZN, ZZLAPIN)      &
!$ACC      COPYIN(PEPSNM, PVOR, PDIV)          &
!$ACC      COPYIN (D,D%MYMS,F,F%RLAPIN,F%RN)   &
!$ACC      COPYOUT(PU, PV)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA                                     &
!$OMP&     MAP(ALLOC:ZZEPSNM, ZZN, ZZLAPIN)      &
!$OMP&     MAP(TO:PEPSNM, PVOR, PDIV)          &
!$OMP&     MAP(TO:D,D%MYMS,F,F%RLAPIN,F%RN)   &
!$OMP&     MAP(FROM:PU, PV)
#endif

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ISMAX = R%NSMAX
DO KMLOC=1,D%NUMP
  ZKM = D%MYMS(KMLOC)
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP
#endif
  DO JN=ZKM-1,ISMAX+2
    IJ = ISMAX+3-JN
    ZZN(IJ) = F%RN(JN)
    ZZLAPIN(IJ) = F%RLAPIN(JN)
    IF( JN >= 0 ) ZZEPSNM(IJ) = PEPSNM(KMLOC,JN)
  ENDDO
  ZZN(0) = F%RN(ISMAX+3)

!*       1.1      U AND V (KM=0) .

IF(ZKM == 0) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(IR)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IR)
#endif
  DO J=1,KFIELD
    DO JI=2,ISMAX+3
      IR = 2*J-1
      PU(IR,JI,KMLOC) = +&
      &ZZN(JI+1)*ZZEPSNM(JI)*ZZLAPIN(JI+1)*PVOR(IR,JI+1,KMLOC)-&
      &ZZN(JI-2)*ZZEPSNM(JI-1)*ZZLAPIN(JI-1)*PVOR(IR,JI-1,KMLOC)
      PV(IR,JI,KMLOC) = -&
      &ZZN(JI+1)*ZZEPSNM(JI)*ZZLAPIN(JI+1)*PDIV(IR,JI+1,KMLOC)+&
      &ZZN(JI-2)*ZZEPSNM(JI-1)*ZZLAPIN(JI-1)*PDIV(IR,JI-1,KMLOC)
    ENDDO
  ENDDO
ELSE
!*       1.2      U AND V (KM!=0) .

#ifdef OMPGPU
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(IR,II)
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IR,II)
#endif
    DO J=1,KFIELD
      DO JI=2,ISMAX+3-ZKM
        !ZKM = D_MYMS(KMLOC)
        IR = 2*J-1
        II = IR+1
        !IF (ZKM>0 .AND. JI<=ISMAX+3-ZKM) THEN
          PU(IR,JI,KMLOC) = -ZKM*ZZLAPIN(JI)*PDIV(II,JI,KMLOC)+&
          &ZZN(JI+1)*ZZEPSNM(JI)*ZZLAPIN(JI+1)*PVOR(IR,JI+1,KMLOC)-&
          &ZZN(JI-2)*ZZEPSNM(JI-1)*ZZLAPIN(JI-1)*PVOR(IR,JI-1,KMLOC)
          PU(II,JI,KMLOC) = +ZKM*ZZLAPIN(JI)*PDIV(IR,JI,KMLOC)+&
          &ZZN(JI+1)*ZZEPSNM(JI)*ZZLAPIN(JI+1)*PVOR(II,JI+1,KMLOC)-&
          &ZZN(JI-2)*ZZEPSNM(JI-1)*ZZLAPIN(JI-1)*PVOR(II,JI-1,KMLOC)
          PV(IR,JI,KMLOC) = -ZKM*ZZLAPIN(JI)*PVOR(II,JI,KMLOC)-&
          &ZZN(JI+1)*ZZEPSNM(JI)*ZZLAPIN(JI+1)*PDIV(IR,JI+1,KMLOC)+&
          &ZZN(JI-2)*ZZEPSNM(JI-1)*ZZLAPIN(JI-1)*PDIV(IR,JI-1,KMLOC)
          PV(II,JI,KMLOC) = +ZKM*ZZLAPIN(JI)*PVOR(IR,JI,KMLOC)-&
          &ZZN(JI+1)*ZZEPSNM(JI)*ZZLAPIN(JI+1)*PDIV(II,JI+1,KMLOC)+&
          &ZZN(JI-2)*ZZEPSNM(JI-1)*ZZLAPIN(JI-1)*PDIV(II,JI-1,KMLOC)
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
