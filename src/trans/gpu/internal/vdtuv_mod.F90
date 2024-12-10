! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
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

USE PARKIND_ECTRANS, ONLY: JPIM, JPRB, JPRBT
USE TPM_DIM,         ONLY: R
USE TPM_FIELDS,      ONLY: F
USE TPM_DISTR,       ONLY: D

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

INTEGER(KIND=JPIM) :: KM, KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRBT), INTENT(IN)   :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
REAL(KIND=JPRB), INTENT(INOUT) :: PVOR(:,:,:),PDIV(:,:,:)
REAL(KIND=JPRB), INTENT(OUT)   :: PU  (:,:,:),PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, JI

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM

ASSOCIATE(D_NUMP=>D%NUMP, D_MYMS=>D%MYMS, R_NTMAX=>R%NTMAX, F_RLAPIN=>F%RLAPIN)

#ifdef ACCGPU
!$ACC DATA                                                       &
!$ACC&      PRESENT(R,R_NTMAX,D,D_MYMS,D_NUMP,F,F_RLAPIN) &
!$ACC&      PRESENT(PEPSNM, PVOR, PDIV)                          &
!$ACC&      PRESENT(PU, PV)
#endif

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
!$OMP& PRIVATE(IR,II,KM,ZKM,JI) SHARED(D,R,F,PEPSNM,PVOR,PDIV,PU,PV,KFIELD) &
!$OMP& MAP(TO:KFIELD)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) PRIVATE(IR,II,KM,ZKM,JI) FIRSTPRIVATE(KFIELD,KMLOC) &
#ifndef _CRAYFTN
!$ACC& ASYNC(1)
#else
!$ACC&
#endif
#endif
DO KMLOC=1,D_NUMP
  DO JN=0,R_NTMAX+1
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1
      KM = D_MYMS(KMLOC)
      ZKM = REAL(KM,JPRBT)

      IF(KM /= 0 .AND. JN >= KM) THEN
        ! (DO JN=KN,R_NTMAX)
        JI = R_NTMAX+3-JN
        PU(IR,JI,KMLOC) = -ZKM*F_RLAPIN(JN)*PDIV(II,JI,KMLOC)+&
         &(JN-1)*PEPSNM(KMLOC,JN)*F_RLAPIN(JN-1)*PVOR(IR,JI+1,KMLOC)-&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F_RLAPIN(JN+1)*PVOR(IR,JI-1,KMLOC)
        PU(II,JI,KMLOC) = +ZKM*F_RLAPIN(JN)*PDIV(IR,JI,KMLOC)+&
         &(JN-1)*PEPSNM(KMLOC,JN)*F_RLAPIN(JN-1)*PVOR(II,JI+1,KMLOC)-&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F_RLAPIN(JN+1)*PVOR(II,JI-1,KMLOC)
        PV(IR,JI,KMLOC) = -ZKM*F_RLAPIN(JN)*PVOR(II,JI,KMLOC)-&
         &(JN-1)*PEPSNM(KMLOC,JN)*F_RLAPIN(JN-1)*PDIV(IR,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F_RLAPIN(JN+1)*PDIV(IR,JI-1,KMLOC)
        PV(II,JI,KMLOC) = +ZKM*F_RLAPIN(JN)*PVOR(IR,JI,KMLOC)-&
         &(JN-1)*PEPSNM(KMLOC,JN)*F_RLAPIN(JN-1)*PDIV(II,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F_RLAPIN(JN+1)*PDIV(II,JI-1,KMLOC)

      ELSEIF(KM == 0) THEN
        ! (DO JN=0,R_NTMAX)
        JI = R_NTMAX+3-JN

        PU(IR,JI,KMLOC) = +&
         &(JN-1)*PEPSNM(KMLOC,JN)*F_RLAPIN(JN-1)*PVOR(IR,JI+1,KMLOC)-&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F_RLAPIN(JN+1)*PVOR(IR,JI-1,KMLOC)
        PV(IR,JI,KMLOC) = -&
         &(JN-1)*PEPSNM(KMLOC,JN)*F_RLAPIN(JN-1)*PDIV(IR,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F_RLAPIN(JN+1)*PDIV(IR,JI-1,KMLOC)
      ENDIF
    ENDDO
  ENDDO
ENDDO

#ifdef ACCGPU
!$ACC END DATA
#endif
!     ------------------------------------------------------------------
END ASSOCIATE

END SUBROUTINE VDTUV
END MODULE VDTUV_MOD

