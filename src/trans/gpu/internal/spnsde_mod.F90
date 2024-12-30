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

MODULE SPNSDE_MOD
CONTAINS
SUBROUTINE SPNSDE(KF_SCALARS,PEPSNM,PF,PNSD)

USE PARKIND_ECTRANS, ONLY: JPIM, JPRB, JPRBT
USE TPM_DIM,         ONLY: R
USE TPM_DISTR,       ONLY: D

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

IMPLICIT NONE

INTEGER(KIND=JPIM)  :: KM, KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_SCALARS
!REAL(KIND=JPRBT),    INTENT(IN)  :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRBT),    INTENT(IN)  :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
REAL(KIND=JPRB),    INTENT(IN)  :: PF(:,:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PNSD(:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IJ, ISKIP, J, JN, JI, IR, II

ASSOCIATE(D_NUMP=>D%NUMP, R_NTMAX=>R%NTMAX, D_MYMS=>D%MYMS)

#ifdef ACCGPU
!$ACC DATA                                  &
!$ACC&      PRESENT (R,R_NTMAX, D,D_MYMS)       &
!$ACC&      PRESENT (D_NUMP,PEPSNM, PF, PNSD) ASYNC(1)
#endif

!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------


!*       1.1      COMPUTE

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,IR,II,JI) MAP(TO:KF_SCALARS) &
!$OMP& SHARED(D,R)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(KM,IR,II,JI) FIRSTPRIVATE(KF_SCALARS) &
#ifndef _CRAYFTN
!$ACC& ASYNC(1)
#else
!$ACC&
#endif
#endif
DO KMLOC=1,D_NUMP
  DO JN=0,R_NTMAX+1
    DO J=1,KF_SCALARS
      IR = 2*J-1
      II = IR+1
      KM = D_MYMS(KMLOC)

      IF(KM /= 0 .AND. JN >= KM) THEN
        ! (DO JN=KN,R_NTMAX+1)
        JI = R_NTMAX+3-JN
        PNSD(IR,JI,KMLOC) = -(JN-1)*PEPSNM(KMLOC,JN)*PF(IR,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*PF(IR,JI-1,KMLOC)
        PNSD(II,JI,KMLOC) = -(JN-1)*PEPSNM(KMLOC,JN)*PF(II,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*PF(II,JI-1,KMLOC)

      ELSEIF(KM == 0) THEN
        ! (DO JN=0,R_NTMAX+1)
        JI = R_NTMAX+3-JN
        PNSD(IR,JI,KMLOC) = -(JN-1)*PEPSNM(KMLOC,JN)*PF(IR,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*PF(IR,JI-1,KMLOC)
      ENDIF
    ENDDO
  ENDDO
END DO

#ifdef ACCGPU
!$ACC END DATA
#endif

!     ------------------------------------------------------------------
END ASSOCIATE

END SUBROUTINE SPNSDE
END MODULE SPNSDE_MOD
