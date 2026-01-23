! (C) Copyright 1991- ECMWF.
! (C) Copyright 1991- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UVTVDAD_MOD
CONTAINS
SUBROUTINE UVTVDAD(KF_UV,PU,PV,PVOR,PDIV)

!**** *UVTVDAD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX.

!**   Interface.
!     ----------
!        CALL UVTVDAD(KM,KF_UV,PEPSNM,PU,PV,PVOR,PDIV)

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KF_UV - number of fields (levels)
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

USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT
USE TPM_DIM,         ONLY: R
USE TPM_DISTR,       ONLY: D
USE TPM_FIELDS_GPU,  ONLY: FG
!

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_UV
REAL(KIND=JPRBT), INTENT(INOUT)    :: PVOR(:,:,:),PDIV(:,:,:)
REAL(KIND=JPRBT), INTENT(INOUT) :: PU  (:,:,:),PV  (:,:,:)
INTEGER(KIND=JPIM)  :: KM, KMLOC

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IN, IR, J, JN

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM,ZJN

!     ------------------------------------------------------------------
ASSOCIATE(D_NUMP=>D%NUMP, R_NTMAX=>R%NTMAX, D_MYMS=>D%MYMS, ZEPSNM=>FG%ZEPSNM)

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

#ifdef ACCGPU
!$ACC DATA &
!$ACC& PRESENT(D,D_MYMS,D_NUMP,R,R_NTMAX) &
!$ACC& PRESENT(FG,ZEPSNM,PU,PV,PVOR,PDIV) ASYNC(1)
#endif

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(IR,II,IN,KM,ZKM,ZJN) &
!$OMP& SHARED(D,R,KF_UV,FG,PVOR,PV,PU,PDIV) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IR,II,IN,KM,ZKM,ZJN) FIRSTPRIVATE(KF_UV) DEFAULT(NONE) &
#ifndef _CRAYFTN
!$ACC& ASYNC(1)
#else
!$ACC&
#endif
#endif
DO KMLOC=1,D_NUMP
  DO JN=-1,R_NTMAX+1
    DO J=1,KF_UV
      IR = 2*J-1
      II = IR+1
      KM = D_MYMS(KMLOC)
      ZKM = REAL(KM,JPRBT)

      IF(KM /= 0 .AND. JN >= (KM-1)) THEN
        ! (DO JN=KN,R_NTMAX)
        IN = R_NTMAX+3-JN
        ZJN = JN

        PU(IR,IN,KMLOC) = 0
        PU(II,IN,KMLOC) = 0
        PV(IR,IN,KMLOC) = 0
        PV(II,IN,KMLOC) = 0

        IF (2 <= IN .AND. IN <= R_NTMAX + 2) THEN
          PU(IR,IN,KMLOC) = PU(IR,IN,KMLOC) - (ZJN-1)*ZEPSNM(KMLOC,JN)*PVOR(IR,IN+1,KMLOC)
          PU(II,IN,KMLOC) = PU(II,IN,KMLOC) - (ZJN-1)*ZEPSNM(KMLOC,JN)*PVOR(II,IN+1,KMLOC)
          PV(IR,IN,KMLOC) = PV(IR,IN,KMLOC) + (ZJN-1)*ZEPSNM(KMLOC,JN)*PDIV(IR,IN+1,KMLOC)
          PV(II,IN,KMLOC) = PV(II,IN,KMLOC) + (ZJN-1)*ZEPSNM(KMLOC,JN)*PDIV(II,IN+1,KMLOC)
        ENDIF
        IF (3 <= IN .AND. IN <= R_NTMAX + 3) THEN
          PU(IR,IN,KMLOC) = PU(IR,IN,KMLOC) + ZKM*PDIV(II,IN,KMLOC)
          PU(II,IN,KMLOC) = PU(II,IN,KMLOC) - ZKM*PDIV(IR,IN,KMLOC)
          PV(IR,IN,KMLOC) = PV(IR,IN,KMLOC) + ZKM*PVOR(II,IN,KMLOC)
          PV(II,IN,KMLOC) = PV(II,IN,KMLOC) - ZKM*PVOR(IR,IN,KMLOC)
        ENDIF
        IF (4 <= IN .AND. IN <= R_NTMAX + 4) THEN
          PU(IR,IN,KMLOC) = PU(IR,IN,KMLOC) + (ZJN+2)*ZEPSNM(KMLOC,JN+1)*PVOR(IR,IN-1,KMLOC)
          PU(II,IN,KMLOC) = PU(II,IN,KMLOC) + (ZJN+2)*ZEPSNM(KMLOC,JN+1)*PVOR(II,IN-1,KMLOC)
          PV(IR,IN,KMLOC) = PV(IR,IN,KMLOC) - (ZJN+2)*ZEPSNM(KMLOC,JN+1)*PDIV(IR,IN-1,KMLOC)
          PV(II,IN,KMLOC) = PV(II,IN,KMLOC) - (ZJN+2)*ZEPSNM(KMLOC,JN+1)*PDIV(II,IN-1,KMLOC)
        ENDIF

      ELSEIF(KM == 0) THEN
        ! (DO JN=0,R_NTMAX)
        IN = R_NTMAX+3-JN
        ZJN = JN

        PU(IR,IN,KMLOC) = 0
        PU(II,IN,KMLOC) = 0
        PV(IR,IN,KMLOC) = 0
        PV(II,IN,KMLOC) = 0
        
        IF (2 <= IN .AND. IN <= R_NTMAX + 2) THEN
          PU(IR,IN,KMLOC) = PU(IR,IN,KMLOC) - (ZJN-1)*ZEPSNM(KMLOC,JN)*PVOR(IR,IN+1,KMLOC)
          PV(IR,IN,KMLOC) = PV(IR,IN,KMLOC) + (ZJN-1)*ZEPSNM(KMLOC,JN)*PDIV(IR,IN+1,KMLOC)
        ENDIF
        IF (4 <= IN .AND. IN <= R_NTMAX + 4) THEN
          PU(IR,IN,KMLOC) = PU(IR,IN,KMLOC) + (ZJN+2)*ZEPSNM(KMLOC,JN+1)*PVOR(IR,IN-1,KMLOC)
          PV(IR,IN,KMLOC) = PV(IR,IN,KMLOC) - (ZJN+2)*ZEPSNM(KMLOC,JN+1)*PDIV(IR,IN-1,KMLOC)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO
#ifdef ACCGPU
!$ACC END DATA
#endif
!     ------------------------------------------------------------------
END ASSOCIATE

END SUBROUTINE UVTVDAD
END MODULE UVTVDAD_MOD
