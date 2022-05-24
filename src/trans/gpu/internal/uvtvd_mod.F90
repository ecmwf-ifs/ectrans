! (C) Copyright 1991- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UVTVD_MOD
CONTAINS
SUBROUTINE UVTVD(KF_UV,PU,PV,PVOR,PDIV)

!**** *UVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX.

!**   Interface.
!     ----------
!        CALL UVTVD(KM,KF_UV,PEPSNM,PU,PV,PVOR,PDIV)

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

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT

USE TPM_DIM         ,ONLY : R, R_NTMAX
USE TPM_FIELDS      ,ONLY : F
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_FIELDS      ,ONLY : ZEPSNM
!

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_UV
INTEGER(KIND=JPIM)  :: KM, KMLOC

REAL(KIND=JPRBT), INTENT(OUT)    :: PVOR(:,:,:),PDIV(:,:,:)
REAL(KIND=JPRBT), INTENT(INOUT)  :: PU  (:,:,:),PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IN, IR, J, JN, ITMAX

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM,ZJN

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

!$ACC DATA&
!$ACC& PRESENT(D_MYMS,D_NUMP,R_NTMAX) &
!$ACC& PRESENT(F,F%RN,F%NLTN) &
!$ACC& PRESENT(ZEPSNM,PU,PV,PVOR,PDIV)

!*       1.1      SET N=KM-1 COMPONENT TO 0 FOR U AND V

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO J=1,2*KF_UV
    KM = D_MYMS(KMLOC)
    PU(J,R_NTMAX+4-KM,KMLOC) = 0.0_JPRBT
    PV(J,R_NTMAX+4-KM,KMLOC) = 0.0_JPRBT
  ENDDO
ENDDO

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IR,II,IN,KM,ZKM,JN,ZJN) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO J=1,KF_UV
    IR = 2*J-1
    II = IR+1
    KM = D_MYMS(KMLOC)
    ZKM = REAL(KM,JPRBT)

    IF(KM /= 0) THEN
      !$ACC LOOP SEQ
      DO JN=KM,R_NTMAX
        IN = R_NTMAX+3-JN
        ZJN = JN

        PVOR(IR,IN,kmloc) = -ZKM*PV(II,IN,kmloc)-&
         &ZJN*ZEPSNM(KMLOC,JN+1)*PU(IR,IN-1,kmloc)+&
         &(ZJN+1)*ZEPSNM(KMLOC,JN)*PU(IR,IN+1,kmloc)
        PVOR(II,IN,kmloc) = +ZKM*PV(IR,IN,kmloc)-&
         &ZJN*ZEPSNM(KMLOC,JN+1)*PU(II,IN-1,kmloc)+&
         &(ZJN+1)*ZEPSNM(KMLOC,JN)*PU(II,IN+1,kmloc)
        PDIV(IR,IN,kmloc) = -ZKM*PU(II,IN,kmloc)+&
         &ZJN*ZEPSNM(KMLOC,JN+1)*PV(IR,IN-1,kmloc)-&
         &(ZJN+1)*ZEPSNM(KMLOC,JN)*PV(IR,IN+1,kmloc)
        PDIV(II,IN,kmloc) = +ZKM*PU(IR,IN,kmloc)+&
         &ZJN*ZEPSNM(KMLOC,JN+1)*PV(II,IN-1,kmloc)-&
         &(ZJN+1)*ZEPSNM(KMLOC,JN)*PV(II,IN+1,kmloc)
      ENDDO
    ELSE
      !$ACC LOOP SEQ
      DO JN=0,R_NTMAX
        IN = R_NTMAX+3-JN
        ZJN = JN

        PVOR(IR,IN,kmloc) = -&
         &ZJN*ZEPSNM(KMLOC,JN+1)*PU(IR,IN-1,kmloc)+&
         &(ZJN+1)*ZEPSNM(KMLOC,JN)*PU(IR,IN+1,kmloc)
        PDIV(IR,IN,kmloc) = &
         &ZJN*ZEPSNM(KMLOC,JN+1)*PV(IR,IN-1,kmloc)-&
         &(ZJN+1)*ZEPSNM(KMLOC,JN)*PV(IR,IN+1,kmloc)
      ENDDO
    ENDIF
  ENDDO
ENDDO
!$acc end data
!     ------------------------------------------------------------------

END SUBROUTINE UVTVD
END MODULE UVTVD_MOD
