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
USE TPM_FIELDS      ,ONLY : F_RN
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_FIELDS      ,ONLY : ZOA1,ZOA2,ZEPSNM
!

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELD
INTEGER(KIND=JPIM)  :: KM, KMLOC

!REAL(KIND=JPRBT), INTENT(IN)     :: PEPSNM(1:d%nump,0:R%NTMAX+2)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IN, IR, J, JN, ITMAX
INTEGER(KIND=JPIM) :: I_DIV_OFFSET

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM
REAL(KIND=JPRBT) :: ZN(-1:R%NTMAX+3)

! ZOA1 and ZOA2 are arranged with vorticity/U from 1 to 2 * KFIELD and divergence/V from
! 2 * KFIELD + 1 to 4 * KFIELD
I_DIV_OFFSET = 2 * KFIELD

!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

#ifdef ACCGPU
!$ACC DATA &
!$ACC& CREATE(ZN) COPYIN(I_DIV_OFFSET) &
!$ACC& PRESENT(D_MYMS,D_NUMP,R_NTMAX,F_RN,ZEPSNM,ZOA1,ZOA2)
#endif
#ifdef OMPGPU
!WARNING: following line should be PRESENT,ALLOC but causes issues with AMD compiler!
!$OMP TARGET DATA&
!$OMP& MAP(ALLOC:ZN) &
!$OMP& MAP(TO:D_MYMS,D_NUMP,R_NTMAX) &
!$OMP& MAP(TO:F_RN) &
!$OMP& MAP(ALLOC:ZEPSNM,ZOA1,ZOA2)
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!! DEFAULT(NONE) SHARED(R_NTMAX,ZN,F_RN)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP DEFAULT(NONE) PRESENT(R_NTMAX,ZN,F_RN)
#endif
DO J=-1,R_NTMAX+3
  ZN(j) = F_RN(j)
ENDDO
!*       1.1      SET N=KM-1 COMPONENT TO 0 FOR U AND V

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
!! PRIVATE(KM,IN) DEFAULT(NONE) &
!!$OMP& SHARED(D_NUMP,KFIELD,D_MYMS)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IN) DEFAULT(NONE) &
!$ACC& COPYIN(KFIELD) &
!$ACC& PRESENT(D_NUMP,D_MYMS,F_RN,R_NTMAX,ZOA1)
#endif
DO KMLOC=1,D_NUMP
  DO J=1,2*KFIELD
    KM = D_MYMS(KMLOC)
    !IN = F%NLTN(KM-1)
    IN=R_NTMAX+3-KM
    ZOA1(J,IN,KMLOC) = 0.0_JPRBT
    ZOA1(I_DIV_OFFSET+J,IN,KMLOC) = 0.0_JPRBT
  ENDDO
ENDDO

!*       1.2      COMPUTE VORTICITY AND DIVERGENCE.

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(IR,II,IN,KM,ZKM) DEFAULT(NONE) &
!$OMP& SHARED(D_NUMP,R_NTMAX,KFIELD,D_MYMS,ZN,ZEPSNM,I_DIV_OFFSET,ZOA1,ZOA2)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IR,II,IN,KM,ZKM) DEFAULT(NONE) &
!$ACC& COPYIN(KFIELD) &
!$ACC& PRESENT(D_NUMP,R_NTMAX,D_MYMS,ZN,ZEPSNM)
#endif
DO KMLOC=1,D_NUMP
  DO JN=0,R_NTMAX
    DO J=1,KFIELD
      IR = 2*J-1
      II = IR+1
      KM = D_MYMS(KMLOC)
      ZKM = REAL(KM,JPRBT)
      IN = R_NTMAX+2-JN

      IF (KM /= 0 .AND. JN >= KM) THEN
        ! Vorticity
        ZOA2(IR,IN,KMLOC) = -ZKM * ZOA1(I_DIV_OFFSET+II,IN,KMLOC) &
          & - ZN(JN)   * ZEPSNM(KMLOC,JN+1) * ZOA1(IR,IN-1,KMLOC) &
          & + ZN(JN+1) * ZEPSNM(KMLOC,JN)   * ZOA1(IR,IN+1,KMLOC)
        ZOA2(II,IN,KMLOC) = ZKM * ZOA1(I_DIV_OFFSET+IR,IN,KMLOC) &
          & - ZN(JN)   * ZEPSNM(KMLOC,JN+1) * ZOA1(II,IN-1,KMLOC) &
          & + ZN(JN+1) * ZEPSNM(KMLOC,JN)   * ZOA1(II,IN+1,KMLOC)

        ! Divergence
        ZOA2(I_DIV_OFFSET+IR,IN,KMLOC) = -ZKM * ZOA1(II,IN,KMLOC) &
          & + ZN(JN)   * ZEPSNM(KMLOC,JN+1) * ZOA1(I_DIV_OFFSET+IR,IN-1,KMLOC) &
          & - ZN(JN+1) * ZEPSNM(KMLOC,JN)   * ZOA1(I_DIV_OFFSET+IR,IN+1,KMLOC)
        ZOA2(I_DIV_OFFSET+II,IN,KMLOC) = ZKM * ZOA1(IR,IN,KMLOC) &
          & + ZN(JN)   * ZEPSNM(KMLOC,JN+1) * ZOA1(I_DIV_OFFSET+II,IN-1,KMLOC) &
          & - ZN(JN+1) * ZEPSNM(KMLOC,JN)   * ZOA1(I_DIV_OFFSET+II,IN+1,KMLOC)
      ELSEIF (KM == 0) THEN
        ! Vorticity
        ZOA2(IR,IN,KMLOC) = -ZN(JN) * ZEPSNM(KMLOC,JN+1) * ZOA1(IR,IN-1,KMLOC) &
          & + ZN(JN+1) * ZEPSNM(KMLOC,JN) * ZOA1(IR,IN+1,KMLOC)

        ! Divergence
        ZOA2(I_DIV_OFFSET+IR,IN,KMLOC) = &
          &   ZN(JN) * ZEPSNM(KMLOC,JN+1) * ZOA1(I_DIV_OFFSET+IR,IN-1,KMLOC) &
          & - ZN(JN+1) * ZEPSNM(KMLOC,JN) * ZOA1(I_DIV_OFFSET+IR,IN+1,KMLOC)
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
