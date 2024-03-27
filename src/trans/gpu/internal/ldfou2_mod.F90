! (C) Copyright 1991- ECMWF.
! (C) Copyright 1991- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LDFOU2_MOD
CONTAINS
SUBROUTINE LDFOU2(KF_UV, P_FOURIER_FIELDS)

!**** *LDFOU2* - Division by a*cos(theta) of u and v

!     Purpose.
!     --------
!        In Fourier space divide u and v by  a*cos(theta).

!**   Interface.
!     ----------
!        CALL LDFOU2(KM,P_FOURIER_FIELDS)

!        Explicit arguments :
!        --------------------  KM - zonal wavenumber
!                              P_FOURIER_FIELDS - Fourier fields

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
!        Original : 91-07-01
!        Modified : 94-04-06 R. El Khatib - Full-POS configuration 'P'
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                    instead of u,v->vor,div
!        MPP Group: 95-10-01 Message Passing option added
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS, ONLY : JPIM, JPRBT
USE TPM_FIELDS,      ONLY : F_RACTHE
USE TPM_DISTR,       ONLY : D, D_NUMP, D_MYMS
USE TPM_DIM,         ONLY : R, R_NDGNH, R_NDGL
USE TPM_GEOMETRY,    ONLY : G, G_NDGLU
USE ABORT_TRANS_MOD  ,ONLY : ABORT_TRANS

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
REAL(KIND=JPRBT), INTENT(INOUT) :: P_FOURIER_FIELDS(:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: J, JGL, ISL, KM, KMLOC

!     ------------------------------------------------------------------

!*       1.    DIVIDE U V BY A*COS(THETA)
!              --------------------------

IF (KF_UV > 0) THEN

#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
  !$OMP& PRIVATE(KMLOC,JGL,J,KM,ISL) MAP(TO:KF_UV,F_RACTHE,D_MYMS,G_NDGLU,P_FOURIER_FIELDS) &
  !$OMP& SHARED(D_NUMP,R_NDGNH,KF_UV,F_RACTHE,D_MYMS,G_NDGLU,P_FOURIER_FIELDS)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) &
  !$ACC& PRIVATE(KMLOC,JGL,J,KM,ISL) COPYIN(KF_UV,F_RACTHE) &
  !$ACC& PRESENT(D_NUMP,D_MYMS,R_NDGNH,R_NDGL,G_NDGLU,P_FOURIER_FIELDS)
#endif
  DO KMLOC = 1, D_NUMP
    DO JGL = 1, R_NDGNH
      DO J = 1, 4 * KF_UV
         KM = D_MYMS(KMLOC)
         ISL = MAX(R_NDGNH - G_NDGLU(KM) + 1, 1)
         IF (JGL >= ISL) THEN
           P_FOURIER_FIELDS(J,JGL,KMLOC) = P_FOURIER_FIELDS(J,JGL,KMLOC) * F_RACTHE(JGL)
         ENDIF
       ENDDO
     ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LDFOU2
END MODULE LDFOU2_MOD
