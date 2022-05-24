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

MODULE LDFOU2_MOD
CONTAINS
SUBROUTINE LDFOU2(KF_UV,PAIA)

!**** *LDFOU2* - Division by a*cos(theta) of u and v

!     Purpose.
!     --------
!        In Fourier space divide u and v by  a*cos(theta).

!**   Interface.
!     ----------
!        CALL LDFOU2(KM,PAIA,PSIA)

!        Explicit arguments :
!        --------------------  KM - zonal wavenumber
!                              PAIA - antisymmetric fourier fields
!                              PSIA - symmetric fourierfields

!        Implicit arguments :  RACTHE - 1./(a*cos(theta))
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
!        Original : 91-07-01
!        Modified : 94-04-06 R. El Khatib - Full-POS configuration 'P'
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                    instead of u,v->vor,div
!        MPP Group: 95-10-01 Message Passing option added
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRBT

USE TPM_FIELDS      ,ONLY : F_RACTHE
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_DIM         ,ONLY : R, R_NDGNH, R_NDGL
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU

!

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM) :: KM,KMLOC

REAL(KIND=JPRBT) ,INTENT(INOUT) :: PAIA(:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: J, JGL, ISL

!     ------------------------------------------------------------------

!*       1.    DIVIDE U V BY A*COS(THETA)
!              --------------------------

IF( KF_UV > 0 ) THEN

#ifdef ACCGPU
!$ACC DATA &
!$ACC& COPYIN(KF_UV,F_RACTHE) &
!$ACC& PRESENT(D_NUMP,D_MYMS,R_NDGNH,R_NDGL,G_NDGLU) &
!$ACC& PRESENT(PAIA)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA &
!$OMP& MAP(TO:F_RACTHE,D,D_NUMP,D_MYMS,R_NDGNH,R_NDGL,G_NDGLU) &
!$OMP& MAP(TOFROM:PAIA)
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
!! DEFAULT(NONE) PRIVATE(KM,ISL) &
!!$OMP&      SHARED(D_NUMP,R_NDGNH,KF_UV,D_MYMS,R_NDGNH,G_NDGLU,R_NDGL,PAIA,F_RACTHE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISL) DEFAULT(NONE) &
!$ACC&      PRESENT(D_NUMP,R_NDGNH,KF_UV,D_MYMS,R_NDGNH,G_NDGLU,R_NDGL,PAIA,F_RACTHE)
#endif
DO KMLOC=1,D_NUMP
  DO J=1,4*KF_UV ! (real+complex) * (U+V)
    KM = D_MYMS(KMLOC)
    ISL  = R_NDGNH-G_NDGLU(KM)+1
    !$ACC LOOP SEQ
    DO JGL=ISL,R_NDGNH
!*       1.1      U AND V
      PAIA(J,JGL,KMLOC) = PAIA(J,JGL,KMLOC)*F_RACTHE(JGL)
    ENDDO
  ENDDO
ENDDO
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LDFOU2
END MODULE LDFOU2_MOD
