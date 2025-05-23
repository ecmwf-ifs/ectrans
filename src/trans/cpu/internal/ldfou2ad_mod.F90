! (C) Copyright 1991- ECMWF.
! (C) Copyright 1991- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LDFOU2AD_MOD
CONTAINS
SUBROUTINE LDFOU2AD(KM,KF_UV,PAIA,PSIA)

!**** *LDFOU2AD* - Division by a*cos(theta) of u and v

!     Purpose.
!     --------
!        In Fourier space divide u and v by  a*cos(theta).

!**   Interface.
!     ----------
!        CALL LDFOU2AD(KM,PAIA,PSIA)

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

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FIELDS      ,ONLY : F
!

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KM,KF_UV

REAL(KIND=JPRB) ,INTENT(INOUT) :: PSIA(:,:),   PAIA(:,:)

!     LOCAL REAL SCALARS
REAL(KIND=JPRB)   :: ZACTHE

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: J, JGL ,IFLD ,ISL


!     ------------------------------------------------------------------

!*       1.    DIVIDE U V BY A*COS(THETA)
!              --------------------------

ISL  = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IFLD = 4*KF_UV

!*       1.1      U AND V

DO JGL=ISL,R%NDGNH
  ZACTHE = REAL(F%RACTHE(JGL),JPRB)
  DO J=1,IFLD
    PAIA(J,JGL) = PAIA(J,JGL)*ZACTHE
    PSIA(J,JGL) = PSIA(J,JGL)*ZACTHE
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE LDFOU2AD
END MODULE LDFOU2AD_MOD
