! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_CTLAD_MOD
CONTAINS
SUBROUTINE FTINV_CTLAD(KF_UV_G,KF_SCALARS_G,&
 & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,KVSETUV,KVSETSC,KPTRGP, &
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *FTINV_CTLAD - Inverse Fourier transform control - adjoint

!     Purpose. Control routine for Fourier to Gridpoint transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV_CTLAD(..)

!        Explicit arguments :
!        --------------------
!        PGP     -  gridpoint array
!        KF_UV_G      - global number of spectral u-v fields
!        KF_SCALARS_G - global number of scalar spectral fields
!        KF_UV        - local number of spectral u-v fields
!        KF_SCALARS   - local number of scalar spectral fields
!        KF_SCDERS    - local number of derivatives of scalar spectral fields
!        KF_GP        - total number of output gridpoint fields
!        KF_FS        - total number of fields in fourier space
!        KF_OUT_LT    - total number of fields coming out from inverse LT
!        KVSETUV - "B"  set in spectral/fourier space for
!                   u and v variables
!        KVSETSC - "B" set in spectral/fourier space for
!                  scalar variables
!        KPTRGP - pointer array to fi3elds in gridpoint space


!     Method.
!     -------

!     Externals.  TRLTOG      - transposition routine
!     ----------  FOURIER_IN  - copy fourier data from Fourier buffer
!                 FTINV       - fourier transform
!                 FSC         - Fourier space computations

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_SCDERS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)

!     ------------------------------------------------------------------

CALL ABORT_TRANS("FTINV_CTLAD: Not supported for overlap prototype")

!     ------------------------------------------------------------------

END SUBROUTINE FTINV_CTLAD
END MODULE FTINV_CTLAD_MOD



