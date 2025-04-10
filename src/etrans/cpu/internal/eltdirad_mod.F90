! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ELTDIRAD_MOD
CONTAINS
SUBROUTINE ELTDIRAD(KM,KMLOC,KF_FS,KF_UV,KF_SCALARS,KLED2,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD

USE EPRFI2AD_MOD    ,ONLY : EPRFI2AD
USE ELEDIRAD_MOD    ,ONLY : ELEDIRAD
USE EUVTVDAD_MOD
USE EUPDSPAD_MOD    ,ONLY : EUPDSPAD
 

!**** *ELTDIRAD* - Control of Direct Legendre transform step - adjoint

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *ELTDIRAD(...)*

!        Explicit arguments :
!        --------------------  KM     - zonal wavenumber
!                              KMLOC  - local zonal wavenumber
!                              PSPVOR - spectral vorticity
!                              PSPDIV - spectral divergence
!                              PSPSCALAR - spectral scalar variables

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         EPRFI2AD   - prepares the Fourier work arrays for model variables.
!         ELEDIRAD   - direct Legendre transform
!         EUVTVDAD   -
!         EUPDSPAD   - updating of spectral arrays (fields)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996):
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement +
!        thread-safety
!     ------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(INOUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(INOUT) :: PSPMEANU(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(INOUT) :: PSPMEANV(:)

INTEGER(KIND=JPIM) :: IFC
INTEGER(KIND=JPIM) :: IUS,IUE,IVS,IVE,IVORS,IVORE,IDIVS,IDIVE

REAL(KIND=JPRB) :: ZFFT(RALD%NDGLSUR+R%NNOEXTZG,KLED2)
REAL(KIND=JPRB) :: ZVODI(RALD%NDGLSUR+R%NNOEXTZG,MAX(4*KF_UV,1))
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
!              --------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTDIRAD_MOD:ELTDIRAD',0,ZHOOK_HANDLE)
ZFFT=0.0_JPRB
ZVODI=0.0_JPRB

!     ------------------------------------------------------------------

!*       6.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

CALL EUPDSPAD(KM,KF_UV,KF_SCALARS,ZFFT,ZVODI, &
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC)

!     ------------------------------------------------------------------

!*       5.    COMPUTE VORTICITY AND DIVERGENCE.
!              ---------------------------------
IF( KF_UV > 0 ) THEN
  IUS = 1
  IUE = 2*KF_UV
  IVS = 2*KF_UV+1
  IVE = 4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
! SET PART OF ZFFT CONTAINING U AND V TO 0.
  ZFFT(:,IUS:IVE) = 0.0_JPRB
  CALL EUVTVDAD(KM,KMLOC,KF_UV,KFLDPTRUV,ZFFT(:,IUS:IUE),ZFFT(:,IVS:IVE),&
   & ZVODI(:,IVORS:IVORE),ZVODI(:,IDIVS:IDIVE),PSPMEANU,PSPMEANV)
ENDIF

!     ------------------------------------------------------------------

!*       4.    DIRECT LEGENDRE TRANSFORM.
!              --------------------------
IFC = 2*KF_FS
CALL ELEDIRAD(KM,IFC,KLED2,ZFFT)

!     ------------------------------------------------------------------

!*       3.    FOURIER SPACE COMPUTATIONS.
!              ---------------------------

!     ------------------------------------------------------------------

!*       2.    PREPARE WORK ARRAYS.
!              --------------------

CALL EPRFI2AD(KM,KMLOC,KF_FS,ZFFT)
IF (LHOOK) CALL DR_HOOK('ELTDIRAD_MOD:ELTDIRAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELTDIRAD
END MODULE ELTDIRAD_MOD

