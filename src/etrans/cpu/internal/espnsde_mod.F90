! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ESPNSDE_MOD
CONTAINS
SUBROUTINE ESPNSDE(KM,KF_SCALARS,PF,PNSD)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_GEN
!USE TPM_DIM
!USE TPM_FIELDS
!USE TPM_TRANS
USE TPMALD_DISTR    ,ONLY : DALD
USE TPMALD_GEO      ,ONLY : GALD


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
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_SCALARS
REAL(KIND=JPRB),    INTENT(IN)  :: PF(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PNSD(:,:)

INTEGER(KIND=JPIM) ::   J, JN,IN
REAL(KIND=JPRB)    :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------

!*       1.1      COMPUTE

IF (LHOOK) CALL DR_HOOK('ESPNSDE_MOD:ESPNSDE',0,ZHOOK_HANDLE)
DO JN=1,DALD%NCPL2M(KM),2
  IN =(JN-1)/2
  ZIN = REAL(IN,JPRB)*GALD%EYWN
  DO J=1,2*KF_SCALARS
    PNSD(JN  ,J) = -ZIN*PF(JN+1,J)
    PNSD(JN+1,J) =  ZIN*PF(JN,J)
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('ESPNSDE_MOD:ESPNSDE',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ESPNSDE
END MODULE ESPNSDE_MOD
