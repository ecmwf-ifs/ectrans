! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EUVTVDAD_MOD
CONTAINS
SUBROUTINE EUVTVDAD(KM,KMLOC,KFIELD,KFLDPTR,PU,PV,PVOR,PDIV,PSPMEANU,PSPMEANV)

!**** *EUVTVDAD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX.

!**   Interface.
!     ----------
!        CALL EUVTVDAD()

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KFIELD - number of fields (levels)
!                              KFLDPTR - fields pointers
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
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        03-03-03   G. Radnoti: b-level conform mean wind distribution
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        01-Dec-2004   A. Deckmyn    removed erasing of mean wind 
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
!USE TPM_FIELDS

USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD, KM, KMLOC
REAL(KIND=JPRB), INTENT(IN)    :: PVOR(:,:),PDIV(:,:)
REAL(KIND=JPRB), INTENT(INOUT) :: PU  (:,:),PV  (:,:)

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KFLDPTR(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(INOUT) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IN, IR, J, JN, IFLD

REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVDAD_MOD:EUVTVDAD',0,ZHOOK_HANDLE)

IF (KM == 0) THEN
  IF (PRESENT(KFLDPTR)) THEN
    DO J=1,KFIELD
      IR=2*J-1
      IFLD=KFLDPTR(J)
      PU(1,IR)=PSPMEANU(IFLD)
      PV(1,IR)=PSPMEANV(IFLD)
    ENDDO
  ELSE
    DO J=1,KFIELD
      IR=2*J-1
      PU(1,IR)=PSPMEANU(J)
      PV(1,IR)=PSPMEANV(J)
    ENDDO
  ENDIF
ENDIF

DO J=1,2*KFIELD
  DO JN=1,DALD%NCPL2M(KM),2
    IN=(JN-1)/2
    ZIN=REAL(IN,JPRB)*GALD%EYWN
    PU(JN+1,J) =  PU(JN+1,J) + ZIN * PVOR(JN  ,J)
    PU(JN  ,J) =  PU(JN  ,J) - ZIN * PVOR(JN+1,J)
    PV(JN+1,J) =  PV(JN+1,J) - ZIN * PDIV(JN  ,J)
    PV(JN  ,J) =  PV(JN  ,J) + ZIN * PDIV(JN+1,J)
  ENDDO
ENDDO

ZKM=REAL(KM,JPRB)*GALD%EXWN
DO J=1,KFIELD
  IR=2*J-1
  II=IR+1
  DO JN=1,R%NDGL+R%NNOEXTZG
    PU(JN,II) = PU(JN,II) - ZKM * PDIV(JN,IR)
    PU(JN,IR) = PU(JN,IR) + ZKM * PDIV(JN,II)
    PV(JN,II) = PV(JN,II) - ZKM * PVOR(JN,IR)
    PV(JN,IR) = PV(JN,IR) + ZKM * PVOR(JN,II)
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('EUVTVDAD_MOD:EUVTVDAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EUVTVDAD
END MODULE EUVTVDAD_MOD
