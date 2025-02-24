! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EUVTVD_MOD
CONTAINS
SUBROUTINE EUVTVD(KM,KMLOC,KFIELD,PU,PV,PVOR,PDIV)

!**** *EUVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX - calculation part.

!**   Interface.
!     ----------
!        CALL EUVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

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
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        03-03-03 : G. Radnoti: b-level conform mean-wind distribution
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!        R. El Khatib 02-Jun-2022 Optimization/Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRB),    INTENT(IN) :: PU(:,:)
REAL(KIND=JPRB),    INTENT(IN) :: PV(:,:)
REAL(KIND=JPRB),    INTENT(OUT):: PVOR(:,:)
REAL(KIND=JPRB),    INTENT(OUT):: PDIV(:,:)

INTEGER(KIND=JPIM) :: II, IN, IR, J, JN

REAL(KIND=JPRB) :: ZKM, ZIN

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',0,ZHOOK_HANDLE)

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ZKM=REAL(KM,JPRB)*GALD%EXWN
DO J=1,KFIELD
  IR=2*J-1
  II=IR+1
  DO JN=1,R%NDGL+R%NNOEXTZG
    PDIV(JN,IR)=-ZKM*PU(JN,II)
    PDIV(JN,II)= ZKM*PU(JN,IR)
    PVOR(JN,IR)=-ZKM*PV(JN,II)
    PVOR(JN,II)= ZKM*PV(JN,IR)
  ENDDO
ENDDO
DO J=1,2*KFIELD
  DO JN=1,DALD%NCPL2M(KM),2
    IN=(JN-1)/2
    ZIN=REAL(IN,JPRB)*GALD%EYWN
    PVOR(JN,J  )=PVOR(JN  ,J)+ZIN*PU(JN+1,J)
    PVOR(JN+1,J)=PVOR(JN+1,J)-ZIN*PU(JN  ,J)
    PDIV(JN,J  )=PDIV(JN  ,J)-ZIN*PV(JN+1,J)
    PDIV(JN+1,J)=PDIV(JN+1,J)+ZIN*PV(JN  ,J)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',1,ZHOOK_HANDLE)

END SUBROUTINE EUVTVD
END MODULE EUVTVD_MOD
