! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EVDTUV_MOD
CONTAINS
SUBROUTINE EVDTUV(KM,KFIELD,KFLDPTR,PVOR,PDIV,PU,PV,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
!USE TPM_FIELDS
USE TPMALD_FIELDS   ,ONLY : FALD
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD

!**** *VDTUV* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL VDTUV(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              KFLDPTR - fields pointers
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!                              PVOR(NLEI1,2*KFIELD) - vorticity (input)
!                              PDIV(NLEI1,2*KFIELD) - divergence (input)
!                              PU(NLEI1,2*KFIELD)   - u wind (output)
!                              PV(NLEI1,2*KFIELD)   - v wind (output)
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

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

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
!        Original : 00-02-01 From VDTUV in IFS CY22R1
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM, KFIELD
REAL(KIND=JPRB),    INTENT(IN)  :: PVOR(:,:),PDIV(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PU  (:,:),PV  (:,:)

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KFLDPTR(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(IN) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, IN, IFLD

REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('EVDTUV_MOD:EVDTUV',0,ZHOOK_HANDLE)
ZKM=REAL(KM,JPRB)*GALD%EXWN
DO J=1,2*KFIELD
  DO JN=1,DALD%NCPL2M(KM),2
    IN = (JN-1)/2
    ZIN = REAL(IN,JPRB)*GALD%EYWN
    PU(JN  ,J) = -ZIN*PVOR(JN+1,J)
    PU(JN+1,J) =  ZIN*PVOR(JN,J)
    PV(JN  ,J) = -ZIN*PDIV(JN+1,J)
    PV(JN+1,J) =  ZIN*PDIV(JN,J)
  ENDDO
ENDDO
DO J=1,KFIELD
  IR = 2*J-1
  II = IR+1
  DO JN=1,DALD%NCPL2M(KM)
    IJ=(JN-1)/2
    PU(JN,IR)= FALD%RLEPINM(DALD%NPME(KM)+IJ)*(-ZKM*PDIV(JN,II)-PU(JN,IR))
    PU(JN,II)= FALD%RLEPINM(DALD%NPME(KM)+IJ)*( ZKM*PDIV(JN,IR)-PU(JN,II))
    PV(JN,IR)= FALD%RLEPINM(DALD%NPME(KM)+IJ)*(-ZKM*PVOR(JN,II)+PV(JN,IR))
    PV(JN,II)= FALD%RLEPINM(DALD%NPME(KM)+IJ)*( ZKM*PVOR(JN,IR)+PV(JN,II))
  ENDDO
ENDDO
IF (KM == 0) THEN
  IF (PRESENT(KFLDPTR)) THEN
    DO J = 1, KFIELD
      IR = 2*J-1
      IFLD=KFLDPTR(J)
      PU(1,IR)=PSPMEANU(IFLD)
      PV(1,IR)=PSPMEANV(IFLD)
    ENDDO
  ELSE
    DO J = 1, KFIELD
      IR = 2*J-1
      PU(1,IR)=PSPMEANU(J)
      PV(1,IR)=PSPMEANV(J)
    ENDDO
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('EVDTUV_MOD:EVDTUV',1,ZHOOK_HANDLE)

END SUBROUTINE EVDTUV
END MODULE EVDTUV_MOD
