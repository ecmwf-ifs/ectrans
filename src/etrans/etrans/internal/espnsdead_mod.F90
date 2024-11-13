MODULE ESPNSDEAD_MOD
CONTAINS
SUBROUTINE ESPNSDEAD(KM,KF_SCALARS,PF,PNSD)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_GEN
!USE TPM_DIM
!USE TPM_FIELDS
!USE TPM_TRANS

USE TPMALD_DISTR    ,ONLY : DALD
USE TPMALD_GEO      ,ONLY : GALD


!**** *ESPNSDEAD* - Compute North-South derivative in spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the North-south derivative

!**   Interface.
!     ----------
!        CALL ESPNSDEAD(...)

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
!        Original : 00-02-01 From SPNSDEAD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KM
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_SCALARS
REAL(KIND=JPRB),    INTENT(INOUT) :: PF(:,:)
REAL(KIND=JPRB),    INTENT(IN)    :: PNSD(:,:)
INTEGER(KIND=JPIM) ::  ISKIP, J, JN
INTEGER(KIND=JPIM) :: IN
REAL(KIND=JPRB)    :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------

!*       1.1      COMPUTE

IF (LHOOK) CALL DR_HOOK('ESPNSDEAD_MOD:ESPNSDEAD',0,ZHOOK_HANDLE)
IF(KM == 0) THEN
  ISKIP = 1
ELSE
  ISKIP = 1
ENDIF

DO JN=1,DALD%NCPL2M(KM),2

  IN = (JN-1)/2
  ZIN = REAL(IN,JPRB)*GALD%EYWN

  DO J=1,2*KF_SCALARS,ISKIP

    PF(JN+1,J) = PF(JN+1,J)-ZIN*PNSD(JN  ,J)
    PF(JN  ,J) = PF(JN  ,J)+ZIN*PNSD(JN+1,J)

  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('ESPNSDEAD_MOD:ESPNSDEAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ESPNSDEAD
END MODULE ESPNSDEAD_MOD
