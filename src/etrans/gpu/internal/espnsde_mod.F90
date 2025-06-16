MODULE ESPNSDE_MOD
CONTAINS
SUBROUTINE ESPNSDE(KF_SCALARS,PF,PNSD)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, D_NUMP, D_MYMS
USE TPMALD_DISTR    ,ONLY : DALD, DALD_NCPL2M
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

INTEGER(KIND=JPIM), INTENT(IN)  :: KF_SCALARS
REAL(KIND=JPRB),    INTENT(IN)  :: PF(:,:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PNSD(:,:,:)

INTEGER(KIND=JPIM) :: J, JN,IN, JM, IM, JNMAX
REAL(KIND=JPRB)    :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------

!*       1.1      COMPUTE

IF (LHOOK) CALL DR_HOOK('ESPNSDE_MOD:ESPNSDE',0,ZHOOK_HANDLE)

JNMAX = MAXVAL (DALD%NCPL2M)

!$acc parallel loop collapse (3) private (JM, J, JN, IM, IN, ZIN) &
!$acc & present (D_NUMP, D_MYMS, DALD_NCPL2M, PNSD, PF)
DO J=1,2*KF_SCALARS
  DO JM = 1, D_NUMP
    DO JN=1,JNMAX,2
      IM = D_MYMS(JM)
      IF (JN <= DALD_NCPL2M(IM)) THEN
        IN =(JN-1)/2
        ZIN = REAL(IN,JPRB)*GALD%EYWN
        PNSD(JN  ,JM,J) = -ZIN*PF(JN+1,JM,J)
        PNSD(JN+1,JM,J) =  ZIN*PF(JN  ,JM,J)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

IF (LHOOK) CALL DR_HOOK('ESPNSDE_MOD:ESPNSDE',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ESPNSDE
END MODULE ESPNSDE_MOD