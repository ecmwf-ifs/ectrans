MODULE EVDTUV_MOD
CONTAINS
SUBROUTINE EVDTUV(KFIELD,KFLDPTR,PVOR,PDIV,PU,PV,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPMALD_FIELDS   ,ONLY : FALD
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD
USE TPM_DISTR       ,ONLY : D
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

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

INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELD
REAL(KIND=JPRB),    INTENT(IN)  :: PVOR(:,:,:),PDIV(:,:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PU  (:,:,:),PV  (:,:,:)

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KFLDPTR(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(IN) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, IN, IFLD
INTEGER(KIND=JPIM) :: JM, IM
INTEGER(KIND=JPIM) :: JNMAX

REAL(KIND=JPRB) :: ZLEPINM
REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('EVDTUV_MOD:EVDTUV',0,ZHOOK_HANDLE)

JNMAX = MAXVAL (DALD%NCPL2M)

!$acc parallel loop collapse (3) private (JM, J, JN, IM, IN, ZIN) &
!$acc & present (D%NUMP, D%MYMS, DALD%NCPL2M, PU, PV, PVOR, PDIV)
DO J=1,2*KFIELD
  DO JM = 1, D%NUMP
    DO JN=1,JNMAX,2
      IM = D%MYMS (JM)
      IF (JN <= DALD%NCPL2M(IM)) THEN
        IN = (JN-1)/2
        ZIN = REAL(IN,JPRB)*GALD%EYWN
        PU(JN  ,JM,J) = -ZIN*PVOR(JN+1,JM,J)
        PU(JN+1,JM,J) =  ZIN*PVOR(JN  ,JM,J)
        PV(JN  ,JM,J) = -ZIN*PDIV(JN+1,JM,J)
        PV(JN+1,JM,J) =  ZIN*PDIV(JN  ,JM,J)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

!$acc parallel loop collapse (3) private (JM, J, JN, IM, ZKM, IR, II, IJ, ZLEPINM) &
!$acc & present (D%NUMP, D%MYMS, DALD%NCPL2M, FALD%RLEPINM, PU, PV, PDIV, PVOR)
DO J=1,KFIELD
  DO JM = 1, D%NUMP
    DO JN=1,JNMAX
      IM = D%MYMS (JM)
      ZKM=REAL(IM,JPRB)*GALD%EXWN
      IR = 2*J-1
      II = IR+1
      IF (JN <= DALD%NCPL2M(IM)) THEN
        IJ=(JN-1)/2
        ZLEPINM = FALD%RLEPINM(DALD%NPME(IM)+IJ)
        PU(JN,JM,IR)= ZLEPINM*(-ZKM*PDIV(JN,JM,II)-PU(JN,JM,IR))
        PU(JN,JM,II)= ZLEPINM*( ZKM*PDIV(JN,JM,IR)-PU(JN,JM,II))
        PV(JN,JM,IR)= ZLEPINM*(-ZKM*PVOR(JN,JM,II)+PV(JN,JM,IR))
        PV(JN,JM,II)= ZLEPINM*( ZKM*PVOR(JN,JM,IR)+PV(JN,JM,II))
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

IF (PRESENT(KFLDPTR)) THEN
!$acc parallel loop collapse (2) private (J, JM, IM, IR, IFLD) &
!$acc & present (D%NUMP, D%MYMS, PU, PV, PSPMEANU, PSPMEANV) copyin (KFLDPTR)
  DO J = 1, KFIELD
    DO JM = 1, D%NUMP
      IM = D%MYMS (JM)
      IF (IM == 0) THEN
        IR = 2*J-1
        IFLD=KFLDPTR(J)
        PU(1,JM,IR)=PSPMEANU(IFLD)
        PV(1,JM,IR)=PSPMEANV(IFLD)
      ENDIF
    ENDDO
  ENDDO
!$acc end parallel loop
ELSE
!$acc parallel loop collapse (2) private (J, JM, IM, IR) &
!$acc & present (D%NUMP, D%MYMS, PU, PV, PSPMEANU, PSPMEANV)
  DO J = 1, KFIELD
    DO JM = 1, D%NUMP
      IM = D%MYMS (JM)
      IF (IM == 0) THEN
        IR = 2*J-1
        PU(1,JM,IR)=PSPMEANU(J)
        PV(1,JM,IR)=PSPMEANV(J)
      ENDIF
    ENDDO
  ENDDO
!$acc end parallel loop
ENDIF

IF (LHOOK) CALL DR_HOOK('EVDTUV_MOD:EVDTUV',1,ZHOOK_HANDLE)

END SUBROUTINE EVDTUV
END MODULE EVDTUV_MOD