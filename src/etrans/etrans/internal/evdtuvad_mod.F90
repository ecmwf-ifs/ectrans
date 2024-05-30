MODULE EVDTUVAD_MOD
CONTAINS
SUBROUTINE EVDTUVAD(KM,KMLOC,KFIELD,KFLDPTR,PVOR,PDIV,PU,PV,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
!USE TPM_FIELDS
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRW, MYSETV, MYSETW, MYPROC, NPROC

USE TPMALD_FIELDS   ,ONLY : FALD
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD

!**** *EVDTUVAD* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL EVDTUVAD(...)

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
!        Original : 00-02-01 From VDTUVAD in IFS CY22R1
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        01-Dec-2004   A. Deckmyn    Fix mean wind for NPRTRW > 1
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement +
!        thread-safety
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KM, KFIELD, KMLOC
REAL(KIND=JPRB),    INTENT(INOUT) :: PVOR(:,:),PDIV(:,:)
REAL(KIND=JPRB),    INTENT(INOUT) :: PU  (:,:),PV  (:,:)

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KFLDPTR(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(OUT) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, IFLD

INTEGER(KIND=JPIM) :: IN
INTEGER(KIND=JPIM) :: ISND, JA, ITAG, ILEN

REAL(KIND=JPRB) :: ZSPU(2*KFIELD)
REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EVDTUVAD_MOD:EVDTUVAD',0,ZHOOK_HANDLE)

IF (KM == 0) THEN
  IF (PRESENT(KFLDPTR)) THEN
    DO J = 1, KFIELD
      IR = 2*J-1
      IFLD=KFLDPTR(J)
      PSPMEANU(IFLD)=PU(1,IR)
      PSPMEANV(IFLD)=PV(1,IR)
    ENDDO
  ELSE
    DO J = 1, KFIELD
      IR = 2*J-1
      PSPMEANU(J)=PU(1,IR)
      PSPMEANV(J)=PV(1,IR)
    ENDDO
  ENDIF
ENDIF

ZKM=REAL(KM,JPRB)*GALD%EXWN
DO J=1,KFIELD
  IR = 2*J-1
  II = IR+1
  DO JN=1,DALD%NCPL2M(KM)
    IJ=(JN-1)/2
    PDIV(JN,II)=PDIV(JN,II)-ZKM*FALD%RLEPINM(DALD%NPME(KM)+IJ)*PU(JN,IR)
    PU(JN,IR)=-FALD%RLEPINM(DALD%NPME(KM)+IJ)*PU(JN,IR)

    PDIV(JN,IR)=PDIV(JN,IR)+ZKM*FALD%RLEPINM(DALD%NPME(KM)+IJ)*PU(JN,II)
    PU(JN,II)=-FALD%RLEPINM(DALD%NPME(KM)+IJ)*PU(JN,II)

    PVOR(JN,II)=PVOR(JN,II)-ZKM*FALD%RLEPINM(DALD%NPME(KM)+IJ)*PV(JN,IR)
    PV(JN,IR)=FALD%RLEPINM(DALD%NPME(KM)+IJ)*PV(JN,IR)

    PVOR(JN,IR)=PVOR(JN,IR)+ZKM*FALD%RLEPINM(DALD%NPME(KM)+IJ)*PV(JN,II)
    PV(JN,II)=FALD%RLEPINM(DALD%NPME(KM)+IJ)*PV(JN,II)

  ENDDO
ENDDO

DO J=1,2*KFIELD
  DO JN=1,DALD%NCPL2M(KM),2
    IN = (JN-1)/2
    ZIN = REAL(IN,JPRB)*GALD%EYWN
    PVOR(JN+1,J) = PVOR(JN+1,J)-ZIN*PU(JN  ,J)
    PVOR(JN  ,J) = PVOR(JN  ,J)+ZIN*PU(JN+1,J)
    PDIV(JN+1,J) = PDIV(JN+1,J)-ZIN*PV(JN  ,J)
    PDIV(JN  ,J) = PDIV(JN  ,J)+ZIN*PV(JN+1,J)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('EVDTUVAD_MOD:EVDTUVAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EVDTUVAD
END MODULE EVDTUVAD_MOD
