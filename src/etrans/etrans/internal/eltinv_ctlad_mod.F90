MODULE ELTINV_CTLAD_MOD
CONTAINS
SUBROUTINE ELTINV_CTLAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC,PSPMEANU,PSPMEANV)

!**** *ELTINV_CTLAD* - Control routine for inverse Legandre transform - adj.

!     Purpose.
!     --------
!     Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL EINV_TRANS_CTL(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KFLDPTRUV(:) - field pointer array for vor./div.
!     KFLDPTRSC(:) - field pointer array for PSPSCALAR
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition

!     Method.
!     -------

!     Externals.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-06-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        O.Spaniel     Oct-2004 phasing for AL29
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement +
!        thread-safety
!   R. El Khatib 02-Jun-2022 Optimization/Cleaning
!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_GEN         ,ONLY : LALLOPERM
!USE TPM_DIM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D
USE ELTINVAD_MOD    ,ONLY : ELTINVAD
USE TRLTOM_MOD      ,ONLY : TRLTOM
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCDERS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT)  :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT)  :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT)  :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPMEANU(:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPMEANV(:)

EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

INTEGER(KIND=JPIM) :: IBLEN, ILEI2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTINV_CTLAD_MOD:ELTINV_CTLAD',0,ZHOOK_HANDLE)

ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IBLEN = D%NLENGT0B*2*KF_OUT_LT
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
    FOUBUF_IN(1)=0._JPRB ! force allocation here
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  FOUBUF_IN(1)=0._JPRB ! force allocation here
ENDIF
CALL GSTATS(180,0)
CALL TRLTOM(FOUBUF,FOUBUF_IN,2*KF_OUT_LT)
CALL GSTATS(180,1)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)

CALL GSTATS(1648,0)
IF(KF_OUT_LT > 0) THEN
  CALL ELTINVAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,&
   & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2, &
   & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC,PSPMEANU,PSPMEANV)
ENDIF
CALL GSTATS(1648,1)

IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)
IF (LHOOK) CALL DR_HOOK('ELTINV_CTLAD_MOD:ELTINV_CTLAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELTINV_CTLAD
END MODULE ELTINV_CTLAD_MOD
