MODULE GATH_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE GATH_SPEC_CONTROL(PSPECG,KFDISTG,KTO,KVSET,PSPEC)

!**** *GATH_SPEC_CONTROL* - Gather global spectral array from processors

!     Purpose.
!     --------
!        Routine for gathering spectral array

!**   Interface.
!     ----------
!     CALL GATH_SPEC_CONTROL(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KTO(:)    - Processor responsible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     PSPEC(:,:)  - Local spectral array

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE
USE YOMGSTATS, ONLY : LSYNCSTATS

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR
USE ABORT_TRANS_MOD

USE SET2PE_MOD

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPEC(:,:)

REAL(KIND=JPRB)    :: ZFLD(R%NSPEC2_G)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,IBSET,ILEN,JA,ISND
INTEGER(KIND=JPIM) :: IRCV,ISP,ILENR,IPROCA,ISTA,ISTP

!     ------------------------------------------------------------------


!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
IF (.NOT.LSYNCSTATS) CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
  DO JM=1,R%NSPEC2_G
    DO JFLD=1,KFDISTG
      PSPECG(JFLD,JM) =PSPEC(JFLD,JM)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
IF (.NOT.LSYNCSTATS) CALL GSTATS(1644,1)
ELSE
  IFLDR = 0
  IFLDS = 0

  DO JFLD=1,KFDISTG

    IBSET = KVSET(JFLD)
  !Send
    IF( IBSET == MYSETV )THEN

      IFLDS = IFLDS+1
      ISND  = KTO(JFLD)
      ITAG  = MTAGDISTSP+JFLD+17
      ZFLD(1:D%NSPEC2)=PSPEC(IFLDS,1:D%NSPEC2)
      CALL MPL_SEND(ZFLD(1:D%NSPEC2),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
       &CDSTRING='GATH_SPEC_CONTROL')
    ENDIF


  ! RECIEVE
    IF(KTO(JFLD) == MYPROC) THEN
      IFLDR = IFLDR+1
      DO JA=1,NPRTRW
        ILEN = D%NPOSSP(JA+1)-D%NPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL SET2PE(IRCV,0,0,JA,IBSET)
          ITAG = MTAGDISTSP+JFLD+17
          ISTA = D%NPOSSP(JA)
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(ZFLD(ISTA:ISTP),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
           &KOUNT=ILENR,CDSTRING='GATH_SPEC_CONTROL')
          IF( ILENR /= ILEN )THEN
            CALL ABORT_TRANS('GATH_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
          ENDIF
          II = 0
          DO JM=0,R%NSMAX
            IPROCA = D%NPROCM(JM)
            DO JN=JM,R%NSMAX
              ISP = D%NDIM0G(JM)+(JN-JM)*2
              II = II+1
              PSPECG(IFLDR,II) = ZFLD(ISP)
              ISP = ISP+1
              II = II+1
              IF(JM /= 0) THEN
                PSPECG(IFLDR,II) = ZFLD(ISP)
              ELSE
                PSPECG(IFLDR,II) = 0.0_JPRB
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF

  !Synchronize processors
    CALL MPL_BARRIER(CDSTRING='GATH_SPEC_CONTROL:')
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC_CONTROL
END MODULE GATH_SPEC_CONTROL_MOD


