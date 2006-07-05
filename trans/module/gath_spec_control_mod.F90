MODULE GATH_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE GATH_SPEC_CONTROL(PSPECG,KFDISTG,KTO,KVSET,PSPEC,LDIM1_IS_FLD,&
 & KSMAX,KSPEC2,KSPEC2_G,KPOSSP,KDIM0G)

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

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR
USE ABORT_TRANS_MOD

USE SET2PE_MOD
USE SUWAVEDI_MOD

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2_G
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPOSSP(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KDIM0G(0:)

REAL(KIND=JPRB)    :: ZFLD(KSPEC2_G),ZRECV(KSPEC2_G)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,IBSET,ILEN,JA,ISND
INTEGER(KIND=JPIM) :: IRCV,ISP,ILENR,ISTA,ISTP,ISENDREQ,IPOS0

!     ------------------------------------------------------------------


!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
  CALL GSTATS(1644,0)
  IF(LDIM1_IS_FLD) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
    DO JM=1,KSPEC2_G
      DO JFLD=1,KFDISTG
        PSPECG(JFLD,JM) =PSPEC(JFLD,JM)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
    DO JFLD=1,KFDISTG
      DO JM=1,KSPEC2_G
        PSPECG(JM,JFLD) =PSPEC(JM,JFLD)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
  CALL GSTATS(1644,1)
ELSE
  IFLDR = 0
  IFLDS = 0

  CALL GSTATS(810,0)
  DO JFLD=1,KFDISTG

    IBSET = KVSET(JFLD)
  !Send
    IF( IBSET == MYSETV )THEN

      IFLDS = IFLDS+1
      ISND  = KTO(JFLD)
      ITAG  = MTAGDISTSP+JFLD+17
      IF(LDIM1_IS_FLD) THEN
        ZFLD(1:KSPEC2)=PSPEC(IFLDS,1:KSPEC2)
        CALL MPL_SEND(ZFLD(1:KSPEC2),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ,&
         &CDSTRING='GATH_SPEC_CONTROL')
      ELSE
        CALL MPL_SEND(PSPEC(1:KSPEC2,IFLDS),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ,&
         &CDSTRING='GATH_SPEC_CONTROL')
      ENDIF
    ENDIF


  ! RECIEVE
    IF(KTO(JFLD) == MYPROC) THEN
      IFLDR = IFLDR+1
      DO JA=1,NPRTRW
        ILEN = KPOSSP(JA+1)-KPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL SET2PE(IRCV,0,0,JA,IBSET)
          ITAG = MTAGDISTSP+JFLD+17
          ISTA = KPOSSP(JA)
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(ZRECV(ISTA:ISTP),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
           &KOUNT=ILENR,CDSTRING='GATH_SPEC_CONTROL')
          IF( ILENR /= ILEN )THEN
            CALL ABORT_TRANS('GATH_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
          ENDIF
          IF(LDIM1_IS_FLD) THEN
            II = 0
            DO JM=0,KSMAX
              DO JN=JM,KSMAX
                ISP = KDIM0G(JM)+(JN-JM)*2
                II = II+1
                PSPECG(IFLDR,II) = ZRECV(ISP)
                ISP = ISP+1
                II = II+1
                IF(JM /= 0) THEN
                  PSPECG(IFLDR,II) = ZRECV(ISP)
                ELSE
                  PSPECG(IFLDR,II) = 0.0_JPRB
                ENDIF
              ENDDO
            ENDDO
          ELSE
            II = 0
            DO JM=0,KSMAX
              DO JN=JM,KSMAX
                ISP = KDIM0G(JM)+(JN-JM)*2
                II = II+1
                PSPECG(II,IFLDR) = ZRECV(ISP)
                ISP = ISP+1
                II = II+1
                IF(JM /= 0) THEN
                  PSPECG(II,IFLDR) = ZRECV(ISP)
                ELSE
                  PSPECG(II,IFLDR) = 0.0_JPRB
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    IF( IBSET == MYSETV )THEN
      CALL MPL_WAIT(ZFLD,KREQUEST=ISENDREQ, &
       & CDSTRING='GATH_GRID_CTL: WAIT')
    ENDIF  
  ENDDO
  CALL GSTATS(810,1)
  !Synchronize processors
  CALL GSTATS(785,0)
  CALL MPL_BARRIER(CDSTRING='GATH_SPEC_CONTROL:')
  CALL GSTATS(785,1)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC_CONTROL
END MODULE GATH_SPEC_CONTROL_MOD


