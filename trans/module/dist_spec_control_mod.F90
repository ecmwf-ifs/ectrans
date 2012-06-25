MODULE DIST_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE DIST_SPEC_CONTROL(PSPECG,KFDISTG,KFROM,KVSET,PSPEC,LDIM1_IS_FLD,&
 & KSMAX,KSPEC2,KSPEC2_G,KPOSSP,KDIM0G)

!**** *DIST_SPEC_CONTROL* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL DIST_SPEC_CONTROL(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KFROM(:)    - Processor resposible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     PSPEC(:,:)  - Local spectral array

!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET2PE_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2_G
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPOSSP(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KDIM0G(0:)
    
INTEGER(KIND=JPIM) :: IDIST(KSPEC2_G)
REAL(KIND=JPRB)    :: ZFLD(KSPEC2_G),ZBUF(KSPEC2_G)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,JNM,IBSET,ILEN,JA,ISND
INTEGER(KIND=JPIM) :: IRCV,ISTA,ISTP,ILENR,ISENDREQ(NPRTRW)
INTEGER(KIND=JPIM) :: ISMAX, ISPEC2, IPOS0

!     ------------------------------------------------------------------


! Compute help array for distribution


IF( NPROC == 1 ) THEN
  CALL GSTATS(1644,0)
  IF(LDIM1_IS_FLD) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
    DO JM=1,KSPEC2_G
      DO JFLD=1,KFDISTG
        PSPEC(JFLD,JM) = PSPECG(JFLD,JM)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
    DO JFLD=1,KFDISTG
      DO JM=1,KSPEC2_G
        PSPEC(JM,JFLD) = PSPECG(JM,JFLD)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
  CALL GSTATS(1644,1)
ELSE
  II = 0
  CALL GSTATS(1804,0)
  DO JM=0,KSMAX
    DO JN=JM,KSMAX
      IDIST(II+1) = KDIM0G(JM)+(JN-JM)*2
      IDIST(II+2) = KDIM0G(JM)+(JN-JM)*2+1
      II = II+2
    ENDDO
  ENDDO
  CALL GSTATS(1804,1)

!Distribute spectral array

  IFLDR = 0
  IFLDS = 0

  CALL GSTATS_BARRIER(790)
  CALL GSTATS(812,0)
  DO JFLD=1,KFDISTG
    
    IBSET = KVSET(JFLD)
    ITAG = MTAGDISTSP+JFLD
  ! Send
    IF(KFROM(JFLD) == MYPROC) THEN
      IFLDS = IFLDS+1
      IF(LDIM1_IS_FLD) THEN
        DO JNM=1,KSPEC2_G
          ZBUF(IDIST(JNM)) = PSPECG(IFLDS,JNM) 
        ENDDO
      ELSE
        DO JNM=1,KSPEC2_G
          ZBUF(IDIST(JNM)) = PSPECG(JNM,IFLDS) 
        ENDDO
      ENDIF
      DO JA=1,NPRTRW
        ILEN = KPOSSP(JA+1)-KPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL SET2PE(ISND,0,0,JA,IBSET)
          ISTA = KPOSSP(JA)
          ISTP = ISTA+ILEN-1
          CALL MPL_SEND(ZBUF(ISTA:ISTP),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
           &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JA),&
           &CDSTRING='DIST_SPEC_CONTROL:')
        ENDIF
      ENDDO
    ENDIF

  !Recieve
    IF( IBSET == MYSETV )THEN

      IF( KSPEC2 > 0 )THEN
        IRCV = KFROM(JFLD)
        CALL MPL_RECV(ZFLD(1:KSPEC2),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
         &KOUNT=ILENR,CDSTRING='DIST_SPEC_CONTROL:')
        IF( ILENR /= KSPEC2 )THEN
          CALL ABORT_TRANS('DIST_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
        ENDIF
        IFLDR = IFLDR+1
        IF(LDIM1_IS_FLD) THEN
          PSPEC(IFLDR,:) = ZFLD(1:KSPEC2)
        ELSE
          PSPEC(:,IFLDR) = ZFLD(1:KSPEC2)
        ENDIF
      ENDIF
    ENDIF

    IF(KFROM(JFLD) == MYPROC) THEN
      DO JA=1,NPRTRW
        ILEN = KPOSSP(JA+1)-KPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL MPL_WAIT(ZBUF,KREQUEST=ISENDREQ(JA), &
           & CDSTRING='DIST_SPEC_CTL: WAIT')
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  CALL GSTATS(812,1)
  CALL GSTATS_BARRIER2(790)

!Synchronize processors
  CALL GSTATS(787,0)
  IF( NPROC > 1 )THEN
    CALL MPL_BARRIER(CDSTRING='DIST_SPEC_CONTROL:')
  ENDIF
  CALL GSTATS(787,1)
ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC_CONTROL
END MODULE DIST_SPEC_CONTROL_MOD


