MODULE DIST_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE DIST_SPEC_CONTROL(PSPECG,KFDISTG,KFROM,KVSET,PSPEC)

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

INTEGER(KIND=JPIM) :: IDIST(R%NSPEC2_G)
REAL(KIND=JPRB)    :: ZFLD(R%NSPEC2_G)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,JNM,IBSET,ILEN,JA,ISND
INTEGER(KIND=JPIM) :: IRCV,ISTA,ISTP,ILENR

!     ------------------------------------------------------------------

! Compute help array for distribution
II = 0
DO JM=0,R%NSMAX
  DO JN=JM,R%NSMAX
    IDIST(II+1) = D%NDIM0G(JM)+(JN-JM)*2
    IDIST(II+2) = D%NDIM0G(JM)+(JN-JM)*2+1
    II = II+2
  ENDDO
ENDDO

!Distribute spectral array

IFLDR = 0
IFLDS = 0

DO JFLD=1,KFDISTG

  IBSET = KVSET(JFLD)

  ! Send
  IF(KFROM(JFLD) == MYPROC) THEN
    IFLDS = IFLDS+1
    DO JNM=1,R%NSPEC2_G
      ZFLD(IDIST(JNM)) = PSPECG(IFLDS,JNM) 
    ENDDO
    DO JA=1,NPRTRW
      ILEN = D%NPOSSP(JA+1)-D%NPOSSP(JA)
      IF( ILEN > 0 )THEN
        CALL SET2PE(ISND,0,0,JA,IBSET)
        IF( NPROC > 1 )THEN
          ITAG = MTAGDISTSP
          ISTA = D%NPOSSP(JA)
          ISTP = ISTA+ILEN-1
          CALL MPL_SEND(ZFLD(ISTA:ISTP),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
           &CDSTRING='DIST_SPEC_CONTROL:')
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  !Recieve
  IF( IBSET == MYSETV )THEN

    IF( D%NSPEC2 > 0 )THEN
      IF( NPROC > 1 )THEN
        IRCV = KFROM(JFLD)
        ITAG = MTAGDISTSP
        CALL MPL_RECV(ZFLD(1:D%NSPEC2),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
         &KOUNT=ILENR,CDSTRING='DIST_SPEC_CONTROL:')
        IF( ILENR /= D%NSPEC2 )THEN
          CALL ABORT_TRANS('DIST_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
        ENDIF
      ENDIF
    ENDIF
    IFLDR = IFLDR+1
    PSPEC(IFLDR,:) = ZFLD(1:D%NSPEC2)
  ENDIF

  !Synchronize processors
  IF( NPROC > 1 )THEN
    CALL MPL_BARRIER(CDSTRING='DIST_SPEC_CONTROL:')
  ENDIF

ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC_CONTROL
END MODULE DIST_SPEC_CONTROL_MOD


