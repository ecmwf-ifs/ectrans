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
!
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET2PE_MOD

IMPLICIT NONE

REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPECG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
INTEGER_M          , INTENT(IN)  :: KVSET(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)

INTEGER_M :: IDIST(R%NSPEC2_G)
REAL_B    :: ZFLD(R%NSPEC2_G)
INTEGER_M :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,JNM,IBSET,ILEN,JA,IERR,ISND
INTEGER_M :: ISENDER,ITAGR,IRCV

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
        ITAG = MTAGDISTSP+JFLD
        IF( NPROC > 1 )THEN
          CALL MPE_SEND(ZFLD(D%NPOSSP(JA)),ILEN,MREALT,&
           &NPRCIDS(ISND),ITAG,0,0,0,IERR)
          IF( IERR < 0 )THEN
            CALL ABOR1(' DIST_SPEC_CONTROL : ERROR IN MPE_SEND (ZFLD)')
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  !Recieve
  IF( IBSET == MYSETV )THEN

    IF( D%NSPEC2 > 0 )THEN
      IRCV = KFROM(JFLD)
      ITAG = MTAGDISTSP+JFLD
      IF( NPROC > 1 )THEN
        CALL MPE_RECV(ZFLD,D%NSPEC2,MREALT,NPRCIDS(IRCV),ITAG,0,&
         &0,0,ILEN,ISENDER,ITAGR,IERR)
        IF( IERR < 0 )THEN
          CALL ABOR1(' DIST_SPEC_CONTROL : ERROR IN MPE_RECV (ZFLD)')
        ENDIF
        IF( ILEN /= D%NSPEC2 )THEN
          CALL ABOR1(' DIST_SPEC_CONTROL: INVALID RECEIVE MESSAGE LENGTH')
        ENDIF
      ENDIF
    ENDIF
    IFLDR = IFLDR+1
    PSPEC(IFLDR,:) = ZFLD(1:D%NSPEC2)
  ENDIF

  !Synchronize processors
  IF( NPROC > 1 )THEN
    CALL MPE_BARRIER(IERR)
    IF( IERR /= 0 )THEN
      CALL ABOR1(' DIST_SPEC_CONTROL: ERROR IN MPE_BARRIER')
    ENDIF
  ENDIF

ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC_CONTROL
END MODULE DIST_SPEC_CONTROL_MOD


