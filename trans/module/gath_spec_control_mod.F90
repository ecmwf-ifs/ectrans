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
!
!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET2PE_MOD

IMPLICIT NONE

REAL_B    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KTO(:)
INTEGER_M          , INTENT(IN)  :: KVSET(:)
REAL_B    ,OPTIONAL, INTENT(IN) :: PSPEC(:,:)

REAL_B    :: ZFLD(R%NSPEC2_G)
INTEGER_M :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,JNM,IBSET,ILEN,JA,IERR,ISND
INTEGER_M :: ISENDER,ITAGR,IRCV,ISP,ILENR,IPROCA

!     ------------------------------------------------------------------


!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
  PSPECG(1:KFDISTG,1:R%NSPEC2_G) =PSPEC(1:KFDISTG,1:R%NSPEC2_G) 
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
      CALL MPE_SEND(ZFLD,D%NSPEC2,MREALT,NPRCIDS(ISND),ITAG,0,&
       &0,0,IERR)
      IF( IERR < 0 )THEN
        CALL ABOR1(' GATH_SPEC_CONTROL : ERROR IN MPE_RECV (ZFLD)')
      ENDIF
    ENDIF


  ! RECIEVE
    IF(KTO(JFLD) == MYPROC) THEN
      IFLDR = IFLDR+1
      DO JA=1,NPRTRW
        ILEN = D%NPOSSP(JA+1)-D%NPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL SET2PE(IRCV,0,0,JA,IBSET)
          ITAG = MTAGDISTSP+JFLD+17
          CALL MPE_RECV(ZFLD(D%NPOSSP(JA)),ILEN,MREALT,&
           &NPRCIDS(IRCV),ITAG,0,0,0,ILENR,ISENDER,ITAGR,IERR)
          IF( IERR < 0 )THEN
            CALL ABOR1(' GATH_SPEC_CONTROL : ERROR IN MPE_SEND (ZFLD)')
          ENDIF
          IF( ILENR /= ILEN )THEN
            CALL ABOR1(' GATH_SPEC_CONTROL: INVALID RECEIVE MESSAGE LENGTH')
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
                PSPECG(IFLDR,II) = _ZERO_
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF

  !Synchronize processors
    CALL MPE_BARRIER(IERR)
    IF( IERR /= 0 )THEN
      CALL ABOR1(' GATH_SPEC_CONTROL: ERROR IN MPE_BARRIER')
    ENDIF

  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC_CONTROL
END MODULE GATH_SPEC_CONTROL_MOD


