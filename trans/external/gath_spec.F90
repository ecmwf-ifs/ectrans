SUBROUTINE GATH_SPEC(PSPECG,KFGATHG,KTO,KVSET,KRESOL,PSPEC)

!**** *GATH_SPEC* - Gather global spectral array from processors

!     Purpose.
!     --------
!        Interface routine for gathering spectral array

!**   Interface.
!     ----------
!     CALL GATH_SPEC(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFGATHG     - Global number of fields to be gathered
!     KTO(:)      - Processor responsible for gathering each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array
!
!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET_RESOL_MOD
USE GATH_SPEC_CONTROL_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER_M          , INTENT(IN)  :: KFGATHG
INTEGER_M          , INTENT(IN)  :: KTO(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)

!ifndef INTERFACE

INTEGER_M :: IVSET(KFGATHG)
INTEGER_M :: IFRECV,IFSEND,J

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

IF(UBOUND(KTO,1) < KFGATHG) THEN
 CALL ABOR1('GATH_SPEC: KTO TOO SHORT!')
ENDIF
 
IFRECV = 0
DO J=1,KFGATHG
  IF(KTO(J) < 1 .OR. KTO(J) > NPROC) THEN
    WRITE(NERR,*) 'GATH_SPEC:ILLEGAL KTO VALUE',KTO(J),J
    CALL ABOR1('GATH_SPEC:ILLEGAL KTO VALUE')
  ENDIF
  IF(KTO(J) == MYPROC) IFRECV = IFRECV+1
ENDDO

IF(IFRECV > 0) THEN
  IF(.NOT.PRESENT(PSPECG)) THEN
    CALL ABOR1('GATH_SPEC:PSPECG MISSING')
  ENDIF
  IF(UBOUND(PSPECG,1) < IFRECV) THEN
    CALL ABOR1('GATH_SPEC:FIRST DIMENSION OF PSPECG TOO SMALL')
  ENDIF 
 IF(UBOUND(PSPECG,2) < R%NSPEC2_G) THEN
    CALL ABOR1('GATH_SPEC:FIRST DIMENSION OF PSPECG TOO SMALL')
  ENDIF
ENDIF

IF(PRESENT(KVSET)) THEN
  IF(UBOUND(KVSET,1) < KFGATHG) THEN
    CALL ABOR1('GATH_SPEC: KVSET TOO SHORT!')
  ENDIF
  DO J=1,KFGATHG
    IF(KVSET(J) > NPRTRV .OR. KVSET(J) < 1) THEN
      WRITE(NERR,*) 'GATH_SPEC:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABOR1('GATH_SPEC:KVSET CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSET(J) == MYSETV) THEN
      IFSEND = IFSEND+1
    ENDIF
  ENDDO
  IVSET(:) = KVSET(1:KFGATHG)
ELSEIF(NPRTRV > 1) THEN
  WRITE(NERR,*) 'GATH_SPEC:KVSET MISSING, NPRTRV ',NPRTRV
  CALL ABOR1('GATH_SPEC:KVSET MISSING, NPRTRV > 1')
ELSE 
  IFSEND   = KFGATHG
  IVSET(:) = 1
ENDIF

IF(IFSEND > 0 ) THEN
  IF(.NOT.PRESENT(PSPEC)) THEN
    CALL ABOR1('GATH_SPEC: FIELDS TO RECIEVE AND PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,1) < IFSEND) THEN
    CALL ABOR1('GATH_SPEC: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,2) < D%NSPEC2 ) THEN
    CALL ABOR1('GATH_SPEC: SECOND DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF

CALL GATH_SPEC_CONTROL(PSPECG,KFGATHG,KTO,IVSET,PSPEC)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC

