SUBROUTINE DIST_SPEC(PSPECG,KFDISTG,KFROM,KVSET,KRESOL,PSPEC)

!**** *DIST_SPEC* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Interface routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL DIST__SPEC(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KFROM(:)    - Processor resposible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array
!
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  DIST_SPEC_CONTROL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET_RESOL_MOD
USE DIST_SPEC_CONTROL_MOD
USE ABORT_TRANS_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPECG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)

!ifndef INTERFACE

INTEGER_M :: IVSET(KFDISTG)
INTEGER_M :: IFSEND,IFRECV,J

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

IF(UBOUND(KFROM,1) < KFDISTG) THEN
 CALL ABORT_TRANS('DIST_SPEC: KFROM TOO SHORT!')
ENDIF
 
IFSEND = 0
IFRECV = 0

DO J=1,KFDISTG
  IF(KFROM(J) < 1 .OR. KFROM(J) > NPROC) THEN
    WRITE(NERR,*) 'DIST_SPEC:ILLEGAL KFROM VALUE',KFROM(J),J
    CALL ABORT_TRANS('DIST_SPEC:ILLEGAL KFROM VALUE')
  ENDIF
  IF(KFROM(J) == MYPROC) IFSEND = IFSEND+1
ENDDO

IF(IFSEND > 0) THEN
  IF(.NOT.PRESENT(PSPECG)) THEN
    CALL ABORT_TRANS('DIST_SPEC:PSPECG MISSING')
  ENDIF
  IF(UBOUND(PSPECG,1) < IFSEND) THEN
    CALL ABORT_TRANS('DIST_SPEC:FIRST DIMENSION OF PSPECG TOO SMALL')
  ENDIF 
 IF(UBOUND(PSPECG,2) < R%NSPEC2_G) THEN
    CALL ABORT_TRANS('DIST_SPEC:SECOND DIMENSION OF PSPECG TOO SMALL')
  ENDIF
ENDIF

IF(PRESENT(KVSET)) THEN
  IF(UBOUND(KVSET,1) < KFDISTG) THEN
    CALL ABORT_TRANS('DIST_SPEC: KVSET TOO SHORT!')
  ENDIF
  DO J=1,KFDISTG
    IF(KVSET(J) > NPRTRV .OR. KVSET(J) < 1) THEN
      WRITE(NERR,*) 'DIST_SPEC:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABORT_TRANS('DIST_SPEC:KVSET CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSET(J) == MYSETV) THEN
      IFRECV = IFRECV+1
    ENDIF
  ENDDO
  IVSET(:) = KVSET(1:KFDISTG)
ELSE
  IFRECV   = KFDISTG
  IVSET(:) = MYSETV
ENDIF

IF(IFRECV > 0 ) THEN
  IF(.NOT.PRESENT(PSPEC)) THEN
    CALL ABORT_TRANS('DIST_SPEC: FIELDS TO RECEIVE AND PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,1) < IFRECV) THEN
    CALL ABORT_TRANS('DIST_SPEC: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,2) < D%NSPEC2 ) THEN
    CALL ABORT_TRANS('DIST_SPEC: SECOND DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF

CALL DIST_SPEC_CONTROL(PSPECG,KFDISTG,KFROM,IVSET,PSPEC)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC

