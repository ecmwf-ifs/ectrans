! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GATH_SPEC(PSPECG,KFGATHG,KTO,KVSET,KRESOL,PSPEC,LDIM1_IS_FLD,KSMAX,LDZA0IP)

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
!     LDZA0IP     - Set to zero imaginary part of first coefficients
!
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  GATH_SPEC_CONTROL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        Modified 03-09-30  Y. Seity, bug correction IFSEND=0
!        Modified 13-10-10  P. Marguinaud add LDZA0IP option
!     ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB

!ifndef INTERFACE

USE TPM_GEN,               ONLY: NERR
USE TPM_DIM,               ONLY: R
USE TPM_DISTR,             ONLY: D, NPRTRV, NPRTRW, MYSETV, MYSETW, MYPROC, NPROC
USE SET_RESOL_MOD,         ONLY: SET_RESOL
USE GATH_SPEC_CONTROL_MOD, ONLY: GATH_SPEC_CONTROL
USE SUWAVEDI_MOD,          ONLY: SUWAVEDI
USE ABORT_TRANS_MOD,       ONLY: ABORT_TRANS
USE YOMHOOK,               ONLY: LHOOK, DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSMAX
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDZA0IP

!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IVSET(KFGATHG)
INTEGER(KIND=JPIM) :: IFRECV,IFSEND,J
INTEGER(KIND=JPIM) :: IFLD,ICOEFF
INTEGER(KIND=JPIM) :: ISMAX, ISPEC2, ISPEC2_G
INTEGER(KIND=JPIM) :: IPOSSP(NPRTRW+1)
INTEGER(KIND=JPIM),ALLOCATABLE :: IDIM0G(:)
LOGICAL :: LLDIM1_IS_FLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GATH_SPEC',0,ZHOOK_HANDLE)
! Set current resolution
CALL SET_RESOL(KRESOL)

LLDIM1_IS_FLD = .TRUE.
IF(PRESENT(LDIM1_IS_FLD)) LLDIM1_IS_FLD = LDIM1_IS_FLD

IF(LLDIM1_IS_FLD) THEN
  IFLD = 1
  ICOEFF = 2
ELSE
  IFLD = 2
  ICOEFF = 1
ENDIF
IF(UBOUND(KTO,1) < KFGATHG) THEN
 CALL ABORT_TRANS('GATH_SPEC: KTO TOO SHORT!')
ENDIF

ISMAX = R%NSMAX
IF(PRESENT(KSMAX)) ISMAX = KSMAX
ALLOCATE(IDIM0G(0:ISMAX))
IF(ISMAX /= R%NSMAX) THEN
  CALL SUWAVEDI(ISMAX,ISMAX,NPRTRW,MYSETW,KPOSSP=IPOSSP,KSPEC2=ISPEC2,&
   & KDIM0G=IDIM0G)
  ISPEC2_G = (ISMAX+1)*(ISMAX+2)
ELSE
  ISPEC2    = D%NSPEC2
  ISPEC2_G  = R%NSPEC2_G
  IPOSSP(:) = D%NPOSSP(:)
  IDIM0G(:) = D%NDIM0G(:)
ENDIF

IFSEND = 0
IFRECV = 0
DO J=1,KFGATHG
  IF(KTO(J) < 1 .OR. KTO(J) > NPROC) THEN
    WRITE(NERR,*) 'GATH_SPEC:ILLEGAL KTO VALUE',KTO(J),J
    CALL ABORT_TRANS('GATH_SPEC:ILLEGAL KTO VALUE')
  ENDIF
  IF(KTO(J) == MYPROC) IFRECV = IFRECV+1
ENDDO

IF(IFRECV > 0) THEN
  IF(.NOT.PRESENT(PSPECG)) THEN
    CALL ABORT_TRANS('GATH_SPEC:PSPECG MISSING')
  ENDIF
  IF(UBOUND(PSPECG,IFLD) < IFRECV) THEN
    WRITE(NERR,*) 'GATH_SPEC: ',IFLD, UBOUND(PSPECG,IFLD), IFRECV
    CALL ABORT_TRANS('GATH_SPEC:FIELD DIMENSION OF PSPECG TOO SMALL')
  ENDIF
 IF(UBOUND(PSPECG,ICOEFF) < ISPEC2_G) THEN
    WRITE(NERR,*) 'GATH_SPEC: ',ICOEFF, UBOUND(PSPECG,ICOEFF), ISPEC2_G
    CALL ABORT_TRANS('GATH_SPEC:COEFF DIMENSION OF PSPECG TOO SMALL')
  ENDIF
ENDIF

IF(PRESENT(KVSET)) THEN
  IF(UBOUND(KVSET,1) < KFGATHG) THEN
    CALL ABORT_TRANS('GATH_SPEC: KVSET TOO SHORT!')
  ENDIF
  DO J=1,KFGATHG
    IF(KVSET(J) > NPRTRV .OR. KVSET(J) < 1) THEN
      WRITE(NERR,*) 'GATH_SPEC:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABORT_TRANS('GATH_SPEC:KVSET CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSET(J) == MYSETV) THEN
      IFSEND = IFSEND+1
    ENDIF
  ENDDO
  IVSET(:) = KVSET(1:KFGATHG)
ELSEIF(NPRTRV > 1) THEN
  WRITE(NERR,*) 'GATH_SPEC:KVSET MISSING, NPRTRV ',NPRTRV
  CALL ABORT_TRANS('GATH_SPEC:KVSET MISSING, NPRTRV > 1')
ELSE
  IFSEND   = KFGATHG
  IVSET(:) = 1
ENDIF

IF(IFSEND > 0 ) THEN
  IF(.NOT.PRESENT(PSPEC)) THEN
    CALL ABORT_TRANS('GATH_SPEC: FIELDS TO RECIEVE AND PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,IFLD) < IFSEND) THEN
    CALL ABORT_TRANS('GATH_SPEC: FIELD DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,ICOEFF) < ISPEC2 ) THEN
    CALL ABORT_TRANS('GATH_SPEC: COEFF DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF

CALL GATH_SPEC_CONTROL(PSPECG,KFGATHG,KTO,IVSET,PSPEC,LLDIM1_IS_FLD,&
 & ISMAX,ISPEC2,ISPEC2_G,IPOSSP,IDIM0G,LDZA0IP)
DEALLOCATE(IDIM0G)

IF (LHOOK) CALL DR_HOOK('GATH_SPEC',1,ZHOOK_HANDLE)
!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC

