! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


SUBROUTINE EDIST_SPEC(PSPECG,KFDISTG,KFROM,KVSET,KRESOL,PSPEC,&
 & LDIM1_IS_FLD,KSORT)

!**** *EDIST_SPEC* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Interface routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL EDIST__SPEC(...)

!     Explicit arguments :
!     --------------------
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KFROM(:)    - Processor resposible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array

!     Method.
!     -------

!     Externals.  ESET_RESOL   - set resolution
!     ----------  DIST_SPEC_CONTROL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Marguinaud  10-Oct-2014 Add KSORT argument (change the order of fields)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR
USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD
USE TPM_DISTR       ,ONLY : D, NPRTRV, NPRTRW, MYSETV, MYPROC, NPROC
USE TPMALD_DISTR    ,ONLY : DALD

USE ESET_RESOL_MOD  ,ONLY : ESET_RESOL
USE DIST_SPEC_CONTROL_MOD ,ONLY : DIST_SPEC_CONTROL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PSPECG(:,:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDISTG
INTEGER(KIND=JPIM),INTENT(IN)    :: KFROM(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KVSET(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KRESOL
LOGICAL           ,OPTIONAL,INTENT(IN)    :: LDIM1_IS_FLD
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT)   :: PSPEC(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KSORT (:)
!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IVSET(KFDISTG)
INTEGER(KIND=JPIM) :: IFSEND,IFRECV,J, IFLD, ICOEFF
INTEGER(KIND=JPIM) :: ISMAX, ISPEC2, ISPEC2_G, ISPEC2MX
INTEGER(KIND=JPIM) :: IPOSSP(NPRTRW+1)
INTEGER(KIND=JPIM) :: IUMPP(NPRTRW)
INTEGER(KIND=JPIM) :: IPTRMS(NPRTRW)
INTEGER(KIND=JPIM),ALLOCATABLE :: IDIM0G(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IALLMS(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IKN(:)
LOGICAL :: LLDIM1_IS_FLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

! Set current resolution
IF (LHOOK) CALL DR_HOOK('EDIST_SPEC',0,ZHOOK_HANDLE)
CALL ESET_RESOL(KRESOL)

LLDIM1_IS_FLD=.TRUE.
IF(PRESENT(LDIM1_IS_FLD)) LLDIM1_IS_FLD=LDIM1_IS_FLD
IF(LLDIM1_IS_FLD) THEN
  IFLD=1
  ICOEFF=2
ELSE
  IFLD=2
  ICOEFF=1
ENDIF

ISMAX = RALD%NMSMAX
ALLOCATE(IDIM0G(0:ISMAX))
ALLOCATE(IALLMS(ISMAX+1))
ALLOCATE(IKN(0:ISMAX))
ISPEC2    = D%NSPEC2
ISPEC2_G  = R%NSPEC2_G
IPOSSP(:) = D%NPOSSP(:)
IDIM0G(:) = D%NDIM0G(:)
ISPEC2MX  = D%NSPEC2MX
IUMPP(:)  = D%NUMPP(:)
IALLMS(:) = D%NALLMS(:)
IPTRMS(:) = D%NPTRMS(:)
DO J=0,ISMAX
  IKN(J)=2*DALD%NCPL2M(J)
ENDDO

IF(UBOUND(KFROM,1) < KFDISTG) THEN
  CALL ABORT_TRANS('EDIST_SPEC: KFROM TOO SHORT!')
ENDIF
 
IFSEND = 0
IFRECV = 0

DO J=1,KFDISTG
  IF(KFROM(J) < 1 .OR. KFROM(J) > NPROC) THEN
    WRITE(NERR,*) 'EDIST_SPEC:ILLEGAL KFROM VALUE',KFROM(J),J
    CALL ABORT_TRANS('EDIST_SPEC:ILLEGAL KFROM VALUE')
  ENDIF
  IF(KFROM(J) == MYPROC) IFSEND = IFSEND+1
ENDDO

IF(IFSEND > 0) THEN
  IF(.NOT.PRESENT(PSPECG)) THEN
    CALL ABORT_TRANS('EDIST_SPEC:PSPECG MISSING')
  ENDIF
  IF(UBOUND(PSPECG,IFLD) < IFSEND) THEN
    WRITE(NERR,*)'EDIST_SPEC: ',IFLD, UBOUND(PSPECG,IFLD), IFSEND
    CALL ABORT_TRANS('EDIST_SPEC:FIELD DIMENSION OF PSPECG TOO SMALL')
  ENDIF
  IF(UBOUND(PSPECG,ICOEFF) < ISPEC2_G) THEN
    WRITE(NERR,*)'EDIST_SPEC: ',ICOEFF, UBOUND(PSPECG,ICOEFF), ISPEC2_G
    CALL ABORT_TRANS('EDIST_SPEC: COEFF DIMENSION OF PSPECG TOO SMALL')
  ENDIF
ENDIF

IF(PRESENT(KVSET)) THEN
  IF(UBOUND(KVSET,1) < KFDISTG) THEN
    CALL ABORT_TRANS('EDIST_SPEC: KVSET TOO SHORT!')
  ENDIF
  DO J=1,KFDISTG
    IF(KVSET(J) > NPRTRV .OR. KVSET(J) < 1) THEN
      WRITE(NERR,*) 'EDIST_SPEC:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABORT_TRANS('EDIST_SPEC:KVSET CONTAINS VALUES OUTSIDE RANGE')
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
    CALL ABORT_TRANS('EDIST_SPEC: FIELDS TO RECEIVE AND PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,IFLD) < IFRECV) THEN
    CALL ABORT_TRANS('EDIST_SPEC: FIELD DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,ICOEFF) < ISPEC2 ) THEN
    CALL ABORT_TRANS('EDIST_SPEC: COEFF DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF

IF (PRESENT (KSORT)) THEN
  IF (.NOT. PRESENT (PSPEC)) THEN
    CALL ABORT_TRANS('EDIST_SPEC: KSORT REQUIRES PSPEC')
  ENDIF
  IF (UBOUND (KSORT, 1) /= UBOUND (PSPEC, IFLD)) THEN
    CALL ABORT_TRANS('EDIST_SPEC: DIMENSION MISMATCH KSORT, PSPEC')
  ENDIF
ENDIF

CALL DIST_SPEC_CONTROL(PSPECG,KFDISTG,KFROM,IVSET,PSPEC,LLDIM1_IS_FLD,&
 & ISMAX,ISPEC2,ISPEC2MX,ISPEC2_G,IPOSSP,IDIM0G,IUMPP,IALLMS,IPTRMS,IKN,KSORT)
DEALLOCATE(IDIM0G)
IF (LHOOK) CALL DR_HOOK('EDIST_SPEC',1,ZHOOK_HANDLE)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE EDIST_SPEC

