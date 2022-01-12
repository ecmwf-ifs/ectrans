! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIST_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE DIST_SPEC_CONTROL(PSPECG,KFDISTG,KFROM,KVSET,PSPEC,LDIM1_IS_FLD,&
 & KSMAX,KSPEC2,KSPEC2_G,KPOSSP,KDIM0G,KSORT)

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
!     KSORT(:)   - Re-order fields on output

!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01
!    P.Marguinaud : 2014-10-10

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_BARRIER, MPL_WAIT, &
     &                  JP_NON_BLOCKING_STANDARD

!USE TPM_GEN
!USE TPM_DIM
USE TPM_DISTR       ,ONLY : MTAGDISTSP, MYSETV, NPRCIDS, NPRTRW, MYPROC, NPROC

USE SET2PE_MOD      ,ONLY : SET2PE
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

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
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), TARGET :: KSORT (:)
    
INTEGER(KIND=JPIM) :: IDIST(KSPEC2_G)
REAL(KIND=JPRB)    :: ZFLD(KSPEC2)
REAL(KIND=JPRB),ALLOCATABLE  :: ZBUF(:,:)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,JNM,IBSET,ILEN,JA,ISND
INTEGER(KIND=JPIM) :: IRCV,ISTA,ISTP,ILENR,ISENDREQ(NPRTRW*KFDISTG)
INTEGER(KIND=JPIM) :: ISMAX, ISPEC2, IPOS0,ISENT
INTEGER(KIND=JPIM), POINTER :: ISORT (:)

!     ------------------------------------------------------------------


! Compute help array for distribution

IF (PRESENT (KSORT)) THEN
  ISORT => KSORT
ELSE
  ALLOCATE (ISORT (KFDISTG))
  DO JFLD = 1, KFDISTG
    ISORT (JFLD) = JFLD
  ENDDO
ENDIF

IF( NPROC == 1 ) THEN
  CALL GSTATS(1644,0)
  IF(LDIM1_IS_FLD) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JNM,JFLD)
    DO JNM=1,KSPEC2_G
      DO JFLD=1,KFDISTG
        PSPEC(ISORT (JFLD),JNM) = PSPECG(JFLD,JNM)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JNM,JFLD)
    DO JFLD=1,KFDISTG
      DO JNM=1,KSPEC2_G
        PSPEC(JNM,ISORT (JFLD)) = PSPECG(JNM,JFLD)
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

  IFLDS = 0
  DO JFLD=1,KFDISTG
    IF(KFROM(JFLD) == MYPROC) THEN
      IFLDS = IFLDS+1
    ENDIF
  ENDDO
  ALLOCATE(ZBUF(KSPEC2_G,IFLDS))

  CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JNM,JFLD)
  DO JFLD=1,IFLDS
    IF(LDIM1_IS_FLD) THEN
      DO JNM=1,KSPEC2_G
        ZBUF(IDIST(JNM),JFLD) = PSPECG(JFLD,JNM)
      ENDDO
    ELSE
      DO JNM=1,KSPEC2_G
        ZBUF(IDIST(JNM),JFLD) = PSPECG(JNM,JFLD)
      ENDDO
    ENDIF
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1644,1)

  IFLDR = 0
  IFLDS = 0
  ISENT = 0

  CALL GSTATS_BARRIER(790)
  CALL GSTATS(812,0)
  DO JFLD=1,KFDISTG
    
  ! Send
    IF(KFROM(JFLD) == MYPROC) THEN
      IFLDS = IFLDS+1
      IBSET = KVSET(JFLD)
      ITAG  = MTAGDISTSP+JFLD

      DO JA=1,NPRTRW
        ILEN = KPOSSP(JA+1)-KPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL SET2PE(ISND,0,0,JA,IBSET)
          ISTA  = KPOSSP(JA)
          ISTP  = ISTA+ILEN-1
          ISENT = ISENT+1
          CALL MPL_SEND(ZBUF(ISTA:ISTP,IFLDS),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
           &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISENT),&
           &CDSTRING='DIST_SPEC_CONTROL:')
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  !Recieve
  DO JFLD=1,KFDISTG
    IBSET = KVSET(JFLD)
    IF( IBSET == MYSETV )THEN
      ITAG = MTAGDISTSP+JFLD
      IF( KSPEC2 > 0 )THEN
        IRCV  = KFROM(JFLD)
        IFLDR = IFLDR+1
        IF(LDIM1_IS_FLD) THEN
          CALL MPL_RECV(ZFLD,KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
           &KOUNT=ILENR,CDSTRING='DIST_SPEC_CONTROL:')
          PSPEC(ISORT (IFLDR),1:KSPEC2) = ZFLD(:)
        ELSE
          CALL MPL_RECV(PSPEC(:,ISORT (IFLDR)),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
           &KOUNT=ILENR,CDSTRING='DIST_SPEC_CONTROL:')
        ENDIF
        IF( ILENR /= KSPEC2 )THEN
          CALL ABORT_TRANS('DIST_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  DO JA=1,ISENT
    CALL MPL_WAIT(KREQUEST=ISENDREQ(JA), &
     & CDSTRING='DIST_SPEC_CTL: WAIT')
  ENDDO

  CALL GSTATS(812,1)
  CALL GSTATS_BARRIER2(790)

!Synchronize processors
  CALL GSTATS(787,0)
  IF( NPROC > 1 )THEN
    CALL MPL_BARRIER(CDSTRING='DIST_SPEC_CONTROL:')
  ENDIF
  CALL GSTATS(787,1)
  IF(ALLOCATED(ZBUF)) DEALLOCATE(ZBUF)
ENDIF

IF (.NOT. PRESENT (KSORT)) THEN
  DEALLOCATE (ISORT)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC_CONTROL
END MODULE DIST_SPEC_CONTROL_MOD
