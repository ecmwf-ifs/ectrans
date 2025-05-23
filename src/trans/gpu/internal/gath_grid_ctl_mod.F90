! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE GATH_GRID_CTL_MOD
CONTAINS
SUBROUTINE GATH_GRID_CTL(PGPG,KFGATHG,KPROMA,KTO,PGP)

!**** *GATH_GRID_CTL* - Gather global gridpoint array from processors

!     Purpose.
!     --------
!        Routine for gathering gridpoint array

!**   Interface.
!     ----------
!     CALL GATH_GRID_CTL(...)

!     Explicit arguments :
!     --------------------
!     PGPG(:,:)   - Global gridpoint array
!     KFGATHG     - Global number of fields to be gathered
!     KPROMA      - blocking factor for gridpoint input
!     KTO(:)      - Processor responsible for gathering each field
!     PGP(:,:,:)  - Local gridpoint array
!
!     ------------------------------------------------------------------


USE PARKIND1,       ONLY: JPIM, JPRB
USE MPL_MODULE,     ONLY: MPL_ALLTOALLV, MPL_RECV, MPL_SEND, MPL_WAIT, JP_BLOCKING_STANDARD, &
  &                       JP_NON_BLOCKING_STANDARD
USE TPM_GEOMETRY,   ONLY: G
USE TPM_DISTR,      ONLY: D, MTAGDISTSP, NPRCIDS, MYPROC, NPROC
USE SET2PE_MOD,     ONLY: SET2PE
USE EQ_REGIONS_MOD, ONLY: N_REGIONS, N_REGIONS_NS

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
REAL(KIND=JPRB)             , INTENT(IN)  :: PGP(:,:,:)

! Declaration of local variables

REAL(KIND=JPRB)    :: ZFLD(D%NGPTOTMX*KFGATHG),ZDUM(D%NGPTOTMX)
REAL(KIND=JPRB),ALLOCATABLE :: ZBUF(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IREQ(:)
INTEGER(KIND=JPIM) :: IFLDR,JFLD,ITAG,ILEN,JA,JB,ISND,JGL,JLON,ILOFF,ILENB
INTEGER(KIND=JPIM) :: IRCV,IOFF,ILAST,IGL1,IGL2,IGLOFF,IR
INTEGER(KIND=JPIM) :: JKGLO,JROF,IEND,J,IBL,IPROC,JROC,IMYFIELDS,ILRECV
INTEGER(KIND=JPIM) :: ISENDREQ(KFGATHG),ITO
INTEGER(KIND=JPIM) :: ILENS(NPROC),IOFFS(NPROC),ILENR(NPROC),IOFFR(NPROC)
INTEGER(KIND=JPIM) :: IFLDL,IFLDS
LOGICAL :: LLSAME
!     ------------------------------------------------------------------


!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
  CALL GSTATS(1643,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IEND,IOFF,IBL,JFLD,JROF)
  DO JKGLO=1,D%NGPTOT,KPROMA
    IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
    IOFF = JKGLO-1
    IBL  = (JKGLO-1)/KPROMA+1
    DO JFLD=1,KFGATHG
      DO JROF=1,IEND
        PGPG(IOFF+JROF,JFLD) = PGP(JROF,JFLD,IBL)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1643,1)

ELSE
! test if values in KTO are all the same
  LLSAME=.TRUE.
  ITO=KTO(1)
  DO JFLD=2,KFGATHG
    IF(KTO(JFLD) /= ITO) THEN
      LLSAME=.FALSE.
      EXIT
    ENDIF
  ENDDO

  IFLDL=D%NGPTOTMX
  IF(LLSAME) THEN
    CALL GSTATS(1643,0)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IEND,IOFF,IBL,JFLD,JROF)
    DO JFLD=1,KFGATHG
      DO JKGLO=1,D%NGPTOT,KPROMA
        IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
        IOFF = JKGLO-1
        IBL  = (JKGLO-1)/KPROMA+1
        DO JROF=1,IEND
          ZFLD(IOFF+JROF+(JFLD-1)*IFLDL) = PGP(JROF,JFLD,IBL)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    CALL GSTATS(1643,1)
  ELSE
    ILENS(:)=0
    IOFFS(:)=0
    ILENR(:)=0
    IOFFR(:)=0
    DO JFLD=1,KFGATHG
      ILENS(KTO(JFLD))=ILENS(KTO(JFLD))+IFLDL
      IF(KTO(JFLD) == MYPROC) THEN
        ILENR(:)=ILENR(:)+IFLDL
      ENDIF
    ENDDO
    DO JROC=2,NPROC
      IOFFR(JROC)=IOFFR(JROC-1)+ ILENR(JROC-1)
      IOFFS(JROC)=IOFFS(JROC-1)+ ILENS(JROC-1)
    ENDDO
    IFLDS=0
    DO JROC=1,NPROC
      IF(ILENS(JROC) > 0) THEN
        DO JFLD=1,KFGATHG
          IF(KTO(JFLD) == JROC) THEN
            DO JKGLO=1,D%NGPTOT,KPROMA
              IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
              IOFF = JKGLO-1
              IBL  = (JKGLO-1)/KPROMA+1
              DO JROF=1,IEND
                ZFLD(IOFF+JROF+IFLDS*IFLDL) = PGP(JROF,JFLD,IBL)
              ENDDO
            ENDDO
            IFLDS=IFLDS+1
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
          
  IMYFIELDS = 0
  DO JFLD=1,KFGATHG
    IF(KTO(JFLD) == MYPROC) THEN
      IMYFIELDS = IMYFIELDS+1
    ENDIF
  ENDDO

  IF(IMYFIELDS > 0) THEN
    ALLOCATE(ZBUF(D%NGPTOTMX*IMYFIELDS*NPROC))
  ELSE
    ALLOCATE(ZBUF(1))
  ENDIF
  IFLDR = 0
  CALL GSTATS_BARRIER(789)
  CALL GSTATS(809,0)

  IF( LLSAME )THEN
    !Send
    ISND  = KTO(1)
    ITAG  = MTAGDISTSP+1+17
    CALL MPL_SEND(ZFLD,KDEST=NPRCIDS(ISND),KTAG=ITAG,&
     &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(1),&
     &CDSTRING='GATH_GRID_CTL:')
    ! RECIEVE
    IF(KTO(1) == MYPROC) THEN
      IFLDR = KFGATHG
      DO JROC=1,NPROC
        ITAG  = MTAGDISTSP+1+17
        IRCV  = JROC
        IOFF=IFLDL*KFGATHG*(JROC-1)
        CALL MPL_RECV(ZBUF(IOFF+1:IOFF+IFLDL*KFGATHG),KSOURCE=NPRCIDS(IRCV),&
         &KMP_TYPE=JP_BLOCKING_STANDARD,KOUNT=ILRECV,&
         &KTAG=ITAG,CDSTRING='GATH_GRID_CTL:')
      ENDDO
    ENDIF
    CALL MPL_WAIT(KREQUEST=ISENDREQ(1), &
     & CDSTRING='GATH_GRID_CTL: WAIT')
  ELSE
    IFLDR=IMYFIELDS

!   ALLTOALLV performance is really slow when number of fields (KFGATHG) is << NPROC
!   This was for IBM - and RECV/SEND alternative causes problems for large number of MPI tasks.

!   IF( KFGATHG >= NPROC/8 )THEN
    IF( .TRUE. )THEN
      CALL MPL_ALLTOALLV(PSENDBUF=ZFLD,KSENDCOUNTS=ILENS,&
       & PRECVBUF=ZBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
       & CDSTRING='GATH_GRID_CTL:')
    ELSE
      IR=0
      DO JFLD=1,KFGATHG
        IF(KTO(JFLD) == MYPROC) THEN
          IR=IR+NPROC
        ENDIF
      ENDDO
      IR=IR+KFGATHG
      ALLOCATE(IREQ(IR))
      IR=0
      ITAG  = MTAGDISTSP+1+17
      DO JROC=1,NPROC
        DO JFLD=1,KFGATHG
          IF(KTO(JFLD) == MYPROC) THEN
            IRCV  = JROC
            IR=IR+1
            CALL MPL_RECV(ZBUF(1+IOFFR(IRCV):IOFFR(IRCV)+ILENR(IRCV)),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
             &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),&
             &CDSTRING='GATH_GRID_CTL:')
          ENDIF
        ENDDO
      ENDDO
      DO JFLD=1,KFGATHG
        ISND  = KTO(JFLD)
        IR=IR+1
        CALL MPL_SEND(ZFLD(1+IOFFS(ISND):IOFFS(ISND)+ILENS(ISND)),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),&
         &CDSTRING='GATH_GRID_CTL:')
      ENDDO
      CALL MPL_WAIT(KREQUEST=IREQ(1:IR), &
       & CDSTRING='GATH_GRID_CTL: WAIT')
      DEALLOCATE(IREQ)
    ENDIF
  ENDIF
  
  CALL GSTATS(809,1)
  CALL GSTATS_BARRIER2(789) 
  CALL GSTATS(1643,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)&
!$OMP&PRIVATE(JA,JB,IPROC,IGLOFF,IGL1,IGL2,IOFF,ILAST,J,&
!$OMP&ILEN,ILOFF,JGL,JLON,JFLD)
  DO JFLD=1,IFLDR
    DO JA=1,N_REGIONS_NS
      DO JB=1,N_REGIONS(JA)
        CALL SET2PE(IPROC,JA,JB,0,0)
        IGLOFF = D%NPTRFRSTLAT(JA)
        IGL1 = D%NFRSTLAT(JA)
        IGL2 = D%NLSTLAT(JA)
        IOFF = 0
        IF(JA > 1) THEN
          IF( D%NLSTLAT(JA-1) == D%NFRSTLAT(JA) )THEN
            ILAST = D%NLSTLAT(JA-1)-1
          ELSE
            ILAST = D%NLSTLAT(JA-1)
          ENDIF
          DO J=D%NFRSTLAT(1),ILAST
            IOFF = IOFF+G%NLOEN(J)
          ENDDO
        ENDIF

        ILEN = 0
        ILOFF = 0
        DO JGL=IGL1,IGL2
          DO JLON=1,D%NONL(IGLOFF+JGL-IGL1,JB)
            PGPG(IOFF+ILOFF+D%NSTA(IGLOFF+JGL-IGL1,JB)+JLON-1,JFLD) = &
             & ZBUF(ILEN+JLON+(JFLD-1)*IFLDL+(IPROC-1)*IFLDL*IMYFIELDS)
          ENDDO
          ILEN = ILEN + D%NONL(IGLOFF+JGL-IGL1,JB)
          ILOFF = ILOFF + G%NLOEN(JGL)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  CALL GSTATS(1643,1)
! Synhronize processors
! Should not be necessary
!!$  CALL GSTATS(784,0)
!!$  CALL MPL_BARRIER(CDSTRING='GATH_GRID_CTL:')
!!$  CALL GSTATS(784,1)
  IF(ALLOCATED(ZBUF)) DEALLOCATE(ZBUF)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID_CTL
END MODULE GATH_GRID_CTL_MOD


