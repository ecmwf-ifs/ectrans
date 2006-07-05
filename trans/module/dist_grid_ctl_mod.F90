MODULE DIST_GRID_CTL_MOD
CONTAINS
SUBROUTINE DIST_GRID_CTL(PGPG,KFDISTG,KPROMA,KFROM,PGP)

!**** *DIST_GRID_CTL* - Distributing global gridpoint array to processors

!     Purpose.
!     --------
!        Routine for distributing gridpoint array

!**   Interface.
!     ----------
!     CALL DIST_GRID_CTL(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:)   - Global gridpoint array
!     KFDISTG     - Global number of fields to be distributed
!     KPROMA      - required blocking factor for gridpoint output
!     KFROM(:)    - Processor responsible for distributing each field
!     PGP(:,:,:)  - Local spectral array

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

USE TPM_DISTR
USE TPM_GEOMETRY

USE SET2PE_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
REAL(KIND=JPRB)             , INTENT(OUT) :: PGP(:,:,:)

! Declaration of local variables

REAL(KIND=JPRB)    :: ZFLD(D%NGPTOTMX,NPROC),ZDUM(D%NGPTOTMX)
REAL(KIND=JPRB)    :: ZRCV(D%NGPTOT)
INTEGER(KIND=JPIM) :: JFLD,JB,JA,IGLOFF,IGL1,IGL2,IOFF,ILAST,ILOFF,ILENR
INTEGER(KIND=JPIM) :: JGL,JLON,ISND,ITAG,J,IRCV,IFLDS
INTEGER(KIND=JPIM) :: JKGLO,IEND,JROF,IBL,JROC
INTEGER(KIND=JPIM) :: ISENDREQ(NPROC),ILEN(NPROC)

!     ------------------------------------------------------------------

! Copy for single PE

IF(NPROC == 1) THEN
  DO JFLD=1,KFDISTG
    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      IBL  = (JKGLO-1)/KPROMA+1
      DO JROF=1,IEND
        PGP(JROF,JFLD,IBL) = PGPG(IOFF+JROF,JFLD) 
      ENDDO
    ENDDO
  ENDDO
ELSE

  ! Message passing

  IFLDS = 0
  DO JFLD=1,KFDISTG

    ! Send

    IF(KFROM(JFLD) == MYPROC) THEN
      IFLDS = IFLDS+1
      CALL GSTATS(1663,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)&
!$OMP&PRIVATE(JB,JA,ISND,IGLOFF,IGL1,IGL2,IOFF,ILAST,J,&
!$OMP&ILOFF,JGL,JLON)
      DO JB=1,NPRGPEW
        DO JA=1,NPRGPNS
          CALL SET2PE(ISND,JA,JB,0,0)

          IGLOFF = D%NPTRFRSTLAT(JA)
          IGL1   = D%NFRSTLAT(JA)
          IGL2   = D%NLSTLAT(JA)
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

          ILEN(ISND) = 0
          ILOFF = 0
          DO JGL=IGL1,IGL2
            DO JLON=1,D%NONL(IGLOFF+JGL-IGL1,JB)
              ZFLD(ILEN(ISND)+JLON,ISND) = &
               & PGPG(IOFF+ILOFF+D%NSTA(IGLOFF+JGL-IGL1,JB)+JLON-1,IFLDS)
            ENDDO
            ILEN(ISND) = ILEN(ISND) + D%NONL(IGLOFF+JGL-IGL1,JB)
            ILOFF = ILOFF + G%NLOEN(JGL)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1663,1)
      CALL GSTATS(811,0)

      ITAG = MTAGDISTGP+JFLD
      DO JROC=1,NPROC
        CALL MPL_SEND(ZFLD(1:ILEN(JROC),JROC),KDEST=NPRCIDS(JROC),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JROC),&
         &CDSTRING='DIST_GRID_CTL')
      ENDDO
      CALL GSTATS(811,1)
    ENDIF

    ! Receive

    CALL GSTATS(811,0)
    IRCV = KFROM(JFLD)
    ITAG = MTAGDISTGP+JFLD
    CALL MPL_RECV(ZRCV,KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
     &KMP_TYPE=JP_BLOCKING_STANDARD,KOUNT=ILENR,CDSTRING='DIST_GRID_CTL:')
    IF( ILENR /= D%NGPTOT )THEN
      CALL ABORT_TRANS(' DIST_GRID_CTL: INVALID RECEIVE MESSAGE LENGTH')
    ENDIF

    IF(KFROM(JFLD) == MYPROC) THEN
      CALL MPL_WAIT(ZDUM,KREQUEST=ISENDREQ(:), &
       & CDSTRING='DIST_GRID_CTL: WAIT')
    ENDIF
    CALL GSTATS(811,1)

    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      IBL  = (JKGLO-1)/KPROMA+1
      DO JROF=1,IEND
        PGP(JROF,JFLD,IBL) = ZRCV(IOFF+JROF) 
      ENDDO
    ENDDO
  ENDDO
  !Synchronize processors
  CALL GSTATS(786,0)
  CALL MPL_BARRIER(CDSTRING='DIST_GRID_CTL:')
  CALL GSTATS(786,1)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID_CTL
END MODULE DIST_GRID_CTL_MOD




