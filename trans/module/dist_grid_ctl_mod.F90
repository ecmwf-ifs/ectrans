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
!
!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

#include "tsmbkind.h"
USE MPL_MODULE
USE YOMGSTATS, ONLY : LSYNCSTATS

USE TPM_DISTR
USE TPM_GEOMETRY

USE SET2PE_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KFROM(:)
REAL_B             , INTENT(OUT) :: PGP(:,:,:)

! Declaration of local variables

REAL_B    :: ZFLD(D%NGPTOTMX)
INTEGER_M :: JFLD,JB,JA,IGLOFF,IGL1,IGL2,IOFF,ILAST,ILEN,ILOFF,ILENR
INTEGER_M :: JGL,JLON,ISND,ITAG,IERR,ISENDER,ITAGR,J,IRCV
INTEGER_M :: JKGLO,IEND,JROF,IBL

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

  DO JFLD=1,KFDISTG

    ! Send

    IF(KFROM(JFLD) == MYPROC) THEN
      DO JB=1,NPRGPEW
        DO JA=1,NPRGPNS

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

          ILEN = 0
          ILOFF = 0
          DO JGL=IGL1,IGL2
            DO JLON=1,D%NONL(IGLOFF+JGL-IGL1,JB)
              ZFLD(ILEN+JLON) = &
               & PGPG(IOFF+ILOFF+D%NSTA(IGLOFF+JGL-IGL1,JB)+JLON-1,JFLD)
            ENDDO
            ILEN = ILEN + D%NONL(IGLOFF+JGL-IGL1,JB)
            ILOFF = ILOFF + G%NLOEN(JGL)
          ENDDO

          CALL SET2PE(ISND,JA,JB,0,0)
          ITAG = MTAGDISTGP
          CALL MPL_SEND(ZFLD(1:ILEN),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
           &CDSTRING='DIST_GRID_CTL')
        ENDDO
      ENDDO
    ENDIF

    ! Recieve

    IRCV = KFROM(JFLD)
    ITAG = MTAGDISTGP
    CALL MPL_RECV(ZFLD(1:D%NGPTOT),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
     &KOUNT=ILENR,CDSTRING='DIST_GRID_CTL:')
    IF( ILENR /= D%NGPTOT )THEN
      CALL ABORT_TRANS(' DIST_GRID_CTL: INVALID RECEIVE MESSAGE LENGTH')
    ENDIF

    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      IBL  = (JKGLO-1)/KPROMA+1
      DO JROF=1,IEND
        PGP(JROF,JFLD,IBL) = ZFLD(IOFF+JROF) 
      ENDDO
    ENDDO

    !Synchronize processors
    CALL MPL_BARRIER(CDSTRING='DIST_GRID_CTL:')

  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID_CTL
END MODULE DIST_GRID_CTL_MOD




