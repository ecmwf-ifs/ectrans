MODULE DIST_GRID_CTL_MOD
CONTAINS
SUBROUTINE DIST_GRID_CTL(PGPG,KFDISTG,KPROMA,KFROM,PGP)

!**** *GATH_GRID_CTL* - Gather global gridpoint array from processors

!     Purpose.
!     --------
!        Routine for gathering gridpoint array

!**   Interface.
!     ----------
!     CALL GATH_GRID_CTL(...)

!     Explicit arguments : 
!**** *GATH_GRID_CTL* - Gather global gridpoint array from processors

!     Purpose.
!     --------
!        Routine for gathering gridpoint array

!**   Interface.
!     ----------
!     CALL GATH_GRID_CTL(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:) - Global gridpoint array
!     KFDISTG     - Global number of fields to be gathered
!     KTO(:)    - Processor responsible for gathering each field
!     PGP(:,:)  - Local spectral array
!

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DISTR
USE TPM_GEOMETRY

USE SET2PE_MOD

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KFROM(:)
REAL_B             , INTENT(OUT) :: PGP(:,:,:)

REAL_B :: ZFLD(D%NGPTOTMX)
INTEGER_M :: JFLD,JB,JA,IGLOFF,IGL1,IGL2,IOFF,ILAST,ILEN,ILOFF
INTEGER_M :: JGL,JLON,ISND,ITAG,IERR,ISENDER,ITAGR,J,IRCV
INTEGER_M :: JKGLO,IEND,JROF
!     ------------------------------------------------------------------

IF(NPROC == 1) THEN
  DO JFLD=1,KFDISTG
    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      DO JROF=1,IEND
        PGP(JROF,JFLD,JKGLO) = PGPG(IOFF+JROF,JFLD) 
      ENDDO
    ENDDO
  ENDDO
ELSE
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
          CALL MPE_SEND(ZFLD,ILEN,MREALT,NPRCIDS(ISND),ITAG,0,0,0,IERR)
          IF( IERR < 0 )THEN
            CALL ABOR1(' DIST_GRID_CTL:ERROR IN MPE_SEND')
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Recieve

    IRCV = KFROM(JFLD)
    ITAG = MTAGDISTGP
    CALL MPE_RECV(ZFLD,D%NGPTOT,MREALT,NPRCIDS(IRCV),ITAG,&
     &0,0,0,ILEN,ISENDER,ITAGR,IERR)
    IF( IERR < 0 ) CALL ABOR1(' DIST_GRID_CTL:ERROR IN MPE_RECV')
    IF( ILEN /= D%NGPTOT )THEN
      CALL ABOR1(' DIST_GRID_CTL: INVALID RECEIVE MESSAGE LENGTH')
    ENDIF

    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      DO JROF=1,IEND
        PGP(JROF,JFLD,JKGLO) = ZFLD(IOFF+JROF) 
      ENDDO
    ENDDO

    !Synchronize processors
    CALL MPE_BARRIER(IERR)
    IF( IERR /= 0 )THEN
      CALL ABOR1(' DIST_GRID_CTL: ERROR IN MPE_BARRIER')
    ENDIF

  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID_CTL
END MODULE DIST_GRID_CTL_MOD




