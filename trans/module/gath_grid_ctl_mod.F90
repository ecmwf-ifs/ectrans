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
!     PGP(:,:,:)  - Local spectral array
!
!     ------------------------------------------------------------------


#include "tsmbkind.h"
USE MPL_MODULE
USE YOMGSTATS, ONLY : LSYNCSTATS

USE TPM_GEN
USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_DISTR

USE SET2PE_MOD

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER_M          , INTENT(IN)  :: KFGATHG
INTEGER_M          , INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KTO(:)
REAL_B             , INTENT(IN)  :: PGP(:,:,:)

! Declaration of local variables

REAL_B    :: ZFLD(D%NGPTOTMX)
INTEGER_M :: IFLDR,JFLD,ITAG,ILEN,JA,JB,IERR,ISND,JGL,JLON,ILOFF,ILENB
INTEGER_M :: IRCV,IOFF,ILAST,IGL1,IGL2,IGLOFF
INTEGER_M :: JKGLO,JROF,IEND,J,IBL
!     ------------------------------------------------------------------


!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
IF (.NOT.LSYNCSTATS) CALL GSTATS(1643,0)
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
IF (.NOT.LSYNCSTATS) CALL GSTATS(1643,1)

ELSE
  IFLDR = 0

  DO JFLD=1,KFGATHG

  !Send

    ISND  = KTO(JFLD)
    ITAG  = MTAGDISTSP+JFLD+17
    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      IBL  = (JKGLO-1)/KPROMA+1
      DO JROF=1,IEND
        ZFLD(IOFF+JROF) = PGP(JROF,JFLD,IBL) 
      ENDDO
    ENDDO
    CALL MPL_SEND(ZFLD(1:D%NGPTOT),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
     &CDSTRING='GATH_GRID_CTL:')
     
  ! RECIEVE
    IF(KTO(JFLD) == MYPROC) THEN
      IFLDR = IFLDR+1
      DO JB=1,NPRGPEW
        DO JA=1,NPRGPNS
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
          CALL SET2PE(IRCV,JA,JB,0,0)
          ITAG  = MTAGDISTSP+JFLD+17
          ILENB = D%NGPTOTMX
          CALL MPL_RECV(ZFLD(1:ILENB),KSOURCE=NPRCIDS(IRCV),&
           &KTAG=ITAG,CDSTRING='GATH_GRID_CTL:')

          ILEN = 0
          ILOFF = 0
          DO JGL=IGL1,IGL2
            DO JLON=1,D%NONL(IGLOFF+JGL-IGL1,JB)
              PGPG(IOFF+ILOFF+D%NSTA(IGLOFF+JGL-IGL1,JB)+JLON-1,IFLDR) = &
               & ZFLD(ILEN+JLON)
            ENDDO
            ILEN = ILEN + D%NONL(IGLOFF+JGL-IGL1,JB)
            ILOFF = ILOFF + G%NLOEN(JGL)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  !Synchronize processors
    CALL MPL_BARRIER(CDSTRING='GATH_GRID_CTL:')

  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID_CTL
END MODULE GATH_GRID_CTL_MOD


