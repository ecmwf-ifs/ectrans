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

USE TPM_GEN
USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_DISTR

USE SET2PE_MOD

IMPLICIT NONE

REAL_B    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER_M          , INTENT(IN)  :: KFGATHG
INTEGER_M          , INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KTO(:)
REAL_B             , INTENT(IN)  :: PGP(:,:,:)

REAL_B    :: ZFLD(D%NGPTOTMX)
INTEGER_M :: IFLDR,JFLD,ITAG,ILEN,JA,JB,IERR,ISND,JGL,JLON,ILOFF,ILENB
INTEGER_M :: ISENDER,ITAGR,IRCV,ILENR,IOFF,ILAST,IGL1,IGL2,IGLOFF
INTEGER_M :: JKGLO,JROF,IEND,J
!     ------------------------------------------------------------------


!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
  DO JFLD=1,KFGATHG
    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      DO JROF=1,IEND
        PGPG(IOFF+JROF,JFLD) = PGP(JROF,JFLD,JKGLO) 
      ENDDO
    ENDDO
  ENDDO
ELSE
  IFLDR = 0

  DO JFLD=1,KFGATHG

  !Send

    ISND  = KTO(JFLD)
    ITAG  = MTAGDISTSP+JFLD+17
    DO JKGLO=1,D%NGPTOT,KPROMA
      IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
      IOFF = JKGLO-1
      DO JROF=1,IEND
        ZFLD(IOFF+JROF) = PGP(JROF,JFLD,JKGLO) 
      ENDDO
    ENDDO
    CALL MPE_SEND(ZFLD,D%NGPTOT,MREALT,NPRCIDS(ISND),ITAG,0,&
     &0,0,IERR)
    IF( IERR < 0 )THEN
      CALL ABOR1(' GATH_GRID_CTL : ERROR IN MPE_SEND (ZFLD)')
    ENDIF
     
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
          CALL MPE_RECV(ZFLD,ILENB,MREALT,NPRCIDS(IRCV),&
           &ITAG,0,0,0,ILENR,ISENDER,ITAGR,IERR)
          IF( IERR < 0 )THEN
            CALL ABOR1(' GATH_GRID_CTL : GG ERROR IN MPE_RECV')
          ENDIF

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
    CALL MPE_BARRIER(IERR)
    IF( IERR /= 0 )THEN
      CALL ABOR1(' GATH_GRID_CTL: ERROR IN MPE_BARRIER')
    ENDIF

  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID_CTL
END MODULE GATH_GRID_CTL_MOD


