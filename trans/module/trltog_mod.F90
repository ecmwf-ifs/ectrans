MODULE TRLTOG_MOD
CONTAINS
SUBROUTINE TRLTOG(PGLAT,PGCOL,KVSET)

!**** *trltog * - transposition of grid point data from latitudinal
!                 to column structure. This takes place between inverse
!                 FFT and grid point calculations. 
!                 TRLTOG is the inverse of TRGTOL

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *call* *trltog(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGCOL    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original  : 95-10-01
!        D.Dent    : 97-08-04 Reorganisation to allow NPRTRV
!                             to differ from NPRGPEW
!        =99-03-29= Mats Hamrud and Deborah Salmond
!                   JUMP in FFT's changed to 1
!                   INDEX introduced and ZCOMBUF not used for same PE
!     ------------------------------------------------------------------



#include "tsmbkind.h"
USE MPL_MODULE

USE TPM_GEN
USE TPM_DISTR
USE TPM_TRANS

USE INIGPTR_MOD
USE PE2SET_MOD


IMPLICIT NONE



REAL_B,INTENT(IN)     :: PGLAT(:,:)
REAL_B,INTENT(OUT)    :: PGCOL(:,:,:)
INTEGER_M,INTENT(IN)  :: KVSET(:)

REAL_B,ALLOCATABLE :: ZCOMBUFS(:),ZCOMBUFR(:)

INTEGER_M :: ISENT    (NPROC)
INTEGER_M :: IRCVD    (NPROC)
INTEGER_M :: ISENDTOT (NPROC)
INTEGER_M :: IRECVTOT (NPROC)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IERR, IFIRST, IFIRSTLAT, IFLD, IGL, IGLL,&
             &ILAST, ILASTLAT, ILRECV, IPE, IPOS, ISETA, &
             &ISETB, IRCVTG, IRECV, IRECVD, IRECVSET, &
             &ISETV, ISEND, ISENDSET, ITAG,  JBLK, JFLD, &
             &JGL, JK, JL, JLOOP, ISETW, IFLDS, IPROC,JROC
INTEGER_M :: II,INDOFFX,ILEN,IBUFLENS,IBUFLENR,INRECV

!     LOCAL LOGICAL SCALARS
LOGICAL   :: LLDONE, LLEXIST, LLRECV
INTEGER_M :: INDEX(D%NLENGTF),INDOFF(NPROC),IFLDOFF(NF_FS)
INTEGER_M :: IRECV_FLD_START,IRECV_FLD_END
INTEGER_M :: ISEND_FLD_START(NPROC),ISEND_FLD_END
INTEGER_M :: ICOMBFLEN
INTEGER_M :: INUMFLDS
INTEGER_M :: IPROCS_COMM_MAX
INTEGER_M :: IGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER_M :: IGPTRRECV(NPRTRNS)
INTEGER_M :: IGPTROFF(NGPBLKS)


!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------

CALL INIGPTR(IGPTRSEND,IGPTRRECV)
LLDONE = .FALSE.
ITAG   = MTAGLG

INDOFFX  = 0
IBUFLENS = 0
IBUFLENR = 0
INRECV   = 0

DO JROC=1,NPROC

  CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)
  ISEND      = JROC
  ISENT(JROC) = 0
  IRCVD(JROC) = 0

!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,NF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1) IPOS = IPOS+1
  ENDDO
  IRECVTOT(JROC) = IGPTRRECV(ISETW)*IPOS

  IF(IRECVTOT(JROC) > 0 .AND. MYPROC /= JROC) INRECV = INRECV + 1

  IF( JROC /= MYPROC) IBUFLENR = MAX(IBUFLENR,IRECVTOT(JROC))

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(ISETA))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(ISETA))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO

  ISENDTOT(JROC) = IPOS*NF_FS
  IF( JROC /= MYPROC) IBUFLENS = MAX(IBUFLENS,ISENDTOT(JROC))

  IF(IPOS > 0) THEN
    INDOFF(JROC) = INDOFFX
    INDOFFX = INDOFFX+IPOS
    IPOS = 0
    DO JGL=IFIRSTLAT,ILASTLAT
      IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
      IGLL = JGL-D%NPTRLS(MYSETW)+1
      DO JL=D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL),&
       &D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL)+D%NONL(IGL,ISETB)-1
        IPOS = IPOS+1
        INDEX(IPOS+INDOFF(JROC)) = JL
      ENDDO
    ENDDO
  ENDIF
ENDDO

IPROCS_COMM_MAX = (NPRTRNS/NPRGPNS+1)*NPRTRV

ICOMBFLEN = NCOMBFLEN/IPROCS_COMM_MAX

IF (IBUFLENS > 0) ALLOCATE(ZCOMBUFS(-1:ICOMBFLEN))
IF (IBUFLENR > 0) ALLOCATE(ZCOMBUFR(-1:ICOMBFLEN))


! Send loop.............................................................

ISEND_FLD_START(:) = 1

DO WHILE( .NOT.LLDONE ) 
LLDONE = .TRUE.  

! 
! loop over the number of processors we need to communicate with. 
! NOT MYPROC 
! 

IPROC = MYPROC 
DO JLOOP=1,NPROC-1 
  IPROC = IPROC+1 
  IF(IPROC > NPROC) THEN
    IPROC = IPROC-NPROC
  ENDIF

! Check if we have data to send

  IF(ISENDTOT(IPROC)-ISENT(IPROC) > 0 )THEN

    CALL PE2SET(IPROC,ISETA,ISETB,ISETW,ISETV)
    ISEND = IPROC
    ISENDSET = ISETA

    ILEN = ISENDTOT(ISEND)/NF_FS

    INUMFLDS = ICOMBFLEN/ILEN

    ISEND_FLD_END = MIN(NF_FS,ISEND_FLD_START(IPROC)+INUMFLDS-1)

    DO JL=1,ILEN
      II = INDEX(INDOFF(ISEND)+JL)
      DO JFLD=ISEND_FLD_START(IPROC),ISEND_FLD_END
        ZCOMBUFS(JL+(JFLD-ISEND_FLD_START(IPROC))*ILEN) = PGLAT(JFLD,II)
      ENDDO
    ENDDO

    ZCOMBUFS(-1) = ISEND_FLD_START(IPROC)
    ZCOMBUFS(0)  = ISEND_FLD_END

    IPOS = ILEN*(ISEND_FLD_END-ISEND_FLD_START(IPROC)+1)

    ISEND_FLD_START(IPROC) = ISEND_FLD_END+1

    LLDONE = .FALSE.
    CALL MPL_SEND(ZCOMBUFS(-1:IPOS),KDEST=NPRCIDS(ISEND),&
     & KTAG=ITAG,CDSTRING='TRLTOG:')
    ISENT(IPROC) = ISENT(IPROC)+IPOS
  ENDIF
ENDDO


!  Receive loop.........................................................


! Check if there is data to receive

  LLRECV = .FALSE.
  DO JROC=1,NPROC
    IF( IRECVTOT(JROC)-IRCVD(JROC) > 0 )THEN
      LLRECV = .TRUE.
    ENDIF
  ENDDO

  IF( LLRECV )THEN

! Probe and wait for a message from anybody

    IF(INRECV > 1) THEN
      CALL MPL_PROBE(KTAG=ITAG,LDWAIT=.TRUE.,LDFLAG=LLEXIST, &
       &    CDSTRING='TRLTOG: PROBEING FOR ANY MESSAGE ' )
    ENDIF

  ENDIF

! A message should now be ready to be received
!
! loop over the number of processors we need to communicate with
! For MYPROC copy straight from PGLAT to PGCOL
!
  IPROC = MYPROC-1
  DO JLOOP=1,NPROC
    IPROC = IPROC+1
    IF(IPROC > NPROC) THEN
      IPROC = IPROC-NPROC
    ENDIF

! Check if there is data to receive

    IF( IRECVTOT(IPROC)-IRCVD(IPROC) > 0 )THEN
      CALL PE2SET(IPROC,ISETA,ISETB,ISETW,ISETV)
      IRECVSET = ISETV

! There is data to receive, probe for a message

      IRECV = IPROC

      IF(INRECV > 1 .AND. MYPROC /= IRECV) THEN
        CALL MPL_PROBE(KSOURCE=NPRCIDS(IRECV),KTAG=ITAG,LDWAIT=.FALSE., &
         &   LDFLAG=LLEXIST,CDSTRING='TRLTOG: PROBEING FOR MESSAGE ' )
      ELSE
        LLEXIST = .TRUE.
      ENDIF

! If a message exists, receive it, otherwise flag an outstanding receive

      IF( LLEXIST )THEN

!*   copy data from buffer to column structure

        IF (IRECV /= MYPROC) THEN

!*   receive message

          LLDONE = .FALSE.

          CALL MPL_RECV(ZCOMBUFR(-1:ICOMBFLEN),KSOURCE=NPRCIDS(IRECV), &
           & KTAG=ITAG,KOUNT=ILRECV,CDSTRING='TRLTOG:' )
  
          IRECV_FLD_START = ZCOMBUFR(-1)
          IRECV_FLD_END   = ZCOMBUFR(0)
          IFLD = 0
          IPOS = 0
          DO JFLD=1,NF_GP
            IF(KVSET(JFLD) == IRECVSET .OR. KVSET(JFLD) == -1 ) THEN
              IF (IFLD == IRECV_FLD_END) EXIT
              IFLD = IFLD+1
              IF (IFLD >= IRECV_FLD_START) THEN
                DO JBLK=1,NGPBLKS
                  IFIRST = IGPTRSEND(1,JBLK,ISETW)
                  IF(IFIRST > 0) THEN
                    ILAST = IGPTRSEND(2,JBLK,ISETW)
                    DO JK=IFIRST,ILAST
                      IPOS = IPOS+1
                      PGCOL(JK,JFLD,JBLK) = ZCOMBUFR(IPOS)
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO
    
          INRECV = INRECV -1
    
          IF( ILRECV /= IPOS+2 )THEN
            WRITE(NOUT,*) MYPROC,' receiving ',ILRECV,'expected ',IPOS
            WRITE(NOUT,*) 'IRECV=',IRECV,'IRECVSET=',IRECVSET
            WRITE(NOUT,*) 'NF_GP=',NF_GP,'KVSET=',KVSET(1:NF_GP)
            CALL ABOR1('TRLTOG: RECEIVED MESSAGE OF INCORRECT LENGTH')
          ENDIF
  
          IRCVD(IPROC) = IRCVD(IPROC)+IPOS
  
        ELSE

          IFLDS = 0
          DO JFLD=1,NF_GP
            IF(KVSET(JFLD) == IRECVSET .OR. KVSET(JFLD) == -1) THEN
              IFLDS = IFLDS+1
              IFLDOFF(IFLDS) = JFLD
            ENDIF
          ENDDO

          IPOS=0
          DO JBLK=1,NGPBLKS
            IGPTROFF(JBLK)=IPOS
            IFIRST = IGPTRSEND(1,JBLK,ISETW)
            IF(IFIRST > 0) THEN
              ILAST = IGPTRSEND(2,JBLK,ISETW)
              IPOS=IPOS+ILAST-IFIRST+1
            ENDIF
          ENDDO
!$OMP PARALLEL DO PRIVATE(JFLD,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
          DO JBLK=1,NGPBLKS
            IFIRST = IGPTRSEND(1,JBLK,ISETW)
            IF(IFIRST > 0) THEN
              ILAST = IGPTRSEND(2,JBLK,ISETW)
              DO JFLD=1,IFLDS
                IFLD = IFLDOFF(JFLD)
                DO JK=IFIRST,ILAST
                  IPOS = INDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
                  PGCOL(JK,IFLD,JBLK) = PGLAT(JFLD,INDEX(IPOS))
                ENDDO
              ENDDO
            ENDIF
          ENDDO
!$OMP END PARALLEL DO


          IRCVD(IPROC) = IRECVTOT(IPROC)

        ENDIF

      ELSE

! Expecting data but none found, set outstanding receive flag

        LLDONE = .FALSE.

      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (IBUFLENS > 0) DEALLOCATE(ZCOMBUFS)
IF (IBUFLENR > 0) DEALLOCATE(ZCOMBUFR)

! Perform barrier synchronisation to guarantee all processors have
! completed communication

IF( NPROC > 1 .AND. (D%LSPLIT .OR. NPRTRV > 1 .OR. NPRGPEW > 1 ))THEN
  CALL MPL_BARRIER(CDSTRING='TRLTOG:')
ENDIF


RETURN
END SUBROUTINE TRLTOG
END MODULE TRLTOG_MOD






