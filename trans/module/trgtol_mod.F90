MODULE TRGTOL_MOD
CONTAINS
SUBROUTINE TRGTOL(PGLAT,PGCOL,KVSET)

!**** *TRGTOL * - transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *call* *trgtol(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (output)
!           PGCOL    -  Blocked grid point data    (input)

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
!        Original: 95-10-01 
!        D.Dent  : 97-08-04   Reorganisation to allow
!                             NPRTRV to differ from NPRGPEW 
!                : 98-06-17   add mailbox control logic (from TRLTOM)
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

REAL_B,INTENT(OUT)   :: PGLAT(:,:)
REAL_B,INTENT(IN)    :: PGCOL(:,:,:)
INTEGER_M,INTENT(IN) :: KVSET(:)

REAL_B,ALLOCATABLE :: ZCOMBUFS(:),ZCOMBUFR(:)

INTEGER_M :: ISENT    (NPROC)
INTEGER_M :: IRCVD    (NPROC)
INTEGER_M :: ISENDTOT (NPROC)
INTEGER_M :: IRECVTOT (NPROC)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IERR, IFIRST, IFIRSTLAT, IGL, IGLL, ILAST,&
             &ILASTLAT, ILEN, ILRECV, JROC, IPOS, ISETA, &
             &ISETB, IRCVTG, IRECV, IRECVD, IRECVSET, &
             &ISETV, ISEND, ISENDSET, ITAG, J, JBLK, JFLD, &
             &JGL, JK, JL, JLOOP, ISETW,  JVAR,IFLD, &
             &II,INDOFFX,IBUFLENS,IBUFLENR,INRECV, IPROC,IFLDS


!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE, LLEXIST,  LLRECV
INTEGER_M :: INDEX(D%NLENGTF),INDOFF(NPROC),IFLDOFF(NF_FS)
INTEGER_M :: ISEND_FLD_TOTAL(NPROC),ISEND_FLD_START(NPROC),ISEND_FLD_END
INTEGER_M :: IRECV_FLD_START,IRECV_FLD_END
INTEGER_M :: ICOMBFLEN
INTEGER_M :: INUMFLDS
INTEGER_M :: IPROCS_COMM_MAX
INTEGER_M :: IGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER_M :: IGPTRRECV(NPRTRNS)

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------

CALL INIGPTR(IGPTRSEND,IGPTRRECV)
LLDONE = .FALSE.

ITAG = MTAGGL

INDOFFX  = 0
IBUFLENS = 0
IBUFLENR = 0
INRECV   = 0

DO JROC=1,NPROC

  CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)
  IRECVSET = ISETA
  ISEND = JROC
  ISENDSET = ISETV
  ISENT(JROC) = 0
  IRCVD(JROC) = 0

!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,NF_GP
    IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1) IPOS = IPOS+1
  ENDDO
  ISEND_FLD_TOTAL(JROC) = IPOS
  ISENDTOT(JROC) = IGPTRRECV(ISETW)*IPOS

  IF( JROC /= MYPROC) IBUFLENS = MAX(IBUFLENS,ISENDTOT(JROC))

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(IRECVSET))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(IRECVSET))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(IRECVSET)+JGL-D%NFRSTLAT(IRECVSET)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO
  
  IRECVTOT(JROC) = IPOS*NF_FS

  IF(IRECVTOT(JROC) > 0 .AND. MYPROC /= JROC) INRECV = INRECV + 1

  IF( JROC /= MYPROC) IBUFLENR = MAX(IBUFLENR,IRECVTOT(JROC))

  IF(IPOS > 0) THEN
    INDOFF(JROC) = INDOFFX
    INDOFFX = INDOFFX+IPOS
    IPOS = 0
    DO JGL=IFIRSTLAT,ILASTLAT
      IGL  = D%NPTRFRSTLAT(IRECVSET)+JGL-D%NFRSTLAT(IRECVSET)
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

IF (IBUFLENR > 0) ALLOCATE(ZCOMBUFR(-1:ICOMBFLEN))
IF (IBUFLENS > 0) ALLOCATE(ZCOMBUFS(-1:ICOMBFLEN))


! Send loop.............................................................

ISEND_FLD_START(:) = 1

DO WHILE( .NOT.LLDONE )

  LLDONE = .TRUE.
!
! loop over the number of processors we need to communicate with
! For MYPROC copy straight from PGCOL to PGLAT
!

  IPROC = MYPROC
  DO JLOOP=1,NPROC
    IPROC = IPROC+1
    IF(IPROC > NPROC) THEN
      IPROC = IPROC-NPROC
    ENDIF

! Check if we have data to send

    IF(ISENDTOT(IPROC)-ISENT(IPROC) > 0 )THEN

      CALL PE2SET(IPROC,ISETA,ISETB,ISETW,ISETV)
      ISEND = IPROC
      ISENDSET = ISETV


      IF (ISEND /= MYPROC) THEN

        INUMFLDS = ICOMBFLEN/IGPTRRECV(ISETW)

        ISEND_FLD_END = MIN(ISEND_FLD_TOTAL(IPROC),&
         &ISEND_FLD_START(IPROC)+INUMFLDS-1)

        ZCOMBUFS(-1) = ISEND_FLD_START(IPROC)
        IFLD = 0
        IPOS = 0
        DO JFLD=1,NF_GP
          IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1 ) THEN

            IF (IFLD == ISEND_FLD_END) EXIT
            IFLD = IFLD+1
            IF (IFLD >= ISEND_FLD_START(IPROC)) THEN
            
              DO JBLK=1,NGPBLKS
                IFIRST = IGPTRSEND(1,JBLK,ISETW)
                IF(IFIRST > 0) THEN
                  ILAST = IGPTRSEND(2,JBLK,ISETW)
                  DO JK=IFIRST,ILAST
                    IPOS = IPOS+1
                    ZCOMBUFS(IPOS) = PGCOL(JK,JFLD,JBLK)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO

        ZCOMBUFS(0) = IFLD

        ISEND_FLD_START(IPROC) = IFLD+1

        LLDONE = .FALSE.
        CALL MPL_SEND(ZCOMBUFS(-1:IPOS),KDEST=NPRCIDS(ISEND), &
         &    KTAG=ITAG,CDSTRING='TRGTOL:' )
        ISENT(IPROC) = ISENT(IPROC)+IPOS

      ELSE

        IFLDS = 0
        DO JFLD=1,NF_GP
          IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1) THEN
            IFLDS = IFLDS+1
            IFLDOFF(IFLDS) = JFLD
          ENDIF
        ENDDO

!$OMP PARALLEL DO PRIVATE(JFLD,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
        DO JFLD=1,IFLDS
          IFLD = IFLDOFF(JFLD)
          IPOS = INDOFF(MYPROC)
          DO JBLK=1,NGPBLKS
            IFIRST = IGPTRSEND(1,JBLK,ISETW)
            IF(IFIRST > 0) THEN
              ILAST = IGPTRSEND(2,JBLK,ISETW)
              DO JK=IFIRST,ILAST
                IPOS = IPOS+1
                PGLAT(JFLD,INDEX(IPOS)) = PGCOL(JK,IFLD,JBLK)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        ISENT(IPROC) = ISENDTOT(IPROC)

      ENDIF
    ENDIF
  ENDDO

!  Receive loop.........................................................



! Check if there is data to receive

  LLRECV = .FALSE.
  DO J=1,NPROC
    IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
      LLRECV = .TRUE.
    ENDIF
  ENDDO

  IF( LLRECV )THEN

! Probe and wait for a message from anybody

    IF(INRECV > 1) THEN
      CALL MPL_PROBE(KTAG=ITAG,LDWAIT=.TRUE.,LDFLAG=LLEXIST, &
       &    CDSTRING='TRGTOL: PROBEING FOR ANY MESSAGE')
    ENDIF

  ENDIF

! A message should now be ready to be received

  IPROC = MYPROC
  DO JLOOP=1,NPROC-1
    IPROC = IPROC+1
    IF(IPROC > NPROC) THEN
      IPROC = IPROC-NPROC
    ENDIF
    IF( IRECVTOT(IPROC)-IRCVD(IPROC) > 0 )THEN
      CALL PE2SET(IPROC,ISETA,ISETB,ISETW,ISETV)
      IRECVSET = ISETA

! There is data to receive, probe for a message

      IRECV = IPROC

      IF(INRECV > 1) THEN
        CALL MPL_PROBE(KSOURCE=NPRCIDS(IRECV),KTAG=ITAG,LDWAIT=.FALSE., &
         &   LDFLAG=LLEXIST,CDSTRING='TRGTOL: PROBEING FOR MESSAGE' )
      ELSE
        LLEXIST = .TRUE.
      ENDIF

! If a message exists, receive it, otherwise flag an outstanding receive

      IF( LLEXIST )THEN

!*   receive message

        LLDONE = .FALSE.

        CALL MPL_RECV(ZCOMBUFR(-1:ICOMBFLEN),KSOURCE=NPRCIDS(IRECV), &
         & KTAG=ITAG,KOUNT=ILRECV,CDSTRING='TRGTOL:' )

        IRECV_FLD_START = ZCOMBUFR(-1)
        IRECV_FLD_END   = ZCOMBUFR(0)


!*   store data in Fourier array

        ILEN = IRECVTOT(IRECV)/NF_FS
        DO JL=1,ILEN
          II = INDEX(INDOFF(IRECV)+JL)
          DO JFLD=IRECV_FLD_START,IRECV_FLD_END
            PGLAT(JFLD,II) = ZCOMBUFR(JL+(JFLD-IRECV_FLD_START)*ILEN)
          ENDDO
        ENDDO
        IPOS = ILEN*(IRECV_FLD_END-IRECV_FLD_START+1)
        IF (ILRECV /= IPOS+2) THEN
          WRITE(NOUT,*) MYPROC,MYSETV,' receiving ',ILRECV,'expected',IPOS
          CALL ABOR1('TRGTOL: RECEIVED MESSAGE OF INCORRECT LENGTH')
        ENDIF
        IRCVD(IPROC) = IRCVD(IPROC)+IPOS
        INRECV = INRECV-1

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

IF( NPROC > 1.AND.(D%LSPLIT.OR.NPRTRV > 1))THEN
  CALL MPL_BARRIER(CDSTRING='TRGTOL:')
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE TRGTOL
END MODULE TRGTOL_MOD
