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
             &ILASTLAT, ILEN, ILRECV, IPE, IPOS, IPROCA, &
             &IPROCB, IRCVTG, IRECV, IRECVD, IRECVSET, &
             &IROCV, ISEND, ISENDSET, ITAG, J, JBLK, JFLD, &
             &JGL, JK, JL, JLOOP, JPRTRNS, JROCV, JVAR,IFLD, &
             &II,INDOFFX,IBUFLENS,IBUFLENR,INRECV


!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE, LLEXIST, LLWAIT, LLRECV
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
LLDONE=.FALSE.

ITAG=MTAGGL

INDOFFX = 0
IBUFLENS = 0
IBUFLENR = 0
INRECV = 0

DO IPE=1,NPROC

  CALL PE2SET(IPE,IPROCA,IPROCB,JPRTRNS,JROCV)
  IRECVSET = IPROCA
  ISEND = IPE
  ISENDSET = JROCV
  ISENT(IPE) = 0
  IRCVD(IPE) = 0

!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,NF_GP
    IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1) IPOS=IPOS+1
  ENDDO
  ISEND_FLD_TOTAL(IPE) = IPOS
  ISENDTOT(IPE) = IGPTRRECV(JPRTRNS)*IPOS

  IF( IPE /= MYPROC) IBUFLENS = MAX(IBUFLENS,ISENDTOT(IPE))

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(IRECVSET))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(IRECVSET))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(IRECVSET)+JGL-D%NFRSTLAT(IRECVSET)
    IPOS = IPOS+D%NONL(IGL,IPROCB)
  ENDDO
  
  IRECVTOT(IPE) = IPOS*NF_FS

  IF(IRECVTOT(IPE) > 0 .AND. MYPROC /= IPE) INRECV = INRECV + 1

  IF( IPE /= MYPROC) IBUFLENR=MAX(IBUFLENR,IRECVTOT(IPE))

  IF(IPOS > 0) THEN
    INDOFF(IPE) = INDOFFX
    INDOFFX = INDOFFX+IPOS
    IPOS = 0
    DO JGL=IFIRSTLAT,ILASTLAT
      IGL  = D%NPTRFRSTLAT(IRECVSET)+JGL-D%NFRSTLAT(IRECVSET)
      IGLL = JGL-D%NPTRLS(MYSETW)+1
      DO JL=D%NSTA(IGL,IPROCB)+D%NSTAGTF(IGLL),&
       &D%NSTA(IGL,IPROCB)+D%NSTAGTF(IGLL)+D%NONL(IGL,IPROCB)-1
        IPOS = IPOS+1
        INDEX(IPOS+INDOFF(IPE)) = JL
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

  J = MYPROC
  DO JLOOP=1,NPROC
    J = J+1
    IF(J > NPROC) THEN
      J = J-NPROC
    ENDIF

! Check if we have data to send

    IF(ISENDTOT(J)-ISENT(J) > 0 )THEN

      CALL PE2SET(J,IPROCA,IPROCB,JPRTRNS,JROCV)
      ISEND = J
      ISENDSET = JROCV


      IF (ISEND /= MYPROC) THEN

        INUMFLDS = ICOMBFLEN/IGPTRRECV(JPRTRNS)

        ISEND_FLD_END = MIN(ISEND_FLD_TOTAL(J),ISEND_FLD_START(J)+INUMFLDS-1)

        ZCOMBUFS(-1) = ISEND_FLD_START(J)
        IFLD = 0
        IPOS = 0
        DO JFLD=1,NF_GP
          IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1 ) THEN

            IF (IFLD == ISEND_FLD_END) EXIT
            IFLD = IFLD+1
            IF (IFLD >= ISEND_FLD_START(J)) THEN
            
              DO JBLK=1,NGPBLKS
                IFIRST = IGPTRSEND(1,JBLK,JPRTRNS)
                IF(IFIRST > 0) THEN
                  ILAST = IGPTRSEND(2,JBLK,JPRTRNS)
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

        ISEND_FLD_START(J) = IFLD+1

        LLDONE = .FALSE.
        CALL MPE_SEND(ZCOMBUFS(-1),IPOS+2,MREALT,NPRCIDS(ISEND),&
                     &ITAG,0,0,0,IERR)
        IF( IERR < 0 )THEN
          CALL ABOR1('TRGTOL : ERROR IN MPE_SEND')
        ENDIF
        ISENT(J) = ISENT(J)+IPOS

      ELSE

        IFLD = 0
        DO JFLD=1,NF_GP
          IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1) THEN
            IFLD = IFLD+1
            IFLDOFF(IFLD) = JFLD
          ENDIF
        ENDDO

!$OMP PARALLEL DO PRIVATE(JFLD,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
        DO JFLD=1,IFLD
          IFLD = IFLDOFF(JFLD)
          IPOS = INDOFF(MYPROC)
          DO JBLK=1,NGPBLKS
            IFIRST = IGPTRSEND(1,JBLK,JPRTRNS)
            IF(IFIRST > 0) THEN
              ILAST = IGPTRSEND(2,JBLK,JPRTRNS)
              DO JK=IFIRST,ILAST
                IPOS = IPOS+1
                PGLAT(JFLD,INDEX(IPOS)) = PGCOL(JK,IFLD,JBLK)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        ISENT(J) = ISENDTOT(J)

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

    IRECV = -1
    LLWAIT = .TRUE.
    IF(INRECV > 1) THEN
      CALL MPE_PROBE(IRECV,ITAG,0,LLWAIT,LLEXIST,IERR)
      IF( IERR < 0 )THEN
        CALL ABOR1('TRGTOL: ERROR IN MPE_PROBE WAIT=T')
      ENDIF
    ENDIF

  ENDIF

! A message should now be ready to be received

  J = MYPROC
  DO JLOOP=1,NPROC-1
    J = J+1
    IF(J > NPROC) THEN
      J = J-NPROC
    ENDIF
    IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
      CALL PE2SET(J,IPROCA,IPROCB,JPRTRNS,IROCV)
      IRECVSET = IPROCA

! There is data to receive, probe for a message

      IRECV = J
      LLWAIT = .FALSE.

      IF(INRECV > 1) THEN
        CALL MPE_PROBE(NPRCIDS(IRECV),ITAG,0,LLWAIT,LLEXIST,IERR)
        IF( IERR < 0 )CALL ABOR1('TRGTOL: ERROR IN MPE_PROBE WAIT=F')
      ELSE
        LLEXIST = .TRUE.
      ENDIF

! If a message exists, receive it, otherwise flag an outstanding receive

      IF( LLEXIST )THEN

!*   receive message

        LLDONE = .FALSE.

        CALL MPE_RECV(ZCOMBUFR(-1),ICOMBFLEN+2,MREALT,&
         &NPRCIDS(IRECV),ITAG,0,0,0,ILRECV,IRECVD,IRCVTG,IERR)
        IF( IERR < 0 ) CALL ABOR1('TRGTOL: ERROR IN MPE_RECV')

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
        IRCVD(J) = IRCVD(J)+IPOS
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
  CALL MPE_BARRIER(IERR)
  IF( IERR /= 0 )THEN
    CALL ABOR1('TRGTOL: ERROR IN MPE_BARRIER')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE TRGTOL
END MODULE TRGTOL_MOD
