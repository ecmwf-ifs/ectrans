MODULE TRGTOL_MOD
CONTAINS
SUBROUTINE TRGTOL(PGLAT,PGCOL,KF_FS,KF_GP,KVSET,KPTRGP)

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
INTEGER_M,INTENT(IN) :: KF_FS,KF_GP
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPTRGP(:)

REAL_B,ALLOCATABLE :: ZCOMBUFS(:,:),ZCOMBUFR(:,:)

INTEGER_M :: ISENT    (NPROC)
INTEGER_M :: IRCVD    (NPROC)
INTEGER_M :: ISENDTOT (NPROC)
INTEGER_M :: IRECVTOT (NPROC)
INTEGER_M :: ISENDREQ (NPROC)
INTEGER_M :: IRECVREQ (NPROC)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IERR, IFIRST, IFIRSTLAT, IGL, IGLL, ILAST,&
             &ILASTLAT, ILEN, ILRECV, JROC, IPOS, ISETA, &
             &ISETB, IRCVTG, IRECV, IRECVD, IRECVSET, &
             &ISETV, ISEND, ISENDSET, ITAG, J, JBLK, JFLD, &
             &JGL, JK, JL, JLOOP, ISETW,  JVAR,IFLD, &
             &II,INDOFFX,IBUFLENS,IBUFLENR,INRECV, IPROC,IFLDS, &
             &INSEND,INS,INR


!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE, LLEXIST,  LLRECV,LLINDER
INTEGER_M :: INDEX(D%NLENGTF),INDOFF(NPROC),IFLDOFF(KF_FS)
INTEGER_M :: ISEND_FLD_TOTAL(NPROC),ISEND_FLD_START(NPROC),ISEND_FLD_END
INTEGER_M :: IRECV_FLD_START,IRECV_FLD_END
INTEGER_M :: ICOMBFLEN
INTEGER_M :: INUMFLDS
INTEGER_M :: IPROCS_COMM_MAX
INTEGER_M :: IGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER_M :: IGPTRRECV(NPRTRNS)
INTEGER_M :: IGPTROFF(NGPBLKS)
!     INTEGER FUNCTIONS
INTEGER_M :: MYSENDSET,MYRECVSET

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------
IF(PRESENT(KPTRGP)) THEN
  LLINDER = .TRUE.
ELSE
  LLINDER = .FALSE.
ENDIF

CALL INIGPTR(IGPTRSEND,IGPTRRECV)
LLDONE = .FALSE.

ITAG = MTAGGL

INDOFFX  = 0
IBUFLENS = 0
IBUFLENR = 0
INRECV   = 0
INSEND   = 0

DO JROC=1,NPROC

  CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)
  IRECVSET = ISETA
  ISEND = JROC
  ISENDSET = ISETV
  ISENT(JROC) = 0
  IRCVD(JROC) = 0

!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1) IPOS = IPOS+1
  ENDDO
  ISEND_FLD_TOTAL(JROC) = IPOS
  ISENDTOT(JROC) = IGPTRRECV(ISETW)*IPOS

  IF( JROC /= MYPROC) THEN
    IBUFLENS = MAX(IBUFLENS,ISENDTOT(JROC))
    IF(ISENDTOT(JROC) > 0) INSEND = INSEND+1
  ENDIF

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(IRECVSET))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(IRECVSET))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(IRECVSET)+JGL-D%NFRSTLAT(IRECVSET)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO
  
  IRECVTOT(JROC) = IPOS*KF_FS

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

IF( .NOT.LIMP) THEN
  INSEND = 1
  INRECV = 1
  INS = 1
  INR = 1
ENDIF

IF (IBUFLENR > 0) ALLOCATE(ZCOMBUFR(-1:ICOMBFLEN,INRECV))
IF (IBUFLENS > 0) ALLOCATE(ZCOMBUFS(-1:ICOMBFLEN,INSEND))


! Send loop.............................................................

ISEND_FLD_START(:) = 1

DO WHILE( .NOT.LLDONE )

  LLDONE = .TRUE.
!  For immediate send/recv, post receives first

  IF(LIMP) THEN
    INR = 0
    DO JLOOP=1,NPROC-1
      IPROC = MYRECVSET(NPROC,MYPROC,JLOOP)

! Check if there is data to receive

      IF( IRECVTOT(IPROC)-IRCVD(IPROC) > 0 )THEN
        IRECV = IPROC
        INR = INR+1

!*   receive message

        CALL MPL_RECV(ZCOMBUFR(-1:ICOMBFLEN,INR),KSOURCE=NPRCIDS(IRECV), &
         & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(INR), &
         & KTAG=ITAG,CDSTRING='TRLTOG:' )
      ENDIF
    ENDDO
  ENDIF
!
! loop over the number of processors we need to communicate with
!
  IF(LIMP) INS = 0

  DO JLOOP=1,NPROC-1
    IPROC = MYSENDSET(NPROC,MYPROC,JLOOP)

! Check if we have data to send

    IF(ISENDTOT(IPROC)-ISENT(IPROC) > 0 )THEN
      LLDONE = .FALSE.

      CALL PE2SET(IPROC,ISETA,ISETB,ISETW,ISETV)
      ISEND = IPROC
      ISENDSET = ISETV
      IF(LIMP) INS = INS+1

      INUMFLDS = ICOMBFLEN/IGPTRRECV(ISETW)

      ISEND_FLD_END = MIN(ISEND_FLD_TOTAL(IPROC),&
       &ISEND_FLD_START(IPROC)+INUMFLDS-1)

      ZCOMBUFS(-1,INS) = ISEND_FLD_START(IPROC)
      IFLD = 0
      IPOS = 0
      DO JFLD=1,KF_GP
        IF(KVSET(JFLD) == ISENDSET .OR. KVSET(JFLD) == -1 ) THEN

          IF (IFLD == ISEND_FLD_END) EXIT
          IFLD = IFLD+1
          IF (IFLD >= ISEND_FLD_START(IPROC)) THEN
            
            DO JBLK=1,NGPBLKS
              IFIRST = IGPTRSEND(1,JBLK,ISETW)
              IF(IFIRST > 0) THEN
                ILAST = IGPTRSEND(2,JBLK,ISETW)
                IF(LLINDER) THEN
                  DO JK=IFIRST,ILAST
                    IPOS = IPOS+1
                    ZCOMBUFS(IPOS,INS) = PGCOL(JK,KPTRGP(JFLD),JBLK)
                  ENDDO
                ELSE 
                  DO JK=IFIRST,ILAST
                    IPOS = IPOS+1
                    ZCOMBUFS(IPOS,INS) = PGCOL(JK,JFLD,JBLK)
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      ZCOMBUFS(0,INS) = IFLD
      
      ISEND_FLD_START(IPROC) = IFLD+1

      IF(LIMP) THEN
        CALL MPL_SEND(ZCOMBUFS(-1:IPOS,INS),KDEST=NPRCIDS(ISEND), &
         & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(INS), &
         & KTAG=ITAG,CDSTRING='TRGTOL:' )
      ELSE
        CALL MPL_SEND(ZCOMBUFS(-1:IPOS,INS),KDEST=NPRCIDS(ISEND), &
         &    KTAG=ITAG,CDSTRING='TRGTOL:' )
      ENDIF
      ISENT(IPROC) = ISENT(IPROC)+IPOS

    ENDIF
  ENDDO

  ! Copy local data

  IF(ISENDTOT(MYPROC)-ISENT(MYPROC) > 0 )THEN
    
    IFLDS = 0
    DO JFLD=1,KF_GP
      IF(KVSET(JFLD) == MYSETV .OR. KVSET(JFLD) == -1) THEN
        IFLDS = IFLDS+1
        IF(LLINDER) THEN
          IFLDOFF(IFLDS) = KPTRGP(JFLD)
        ELSE
          IFLDOFF(IFLDS) = JFLD
        ENDIF
      ENDIF
    ENDDO

    IPOS=0
    DO JBLK=1,NGPBLKS
      IGPTROFF(JBLK)=IPOS
      IFIRST = IGPTRSEND(1,JBLK,MYSETW)
      IF(IFIRST > 0) THEN
        ILAST = IGPTRSEND(2,JBLK,MYSETW)
        IPOS=IPOS+ILAST-IFIRST+1
      ENDIF
    ENDDO
!$OMP PARALLEL DO PRIVATE(JFLD,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
    DO JBLK=1,NGPBLKS
      IFIRST = IGPTRSEND(1,JBLK,MYSETW)
      IF(IFIRST > 0) THEN
        ILAST = IGPTRSEND(2,JBLK,MYSETW)
        DO JFLD=1,IFLDS
          IFLD = IFLDOFF(JFLD)
          DO JK=IFIRST,ILAST
            IPOS = INDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
            PGLAT(JFLD,INDEX(IPOS)) = PGCOL(JK,IFLD,JBLK)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
!$OMP END PARALLEL DO

    ISENT(MYPROC) = ISENDTOT(MYPROC)

  ENDIF

!  Receive loop.........................................................

  IF (LIMP) INR = 0

  DO JLOOP=1,NPROC-1
    IPROC = MYRECVSET(NPROC,MYPROC,JLOOP)
    IF( IRECVTOT(IPROC)-IRCVD(IPROC) > 0 )THEN
      LLDONE = .FALSE.
      CALL PE2SET(IPROC,ISETA,ISETB,ISETW,ISETV)
      IRECVSET = ISETA
      IRECV = IPROC
!*   receive message


      IF (LIMP) THEN
        INR = INR+1
        CALL MPL_WAIT(ZCOMBUFR(-1:ICOMBFLEN,INR),KREQUEST=IRECVREQ(INR), &
         & KCOUNT=ILRECV,CDSTRING='TRGTOL: LIMP WAITS ' )
      ELSE
        CALL MPL_RECV(ZCOMBUFR(-1:ICOMBFLEN,INR),KSOURCE=NPRCIDS(IRECV), &
         & KTAG=ITAG,KOUNT=ILRECV,CDSTRING='TRGTOL:' )
      ENDIF
      IRECV_FLD_START = ZCOMBUFR(-1,INR)
      IRECV_FLD_END   = ZCOMBUFR(0,INR)


!*   store data in Fourier array

      ILEN = IRECVTOT(IRECV)/KF_FS
      DO JL=1,ILEN
        II = INDEX(INDOFF(IRECV)+JL)
        DO JFLD=IRECV_FLD_START,IRECV_FLD_END
          PGLAT(JFLD,II) = ZCOMBUFR(JL+(JFLD-IRECV_FLD_START)*ILEN,INR)
        ENDDO
      ENDDO
      IPOS = ILEN*(IRECV_FLD_END-IRECV_FLD_START+1)
      IF (ILRECV /= IPOS+2) THEN
        WRITE(NOUT,*) MYPROC,MYSETV,' receiving ',ILRECV,'expected',IPOS
        CALL ABOR1('TRGTOL: RECEIVED MESSAGE OF INCORRECT LENGTH')
      ENDIF
      IRCVD(IPROC) = IRCVD(IPROC)+IPOS

    ENDIF
  ENDDO

  IF(LIMP)THEN
    IF(INS > 0) THEN
      CALL MPL_WAIT(ZCOMBUFS(-1:ICOMBFLEN,1),KREQUEST=ISENDREQ(1:INS), &
      & CDSTRING='TRGTOL: ERROR IN MPL_WAIT FOR SENDS')
    ENDIF
  ENDIF

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
