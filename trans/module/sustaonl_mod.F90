MODULE SUSTAONL_MOD
CONTAINS
SUBROUTINE SUSTAONL(KMEDIAP,KRESTM)

!**** *SUSTAONL * - Routine to initialize parallel environment

!     Purpose.
!     --------
!           Initialize D%NSTA and D%NONL.
!           Calculation of distribution of grid points to processors :
!           Splitting of grid in B direction

!**   Interface.
!     ----------
!        *CALL* *SUSTAONL *

!        Explicit arguments : KMEDIAP - mean number of grid points per PE
!        -------------------- KRESTM  - number of PEs with one extra point

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option.
!          - removal of LRPOLE in YOMCT0.
!          - removal of code under LRPOLE.
!        Modified 98-12-04 C. Fischer: merge with SUESTAONL (Aladin)
!     ------------------------------------------------------------------

#include "tsmbkind.h"
USE MPL_MODULE

USE TPM_GEN
USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_DISTR

USE SET2PE_MOD

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KMEDIAP
INTEGER_M :: KRESTM


INTEGER_M :: IXPTLAT(R%NDGL), ILSTPTLAT(R%NDGL), INTOTRE(NPRGPEW)
INTEGER_M :: ICHK(R%NDLON,R%NDGL), ICOMBUF(R%NDGL*NPRGPEW*2)

!     LOCAL INTEGER SCALARS
INTEGER_M :: I1, I2, IB, IBS, IBUFLEN, IDGLG, IDWIDE,&
             &IERR, IGL, IGL1, IGL2, IGLOFF, IGPTA, IGPTOT, &
             &IGPTPRSETS, IGPTS, IGPTSP, ILEN, ILRECV, &
             &ILSEND, INPLAT, INXLAT, IPART, IPONL, IPOS, &
             &IPROCB, IPSTA, IPTSRE, IRCVID, IRCVTG, IRECV, &
             &IREST, ISEND, ITAG, JA, JB, JGL, JL, JNPTSRE, &
             &JS

!     LOCAL LOGICAL SCALARS
LOGICAL :: LLABORT, LLALLAT
LOGICAL :: LLP1,LLP2

!     LOCAL REAL SCALARS
REAL_B ::  ZLAT, ZLAT1

REAL_B :: ZDIVID(R%NDGL),ZXPTLAT(R%NDGL)
!      -----------------------------------------------------------------

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1


IDWIDE  = R%NDGL/2
IBUFLEN = R%NDGL*NPRGPEW*2
IDGLG   = R%NDGL


I1 = MAX(   1,D%NFRSTLAT(MYSETNS)-D%NFRSTLOFF)
I2 = MIN(IDGLG,D%NLSTLAT (MYSETNS)-D%NFRSTLOFF)

ILEN = D%NLSTLAT(MYSETNS) - D%NFRSTLAT(MYSETNS)+1

IGPTPRSETS = SUM(G%NLOEN(1:D%NFRSTLAT(MYSETNS)-1))

IGPTOT = SUM(G%NLOEN(1:R%NDGL))

IF (D%LSPLIT) THEN
  IF (MYSETNS <= KRESTM.OR.KRESTM == 0) THEN
    IGPTS = KMEDIAP
    IGPTA = KMEDIAP*(MYSETNS-1)
  ELSE
    IGPTS = KMEDIAP-1
    IGPTA = KMEDIAP*KRESTM+IGPTS*(MYSETNS-1-KRESTM)
  ENDIF
ELSE
  IGPTA = IGPTPRSETS
  IGPTS = SUM(G%NLOEN(D%NFRSTLAT(MYSETNS):D%NLSTLAT(MYSETNS)))
ENDIF
IGPTSP = IGPTS/NPRGPEW
IREST = IGPTS-NPRGPEW*IGPTSP
IXPTLAT(1) = IGPTA-IGPTPRSETS+1
ZXPTLAT(1) = REAL(IXPTLAT(1))
ILSTPTLAT(1) = G%NLOEN(D%NFRSTLAT(MYSETNS))
INPLAT = G%NLOEN(D%NFRSTLAT(MYSETNS))-IXPTLAT(1)+1
DO JGL=2,ILEN
  IXPTLAT(JGL) = 1
  ZXPTLAT(JGL) = _ONE_
  ILSTPTLAT(JGL) =  G%NLOEN(D%NFRSTLAT(MYSETNS)+JGL-1)
  INPLAT = INPLAT+G%NLOEN(D%NFRSTLAT(MYSETNS)+JGL-1)
ENDDO
ILSTPTLAT(ILEN) = G%NLOEN(D%NLSTLAT(MYSETNS))-INPLAT+IGPTS

DO JB=1,NPRGPEW
  DO JGL=1,R%NDGL+NPRGPNS-1
    D%NSTA(JGL,JB) = 0
    D%NONL(JGL,JB) = 0
  ENDDO
ENDDO


!  grid point decomposition
!  ---------------------------------------
LLALLAT = (NPRGPNS == 1)
DO JGL=1,ILEN
  ZDIVID(JGL)=REAL(G%NLOEN(D%NFRSTLAT(MYSETNS)+JGL-1),JPRB)
ENDDO
DO JB=1,NPRGPEW

  IF (JB <= IREST) THEN
    IPTSRE = IGPTSP+1
  ELSE
    IPTSRE = IGPTSP
  ENDIF

  IPART=0
  DO JNPTSRE=1,IPTSRE
    ZLAT  = 1._JPRB
    ZLAT1 = 1._JPRB
    IF (MYSETNS <= D%NAPSETS .AND.(IPART /= 2.OR.LLALLAT)) THEN
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1  = (ZXPTLAT(JGL)-_ONE_)/ZDIVID(JGL)
          ZLAT   = MIN(ZLAT1,ZLAT)
          INXLAT = JGL
          IPART  = 1
          EXIT
        ENDIF
      ENDDO
    ELSEIF (MYSETNS > NPRGPNS-D%NAPSETS.AND.(IPART /= 1.OR.LLALLAT)) THEN
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1  = (ZXPTLAT(JGL)-_ONE_)/ZDIVID(JGL)
          ZLAT   = MIN(ZLAT1,ZLAT)
          INXLAT = JGL
          IPART  = 2
          EXIT
        ENDIF
      ENDDO
    ELSE
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          ZLAT1 = (ZXPTLAT(JGL)-_ONE_)/ZDIVID(JGL)
          IF (ZLAT1 < ZLAT) THEN
            ZLAT   = ZLAT1
            INXLAT = JGL
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF (INXLAT >= I1 .AND. INXLAT <= I2) THEN
      IF (D%NSTA(D%NPTRFLOFF+INXLAT,JB) == 0) THEN
        D%NSTA(D%NPTRFLOFF+INXLAT,JB) = IXPTLAT(INXLAT)
      ENDIF
      D%NONL(D%NPTRFLOFF+INXLAT,JB) = D%NONL(D%NPTRFLOFF+INXLAT,JB)+1
    ENDIF
    IXPTLAT(INXLAT) = IXPTLAT(INXLAT)+1
    ZXPTLAT(INXLAT) = REAL(IXPTLAT(INXLAT),JPRB)
  ENDDO
ENDDO


! Exchange local partitioning info to produce global view
!

IF( NPROC > 1 )THEN
  ITAG = MTAGPART
  IPOS = 0
  DO JB=1,NPRGPEW
    DO JGL=1,D%NLSTLAT(MYSETNS)-D%NFRSTLAT(MYSETNS)+1
      IPOS = IPOS+1
      ICOMBUF(IPOS) = D%NSTA(D%NPTRFLOFF+JGL,JB)
      IPOS = IPOS+1
      ICOMBUF(IPOS) = D%NONL(D%NPTRFLOFF+JGL,JB)
    ENDDO
  ENDDO
  IF( IPOS > IBUFLEN )THEN
    CALL ABOR1(' SUSTAONL: SEND BUFFER TOO SMALL FOR GLOBAL INFO')
  ENDIF
  ILSEND = IPOS

  DO JA=1,NPRGPNS
    CALL SET2PE(ISEND,JA,MYSETEW,0,0)
    IF(ISEND /= MYPROC) THEN
      CALL MPL_SEND(ICOMBUF(1:ILSEND),KDEST=NPRCIDS(ISEND),KTAG=ITAG, &
       &   CDSTRING='SUSTAONL:')
    ENDIF
  ENDDO

  DO JA=1,NPRGPNS
    ILEN = (D%NLSTLAT(JA)-D%NFRSTLAT(JA)+1)*NPRGPEW*2
    CALL SET2PE(IRECV,JA,MYSETEW,0,0)
    IF(IRECV /= MYPROC) THEN
      CALL MPL_RECV(ICOMBUF(1:ILEN),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
       & KOUNT=ILRECV,CDSTRING='SUSTAONL:')
      IGL1 = D%NFRSTLAT(JA)
      IGL2 = D%NLSTLAT(JA)
      IPOS = 0
      DO JB=1,NPRGPEW
        DO JGL=IGL1,IGL2
          IGL = D%NPTRFRSTLAT(JA)+JGL-IGL1
          IPOS = IPOS+1
          D%NSTA(IGL,JB) = ICOMBUF(IPOS)
          IPOS = IPOS+1
          D%NONL(IGL,JB) = ICOMBUF(IPOS)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDIF

! Confirm consistency of global partitioning, specifically testing for
! multiple assignments of same grid point and unassigned grid points

LLABORT = .FALSE.
DO JGL=1,R%NDGL
  DO JL=1,G%NLOEN(JGL)
    ICHK(JL,JGL) = 1
  ENDDO
ENDDO
DO JA=1,NPRGPNS
  IGLOFF = D%NPTRFRSTLAT(JA)
  DO JB=1,NPRGPEW
    IGL1 = D%NFRSTLAT(JA)
    IGL2 = D%NLSTLAT(JA)
    DO JGL=IGL1,IGL2
      IGL = IGLOFF+JGL-IGL1
      DO JL=D%NSTA(IGL,JB),D%NSTA(IGL,JB)+D%NONL(IGL,JB)-1
        IF( ICHK(JL,JGL) /= 1 )THEN
          WRITE(NOUT,'(" SUSTAONL : seta=",i4," setb=",i4,&
           &" row=",I4," sta=",I4," INVALID GRID POINT")')&
           &JA,JB,JGL,JL
          WRITE(0,'(" SUSTAONL : seta=",i4," setb=",i4,&
           &" ROW=",I4," sta=",I4," INVALID GRID POINT")')&
           &JA,JB,JGL,JL
          LLABORT = .TRUE.
        ENDIF
        ICHK(JL,JGL) = 2
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JGL=1,R%NDGL
  DO JL=1,G%NLOEN(JGL)
    IF( ICHK(JL,JGL) /= 2 )THEN
      WRITE(NOUT,'(" SUSTAONL : row=",i4," sta=",i4,&
       &" GRID POINT NOT ASSIGNED")') JGL,JL
      LLABORT = .TRUE.
    ENDIF
  ENDDO
ENDDO
IF( LLABORT )THEN
  WRITE(NOUT,'(" SUSTAONL : inconsistent partitioning")')
  CALL ABOR1(' SUSTAONL: inconsistent partitioning')
ENDIF


IF (LLP1) THEN
  WRITE(UNIT=NOUT,FMT='('' OUTPUT FROM ROUTINE SUSTAONL '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
  WRITE(UNIT=NOUT,FMT='('' PARTITIONING INFORMATION '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
  IPROCB = MIN(32,NPRGPEW)
  WRITE(UNIT=NOUT,FMT='(17X," SETB=",32(1X,I3))') (JB,JB=1,IPROCB)
  DO JA=1,NPRGPNS
    WRITE(UNIT=NOUT,FMT='('' '')')
    IGLOFF = D%NPTRFRSTLAT(JA)
    IGL1 = D%NFRSTLAT(JA)
    IGL2 = D%NLSTLAT(JA)
    DO JGL=IGL1,IGL2
      IGL=IGLOFF+JGL-IGL1
      WRITE(UNIT=NOUT,FMT='(" SETA=",I3," LAT=",I3," NSTA=",&
       &32(1X,I3))') JA,JGL,(D%NSTA(IGL,JB),JB=1,IPROCB)
      WRITE(UNIT=NOUT,FMT='(" SETA=",I3," LAT=",I3," D%NONL=",&
       &32(1X,I3))') JA,JGL,(D%NONL(IGL,JB),JB=1,IPROCB)
      WRITE(UNIT=NOUT,FMT='('' '')')
    ENDDO
    WRITE(UNIT=NOUT,FMT='('' '')')
  ENDDO
  WRITE(UNIT=NOUT,FMT='('' '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUSTAONL
END MODULE SUSTAONL_MOD

