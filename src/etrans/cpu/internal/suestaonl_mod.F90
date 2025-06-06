! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE SUESTAONL_MOD
CONTAINS
SUBROUTINE SUESTAONL(KMEDIAP,KRESTM,LDWEIGHTED_DISTR,PWEIGHT,PMEDIAP,KPROCAGP)

!**** *SUESTAONL * - Routine to initialize parallel environment, TAL

!     Purpose.
!     --------
!           Initialize D%NSTA and D%NONL.
!           Calculation of distribution of grid points to processors :
!           Splitting of grid in B direction

!**   Interface.
!     ----------
!        *CALL* *SUESTAONL *

!        Explicit arguments :
!        --------------------
!                     KMEDIAP - mean number of grid points per PE
!                     KRESTM  - number of PEs with one extra point
!                     LDWEIGHTED_DISTR -true if weighted distribution
!                     PWEIGHT    -weight per grid-point if weighted
!                                   distribution
!                     PMEDIAP    -mean weight per PE if weighted
!                                   distribution
!                     KPROCAGP   -number of grid points per A set
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
!                 03-03-03 G. Radnoti: no merge: only difference with
!                                      sustaonl: ezone added to last a-set
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        O.Spaniel     Oct-2004 phasing for AL29
!        A.Bogatchev   Sep-2010 phasing for AL37
!        R. El Khatib 09-Aug-2013 Allow LEQ_REGIONS
!        R. El Khatib 26-Apr-2018 vectorization
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MPL_MODULE  ,ONLY : MPL_ALLGATHERV, MPL_RECV, MPL_SEND

USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_DISTR       ,ONLY : D, LEQ_REGIONS, MTAGPART, NPRCIDS, MYPROC, NPROC
USE TPMALD_DIM      ,ONLY : RALD
USE SET2PE_MOD      ,ONLY : SET2PE
USE EQ_REGIONS_MOD  ,ONLY : MY_REGION_EW, MY_REGION_NS,           &
     &                      N_REGIONS, N_REGIONS_NS, N_REGIONS_EW
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(IN) :: KRESTM
REAL(KIND=JPRD),INTENT(IN)    :: PWEIGHT(:)
LOGICAL,INTENT(IN)            :: LDWEIGHTED_DISTR
REAL(KIND=JPRD),INTENT(IN)    :: PMEDIAP
INTEGER(KIND=JPIM),INTENT(IN) :: KPROCAGP(:)

INTEGER(KIND=JPIM) :: IXPTLAT(R%NDGL), ILSTPTLAT(R%NDGL)
INTEGER(KIND=JPIM) :: ICHK(R%NDLON,R%NDGL), ICOMBUF(R%NDGL*N_REGIONS_EW*2)

INTEGER(KIND=JPIM) :: I1, I2, IBUFLEN, IDGLG, IDWIDE, &
             & IGL, IGL1, IGL2, IGLOFF, IGPTA, &
             & IGPTPRSETS, IGPTS, IGPTSP, ILEN, ILRECV, &
             & ILSEND, INPLAT, INXLAT, IPOS, &
             & IPROCB, IPTSRE, IRECV, &
             & IREST, ISEND, ITAG, JA, JB, JGL, JL, JNPTSRE, &
             & ILAT, ILON, ILOEN  
INTEGER(KIND=JPIM),ALLOCATABLE :: ICOMBUFG(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZWEIGHT(:,:)
INTEGER(KIND=JPIM) :: JJ, ILENG(NPROC), IOFF(NPROC)

LOGICAL :: LLABORT
LOGICAL :: LLP1,LLP2

REAL(KIND=JPRB) ::  ZLAT, ZLAT1(R%NDGL), ZCOMP
REAL(KIND=JPRB) :: ZDIVID(R%NDGL),ZXPTLAT(R%NDGL)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!      -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUESTAONL_MOD:SUESTAONL',0,ZHOOK_HANDLE)
IXPTLAT  (:)=999999
ILSTPTLAT(:)=999999
LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1

IDWIDE  = R%NDGL/2
IBUFLEN = R%NDGL*N_REGIONS_EW*2
IDGLG   = R%NDGL

I1 = MAX(   1,D%NFRSTLAT(MY_REGION_NS)-D%NFRSTLOFF)
I2 = MIN(IDGLG,D%NLSTLAT (MY_REGION_NS)-D%NFRSTLOFF)

ILEN = D%NLSTLAT(MY_REGION_NS) - D%NFRSTLAT(MY_REGION_NS)+1

IGPTPRSETS = SUM(G%NLOEN(1:D%NFRSTLAT(MY_REGION_NS)-1))


IF (D%LSPLIT) THEN
  IF( LEQ_REGIONS )THEN
    IGPTA=0
    DO JA=1,MY_REGION_NS-1
      IGPTA = IGPTA + KPROCAGP(JA)
    ENDDO
    IGPTS = KPROCAGP(MY_REGION_NS)
  ELSE
    IF (MY_REGION_NS <= KRESTM.OR.KRESTM == 0) THEN
      IF (MY_REGION_NS < N_REGIONS_NS) THEN
        IGPTS = KMEDIAP
        IGPTA = KMEDIAP*(MY_REGION_NS-1)
      ELSE
        IGPTS = KMEDIAP+SUM(G%NLOEN(RALD%NDGUX+1:R%NDGL))
        IGPTA = KMEDIAP*(MY_REGION_NS-1)
      ENDIF
    ELSE
      IF (MY_REGION_NS < N_REGIONS_NS) THEN
        IGPTS = KMEDIAP-1
        IGPTA = KMEDIAP*KRESTM+IGPTS*(MY_REGION_NS-1-KRESTM)
      ELSE
        IGPTS = KMEDIAP-1+SUM(G%NLOEN(RALD%NDGUX+1:R%NDGL))
        IGPTA = KMEDIAP*KRESTM+(KMEDIAP-1)*(MY_REGION_NS-1-KRESTM)
      ENDIF
    ENDIF
  ENDIF
ELSE
  IGPTA = IGPTPRSETS
  IGPTS = SUM(G%NLOEN(D%NFRSTLAT(MY_REGION_NS):D%NLSTLAT(MY_REGION_NS)))
ENDIF
IGPTSP = IGPTS/N_REGIONS(MY_REGION_NS)
IREST = IGPTS-N_REGIONS(MY_REGION_NS)*IGPTSP
IXPTLAT(1) = IGPTA-IGPTPRSETS+1
ZXPTLAT(1) = REAL(IXPTLAT(1))
ILSTPTLAT(1) = G%NLOEN(D%NFRSTLAT(MY_REGION_NS))
INPLAT = G%NLOEN(D%NFRSTLAT(MY_REGION_NS))-IXPTLAT(1)+1
DO JGL=2,ILEN
  IXPTLAT(JGL) = 1
  ZXPTLAT(JGL) = 1.0_JPRB
  ILSTPTLAT(JGL) =  G%NLOEN(D%NFRSTLAT(MY_REGION_NS)+JGL-1)
  INPLAT = INPLAT+G%NLOEN(D%NFRSTLAT(MY_REGION_NS)+JGL-1)
ENDDO
ILSTPTLAT(ILEN) = G%NLOEN(D%NLSTLAT(MY_REGION_NS))-INPLAT+IGPTS

DO JB=1,N_REGIONS_EW
  DO JGL=1,R%NDGL+N_REGIONS_NS-1
    D%NSTA(JGL,JB) = 0
    D%NONL(JGL,JB) = 0
  ENDDO
ENDDO

!  grid point decomposition
!  ---------------------------------------
DO JGL=1,ILEN
  ZDIVID(JGL)=1._JPRB/REAL(G%NLOEN(D%NFRSTLAT(MY_REGION_NS)+JGL-1),JPRB)
ENDDO
IF( LDWEIGHTED_DISTR )THEN
  ALLOCATE(ZWEIGHT(G%NLOEN(R%NDGL/2),R%NDGL))
  IGL=0
  DO JGL=1,R%NDGL
    DO JL=1,G%NLOEN(JGL)
      IGL=IGL+1
      ZWEIGHT(JL,JGL)=PWEIGHT(IGL)
    ENDDO
  ENDDO
  ZCOMP=0
  IGPTS=0
ENDIF
DO JB=1,N_REGIONS(MY_REGION_NS)

 IF( .NOT.LDWEIGHTED_DISTR )THEN

  IF (JB <= IREST) THEN
    IPTSRE = IGPTSP+1
  ELSE
    IPTSRE = IGPTSP
  ENDIF

  DO JNPTSRE=1,IPTSRE
    ZLAT  = 1._JPRB
    DO JGL=1,ILEN
      ZLAT1(JGL)  = (ZXPTLAT(JGL)-1.0_JPRB)*ZDIVID(JGL)
    ENDDO
    DO JGL=1,ILEN
      IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
        IF (ZLAT1(JGL) < ZLAT) THEN
         ZLAT=ZLAT1(JGL)
         INXLAT = JGL
        ENDIF
      ENDIF
    ENDDO
    IF (INXLAT >= I1 .AND. INXLAT <= I2) THEN
      IGL=D%NPTRFLOFF+INXLAT
      IF (D%NSTA(IGL,JB) == 0) THEN
        D%NSTA(IGL,JB) = IXPTLAT(INXLAT)
      ENDIF
      D%NONL(IGL,JB) = D%NONL(IGL,JB)+1
    ENDIF
    IXPTLAT(INXLAT) = IXPTLAT(INXLAT)+1
    ZXPTLAT(INXLAT) = REAL(IXPTLAT(INXLAT),JPRB)
  ENDDO

 ELSE
    DO WHILE ( (JB <  N_REGIONS(MY_REGION_NS) .AND. ZCOMP < PMEDIAP) &
        & .OR. (JB == N_REGIONS(MY_REGION_NS) .AND. IGPTS < KPROCAGP(MY_REGION_NS)) )

      IGPTS = IGPTS + 1
      ZLAT  = 1._JPRB
      DO JGL=1,ILEN
        ZLAT1(JGL) = (ZXPTLAT(JGL)-1.0_JPRB)*ZDIVID(JGL)
      ENDDO
      DO JGL=1,ILEN
        IF (IXPTLAT(JGL)  <=  ILSTPTLAT(JGL)) THEN
          IF (ZLAT1(JGL) < ZLAT) THEN
            ZLAT   = ZLAT1(JGL)
            INXLAT = JGL
          ENDIF
        ENDIF
      ENDDO
  
      IF (INXLAT >= I1 .AND. INXLAT <= I2) THEN
        IGL=D%NPTRFLOFF+INXLAT
        IF (D%NSTA(IGL,JB) == 0) THEN
          D%NSTA(IGL,JB) = IXPTLAT(INXLAT)
        ENDIF
        D%NONL(IGL,JB) = D%NONL(IGL,JB)+1
        IF(IGL<1.OR.IGL>R%NDGL+N_REGIONS_NS-1)THEN
          CALL ABORT_TRANS(' SUSTAONL: IGL<1.OR.IGL>R%NDGL+N_REGIONS_NS-1')
        ENDIF
        ILON=D%NSTA(IGL,JB)+D%NONL(IGL,JB)-1
        ILAT=D%NFRSTLAT(MY_REGION_NS)+INXLAT-1
        ILOEN=G%NLOEN(ILAT)
        IF(ILON<1.OR.ILON>ILOEN)THEN
          CALL ABORT_TRANS(' SUSTAONL: ILON<1.OR.ILON>ILOEN')
        ENDIF
        ZCOMP = ZCOMP + ZWEIGHT(ILON,ILAT)
      ENDIF
      IXPTLAT(INXLAT) = IXPTLAT(INXLAT)+1
      ZXPTLAT(INXLAT) = REAL(IXPTLAT(INXLAT),JPRB)
    ENDDO

    ZCOMP = ZCOMP - PMEDIAP

  ENDIF

ENDDO

IF( LDWEIGHTED_DISTR )THEN
  DEALLOCATE(ZWEIGHT)
ENDIF
! Exchange local partitioning info to produce global view

IF( NPROC > 1 )THEN
  IF( LEQ_REGIONS )THEN

    ITAG = MTAGPART
    IPOS = 0
    DO JGL=1,D%NLSTLAT(MY_REGION_NS)-D%NFRSTLAT(MY_REGION_NS)+1
      IPOS = IPOS+1
      ICOMBUF(IPOS) = D%NSTA(D%NPTRFLOFF+JGL,MY_REGION_EW)
      IPOS = IPOS+1
      ICOMBUF(IPOS) = D%NONL(D%NPTRFLOFF+JGL,MY_REGION_EW)
    ENDDO
    IF( IPOS > IBUFLEN )THEN
      CALL ABORT_TRANS(' SUSTAONL: SEND BUFFER TOO SMALL FOR GLOBAL INFO')
    ENDIF
    ILSEND = IPOS

    DO JA=1,N_REGIONS_NS
      DO JB=1,N_REGIONS(JA)
        CALL SET2PE(IRECV,JA,JB,0,0)
        ILEN = (D%NLSTLAT(JA)-D%NFRSTLAT(JA)+1)*2
        ILENG(NPRCIDS(IRECV))=ILEN
      ENDDO
    ENDDO
    IOFF(1)=0
    DO JJ=2,NPROC
      IOFF(JJ)=IOFF(JJ-1)+ILENG(JJ-1)
    ENDDO
    ALLOCATE(ICOMBUFG(SUM(ILENG(:))))
    CALL MPL_ALLGATHERV(ICOMBUF(1:ILSEND),ICOMBUFG,ILENG,CDSTRING='SUSTAONL')
    DO JA=1,N_REGIONS_NS
      IGL1 = D%NFRSTLAT(JA)
      IGL2 = D%NLSTLAT(JA)
      DO JB=1,N_REGIONS(JA)
        CALL SET2PE(IRECV,JA,JB,0,0)
        IF(IRECV /= MYPROC) THEN
          ILEN = (D%NLSTLAT(JA)-D%NFRSTLAT(JA)+1)*2
          IPOS = IOFF(NPRCIDS(IRECV))
          DO JGL=IGL1,IGL2
            IGL = D%NPTRFRSTLAT(JA)+JGL-IGL1
            IPOS = IPOS+1
            D%NSTA(IGL,JB) = ICOMBUFG(IPOS)
            IPOS = IPOS+1
            D%NONL(IGL,JB) = ICOMBUFG(IPOS)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(ICOMBUFG)

  ELSE

    ITAG = MTAGPART
    IPOS = 0
    DO JB=1,N_REGIONS(MY_REGION_NS)
      DO JGL=1,D%NLSTLAT(MY_REGION_NS)-D%NFRSTLAT(MY_REGION_NS)+1
        IPOS = IPOS+1
        ICOMBUF(IPOS) = D%NSTA(D%NPTRFLOFF+JGL,JB)
        IPOS = IPOS+1
        ICOMBUF(IPOS) = D%NONL(D%NPTRFLOFF+JGL,JB)
      ENDDO
    ENDDO
    IF( IPOS > IBUFLEN )THEN
      CALL ABORT_TRANS(' SUESTAONL: SEND BUFFER TOO SMALL FOR GLOBAL INFO')
    ENDIF
    ILSEND = IPOS

    DO JA=1,N_REGIONS_NS
      CALL SET2PE(ISEND,JA,MY_REGION_EW,0,0)
      IF(ISEND /= MYPROC) THEN
        CALL MPL_SEND(ICOMBUF(1:ILSEND),KDEST=NPRCIDS(ISEND),KTAG=ITAG, &
         & CDSTRING='SUESTAONL:') 
      ENDIF
    ENDDO
    DO JA=1,N_REGIONS_NS
      CALL SET2PE(IRECV,JA,MY_REGION_EW,0,0)
      IF(IRECV /= MYPROC) THEN
        ILEN = (D%NLSTLAT(JA)-D%NFRSTLAT(JA)+1)*N_REGIONS(JA)*2
        CALL MPL_RECV(ICOMBUF(1:ILEN),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
          & KOUNT=ILRECV,CDSTRING='SUESTAONL:')  
        IGL1 = D%NFRSTLAT(JA)
        IGL2 = D%NLSTLAT(JA)
        IPOS = 0
        DO JB=1,N_REGIONS(JA)
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
ENDIF

! Confirm consistency of global partitioning, specifically testing for
! multiple assignments of same grid point and unassigned grid points

LLABORT = .FALSE.
DO JGL=1,R%NDGL
  DO JL=1,G%NLOEN(JGL)
    ICHK(JL,JGL) = 1
  ENDDO
ENDDO
DO JA=1,N_REGIONS_NS
  IGLOFF = D%NPTRFRSTLAT(JA)
  DO JB=1,N_REGIONS(JA)
    IGL1 = D%NFRSTLAT(JA)
    IGL2 = D%NLSTLAT(JA)
    DO JGL=IGL1,IGL2
      IGL = IGLOFF+JGL-IGL1
      DO JL=D%NSTA(IGL,JB),D%NSTA(IGL,JB)+D%NONL(IGL,JB)-1
        IF( ICHK(JL,JGL) /= 1 )THEN
          WRITE(NOUT,'(" SUESTAONL : seta=",i4," setb=",i4,&
           & " row=",I4," sta=",I4," INVALID GRID POINT")')&
           & JA,JB,JGL,JL  
          WRITE(0,'(" SUESTAONL : seta=",i4," setb=",i4,&
           & " ROW=",I4," sta=",I4," INVALID GRID POINT")')&
           & JA,JB,JGL,JL  
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
      WRITE(NOUT,'(" SUESTAONL : row=",i4," sta=",i4,&
       & " GRID POINT NOT ASSIGNED")') JGL,JL  
      LLABORT = .TRUE.
    ENDIF
  ENDDO
ENDDO
IF( LLABORT )THEN
  WRITE(NOUT,'(" SUESTAONL : inconsistent partitioning")')
  CALL ABORT_TRANS(' SUESTAONL: inconsistent partitioning')
ENDIF

IF (LLP1) THEN
  WRITE(UNIT=NOUT,FMT='('' OUTPUT FROM ROUTINE SUESTAONL '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
  WRITE(UNIT=NOUT,FMT='('' PARTITIONING INFORMATION '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
  IPROCB = MIN(32,N_REGIONS_EW)
  WRITE(UNIT=NOUT,FMT='(17X," SETB=",32(1X,I3))') (JB,JB=1,IPROCB)
  DO JA=1,N_REGIONS_NS
    IPROCB = MIN(32,N_REGIONS(JA))
    WRITE(UNIT=NOUT,FMT='('' '')')
    IGLOFF = D%NPTRFRSTLAT(JA)
    IGL1 = D%NFRSTLAT(JA)
    IGL2 = D%NLSTLAT(JA)
    DO JGL=IGL1,IGL2
      IGL=IGLOFF+JGL-IGL1
      WRITE(UNIT=NOUT,FMT='(" SETA=",I3," LAT=",I3," NSTA=",&
       & 32(1X,I3))') JA,JGL,(D%NSTA(IGL,JB),JB=1,IPROCB)  
      WRITE(UNIT=NOUT,FMT='(" SETA=",I3," LAT=",I3," D%NONL=",&
       & 32(1X,I3))') JA,JGL,(D%NONL(IGL,JB),JB=1,IPROCB)  
      WRITE(UNIT=NOUT,FMT='('' '')')
    ENDDO
    WRITE(UNIT=NOUT,FMT='('' '')')
  ENDDO
  WRITE(UNIT=NOUT,FMT='('' '')')
  WRITE(UNIT=NOUT,FMT='('' '')')
ENDIF
IF (LHOOK) CALL DR_HOOK('SUESTAONL_MOD:SUESTAONL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUESTAONL
END MODULE SUESTAONL_MOD
