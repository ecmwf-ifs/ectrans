! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ESPLINE_MOD
CONTAINS
SUBROUTINE ESPLINE(KDLUN,KDLON,KDGUN,KDGL,KDLUX,KDGUX,KSTART,&
 & KDLSM,KDGSA,KDGEN,KNUBI,PWORK,PALFA,LDBIX,LDBIY,KDAD)  

!   purpose  :
!   --------
!     Make spline extension.

!    *CALL* *ESPLINE*(...)

!      externals :
!      ----------
!             None

!      explicit arguments :
!      ------------------
!     KDLUN : lower bound for the x (or longitude) dimension
!             of the gridpoint array
!     KDLON : upper bound for the x (or longitude) dimension
!             of the gridpoint array on C U I U E
!     KDGUN : lower bound for the y (or latitude) dimension
!             of the gridpoint array
!     KDGL : upper bound for the y (or latitude) dimension
!             of the gridpoint array on C U I U E
!     KDLUX : upper bound for the x (or longitude) dimension
!             of  C U I.
!     KDGUX : upper bound for the y (or latitude) dimension
!             of  C U I.
!     KSTART : first dimension in x direction of g-p array
!     KDLSM  : last dimension in x direction of g-p array
!     KDGSA  : first dimension in y of g-p array
!     KDGEN  :  last dimension in y of g-p array
!     KNUBI  : number of levels to biperiodicise
!     PWORK : gridpoint array on C U I U E.
!     PALFA : boundary condition of a spline:
!             = 0. ... natural spline
!             = 1. ... boundary condition computed differentially
!                      (additional option)
!     LDBIX : .TRUE. biperiodicisation in x  ( and vice versa )
!     LDBIY : .TRUE. biperiodicisation in y  ( and vice versa )
!     KDAD  : 1 for test of biperiodic.

!      references :
!      ----------

!      author :
!      ------
!       Michal Batka and Radmila Bubnova ( B & B )

!      modifications :
!      -------------
!      J.Vivoda 03-2002 2D model fix
!      A. Stanesic  : 28-03-08: KDADD - test of externalized biper.
!      -------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!      -------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLSM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUBI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWORK(KSTART:KDLSM,KNUBI,KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALFA 
LOGICAL           ,INTENT(IN)    :: LDBIX 
LOGICAL           ,INTENT(IN)    :: LDBIY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDAD 

!      -------------------------------------------------------------

LOGICAL :: LLBIX
LOGICAL :: LLBIY
INTEGER(KIND=JPIM) :: IENDX, IENDY, JFL, JLAT, JLON, IA
REAL(KIND=JPRB) :: ZA, ZB, ZC, ZD, ZEPSA, ZEPSB, ZJ, ZK, ZKP1,&
 & ZLAM, ZLAMB, ZM1, ZM2, ZMM, ZNY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ESPLINE',0,ZHOOK_HANDLE)
!      -------------------------------------------------------------

!*         1. Spline   Extension.
!             -------------------

LLBIX=LDBIX
LLBIY=LDBIY

IF( KDLUN==1.AND.KDLUX==1 ) LLBIX=.FALSE.
IF( KDGUN==1.AND.KDGUX==1 ) LLBIY=.FALSE.

IENDX = KDGUX 
IENDY = KDLON 

IF(LLBIX.AND.(.NOT.LLBIY)) THEN
  IENDY = KDLUN - 1

ELSEIF((.NOT.LLBIX).AND.LLBIY) THEN
  IENDX = KDGUN - 1
  IENDY = KDLUX

ELSEIF((.NOT.LLBIX).AND.(.NOT.LLBIY)) THEN
  IF (LHOOK) CALL DR_HOOK('ESPLINE',1,ZHOOK_HANDLE)
  RETURN
ENDIF
DO JFL = 1, KNUBI

  ZK    = REAL(KDLON-KDLUX+1,JPRB)
  ZKP1  = ZK + 1.0_JPRB
  ZLAMB = ZK/ZKP1
  ZNY   = PALFA/ZKP1

  DO JLAT=KDGUN,IENDX

    ZEPSA = ((PWORK(KDLUN,JFL,JLAT)-PWORK(KDLUX,JFL,JLAT))/ZK -&
     & PWORK(KDLUX,JFL,JLAT)+PWORK(KDLUX-1,JFL,JLAT))*6._JPRB/ZKP1 -&
     & ZNY*(PWORK(KDLUX,JFL,JLAT)-2.0_JPRB* PWORK(KDLUX-1,JFL,JLAT)+&
     & PWORK(KDLUX-2,JFL,JLAT))  

    ZEPSB = (PWORK(KDLUN+1,JFL,JLAT)-PWORK(KDLUN,JFL,JLAT) -&
     & (PWORK(KDLUN,JFL,JLAT)-PWORK(KDLUX,JFL,JLAT))/ZK)*6._JPRB/ZKP1-&
     & ZNY*(PWORK(KDLUN+2,JFL,JLAT)-2.0_JPRB* PWORK(KDLUN+1,JFL,JLAT)+&
     & PWORK(KDLUN,JFL,JLAT))  

    ZMM = 4._JPRB - ZLAMB*ZLAMB
    ZM1 = (2.0_JPRB*ZEPSA - ZLAMB*ZEPSB)/ZMM
    ZM2 = (2.0_JPRB*ZEPSB - ZLAMB*ZEPSA)/ZMM
    ZA  = PWORK(KDLUX,JFL,JLAT)
    ZB  = (PWORK(KDLUN,JFL,JLAT)-PWORK(KDLUX,JFL,JLAT))/ZK-&
     & (2.0_JPRB*ZM1 + ZM2) * ZK/6._JPRB  
    ZC  = 0.5_JPRB * ZM1
    ZD  = (ZM2 - ZM1)/(6._JPRB*ZK)

    DO JLON=KDLUX+1,KDLON+KDAD
      ZJ  = REAL(JLON - KDLUX,JPRB)
      PWORK(JLON,JFL,JLAT) = ZA + ZJ * (ZB + ZJ * (ZC + ZD * ZJ))
    ENDDO 
  ENDDO 

  ZK   = REAL(KDGL - KDGUX + 1,JPRB)
  ZKP1 = ZK + 1
  ZLAM = ZK/ZKP1
  ZNY  = PALFA/ZKP1

  DO JLON=KDLUN,IENDY+KDAD

    ZEPSA = ((PWORK(JLON,JFL,KDGUN)-PWORK(JLON,JFL,KDGUX))/ZK -&
     & PWORK(JLON,JFL,KDGUX)+PWORK(JLON,JFL,KDGUX-1))*6._JPRB/ZKP1-&
     & ZNY*(PWORK(JLON,JFL,KDGUX)-2.0_JPRB*PWORK(JLON,JFL,KDGUX-1)+&
     & PWORK(JLON,JFL,KDGUX-2))  

    ZEPSB = (PWORK(JLON,JFL,KDGUN+1)-PWORK(JLON,JFL,KDGUN) -&
     & (PWORK(JLON,JFL,KDGUN)-PWORK(JLON,JFL,KDGUX))/ZK)*6._JPRB/ZKP1-&
     & ZNY*(PWORK(JLON,JFL,KDGUN+2)-2.0_JPRB*PWORK(JLON,JFL,KDGUN+1) +&
     & PWORK(JLON,JFL,KDGUN))  

    ZMM = 4._JPRB - ZLAMB*ZLAMB
    ZM1 = (2.0_JPRB*ZEPSA - ZLAMB*ZEPSB)/ ZMM
    ZM2 = (2.0_JPRB*ZEPSB - ZLAMB*ZEPSA)/ ZMM
    ZA  = PWORK(JLON,JFL,KDGUX)
    ZB  = (PWORK(JLON,JFL,KDGUN)-PWORK(JLON,JFL,KDGUX))/ZK - (2.0_JPRB*&
     & ZM1 &
     & + ZM2) * ZK/6._JPRB  
    ZC  = 0.5_JPRB * ZM1
    ZD  = (ZM2 - ZM1)/(6._JPRB*ZK)

    DO JLAT=KDGUX+1,KDGL+KDAD
      ZJ = REAL(JLAT - KDGUX,JPRB)
      PWORK(JLON,JFL,JLAT) = ZA +ZJ*(ZB +ZJ*(ZC + ZJ * ZD))
    ENDDO 
  ENDDO 

ENDDO

!      -------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ESPLINE',1,ZHOOK_HANDLE)
END SUBROUTINE ESPLINE
END MODULE ESPLINE_MOD
