! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ESMOOTHE_MOD
CONTAINS
SUBROUTINE ESMOOTHE(KDLUN,KDLON,KDGUN,KDGL,KDLUX,KDGUX,KSTART,&
 & KDLSM,KDGSA,KDGEN,KNUBI,PWORK,LDBIX,LDBIY)  

!   purpose  :
!   --------
!     To smooth the fields over the extension zone.

!*    *CALL* *ESMOOTHE*(...)

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
!     KDLSM  : dimension in x direction of g-p array
!     KDGSA  : first dimension index in y of g-p array
!     KDGEN  : last dimension index in y of g-p array
!     KSTART :  first dimension index in x of g-p array
!     KDLSM  :  last dimension index in x of g-p array
!     KNUBI  : number of levels to biperiodicise

!     PWORK : gridpoint array on C U I U E.

!     LDBIX  : .TRUE.: biperiodicise in x direction (and vice versa)
!     LDBIY  : .TRUE.: biperiodicise in y direction (and vice versa)

!      references :
!      ----------

!      author :
!      ------
!       Michal Batka and Radmila Bubnova ( B & B )

!      modifications :
!      -------------
!      R. El Khatib 03-05-05 Optimizations
!      --------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!      --------------------------------------------------------------

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
LOGICAL           ,INTENT(IN)    :: LDBIX 
LOGICAL           ,INTENT(IN)    :: LDBIY 

!      --------------------------------------------------------------

REAL(KIND=JPRB) :: ZPRAC(KDLUN-1:KDLON+1,KDGUN-1:KDGL+1)
INTEGER(KIND=JPIM) :: IEND, IENX1, IENX2, IENY1, IENY2, JFL, JLAT, JLL, JLON
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      --------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ESMOOTHE',0,ZHOOK_HANDLE)
!      --------------------------------------------------------------

!*         1. Calculation.
!             ------------

IEND = MAX(KDLON-KDLUX,KDGL-KDGUX)
IEND = (IEND+1)/2
IENX1= KDLON
IENX2= KDGL
IENY1= KDGL
IENY2= KDLON
IF(LDBIX.AND.(.NOT.LDBIY)) THEN
  IENX2 = KDGUX
  IENY1 = KDGUX
ELSEIF((.NOT.LDBIX).AND.LDBIY) THEN
  IENX1 = KDLUX
  IENY2 = KDLUX
ELSEIF((.NOT.LDBIX).AND.(.NOT.LDBIY)) THEN
  IF (LHOOK) CALL DR_HOOK('ESMOOTHE',1,ZHOOK_HANDLE)
  RETURN
ENDIF

DO JFL = 1, KNUBI

  DO JLL = 1, IEND

    DO JLON = KDLUX,KDLON
      DO JLAT = KDGUN,KDGL
        ZPRAC(JLON,JLAT) = PWORK(JLON,JFL,JLAT)
      ENDDO
    ENDDO

    DO JLON = KDLUX,KDLON
      ZPRAC(JLON,KDGUN-1) = PWORK(JLON,JFL,KDGL)
      ZPRAC(JLON,KDGL +1) = PWORK(JLON,JFL,KDGUN)
    ENDDO
    DO JLAT = KDGUN,KDGL
      ZPRAC(KDLON+1,JLAT) = PWORK(KDLUN,JFL,JLAT)
    ENDDO
    ZPRAC(KDLON+1,KDGUN-1) = PWORK(KDLUN,JFL,KDGL)
    ZPRAC(KDLON+1,KDGL +1) = PWORK(KDLUN,JFL,KDGUN)

    DO JLON = KDLUX + JLL,IENX1 - JLL + 1
      DO JLAT = KDGUN, IENX2  
        PWORK(JLON,JFL,JLAT)=(4._JPRB*ZPRAC(JLON,JLAT)+2.0_JPRB*(ZPRAC(JLON+&
         & 1,JLAT)+&
         & ZPRAC(JLON-1,JLAT) + ZPRAC(JLON,JLAT+1) +&
         & ZPRAC(JLON,JLAT-1)) + ZPRAC(JLON+1,JLAT+1) +&
         & ZPRAC(JLON-1,JLAT+1) + ZPRAC(JLON+1,JLAT-1)+&
         & ZPRAC(JLON-1,JLAT-1))/16._JPRB  
      ENDDO
    ENDDO

    DO JLAT = KDGUX,KDGL
      DO JLON = KDLUN,KDLON
        ZPRAC(JLON,JLAT) = PWORK(JLON,JFL,JLAT)
      ENDDO
    ENDDO

    DO JLAT = KDGUX,KDGL
      ZPRAC(KDLUN-1,JLAT) = PWORK(KDLON,JFL,JLAT)
      ZPRAC(KDLON+1,JLAT) = PWORK(KDLUN,JFL,JLAT)
    ENDDO
    DO JLON = KDLUN,KDLON
      ZPRAC(JLON,KDGL +1) = PWORK(JLON,JFL,KDGUN)
    ENDDO
    ZPRAC(KDLUN-1,KDGL +1) = PWORK(KDLON,JFL,KDGUN)
    ZPRAC(KDLON+1,KDGL +1) = PWORK(KDLUN,JFL,KDGUN)

    DO JLAT = KDGUX + JLL, IENY1 - JLL + 1
      DO JLON = KDLUN,IENY2
        PWORK(JLON,JFL,JLAT)=(4._JPRB*ZPRAC(JLON,JLAT)+2.0_JPRB*(ZPRAC(JLON+&
         & 1,JLAT)+&
         & ZPRAC(JLON-1,JLAT) + ZPRAC(JLON,JLAT+1) +&
         & ZPRAC(JLON,JLAT-1)) + ZPRAC(JLON+1,JLAT+1) +&
         & ZPRAC(JLON-1,JLAT+1) + ZPRAC(JLON+1,JLAT-1)+&
         & ZPRAC(JLON-1,JLAT-1))/16._JPRB  
      ENDDO
    ENDDO

  ENDDO

ENDDO

!      --------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ESMOOTHE',1,ZHOOK_HANDLE)
END SUBROUTINE ESMOOTHE
END MODULE ESMOOTHE_MOD
