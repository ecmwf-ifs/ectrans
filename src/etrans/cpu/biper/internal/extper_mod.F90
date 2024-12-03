! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EXTPER_MOD
CONTAINS
SUBROUTINE EXTPER(PWORK,KDIM,KPSTA,KPOINTS,KFLDS,KUNITS,&
  & KPOINTERS,KALFA) 

!   purpose  :
!   --------
!     Make spline extension.

!    *CALL* *EXTPER(PWORK,KDIM,KPSTA,KPOINTS,KFLDS,KUNITS,&
!         & KPOINTERS,KALFA)

!      externals :
!      ----------
!             None

!      explicit arguments :
!      ------------------
!     PWORK   : Input:  values in C U I area
!             : Output: input+(spline extension on the E area)
!     KDIM    : Dimension of the C U I U E unit of work (one row or one m)
!     KPSTA   : Position where the unit of work starts
!     KPOINTS : Position where the unit of work ends       
!     KFLDS   : number of 2D fields
!     KUNITS  : Number of units of work
!     KPOINTERS : Array of pointers for the units of work
!     KALFA : boundary condition of a spline:
!             = 0 ... natural spline
!             = 1 ... boundary condition computed differentially
!                      (additional option)
!      references :
!      ----------

!      author :
!      ------
!         M. Hortal 03-11-2009
!      -----------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN
USE TPM_DISTR

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWORK(:,:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPOINTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KUNITS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPOINTERS(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KALFA 

!      arrays :
!     --------
INTEGER(KIND=JPIM) :: IENDX, IENDY, JFL, JLAT, JLON, IA

REAL(KIND=JPRB) :: ZA, ZB, ZC, ZD, ZEPSA, ZEPSB, ZJ, ZK, ZKP1,&
 & ZLAM, ZLAMB, ZM1, ZM2, ZMM, ZNY  
REAL(KIND=JPRB) :: ZMAX(KUNITS), ZMIN(KUNITS)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
 
!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EXTPER',0,ZHOOK_HANDLE)

!*         0. Security
!             --------

IF(UBOUND(PWORK,1) < KFLDS) THEN
  CALL ABOR1(' EXTPER, PWORK first dimension too small')
ENDIF
IF(UBOUND(PWORK,2) < KDIM+2) THEN
  WRITE(NOUT,*) ' UBOUND(PWORK,2)=',UBOUND(PWORK,2),' KDIM=',KDIM,' KUNITS=',&
    &KUNITS
  CALL ABOR1(' EXTPER, PWORK second dimension too small')
ENDIF
IF(UBOUND(KPOINTERS,1) < KUNITS) THEN
  CALL ABOR1(' EXTPER, KPOINTERS too small')
ENDIF
IF(UBOUND(PWORK,2) < KPOINTERS(KUNITS)+KDIM) THEN
  WRITE(NERR,*) ' EXTPER, KUNITS=',KUNITS,' KPOINTERS=',KPOINTERS(1:KUNITS),&
    &' KDIM=',KDIM,' UBOUND(PWORK,2)=',UBOUND(PWORK,2)
  CALL ABOR1(' EXTPER, value of KPOINTERS too large')
ENDIF

!*         1. Spline   Extension.
!             -------------------

DO JFL = 1, KFLDS

  ZK    = REAL(KDIM-KPOINTS+1,JPRB)
  ZKP1  = ZK + 1.0_JPRB
  ZLAMB = ZK/ZKP1
  ZNY   = REAL(KALFA,JPRB)/ZKP1

  DO JLAT=1,KUNITS
    ZEPSA = &
     &((PWORK(JFL,KPOINTERS(JLAT)+KPSTA)-&
     & PWORK(JFL,KPOINTERS(JLAT)+KPOINTS))/ZK -&
     & PWORK(JFL,KPOINTERS(JLAT)+KPOINTS)+&
     & PWORK(JFL,KPOINTERS(JLAT)+KPOINTS-1))*6._JPRB/ZKP1 -&
     & ZNY*(PWORK(JFL,KPOINTERS(JLAT)+KPOINTS)-&
     & 2.0_JPRB* PWORK(JFL,KPOINTERS(JLAT)+KPOINTS-1)+&
     & PWORK(JFL,KPOINTERS(JLAT)+KPOINTS-2)) 

    ZEPSB = (PWORK(JFL,KPOINTERS(JLAT)+KPSTA+1)-&
     & PWORK(JFL,KPOINTERS(JLAT)+KPSTA) -&
     & (PWORK(JFL,KPOINTERS(JLAT)+KPSTA)-&
     & PWORK(JFL,KPOINTERS(JLAT)+KPOINTS))/ZK)*6._JPRB/ZKP1-&
     & ZNY*(PWORK(JFL,KPOINTERS(JLAT)+KPSTA+2)-&
     & 2.0_JPRB* PWORK(JFL,KPOINTERS(JLAT)+KPSTA+1)+&
     & PWORK(JFL,KPOINTERS(JLAT)+KPSTA))  

    ZMM = 4._JPRB - ZLAMB*ZLAMB
    ZM1 = (2.0_JPRB*ZEPSA - ZLAMB*ZEPSB)/ZMM
    ZM2 = (2.0_JPRB*ZEPSB - ZLAMB*ZEPSA)/ZMM
    ZA  = PWORK(JFL,KPOINTERS(JLAT)+KPOINTS)
    ZB  = (PWORK(JFL,KPOINTERS(JLAT)+KPSTA)-&
        & PWORK(JFL,KPOINTERS(JLAT)+KPOINTS))/ZK-&
     & (2.0_JPRB*ZM1 + ZM2) * ZK/6._JPRB  
    ZC  = 0.5_JPRB * ZM1
    ZD  = (ZM2 - ZM1)/(6._JPRB*ZK)


    DO JLON=KPOINTERS(JLAT)+KPOINTS+1,KPOINTERS(JLAT)+KDIM

      ZJ  = REAL(JLON - (KPOINTERS(JLAT)+KPOINTS),JPRB)
      PWORK(JFL,JLON) = ZA + ZJ * (ZB + ZJ * (ZC + ZD * ZJ))
    ENDDO 
  ENDDO 


ENDDO

IF (LHOOK) CALL DR_HOOK('EXTPER',1,ZHOOK_HANDLE)
END SUBROUTINE EXTPER
END MODULE EXTPER_MOD
