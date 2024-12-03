! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EWINDOWE_MOD

CONTAINS

SUBROUTINE EWINDOWE(KDLON,KDLUX,KBWX,KDGL,KDGUX,KBWY,KFLD,PGPIN,PSCAL,LDBIX,LDBIY)

!   purpose  :
!   --------
!     Make boyd periodic extension.

!      externals :
!      ----------
!             None

!      explicit arguments :
!      ------------------
!     KDLON  : upper bound for the x (or longitude) dimension
!              of  C U I U P.  
!     KDGL  : upper bound for the y (or latitude) dimension
!              of the gridpoint array on C U I U P
!     PGPIN  : gridpoint array on C U I U P (gp:fields).
!     PSCAL  : window function scaling  parameter
!     LDBIX  : .TRUE. windowing  in x direction  ( and vice versa )
!     LDBIY  : .TRUE. windowing  in y direction  ( and vice versa )


!      references :
!      ----------

!      author : Fabrice Voitus and Piet Termonia, 07/2009
!      ------
!      
!      modification :
!         Daan Degrauwe    02/2012    Cleaned and generalized
!         S. Martinez      03/2012    Calls to ERF under CPP key __PGI
!                                     (ERF function is not intrinsic with PGI)
!         R. El Khatib 27-Sep-2013 implicit sized PGPIN 
!         R. El Khatib 04-Aug-2016 new interface
!      -----------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KBWX
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KBWY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLD
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPIN((KDLUX+2*KBWX+2*(KDLON-KDLUX))*(KDGUX+2*KBWY+2*(KDGL-KDGUX)),KFLD)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCAL 
LOGICAL           ,INTENT(IN)    :: LDBIX 
LOGICAL           ,INTENT(IN)    :: LDBIY  

!      FERF function
!      -------------

#ifdef __PGI
REAL(KIND=JPRB), EXTERNAL :: ERF
#endif

!      scalars
!     --------

INTEGER(KIND=JPIM) :: JFL, JGL, JLON, IOFF, IDLW, IDGW
INTEGER(KIND=JPIM) :: IWX, ILWX, IRWX, IWY, ILWY, IRWY, IBWXO, IBWYO
INTEGER(KIND=JPIM) :: ILATF, ILONF, IND1, IND, IOFF_LEFT,IOFF_RIGHT,IOFF_BOTTOM,IOFF_TOP
REAL(KIND=JPRB) :: ZI, ZJ, ZK, ZL  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     local arrays :
!     ------------

REAL(KIND=JPRB)  :: ZBELX(2*KBWX+(KDLON-KDLUX))
REAL(KIND=JPRB)  :: ZBELY(2*KBWY+(KDGL -KDGUX))

!*         1. Boyd Bi-periodic Extension Method.
!             ---------------------------------

IF (LHOOK) CALL DR_HOOK('EWINDOWE',0,ZHOOK_HANDLE)

IF ((.NOT.LDBIX).AND.(.NOT.LDBIY)) THEN
  IF (LHOOK) CALL DR_HOOK('EWINDOWE',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IDGW=SIZE(ZBELY)
IDLW=SIZE(ZBELX)

!   Bell window functions :
!   ---------------------

IF (LDBIX) THEN
  DO JLON=1,IDLW
    ! variable between -1 and 1
    ZJ=REAL(-IDLW-1+2*JLON,JPRB)/(IDLW+1)
    ZL=ZJ/SQRT(1.0_JPRB-(ZJ*ZJ))
#ifdef __PGI
    ZBELX(JLON)=(1.0_JPRB+ERF(REAL(PSCAL*ZL)))/2.0_JPRB
#else
    ZBELX(JLON)=(1.0_JPRB+ERF(PSCAL*ZL))/2.0_JPRB
#endif
  ENDDO
ENDIF

IF (LDBIY) THEN
  DO JGL=1,IDGW
    ! variable between -1 and 1
    ZJ=REAL(-IDGW-1+2*JGL,JPRB)/(IDGW+1)
    ZL=ZJ/SQRT(1.0_JPRB-(ZJ*ZJ))
#ifdef __PGI
    ZBELY(JGL)=(1.0_JPRB+ERF(REAL(PSCAL*ZL)))/2.0_JPRB
#else
    ZBELY(JGL)=(1.0_JPRB+ERF(PSCAL*ZL))/2.0_JPRB
#endif
  ENDDO
ENDIF


!    Windowing on P+G-zone :
!    --------------------

IOFF=(KDLUX+2*(KBWX+KDGL-KDGUX))
IBWXO=KBWX+(KDLON-KDLUX)
IBWYO=KBWY+(KDGL-KDGUX)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFL,JGL,JLON,ILONF,ILATF,IND1,IND,IOFF_LEFT,IOFF_RIGHT,IOFF_BOTTOM,IOFF_TOP)
DO JFL=1,KFLD
  IF (LDBIX) THEN
    ! X-direction
    DO JGL=1,KDGL+IDGW
      IOFF_LEFT=(JGL-1)*IOFF
      IOFF_RIGHT=IOFF_LEFT+KDLON
      DO JLON=1,IDLW
        PGPIN(IOFF_RIGHT+JLON,JFL) = ZBELX(JLON)*PGPIN(IOFF_LEFT+JLON,JFL) +&
         & (1.0_JPRB-ZBELX(JLON))*PGPIN(IOFF_RIGHT+JLON,JFL)
      ENDDO
    ENDDO 
  ENDIF
  IF (LDBIY) THEN
    ! Y-direction  
    DO JGL=1,IDGW
      IOFF_BOTTOM=(JGL-1)*IOFF
      IOFF_TOP=(KDGL+JGL-1)*IOFF
!DIR$ IVDEP
      DO JLON=1,KDLON+IDLW
        PGPIN(IOFF_TOP+JLON,JFL) = ZBELY(JGL)*PGPIN(IOFF_BOTTOM+JLON,JFL) +&
         & (1.0_JPRB-ZBELY(JGL))*PGPIN(IOFF_TOP+JLON,JFL)
      ENDDO 
    ENDDO 
  ENDIF
ENDDO
!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('EWINDOWE',1,ZHOOK_HANDLE)

END SUBROUTINE EWINDOWE

END MODULE EWINDOWE_MOD
