! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


SUBROUTINE ETIBIHIE(KDLON,KDGL,KNUBI,KDLUX,KDGUX,&
 & KSTART,KDLSM,PGPBI,LDBIX,LDBIY,KDADD)  

!****   tool  ETIBIHIE : Doubly-periodicisation : isotropic spline
!       -------------   method.

!   purpose  :
!   --------
!     KNUBI  horizontal fields which are known on C U I,
!     are extended over E, in order to obtain  doubly-periodic
!     fields.
!     IF LDBIX is equal .TRUE. , then the fields are periodicise
!     in the x ( or longitude ) direction.  If it is not the case,
!     KDLUX must be equal to KDLON.
!     IF LDBIY is equal .TRUE. , then the fields are periodicise
!     in the y ( or latitude ) direction.   If it is not the case,
!     KDGUX must be equal to KDGL.

!*     *CALL* *ETIBIHIE*(...)

!      externals :
!      ----------
!      ESPLIN  spline extension
!      ESMOOTH smoothing across to get isotropy.

!      explicit arguments :
!      ------------------
!     KDLON : upper bound for the x (or longitude) dimension
!             of the gridpoint array on C U I U E
!     KDGL : upper bound for the y (or latitude) dimension
!             of the gridpoint array on C U I U E
!     KNUBI : number of horizontal fields to doubly-periodicise.
!     KDLUX : upper bound for the x (or longitude) dimension
!             of  C U I.
!     KDGUX : upper bound for the y (or latitude) dimension
!             of  C U I.
!     KSTART : first  dimension in x direction of g-p array
!     KDLSM  : second dimension in x direction of g-p array
!     PGPBI : gridpoint array on C U I U E.
!     LDBIX  : logical to periodicize or not
!             in the x ( or longitude ) direction.
!     LDBIY  : logical to periodicize  or not
!             in the y ( or latitude ) direction.
!     KDADD : 1 to test biperiodiz.

!      references :
!      ----------

!      author :
!      ------
!          V. Ducrocq

!      modification :
!      ------------
!          A. Stanesic  28/03/2008: KDADD - test of externalized biper.
! -------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE ESPLINE_MOD
USE ESMOOTHE_MOD

! -------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)               :: KNUBI
INTEGER(KIND=JPIM),INTENT(IN)               :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)               :: KDLSM 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDLUX 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDGUX 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDADD
REAL(KIND=JPRB),INTENT(INOUT)               :: PGPBI(KSTART:KDLSM+KDADD,KNUBI,1:KDGL+KDADD) 
LOGICAL,INTENT(IN)                          :: LDBIX 
LOGICAL,INTENT(IN)                          :: LDBIY 

! -------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZALFA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ETIBIHIE',0,ZHOOK_HANDLE)
! -------------------------------------------------------------------------

!*         1. DOUBLY-PERIODICISE :
!             ------------------

ZALFA = 0.0_JPRB

CALL ESPLINE(1,KDLON,1,KDGL,KDLUX,KDGUX,KSTART,&
 & KDLSM+KDADD,1,KDGL+KDADD,KNUBI,PGPBI,ZALFA,LDBIX,LDBIY,KDADD)
CALL ESMOOTHE(1,KDLON,1,KDGL,KDLUX,KDGUX,KSTART,&
 & KDLSM+KDADD,1,KDGL+KDADD,KNUBI,PGPBI,LDBIX,LDBIY)  

! -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ETIBIHIE',1,ZHOOK_HANDLE)
END SUBROUTINE ETIBIHIE
