! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


SUBROUTINE FPBIPERE(KDLUX,KDGUX,KDLON,KDGL,KNUBI,KD1,PGPBI,KDADD,LDZON, &
 & LDBOYD, KDBOYD, PLBOYD, LDSEAM, KDSEAM, PLSEAM)

!****   *FPBIPERE*  - Full-POS interface for double periodicisation

!   purpose  :
!   --------
!     To bi-periodicise the post-processed fields, or just fill the extension zone
!      with the mean value of C+I area

!**   INTERFACE.
!     ----------
!         *CALL*   *FPBIPERE*(...)

!        EXPLICIT ARGUMENTS
!        --------------------
!         KDLUX : upper bound for the x (or longitude) dimension of  C U I.
!         KDGUX : upper bound for the y (or latitude) dimension of  C U I.
!         KDLON : upper bound for the x (or longitude) dimension of the gridpoint array on C U I U E
!         KDGL  : upper bound for the y (or latitude) dimension of the gridpoint array on C U I U E
!         KNUBI : number of horizontal fields to doubly-periodicise.
!         KD1   : dimension of input/output array
!         PGPBI : input/output gridpoint array on C U I U E.
!         LDZON : .true. if input grid on C U I U E (.false. if C U I)
!         KDADD : 1 to test biperiodiz.
!         LDBOYD: perform boyd periodization (inside C U I)
!         KDBOYD: array containing dimensions of boyd domain
!         PLBOYD: scalar parameter for boyd (variable L in paper)
!         LDSEAM: perform SEAM periodization
!         KDSEAM: size of the SEAM mixing zone
!         PLSEAM: scalar parameter for SEAM (Lmix and Lboyd)

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        ESMOOTHE smoothing across to get isotropy.
!        ESPLINE  spline extension
!        ESEAM    Symmetric Extension And Mixing method

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!      F. Taillefer : 04-10-21 Add LDZON
!      A. Stanesic  : 28-03-08: KDADD - test of externalized biper.
!      D. Degrauwe  : feb 2012 Boyd periodization
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2
!      R. El Khatib 04-Aug-2016 new interface to ewindowe + cleaning
!      B. Menetrier : 02-06-2025 new SEAM method
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK,  JPHOOK
USE ESEAM_MOD
USE ESPLINE_MOD
USE ESMOOTHE_MOD
USE EWINDOWE_MOD
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KNUBI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KD1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPBI(KD1,KNUBI) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDADD
LOGICAL, OPTIONAL ,INTENT(IN)    :: LDZON
LOGICAL           ,INTENT(IN) ,OPTIONAL :: LDBOYD
INTEGER(KIND=JPIM),INTENT(IN) ,OPTIONAL :: KDBOYD(6)
REAL(KIND=JPRB)   ,INTENT(IN) ,OPTIONAL :: PLBOYD
LOGICAL           ,INTENT(IN) ,OPTIONAL :: LDSEAM
INTEGER(KIND=JPIM),INTENT(IN) ,OPTIONAL :: KDSEAM
REAL(KIND=JPRB)   ,INTENT(IN) ,OPTIONAL :: PLSEAM(2)

!     ------------------------------------------------------------------

REAL(KIND=JPRB), ALLOCATABLE :: ZGPBI(:,:,:)
INTEGER(KIND=JPIM) :: IND, ISTAE, JGL, JLON, JNUBI, IBWX, IBWY
LOGICAL         :: LLZON, LLBOYD, LLSEAM
REAL(KIND=JPRB) :: ZALFA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPBIPERE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLBOYD=.FALSE.
IF (PRESENT(LDBOYD)) LLBOYD=LDBOYD
LLSEAM=.FALSE.
IF (PRESENT(LDSEAM)) LLSEAM=LDSEAM

!*         2. DOUBLY-PERIODICISE
!             ------------------

IF (LLBOYD) THEN
  IF (.NOT.PRESENT(KDBOYD)) CALL ABOR1('FPBIPERE: Boyd periodization requires KDBOYD argument')
  IF (.NOT.PRESENT(PLBOYD)) CALL ABOR1('FPBIPERE: Boyd periodization requires PLBOYD argument')
  IBWX=KDBOYD(3)
  IBWY=KDBOYD(6)
  CALL EWINDOWE(KDLON,KDLUX,IBWX,KDGL,KDGUX,IBWY,KNUBI,PGPBI,PLBOYD,.TRUE.,.TRUE.)
ELSE
  LLZON=.FALSE.
  IF(PRESENT(LDZON)) LLZON=LDZON
  ALLOCATE(ZGPBI(KDLON+KDADD,KNUBI,KDGL+KDADD))
  IF(LLZON) THEN
!   Copy C+I+E
    IND=KDLON
  ELSE
!   Copy C+I
    IND=KDLUX
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JNUBI,ISTAE,JGL,JLON)
  DO JNUBI=1,KNUBI
    ISTAE=0
    DO JGL=1,KDGUX
      DO JLON=1,KDLUX
        ZGPBI(JLON,JNUBI,JGL)=PGPBI(ISTAE+JLON,JNUBI)
      ENDDO
      ISTAE=ISTAE+IND
    ENDDO
  ENDDO 
!$OMP END PARALLEL DO
  IF (LLSEAM) THEN
    IF (.NOT.PRESENT(KDSEAM)) CALL ABOR1('FPBIPERE: SEAM periodization requires KDSEAM argument')
    IF (.NOT.PRESENT(PLSEAM)) CALL ABOR1('FPBIPERE: SEAM periodization requires PLSEAM argument')
    IF (KDADD>0) CALL ABOR1('FPBIPERE: SEAM periodization requires KDADD=0')
    CALL ESEAM(1,KDLON,1,KDGL,KDLUX,KDGUX,KNUBI,ZGPBI,KDSEAM,PLSEAM,.TRUE.,.TRUE.)
  ELSE
    ZALFA = 0.0_JPRB
    CALL ESPLINE(1,KDLON,1,KDGL,KDLUX,KDGUX,1,KDLON+KDADD,1,KDGL+KDADD,KNUBI,ZGPBI,&
     & ZALFA,.TRUE.,.TRUE.,KDADD)
    CALL ESMOOTHE(1,KDLON,1,KDGL,KDLUX,KDGUX,1,KDLON+KDADD,1,KDGL+KDADD,KNUBI,ZGPBI,&
     & .TRUE.,.TRUE.)
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JNUBI,ISTAE,JGL,JLON)
  DO JNUBI=1,KNUBI 
    ISTAE=0
    DO JGL=1,KDGL
      DO JLON=1,KDLON
        PGPBI(ISTAE+JLON,JNUBI)=ZGPBI(JLON,JNUBI,JGL)
      ENDDO
      ISTAE=ISTAE+KDLON
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  DEALLOCATE(ZGPBI)
ENDIF


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPBIPERE',1,ZHOOK_HANDLE)
END SUBROUTINE FPBIPERE
