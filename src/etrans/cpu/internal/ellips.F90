! Jan-2011 P. Marguinaud Interface to thread-safe FA
SUBROUTINE ELLIPS (KSMAX,KMSMAX,KNTMP,KMTMP)
USE PARKIND1, ONLY : JPRD, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
!
! ***ELLIPS*** - General routine for computing elliptic truncation
!
!    Purpose.
!    --------
!       Computation of zonal and meridional limit wavenumbers within the ellipse
!    Interface:
!    ----------
!                   *CALL* *ELLIPS *
!
!        Explicit arguments :
!        --------------------
!
!        Implicit arguments :
!        --------------------
!
!
!     Method.
!     -------
!        See documentation
!
!     Externals.   NONE.
!     ----------
!
!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation
!
!     Author.
!     -------
!        G. Radnoti LACE 97/04/04
!
!     Modifications.
!
!-------------------------------------------------------------
!        J.Vivoda, 99/05/19  treating NSMAX=0 and NMSMAX=0
!        O.Nuissier, 23/09/01 Change type of real (simple --> 
!        double precision)        
!
!
INTEGER (KIND=JPIM) KSMAX, KMSMAX
INTEGER (KIND=JPIM) KNTMP(0:KMSMAX),KMTMP(0:KSMAX)
!
INTEGER (KIND=JPIM) JM, JN
!
REAL (KIND=JPRD) ZEPS, ZKN, ZKM, ZAUXIL
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ELLIPS',0,ZHOOK_HANDLE)
ZEPS=1.E-10
ZAUXIL=0.
!
! 1. Computing meridional limit wavenumbers along zonal wavenumbers
!
DO JM=1,KMSMAX-1
ZKN = REAL(KSMAX,JPRD)/REAL(KMSMAX,JPRD)* &
& SQRT(MAX(ZAUXIL,REAL(KMSMAX**2-JM**2,JPRD)))
  KNTMP(JM)=INT(ZKN+ZEPS, JPIM)
ENDDO

IF( KMSMAX.EQ.0 )THEN
   KNTMP(0)=KSMAX
ELSE
   KNTMP(0)=KSMAX
   KNTMP(KMSMAX)=0
ENDIF
!
! 2. Computing zonal limit wavenumbers along meridional wavenumbers
!             
DO JN=1,KSMAX-1
ZKM = REAL(KMSMAX,JPRD)/REAL(KSMAX,JPRD)* &
     & SQRT(MAX(ZAUXIL,REAL(KSMAX**2-JN**2,JPRD)))
  KMTMP(JN)=INT(ZKM+ZEPS, JPIM)
ENDDO   

IF( KSMAX.EQ.0 )THEN
   KMTMP(0)=KMSMAX
ELSE
   KMTMP(0)=KMSMAX
   KMTMP(KSMAX)=0
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ELLIPS',1,ZHOOK_HANDLE)
END SUBROUTINE ELLIPS
