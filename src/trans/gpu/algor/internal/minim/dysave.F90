SUBROUTINE DYSAVE (K_N,YD_Y,YD_S,P_YS,K_M,K_NYS,K_JMIN,K_JMAX,YD_YBAR,YD_SBAR,K_SELECT,&
 & K_IITER,P_OL,K_JOL,P_EPS,P_SIZE,K_MODE,K_PLEV,K_IO)  

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND2  ,ONLY : JPRH

USE CONTROL_VECTORS_MOD

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: K_M 
INTEGER(KIND=JPIM)               :: K_N ! Argument NOT used
TYPE(CONTROL_VECTOR),INTENT(IN)    :: YD_Y 
TYPE(CONTROL_VECTOR),INTENT(IN)    :: YD_S 
REAL(KIND=JPRH)   ,INTENT(IN)    :: P_YS 
INTEGER(KIND=JPIM),INTENT(INOUT) :: K_NYS 
INTEGER(KIND=JPIM),INTENT(INOUT) :: K_JMIN 
INTEGER(KIND=JPIM),INTENT(INOUT) :: K_JMAX 
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_YBAR(K_M) 
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_SBAR(K_M) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_SELECT 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IITER 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: P_OL(K_M) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: K_JOL(K_M) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_EPS 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SIZE 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_MODE 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_PLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IO 
!----

!     Save a pair (y,s)/sqrt(ys) in YBAR and SBAR and update the
!     Oren-Luenberger factor OLFACT.

!     Input:

!       ys: Euclidean scalar product (y,s). Must be positive.
!       nys: number of (y,s) pairs having been received. Must be
!         initialized to 0 before calling dysave for the first time.
!       select monitors the selection of the pairs (y,s) to build the
!         l-BFGS preconditioning matrix, with the following meanings:
!         =0: no particular selection, FIFO policy
!         =1: mexican selection (the pairs are distributed uniformely
!             on the iteration counter iiter for a particular CG run,
!             at each new run the previous pairs are discarded and a
!             new l-bfgs matrix is build from scratch)
!         =2: selection by the Rayleigh quotient; the aim is to obtain
!             uniformely distributed OL factors in the logarithmic
!             scale
!       iiter: integer giving the index of the current inner iteration
!         when dysave is called.
!       jol(m): jol(i) is the index of the pair (y,s) with the i-th
!         smallest OL factor.
!       eps = epsilon for detecting negative curvature directions,
!         whose square root is used for selecting good pairs.
!       size: scalar preconditioner.

!     Output:

!       mode: output mode
!         =-1: the (y,s) pair is not selected
!         =0 : good call,
!         =1 : the number m of pairs that can be stored in (ybar,sbar)
!              is <= 0.
!         =2 : ys <= 0.
!         =3 : unknown value of select

!----

! --- local variables

INTEGER(KIND=JPIM) :: I,I2,J,JCOUR,I_ODIFF
REAL(KIND=JPRB) :: Z_RMIN,Z_R,Z_OLMED,Z_OLDQ,Z_NEWQ,Z_OL0,Z_OL1,Z_OL2,Z_OL3,Z_OLX,Z_OLVAR,Z_OLVAR1,Z_OLVAR2,Z_EXPO
REAL(KIND=JPRH) :: Z_OLF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "dpseuclid.h"

! --- check the input

IF (LHOOK) CALL DR_HOOK('DYSAVE',0,ZHOOK_HANDLE)
IF (K_M <= 0) THEN
  K_MODE=1
  IF (K_PLEV >= 1) WRITE (K_IO,900) K_M
  IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
  RETURN
ENDIF
900 FORMAT (/" dysave-ERROR: non positive number of updates K_M = ",I5)

IF (P_YS <= 0.0_JPRB) THEN
  K_MODE=2
  IF (K_PLEV >= 1) WRITE (K_IO,901) P_YS
  901 FORMAT (/" dysave-ERROR: non positive scalar product",&
   & " (YD_Y,YD_S) = ",1PD12.5)  
  IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! --- initialization

K_MODE=0
Z_RMIN=1.E-20_JPRB

! --- Oren-Luenberger factor

!     call dpseuclid (n,y,y,olf)
!     olf=ys/olf

!     for a quadratic function, I prefer

CALL DPSEUCLID (K_N,YD_S,YD_S,Z_OLF)
Z_OLF=Z_OLF/P_YS
!>          print *,"olf=",olf
IF (K_PLEV >= 5) WRITE (K_IO,902) Z_OLF
902 FORMAT (/4X,"dysave: inverse Oren-Luenberger factor: ",&
 & "|YD_S|^2/(YD_Y,YD_S) = ",1PE10.3)  

! --- select the pair to discard if any and update the pointers

!     --- is this a good pair to save ?

IF (K_SELECT >= 0.AND. Z_OLF < SQRT(P_EPS)) THEN
  K_MODE=-1
  IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
  RETURN
ENDIF

!     --- accept any pair, discard the oldest one

IF (K_SELECT == 0) THEN
  IF (K_NYS == 0) THEN
    K_JMIN=1
    K_JMAX=0
  ENDIF
  K_NYS=K_NYS+1
  K_JMAX=K_JMAX+1
  IF (K_JMAX > K_M) K_JMAX=K_JMAX-K_M
  IF (K_NYS > K_M) THEN
    K_JMIN=K_JMIN+1
    IF (K_JMIN > K_M) K_JMIN=K_JMIN-K_M
  ENDIF
  JCOUR=K_JMAX
  P_OL(JCOUR)=Z_OLF

!     --- mexican selection

ELSEIF (K_SELECT == 1) THEN
  IF (K_IITER == 1) K_NYS=0
  IF (K_NYS == 0) THEN
    K_JMIN=1
    K_JMAX=0
    DO I=1,K_M
      K_JOL(I)=0
    ENDDO
  ENDIF
  K_NYS=K_NYS+1
  IF (K_NYS <= K_M) THEN
    K_JMAX=K_JMAX+1
    K_JOL(K_JMAX)=K_IITER
    JCOUR=K_JMAX
  ELSEIF (K_M == 2) THEN
    K_JOL(2)=K_IITER
    JCOUR=2
  ELSEIF (K_M > 2) THEN

!         --- reject current pair if iiter is too close to jol(m-1)
!             jol gives the iteration number of the selected pairs
!             they are supposed to be in order of iteration
!             here jmin=1 and jmax=m (always)

    I_ODIFF=K_JOL(K_M)-K_JOL(K_M-1)
    IF (K_IITER-K_JOL(K_M)  <  I_ODIFF) THEN
      K_MODE=-1
      IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
      RETURN
    ENDIF

!         --- select the pair jcour to discard

    DO J=1,K_M-1
      JCOUR=J+1
      IF (K_JOL(JCOUR)-K_JOL(J)  ==  I_ODIFF) GOTO 10
    ENDDO
    10 CONTINUE

!         --- shift the other (ybar,sbar) and jol
!>>>>>>>>>>>> playing with pointers instead of with the actual pairs
!             would yield a saving in time (TO DO if the selection is
!             efficient)

    DO J=JCOUR+1,K_M
      K_JOL(J-1)=K_JOL(J)
      YD_YBAR(J-1) = YD_YBAR(J)
      YD_SBAR(J-1) = YD_SBAR(J)
    ENDDO
    K_JOL(K_M)=K_IITER
    JCOUR=K_M
  ENDIF
  P_OL(JCOUR)=Z_OLF
!         print *,"jol=",(jol(i),i=1,m)

!     --- selection by the Rayleigh quotient

ELSEIF (K_SELECT == 2) THEN

!>>>>>>>> the following should be changed in the futur so that a
!         matrix build in a previous CG run could be updated (TO DO)

  IF (K_IITER == 1) K_NYS=0
  IF (K_NYS == 0) THEN
    K_JMIN=1
    K_JMAX=0
    DO I=1,K_M
      K_JOL(I)=0
    ENDDO
  ENDIF

!         --- it remains to set nys, jmax, jcour, ol(), and jol()

  IF (K_M == 1) THEN
    K_NYS=1
    K_JMAX=1
    JCOUR=1
    P_OL(1)=Z_OLF
    K_JOL(1)=1
  ELSEIF (K_NYS < K_M) THEN

!             --- select the current pair, since there is enough space
!                 to save it

    K_JMAX=K_JMAX+1
    JCOUR=K_JMAX
    P_OL(JCOUR)=Z_OLF

!             --- it remains to set jol(); for this, sort the saved
!                 pairs according to their OL factor

    IF (K_NYS == 0) THEN
      K_JOL(1)=1
    ELSE

!                 --- find the smallest index i2 such that
!                     ol(jol(i2)) >= olf, this is the position of the
!                     current OLF

      I2=1
      101 CONTINUE
      IF (P_OL(K_JOL(I2)) >= Z_OLF) GOTO 102
      I2=I2+1
      IF (I2 <= K_NYS) GOTO 101
      102 CONTINUE

!                 --- shift jol(i), for i>=i2

      IF (I2 <= K_NYS) THEN
        DO I=K_NYS,I2,-1
          K_JOL(I+1)=K_JOL(I)
        ENDDO
      ENDIF
      K_JOL(I2)=JCOUR
    ENDIF
    K_NYS=K_NYS+1

  ELSE

!             --- here, if the current (y,s) is selected, an older one
!                 must be discarded;

!                 here, nys>=m (updated at the end of the ELSE), jmax=m
!                 (no longer udated); it remains to set jcour (inside
!                 the long if-then-else below and ol(jcour)=olf (at the
!                 end of the ELSE);

!                 start by finding the smallest index i2 such that
!                 ol(jol(i2)) >= olf, this is the position of the
!                 current OLF

    I2=1
    103 CONTINUE
    IF (P_OL(K_JOL(I2)) >= Z_OLF) GOTO 104
    I2=I2+1
    IF (I2 <= K_M) GOTO 103
    104 CONTINUE
    IF (I2 == 1) THEN

!                 --- olf is smaller than all the other OL factors
!                     ==> save the (y,s) pair

      JCOUR=K_JOL(1)
    ELSEIF (I2 > K_M) THEN

!                 --- olf is greater than all the other OL factors
!                     ==> save the (y,s) pair

      JCOUR=K_JOL(K_M)
    ELSEIF (K_M == 2) THEN

!                 --- reject the (y,s) pair, since when m=2, we want to
!                     save the pairs with the extreme OL factors

      K_MODE=-1
      IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
      RETURN
    ELSEIF (I2 == 2) THEN

!                 --- check whether to replace the pair jol(2), by
!                     looking at jol(1), jol(2), and jol(3)

      Z_OLMED=SQRT(P_OL(K_JOL(1))*P_OL(K_JOL(3)))
      Z_OLDQ=P_OL(K_JOL(2))/Z_OLMED
      IF (Z_OLDQ < 1.E0_JPRB) Z_OLDQ=1.E0_JPRB/Z_OLDQ
      Z_NEWQ=Z_OLF/Z_OLMED
      IF (Z_NEWQ < 1.E0_JPRB) Z_NEWQ=1.E0_JPRB/Z_NEWQ
      IF (Z_NEWQ >= Z_OLDQ) THEN
        K_MODE=-1
        IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      JCOUR=K_JOL(2)
    ELSEIF (I2 == K_M) THEN

!                 --- check whether to replace the pair jol(m-1), by
!                     looking at jol(m-2), jol(m-1), and jol(m)

      Z_OLMED=SQRT(P_OL(K_JOL(K_M-2))*P_OL(K_JOL(K_M)))
      Z_OLDQ=P_OL(K_JOL(K_M-1))/Z_OLMED
      IF (Z_OLDQ < 1.E0_JPRB) Z_OLDQ=1.E0_JPRB/Z_OLDQ
      Z_NEWQ=Z_OLF/Z_OLMED
      IF (Z_NEWQ < 1.E0_JPRB) Z_NEWQ=1.E0_JPRB/Z_NEWQ
      IF (Z_NEWQ >= Z_OLDQ) THEN
        K_MODE=-1
        IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      JCOUR=K_JOL(K_M-1)
    ELSE

!                 --- here m>=4;
!                     check whether to replace the pair jol(i2-1) or
!                     jol(i2), by looking at jol(i2-2), jol(i2-1),
!                     jol(i2), and jol(i2+1)

      Z_OL0=LOG(P_OL(K_JOL(I2-2)))
      Z_OL1=LOG(P_OL(K_JOL(I2-1)))
      Z_OL2=LOG(P_OL(K_JOL(I2)))
      Z_OL3=LOG(P_OL(K_JOL(I2+1)))
      Z_OLX=LOG(Z_OLF)
      Z_OLVAR =MAX(Z_OL1-Z_OL0,Z_OL2-Z_OL1,Z_OL3-Z_OL2)
      Z_OLVAR1=MAX(Z_OLX-Z_OL0,Z_OL2-Z_OLX,Z_OL3-Z_OL2)
      Z_OLVAR2=MAX(Z_OL1-Z_OL0,Z_OLX-Z_OL1,Z_OL3-Z_OLX)
      IF (Z_OLVAR <= MIN(Z_OLVAR1,Z_OLVAR2)) THEN
        K_MODE=-1
        IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF (Z_OLVAR1 < Z_OLVAR2) THEN
        JCOUR=K_JOL(I2-1)
      ELSE
        JCOUR=K_JOL(I2)
      ENDIF
    ENDIF

!             --- save the pair

    K_NYS=K_NYS+1
    P_OL(JCOUR)=Z_OLF
  ENDIF
!>          print *," "
!>          print *,(ol(i),i=1,min(nys,m))
!>          print *,"jol=",(jol(i),i=1,min(nys,m))
!>          write(6,999) (ol(jol(i)),i=1,min(nys,m))
!>  999     format (1pe12.6,2x,50(e12.6,2x))
ELSE

!     --- unknown value of select

  K_MODE=3
  IF (K_PLEV >= 1) WRITE (K_IO,903) K_SELECT
  903 FORMAT (/" dysave-ERROR: unknown selection procedure,"" K_SELECT = ",I5)
  IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! --- store ybar et sbar

Z_R=SQRT(1.E0_JPRB/P_YS)
YD_YBAR(JCOUR)%DATA = Z_R*YD_Y%DATA
YD_SBAR(JCOUR)%DATA = Z_R*YD_S%DATA

! --- save the OL factor

!     --- the simplest strategy is to take the OL factor at the
!         first CG iteration (the function is supposed to be
!         quadratic)

IF (K_NYS == 1) P_SIZE=Z_OLF

!     --- another strategy is to take the geometric means of the
!         previous OL factors (taking the d'OL factor, the closest
!         in the logarithmic scale to this geometric means give
!         approximately the same results)

Z_OLMED=1.E0_JPRB
Z_EXPO=1.E0_JPRB/REAL(MIN(K_NYS,K_M),JPRB)
DO J=1,MIN(K_NYS,K_M)
  Z_OLMED=Z_OLMED*(P_OL(J)**Z_EXPO)
ENDDO
P_SIZE = Z_OLMED
IF (LHOOK) CALL DR_HOOK('DYSAVE',1,ZHOOK_HANDLE)
!      print *,"ol=",(ol(j),j=1,min(nys,m))
!      print *,"olmed=",olmed

END SUBROUTINE DYSAVE
