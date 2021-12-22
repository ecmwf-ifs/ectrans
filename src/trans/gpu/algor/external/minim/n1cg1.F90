SUBROUTINE N1CG1 (SIMUL,K_N,YD_X,P_EPSNEG,P_EPS,K_ITER,K_IMP,K_IO,K_MODE,&
 & K_PRECO,K_M0,K_ILM0,K_NILM0,YD_YBAR0,YD_SBAR0,P_SIZE0,&
 & K_BFGSB,K_M1,K_ILM1,K_NILM1,YD_YBAR1,YD_SBAR1,P_SIZE1,&
 & K_SELECT,K_NUPTRA,P_F,YD_R,YD_YS,YD_YRS)  

!-----------------------------------------------------------------------

!     This routine solves for x the linear system (LS)

!         A x = b,

!     by a possibly preconditioned Fletcher-Reeves conjugate gradient
!     (CG) algorithm. The matrix A and the vector b are supposed given
!     The matrix A must be symmetric positive definite. The way A is
!     stored is irrelevant for N1CG1, since the algorithm requires
!     only products Av of A by various vectors v. This operation is
!     supposed to be made in the routine MVPROD.

!     The CG iterations can be preconditioned by a full positive
!     definite matrix stored in PMAT0 in the standard way (when
!     preco=1) or by a limited memory BFGS (l-BFGS) matrix stored in
!     the structure (M0, ILM0, NILM0, WLM0, NWLM0) (when preco=2). This
!     particularly appropriate when the dimension N of the LS is very
!     large. The l-BFGS structure, or even PMAT0, may be filled in by a
!     previous call to N1CG1. This routine has indeed the feature
!     (when BFGSB is nonzero) to build a preconditioning matrix using
!     pairs (Av,v). This preconditioning matrix can then be used to
!     speed-up another CG run for a system having a matrix not too
!     different from A.

!     Input
!     ^^^^^

!       mvprod, external: routine making the matrix A times vector v
!         product; it must be defined as follows:
!           subroutine mvprod (n,a,v,av,izs,rzs,dzs)
!         where, n (I), a (I), and izs (I), rzs (I), dzs (I) have the
!         same meanings as below, v(n) (I) is a double precision vector
!         and av(n) (O) is a double precision vector containing the
!         product of the matrix A times the vector V.
!       n, integer: dimension of the LS to solve
!       x(n), double precision: initial guess for x
!       b(n), double precision: RHS of the LS
!       a(), double precision: address of the matrix A, it is passed
!         to MVPROD
!       epsneg, double precision: if during CG iterations, a Rayleigh
!         quotient v'Av/v'v less than epsneg is encountered, the
!         algorithm stops, since the matrix A is not considered to be
!         positive definite
!       eps, double precision: if during CG iterations, the relative
!         norm of the residual |Ax-b|/|Ax_0-b| is less than eps, the
!         algorithm stops, since the LS is considered to be solved by
!         the current x
!       iter, integer: max number of iterations accepted
!       imp, integer: printing level
!       io, integer: output channel for printings
!       w(nw), double precision: working area
!       nw, integer: dimension of w, must be >= 5*n+m0
!       preco, integer: specifies the type of preconditioning matrix
!         made available to precondition the current CG run
!         =0: no preconditioning
!         =1: a full symmetric positive definite preconditioning matrix
!             of order n is available in the structure (pmat0,npmat0)
!         =2: a limited memory BFGS preconditioning matrix is available
!             in the l-BFGS structure (m0,ilm0,nilm0,wlm0,nwlm0), built
!             up by a previous called to n1cg1
!       (pmat0,npmat0): if preco=1, this structure must contain
!         preconditioning information (not used otherwise)
!         . pmat0(npmat0), double precision: must contain an order n
!           symmetric positive definite matrix approximating the
!           inverse of A
!         . npmat0, integer: dimension of pmat0, must be >= n**2
!       (m0,ilm0,nilm0,wlm0,nwlm0): if preco=2, must contain an l-BFGS
!         preconditioning structure (not used otherwise)
!         . m0, integer: number of (y,s) pairs stored in the structure,
!         . ilm0(nilm0), integer
!         . nilm0, integer: dimension of ilm0
!         . wlm0(nwlm0), double precision: contains the (y,s) pairs
!         . nwlm0, integer: dimension of wlm0, must be >= m0*(2*n+1)+1
!       bfgsb, integer: specifies the type of BFGS preconditioning
!         matrix to build during the CG iterations,
!         =0: don't build any preconditioning matrix
!         =1: build a full BFGS preconditioning matrix in the BFGS
!             structure (pmat1,npmat1); if pmat1(1)=0 the updates start
!             from scratch, otherwise pmat1 is updated
!         =2: build a limited memory BFGS preconditioning matrix in the
!             l-BFGS structure (m1,nys1,jmin1,jmax1,ybar1,sbar1,size1)
!>>>>>>       (FOR THE WHILE, THE UPDATES START FROM SCRATCH AT EACH NEW
!>>>>>>       CG RUN)
!       select, in the case when bfgsb=2, it monitors the selection of
!         the pairs (Av,v) to build the l-BFGS preconditioning matrix,
!         with the following meanings:
!         =0: no particular selection, FIFO policy
!>>>>>>       (TO DATE, THIS SEEMS TO BE THE BEST CHOICE)
!         =1: mexican selection (the pairs are distributed uniformely
!             according to the iteration counter for a particular CG
!             run)
!         =2: by the Rayleigh quotient (it is tried that the Rayleigh
!             quotient of the pairs be distributed as uniformely as
!             possible on the logarithmic scale)

!     Output
!     ^^^^^^

!       iter = number of iterations
!       mode:
!         =0: normal terminaison,
!         =1: bad dimensions
!         =2: A is probably indefinite,
!         =3: maximum number of iterations reached.
!       (pmat1,npmat1): if bfgsb=1, will contain the plain BFGS
!         preconditioning structure (not used otherwise)
!         . pmat1(npmat1), double precision: will contain the matrix
!         . npmat1, integer: dimension of pmat1, must be >= n**2
!       (m1,ilm1,nilm1,wlm1,nwlm1): if bfgsb=2, will contain the l-BFGS
!         preconditioning structure (not used otherwise)
!         . m1, integer: number of (y,s) pairs to store,
!         . ilm1(nilm1), integer
!         . nilm1, integer: dimension of ilm1
!         . wlm1(nwlm1), double precision: contains the (y,s) pairs
!         . nwlm1, integer: dimension of wlm1, must be >= m1*(2*n+1)+1

!     Other arguments
!     ^^^^^^^^^^^^^^^

!       izs() [integer], rzs() [real], and dzs() [double precision]:
!         these variables are not used by N1CG1, they are just passed
!         to the user-supplied routine mvprod

!     Author
!     ^^^^^^

!       Version 1.0: J.Ch. Gilbert, Inria, July 1998

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND2  ,ONLY : JPRH

USE CONTROL_VECTORS_MOD

!-----------------------------------------------------------------------

IMPLICIT NONE

EXTERNAL SIMUL
INTEGER(KIND=JPIM),INTENT(IN)       :: K_NILM0 
INTEGER(KIND=JPIM),INTENT(IN)       :: K_NILM1  
INTEGER(KIND=JPIM),INTENT(IN)       :: K_N
TYPE(CONTROL_VECTOR),INTENT(INOUT)  :: YD_X
REAL(KIND=JPRB)   ,INTENT(IN)       :: P_EPSNEG 
REAL(KIND=JPRB)   ,INTENT(IN)       :: P_EPS 
INTEGER(KIND=JPIM),INTENT(INOUT)    :: K_ITER 
INTEGER(KIND=JPIM),INTENT(IN)       :: K_IMP
INTEGER(KIND=JPIM),INTENT(IN)       :: K_IO 
INTEGER(KIND=JPIM),INTENT(OUT)      :: K_MODE
INTEGER(KIND=JPIM),INTENT(IN)       :: K_PRECO 
INTEGER(KIND=JPIM),INTENT(IN)       :: K_M0 
INTEGER(KIND=JPIM),INTENT(IN)       :: K_ILM0(K_NILM0) 
TYPE(CONTROL_VECTOR),INTENT(IN)     :: YD_YBAR0(K_M0) 
TYPE(CONTROL_VECTOR),INTENT(IN)     :: YD_SBAR0(K_M0) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: P_SIZE0
INTEGER(KIND=JPIM),INTENT(IN)       :: K_BFGSB 
INTEGER(KIND=JPIM),INTENT(INOUT)    :: K_M1 
INTEGER(KIND=JPIM),INTENT(INOUT)    :: K_ILM1(K_NILM1) 
TYPE(CONTROL_VECTOR),INTENT(INOUT)  :: YD_YBAR1(K_M1) 
TYPE(CONTROL_VECTOR),INTENT(INOUT)  :: YD_SBAR1(K_M1)
REAL(KIND=JPRB)   ,INTENT(OUT)      :: P_SIZE1 
INTEGER(KIND=JPIM),INTENT(IN)       :: K_SELECT 
INTEGER(KIND=JPIM),INTENT(IN)       :: K_NUPTRA 
REAL(KIND=JPRB)   ,INTENT(IN)       :: P_F 
TYPE(CONTROL_VECTOR),INTENT(INOUT)  :: YD_R
TYPE(CONTROL_VECTOR) ,INTENT(IN)    :: YD_YS
TYPE(CONTROL_VECTOR) ,INTENT(IN)    :: YD_YRS

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INDIC
REAL(KIND=JPRB) :: Z_OL1(K_M1)
REAL(KIND=JPRB) :: Z_EPS2,Z_F0,Z_RHO(K_M0)
TYPE(CONTROL_VECTOR) :: YL_B, YL_Q, YL_V, YL_RM, YL_NUL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "n1cga.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('N1CG1',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------

CALL ALLOCATE_CTLVEC(YL_B)
CALL ALLOCATE_CTLVEC(YL_Q)
CALL ALLOCATE_CTLVEC(YL_V)
CALL ALLOCATE_CTLVEC(YL_RM)

!---- initial printings and check the arguments

IF (K_IMP >= 1) THEN
  WRITE (K_IO,900) K_N,P_EPSNEG,P_EPS,K_ITER,K_IMP
  IF (K_PRECO == 0) THEN
    WRITE (K_IO,901)
  ELSEIF (K_PRECO == 1) THEN
    WRITE (K_IO,902)
  ELSEIF (K_PRECO == 2) THEN
    WRITE (K_IO,903) K_M0
  ENDIF
  IF (K_BFGSB == 0) THEN
    WRITE (K_IO,904)
  ELSEIF (K_BFGSB == 1) THEN
    WRITE (K_IO,905)
  ELSEIF (K_BFGSB == 2) THEN
    WRITE (K_IO,906) K_M1
  ENDIF
  IF (K_IMP >= 2) WRITE (K_IO,910)
ENDIF
900 FORMAT (/" n1cg1 (Version 1.0_JPRB, July 1998): entry point"/&
 & 4X,"dimension of the problem (K_N):",I14/&
 & 4X,"Rayleigh quotient threshold (P_EPSNEG):",4X,1PE9.2/&
 & 4X,"relative precision on Z_G (P_EPS):",11X,E9.2/&
 & 4X,"maximal number of iterations (K_ITER):",I7/&
 & 4X,"printing level (K_IMP):",I22)  
901 FORMAT (/4X,"Plain CG iterations")
902 FORMAT (/4X,"Preconditionned CG iterations")
903 FORMAT (/4X,"l-BFGS preconditioned CG iterations, with",I3," pairs")
904 FORMAT (4X,"No preconditioning matrix is built")
905 FORMAT (4X,"A full BFGS preconditioning matrix is built")
906 FORMAT (4X,"An l-BFGS preconditioning matrix is built, with",I3," pairs")
910 FORMAT (/1X,79("-"))

IF (K_N <= 0.OR. K_ITER <= 0.OR.P_EPS <= 0.E+0_JPRB .OR. P_EPS > 1.E+0_JPRB) THEN
  K_MODE=1
  IF (K_IMP >= 1) WRITE (K_IO,920)
  IF (LHOOK) CALL DR_HOOK('N1CG1',1,ZHOOK_HANDLE)
  RETURN
ENDIF
920 FORMAT (/" >>> n1cg1: inconsistent call")

!      nq=1
!      nr=nq+n
!      nv=nr+n
!      nrm=nv+n
!      naux=nrm+n
!      nrho=naux+n
!      nybar0=1
!      nsbar0=nybar0+n*m0
!      nsize0=nsbar0+n*m0
!      nybar1=1
!      nsbar1=nybar1+n*m1
!      nsize1=nsbar1+n*m1
!      nol1=nsize1+1
Z_EPS2=P_EPS*P_EPS

!------- computation of b (RHS)

IF (K_NUPTRA == 0) THEN
  Z_F0 = P_F
  YL_B%DATA = -YD_R%DATA
ELSE
  CALL ALLOCATE_CTLVEC(YL_NUL)
  YL_NUL = 0.0_JPRB
  CALL SIMUL (INDIC,K_N,YL_NUL,Z_F0,YL_B,0)
  YL_B%DATA = - YL_B%DATA
  CALL DEALLOCATE_CTLVEC(YL_NUL)
ENDIF

CALL N1CGA (SIMUL,K_N,YD_X,YL_B,YL_Q,YD_R,YL_V,YL_RM,P_EPSNEG,Z_EPS2,&
 & K_ITER,K_IMP,K_IO,K_MODE,&
 & K_PRECO,&
 & K_M0,K_ILM0(1),K_ILM0(2),K_ILM0(3),YD_YBAR0,&
 & YD_SBAR0,P_SIZE0,&
 & K_BFGSB,&
 & K_M1,K_ILM1(1),K_ILM1(2),K_ILM1(3),K_ILM1(4),YD_YBAR1,&
 & YD_SBAR1,P_SIZE1,Z_OL1,K_SELECT,&
 & Z_RHO,Z_F0,&
 & YD_YS,YD_YRS)  

!---- final printings

!      if (imp.ge.1) then
!          call simul (indic,n,x,f,aux,iter)
!          call dpseuclid (n,aux,aux,g)
!          g=dsqrt(g)
!          if (imp.ge.2) write (io,910)
!          write (io,930) mode,iter,f,sqrt(eps2),g
!      endif
!  930 format (/" n1cg1: output mode is ",i8
!     &        /4x,"number of iterations = ",i4
!     &        /4x,"cost function f =",8x,1pe22.15
!     &        /4x,"relative |g| = ",10x,e22.15
!     &        /4x,"|g| = ",19x,d22.15)

CALL DEALLOCATE_CTLVEC(YL_B)
CALL DEALLOCATE_CTLVEC(YL_Q)
CALL DEALLOCATE_CTLVEC(YL_V)
CALL DEALLOCATE_CTLVEC(YL_RM)

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('N1CG1',1,ZHOOK_HANDLE)
END SUBROUTINE N1CG1
