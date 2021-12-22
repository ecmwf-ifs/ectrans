SUBROUTINE N1CGA (SIMUL,K_N,YD_X,YD_B,YD_Q,YD_R,YD_V,YD_RM,P_EPSNEG,P_EPS2,K_ITER,K_IMP,K_IO,&
 & K_MODE,&
 & K_BFGSP,&
 & K_M0,K_NYS0,K_JMIN0,K_JMAX0,YD_YBAR0,YD_SBAR0,P_SIZE0,&
 & K_BFGSB,&
 & K_M1,K_NYS1,K_JMIN1,K_JMAX1,K_JOL,YD_YBAR1,YD_SBAR1,P_SIZE1,P_OL,&
 & K_SELECT,&
 & P_RHO,P_F0,&
 & YD_YS,YD_YRS)  

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND2  ,ONLY : JPRH

USE CONTROL_VECTORS_MOD

!-----------------------------------------------------------------------------

!     This routine solves the linear system (LS) Ax=b, for x,
!     by the possibly preconditioned Fletcher-Reeves CG algorithm.
!     The preconditioning can come from l-BFGS updates stored in ybar0,
!     sbar0 and rho.

!     The residual is denoted by r = Ax-b and the preconditioning
!     matrix by P (=+- inverse of A).

!     The first iterate is supposed to be in x.

!     Input:

!       k_n (integer): dimension of the LS to solve.
!       ??? (double precision vector): name of the vector storing the
!           matrix A, when necessary.
!       p_epsneg = epsilon for detecting negative curvature directions.
!       k_iter = max number of iterations
!       k_imp: printing levels
!         <= 0: nothing is printed,
!          = 1: error, terminal, and stopping messages,
!          = 2: also, one line per iteration.
!       k_bfgsb:
!         =0: don't build a BFGS preconditioning matrix
!         =1: build a full BFGS preconditioning matrix in pmat1, if
!             pmat1(1,1) <= 0., the matrix is re-initialized
!         =2: build a limited memory BFGS preconditioning matrix
!       k_bfgsp:
!         =0: no BFGS preconditioning
!         =1: full BFGS preconditioning with pmat0
!         =2: limited memory BFGS preconditioning

!     Output:

!       yd_x (double precision vector of dim n): approximate solution of
!           the (LS).
!       yd_b (double precision vector of dim n): RHS of the (LS).
!       yd_q (double precision vector of dim n): auxiliary vector for
!           storing Av.
!       yd_r (double precision vector of dim n): updated residual.
!       yd_v (double precision vector of dim n): inner CG directions.
!       yd_rm (double precision vector of dim n): auxiliary vector for
!           storing Pr.
!       k_iter = number of iterations
!       k_mode:
!         =0: normal terminaison (on eps2),
!         =3: A is probably indefinite (Rayleigh quotient less than
!             epsneg),
!         =4: maximum number of iterations reached,
!         =5: number of iterations exceeds 2n.
!         =6: the initial residual is zero.
!         =7: bad call to dysave.
!         =8: stop on a cost increase

!-----------------------------------------------------------------------------

IMPLICIT NONE
EXTERNAL SIMUL
INTEGER(KIND=JPIM),INTENT(IN)      :: K_M0 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: K_M1 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_N
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_X
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_B 
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_Q
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_R
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_V
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_RM 
REAL(KIND=JPRB)   ,INTENT(IN)      :: P_EPSNEG 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: P_EPS2 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: K_ITER
INTEGER(KIND=JPIM),INTENT(IN)      :: K_IMP 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_IO
INTEGER(KIND=JPIM),INTENT(OUT)     :: K_MODE 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_BFGSP 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_NYS0 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_JMIN0 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_JMAX0 
TYPE(CONTROL_VECTOR),INTENT(IN)    :: YD_YBAR0(K_M0) 
TYPE(CONTROL_VECTOR),INTENT(IN)    :: YD_SBAR0(K_M0) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: P_SIZE0 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_BFGSB 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: K_NYS1 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: K_JMIN1 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: K_JMAX1 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: K_JOL(K_M1) 
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_YBAR1(K_M1) 
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_SBAR1(K_M1) 
REAL(KIND=JPRB)   ,INTENT(OUT)     :: P_SIZE1 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: P_OL(K_M1) 
INTEGER(KIND=JPIM),INTENT(IN)      :: K_SELECT 
REAL(KIND=JPRB)   ,INTENT(OUT)     :: P_RHO(K_M0) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: P_F0 
TYPE(CONTROL_VECTOR) ,INTENT(IN)   :: YD_YS
TYPE(CONTROL_VECTOR) ,INTENT(IN)   :: YD_YRS

!-----------------------------------------------------------------------------

INTEGER(KIND=JPIM) ::INDIC
REAL(KIND=JPRB) ::Z_F

! --- parameter

REAL(KIND=JPRB) :: Z_PI
PARAMETER (Z_PI=3.1415927E0_JPRB)

! --- local variables
!     restart=.true.:
!           restarts are allowed. A restart of the CG iterations
!           occurs at iterations n, 3n/2, 9n/4, ... or when
!           (Pr_,r)>=0.2|Pr_||r| (Powell's criterion).
!     r2i = |initial residual|^2
!     v2 = |v|^2
!     pprr = previous prr
!     prr = <Pr,r>
!     r2 = |r|^2
!     avv = <Av,v>
!     rr=(|r|/|r0|)**2: relative r
!     cost = (Ax,x)/2 - (b,x)
!     x2 = |x|^2

LOGICAL :: LL_RESTART
INTEGER(KIND=JPIM) :: I_MITER,ITERR,I_NITER,I_BITER,I_YSMODE
CHARACTER :: CL_SC
REAL(KIND=JPRB) :: Z_DANGLE,Z_AUX1,Z_COST,Z_MCOST,Z_COST0,Z_RR,Z_BETA
REAL(KIND=JPRH) :: Z_R2I,Z_R2,Z_RM2,Z_PRR,PPRR,Z_V2,Z_AVV,Z_X2,Z_BTX,Z_RTX,Z_GPERP,Z_SP,Z_RAYL,Z_ALPHA,Z_T
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------

#include "dbfgsl.h"
#include "dpseuclid.h"
#include "dysave.h"

!-----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('N1CGA',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------------

! --- initialisation
!     i_miter = max iterations
!     i_niter = iteration nbr of the next forced restart
!     iterr = iteration nbr of the previous restart
!     i_biter = bad iteration nbr (will produce failure)

I_MITER=K_ITER
I_NITER=K_N
ITERR=1
I_BITER=INT(SQRT(REAL(K_N*K_N*K_N,JPRB)))
K_ITER=0
Z_RR=1.0_JPRB
Z_GPERP=0.0_JPRB
IF (K_SELECT == 0) CL_SC='S'

!     --- In my experience, it is better no doing restart

LL_RESTART=.FALSE.

!     --- force initialize the l-BFGS matrix if appropriate

K_NYS1=0

!     --- initial residual r = Ax-b
!         r2i = |initial residual|^2

!      call simul (indic,n,x,f,r,iter)
!      appel au simul. inutile car deja effectue a l'exterieur !

!     --- change of variable for x, r and b (RHS)

YD_X%DATA = YD_X%DATA * YD_YS%DATA
YD_R%DATA = 0.5_JPRB * YD_R%DATA * YD_YS%DATA
YD_B%DATA = 0.5_JPRB * YD_B%DATA * YD_YS%DATA
YD_RM = YD_R
CALL DPSEUCLID (K_N,YD_R,YD_R,Z_R2I)
Z_R2=Z_R2I
IF (Z_R2I == 0.0_JPRB) THEN
  K_MODE=6
  IF (K_IMP > 0) WRITE (K_IO,900)
  GOTO 1000
ENDIF
900 FORMAT(/4X,"n1cga-WARNING: zero initial residual")

!     --- rm = Pr
!         prr = <Pr,r>

!      rm = r     ! Inutile a mon avis

IF (K_BFGSP == 2) THEN
  CALL DBFGSL (K_N,YD_RM,K_M0,K_NYS0,K_JMIN0,K_JMAX0,YD_YBAR0,YD_SBAR0,P_RHO,P_SIZE0)
ENDIF
CALL DPSEUCLID (K_N,YD_RM,YD_RM,Z_RM2)
CALL DPSEUCLID (K_N,YD_R,YD_RM,Z_PRR)

!     --- initial search direction v
!         v2 = |v|^2

YD_V%DATA = -YD_RM%DATA
Z_V2=Z_RM2

! --- cost
!     cost = current cost
!     mcost = min cost encountered

Z_COST=0.E0_JPRB
Z_MCOST=0.E0_JPRB

! --- printings

IF (K_IMP >= 2) THEN
  WRITE (K_IO,901) SQRT(Z_R2I)
ENDIF
901 FORMAT (/4X,"n1cga: initial residual |r0| =",1PD13.6)

! --- start the CG iterations at 100, end at 1000

100 CONTINUE
K_ITER=K_ITER+1

!     --- check iter (force convergence, however)

IF (K_ITER > I_MITER) THEN
  K_MODE=4
  K_ITER=K_ITER-1
  IF (K_IMP > 0) WRITE (K_IO,910)
  GOTO 1000
ENDIF
910 FORMAT(/4X,"n1cga: max inner iterations")

!     --- too many iterations ?

!     if (iter > biter) then
!         mode=5
!         if (imp > 0) write (io,911) biter
!         goto 1000
!     endif
! 911 format(/4x,"n1cga-ERROR: number of inner iterations exceeds",i6)

!     --- restart ?

IF (LL_RESTART.AND.(Z_GPERP > 0.2E0_JPRB.OR.K_ITER > I_NITER)) THEN
  IF (K_ITER > I_NITER) I_NITER=I_NITER*3/2
  ITERR=K_ITER
  IF (K_IMP >= 2) WRITE (K_IO,912)
  912 FORMAT(4X,"n1cga: LL_RESTART")

!         --- recompute the residuals r and rm

!         --- inverse change of variable

  YD_X%DATA = YD_X%DATA * YD_YRS%DATA

  CALL SIMUL (INDIC,K_N,YD_X,Z_F,YD_R,K_ITER)

!         --- change of variable for x and r

  YD_X%DATA = YD_X%DATA * YD_YS%DATA
  YD_R%DATA = 0.5_JPRB * YD_R%DATA * YD_YS%DATA
  YD_RM = YD_R
  CALL DPSEUCLID (K_N,YD_R,YD_R,Z_R2)

  IF (K_BFGSP == 2) THEN
    CALL DBFGSL (K_N,YD_RM,K_M0,K_NYS0,K_JMIN0,K_JMAX0,YD_YBAR0,YD_SBAR0,P_RHO,P_SIZE0)
  ENDIF
  CALL DPSEUCLID (K_N,YD_RM,YD_RM,Z_RM2)
  CALL DPSEUCLID (K_N,YD_R,YD_RM,Z_PRR)

!         --- and the search direction v

  YD_V%DATA = -YD_RM%DATA
  Z_V2=Z_RM2
ENDIF

!     --- q = Av
!         avv = <v,q> = <Av,v>

!     --- inverse of change of variable for v

YD_V%DATA = YD_V%DATA * YD_YRS%DATA

CALL SIMUL (INDIC,K_N,YD_V,Z_F,YD_Q,K_ITER)

!     --- change of variable for v and q

YD_V%DATA = YD_V%DATA * YD_YS%DATA
YD_Q%DATA = 0.5_JPRB * YD_YS%DATA * YD_Q%DATA
YD_Q%DATA = YD_Q%DATA + YD_B%DATA

CALL DPSEUCLID (K_N,YD_V,YD_Q,Z_AVV)

!     --- check positive definiteness

IF (K_IMP > 0.AND.Z_V2 /= 0.0_JPRB) Z_RAYL=Z_AVV/Z_V2
IF (Z_AVV <= P_EPSNEG*Z_V2) THEN
  K_MODE=3
  IF (K_IMP > 0.AND.Z_V2 /= 0.0_JPRB) WRITE (K_IO,913) Z_RAYL
  GOTO 1000
ENDIF
913 FORMAT(/4X,"n1cga-ERROR: non positive definite matrix,",&
 & " <Av,YD_V>/|YD_V|^2 = ",1PE9.2)  

!     --- update the preconditioning matrix when appropriate
!         y=q=Av, s=v, ys=avv

IF (K_BFGSB == 2) THEN

!         --- limited memory inverse BFGS update

  CL_SC='-'
  CALL DYSAVE (K_N,YD_Q,YD_V,Z_AVV,K_M1,K_NYS1,K_JMIN1,K_JMAX1,YD_YBAR1,YD_SBAR1,K_SELECT,&
   & K_ITER,P_OL,K_JOL,P_EPSNEG,P_SIZE1,I_YSMODE,4,K_IO)  
  IF (I_YSMODE == 0) CL_SC='S'
  IF (I_YSMODE > 0) THEN
    K_MODE=7
    IF (LHOOK) CALL DR_HOOK('N1CGA',1,ZHOOK_HANDLE)
    RETURN
  ENDIF
ENDIF
914 FORMAT(4X,"n1cga: initialization of the BFGS matrix with P_OL",&
 & " factor =",1PE13.6)  

!     --- new iterate x and residual r

Z_ALPHA=Z_PRR/Z_AVV
YD_X%DATA = YD_X%DATA + Z_ALPHA*YD_V%DATA
CALL DPSEUCLID (K_N,YD_X,YD_X,Z_X2)
YD_R%DATA = YD_R%DATA + Z_ALPHA*YD_Q%DATA

! --- stop if the cost increases

Z_MCOST=MIN(Z_COST,Z_MCOST)
CALL DPSEUCLID (K_N,YD_R,YD_X,Z_RTX)
CALL DPSEUCLID (K_N,YD_B,YD_X,Z_BTX)
Z_COST = 0.5E0_JPRB*(Z_RTX - Z_BTX)
Z_COST0 = P_F0 + 2*Z_COST
IF (Z_COST > 0.999999E0_JPRB*Z_MCOST) THEN
  K_MODE=8
  IF (K_IMP > 0) WRITE (K_IO,915)
  GOTO 1000
ENDIF
915 FORMAT(/4X,"n1cga-ERROR: Z_COST increase (rounding error)")

!     --- r2=|r|^2

CALL DPSEUCLID (K_N,YD_R,YD_R,Z_R2)

!     --- prepare to check conjugacy
!         gperp = (rmo,r)/|rmo|/|r| == 0 ?

CALL DPSEUCLID (K_N,YD_RM,YD_R,Z_GPERP)
Z_T=SQRT(Z_R2*Z_RM2)
IF (Z_T /= 0.0_JPRB) THEN
  Z_GPERP=Z_GPERP/Z_T
ELSE
  Z_GPERP=0.0_JPRB
ENDIF

!     --- new rm
!         rm2 = |rm|^2
!         new pprr=prr
!         prr=<Pr,r>

YD_RM = YD_R

IF (K_BFGSP == 2) THEN
  CALL DBFGSL (K_N,YD_RM,K_M0,K_NYS0,K_JMIN0,K_JMAX0,YD_YBAR0,YD_SBAR0,P_RHO,P_SIZE0)
ENDIF
CALL DPSEUCLID (K_N,YD_RM,YD_RM,Z_RM2)
PPRR=Z_PRR
CALL DPSEUCLID (K_N,YD_R,YD_RM,Z_PRR)

!     --- some printings

IF (K_IMP >= 2) THEN
  CALL DPSEUCLID (K_N,YD_RM,YD_RM,Z_SP)
  Z_DANGLE=0.E0_JPRB
  IF (SQRT(Z_SP*Z_R2) /= 0.0_JPRB) THEN
    Z_AUX1 = Z_PRR/SQRT(Z_SP*Z_R2)
    IF (Z_AUX1 < -1.E0_JPRB) Z_AUX1 = -1.E0_JPRB
    IF (Z_AUX1 > 1E0_JPRB) Z_AUX1 = 1.E0_JPRB
    Z_DANGLE=ACOS(Z_AUX1)/Z_PI*180.E0_JPRB
  ENDIF
  IF (K_ITER == 1) WRITE (K_IO,905)
  WRITE (K_IO,902) K_ITER,Z_COST0,SQRT(Z_R2/Z_R2I),Z_DANGLE,Z_RAYL,Z_GPERP,&
   & Z_ALPHA,SQRT(Z_X2),CL_SC  
  CALL FLUSH (K_IO)
ENDIF
905 FORMAT (/4X,"n1cga: K_ITER",6X,"Z_COST",7X,"|YD_R|/|r0|",2X,"<(YD_R,Pr)",&
 & "(Av,YD_V)/|YD_V|^2",4X,"conj",3X,"Z_ALPHA",4X,"|YD_X|",3X,"S")  
902 FORMAT (4X,"n1cga:",I5,1X,1PE13.6,1X,E12.6,1X,0PF5.2,1X,1PE12.6,&
 & 1X,E9.2,1X,E6.0,1X,E7.1,1X,A1)  

!     --- stop on the relative residual (eps2)

Z_RR=Z_R2/Z_R2I
IF (Z_RR <= P_EPS2) THEN
  K_MODE=0
  IF (K_IMP > 0) WRITE (K_IO,916) SQRT(Z_RR),SQRT(P_EPS2)
  GOTO 1000
ENDIF
916 FORMAT(/4X,"n1cga: stopping criterion satisfied, |YD_R|/|r0| = ",&
 & 1PE8.2," <= ",E8.2)  

!     --- FR beta, new v, v2=|v|^2

80 CONTINUE 
Z_BETA=Z_PRR/PPRR
YD_V%DATA = -YD_RM%DATA + Z_BETA*YD_V%DATA
CALL DPSEUCLID (K_N,YD_V,YD_V,Z_V2)

!     --- loop

GOTO 100

! --- terminaison

1000 CONTINUE
P_EPS2=Z_RR
IF (K_IMP == 1) WRITE (K_IO,920) K_ITER,SQRT(Z_R2I),SQRT(Z_R2)
920 FORMAT(4X,"n1cga: I_NITER=",I5,", |YD_R|:",D13.6," -> ",D13.6)

!     ---- inverse change of var. for x

YD_X%DATA = YD_X%DATA * YD_YRS%DATA

!-----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('N1CGA',1,ZHOOK_HANDLE)
END SUBROUTINE N1CGA
