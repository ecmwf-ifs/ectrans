!************************* SUBROUTINE INTAVG ***************************

   SUBROUTINE INTAVG(PVLEV,PVI,KNIDIM,KNI,KNPROF,KNO,PPO,PVO)

!-----------------------------------------------------------------------

!     This software was developed and tested by Y.J. Rochon
!     with contributions from L. Garand, D.S. Turner, and S. Polavarapu
!     of the Atmospheric Science and Technology Directorate,
!     Environment Canada (EC).

!     Copyright 2006, Environment Canada, All Rights Reserved.

!     Please send any improvements to any one of the above authors
!     at firstname.lastname@ec.gc.ca (e.g yves.rochon@ec.gc.ca) and provide
!     appropriate aknowledgement through referencing to the authors or
!     the corresponding journal article (see below).

!     Description:

!     Performs forward interpolation (IFLAG=0) based on piecewise 
!     weighted averaging in log of pressure of one-dimensional vectors.

!     Method:
!     - Piecewise weighted interpolation in ln(P).
!     - Journal reference:
!       Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
!       Jacobian mapping between coordinate systems in data assimilation,
!       Q. J. of the Royal Met. Soc., Accepted pending revisions on 25
!       Sept. 2006.

!     History:

!     Version   Date        Comment
!     -------   ----        -------
!     1         10/2005     Original F77 code by Y.J. Rochon with
!                           contributions from S. Polavarapu. (EC)
!     2         10/2006     Adaptation for F90 and for improved consistency
!                           with RTTOV code format (Y.J. Rochon, EC).
!                           Argument names of original code were not changed.

!     Arguments:

!     Input


!     PVLEV(KNIDIM,KNPROF).....Vertical levels, pressure (source/model)
!     PVI(KNIDIM,KNPROF).......Input vector to be interpolated (model)
!     KNIDIM...................Dimension of input levels (source)
!     KNI......................Number of input levels (source)
!     KNPROF...................Number of profiles (source)
!     KNO......................Number of output levels (model)
!     PPO(KNO).................Vertical levels, pressure (destination)

!     Output

!     PVO(KNO,KNPROF)..........Interpolated vector on RTM levels
!                              
!     External functions and subroutines:
!     
!     LAYERAVG: Performs integration between points ZLNPO(JO-1) and
!               ZLNPO(JO+1) with  piecewise linear weighting having
!               weights of zero at  JO-1 and JO+1 and max weight at JO.
!     
!     Related routines: INTAVGAD, INTAVGTL, INTAVGK, SUBLAYER, LAYERAVG
!     
!     Assumptions:

!     1) PVLEV(i)<PVLEV(i+1) & PPO(i)<PPO(i+1)
!     
!     Comments:
!     
!     1) Types of operations to which routine LAYERAVG contributes

!        Cases  a: Interpolate (NWP to RTM)
!               b: Apply TLM increment operator (NWP to RTM)
!               c: Apply ADJ operator (RTM to NWP)
!               d: Provide K-matrix (NWP to RTM increment operator)

!               a) Application as forward interpolator

!                  Y = sum(H*PVI)
!                    = PVO

!               with ZPZ=H on output of LAYERAVG

!               b) TLM case:

!                  dY = sum(H*PVI) + PPS*sum(PVIG*dH/dPs))
!                     = TOTAL(ZPZ*PVI)+ZPVOPS*PPS
!                     = PVO

!               where

!                     dPZ(i)/dPs = sum_k (dH(i)/dPVLEV(k) * dPVLEV(k)/dPs)
!                                = sum_k (dH(i)/dZLNPI(k) * zpresb(k)/PVLEV(k))

!               with ZPZ(k)=zpresb(k)/PVLEV(k) on input to LAYERAVG and
!                    ZPZ(i)=H(i) on output of LAYERAVG

!               c) Adjoint case:

!                  PVI(i,jn) = PVI(i,jn) + sum_jo (PVO(jo,jn)*H(i))
!                            = PVI(i,jn) + sum_jo (ZPZ(jo,i)*PVO(jo,jn))

!                  PPS(jn) = PPS(jn)
!                            + sum_jo (PVO(jo,jn)*sum_i (PVIG(i,jn)*dH(i)/dPs))
!                          = PPS(jn) + sum_jo (ZPZPS(jo)*PVO(jn,jn))

!               d) K-matrix:

!                  ZK(1:n,jn)=H(1:n)   (COORD='PRESSURE')

!                  or, when COORD is not 'PRESSURE',

!                  ZK(1:n+1,jn)=[H(1:n),sum(PVIG*dH/dPs))]
!               
!     2) The gradient calc w.r.t Ps does not account for possible
!     re-assigment of lowest contributing RTM level in INTAVG* and LAYERAVG
!     to the surface model level when the original values are below
!     the model surface. This will result in differences from gradients
!     calculated using finite differences.

!--------------------------------------------------------------------

      USE PARKIND1  ,ONLY : JPIM     ,JPRB
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      IMPLICIT NONE

! --- Subroutine arguments

      INTEGER(KIND=JPIM), INTENT(IN) :: KNIDIM, KNI, KNO, KNPROF
      REAL(KIND=JPRB), INTENT(IN)    :: PVLEV(KNIDIM,KNPROF)
      REAL(KIND=JPRB), INTENT(IN)    :: PPO(KNO),PVI(KNIDIM,KNPROF)
      REAL(KIND=JPRB), INTENT(INOUT) :: PVO(KNO,KNPROF)

! --- Local scalars and arrays
!                                
      INTEGER(KIND=JPIM) :: JO, JN
      LOGICAL :: LLGRADPS    ! Indicates if gradient w.r.t. Ps is required
      REAL(KIND=JPRB)    :: ZLNPI(KNI),ZPZ(KNI),ZPVI(KNI),ZPS(KNI)
      REAL(KIND=JPRB)    :: ZLNPO(KNO),ZPZPS
      REAL(KIND=JPRB)    :: ZHOOK_HANDLE


#include "layeravg.h"
!-----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('INTAVG',0,ZHOOK_HANDLE)
!               
! --- Initialization
!                  
      LLGRADPS=.FALSE.                            
      ZLNPO(1:KNO)=LOG(PPO(1:KNO))
      ZPZ(1:KNI)=0.0
      ZPZPS=0.0

! --- Loop over profiles

      DO JN = 1, KNPROF
         ZLNPI(1:KNI)=LOG(PVLEV(1:KNI,JN))
!                  
! --     Identify index of bottom RTM level relative to model surface.         
!        Impose model surface as lower boundary of lowest usable RTM level.
! 
         JO=KNO
  100    CONTINUE
         IF (PVLEV(KNI,JN) < PPO(JO)) THEN
            ZLNPO(JO)= ZLNPI(KNI)
            PVO(JO,JN)=0.0
            JO=JO-1
            GO TO 100
         ENDIF

         ZPVI(1:KNI)=PVI(1:KNI,JN)

         DO JO=1,KNO
            CALL LAYERAVG(LLGRADPS,ZLNPO,ZLNPI,ZPVI,KNO,KNI, &
                    &      JO,ZPZ,ZPS,ZPZPS)
            PVO(JO,JN)=SUM(ZPZ(1:KNI)*ZPVI(1:KNI))
         ENDDO
      ENDDO

   IF (LHOOK) CALL DR_HOOK('INTAVG',1,ZHOOK_HANDLE)
   END SUBROUTINE INTAVG
