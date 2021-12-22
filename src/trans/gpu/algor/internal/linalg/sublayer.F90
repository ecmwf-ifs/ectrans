!*********************** SUBROUTINE SUBLAYER *************************

   SUBROUTINE SUBLAYER(PZ1,PZ2,PZ3,PX1,PX2,LDGRADPS, &
                &       PT1,PT2,PW1,PW2,PZS1,PZS2,PZPS)
!     
!-----------------------------------------------------------------------
!     
!     This software was developed and tested by Y.J. Rochon
!     with contributions from L. Garand, D.S. Turner, and S. Polavarapu
!     of the Atmospheric Science and Technology Directorate,
!     Environment Canada (EC).
!                              
!     Copyright 2006, Environment Canada, All Rights Reserved.

!     Please send any improvements to any one of the above authors
!     at firstname.lastname@ec.gc.ca (e.g yves.rochon@ec.gc.ca) and provide
!     appropriate aknowledgement through referencing to the authors or
!     the corresponding journal article (see below).
!     
!     Description:
!       
!     Determine weight coefficient contributions to assign to NWP/model 
!     variables at px1 and px2. Weights are determined from integration over
!     the range (zy1,zy2), which is within the ranges of the
!     NWP/model layer (px1,px2) and of the RTM layer (pz1,pz2). Intergrals
!     are approximated via the trapezoidal rule:

!         integral of f(x) from zy1 to zy2 = (f(zy1)+f(zy2))/2*abs(zy1-zy2)

!     This is synonomous to having an integrand linear in x.

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
!     
!     Arguments:

!     Input

!     pz1.......................Outer boundary of RTM level (above or below pz2)
!     pz2.......................Inner boundary of RTM level
!                              (position of reference RTM level)
!     pz3.......................Second outer boundary
!     px1.......................Upper boundary of NWP layer (px1<px2)
!     px2.......................Lower boundary of NWP layer
!     pt1.......................Variable value at upper NWP level.
!     pt2.......................Variable value at lower NWP level.
!     ldgradps..................Flag indicating if gradient w.r.t. Ps is
!                              required. True for yes.
!     dzs1.....................dlnP/dPs = dx1/dPs (required when ldgradps=.true.)
!     dzs2.....................dlnP/dPs = dx2/dPs (required when ldgradps=.true.)
!     pzps.....................Current gradient contribution for 
!                              weights*variable w.r.t Ps
!                              (required when ldgradps=.true.)
!     pzps2....................Current gradient contribution for 
!                              weights w.r.t Ps
!                              (required when ldgradps=.true.)

!     Output

!     pw1.......................Weight assigned to variable at upper NWP level
!     pw2.......................Weight assigned to variable at lower NWP level
!     pzps.....................Updated gradient contribution for 
!                              weights*variable w.r.t Ps
!                              (provided when ldgradps=.true.)
!     pzps2....................Updated gradient contribution for  
!                              weights w.r.t Ps 
!                              (required when ldgradps=.true.) 

!     Other

!     ztot.....................Evaluated integral
!     zg1......................Gradient of weights*variables w.r.t. px1
!     zg2......................Gradient of weights*variables w.r.t. px2
!     zy1......................Upper boundary of integral range (zy1<zy2)
!     zy2......................Lower boundary of integral range
!                              
!     External functions and subroutines: None

!     Related routines: INTAVGAD, INTAVG, INTAVGTL, INTAVGK, LAYERAVG
!     
!     Assumptions:
!     
!     1) px1<px2

!     Comments: 

!--------------------------------------------------------------------
!     
      USE PARKIND1  ,ONLY : JPIM     ,JPRB
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      IMPLICIT NONE
!     
! --- Subroutine arguments
!                              
      LOGICAL, INTENT(IN) :: LDGRADPS
      REAL(KIND=JPRB), INTENT(IN)    :: PZ1,PZ2,PZ3,PX1,PX2,PT1,PT2,PZS1,PZS2
      REAL(KIND=JPRB), INTENT(OUT)    :: PW1,PW2,PZPS
!     
! --- Local scalars and arrays
!               
      REAL(KIND=JPRB)    :: ZY1,ZY2,ZTOT,ZD,ZD2,ZW10,ZW20,ZDZ,ZDX,ZDY,ZDZD,ZDXD,ZG1,ZG2
      REAL(KIND=JPRB)    :: ZA1,ZA2,ZAA1
      INTEGER(KIND=JPIM) :: IBOT,ITOP
      REAL(KIND=JPRB)    :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!-----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('SUBLAYER',0,ZHOOK_HANDLE)

! --- Identify and set upper and lower boundaries of
!     integration/averaging layers (zy1 and zy2)

      ITOP=0
      IBOT=0
      IF (PZ1 < PZ3) THEN
         ZY1=PZ1
         IF (PX1 > PZ1) THEN
            ZY1=PX1
            ITOP=1
         ENDIF
         ZY2=PZ2
         IF (PX2 < PZ2) THEN
            ZY2=PX2
            IBOT=1
         ENDIF
      ELSE
         ZY1=PZ2
         IF (PX1 > PZ2) THEN
            ZY1=PX1
            ITOP=1
         ENDIF
         ZY2=PZ1
         IF (PX2 < PZ1) THEN
            ZY2=PX2
            IBOT=1
         ENDIF
      ENDIF

! --- Set weights for forward interpolator and linear model contribution to TLM

      ZDY=ZY2-ZY1
      ZDZ=PZ1-PZ2
      IF (ABS(ZDZ) < 1E-15) THEN
         PRINT*,'SUBLAYER: ERROR: zdz is <=0. zdz = ',ZDZ
         PRINT*,'pz1,pz2,pz3 = ',PZ1,PZ2,PZ3
         PRINT*,'px1,px2    = ',PX1,PX2
         PRINT*,'pt1,pt2    = ',PT1,PT2
         CALL ABOR1('SUBLAYER: zdz is <= 0.')
      ELSE
         ZDZD=1.0/ZDZ
      ENDIF
      PW1=(PZ1-ZY1)*ZDZD*ZDY
      PW2=(PZ1-ZY2)*ZDZD*ZDY
      ZW10=PW1
      ZW20=PW2
      ZDX=(PX2-PX1)
      IF (ABS(ZDX) < 1E-10) THEN
         PRINT*,'SUBLAYER: ERROR: zdx is <=0. zdx = ',ZDX
         PRINT*,'pz1,pz2,pz3 = ',PZ1,PZ2,PZ3
         PRINT*,'px1,px2    = ',PX1,PX2
         PRINT*,'pt1,pt2    = ',PT1,PT2
         CALL ABOR1('SUBLAYER: zdx is <= 0.')
      ELSE
         ZDXD=1.0/ZDX
      ENDIF

      ZD=(PX2-PZ2)*ZDXD
      IF (PZ1 < PZ3.AND.IBOT == 0) THEN
         PW1=PW1+PW2*ZD
         PW2=PW2*(1.0-ZD)
      ELSEIF (PZ1 > PZ3.AND.ITOP == 0) THEN
         PW2=PW2+PW1*(1.0-ZD)
         PW1=PW1*ZD
      ENDIF
      ZTOT=PT1*PW1+PT2*PW2
!     
! --- Provide linearized contribution of non-linear model component 
!     to TLM (gradients w.r.t. Ps)

!     Note: The gradient calc w.r.t Ps does not account for possible
!     re-assigment of z* in LAYERAVG and of the bottom RTM level in INTAVG*
!     to the surface model level when the original values are below
!     the model surface. This will result in differences from gradients
!     calculated using finite differences.

      IF (LDGRADPS) THEN

!        Determine gradient of 'ztot' w.r.t. px1

         IF (ITOP == 1) THEN
            ZA1=-(ZDY+(PZ1-ZY1))*ZDZD 
            ZA2=-(PZ1-ZY2)*ZDZD
         ELSE
            ZA1=0.0
            ZA2=0.0
         ENDIF

         ZD2=1.0-ZD
         IF (PZ1 < PZ3.AND.IBOT == 0) THEN
            ZAA1=ZW20*ZD*ZDXD
            ZA1=ZA1+ZA2*ZD+ZAA1
            ZA2=ZA2*ZD2-ZAA1
         ELSEIF (PZ1 > PZ3.AND.ITOP == 0) THEN
            ZAA1=ZW10*ZD*ZDXD
            ZA2=-ZAA1
            ZA1=ZAA1
         ENDIF

         ZG1=ZA1*PT1+ZA2*PT2
!                      
!        Determine gradient of 'ztot' w.r.t. px2

         IF (IBOT == 1) THEN
            ZA1=ZDZD*(PZ1-ZY1)
            ZA2=((PZ1-ZY2)-ZDY)*ZDZD
         ELSE
            ZA1=0.0
            ZA2=0.0
         ENDIF

         IF (PZ1 < PZ3.AND.IBOT == 0) THEN
            ZAA1=ZW20*ZD2*ZDXD
            ZA1=ZAA1
            ZA2=-ZAA1
         ELSEIF (PZ1 > PZ3.AND.ITOP == 0) THEN
            ZAA1=ZW10*ZD2*ZDXD
            ZA2=ZA2+ZA1*ZD2-ZAA1
            ZA1=ZA1*ZD+ZAA1
         ENDIF

         ZG2=ZA1*PT1+ZA2*PT2

!        Accumulate for gradient w.r.t. Ps

         PZPS=PZPS+ZG1*PZS1+ZG2*PZS2

      ENDIF

   IF (LHOOK) CALL DR_HOOK('SUBLAYER',1,ZHOOK_HANDLE)
   END SUBROUTINE SUBLAYER
