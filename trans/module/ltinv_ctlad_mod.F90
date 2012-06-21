MODULE LTINV_CTLAD_MOD
CONTAINS
SUBROUTINE LTINV_CTLAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

!**** *LTINV_CTLAD* - Control routine for inverse Legandre transform - adj.

!     Purpose.
!     --------
!     Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTL(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KFLDPTRUV(:) - field pointer array for vor./div.
!     KFLDPTRSC(:) - field pointer array for PSPSCALAR
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition  

!     Method.
!     -------

!     Externals.  
!     ----------  
!                 

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-06-03

!     ------------------------------------------------------------------
#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE LTINVAD_MOD
USE TRLTOM_MOD
USE YOMGSTATS, ONLY : LSYNCSTATS

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL_B ,OPTIONAL, INTENT(OUT)  :: PSPSC3A(:,:,:)
REAL_B ,OPTIONAL, INTENT(OUT)  :: PSPSC3B(:,:,:)
REAL_B ,OPTIONAL, INTENT(OUT)  :: PSPSC2(:,:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

INTEGER_M :: JM,IM,IBLEN,ILEI2,IDIM1

!     ------------------------------------------------------------------

ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IDIM1 = 2*KF_OUT_LT
IBLEN = D%NLENGT0B*2*KF_OUT_LT
ALLOCATE(FOUBUF_IN(IBLEN))
CALL GSTATS(180,0)
CALL TRLTOM(FOUBUF,FOUBUF_IN,2*KF_OUT_LT)
CALL GSTATS(180,1)
DEALLOCATE(FOUBUF)

IF (.NOT.LSYNCSTATS) CALL GSTATS(1648,0)
IF(KF_OUT_LT > 0) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTINVAD(IM,JM,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,&
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
IF (.NOT.LSYNCSTATS) CALL GSTATS(1648,1)

DEALLOCATE(FOUBUF_IN)

!     ------------------------------------------------------------------

END SUBROUTINE LTINV_CTLAD
END MODULE LTINV_CTLAD_MOD
