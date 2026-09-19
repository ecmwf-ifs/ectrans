! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_VIEW_CTL_MOD
CONTAINS
SUBROUTINE LTDIR_VIEW_CTL(KF_FS, YDSPVVOR, YDSPVDIV, YDSPVSCALAR)

!**** *LTDIR_VIEW_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_VIEW_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : LALLOPERM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

USE LTDIR_VIEW_MOD       ,ONLY : LTDIR_VIEW
USE TRLTOM_MOD      ,ONLY : TRLTOM
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS

TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)


INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IBLEN = D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF

CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
CALL GSTATS(153,1)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)

! Direct Legendre transform

CALL GSTATS(103,0)
ILED2 = 2*KF_FS
CALL GSTATS(1645,0)
IF(KF_FS>0) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIR_VIEW(IM,JM,KF_FS,SIZE(YDSPVVOR),SIZE(YDSPVSCALAR),ILED2, &
                  & YDSPVVOR, YDSPVDIV, YDSPVSCALAR)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1645,1)

IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)
CALL GSTATS(103,1)

!     -----------------------------------------------------------------

END SUBROUTINE LTDIR_VIEW_CTL
END MODULE LTDIR_VIEW_CTL_MOD
