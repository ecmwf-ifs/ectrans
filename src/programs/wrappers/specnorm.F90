! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE SPECNORM(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD,JPRM
#if PRECOPT == 1
USE ECTRANS_MOD_SP, ONLY : SPECNORM_SP => SPECNORM 
#elif PRECOPT == 2
USE ECTRANS_MOD_DP, ONLY : SPECNORM_DP => SPECNORM 
#endif

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PMET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PNORM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
IF(JPRB == JPRM) THEN
#if PRECOPT == 1
CALL SPECNORM_SP(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)
#endif
ELSE IF (JPRB == JPRD) THEN
#if PRECOPT == 2
CALL SPECNORM_DP(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)
#endif
ENDIF
END SUBROUTINE SPECNORM

