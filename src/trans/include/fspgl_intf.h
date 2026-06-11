! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 

ABSTRACT INTERFACE

SUBROUTINE FSPGL_INTF (KM,KSL,KDGL,KFIELDS,PR1MU2,PFIELD,&
                     & KPTRU,KFLDUV,KFLDSC,KFLDPTRUV)

! interface of a procedure to be executed in fourier space by inv_trans before transposition
USE PARKIND1  ,ONLY : JPIM     ,JPRB

INTEGER(KIND=JPIM),INTENT(IN)           :: KM
INTEGER(KIND=JPIM),INTENT(IN)           :: KSL
INTEGER(KIND=JPIM),INTENT(IN)           :: KDGL
REAL(KIND=JPRB)   ,INTENT(IN)           :: PR1MU2(KDGL)
INTEGER(KIND=JPIM),INTENT(IN)           :: KFIELDS
REAL(KIND=JPRB)   ,INTENT(INOUT),TARGET :: PFIELD(2*KFIELDS,0:KDGL+1)
INTEGER(KIND=JPIM),INTENT(IN)           :: KPTRU
INTEGER(KIND=JPIM),INTENT(IN)           :: KFLDUV
INTEGER(KIND=JPIM),INTENT(IN)           :: KFLDSC
INTEGER(KIND=JPIM),INTENT(IN)           :: KFLDPTRUV(KFLDUV)

END SUBROUTINE FSPGL_INTF

END INTERFACE