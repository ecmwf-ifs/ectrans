! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


INTERFACE
SUBROUTINE ETIBIHIE(KDLON,KDGL,KNUBI,KDLUX,KDGUX,&
 & KSTART,KDLSM,PGPBI,LDBIX,LDBIY,KDADD)  

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)               :: KNUBI
INTEGER(KIND=JPIM),INTENT(IN)               :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)               :: KDLSM 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDLUX 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDGUX 
INTEGER(KIND=JPIM),INTENT(IN)               :: KDADD
REAL(KIND=JPRB),INTENT(INOUT)               :: PGPBI(KSTART:KDLSM+KDADD,KNUBI,1:KDGL+KDADD) 
LOGICAL,INTENT(IN)                          :: LDBIX 
LOGICAL,INTENT(IN)                          :: LDBIY 

END SUBROUTINE ETIBIHIE
END INTERFACE
