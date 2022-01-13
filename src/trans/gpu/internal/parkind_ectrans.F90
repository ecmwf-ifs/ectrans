! (C) Copyright 2021- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PARKIND_ECTRANS
!
! Re-export precision-related symbols defined in fiat / parkind1, 
! and add ECTRANS-internal precision-related symbols

USE PARKIND1 
!
IMPLICIT NONE
SAVE
!
!     Real Kind of compile-time precision for internal trans use
!     ----------------------------------------------------------
!
#ifdef PARKINDTRANS_SINGLE
INTEGER, PARAMETER :: JPRBT = SELECTED_REAL_KIND(6,37)
#else
INTEGER, PARAMETER :: JPRBT = SELECTED_REAL_KIND(13,300)
#endif


!
!     Half precision
!     --------------

!!INTEGER, PARAMETER :: JPRL = 2



END MODULE PARKIND_ECTRANS
