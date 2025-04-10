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
SUBROUTINE ETRANS_END(CDMODE)

!**** *ETRANS_END* - Terminate transform package 

!     Purpose.
!     --------
!     Terminate transform package. Release all allocated arrays.

!**   Interface.
!     ----------
!     CALL ETRANS_END

!     Explicit arguments : None
!     -------------------- 

!     Method.
!     -------

!     Externals.  None
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Nmiri       15-Nov-2007 Phasing with TFL 32R3
!        A.Bogatchev   16-Sep-2010 Phasing cy37 after G.Radnoti

!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER*5, OPTIONAL,  INTENT(IN) :: CDMODE

END SUBROUTINE ETRANS_END
END INTERFACE
