! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_INAD_MOD
CONTAINS
SUBROUTINE FOURIER_INAD(PREEL, KFIELDS, KGL)

!**** *FOURIER_INAD* - Copy fourier data from buffer to local array - adjoint

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_INAD(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!                           KGL - local index of latitude we are currently on
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB
USE TPM_DISTR,    ONLY : D, MYSETW
USE TPM_TRANS,    ONLY : FOUBUF
USE TPM_GEOMETRY, ONLY : G

IMPLICIT NONE

REAL(KIND=JPRB),    INTENT(IN) :: PREEL(:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN) :: KGL

INTEGER(KIND=JPIM) :: JM, JF, IGLG, IPROC, IR, II, ISTA

!     ------------------------------------------------------------------

! Determine global latitude index corresponding to local latitude index KGL
IGLG = D%NPTRLS(MYSETW) + KGL - 1

! Loop over all zonal wavenumbers relevant for this latitude
DO JM = 0, G%NMEN(IGLG)
  ! Get the member of the W-set responsible for this zonal wavenumber in the "m" representation
  IPROC = D%NPROCM(JM)

  ! Compute offset in FFT work array PREEL corresponding to wavenumber JM and latitude KGL
  IR = 2 * JM + 1 + D%NSTAGTF(KGL)
  II = 2 * JM + 2 + D%NSTAGTF(KGL)

  ! Compute offset for insertion of the fields in the  m-to-l transposition buffer, FOUBUF
  ISTA = (D%NSTAGT0B(D%MSTABF(IPROC)) + D%NPNTGTB0(JM,KGL)) * 2 * KFIELDS

  ! Copy all fields from FFT work array to m-to-l transposition buffer
  DO JF = 1, KFIELDS
    FOUBUF(ISTA+2*JF-1) = PREEL(JF,IR)
    FOUBUF(ISTA+2*JF)   = PREEL(JF,II)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_INAD
END MODULE FOURIER_INAD_MOD