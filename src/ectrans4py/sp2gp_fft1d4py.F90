! (C) Copyright 2025- Meteo-France.
! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! ===================================
! NOTE: this subroutine is not tested
! ===================================
SUBROUTINE SP2GP_FFT1D4PY(KSIZES, KTRUNC, PSPEC, KSIZEG, PGPT)
! ** PURPOSE
!    Transform spectral coefficients into grid-point values,
!    for a 1D array (vertical section academic model)
!
! ** DUMMY ARGUMENTS
!    KSIZES size of PSPEC
!    KTRUNC: troncature
!    PSPEC: spectral coefficient array
!    KSIZEG: size of grid-point field (with extension zone)
!    PGPT: grid-point field
!
! ** AUTHOR
!    26 March 2015, A. Mary, from utilities/pinuts/module/fa_datas_mod.F90
!
! ** MODIFICATIONS
!
! I. Dummy arguments declaration

USE ISO_FORTRAN_ENV, ONLY: INT64, REAL64
USE TPM_FFTW_DP, ONLY: EXEC_FFTW

IMPLICIT NONE

INTEGER(KIND=INT64), INTENT(IN) :: KSIZES
INTEGER(KIND=INT64), INTENT(IN) :: KTRUNC
REAL(KIND=REAL64), DIMENSION(KSIZES), INTENT(IN) :: PSPEC
INTEGER(KIND=INT64), INTENT(IN) :: KSIZEG
REAL(KIND=REAL64), DIMENSION(KSIZEG), INTENT(OUT) :: PGPT

INTEGER(KIND=INT64) :: NFTM, NDGLSUR
REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: SP2
INTEGER(KIND=INT64), PARAMETER :: NZERO=0

NDGLSUR = KSIZEG+MOD(KSIZEG,2)+2
NFTM    = 2*(KTRUNC+1)
ALLOCATE(SP2(1,NDGLSUR*NFTM))
SP2      = 0.0
SP2(1,:) = CONVRT2FFT(PSPEC,NZERO,KTRUNC,NDGLSUR)
CALL EXEC_FFTW(1, INT(KSIZEG,4), (INT(KSIZEG,4)/2+1)*2, 1, 1, .FALSE., SP2(:,1:KSIZEG))
PGPT(:) = SP2(1,1:KSIZEG)

CONTAINS

! from utilities/pinuts/module/fa_datas_mod.F90
! and utilities/pinuts/module/array_lib_mod.F90

FUNCTION CONVRT2FFT(IN,X,Y,N) RESULT(OU)
REAL(KIND=REAL64),DIMENSION(:),INTENT(IN)      :: IN
INTEGER(KIND=INT64),INTENT(IN)                       :: X, Y, N
REAL(KIND=REAL64),DIMENSION(N*2*(X+1))         :: OU

INTEGER(KIND=INT64),DIMENSION(2*(X+1),(N/2))         :: MINQ
INTEGER(KIND=INT64),DIMENSION((N/2),2*(X+1))         :: TMINQ
REAL(KIND=REAL64),DIMENSION(2*(X+1),(N/2))     :: OMINQ, EMINQ
REAL(KIND=REAL64),DIMENSION((N/2),2*(X+1))     :: TOMINQ, TEMINQ
REAL(KIND=REAL64),DIMENSION(N*(X+1))           :: OINI, EINI
REAL(KIND=REAL64), PARAMETER                   :: ZZERO=0.0

CALL SPLIT_ODEV(IN,OINI,EINI)
MINQ   = MASQ(X,Y,N)
OMINQ  = UNPACK(OINI,MINQ == 1,ZZERO)
TOMINQ = TRANSPOSE(OMINQ)
EMINQ  = UNPACK(EINI,MINQ == 1,ZZERO)
TEMINQ = TRANSPOSE(EMINQ)
TMINQ  = 1
OINI   = PACK(TOMINQ,TMINQ > 0)
EINI   = PACK(TEMINQ,TMINQ > 0)
OU     = MIX_ODEV(OINI,EINI)
END FUNCTION CONVRT2FFT

FUNCTION MASQ(X,Y,N) RESULT(T)
INTEGER(KIND=INT64),INTENT(IN)                       :: X, Y, N
INTEGER(KIND=INT64),DIMENSION(1:2*(X+1),1:(N/2))     :: T

INTEGER(KIND=INT64)                                  :: I, J
INTEGER(KIND=INT64),DIMENSION(0:X)                   :: KM
INTEGER(KIND=INT64),DIMENSION(0:Y)                   :: KN
CALL ELLIPS(INT(X,4),INT(Y,4),INT(KN,4),INT(KM,4))
T = 0
DO I=0,Y
  DO J=0,2*KN(I)+1
    T(J+1,I+1)=1
  END DO
END DO
END FUNCTION MASQ

FUNCTION MIX_ODEV(TO,TE) RESULT(T)
REAL(KIND=REAL64),DIMENSION(:),INTENT(IN)        :: TO,TE
REAL(KIND=REAL64),DIMENSION(SIZE(TO)+SIZE(TE))   :: T

INTEGER(KIND=INT64) :: I

DO I=1,(SIZE(TO)+SIZE(TE))/2
  T((2*I)-1)=TE(I)
  T(2*I)=TO(I)
END DO
END FUNCTION MIX_ODEV

SUBROUTINE SPLIT_ODEV(T,TO,TE)
REAL(KIND=REAL64),DIMENSION(:),INTENT(IN)          :: T
REAL(KIND=REAL64),DIMENSION(SIZE(T)/2),INTENT(OUT) :: TO,TE

INTEGER(KIND=INT64) :: I

DO I=1,SIZE(T)/2
  TO(I)=T(2*I)
  TE(I)=T((2*I)-1)
END DO
END SUBROUTINE SPLIT_ODEV

END SUBROUTINE SP2GP_FFT1D4PY