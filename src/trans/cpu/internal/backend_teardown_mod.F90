! (C) Copyright 2026- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE BACKEND_TEARDOWN_MOD

USE TPM_FLT,   ONLY: S, FLT_RESOL
USE TPM_FFTW,  ONLY: TW, FFTW_RESOL
#ifdef WITH_FFT992
USE TPM_FFT992,ONLY: T992, FFT992_RESOL
#endif
USE TPM_TRANS, ONLY: FOUBUF, FOUBUF_IN

CONTAINS
SUBROUTINE BACKEND_TEARDOWN

IMPLICIT NONE

NULLIFY(S)
IF (ALLOCATED(FLT_RESOL)) DEALLOCATE(FLT_RESOL)

NULLIFY(TW)
IF (ALLOCATED(FFTW_RESOL)) DEALLOCATE(FFTW_RESOL)

#ifdef WITH_FFT992
NULLIFY(T992)
IF (ALLOCATED(FFT992_RESOL)) DEALLOCATE(FFT992_RESOL)
#endif

IF(ALLOCATED(FOUBUF_IN)) DEALLOCATE(FOUBUF_IN)
IF(ALLOCATED(FOUBUF)) DEALLOCATE(FOUBUF)

END SUBROUTINE BACKEND_TEARDOWN
END MODULE BACKEND_TEARDOWN_MOD
