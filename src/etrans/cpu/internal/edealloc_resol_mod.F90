! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EDEALLOC_RESOL_MOD
CONTAINS
SUBROUTINE EDEALLOC_RESOL(KRESOL)

!**** *EDEALLOC_RESOL_MOD* - Deallocations of a resolution

!     Purpose.
!     --------
!     Release allocated arrays for a given resolution

!**   Interface.
!     ----------
!     CALL EDEALLOC_RESOL_MOD

!     Explicit arguments : KRESOL : resolution tag
!     --------------------

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!       Original : 09-Jul-2013 from etrans_end
!       B. Bochenek (Apr 2015): Phasing: update
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : LENABLED, NOUT
USE TPM_DISTR       ,ONLY : D
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FIELDS      ,ONLY : F
#ifdef WITH_FFT992
USE TPM_FFT         ,ONLY : T
#endif
USE TPM_FFTW        ,ONLY : TW,DESTROY_PLANS_FFTW
USE TPM_FLT         ,ONLY : S

USE ESET_RESOL_MOD  ,ONLY : ESET_RESOL

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL

!     ------------------------------------------------------------------

IF (.NOT.LENABLED(KRESOL)) THEN

  WRITE(UNIT=NOUT,FMT='('' EDEALLOC_RESOL WARNING: KRESOL = '',I3,'' ALREADY DISABLED '')') KRESOL

ELSE

  CALL ESET_RESOL(KRESOL)

  !TPM_DISTR
  DEALLOCATE(D%NFRSTLAT,D%NLSTLAT,D%NPTRLAT,D%NPTRFRSTLAT,D%NPTRLSTLAT)
  DEALLOCATE(D%LSPLITLAT,D%NSTA,D%NONL,D%NGPTOTL,D%NPROCA_GP)

  IF(D%LWEIGHTED_DISTR) THEN
    DEALLOCATE(D%RWEIGHT)
  ENDIF

  IF(.NOT.D%LGRIDONLY) THEN

    DEALLOCATE(D%MYMS,D%NUMPP,D%NPOSSP,D%NPROCM,D%NDIM0G,D%NASM0,D%NATM0)
    DEALLOCATE(D%NLATLS,D%NLATLE,D%NPMT,D%NPMS,D%NPMG,D%NULTPP,D%NPROCL)
    DEALLOCATE(D%NPTRLS,D%NALLMS,D%NPTRMS,D%NSTAGT0B,D%NSTAGT1B,D%NPNTGTB0)
    DEALLOCATE(D%NPNTGTB1,D%NLTSFTB,D%NLTSGTB,D%MSTABF)
    DEALLOCATE(D%NSTAGTF)

#ifdef WITH_FFT992
  !TPM_FFT
    DEALLOCATE(T%TRIGS,T%NFAX)
#endif
  !TPM_FFTW
    CALL DESTROY_PLANS_FFTW
  !TPM_GEOMETRY
    DEALLOCATE(G%NMEN,G%NDGLU)

  ELSE

    DEALLOCATE(G%NLOEN)

  ENDIF

  LENABLED(KRESOL)=.FALSE.

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE EDEALLOC_RESOL
END MODULE EDEALLOC_RESOL_MOD
