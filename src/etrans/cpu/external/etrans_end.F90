! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


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
!        R. El Khatib 02-Mar-2012 Support for mixed multi-resolutions
!        R. El Khatib 09-Jul-2013 LENABLED
!        R. El Khatib 01-Set-2015 Support for FFTW
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : MSETUP0, NCUR_RESOL, NDEF_RESOL, NMAX_RESOL, LENABLED
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL, NPRCIDS
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : F, FIELDS_RESOL
#ifdef WITH_FFT992
USE TPM_FFT         ,ONLY : T, FFT_RESOL
USE TPMALD_FFT      ,ONLY : TALD, ALDFFT_RESOL
#endif
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
USE TPM_FLT         ,ONLY : S, FLT_RESOL
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPMALD_DIM      ,ONLY : RALD, ALDDIM_RESOL
USE TPMALD_DISTR    ,ONLY : DALD, ALDDISTR_RESOL
USE TPMALD_FIELDS   ,ONLY : FALD, ALDFIELDS_RESOL
USE TPMALD_GEO      ,ONLY : GALD, ALDGEO_RESOL

USE ESET_RESOL_MOD  ,ONLY : ESET_RESOL
USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS
USE EDEALLOC_RESOL_MOD   ,ONLY : EDEALLOC_RESOL

IMPLICIT NONE

CHARACTER*5, OPTIONAL,  INTENT(IN) :: CDMODE
! Local variables
CHARACTER*5 :: CLMODE
INTEGER(KIND=JPIM) :: JRES
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ETRANS_END',0,ZHOOK_HANDLE)
CLMODE='FINAL'
IF (PRESENT(CDMODE)) CLMODE=CDMODE
IF (CLMODE == 'FINAL') THEN
 DO JRES=1,NDEF_RESOL
   CALL EDEALLOC_RESOL(JRES)
 ENDDO
 NULLIFY(R)
 IF (ALLOCATED(DIM_RESOL)) DEALLOCATE(DIM_RESOL)
 NULLIFY(RALD)
 IF (ALLOCATED(ALDDIM_RESOL)) DEALLOCATE(ALDDIM_RESOL)
!EQ_REGIONS
 IF (ASSOCIATED(N_REGIONS)) THEN
   DEALLOCATE(N_REGIONS)
   NULLIFY (N_REGIONS)
 ENDIF
!TPM_DISTR
 NULLIFY(D)
 IF (ALLOCATED(DISTR_RESOL)) DEALLOCATE(DISTR_RESOL)
 NULLIFY(DALD)
 IF (ALLOCATED(ALDDISTR_RESOL)) DEALLOCATE(ALDDISTR_RESOL)
#ifdef WITH_FFT992
!TPM_FFT
 NULLIFY(T)
 IF (ALLOCATED(FFT_RESOL)) DEALLOCATE(FFT_RESOL)
#endif
 !TPM_FFTW
 NULLIFY(TW)
 DEALLOCATE(FFTW_RESOL)
!TPM_FLT
 NULLIFY(S)
 IF (ALLOCATED(FLT_RESOL)) DEALLOCATE(FLT_RESOL)
#ifdef WITH_FFT992
 NULLIFY(TALD)
 IF (ALLOCATED(ALDFFT_RESOL)) DEALLOCATE(ALDFFT_RESOL)
#endif

!TPM_FIELDS
 NULLIFY(F)
 IF (ALLOCATED(FIELDS_RESOL)) DEALLOCATE(FIELDS_RESOL)
 NULLIFY(FALD)
 IF (ALLOCATED(ALDFIELDS_RESOL)) DEALLOCATE(ALDFIELDS_RESOL)

!TPM_GEOMETRY
 NULLIFY(G)
 IF(ALLOCATED(GEOM_RESOL)) DEALLOCATE(GEOM_RESOL)
 NULLIFY(GALD)
 IF(ALLOCATED(ALDGEO_RESOL)) DEALLOCATE(ALDGEO_RESOL)
!TPM_TRANS
 IF(ALLOCATED(FOUBUF_IN)) DEALLOCATE(FOUBUF_IN)
 IF(ALLOCATED(FOUBUF)) DEALLOCATE(FOUBUF)

 IF (ALLOCATED(LENABLED)) DEALLOCATE(LENABLED)
 MSETUP0 = 0
 NMAX_RESOL = 0
 NCUR_RESOL = 0
 NDEF_RESOL = 0
ENDIF

IF (CLMODE == 'FINAL' .OR. CLMODE == 'INTER') THEN
  !EQ_REGIONS
  IF (ASSOCIATED(N_REGIONS)) THEN
    DEALLOCATE(N_REGIONS)
    NULLIFY (N_REGIONS)
  ENDIF
 !TPM_DISTR
  IF (ALLOCATED(NPRCIDS)) DEALLOCATE(NPRCIDS)
ENDIF
IF (LHOOK) CALL DR_HOOK('ETRANS_END',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE ETRANS_END

