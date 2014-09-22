SUBROUTINE TRANS_END(CDMODE)

!**** *TRANS_END* - Terminate transform package

!     Purpose.
!     --------
!     Terminate transform package. Release all allocated arrays.

!**   Interface.
!     ----------
!     CALL TRANS_END

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
!          G. Radnoti: 19-03-2009: intermediate end of transf to allow to switch to mono-task transforms
!        R. El Khatib 09-Jul-2013 LENABLED

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : MSETUP0, NCUR_RESOL, NDEF_RESOL, NMAX_RESOL, LENABLED
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL, NPRCIDS
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : F, FIELDS_RESOL
USE TPM_FFT         ,ONLY : T, FFT_RESOL
USE TPM_FLT
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN

USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE DEALLOC_RESOL_MOD   ,ONLY : DEALLOC_RESOL
!

IMPLICIT NONE
CHARACTER*5, OPTIONAL,  INTENT(IN) :: CDMODE
! Local variables
INTEGER(KIND=JPIM) :: JRES
CHARACTER*5 :: CLMODE
!     ------------------------------------------------------------------
CLMODE='FINAL'
IF (PRESENT(CDMODE)) CLMODE=CDMODE
IF (CLMODE == 'FINAL') THEN
  DO JRES=1,NDEF_RESOL
    CALL DEALLOC_RESOL(JRES)
  ENDDO

  NULLIFY(R)
  DEALLOCATE(DIM_RESOL)

  NULLIFY(D)
  DEALLOCATE(DISTR_RESOL)

  !TPM_FFT
  NULLIFY(T)
  DEALLOCATE(FFT_RESOL)

  !TPM_FLT
  NULLIFY(S)
  DEALLOCATE(FLT_RESOL)

  !TPM_FIELDS
  NULLIFY(F)
  DEALLOCATE(FIELDS_RESOL)


  !TPM_GEOMETRY
  NULLIFY(G)
  DEALLOCATE(GEOM_RESOL)

  !TPM_TRANS
  IF(ALLOCATED(FOUBUF_IN)) DEALLOCATE(FOUBUF_IN)
  IF(ALLOCATED(FOUBUF)) DEALLOCATE(FOUBUF)

  DEALLOCATE(LENABLED)
  MSETUP0 = 0
  NMAX_RESOL = 0
  NCUR_RESOL = 0
  NDEF_RESOL = 0
ENDIF
IF (CLMODE == 'FINAL' .OR. CLMODE == 'INTER') THEN
  !EQ_REGIONS
  DEALLOCATE(N_REGIONS)
  !TPM_DISTR
  DEALLOCATE(NPRCIDS)
ENDIF

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE TRANS_END
