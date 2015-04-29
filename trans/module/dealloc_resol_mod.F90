MODULE DEALLOC_RESOL_MOD
CONTAINS
SUBROUTINE DEALLOC_RESOL(KRESOL)

!**** *DEALLOC_RESOL* - Deallocations of a resolution

!     Purpose.
!     --------
!     Release allocated arrays for a given resolution

!**   Interface.
!     ----------
!     CALL DEALLOC_RESOL

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
!        Original : 09-Jul-2013 from trans_end

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : LENABLED, NOUT
USE TPM_DISTR       ,ONLY : D
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FIELDS      ,ONLY : F
USE TPM_FFT         ,ONLY : T
USE TPM_FFTW        ,ONLY : TW,DESTROY_PLANS_FFTW
USE TPM_FLT         ,ONLY : S

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),  INTENT(IN) :: KRESOL

!     ------------------------------------------------------------------

IF (.NOT.LENABLED(KRESOL)) THEN

  WRITE(UNIT=NOUT,FMT='('' DEALLOC_RESOL WARNING : KRESOL = '',I3,'' ALREADY DISABLED '')') KRESOL

ELSE

  CALL SET_RESOL(KRESOL)

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

    !TPM_FFT
     IF (.NOT.D%LCPNMONLY) THEN
       IF( ASSOCIATED(T) ) THEN
         IF( ALLOCATED(T%TRIGS) ) DEALLOCATE(T%TRIGS)
         IF( ALLOCATED(T%NFAX) )  DEALLOCATE(T%NFAX)
        ENDIF
     ENDIF
    !TPM_FFTW

    IF( TW%LFFTW )THEN
      CALL DESTROY_PLANS_FFTW
    ENDIF


    !TPM_FIELDS
    DEALLOCATE(F%RMU,F%RW,F%R1MU2,F%RACTHE)
    DEALLOCATE(F%REPSNM,F%RN,F%RLAPIN,F%NLTN)
    IF( S%LKEEPRPNM ) THEN
      DEALLOCATE(F%RPNM)
    ENDIF

    !TPM_GEOMETRY
    DEALLOCATE(G%NMEN,G%NDGLU)

  ELSE

    DEALLOCATE(G%NLOEN)

  ENDIF

  LENABLED(KRESOL)=.FALSE.

ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE DEALLOC_RESOL
END MODULE DEALLOC_RESOL_MOD
