MODULE EPRFI2AD_MOD
CONTAINS
SUBROUTINE EPRFI2AD(KM,KMLOC,KF_FS,PFFT)

!**** *EPRFI2AD* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!        *CALL* *EPRFI2AD(..)

!        Explicit arguments :
!        --------------------  KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAIA - antisymmetric part of Fourier
!                              components for KM  (output)
!                              PSIA - symmetric part of Fourier
!                              components for KM  (output)

!        Implicit arguments :  The Grid point arrays of the model.
!        --------------------

!     Method.
!     -------

!     Externals.   EPRFI2BAD   - basic copying routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-25
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - CDCONF='T'
!        Modified : 93-??-?? ???????? - CDCONF='C'--> bug if CDCONF='T'
!        Modified : 93-05-13 D. Giard - correction of the previous bug
!        Modified : 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified : 94-04-06 R. El Khatib - Full-POS configuration 'P'
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE EPRFI2BAD_MOD   ,ONLY : EPRFI2BAD
!

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(IN) :: KM
INTEGER(KIND=JPIM) , INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM) , INTENT(IN) :: KF_FS

REAL(KIND=JPRB) , INTENT(IN) :: PFFT(:,:)

!     ------------------------------------------------------------------

!*       2.    EXTRACT SYM./ANTISYM. FIELDS FROM TIME T+1.
!              -------------------------------------------

CALL EPRFI2BAD(KF_FS,KM,KMLOC,PFFT)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI2AD
END MODULE EPRFI2AD_MOD
