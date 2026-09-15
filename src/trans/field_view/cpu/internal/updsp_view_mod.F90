! (C) Copyright 1988- ECMWF.
! (C) Copyright 1988- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSP_VIEW_MOD
CONTAINS
SUBROUTINE UPDSP_VIEW(KM,POA1,POA2, YDSPVVOR, YDSPVDIV, YDSPVSCALAR)

!**** *UPDSP_VIEW* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update the spectral arrays for a fixed zonal wave-number
!        from values in POA1 and POA2.

!**   Interface.
!     ----------
!        CALL UPDSP_VIEW(...)

!        Explicit arguments :
!        --------------------
!        KM - zonal wave-number
!        POA1 - spectral fields for zonal wavenumber KM (basic var.)
!        POA2 - spectral fields for zonal wavenumber KM (vor. div.)
!        PSPVOR - spectral vorticity
!        PSPDIV - spectral divergence
!        PSPSCALAR - spectral scalar variables

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.  UPDSP_VIEWB - basic transfer routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-02
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 94-08-02 R. El Khatib - interface to UPDSP_VIEWB
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : NF_SC2, NF_SC3A, NF_SC3B
USE TPM_DISTR       ,ONLY : D

USE UPDSPB_VIEW_MOD      ,ONLY : UPDSPB_VIEW
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
REAL(KIND=JPRB) , INTENT(IN)  :: POA1(:,:)
REAL(KIND=JPRB) , INTENT(IN)  :: POA2(:,:)
TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IST, IEND, IVORS, IVORE, IDIVS, IDIVE, JFLD
INTEGER(KIND=JPIM) :: IF_UV, IF_SCALARS

!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------


!*       1.1      VORTICITY AND DIVERGENCE.

IF_UV = SIZE(YDSPVVOR)

IST = 1
IF (IF_UV > 0) THEN
  IST = IST+4*IF_UV
  IVORS = 1
  IVORE = 2*IF_UV
  IDIVS = 2*IF_UV+1
  IDIVE = 4*IF_UV
  CALL UPDSPB_VIEW(KM,POA2(:,IVORS:IVORE),YDSPVVOR)
  CALL UPDSPB_VIEW(KM,POA2(:,IDIVS:IDIVE),YDSPVDIV)
  IF (KM == 0) THEN
    DO JFLD=1,IF_UV
      YDSPVVOR(JFLD)%P(D%NASM0(0)) = 0.0_JPRB
      YDSPVDIV(JFLD)%P(D%NASM0(0)) = 0.0_JPRB
    ENDDO
  ENDIF
ENDIF

!*       1.2   SCALARS

IF_SCALARS = SIZE(YDSPVSCALAR)
IF (IF_SCALARS > 0) THEN
    IEND = IST+2*IF_SCALARS-1
    CALL UPDSPB_VIEW(KM,POA1(:,IST:IEND),YDSPVSCALAR)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UPDSP_VIEW
END MODULE UPDSP_VIEW_MOD
