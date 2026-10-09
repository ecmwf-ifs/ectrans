! (C) Copyright 1988- ECMWF.
! (C) Copyright 1988- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSP_VIEW_MOD
CONTAINS
SUBROUTINE UPDSP_VIEW(KF_UV,KF_SCALARS,POA1, YDSPVSCALAR)

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
!        YDSPVSCALAR - scalar field

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.  UPDSPB_VIEW - basic transfer routine
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
!        Modified : 94-08-02 R. El Khatib - interface to UPDSPB_VIEW
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS, ONLY: JPIM ,JPRB, JPRBT
USE TPM_TRANS,       ONLY: NF_SC2, NF_SC3A, NF_SC3B
USE TPM_DISTR,       ONLY: D
USE UPDSPB_VIEW_MOD,      ONLY: UPDSPB_VIEW
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_UV,KF_SCALARS
REAL(KIND=JPRBT) , INTENT(IN)  :: POA1(:,:,:)
TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IST ,IEND,IDIM1,IDIM3,J3


!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------

!*       1.1      VORTICITY AND DIVERGENCE.



IST = 1
IST = IST+4*KF_UV


!*       1.2   SCALARS

IF (KF_SCALARS > 0) THEN
  IEND = IST+2*KF_SCALARS-1
  CALL UPDSPB_VIEW(POA1(IST:IEND,:,:),YDSPVSCALAR)
ENDIF


!     ------------------------------------------------------------------

END SUBROUTINE UPDSP_VIEW
END MODULE UPDSP_VIEW_MOD
