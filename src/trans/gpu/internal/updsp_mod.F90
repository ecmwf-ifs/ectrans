! (C) Copyright 1988- ECMWF.
! (C) Copyright 1988- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSP_MOD
CONTAINS
SUBROUTINE UPDSP(KF_UV,KF_SCALARS,POA1,POA2, &
 &               PSPVOR,PSPDIV,PSPSCALAR,&
 &               PSPSC3A,PSPSC3B,PSPSC2 , &
 &               KFLDPTRUV,KFLDPTRSC)

!**** *UPDSP* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update the spectral arrays for a fixed zonal wave-number
!        from values in POA1 and POA2.

!**   Interface.
!     ----------
!        CALL UPDSP(...)

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

!     Externals.  UPDSPB - basic transfer routine
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
!        Modified : 94-08-02 R. El Khatib - interface to UPDSPB
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRB,  JPRBT

USE TPM_TRANS       ,ONLY : NF_SC2, NF_SC3A, NF_SC3B
USE TPM_DISTR       ,ONLY : D

USE UPDSPB_MOD      ,ONLY : UPDSPB

IMPLICIT NONE


!     DUMMY INTEGER SCALARS

INTEGER(KIND=JPIM), INTENT(IN)  :: KF_UV,KF_SCALARS

REAL(KIND=JPRBT) , INTENT(IN)  :: POA1(:,:,:)
REAL(KIND=JPRBT) , INTENT(IN)  :: POA2(:,:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IVORS, IVORE, IDIVS, IDIVE, IST ,IEND,JFLD,IFLD,IDIM1,IDIM3,J3


!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------

!*       1.1      VORTICITY AND DIVERGENCE.

#ifdef ACCGPU
!$ACC DATA PRESENT(PSPVOR,PSPDIV) IF(KF_UV > 0)
!$ACC DATA PRESENT(PSPSCALAR) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSCALAR))
!$ACC DATA PRESENT(PSPSC2)  IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC2) .AND. NF_SC2 > 0)
!$ACC DATA PRESENT(PSPSC3A) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC3A) .AND. NF_SC3A > 0)
!$ACC DATA PRESENT(PSPSC3B) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC3B) .AND. NF_SC3B > 0)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:PSPVOR,PSPDIV) IF(KF_UV > 0)
!$OMP TARGET DATA MAP(ALLOC:PSPSCALAR) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSCALAR))
!$OMP TARGET DATA MAP(ALLOC:PSPSC2)  IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC2) .AND. NF_SC2 > 0)
!$OMP TARGET DATA MAP(ALLOC:PSPSC3A) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC3A) .AND. NF_SC3A > 0)
!$OMP TARGET DATA MAP(ALLOC:PSPSC3B) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC3B) .AND. NF_SC3B > 0)
#endif

IST = 1
IF (KF_UV > 0) THEN
  !stop 'Error: code path not (yet) supported in GPU version'

   IST = IST+4*KF_UV
   IVORS = 1
   IVORE = 2*KF_UV
   IDIVS = 2*KF_UV+1
   IDIVE = 4*KF_UV
   CALL UPDSPB(KF_UV,IVORS,POA2,PSPVOR,KFLDPTRUV)
   CALL UPDSPB(KF_UV,IDIVS,POA2,PSPDIV,KFLDPTRUV)
ENDIF

!*       1.2   SCALARS

IF (KF_SCALARS > 0) THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IEND = IST+2*KF_SCALARS-1
    CALL UPDSPB(KF_SCALARS,IST,POA1,PSPSCALAR,KFLDPTRSC)
  ELSE
    IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
      IDIM1 = NF_SC2
      IEND  = IST+2*IDIM1-1
      CALL UPDSPB(IDIM1,IST,POA1,PSPSC2)
      IST=IST+2*IDIM1
    ENDIF
    IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
      IDIM1=NF_SC3A
      IDIM3=UBOUND(PSPSC3A,3)
      DO J3=1,IDIM3
        IEND = IST+2*IDIM1-1
        CALL UPDSPB(IDIM1,IST,POA1,PSPSC3A(:,:,J3))
        IST=IST+2*IDIM1
      ENDDO
    ENDIF
    IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
      IDIM1=NF_SC3B
      IDIM3=UBOUND(PSPSC3B,3)
      DO J3=1,IDIM3
        IEND = IST+2*IDIM1-1
        CALL UPDSPB(IDIM1,IST,POA1,PSPSC3B(:,:,J3))
        IST=IST+2*IDIM1
      ENDDO
    ENDIF
  ENDIF
ENDIF

#ifdef OMPGPU
!$OMP END TARGET DATA
!$OMP END TARGET DATA
!$OMP END TARGET DATA
!$OMP END TARGET DATA
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
#endif

!     ------------------------------------------------------------------

END SUBROUTINE UPDSP
END MODULE UPDSP_MOD
