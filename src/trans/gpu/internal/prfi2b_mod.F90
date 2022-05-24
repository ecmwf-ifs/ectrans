! (C) Copyright 1990- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI2B_MOD
  CONTAINS
  SUBROUTINE PRFI2B(KF_FS,PAIA,KMODE)
  
  !**** *PRFI2B* - Prepare input work arrays for direct transform
  
  !     Purpose.
  !     --------
  !        To extract the Fourier fields for a specific zonal wavenumber
  !        and put them in an order suitable for the direct Legendre
  !        tranforms, i.e. split into symmetric and anti-symmetric part.
  
  !**   Interface.
  !     ----------
  !     *CALL* *PRFI2B(..)
  
  !        Explicit arguments :
  !        -------------------   KF_FS - number of fields
  !                              KM - zonal wavenumber
  !                              KMLOC - local zonal wavenumber
  !                              PAOA - antisymmetric part of Fourier
  !                              fields for zonal wavenumber KM
  !                              PSOA - symmetric part of Fourier
  !                              fields for zonal wavenumber KM
  
  !        Implicit arguments :  FOUBUF in TPM_TRANS
  !        --------------------
  
  !     Method.
  !     -------
  
  !     Externals.   None.
  !     ----------
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !        Mats Hamrud and Philippe Courtier  *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 90-07-01
  !        MPP Group: 95-10-01 Support for Distributed Memory version
  !        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
  !     ------------------------------------------------------------------
  
  USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT
  
  USE TPM_DIM         ,ONLY : R, R_NDGNH, R_NDGL
  USE TPM_TRANS       ,ONLY : FOUBUF
  USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
  USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1,MYPROC
  !
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),INTENT(IN)  :: KF_FS
  INTEGER(KIND=JPIM),INTENT(IN)  :: KMODE
  INTEGER(KIND=JPIM)  :: KM,KMLOC
  REAL(KIND=JPRBT)  , INTENT(OUT) :: PAIA(:,:,:)
!!  REAL(KIND=JPRBT)  , INTENT(OUT) :: PSIA(:,:,:),   PAIA(:,:,:)
  
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: IGLS,  ISL, JF, JGL, iunit
  
  INTEGER(KIND=JPIM) :: OFFSET1, OFFSET2
  
  
  !     ------------------------------------------------------------------
  
  !*       1.    EXTRACT SYM./ANTISYM. FIELDS FROM FOURIER ARRAY.
  !              ------------------------------------------------
 

!$ACC DATA PRESENT(PAIA,FOUBUF, D_NPNTGTB1,D_NSTAGT1B,D_MYMS,R_NDGL,R_NDGNH,G_NDGLU,D_NPROCL)

!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(2) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2,JGL)
DO KMLOC=1,D_NUMP
  DO JF=1,KF_FS*2
    KM = D_MYMS(KMLOC)
    ISL = R_NDGNH-G_NDGLU(KM)+1
    !$ACC LOOP SEQ
    DO JGL=ISL,R_NDGNH
      IGLS = R_NDGL+1-JGL
      OFFSET1 = (D_NSTAGT1B(D_NPROCL(JGL) )+D_NPNTGTB1(KMLOC,JGL ))*2*KF_FS
      OFFSET2 = (D_NSTAGT1B(D_NPROCL(IGLS))+D_NPNTGTB1(KMLOC,IGLS))*2*KF_FS
      IF( KMODE == -1 ) THEN
        PAIA(JF,JGL,KMLOC) = FOUBUF(OFFSET1+JF)-FOUBUF(OFFSET2+JF)
      ELSE
        PAIA(JF,JGL,KMLOC) = FOUBUF(OFFSET1+JF)+FOUBUF(OFFSET2+JF)
      ENDIF
    ENDDO
  ENDDO
END DO

!$ACC END DATA
  
  !     ------------------------------------------------------------------
  
  END SUBROUTINE PRFI2B
  END MODULE PRFI2B_MOD  
