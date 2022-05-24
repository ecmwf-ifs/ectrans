! (C) Copyright 1988- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSPB_MOD
  CONTAINS
  SUBROUTINE UPDSPB(KFIELD,POA,PSPEC,KFLDPTR)
  
  
  !**** *UPDSPB* - Update spectral arrays after direct Legendre transform
  
  !     Purpose.
  !     --------
  !        To update spectral arrays for a fixed zonal wave-number
  !         from values in POA.
  
  !**   Interface.
  !     ----------
  !        CALL UPDSPB(....)
  
  !        Explicit arguments :  KM - zonal wavenumber
  !        --------------------  KFIELD  - number of fields
  !                              POA - work array
  !                              PSPEC - spectral array
  
  !        Implicit arguments :  None
  !        --------------------
  
  !     Method.
  !     -------
  
  !     Externals.
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
  !        D. Giard : 93-03-19 truncations NSMAX and NTMAX (see NOTE)
  !        R. El Khatib : 94-08-02 Replace number of fields by indexes of the
  !                       first and last field
  !        L. Isaksen : 95-06-06 Reordering of spectral arrays
  !     ------------------------------------------------------------------
  
USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT
  
  USE TPM_DIM         ,ONLY : R,R_NSMAX,R_NTMAX
  !USE TPM_FIELDS
  USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NASM0
  !
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
  INTEGER(KIND=JPIM)  :: KM,KMLOC
  REAL(KIND=JPRBT)   ,INTENT(IN)  :: POA(:,:,:)
  REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPEC(:,:)
  INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN, ISMAX, ITMAX, IASM0,IFLD
  
  
  !     ------------------------------------------------------------------
  
  !*       0.    NOTE.
  !              -----
  
  ! The following transfer reads :
  ! SPEC(k,NASM0(m)+NLTN(n)*2)  =POA(nn,2*k-1) (real part)
  ! SPEC(k,NASM0(m)+NLTN(n)*2+1)=POA(nn,2*k  ) (imaginary part)
  ! with n from m to NSMAX
  ! and nn=NTMAX+2-n from NTMAX+2-m to NTMAX+2-NSMAX.
  ! NLTN(m)=NTMAX+2-m : n=NLTN(nn),nn=NLTN(n)
  ! nn is the loop index.
  
    IF(PRESENT(KFLDPTR)) THEN
       stop 'Error: code path not (yet) supported in GPU version'
    ENDIF
  
  !*       1.    UPDATE SPECTRAL FIELDS.
  !              -----------------------

  !loop over wavenumber
  !$ACC DATA PRESENT(PSPEC,POA,R,D)
  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IASM0,INM,IR,II) DEFAULT(NONE)
  DO KMLOC=1,D%NUMP
    DO JFLD=1,KFIELD
      KM = D%MYMS(KMLOC)
      IASM0 = D%NASM0(KM)

      IF(KM == 0) THEN
        !$ACC LOOP SEQ
        DO JN=3,R%NTMAX+3
          INM = IASM0+(R%NTMAX+3-JN)*2
          PSPEC(JFLD,INM)   = POA(2*JFLD-1,JN,KMLOC)
          PSPEC(JFLD,INM+1) = 0.0_JPRBT
        ENDDO
      ELSE
        !$ACC LOOP SEQ
        DO JN=3,R%NTMAX+3-KM
          INM = IASM0+((R%NTMAX+3-JN)-KM)*2
          PSPEC(JFLD,INM)   = POA(2*JFLD-1,JN,KMLOC)
          PSPEC(JFLD,INM+1) = POA(2*JFLD  ,JN,KMLOC)
        ENDDO
      END IF
    ENDDO
  ENDDO
!$ACC END PARALLEL
!$ACC END DATA
    
  !     ------------------------------------------------------------------
  
  END SUBROUTINE UPDSPB
  END MODULE UPDSPB_MOD
