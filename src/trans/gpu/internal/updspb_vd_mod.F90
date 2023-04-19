! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSPB_VD_MOD
CONTAINS
SUBROUTINE UPDSPB_VD(KFIELD,PSPVOR,PSPDIV,KFLDPTR)
  
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
 
  USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRB,  JPRBT

  USE TPM_DIM         ,ONLY : R,R_NSMAX,R_NTMAX
  USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NASM0
  USE TPM_FIELDS      ,ONLY : ZOA2
  !
 
  IMPLICIT NONE
 
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
  INTEGER(KIND=JPIM)  :: KM,KMLOC
  REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPVOR(:,:)
  REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPDIV(:,:)
  INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN, ISMAX, ITMAX, IASM0,IFLD
  INTEGER(KIND=JPIM) :: IVORS, IDIVS


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
  
  IVORS = 1
  IDIVS = 2*KFIELD+1

  !*       1.    UPDATE SPECTRAL FIELDS.
  !              -----------------------
#ifdef ACCGPU
  !$ACC DATA &
  !$ACC& PRESENT(ZOA2) &
  !$ACC& COPY(PSPVOR,PSPDIV) &
  !$ACC& COPY(D,D_NUMP,D_MYMS,R,R_NSMAX,R_NTMAX,D,D_NASM0)
#endif
#ifdef OMPGPU
  !$OMP TARGET DATA &
  !$OMP& MAP(PRESENT,ALLOC:ZOA2) &
  !$OMP& MAP(FROM:PSPVOR,PSPDIV) &
  !$OMP& MAP(TO:D,D_NUMP,D_MYMS,R,R_NSMAX,R_NTMAX,D,D_NASM0)
#endif
 
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,INM,IR,II,IASM0,IFLD)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,INM,IR,II,IASM0,IFLD)
#endif
  DO KMLOC=1,D_NUMP     
       DO JN=R_NTMAX+2-R_NSMAX,R_NTMAX+2
          DO JFLD=1,KFIELD
 
             KM = D_MYMS(KMLOC)
             IASM0 = D_NASM0(KM)
 
             IF(KM == 0) THEN

                IF (JN .LE. R_NTMAX+2-KM) THEN
                   INM = IASM0+(R_NTMAX+2-JN)*2
                   IR = 2*JFLD-1
                   PSPVOR(JFLD,INM)   = ZOA2(IVORS+IR-1,JN,KMLOC)
                   PSPDIV(JFLD,INM)   = ZOA2(IDIVS+IR-1,JN,KMLOC)
                   PSPVOR(JFLD,INM+1) = 0.0_JPRBT
                   PSPDIV(JFLD,INM+1) = 0.0_JPRBT
                END IF
                IF(PRESENT(KFLDPTR)) THEN
                  IFLD = KFLDPTR(JFLD)
                  PSPVOR(IFLD,IASM0) = 0.0_JPRBT
                  PSPDIV(IFLD,IASM0) = 0.0_JPRBT
                ELSE
                  PSPVOR(JFLD,IASM0) = 0.0_JPRBT
                  PSPDIV(JFLD,IASM0) = 0.0_JPRBT
                ENDIF

             ELSE
 
 
                IF (JN .LE. R_NTMAX+2-KM) THEN
                   INM = IASM0+((R_NTMAX+2-JN)-KM)*2
 
                   IR = 2*JFLD-1
                   II = IR+1
                   PSPVOR(JFLD,INM)   = ZOA2(IVORS+IR-1,JN,KMLOC)
                   PSPVOR(JFLD,INM+1) = ZOA2(IVORS+II-1,JN,KMLOC)
                   PSPDIV(JFLD,INM)   = ZOA2(IDIVS+IR-1,JN,KMLOC)
                   PSPDIV(JFLD,INM+1) = ZOA2(IDIVS+II-1,JN,KMLOC)
 
                END IF
             END IF
 
          ENDDO
 
       ENDDO
       !end loop over wavenumber
   END DO
#ifdef OMPGPU
   !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
   !$ACC end data
#endif
 
  !     ------------------------------------------------------------------
 
  END SUBROUTINE UPDSPB_VD
END MODULE UPDSPB_VD_MOD
