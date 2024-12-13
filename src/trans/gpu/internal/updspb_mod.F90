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
  
  USE PARKIND_ECTRANS, ONLY: JPIM, JPRB, JPRBT
  USE TPM_DIM,         ONLY: R
  USE TPM_DISTR,       ONLY: D
  !
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
  REAL(KIND=JPRBT)  ,INTENT(IN)  :: POA(:,:,:)
  REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPEC(:,:)
  INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM)  :: KM,KMLOC
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
  ASSOCIATE(D_NUMP=>D%NUMP, D_MYMS=>D%MYMS, D_NASM0=>D%NASM0, R_NTMAX=>R%NTMAX)

  IF(PRESENT(KFLDPTR)) THEN
      stop 'Error: code path not (yet) supported in GPU version'
  ENDIF
  
  !*       1.    UPDATE SPECTRAL FIELDS.
  !              -----------------------

#ifdef OMPGPU
  !$OMP TARGET DATA MAP(TO:KFIELD)
#endif
#ifdef ACCGPU
  !$ACC DATA PRESENT(PSPEC,POA,R,R_NTMAX,D,D_NUMP,D_MYMS,D_NASM0) ASYNC(1)
#endif

! Directive incomplete -> putting more variables in SHARED() triggers internal compiler error
! ftn-7991: INTERNAL COMPILER ERROR:  "Too few arguments on the stack"
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(KM,IASM0,INM) &
  !$OMP& SHARED(KFIELD)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,IASM0,INM) DEFAULT(NONE) COPYIN(KFIELD) &
#ifndef _CRAYFTN
  !$ACC& ASYNC(1)
#else
  !$ACC&
#endif
#endif
  DO KMLOC=1,D_NUMP
    DO JN=3,R_NTMAX+3
      DO JFLD=1,KFIELD
        KM = D_MYMS(KMLOC)
        IASM0 = D_NASM0(KM)

        IF(KM /= 0 .AND. JN <= R_NTMAX+3-KM) THEN
        !(DO JN=3,R_NTMAX+3-KM)
          INM = IASM0+((R_NTMAX+3-JN)-KM)*2
          PSPEC(JFLD,INM)   = POA(2*JFLD-1,JN,KMLOC)
          PSPEC(JFLD,INM+1) = POA(2*JFLD  ,JN,KMLOC)
        ELSEIF (KM == 0) THEN
          !(DO JN=3,R_NTMAX+3)
          INM = IASM0+(R_NTMAX+3-JN)*2
          PSPEC(JFLD,INM)   = POA(2*JFLD-1,JN,KMLOC)
          PSPEC(JFLD,INM+1) = 0.0_JPRBT
        END IF
      ENDDO
    ENDDO
  ENDDO

#ifdef ACCGPU
  !$ACC END DATA
#endif
#ifdef OMPGPU
  !$OMP END TARGET DATA
#endif
 
  END ASSOCIATE
  !     ------------------------------------------------------------------
 
  END SUBROUTINE UPDSPB
END MODULE UPDSPB_MOD
