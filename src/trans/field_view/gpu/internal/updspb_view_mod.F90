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

MODULE UPDSPB_VIEW_MOD
  CONTAINS
  SUBROUTINE UPDSPB_VIEW(POA,YDSPEC)
  
  
  !**** *UPDSPB_VIEW* - Update spectral arrays after direct Legendre transform
  
  !     Purpose.
  !     --------
  !        To update spectral arrays for a fixed zonal wave-number
  !         from values in POA.
  
  !**   Interface.
  !     ----------
  !        CALL UPDSPB_VIEW(....)
  
  !        Explicit arguments :  
  !        --------------------  
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
  USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
  USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW
  !
  
  IMPLICIT NONE
    
  REAL(KIND=JPRBT)  ,INTENT(IN)  :: POA(:,:,:)
  TYPE(SPEC_VIEW) :: YDSPEC(:)

  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: KM,KMLOC
  INTEGER(KIND=JPIM) :: INM, IR, JFLD, JN, IASM0
  INTEGER(KIND=JPIM) :: IFIELD
  
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

  !*       1.    UPDATE SPECTRAL FIELDS.
  !              -----------------------

 IFIELD = SIZE(YDSPEC)

#ifdef OMPGPU
  !$OMP TARGET DATA MAP(PRESENT,ALLOC:POA,R,R_NTMAX,D,D_NUMP,D_MYMS,D_NASM0,YDSPEC)
#endif
#ifdef ACCGPU
  !$ACC DATA PRESENT(POA,R,R_NTMAX,D,D_NUMP,D_MYMS,D_NASM0) COPYIN(YDSPEC) ASYNC(1)
#endif

! Directive incomplete -> putting more variables in SHARED() triggers internal compiler error
! ftn-7991: INTERNAL COMPILER ERROR:  "Too few arguments on the stack"
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(KM,IASM0,INM) &
  !$OMP& SHARED(D,R,IFIELD,POA,YDSPEC) MAP(TO:IFIELD)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,IASM0,INM) DEFAULT(NONE) COPYIN(IFIELD) &
#ifndef _CRAYFTN
  !$ACC& ASYNC(1)
#else
  !$ACC&
#endif
#endif
  DO KMLOC=1,D_NUMP
    DO JN=3,R_NTMAX+3
      DO JFLD=1,IFIELD
        KM = D_MYMS(KMLOC)
        IASM0 = D_NASM0(KM)

        IF(KM /= 0 .AND. JN <= R_NTMAX+3-KM) THEN
        !(DO JN=3,R_NTMAX+3-KM)
          INM = IASM0+((R_NTMAX+3-JN)-KM)*2
          YDSPEC(JFLD)%P(INM)   = POA(2*JFLD-1,JN,KMLOC)
          YDSPEC(JFLD)%P(INM+1) = POA(2*JFLD  ,JN,KMLOC)
        ELSEIF (KM == 0) THEN
          !(DO JN=3,R_NTMAX+3)
          INM = IASM0+(R_NTMAX+3-JN)*2
          YDSPEC(JFLD)%P(INM)   = POA(2*JFLD-1,JN,KMLOC)
          YDSPEC(JFLD)%P(INM+1) = 0.0_JPRBT
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
 
  END SUBROUTINE UPDSPB_VIEW
END MODULE UPDSPB_VIEW_MOD
