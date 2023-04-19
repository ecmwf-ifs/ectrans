! (C) Copyright 1990- ECMWF.
! (C) Copyright 1990- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI2B_MOD
  CONTAINS
  SUBROUTINE PRFI2B(KFIELD,PAIA,KMODE)
  
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
  !        -------------------   KFIELD - number of fields
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
  
  USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRBT
  
  USE TPM_DIM         ,ONLY : R, R_NDGNH, R_NDGL
  USE TPM_TRANS       ,ONLY : FOUBUF
  USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
  USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1,MYPROC
  !
 
  IMPLICIT NONE
 
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
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


#ifdef ACCGPU
!$ACC DATA PRESENT(PAIA,FOUBUF, D_NPNTGTB1,D_NSTAGT1B,D_MYMS,R_NDGL,R_NDGNH,G_NDGLU,D_NPROCL)
#endif
#ifdef OMPGPU
!WARNING: following line should be PRESENT,ALLOC but causes issues with AMD compiler!
!$OMP TARGET DATA MAP(ALLOC:PAIA,FOUBUF, D_NPNTGTB1,D_NSTAGT1B,D_MYMS,R_NDGL,R_NDGNH,G_NDGLU,D_NPROCL)
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2) &
!$OMP&    SHARED(D_NUMP,R_NDGNH,KFIELD,D_MYMS,G_NDGLU,R_NDGL,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1,KMODE,PAIA,FOUBUF)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2) &
!$ACC&    COPYIN(KFIELD,KMODE) &
!$ACC&    PRESENT(D_NUMP,R_NDGNH,D_MYMS,G_NDGLU,R_NDGL,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1,PAIA,FOUBUF)
#endif
DO KMLOC=1,D_NUMP
     DO JGL=1,R_NDGNH
        DO JF=1,KFIELD*2
           KM = D_MYMS(KMLOC)
           ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
           if (JGL .ge. ISL) then
              IGLS = R_NDGL+1-JGL
              OFFSET1 = (D_NSTAGT1B(D_NPROCL(JGL) )+D_NPNTGTB1(KMLOC,JGL ))*2*KFIELD
              OFFSET2 = (D_NSTAGT1B(D_NPROCL(IGLS))+D_NPNTGTB1(KMLOC,IGLS))*2*KFIELD
              IF( KMODE == -1 ) THEN
                PAIA(JF,JGL,KMLOC) = FOUBUF(OFFSET1+JF)-FOUBUF(OFFSET2+JF)
              ELSE
                PAIA(JF,JGL,KMLOC) = FOUBUF(OFFSET1+JF)+FOUBUF(OFFSET2+JF)
!                PSIA(JF,JGL,KMLOC) = FOUBUF(OFFSET1+JF)+FOUBUF(OFFSET2+JF)
              ENDIF
           end if
        ENDDO
     ENDDO
END DO

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif
 
  !     ------------------------------------------------------------------
 
  END SUBROUTINE PRFI2B
END MODULE PRFI2B_MOD  
