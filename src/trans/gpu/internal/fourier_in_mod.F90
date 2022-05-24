! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_IN_MOD
CONTAINS
SUBROUTINE FOURIER_IN(PREEL,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NSTAGTF,D_MSTABF,D_NSTAGT0B,D_NPNTGTB0,D_NPROCM,D_NPTRLS
USE TPM_TRANS       ,ONLY : FOUBUF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
use tpm_gen, only: nout
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS

INTEGER(KIND=JPIM) :: KGL

REAL(KIND=JPRBT), INTENT(OUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,iimax1,iimax2,iimax3,iunit

!     ------------------------------------------------------------------

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

!$ACC DATA PRESENT(D_NPTRLS,G_NMEN,D_NPROCM,D_NSTAGT0B,D_MSTABF,D_NPNTGTB0,FOUBUF,PREEL,D_NSTAGTF)
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IGLG,IPROC,ISTA) DEFAULT(NONE)
DO KGL=IBEG,IEND,IINC
   DO JM=0,G_NMEN_MAX      
      DO JF=1,KFIELDS     
         
         IGLG = D_NPTRLS(MYSETW)+KGL-1
         
         if ( JM .le. G_NMEN(IGLG)) then
            
            IPROC = D_NPROCM(JM)
            ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*2*KFIELDS
            
            PREEL(2*JF-1,2*JM+1+D_NSTAGTF(KGL)) = FOUBUF(ISTA+2*JF-1)
            PREEL(2*JF,  2*JM+1+D_NSTAGTF(KGL)) = FOUBUF(ISTA+2*JF  )
            !write(nout,*) , 'istart1 ...', KGL, JM, JF, ISTA+2*JF,ISTA,D_NSTAGT0B(D_MSTABF(IPROC)),IPROC,KFIELDS
            !write(nout,*) , 'istart2 ...',D_NPNTGTB0(JM,KGL), FOUBUF(ISTA+2*JF-1), FOUBUF(ISTA+2*JF),2*JM+1+D_NSTAGTF(KGL)
            !if(jf==1 .and. 2*JM+1+D_NSTAGTF(KGL)==7972) write(nout,*) 'fourier_in: fidx=7972, kgl=',kgl,' jm=',jm
            !TODO (Andreas): should be able to remove the factor 2 in the second dimension (in front of jm)
            !and reduce the size of the array. Will need to adapt fsc_mod accordingly! This is actually more
            !difficult: d_nstagtf(kgl) is not necessarily even!
            
         end if
      ENDDO
   ENDDO
ENDDO
!$ACC END DATA

!iimax1=0
!iimax2=0
!iimax3=0
!iunit=myproc+300
!DO KGL=IBEG,IEND,IINC
!   DO JM=0,G_NMEN_MAX      
!      DO JF=1,KFIELDS     
!         
!         IGLG = D_NPTRLS(MYSETW)+KGL-1
!         
!         if ( JM .le. G_NMEN(IGLG)) then
!            
!            IPROC = D_NPROCM(JM)
!            ISTA  = (D_NSTAGT0B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,KGL))*2*KFIELDS
!           
!            iimax1=max(iimax1,2*JF)
!            iimax2=max(iimax2,2*JM+1+D_NSTAGTF(KGL))
!            iimax3=max(iimax3,ISTA+2*JF)
!
!          endif
!      ENDDO 
!    ENDDO 
!ENDDO 
!write(iunit,*) 'max_in ',iimax1,size(PREEL,1),iimax2,size(PREEL,2),iimax3,size(FOUBUF)
!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_IN
END MODULE FOURIER_IN_MOD

