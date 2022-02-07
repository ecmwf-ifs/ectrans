! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUT_MOD
CONTAINS
SUBROUTINE FOURIER_OUT(KFIELDS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

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

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT1B,D_NPNTGTB0,D_NPROCM, D_NPROCL
USE TPM_TRANS       ,ONLY : FOUBUF_IN, ZGTF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
!

IMPLICIT NONE

!REAL(KIND=JPRBT), INTENT(IN) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM) :: KGL

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA, ISTA1,JMMAX, iunit

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC, IOFF,iimax1,iimax2,iimax3

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

!! !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(KGL,JM,JF,IGLG,JMMAX,IPROC,ISTA,IOFF)
!$ACC DATA PRESENT(FOUBUF_IN,ZGTF, D_NPTRLS,G_NMEN,D_NPROCM,D_NSTAGT1B,D_MSTABF,D_NPNTGTB0,D_NSTAGTF) COPYIN(IBEG,IEND,IINC)
!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(IGLG,JMMAX,IPROC,ISTA,IOFF)
DO KGL=IBEG,IEND,IINC
   DO JM=0,G_NMEN_MAX      
      DO JF=1,KFIELDS
         
         IGLG = D_NPTRLS(MYSETW)+KGL-1
         JMMAX = G_NMEN(IGLG)
         if  (JM .le. JMMAX) then
            
            IPROC = D_NPROCM(JM)
            ISTA  = (D_NSTAGT1B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,KGL))*2*KFIELDS
            IOFF  = 1+D_NSTAGTF(KGL)
           
            ! imaginary may be not JM+1 but JM+G_NMEN(IGLG)+1
            FOUBUF_IN(ISTA+2*JF-1) = ZGTF(2*JF-1, 2*JM+IOFF)
            FOUBUF_IN(ISTA+2*JF  ) = ZGTF(2*JF  , 2*JM+IOFF)
            !if( myproc.eq.1 .and. jf.eq.1 .and. JM.eq.0 ) write(*,*) 'fou_oD ',ISTA+2*JF-1,FOUBUF_IN(ISTA+2*JF-1),FOUBUF_IN(ISTA+2*JF  )
            !if( myproc.eq.1 .and. jf.eq.1 .and. JM.eq.0 ) write(*,*) 'fou_oD1 ', ZGTF(JF, 2*JM+IOFF), ZGTF(JF, 2*JM+1+IOFF) 
            !if( myproc.eq.2 .and. jf.eq.1 .and. JM.eq.0 ) write(*,*) 'fou_oD ',ISTA+2*JF-1,FOUBUF_IN(ISTA+2*JF-1),FOUBUF_IN(ISTA+2*JF  )
            !if( myproc.eq.2 .and. jf.eq.1 .and. JM.eq.0 ) write(*,*) 'fou_oD1 ', IOFF, KGL, ZGTF(JF, 2*JM+IOFF), ZGTF(JF, 2*JM+1+IOFF) 
            
         end if         
      ENDDO
   ENDDO   
END DO

!$ACC END DATA

!! !$OMP END PARALLEL DO

!iimax1=0
!iimax2=0
!iimax3=0
!iunit=myproc+300
!DO KGL=IBEG,IEND,IINC
!   DO JM=0,G_NMEN_MAX      
!      DO JF=1,KFIELDS
!         IGLG = D_NPTRLS(MYSETW)+KGL-1
!         JMMAX = G_NMEN(IGLG)
!         if  (JM .le. JMMAX) then
!            IPROC = D_NPROCM(JM)
!            ISTA  = (D_NSTAGT1B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,KGL))*2*KFIELDS
!            IOFF  = 1+D_NSTAGTF(KGL)
!            iimax1 = max(iimax1,2*JF)
!            iimax2 = max(iimax2,2*JM+IOFF)
!            iimax3 = max(iimax3,ISTA+2*JF)
!            !if( jf.eq.(41+137-1) .and. JM.eq.0 ) write(iunit,*) 'fou_o ',ISTA+2*JF-1,FOUBUF_IN(ISTA+2*JF-1),FOUBUF_IN(ISTA+2*JF  )
!            if( jf.eq.1 .and. JM.eq.0 ) write(iunit,*) 'fou_o10 ', IOFF, KGL, ZGTF(2*JF-1, 2*JM+IOFF), ZGTF(2*JF, 2*JM+IOFF)
!            !if( jf.eq.1 .and. JM.eq.1 ) write(iunit,*) 'fou_o11 ', IOFF, KGL, ZGTF(JF, 2*JM+IOFF), ZGTF(JF, 2*JM+1+IOFF)
!            !if( jf.eq.1 .and. JM.eq.2 ) write(iunit,*) 'fou_o12 ', IOFF, KGL, ZGTF(JF, 2*JM+IOFF), ZGTF(JF, 2*JM+1+IOFF)
!            !if( jf.eq.1 ) write(iunit,*) 'fou_o2 ', IOFF, ZGTF(JF, 2*JM-1+IOFF),ZGTF(JF, 2*JM+IOFF),ZGTF(JF, 2*JM+1+IOFF)
 !       end if
!          
!      ENDDO
!   ENDDO
!ENDDO
!write(iunit,*), 'maxes ',iimax1,size(ZGTF,1),iimax2,size(ZGTF,2),iimax3,size(FOUBUF_IN)

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD

