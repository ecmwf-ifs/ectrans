! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INIGPTR_MOD
CONTAINS
SUBROUTINE INIGPTR(KGPTRSEND,KGPTRRECV)

!     Compute tables to assist GP to/from Fourier space transpositions

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NOUT
USE TPM_DISTR       ,ONLY : D, NPRTRNS
USE TPM_TRANS       ,ONLY : NGPBLKS, NPROMA
USE EQ_REGIONS_MOD  ,ONLY : MY_REGION_EW
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(OUT) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM),INTENT(OUT) :: KGPTRRECV(NPRTRNS)

INTEGER(KIND=JPIM) :: IBLOCK,IROF,IBFIRST,IPROCLAST,IPROC,IFIRST,ILAST,IBLAST
INTEGER(KIND=JPIM) :: JGL,JBL,JPRTRNS,JBLKS
!     Compute tables to assist GP to/from Fourier space transpositions


KGPTRSEND(:,:,:)=0
IBLOCK=1
IROF=1
IBFIRST=1
IPROCLAST=D%NPROCL(D%NFRSTLOFF+1)
! for each latitude on this processor
DO JGL=1,D%NDGL_GP
  ! find the processor where this row should be saved in the fourier distribution
  ! this is called the "w-set"
  IPROC=D%NPROCL(D%NFRSTLOFF+JGL)

  ! for each latitude on this processor, find first and last points
  ! for each NPROMA chunk, for each destination processor
  IF(IPROC /= IPROCLAST) THEN
    ! we got onto a new process, we still need to finish the last block of the previous
    ! process
    IF(IROF > 1) THEN
      KGPTRSEND(1,IBLOCK,IPROCLAST)=IBFIRST
      KGPTRSEND(2,IBLOCK,IPROCLAST)=IROF-1
    ENDIF
    IF(IROF <= NPROMA) IBFIRST=IROF
    IPROCLAST=IPROC
  ENDIF
  ! my offset of the first gridpoint in this row (globally, in EW-direction)
  IFIRST=D%NSTA(D%NPTRFLOFF+JGL,MY_REGION_EW)
  ! my offset of the last gridpoint in this row (globally, in EW-direction)
  ILAST =IFIRST + D%NONL(D%NPTRFLOFF+JGL,MY_REGION_EW) -1
  ! now go through all gridpoints on this latitude
  DO JBL=IFIRST,ILAST
    IF(IROF == NPROMA) THEN
      ! this block is full!
      IBLAST=IROF
      KGPTRSEND(1,IBLOCK,IPROC)=IBFIRST
      KGPTRSEND(2,IBLOCK,IPROC)=IBLAST
      IF(IBLOCK < NGPBLKS) IBLOCK=IBLOCK+1
      IROF=0
      IBFIRST=1
    ENDIF
    IROF=IROF+1
  ENDDO
ENDDO
IF(IROF /= 1.AND.IROF /= IBFIRST) THEN
  ! non-empty residual block after last latitude line
  IBLAST=IROF-1
  KGPTRSEND(1,IBLOCK,IPROC)=IBFIRST
  KGPTRSEND(2,IBLOCK,IPROC)=IBLAST
ENDIF
!         sum up over blocks
KGPTRRECV(:)=0
DO JPRTRNS=1,NPRTRNS
  DO JBLKS=1,NGPBLKS
    IF(KGPTRSEND(1,JBLKS,JPRTRNS) > 0) THEN
      KGPTRRECV(JPRTRNS)=KGPTRRECV(JPRTRNS)+&
       &KGPTRSEND(2,JBLKS,JPRTRNS)-KGPTRSEND(1,JBLKS,JPRTRNS)+1
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE INIGPTR
END MODULE INIGPTR_MOD
