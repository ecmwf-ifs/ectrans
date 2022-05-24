! (C) Copyright 2008- ECMWF.
! (C) Copyright 2008- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GPNORM_TRANS(PGP,KFIELDS,KPROMA,PAVE,PMIN,PMAX,LDAVE_ONLY,KRESOL)


!**** *GPNORM_TRANS* - calculate grid-point norms

!     Purpose.
!     --------
!        calculate grid-point norms using a 2 stage (NPRTRV,NPRTRW) communication rather
!        than an approach using a more expensive global gather collective communication

!**   Interface.
!     ----------
!     CALL GPNORM_TRANS(...)

!     Explicit arguments :
!     --------------------
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP is  dimensioned (NPROMA,KFIELDS,NGPBLKS) where
!                  NPROMA is the blocking factor, KFIELDS the total number
!                  of fields and NGPBLKS the number of NPROMA blocks.
!     KFIELDS     - number of fields (input)
!                   (these do not have to be just levels)
!     KPROMA      - required blocking factor (input)
!     PAVE        - average (output)
!     PMIN        - minimum (input/output)
!     PMAX        - maximum (input/output)
!     LDAVE_ONLY  - T : PMIN and PMAX already contain local MIN and MAX
!     KRESOL      -  resolution tag (optional)
!                    default assumes first defined resolution
!

!     Author.
!     -------
!        George Mozdzynski *ECMWF*

!     Modifications.
!     --------------
!        Original : 19th Sept 2008
!        R. El Khatib 07-08-2009 Optimisation directive for NEC

!     ------------------------------------------------------------------

USE PARKIND1        ,ONLY : JPIM     ,JPRB , JPRD
USE PARKIND_ECTRANS ,ONLY : JPRBT

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NOUT
USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : LGPNORM, NGPBLKS, NPROMA
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRV, NPRTRW, MYSETV, MYSETW, NPROC, D_NSTAGTF,D_NPTRLS, MYPROC
USE TPM_GEOMETRY    ,ONLY : G,G_NLOEN,G_NLOEN_MAX
USE TPM_FIELDS      ,ONLY : F_RW
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
!USE TRGTOL_MOD      ,ONLY : TRGTOL_CUDAAWARE
USE SET2PE_MOD      ,ONLY : SET2PE
USE MPL_MODULE      ,ONLY : MPL_RECV, MPL_SEND, JP_BLOCKING_STANDARD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK,  JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVE(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMIN(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMAX(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
LOGICAL           ,INTENT(IN)    :: LDAVE_ONLY
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN)  :: KRESOL

!ifndef INTERFACE

! Local variables
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: IUBOUND(4)
INTEGER(KIND=JPIM) :: IVSET(KFIELDS)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETS(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETG(:,:)
!ON CPU
REAL(KIND=JPRBT),ALLOCATABLE :: ZGTFL(:,:), ZGTF(:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAVE(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZMINGL(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZMAXGL(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZMINGPN(:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZMAXGPN(:)

REAL(KIND=JPRD),ALLOCATABLE :: ZAVEG(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZMING(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZMAXG(:)
REAL(KIND=JPRD),ALLOCATABLE :: ZSND(:)
REAL(KIND=JPRD),ALLOCATABLE :: ZRCV(:)
INTEGER(KIND=JPIM) :: J,JGL,IGL,JL,JF,IF_GP,IF_SCALARS_G,IF_FS,JSETV,JSETW,IWLATS,JMAX
INTEGER(KIND=JPIM) :: IPROC,ITAG,ILEN,ILENR,IBEG,IEND,IND
!INTEGER(KIND=JPIM) :: iunit

!     ------------------------------------------------------------------
CALL ABORT_TRANS("GPNORM NOT IMPLEMENTED")

END SUBROUTINE GPNORM_TRANS
