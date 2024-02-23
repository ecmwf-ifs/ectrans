! (C) Copyright 1995- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
!  Copyright (c) 2020, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions
!  are met:
!
!1 . Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!2 . Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!3 . Neither the name of the copyright holder nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
!  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


!#define SEND_BY_LEVEL
!#define DEBUG
!#define DEBUG_COMM

MODULE TRGTOL_MOD

PUBLIC TRGTOL
PRIVATE TRGTOL_PROLOG, TRGTOL_COMM, TRGTOL_INIT

!============================================================================
integer, allocatable :: sendtype(:,:)   ! MPI Datatype for sending 
!integer, allocatable :: Dhor(:)  ! Number of elements in horizontal dimension for each neighbor 
integer, allocatable :: seg(:,:) ! beginning and ending element (for cases where there is padding)
integer :: max_tot  ! Maximum total number of fields received
integer :: tot_count ! Total number of received messages
integer, allocatable :: rcount(:) ! Number of arrays per neighbor
integer, allocatable :: start_blk(:) ! Beginning jblk when sending arrays
integer, allocatable :: tot_recv(:,:) ! Number of elements per array per neighbor
integer, allocatable :: count_send(:,:) ! Number of variables per array (0, 1 or 2) for sending
integer, allocatable :: count_recv(:,:) ! Number of variables per array (0, 1 or 2) for receiving
integer, allocatable :: startsend(:,:) ! Starting point for sending in field space
integer, allocatable :: start_recv(:) ! Starting point for receiving in field space
integer, allocatable :: recv_reqs(:),send_reqs(:) ! receive and send requests
integer, allocatable :: blocklength(:,:) ! Number of fields in send space
#ifndef SEND_BY_LEVEL
integer, allocatable :: num_fld(:,:) ! Number of fields per neighbor
integer, allocatable :: fstart(:,:) ! Starting index for each field
integer, allocatable :: rc2iv(:,:) ! translation table, from rcount to iv
#endif
!integer :: num_lev

!     LOCAL LOGICAL SCALARS
LOGICAL   :: LLPGPUV,LLPGP3A,LLPGP3B,LLPGP2,LLPGPONLY, LLINDER
LOGICAL, allocatable   :: LLUV(:),LLGP2(:),LLGP3A(:),LLGP3B(:)
!LOGICAL   :: LLUV(KF_GP),LLGP2(KF_GP),LLGP3A(KF_GP),LLGP3B(KF_GP)
INTEGER, allocatable :: IFLDOFF(:),IGPTROFF(:)

integer :: iuvlev,iuvpar,igp3alev,igp3apar,igp2par,igppar,igp3bpar,igp3blev
!INTEGER(KIND=JPIM) :: IUVLEVS(KF_GP),IUVPARS(KF_GP),IGP2PARS(KF_GP)
!INTEGER(KIND=JPIM) :: IGP3APARS(KF_GP),IGP3ALEVS(KF_GP),IGP3BPARS(KF_GP),IGP3BLEVS(KF_GP)
INTEGER, allocatable :: IUVLEVS(:),IUVPARS(:),IGP2PARS(:)
INTEGER, allocatable :: IGP3APARS(:),IGP3ALEVS(:),IGP3BPARS(:),IGP3BLEVS(:)

CONTAINS

SUBROUTINE TRGTOL(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)

!**** *TRGTOL * - head routine for transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform

!**   Interface.
!     ----------
!        *call* *trgtol_prolog(...)

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trgtol
!     ------------------------------------------------------------------
!        D. Pekurovsky (NVIDIA Corp.), 19 OCtober 2023: Introduced MPI Datatypes. 
!     ------------------------------------------------------------------
  
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, NPRTRNS, NPROC
USE TPM_TRANS       ,ONLY : LDIVGP, LGPNORM, LSCDERS, LUVDER, LVORGP, NGPBLKS, NPROMA

!use mpi_f08

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(OUT)   :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL     :: PGPUV(:,:,:,:) !,INTENT(IN)
REAL(KIND=JPRB),OPTIONAL     :: PGP3A(:,:,:,:) !,INTENT(IN)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP2(:,:,:)

INTEGER(KIND=JPIM) :: ISENDCOUNT
INTEGER(KIND=JPIM) :: IRECVCOUNT
INTEGER(KIND=JPIM) :: INSEND
INTEGER(KIND=JPIM) :: INRECV
INTEGER(KIND=JPIM) :: ISENDTOT (NPROC)
INTEGER(KIND=JPIM) :: IRECVTOT (NPROC)
INTEGER(KIND=JPIM) :: ISEND    (NPROC)
INTEGER(KIND=JPIM) :: IRECV    (NPROC)
INTEGER(KIND=JPIM) :: IINDEX(D%NLENGTF)
INTEGER(KIND=JPIM) :: INDOFF(NPROC)
INTEGER(KIND=JPIM) :: IGPTRSEND(2,NGPBLKS,NPRTRNS)
integer :: init = 0

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRGTOL',0,ZHOOK_HANDLE)


CALL TRGTOL_PROLOG(KF_FS,KF_GP,KVSET,&
 & ISENDCOUNT,IRECVCOUNT,INSEND,INRECV,ISENDTOT,IRECVTOT,ISEND,IRECV,IINDEX,INDOFF,IGPTRSEND)

! Transition from assumed shape arrays to explicit shape

if(init .eq. 0) then
LLINDER = .FALSE.
LLPGPUV = .FALSE.
LLPGP3A = .FALSE.
LLPGP3B = .FALSE.
LLPGP2  = .FALSE.
LLPGPONLY = .FALSE.
IF(PRESENT(KPTRGP))  LLINDER = .TRUE.
IF(PRESENT(PGP))  then
   LLPGPONLY = .TRUE.
   igppar = ubound(pgp,2)

endif
IF(PRESENT(PGPUV))  then
   LLPGPUV = .TRUE.
   iuvlev = ubound(pgpuv,2)
   iuvpar = ubound(pgpuv,3)
endif
IF(PRESENT(PGP3A))  then
   LLPGP3A = .TRUE.
   igp3alev = ubound(pgp3a,2)
   igp3apar = ubound(pgp3a,3)
endif
IF(PRESENT(PGP3B)) then
   LLPGP3B = .TRUE.
   IGP3BLEV=UBOUND(PGP3B,2)
  IGP3BPAR=UBOUND(PGP3B,3)
  IF(LSCDERS) IGP3BPAR=IGP3BPAR/3
endif
IF(PRESENT(PGP2)) then
   LLPGP2 = .TRUE.
   igp2par = ubound(pgp2,2)
   IF(LSCDERS) IGP2PAR=IGP2PAR/3
endif

! Initialization routine - call only once
   call trgtol_init(kvset,igptrsend,kf_gp,kf_fs,insend,inrecv,irecv,isend,isendtot,kptrgp,kf_scalars_g)
   init = 1
endif



CALL TRGTOL_COMM(PGLAT,KF_FS,KF_GP,KVSET, &
 & INSEND,INRECV,ISENDTOT,ISEND,IRECV,IINDEX,INDOFF,IGPTRSEND, &
 & KPTRGP,iuvlev,iuvpar,igp3alev,igp3apar,igp2par,igp3blev,igp3bpar,igppar,PGP,PGPUV,PGP3A,PGP3B,PGP2)

IF (LHOOK) CALL DR_HOOK('TRGTOL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRGTOL

SUBROUTINE TRGTOL_PROLOG(KF_FS,KF_GP,KVSET,&
 & KSENDCOUNT,KRECVCOUNT,KNSEND,KNRECV,KSENDTOT,KRECVTOT,KSEND,KRECV,KINDEX,KNDOFF,KGPTRSEND)

!**** *TRGTOL_PROLOG * - prolog for transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform
!                 the purpose is essentially 
!                 to compute the size of communication buffers in order to enable
!                 the use of automatic arrays later.


!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *call* *trgtol_prolog(...)

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trgtol
!     ------------------------------------------------------------------



USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, NPRTRNS, MYSETW, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : NGPBLKS

USE INIGPTR_MOD     ,ONLY : INIGPTR
USE PE2SET_MOD      ,ONLY : PE2SET
!
!use mpi_f08
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)

INTEGER(KIND=JPIM), INTENT(OUT) :: KSENDCOUNT
INTEGER(KIND=JPIM), INTENT(OUT) :: KRECVCOUNT
INTEGER(KIND=JPIM), INTENT(OUT) :: KNSEND
INTEGER(KIND=JPIM), INTENT(OUT) :: KNRECV
INTEGER(KIND=JPIM), INTENT(OUT) :: KSENDTOT (NPROC)
INTEGER(KIND=JPIM), INTENT(OUT) :: KRECVTOT (NPROC)
INTEGER(KIND=JPIM), INTENT(OUT) :: KSEND    (NPROC)
INTEGER(KIND=JPIM), INTENT(OUT) :: KRECV    (NPROC)
INTEGER(KIND=JPIM), INTENT(OUT) :: KINDEX(D%NLENGTF)
INTEGER(KIND=JPIM), INTENT(OUT) :: KNDOFF(NPROC)
INTEGER(KIND=JPIM), INTENT(OUT) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)

INTEGER(KIND=JPIM) :: IGPTRRECV(NPRTRNS)
INTEGER(KIND=JPIM) :: IFIRSTLAT, IGL, IGLL, ILASTLAT, JROC, IPOS, ISETB, ISETA
INTEGER(KIND=JPIM) :: ISETV, J, JFLD, JGL, JL, ISETW,  INDOFFX,IBUFLENS,IBUFLENR

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------

CALL GSTATS(1805,0)

CALL INIGPTR(KGPTRSEND,IGPTRRECV)

INDOFFX  = 0
IBUFLENS = 0
IBUFLENR = 0
KNRECV   = 0
KNSEND   = 0

DO JROC=1,NPROC

  CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)

!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1) IPOS = IPOS+1
  ENDDO
  KSENDTOT(JROC) = IGPTRRECV(ISETW)*IPOS

  IF( JROC /= MYPROC) THEN
    IBUFLENS = MAX(IBUFLENS,KSENDTOT(JROC))
    IF(KSENDTOT(JROC) > 0) THEN
      KNSEND = KNSEND+1
      KSEND(KNSEND)=JROC
    ENDIF
  ENDIF

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(ISETA))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(ISETA))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO

  KRECVTOT(JROC) = IPOS*KF_FS
  IF(KRECVTOT(JROC) > 0 .AND. MYPROC /= JROC) THEN
    KNRECV = KNRECV + 1
    KRECV(KNRECV)=JROC
  ENDIF
  IF( JROC /= MYPROC) IBUFLENR = MAX(IBUFLENR,KRECVTOT(JROC))

  IF(IPOS > 0) THEN
    KNDOFF(JROC) = INDOFFX
    INDOFFX = INDOFFX+IPOS
    IPOS = 0
    DO JGL=IFIRSTLAT,ILASTLAT
      IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
      IGLL = JGL-D%NPTRLS(MYSETW)+1
      DO JL=D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL),&
       &D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL)+D%NONL(IGL,ISETB)-1
        IPOS = IPOS+1
        KINDEX(IPOS+KNDOFF(JROC)) = JL
      ENDDO
    ENDDO
  ENDIF

ENDDO

KSENDCOUNT=0
KRECVCOUNT=0
DO J=1,NPROC
  KSENDCOUNT=MAX(KSENDCOUNT,KSENDTOT(J))
  KRECVCOUNT=MAX(KRECVCOUNT,KRECVTOT(J))
ENDDO

CALL GSTATS(1805,1)

END SUBROUTINE TRGTOL_PROLOG





!---------------------------------------------------------
!        D. Pekurovsky (NVIDIA Corp.), 19 OCtober 2023: Introduced MPI Datatypes. 
!     ------------------------------------------------------------------
subroutine trgtol_init(kvset,kgptrsend,kf_gp,kf_fs,knsend,knrecv,krecv,ksend,ksendtot,kptrgp,kf_scalars_g)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY, &
     & JP_BLOCKING_STANDARD, MPL_BARRIER, JP_BLOCKING_BUFFERED

USE TPM_GEN         ,ONLY : NOUT, NTRANS_SYNC_LEVEL
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRNS, MTAGGL,  &
     &                      MYSETV, MYSETW, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : LDIVGP, LGPNORM, LSCDERS, LUVDER, LVORGP, NGPBLKS, NPROMA

USE PE2SET_MOD      ,ONLY : PE2SET
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

use mpi

INTEGER(KIND=JPIM), INTENT(IN) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM), INTENT(IN) :: KSEND    (NPROC),knsend
INTEGER(KIND=JPIM), INTENT(OUT) :: KRECV    (NPROC)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KSENDTOT (NPROC)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCALARS_G

integer :: kf_gp
integer :: ins,isend,ifld,jfld,jblk,ifirst,ilast
integer :: type_cont,type_tmp1,type_tmp2,req_id,tot_req_id
integer(8) :: x,y
integer :: iblock(6),nblocks,jbeg,jend,val(6),val_prev(6),id,ind
integer :: dims(4,6), stride,mytype,newtype
integer ::  ierr,tot_bl
integer f,ff,ib,n,n0,iv,tot1fld,num_recv,seg_send(2),jdim_send(2),num_present,np,s
integer, allocatable :: bl(:,:),nflds(:),jdim(:,:)
integer, allocatable :: recv_reqs_prelim(:,:),typesz(:,:)
integer prec

INTEGER(KIND=JPIM) :: IOFF,IOFF1,IOFFNS,IOFFEW,J1,J2, JNR
INTEGER(KIND=JPIM) :: ISEND_FLD_START,ISEND_FLD_END

INTEGER(KIND=JPIM) :: ILEN, IPOS, ISETA, ISETB, IRECV, ISETV
INTEGER(KIND=JPIM) :: ISETW(KNSEND)
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM) :: IFLDA(KF_GP,KNSEND)
character (len=25) envalue
INTEGER(KIND = MPI_COUNT_KIND) :: lb, extent
!!!===========================================================================

!call getenv("OMPI_WANT_UMR",envalue)
!print *,myproc-1,': OMPI_WANT_UMR=',TRIM(envalue)

CALL GSTATS(1805,0)

IUVPAR=0
!IUVLEV=0
IOFF1=0
IOFFNS=KF_SCALARS_G
IOFFEW=2*KF_SCALARS_G

allocate(lluv(KF_GP))
allocate(iuvpars(KF_GP))
allocate(iuvlevs(KF_GP))
allocate(llgp2(KF_GP))
allocate(igp2pars(KF_GP))
allocate(llgp3a(KF_GP))
allocate(igp3apars(KF_GP))
allocate(igp3alevs(KF_GP))
allocate(llgp3b(KF_GP))
allocate(igp3bpars(KF_GP))
allocate(igp3blevs(KF_GP))
allocate(IFLDOFF(KF_FS))
allocate(IGPTROFF(NGPBLKS))

LLUV(:) = .FALSE.
IUVPARS(:) = -99
IUVLEVS(:) = -99
IF (LLPGPUV) THEN
  IOFF=0
  IF(LVORGP) THEN
    IUVPAR=IUVPAR+1
    DO J=1,IUVLEV
      IUVLEVS(IOFF+J)=J
      IUVPARS(IOFF+J)=IUVPAR
      LLUV(IOFF+J)=.TRUE.
    ENDDO
    IOFF=IOFF+IUVLEV
  ENDIF
  IF(LDIVGP) THEN
    IUVPAR=IUVPAR+1
    DO J=1,IUVLEV
      IUVLEVS(IOFF+J)=J
      IUVPARS(IOFF+J)=IUVPAR
      LLUV(IOFF+J)=.TRUE.
    ENDDO
    IOFF=IOFF+IUVLEV
  ENDIF
  DO J=1,IUVLEV
    IUVLEVS(IOFF+J)=J
    IUVPARS(IOFF+J)=IUVPAR+1
    IUVLEVS(IOFF+J+IUVLEV)=J
    IUVPARS(IOFF+J+IUVLEV)=IUVPAR+2
  ENDDO
  IUVPAR=IUVPAR+2
  LLUV(IOFF+1:IOFF+2*IUVLEV)=.TRUE.
  IOFF=IOFF+2*IUVLEV
  IOFF1=IOFF
  IOFFNS=IOFFNS+IOFF
  IOFFEW=IOFFEW+IOFF

  IOFF=IUVPAR*IUVLEV+KF_SCALARS_G
  IF(LUVDER) THEN
    IF(LSCDERS) IOFF=IOFF+KF_SCALARS_G
    DO J=1,IUVLEV
      IUVLEVS(IOFF+J)=J
      IUVPARS(IOFF+J)=IUVPAR+1
      LLUV(IOFF+J)=.TRUE.
      IUVLEVS(IOFF+J+IUVLEV)=J
      IUVPARS(IOFF+J+IUVLEV)=IUVPAR+2
      LLUV(IOFF+J+IUVLEV)=.TRUE.
    ENDDO
    IUVPAR=IUVPAR+2
    IOFF=IOFF+2*IUVLEV
    IOFFEW=IOFFEW+2*IUVLEV
  ENDIF
ENDIF

LLGP2(:)=.FALSE.
IF(LLPGP2) THEN
  IOFF=IOFF1
  DO J=1,IGP2PAR
    LLGP2(J+IOFF) = .TRUE.
    IGP2PARS(J+IOFF)=J
  ENDDO
  IOFF1=IOFF1+IGP2PAR
  IF(LSCDERS) THEN
    IOFF=IOFFNS
    DO J=1,IGP2PAR
      LLGP2(J+IOFF) = .TRUE.
      IGP2PARS(J+IOFF)=J+IGP2PAR
    ENDDO
    IOFFNS=IOFF+IGP2PAR
    IOFF=IOFFEW
    DO J=1,IGP2PAR
      LLGP2(J+IOFF) = .TRUE.
      IGP2PARS(J+IOFF)=J+2*IGP2PAR
    ENDDO
    IOFFEW=IOFF+IGP2PAR
  ENDIF
ENDIF

LLGP3A(:) = .FALSE.
IF(LLPGP3A) THEN
  IF(LSCDERS) IGP3APAR=IGP3APAR/3
  IOFF=IOFF1
  DO J1=1,IGP3APAR
    DO J2=1,IGP3ALEV
      LLGP3A(J2+(J1-1)*IGP3ALEV+IOFF) = .TRUE.
      IGP3APARS(J2+(J1-1)*IGP3ALEV+IOFF)=J1
      IGP3ALEVS(J2+(J1-1)*IGP3ALEV+IOFF)=J2
    ENDDO
  ENDDO
  IPAROFF=IGP3APAR
  IOFF1=IOFF1+IGP3APAR*IGP3ALEV
  IF(LSCDERS) THEN
    IOFF=IOFFNS
    DO J1=1,IGP3APAR
      DO J2=1,IGP3ALEV
        LLGP3A(J2+(J1-1)*IGP3ALEV+IOFF) = .TRUE.
        IGP3APARS(J2+(J1-1)*IGP3ALEV+IOFF)=J1+IPAROFF
        IGP3ALEVS(J2+(J1-1)*IGP3ALEV+IOFF)=J2
      ENDDO
    ENDDO
    IPAROFF=IPAROFF+IGP3APAR
    IOFFNS=IOFFNS+IGP3APAR*IGP3ALEV
    IOFF=IOFFEW
    DO J1=1,IGP3APAR
      DO J2=1,IGP3ALEV
        LLGP3A(J2+(J1-1)*IGP3ALEV+IOFF) = .TRUE.
        IGP3APARS(J2+(J1-1)*IGP3ALEV+IOFF)=J1+IPAROFF
        IGP3ALEVS(J2+(J1-1)*IGP3ALEV+IOFF)=J2
      ENDDO
    ENDDO
    IOFFEW=IOFFEW+IGP3APAR*IGP3ALEV
  ENDIF
ENDIF

LLGP3B(:) = .FALSE.
IF(LLPGP3B) THEN
  IOFF=IOFF1
  DO J1=1,IGP3BPAR
    DO J2=1,IGP3BLEV
      LLGP3B(J2+(J1-1)*IGP3BLEV+IOFF) = .TRUE.
      IGP3BPARS(J2+(J1-1)*IGP3BLEV+IOFF)=J1
      IGP3BLEVS(J2+(J1-1)*IGP3BLEV+IOFF)=J2
    ENDDO
  ENDDO
  IPAROFF=IGP3BPAR
  IOFF1=IOFF1+IGP3BPAR*IGP3BLEV
  IF(LSCDERS) THEN
    IOFF=IOFFNS
    DO J1=1,IGP3BPAR
      DO J2=1,IGP3BLEV
        LLGP3B(J2+(J1-1)*IGP3BLEV+IOFF) = .TRUE.
        IGP3BPARS(J2+(J1-1)*IGP3BLEV+IOFF)=J1+IPAROFF
        IGP3BLEVS(J2+(J1-1)*IGP3BLEV+IOFF)=J2
      ENDDO
    ENDDO
    IPAROFF=IPAROFF+IGP3BPAR
    IOFFNS=IOFFNS+IGP3BPAR*IGP3BLEV
    IOFF=IOFFEW
    DO J1=1,IGP3BPAR
      DO J2=1,IGP3BLEV
        LLGP3B(J2+(J1-1)*IGP3BLEV+IOFF) = .TRUE.
        IGP3BPARS(J2+(J1-1)*IGP3BLEV+IOFF)=J1+IPAROFF
        IGP3BLEVS(J2+(J1-1)*IGP3BLEV+IOFF)=J2
      ENDDO
    ENDDO
    IOFFEW=IOFFEW+IGP3BPAR*IGP3BLEV
  ENDIF
ENDIF

IF(KSENDTOT(MYPROC) > 0 )THEN
  IFLDS = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == MYSETV .OR. KVSET(JFLD) == -1) THEN
      IFLDS = IFLDS+1
      IF(LLINDER) THEN
        IFLDOFF(IFLDS) = KPTRGP(JFLD)
      ELSE
        IFLDOFF(IFLDS) = JFLD
      ENDIF
    ENDIF
  ENDDO
endif

CALL GSTATS(1805,1)

allocate(seg(2,knrecv))
allocate(sendtype(6,knsend))
allocate(startsend(6,knsend))
allocate(count_send(6,knsend))
allocate(count_recv(6,knrecv))
!allocate(Dhor(knrecv))
allocate(tot_recv(6,knrecv))
allocate(start_blk(knsend))
allocate(blocklength(6,knsend))

allocate(nflds(knsend))
allocate(bl(6,knrecv))
allocate(jdim(2,knrecv))
allocate(rcount(knrecv))
allocate(recv_reqs_prelim(5,knrecv))
allocate(typesz(6,knsend))
allocate(start_recv(knrecv+1))

#ifndef SEND_BY_LEVEL
allocate(num_fld(12,knrecv))
#endif

! Define array dimensions
do i=1,6
   dims(1,i) = nproma
   dims(4,i) = ngpblks
enddo

dims(2,1) = igppar
dims(2,2) = igppar
dims(2,4) = igp2par
dims(2,3) = iuvlev
dims(2,5) = igp3alev
dims(2,6) = igp3blev

dims(3,1) = 1
dims(3,2) = 1
dims(3,4) = 1
dims(3,3) = iuvpar
dims(3,5) = igp3apar
dims(3,6) = igp3bpar

! Begin bookkeeping by posting receives

DO JNR=1,KNRECV
   IRECV=KRECV(JNR)
   id = nprcids(irecv) -1
   call mpi_irecv(jdim(1,jnr),2,MPI_INTEGER,id,id+1,MPI_COMM_WORLD,recv_reqs_prelim(1,jnr),ierr) ! j beg
   call mpi_irecv(seg(1,jnr),2,MPI_INTEGER,id,id+2,MPI_COMM_WORLD,recv_reqs_prelim(2,jnr),ierr) ! segments start and end
   call mpi_irecv(tot_recv(1,jnr),6,MPI_INTEGER,id,id+3,MPI_COMM_WORLD,recv_reqs_prelim(3,jnr),ierr) ! total volume received
   call mpi_irecv(bl(1,jnr),6,MPI_INTEGER,id,id+4,MPI_COMM_WORLD,recv_reqs_prelim(4,jnr),ierr) ! Block lengths
   call mpi_irecv(count_recv(1,jnr),6,MPI_INTEGER,id,id+5,MPI_COMM_WORLD,recv_reqs_prelim(5,jnr),ierr) 

enddo

DO INS=1,KNSEND
  ISEND=KSEND(INS)
  CALL PE2SET(ISEND,ISETA,ISETB,ISETW(INS),ISETV)
  IFLD = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1 ) THEN
      IFLD = IFLD+1
      IFLDA(IFLD,INS)=JFLD
    ENDIF
  ENDDO

  nflds(ins) = ifld

ENDDO

if(jprb .eq. 4) then
prec = MPI_REAL
else
prec = MPI_DOUBLE_PRECISION
endif

stride = nproma * jprb
tot_bl = 0

DO INS=1,KNSEND
   
   ISEND_FLD_START=1
   ISEND_FLD_END = nflds(ins)

! Find continuous segments (to be used in processing the receive buffer)   
   jbeg = -1  ! Beginning jblk containing nonzero elements, for this neighbor
   jend = -2  ! Ending jblk containing nonzero elements, for this neighbor
   ji = 0
   do jblk = 1,ngpblks
      IFIRST = KGPTRSEND(1,JBLK,ISETW(INS))
      IF(IFIRST > 0) THEN    ! For those jblocks that have relevance
         if(jbeg < 0) then    ! If first such jblock
            jbeg = jblk           ! start count
            seg_send(1) = ifirst   ! start of a segment for processing receive buffer
         endif
         jend = jblk

         ilast = KGPTRSEND(2,JBLK,ISETW(INS))
         seg_send(2) = ilast + ji   ! end of a segment
         ji = ji + nproma
      endif
   enddo

   start_blk(ins) = jbeg
   
! For each array (PGPUV etc) find number of entries (vert. levels) and arrange into strided blocks for MPI datatype   
   do iv=1,6
      blocklength(iv,ins) = 0
      startsend(iv,ins) = 0
   enddo
   
   DO JJ=ISEND_FLD_START,ISEND_FLD_END
      IFLDT=IFLDA(JJ,INS)
      iv = -1
      IF(LLINDER) THEN
         iv = 1
         blocklength(iv,ins) = blocklength(iv,ins)+1
         if(startsend(iv,ins) .eq. 0) then
            startsend(iv,ins) = KPTRGP(IFLDT)
         endif
      ELSE IF(LLPGPONLY) THEN
         iv = 2
         blocklength(iv,ins) = blocklength(iv,ins)+1
         if(startsend(iv,ins) .eq. 0) then
            startsend(iv,ins) = IFLDT
         endif
      ELSEIF(LLUV(IFLDT)) THEN
         iv = 3
         blocklength(iv,ins) = blocklength(iv,ins)+1
         if(startsend(iv,ins) .eq. 0) then
            startsend(iv,ins) = iuvlevs(IFLDT)
         endif
      ELSEIF(LLGP2(IFLDT)) THEN
         iv = 4
         blocklength(iv,ins) = blocklength(iv,ins)+1
         if(startsend(iv,ins) .eq. 0) then
            startsend(iv,ins) = 1
         endif
      ELSEIF(LLGP3A(IFLDT)) THEN
         iv = 5
         blocklength(iv,ins) = blocklength(iv,ins)+1
         if(startsend(iv,ins) .eq. 0) then
            startsend(iv,ins) = igp3alevs(IFLDT)
         endif
      ELSEIF(LLGP3B(IFLDT)) THEN
         iv = 6
         blocklength(iv,ins) = blocklength(iv,ins)+1
         if(startsend(iv,ins) .eq. 0) then
            startsend(iv,ins) = igp3blevs(ifldt)
         endif
      endif
      
   enddo

   blocklength(3,ins) = blocklength(3,ins) /2

   do iv=1,6

      if(jend >= jbeg .and. blocklength(iv,ins) > 0) then

  ! Combine 1st and last dimensions of the array (PGPUV etc) into a single dimension (horizontal data), for one level (include some extra elements, that will be subtracted on the receiving side)
#ifdef SEND_BY_LEVEL
         call mpi_type_vector(jend - jbeg+1,dims(1,iv),dims(1,iv)*dims(2,iv)*dims(3,iv),prec,sendtype(iv,ins),ierr)
         typesz(iv,ins) = dims(1,iv) * (jend-jbeg+1)
#else
         call mpi_type_vector(jend - jbeg+1,dims(1,iv),dims(1,iv)*dims(2,iv)*dims(3,iv),prec,type_tmp1,ierr)
! Optionally, combine all levels into a single datatype
         call mpi_type_hvector(blocklength(iv,ins),1,stride,type_tmp1,sendtype(iv,ins),ierr)
         typesz(iv,ins) = dims(1,iv) * (jend-jbeg+1) * blocklength(iv,ins)
#endif
         
         call mpi_type_commit(sendtype(iv,ins),ierr)
         
         if(iv .ne. 3) then
            count_send(iv,ins) = 1
         else
            count_send(iv,ins) = 2
         endif
! This is the number of elements in this datatype
         tot_bl = tot_bl + blocklength(iv,ins) * count_send(iv,ins)
      else
         count_send(iv,ins) = 0
         sendtype(iv,ins) = -999
         typesz(iv,ins) = 0
      endif
      
   enddo ! iv loop

! Now send the bookkeeping information to neighbors
   jdim_send(1) = jbeg
   jdim_Send(2) = jend
   ! Beginning j of segment for this task
   ISEND=KSEND(INS)
   np = NPRCIDS(ISEND)-1
   ! myproc = id + 1
   
   call mpi_send(jdim_send,2,MPI_INTEGER,np,myproc,mpi_comm_world,ierr)
   call mpi_send(seg_send,2,MPI_INTEGER,np,myproc+1,mpi_comm_world,ierr)
   call mpi_send(typesz(1,ins),6,MPI_INTEGER,np,myproc+2,mpi_comm_world,ierr)
   call mpi_send(blocklength(1,ins),6,MPI_INTEGER,np,myproc+3,mpi_comm_world,ierr)
   call mpi_send(count_send(1,ins),6,MPI_INTEGER,np,myproc+4,mpi_comm_world,ierr)
   
enddo ! INS loop

max_tot = 0

!Complete receives for bookkeeping

do i=1,5
   call mpi_waitall(knrecv,recv_reqs_prelim(i,:),MPI_STATUSES_IGNORE,ierr)
enddo

! Finish setting up bookkeeping variables for receives
!num_lev = 0

#ifndef SEND_BY_LEVEL
allocate(fstart(12,knrecv))
allocate(rc2iv(18,knrecv))
#endif


tot_count = 0
start_recv(1) = 1
DO JNR=1,KNRECV
   IRECV=KRECV(jNR)
   id = nprcids(irecv)
!   Dhor(jnr) = (jdim(2,jnr) - jdim(1,jnr) +1) * nproma

   rcount(jnr) = 1
#ifndef SEND_BY_LEVEL
   fstart(1,jnr) = 1
#endif
   do iv=1,6 
      max_tot = max(max_tot,tot_recv(iv,jnr))
      do i=1,count_recv(iv,jnr)
!         if(num_lev .lt. num_fld(rcount,jnr)) then
!            num_lev = num_fld(rcount,jnr)
!         endif
#ifdef SEND_BY_LEVEL
         rcount(jnr) = rcount(jnr) + bl(iv,jnr)
#else
         rc2iv(rcount(jnr),jnr) = iv
         num_fld(rcount(jnr),jnr) = bl(iv,jnr)
         fstart(rcount(jnr)+1,jnr) = fstart(rcount(jnr),jnr)+num_fld(rcount(jnr),jnr)
         rcount(jnr) = rcount(jnr)+1
#endif
      enddo
      
   enddo
rcount(jnr) = rcount(jnr) -1
tot_count = tot_count + rcount(jnr)
start_recv(jnr+1) = start_recv(jnr) + rcount(jnr)
enddo


allocate(recv_reqs(tot_count))
allocate(send_reqs(tot_bl))

deallocate(bl)
deallocate(jdim)
deallocate(recv_reqs_prelim)
deallocate(nflds)
  
end subroutine trgtol_init


SUBROUTINE TRGTOL_COMM(PGLAT,KF_FS,KF_GP,KVSET,&
 & KNSEND,KNRECV,KSENDTOT,KSEND,KRECV,KINDEX,KNDOFF,KGPTRSEND,&
 & KPTRGP,iuvlev,iuvpar,igp3alev,igp3apar,igp2par,igp3blev,igp3bpar,igppar,PGP,PGPUV,PGP3A,PGP3B,PGP2)

!**** *TRGTOL_COMM * - transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *call* *trgtol(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (output)
!           PGP    -  Blocked grid point data    (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original: 95-10-01
!        D.Dent  : 97-08-04   Reorganisation to allow
!                             NPRTRV to differ from NPRGPEW
!                : 98-06-17   add mailbox control logic (from TRLTOM)
!        =99-03-29= Mats Hamrud and Deborah Salmond
!                   JUMP in FFT's changed to 1
!                   KINDEX introduced and ZCOMBUF not used for same PE
!         01-11-23  Deborah Salmond and John Hague
!                    LIMP_NOOLAP Option for non-overlapping message passing
!                    and buffer packing
!         01-12-18  Peter Towers
!                   Improved vector performance of GTOL_PACK,GTOL_UNPACK
!         03-04-02  G. Radnoti: call barrier always when nproc>1
!         08-01-01  G.Mozdzynski: cleanup
!         09-01-02  G.Mozdzynski: use non-blocking recv and send
!         23-10-19  D. Pekurovsky (NVIDIA Corp.): Introduced MPI Datatypes. 
!     ------------------------------------------------------------------



USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY, &
     & JP_BLOCKING_STANDARD, MPL_BARRIER, JP_BLOCKING_BUFFERED

USE TPM_GEN         ,ONLY : NOUT, NTRANS_SYNC_LEVEL
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRNS, MTAGGL,  &
     &                      MYSETV, MYSETW, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : LDIVGP, LGPNORM, LSCDERS, LUVDER, LVORGP, NGPBLKS, NPROMA

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS


!
use mpi

IMPLICIT NONE

!include 'mpif.h'


INTEGER(KIND=JPIM) :: IUVPAR,IUVLEV,IGP2PAR,IGP3ALEV,IGP3APAR,IGP3BLEV,IGP3BPAR,IPAROFF,igppar
INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(OUT)   :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM), INTENT(IN) :: KNSEND
INTEGER(KIND=JPIM), INTENT(IN) :: KNRECV
INTEGER(KIND=JPIM), INTENT(IN) :: KSENDTOT (NPROC)
INTEGER(KIND=JPIM), INTENT(IN) :: KSEND    (NPROC)
INTEGER(KIND=JPIM), INTENT(IN) :: KRECV    (NPROC)
INTEGER(KIND=JPIM), INTENT(IN) :: KINDEX(D%NLENGTF)
INTEGER(KIND=JPIM), INTENT(IN) :: KNDOFF(NPROC)
INTEGER(KIND=JPIM), INTENT(IN) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP(nproma,igppar,ngpblks)
#ifdef DEBUG_COMM
REAL(KIND=JPRB),OPTIONAL     :: PGPUV(nproma,iuvlev,iuvpar,ngpblks) 
REAL(KIND=JPRB),OPTIONAL     :: PGP3A(nproma,igp3alev,igp3apar,ngpblks) 
REAL(KIND=JPRB),OPTIONAL     :: PGP2(nproma,igp2par,ngpblks) 
#else
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGPUV(nproma,iuvlev,iuvpar,ngpblks) 
!REAL(KIND=JPRB),OPTIONAL     :: PGPUV(:,:,:,:) !,INTENT(IN)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3A(nproma,igp3alev,igp3apar,ngpblks) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP2(nproma,igp2par,ngpblks) 
#endif
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3B(nproma,igp3blev,igp3bpar,ngpblks)

! Send buffer is no longer needed: in this version we use MPI Datatypes to copy data directly from the origina arrays
!REAL(KIND=JPRB), ALLOCATABLE :: ZCOMBUFS(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZCOMBUFR(:,:)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFIRST, ILAST
INTEGER(KIND=JPIM) :: ISEND, JBLK, IPOS,JFLD, JK, JL, IFLD, II, IFLDS, INS, INR
INTEGER(KIND=JPIM) :: JJ,JI,IFLDT, J, i, ierr, jnr
integer tmp,r,i2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR
real(kind=jprb) :: num,num1,num2


!============================================================================
integer :: req_id,tot_req_id
integer :: id,ind,fcount,req_count,prec,irecv
integer f,ff,n,n0,iv,np,s
character(len=25) :: s1,s2,s3,str
integer k,m,lev,l
!============================================================================

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------


#ifdef DEBUG_COMM
do m=1,ngpblks
do k=1,iuvpar
do j=1,iuvlev
   do i=1,nproma
      num1 = real(j-1,jprb)*100.0d0 + real(k-1,jprb)*100000.0d0 + real(m-1,jprb)*1000000.0d0
      num2 = real(i,jprb) + real(myproc-1,jprb)*0.001d0
      num = num1 + num2
      pgpuv(i,j,k,m) = num
enddo
enddo
enddo
do k=1,igp3apar
do j=1,igp3alev
   do i=1,nproma
      num1 = real(j-1,jprb)*100.0d0 + real(k-1,jprb)*100000.0d0 + real(m-1,jprb)*1000000.0d0
      num2 = real(i,jprb) + real(myproc-1,jprb)*0.001d0
      num = num1 + num2
      pgp3a(i,j,k,m) = -num
      pgp2(i,1,m) = real(m-1,jprb)*100.0d0 + num2 +0.5d0
   enddo
enddo
enddo
enddo

#endif

!ITAG = MTAGGL

IF (LHOOK) CALL DR_HOOK('TRGTOL_BAR',0,ZHOOK_HANDLE_BAR)
CALL GSTATS_BARRIER(761)
IF (LHOOK) CALL DR_HOOK('TRGTOL_BAR',1,ZHOOK_HANDLE_BAR)

IF(.NOT.LGPNORM)THEN
   CALL GSTATS(803,0)
ELSE
   CALL GSTATS(804,0)
ENDIF

ALLOCATE(ZCOMBUFR(max_tot,tot_count))

!===============================================================================

! Post receives

recv_reqs = MPI_REQUEST_NULL
req_count = 0

if(jprb .eq. 4) then
prec = MPI_REAL
else
prec = MPI_DOUBLE_PRECISION
endif

DO JNR=1,KNRECV  ! Loop over neighbors
   IRECV=KRECV(JNR)
   id = nprcids(irecv) -1
!   fcount = 0
!   do iv=1,6  ! Loop over arrays

   do i=1,rcount(jnr)
!      do i = 1,count_recv(iv,jnr)  ! Loop over how many variables (U/V) inside one array (0 to 2)
!#ifdef SEND_BY_LEVEL
!         do j=0,num_fld(iv,jnr)-1
!            fcount = fcount+1 ! Number of variables for one neighbor
!            req_count = req_count +1  ! Total number of requests
!            call mpi_irecv(zcombufr(1,fcount,jnr),tot_recv(iv,jnr),prec,&
!                 id,100000000*iv+1000000*j + (id+1)*10+i,mpi_comm_world,recv_reqs(req_count),ierr)
!         enddo
!#else         
!         fcount = fcount+1 ! Number of variables for one neighbor
         req_count = req_count +1  ! Total number of requests
         call mpi_irecv(zcombufr(1,req_count),max_tot,prec,id, &
              (id+1)*10000+i,mpi_comm_world,recv_reqs(req_count),ierr)
!#endif
!      enddo
    enddo
enddo

! Post sends, for each active array, using MPI Datatypes 
!tot_req_id = 0
req_id = 1
DO INS=1,KNSEND
   ISEND=KSEND(INS)
   np = nprcids(isend) -1
   r = 1
   if(LLINDER) then
#ifdef SEND_BY_LEVEL
      do j=0,blocklength(1,ins)-1
#else
         j = 0
#endif
         i2 = startsend(1,ins)+j
         call mpi_isend(pgp(1,i2,start_blk(ins)),1,sendtype(1,ins),np, &
           myproc*10000+r,mpi_comm_world,send_reqs(req_id),ierr)
         req_id = req_id +1
         r = r+1
#ifdef SEND_BY_LEVEL
      enddo
#endif

   elseif(LLPGPONLY) then
#ifdef SEND_BY_LEVEL
      do j=0,blocklength(2,ins)-1
#else
         j=0
#endif
         i2 = startsend(2,ins)+j
         call mpi_isend(pgp(1,i2,start_blk(ins)),1,sendtype(2,ins),np, &
              myproc*10000+r,mpi_comm_world,send_reqs(req_id),ierr)
         req_id = req_id +1
         r = r+1
#ifdef SEND_BY_LEVEL
      enddo
#endif
   else
      do i=1,count_send(3,ins)
#ifdef SEND_BY_LEVEL
         do j=0,blocklength(3,ins)-1
#else
         j=0
#endif
         i2 = startsend(3,ins)+j
         call mpi_isend(pgpuv(1,i2,i,start_blk(ins)),1,sendtype(3,ins),np, &
              myproc*10000+r,mpi_comm_world,send_reqs(req_id),ierr)
         req_id = req_id +1
         r = r+1
#ifdef SEND_BY_LEVEL
      enddo
#endif
      enddo

      do i=1,count_send(4,ins)
#ifdef SEND_BY_LEVEL
         do j=0,blocklength(4,ins)-1
#else
         j=0
#endif
         i2 = startsend(4,ins)+j
         call mpi_isend(pgp2(1,i2,start_blk(ins)),1,sendtype(4,ins),np, &
              myproc*10000+r,mpi_comm_world,send_reqs(req_id),ierr)
         req_id = req_id +1
         r = r+1
#ifdef SEND_BY_LEVEL
         enddo
#endif
      enddo

      do i=1,count_send(5,ins)
#ifdef SEND_BY_LEVEL
         do j=0,blocklength(5,ins)-1
#else
         j=0
#endif
         i2 = startsend(5,ins)+j
         call mpi_isend(pgp3a(1,i2,i,start_blk(ins)),1,sendtype(5,ins),np, &
             myproc*10000+r,mpi_comm_world,send_reqs(req_id),ierr)
         req_id = req_id +1
         r = r+1
#ifdef SEND_BY_LEVEL
         enddo
#endif
      enddo

      do i=1,count_send(6,ins)
#ifdef SEND_BY_LEVEL
         do j=0,blocklength(6,ins)-1
#else
         j=0
#endif
         i2 = startsend(6,ins)+j
         call mpi_isend(pgp3b(1,i2,i,start_blk(ins)),1,sendtype(6,ins),np, &
                myproc*10000+r,mpi_comm_world,send_reqs(req_id),ierr)
         req_id = req_id +1
         r = r+1
#ifdef SEND_BY_LEVEL
         enddo
#endif
      enddo
   endif
  
ENDDO

IF(.NOT.LGPNORM)THEN
   CALL GSTATS(803,1)
ELSE
   CALL GSTATS(804,1)
ENDIF

tot_req_id = req_id-1

! Copy local contribution

IF(KSENDTOT(MYPROC) > 0 )THEN
  IFLDS = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == MYSETV .OR. KVSET(JFLD) == -1) THEN
      IFLDS = IFLDS+1
      IF(LLINDER) THEN
        IFLDOFF(IFLDS) = KPTRGP(JFLD)
      ELSE
        IFLDOFF(IFLDS) = JFLD
      ENDIF
    ENDIF
  ENDDO

  IPOS=0
  DO JBLK=1,NGPBLKS
    IGPTROFF(JBLK)=IPOS
    IFIRST = KGPTRSEND(1,JBLK,MYSETW)
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,MYSETW)
      IPOS=IPOS+ILAST-IFIRST+1
    ENDIF
  ENDDO
  CALL GSTATS(1601,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
  DO JBLK=1,NGPBLKS
    IFIRST = KGPTRSEND(1,JBLK,MYSETW)
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,MYSETW)
      IF(LLPGPONLY) THEN
        DO JFLD=1,IFLDS
          IFLD = IFLDOFF(JFLD)
          DO JK=IFIRST,ILAST
            IPOS = KNDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
            PGLAT(JFLD,KINDEX(IPOS)) = PGP(JK,IFLD,JBLK)
          ENDDO
        ENDDO
      ELSE
        DO JFLD=1,IFLDS
          IFLD = IFLDOFF(JFLD)
          IF(LLUV(IFLD)) THEN
            DO JK=IFIRST,ILAST
              IPOS = KNDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
              PGLAT(JFLD,KINDEX(IPOS)) = PGPUV(JK,IUVLEVS(IFLD),IUVPARS(IFLD),JBLK)
            ENDDO
          ELSEIF(LLGP2(IFLD)) THEN
            DO JK=IFIRST,ILAST
              IPOS = KNDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
              PGLAT(JFLD,KINDEX(IPOS)) = PGP2(JK,IGP2PARS(IFLD),JBLK)
            ENDDO
          ELSEIF(LLGP3A(IFLD)) THEN
            DO JK=IFIRST,ILAST
              IPOS = KNDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
              PGLAT(JFLD,KINDEX(IPOS)) = PGP3A(JK,IGP3ALEVS(IFLD),IGP3APARS(IFLD),JBLK)
            ENDDO
          ELSEIF(LLGP3B(IFLD)) THEN
            DO JK=IFIRST,ILAST
              IPOS = KNDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
              PGLAT(JFLD,KINDEX(IPOS)) = PGP3B(JK,IGP3BLEVS(IFLD),IGP3BPARS(IFLD),JBLK)
            ENDDO
          ELSE
            WRITE(NOUT,*)'TRGTOL_MOD: ERROR',JFLD,IFLD
            CALL ABORT_TRANS('TRGTOL_MOD: ERROR')
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1601,1)

ENDIF


!  Unpack loop.........................................................

! Wait for the data and process the receive buffer (remove redundant elements between segments)

IF(.NOT.LGPNORM)THEN
   CALL GSTATS(803,0)
ELSE
   CALL GSTATS(804,0)
ENDIF

#ifdef DEBUG
s1 = 'waitany.'
write(s2,9) myproc-1
9 format(i0)
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')
#endif

DO JNR=1,tot_count
#ifdef DEBUG
      write(11,12) myproc-1,jnr
      call flush
#endif
12    format(i3,': Waitany loop, count ',i3)
      call mpi_waitany(tot_count,recv_reqs,ind,MPI_STATUS_IGNORE,ierr)
      call find_source(ind,inr,f,start_recv,knrecv)
!      if(f .gt. kf_fs) then
!         print *,myproc-1,': Beyong bounds, f=',f, 'ind=',ind,'inr=',inr
!      endif
      
!      call mpi_wait(recv_reqs(jnr),MPI_STATUS_IGNORE,ierr)

!      call find_source(jnr,inr,f,start_recv,knrecv)
!
!      ind = jnr
!      f = mod(ind-1,rcount) +1   ! Group of fields
!      inr = (ind-1)/rcount +1 

      
      IRECV=KRECV(INR)    ! Source ID

#ifdef SEND_BY_LEVEL
      n0 = 0
      s = f
#else      
      iv = rc2iv(f,inr)
      l = tot_recv(iv,inr)/num_fld(f,inr)
      do ff=0,num_fld(f,inr)-1   ! Fields within the block
!         n0 = ff*Dhor(inr) ! Offset 
         n0 = ff*l  ! Offset 
         s = fstart(f,inr)+ff
#endif
         do i=seg(1,inr),seg(2,inr)  ! Copy one segment (omitting cut elements)
            ii = kindex(kndoff(irecv)+i-seg(1,inr)+1)
            PGLAT(s,ii) = zcombufr(i+n0,ind)
         enddo
#ifndef SEND_BY_LEVEL
      enddo
#endif   
   enddo

#ifdef DEBUG
   write (11,13) myproc-1
13 format(i3,': Loop completed')
   call flush
   close(11)
#endif

   call mpi_waitall(tot_req_id,send_reqs,MPI_STATUSES_IGNORE,ierr)
   

#ifdef DEBUG_COMM
s1 = 'pglat.'
write(s2,9) myproc-1
9 format(i0)
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')

do j=1,D%NLENGTF
   do i=1,KF_FS
      write(11,10) i,j,pglat(i,j)
   enddo
enddo

10 format(i2,i8,F16.3)

close(11)
#endif
      
!IF (NTRANS_SYNC_LEVEL >= 1) THEN
   CALL MPL_BARRIER(CDSTRING='TRGTOL_COMM: BARRIER AT END')
!ENDIF

IF(.NOT.LGPNORM)THEN
  CALL GSTATS(803,1)
ELSE
  CALL GSTATS(804,1)
ENDIF

DEALLOCATE(ZCOMBUFR)

!==============================================

!==============================================

CALL GSTATS_BARRIER2(761)

END SUBROUTINE TRGTOL_COMM

subroutine find_source(ind,src,f,ref,n)
  implicit none
  integer, intent(in) :: ind,n,ref(n)
  integer, intent(out) :: src,f
  integer low,high,mid
  logical conv
  
  ! Do binary search
  low = 1
  high = n
  conv = .false.
  do while(low < high .and. .not. conv)
     mid = (low + high + 1)/2
     
     if(ref(mid) <= ind) then
        low = mid
        if(ref(mid) .eq. ind) then
           conv = .true.
        endif
     else if(ref(mid) > ind) then
        high = mid - 1
     else
        conv = .true.
     endif
     
  enddo
  src = low
  f = ind - ref(src) +1
  
end subroutine find_source

END MODULE TRGTOL_MOD

