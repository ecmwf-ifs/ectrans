! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE transform_driver
use parkind1, only: jpim, jprb, jprd
use oml_mod ,only : oml_max_threads
use yomgstats, only: jpmaxstat
use yomhook, only : dr_hook_init

implicit none
external timef
real(kind=jprd) :: timef
!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "specnorm.h"

!!===================================================================================================
CONTAINS
!===================================================================================================
! ecTrans setup routines
!===================================================================================================
SUBROUTINE ectrans_setup0(kout,kerr,kverbosity,kprgpns,kprgpew,kprtrw,lduse_mpi)
integer(kind=jpim), intent(in) :: kout, kerr
integer, intent(in) :: kverbosity
integer(kind=jpim), intent(in) :: kprgpns, kprgpew, kprtrw
logical, intent(in) :: lduse_mpi 
integer(kind=jpim) :: nmax_resol = 37 ! Max number of resolutions
integer(kind=jpim) :: npromatr = 0 ! nproma for trans lib
integer(kind=jpim) :: kcombflen = 1800000 ! Size of comm buffer
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag
real(kind=jprb) :: zra = 6371229._jprd

call setup_trans0(kout=kout, kerr=kerr, kprintlev=merge(2, 0, kverbosity == 1),                &
  &               kmax_resol=nmax_resol, kpromatr=npromatr, kprgpns=kprgpns, kprgpew=kprgpew, &
  &               kprtrw=kprtrw, kcombflen=kcombflen, ldsync_trans=lsync_trans,               &
  &               ldeq_regions=leq_regions, prad=zra, ldalloperm=.true., ldmpoff=lduse_mpi)

END SUBROUTINE ectrans_setup0

SUBROUTINE ectrans_setup(ksmax,kdgl,kloen, lduseflt)
integer(kind=jpim), intent(in) :: ksmax,kdgl 
integer(kind=jpim), intent(in) :: kloen(:)
logical, intent(in) :: lduseflt
logical :: ldsplit = .true.
logical :: lfftw = .true. ! Use FFTW for Fourier transforms
logical :: luserpnm = .false.
logical :: lkeeprpnm = .false.
call setup_trans(ksmax=ksmax, kdgl=kdgl, kloen=kloen, ldsplit=ldsplit,          &
  &                 ldusefftw=lfftw, lduserpnm=luserpnm, ldkeeprpnm=lkeeprpnm, &
  &                 lduseflt=lduseflt)
END SUBROUTINE ectrans_setup

SUBROUTINE ectrans_trans_inq(kspec2, kspec2g, kgptot, kgptotg)
integer(kind=jpim), intent(out) :: kspec2
integer(kind=jpim), intent(out) :: kspec2g
integer(kind=jpim), intent(out) :: kgptot
integer(kind=jpim), intent(out) :: kgptotg
call trans_inq(kspec2=kspec2, kspec2g=kspec2g, kgptot=kgptot, kgptotg=kgptotg)
END SUBROUTINE ectrans_trans_inq

SUBROUTINE ectrans_allocate_spectral(kflevl,kspec2,kfld,ksmax)
USE transform_driver_data_mod, ONLY : zspvor,zspdiv,zspsc3a,zspsc2
USE transform_driver_data_mod, ONLY : sp3d
integer(kind=jpim), intent(in) :: kflevl,kspec2,kfld,ksmax
! Try to mimick IFS layout as much as possible
nullify(zspvor)
nullify(zspdiv)
nullify(zspsc3a)
allocate(sp3d(kflevl,kspec2,2+kfld))
allocate(zspsc2(1,kspec2))
call initialize_spectral_arrays(ksmax, zspsc2, sp3d)
! Point convenience variables to storage variable sp3d
zspvor  => sp3d(:,:,1)
zspdiv  => sp3d(:,:,2)
zspsc3a => sp3d(:,:,3:3+(kfld-1))
END SUBROUTINE ectrans_allocate_spectral


SUBROUTINE ectrans_allocate_grid(nproma, ngpblks, nfld,nflevg, &
        & ldvordiv, lduvders, ldscders) 
  USE transform_driver_data_mod, ONLY:  zgmv,zgmvs
  USE transform_driver_data_mod, ONLY:  zgpuv, zgp3a, zgp2
  integer(kind=jpim), intent(in) :: nproma
  integer(kind=jpim), intent(in) :: ngpblks
  integer(kind=jpim), intent(in) :: nfld
  integer(kind=jpim), intent(in) :: nflevg
  logical, intent(in) :: ldvordiv, lduvders, ldscders
  integer(kind=jpim) :: ndimgmv  = 0 ! Third dim. of gmv "(nproma,nflevg,ndimgmv,ngpblks)"
  integer(kind=jpim) :: ndimgmvs = 0 ! Second dim. gmvs "(nproma,ndimgmvs,ngpblks)"
  integer(kind=jpim) :: jbegin_uv = 0
  integer(kind=jpim) :: jend_uv   = 0
  integer(kind=jpim) :: jbegin_sc = 0
  integer(kind=jpim) :: jend_sc   = 0
  integer(kind=jpim) :: jbegin_scder_NS = 0
  integer(kind=jpim) :: jend_scder_NS = 0
  integer(kind=jpim) :: jbegin_scder_EW = 0
  integer(kind=jpim) :: jend_scder_EW = 0
  integer(kind=jpim) :: jbegin_uder_EW = 0
  integer(kind=jpim) :: jend_uder_EW = 0
  integer(kind=jpim) :: jbegin_vder_EW = 0
  integer(kind=jpim) :: jend_vder_EW = 0


 if (ldvordiv) then
    jbegin_uv = 1
    jend_uv = 2
  endif
  if (lduvders) then
    jbegin_uder_EW  = jend_uv + 1
    jend_uder_EW    = jbegin_uder_EW + 1
    jbegin_vder_EW  = jend_uder_EW + 1
    jend_vder_EW    = jbegin_vder_EW + 1
  else
    jbegin_uder_EW = jend_uv
    jend_uder_EW   = jend_uv
    jbegin_vder_EW = jend_uv
    jend_vder_EW   = jend_uv
  endif

  jbegin_sc = jbegin_vder_EW + 1
  jend_sc   = jbegin_vder_EW + nfld

  if (ldscders) then
    ndimgmvs = 3
    jbegin_scder_NS = jend_sc + 1
    jend_scder_NS   = jend_sc + nfld
    jbegin_scder_EW = jend_scder_NS + 1
    jend_scder_EW   = jend_scder_NS + nfld
  else
    ndimgmvs = 1
    jbegin_scder_NS = jend_sc
    jend_scder_NS   = jend_sc
    jbegin_scder_EW = jend_sc
    jend_scder_EW   = jend_sc
  endif

  ndimgmv = jend_scder_EW



!===================================================================================================
! Allocate gridpoint arrays
!===================================================================================================
allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
allocate(zgmvs(nproma,ndimgmvs,ngpblks))

zgpuv => zgmv(:,:,1:jend_vder_EW,:)
zgp3a => zgmv(:,:,jbegin_sc:jend_scder_EW,:)
zgp2  => zgmvs(:,:,:)

END SUBROUTINE ectrans_allocate_grid
SUBROUTINE ectrans_deallocate_grid
USE transform_driver_data_mod, ONLY : zgmv, zgmvs 
deallocate(zgmv)
deallocate(zgmvs)
END SUBROUTINE ectrans_deallocate_grid

!===================================================================================================
! Allocate norm arrays
!===================================================================================================

SUBROUTINE ectrans_allocate_normdata(nflevl, nflevg)
USE transform_driver_data_mod, ONLY : zspvor,zspdiv,zspsc3a,zspsc2
USE transform_driver_data_mod, ONLY : znormsp, znormsp1 
USE transform_driver_data_mod, ONLY : znormvor, znormvor1 
USE transform_driver_data_mod, ONLY : znormdiv, znormdiv1 
USE transform_driver_data_mod, ONLY : znormt, znormt1 
  integer(kind=jpim), intent(in) :: nflevl,nflevg
  allocate(znormsp(1))
  allocate(znormsp1(1))
  allocate(znormvor(nflevg))
  allocate(znormvor1(nflevg))
  allocate(znormdiv(nflevg))
  allocate(znormdiv1(nflevg))
  allocate(znormt(nflevg))
  allocate(znormt1(nflevg))

END SUBROUTINE ectrans_allocate_normdata

SUBROUTINE ectrans_calculate_norms(indx, nflevl, nflevg, ivset,ivsetsc)
USE transform_driver_data_mod, ONLY : zspvor,zspdiv,zspsc3a,zspsc2
USE transform_driver_data_mod, ONLY : znormsp, znormsp1 
USE transform_driver_data_mod, ONLY : znormvor, znormvor1 
USE transform_driver_data_mod, ONLY : znormdiv, znormdiv1 
USE transform_driver_data_mod, ONLY : znormt, znormt1 
  integer(kind=jpim), intent(in) :: indx,nflevl,nflevg
  integer(kind=jpim), intent(in) :: ivset(:),ivsetsc(1)
  select case(indx)
    case(0) 
  call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt,   kvset=ivset(1:nflevg))
  call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp,  kvset=ivsetsc)
    case(1) 
  call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor1, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv1, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt1,   kvset=ivset(1:nflevg))
  call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp1,  kvset=ivsetsc)
    case(2) 
  call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor, kvset=ivset)
  call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv, kvset=ivset)
  call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt,   kvset=ivset)
  call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp,  kvset=ivsetsc)
  end select
END SUBROUTINE ectrans_calculate_norms

SUBROUTINE ectrans_print_norms_init(nout,nflevg)
USE transform_driver_data_mod, ONLY : znormsp1 
USE transform_driver_data_mod, ONLY : znormvor1 
USE transform_driver_data_mod, ONLY : znormdiv1 
USE transform_driver_data_mod, ONLY : znormt1 
  integer(kind=jpim), intent(in) :: nout,nflevg
  integer(kind=jpim) :: ifld
    do ifld = 1, nflevg
      write(nout,'("norm zspvor( ",i4,",:)   = ",f20.15)') ifld, znormvor1(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv1(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormt1(ifld)
    enddo
    do ifld = 1, 1
      write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15)') ifld, znormsp1(ifld)
    enddo
END SUBROUTINE ectrans_print_norms_init
!===================================================================================================
! Setup timers
!===================================================================================================


SUBROUTINE ectrans_allocate_timers(iters)
USE transform_driver_data_mod, ONLY : ztstep, ztstep1, ztstep2 
USE transform_driver_data_mod, ONLY : ztstepmax , ztstepmin , ztstepavg  
USE transform_driver_data_mod, ONLY : ztstepmax1, ztstepmin1, ztstepavg1 
USE transform_driver_data_mod, ONLY : ztstepmax2, ztstepmin2, ztstepavg2 
  integer(kind=jpim), intent(in) :: iters
allocate(ztstep(iters))
allocate(ztstep1(iters))
allocate(ztstep2(iters))

ztstepavg  = 0._jprd
ztstepmax  = 0._jprd
ztstepmin  = 9999999999999999._jprd
ztstepavg1 = 0._jprd
ztstepmax1 = 0._jprd
ztstepmin1 = 9999999999999999._jprd
ztstepavg2 = 0._jprd
ztstepmax2 = 0._jprd
ztstepmin2 = 9999999999999999._jprd
END SUBROUTINE ectrans_allocate_timers

SUBROUTINE ectrans_set_ztstep_start(indx,jstep)
USE transform_driver_data_mod, ONLY : ztstep,ztstep1,ztstep2
  integer(kind=jpim), intent(in) :: indx,jstep 
  select case(indx)
     case (0)
       ztstep(jstep) = timef()
     case (1)
       ztstep1(jstep) = timef()
     case (2)
       ztstep2(jstep) = timef()
   end select
END SUBROUTINE ectrans_set_ztstep_start

SUBROUTINE ectrans_set_ztstep_end(indx,jstep)
USE transform_driver_data_mod, ONLY : ztstep,ztstep1,ztstep2
  integer(kind=jpim), intent(in) :: indx, jstep 
  select case(indx)
     case (0)
  ztstep(jstep) = (timef() - ztstep(jstep))/1000.0_jprd
     case (1)
  ztstep1(jstep) = (timef() - ztstep1(jstep))/1000.0_jprd
     case (2)
  ztstep2(jstep) = (timef() - ztstep2(jstep))/1000.0_jprd
   end select
END SUBROUTINE ectrans_set_ztstep_end
  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

SUBROUTINE ectrans_inv_trans(nproma,lscders,luvders,ivset,ivsetsc)
USE transform_driver_data_mod, ONLY : zspvor,zspdiv,zspsc3a,zspsc2
USE transform_driver_data_mod, ONLY:  zgpuv, zgp3a, zgp2
  integer(kind=jpim), intent(in) :: nproma
  integer(kind=jpim), intent(in) :: ivset(:),ivsetsc(1)
  logical, intent(in) :: lscders
  logical, intent(in), optional :: luvders
    if(present(luvders)) then
    call inv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspvor=zspvor,                     & ! spectral vorticity
       & pspdiv=zspdiv,                     & ! spectral divergence
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  &
       & ldvorgp=.false.,                   & ! no gridpoint vorticity
       & lddivgp=.false.,                   & ! no gridpoint divergence
       & lduvder=luvders,                   &
       & kvsetuv=ivset,                     &
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgpuv=zgpuv,                       &
       & pgp3a=zgp3a)
  else
    call inv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  & ! scalar derivatives
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgp3a=zgp3a)
  endif
END SUBROUTINE ectrans_inv_trans

  !=================================================================================================
  ! While in grid point space, dump the values to disk, for debugging only
  !=================================================================================================

SUBROUTINE ectrans_dump(jstep, myproc, nproma, ngpblks, nflevg, noutdump)
USE transform_driver_data_mod, ONLY:  zgpuv, zgp3a, zgp2
  integer(kind=jpim), intent(in) :: jstep, myproc, nproma, ngpblks, nflevg, noutdump 
    ! dump a field to a binary file
    call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgp2(:,1,:),         'S', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgpuv(:,nflevg,1,:), 'U', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgpuv(:,nflevg,2,:), 'V', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgp3a(:,nflevg,1,:), 'T', noutdump)
END SUBROUTINE ectrans_dump

  !=================================================================================================
  ! Do direct transform
  !=================================================================================================


SUBROUTINE ectrans_direct_trans(nproma,nfld,ivset,ivsetsc,lvordiv)
USE transform_driver_data_mod, ONLY : zspvor,zspdiv,zspsc3a,zspsc2
USE transform_driver_data_mod, ONLY:  zgpuv, zgp3a, zgmvs
  integer(kind=jpim), intent(in) :: nproma, nfld
  integer(kind=jpim), intent(in) :: ivset(:), ivsetsc(1)
  logical, intent(in), optional :: lvordiv
    if(present(lvordiv)) then
    call dir_trans(kresol=1, kproma=nproma, &
      & pgp2=zgmvs(:,1:1,:),                &
      & pgpuv=zgpuv(:,:,1:2,:),             &
      & pgp3a=zgp3a(:,:,1:nfld,:),          &
      & pspvor=zspvor,                      &
      & pspdiv=zspdiv,                      &
      & pspsc2=zspsc2,                      &
      & pspsc3a=zspsc3a,                    &
      & kvsetuv=ivset,                      &
      & kvsetsc2=ivsetsc,                   &
      & kvsetsc3a=ivset)
  else
    call dir_trans(kresol=1, kproma=nproma, &
      & pgp2=zgmvs(:,1:1,:),                &
      & pgp3a=zgp3a(:,:,1:nfld,:),          &
      & pspsc2=zspsc2,                      &
      & pspsc3a=zspsc3a,                    &
      & kvsetsc2=ivsetsc,                   &
      & kvsetsc3a=ivset)
  endif
END SUBROUTINE ectrans_direct_trans

  !=================================================================================================
  ! Calculate timings
  !=================================================================================================

SUBROUTINE ectrans_calculate_timings(jstep)
USE transform_driver_data_mod, ONLY: ztstepmax , ztstepmin , ztstepavg 
USE transform_driver_data_mod, ONLY: ztstepmax1, ztstepmin1, ztstepavg1
USE transform_driver_data_mod, ONLY: ztstepmax2, ztstepmin2, ztstepavg2
USE transform_driver_data_mod, ONLY : ztstep, ztstep1, ztstep2
  integer(kind=jpim), intent(in) :: jstep 
  ztstep(jstep) = (timef() - ztstep(jstep))/1000.0_jprd

  ztstepavg = ztstepavg + ztstep(jstep)
  ztstepmin = min(ztstep(jstep), ztstepmin)
  ztstepmax = max(ztstep(jstep), ztstepmax)

  ztstepavg1 = ztstepavg1 + ztstep1(jstep)
  ztstepmin1 = min(ztstep1(jstep), ztstepmin1)
  ztstepmax1 = max(ztstep1(jstep), ztstepmax1)

  ztstepavg2 = ztstepavg2 + ztstep2(jstep)
  ztstepmin2 = min(ztstep2(jstep), ztstepmin2)
  ztstepmax2 = max(ztstep2(jstep), ztstepmax2)
END SUBROUTINE ectrans_calculate_timings

  !=================================================================================================
  ! Print norms
  !=================================================================================================


    SUBROUTINE ectrans_print_timestep(nout,jstep)
USE transform_driver_data_mod, ONLY : ztstep
integer(kind=jpim), intent(in) :: nout,jstep
write(nout,'("Time step ",i6," took", f8.4)') jstep, ztstep(jstep)
END SUBROUTINE ectrans_print_timestep
    SUBROUTINE ectrans_print_norms_calc(nout,jstep,myproc,nflevg)
USE transform_driver_data_mod, ONLY : znormsp,znormsp1 
USE transform_driver_data_mod, ONLY : znormvor,znormvor1 
USE transform_driver_data_mod, ONLY : znormdiv,znormdiv1 
USE transform_driver_data_mod, ONLY : znormt,znormt1 
USE transform_driver_data_mod, ONLY : ztstep
USE transform_driver_data_mod, ONLY : zerr,zmaxerr 
integer(kind=jpim), intent(in) :: nout,jstep,myproc,nflevg
integer(kind=jpim) :: ifld
    ! Surface pressure
    if (myproc == 1) then
      zmaxerr(:) = -999.0
      do ifld = 1, 1
        write(nout,*) "znormsp", znormsp
        call flush(nout)
        zerr(1) = abs(znormsp1(ifld)/znormsp(ifld) - 1.0_jprb)
        zmaxerr(1) = max(zmaxerr(1), zerr(1))
      enddo
      ! Divergence
      do ifld = 1, nflevg
        zerr(2) = abs(znormdiv1(ifld)/znormdiv(ifld) - 1.0_jprb)
        zmaxerr(2) = max(zmaxerr(2), zerr(2))
      enddo
      ! Vorticity
      do ifld = 1, nflevg
        zerr(3) = abs(znormvor1(ifld)/znormvor(ifld) - 1.0_jprb)
        zmaxerr(3) = max(zmaxerr(3),zerr(3))
      enddo
      ! Temperature
      do ifld = 1, nflevg
        zerr(4) = abs(znormt1(ifld)/znormt(ifld) - 1.0_jprb)
        zmaxerr(4) = max(zmaxerr(4), zerr(4))
      enddo
      write(nout,'("time step ",i6," took", f8.4," | zspvor max err="e10.3,&
                  & " | zspdiv max err="e10.3," | zspsc3a max err="e10.3," | zspsc2 max err="e10.3)') &
                  &  jstep, ztstep(jstep), zmaxerr(3), zmaxerr(2), zmaxerr(4), zmaxerr(1)
    endif
    END SUBROUTINE ectrans_print_norms_calc

    SUBROUTINE ectrans_print_norms_fin(nout,nflevg,verbosity)
USE transform_driver_data_mod, ONLY : znormsp,znormsp1 
USE transform_driver_data_mod, ONLY : znormvor,znormvor1 
USE transform_driver_data_mod, ONLY : znormdiv,znormdiv1 
USE transform_driver_data_mod, ONLY : znormt,znormt1 
USE transform_driver_data_mod, ONLY : zerr,zmaxerr,zmaxerrg 
integer(kind=jpim), intent(in) :: nout,nflevg
integer(kind=jpim) :: ifld
integer, intent(in)  :: verbosity
    zmaxerr(:) = -999.0
    do ifld = 1, nflevg
      zerr(3) = abs(real(znormvor1(ifld),kind=jprd)/real(znormvor(ifld),kind=jprd) - 1.0_jprd)
      zmaxerr(3) = max(zmaxerr(3), zerr(3))
      if (verbosity >= 1) then
        write(nout,'("norm zspvor( ",i4,")     = ",f20.15,"        error = ",e10.3)') ifld, znormvor1(ifld), zerr(3)
      endif
    enddo
    do ifld = 1, nflevg
      zerr(2) = abs(real(znormdiv1(ifld),kind=jprd)/real(znormdiv(ifld),kind=jprd) - 1.0d0)
      zmaxerr(2) = max(zmaxerr(2),zerr(2))
      if (verbosity >= 1) then
        write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15,"        error = ",e10.3)') ifld, znormdiv1(ifld), zerr(2)
      endif
    enddo
    do ifld = 1, nflevg
      zerr(4) = abs(real(znormt1(ifld),kind=jprd)/real(znormt(ifld),kind=jprd) - 1.0d0)
      zmaxerr(4) = max(zmaxerr(4), zerr(4))
      if (verbosity >= 1) then
        write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15,"        error = ",e10.3)') ifld, znormt1(ifld), zerr(4)
      endif
    enddo
    do ifld = 1, 1
      zerr(1) = abs(real(znormsp1(ifld),kind=jprd)/real(znormsp(ifld),kind=jprd) - 1.0d0)
      zmaxerr(1) = max(zmaxerr(1), zerr(1))
      if (verbosity >= 1) then
        write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15,"        error = ",e10.3)') ifld, znormsp1(ifld), zerr(1)
      endif
    enddo

    ! maximum error across all fields
    zmaxerrg = max(max(zmaxerr(1),zmaxerr(2)), max(zmaxerr(2), zmaxerr(3)))

    if (verbosity >= 1) write(nout,*)
    write(nout,'("max error zspvor(1:nlev,:)    = ",e10.3)') zmaxerr(3)
    write(nout,'("max error zspdiv(1:nlev,:)    = ",e10.3)') zmaxerr(2)
    write(nout,'("max error zspsc3a(1:nlev,:,1) = ",e10.3)') zmaxerr(4)
    write(nout,'("max error zspsc2(1:1,:)       = ",e10.3)') zmaxerr(1)
    write(nout,*)
    write(nout,'("max error combined =          = ",e10.3)') zmaxerrg
    write(nout,*)
    END SUBROUTINE ectrans_print_norms_fin
      ! If the maximum spectral norm error across all fields is greater than 100 times the machine
      ! epsilon, fail the test
    FUNCTION ectrans_print_norms_fails(nout,ncheck) result (ierr)
USE transform_driver_data_mod, ONLY : zmaxerrg 
      integer(kind=jpim), intent(in) :: nout,ncheck
      integer(kind=jpim) :: ierr 
      ierr = 0
      if (zmaxerrg > real(ncheck, jprb) * epsilon(1.0_jprb)) then
        write(nout, '(a)') '*******************************'
        write(nout, '(a)') 'Correctness test failed'
        write(nout, '(a,1e7.2)') 'Maximum spectral norm error = ', zmaxerrg
        write(nout, '(a,1e7.2)') 'Error tolerance = ', real(ncheck, jprb) * epsilon(1.0_jprb)
        write(nout, '(a)') '*******************************'
        ierr = 1
      endif
    END FUNCTION ectrans_print_norms_fails

    SUBROUTINE ectrans_norms_reduce(ztloop)
use mpl_module
USE transform_driver_data_mod, ONLY: ztstepmax , ztstepmin , ztstepavg 
USE transform_driver_data_mod, ONLY: ztstepmax1, ztstepmin1, ztstepavg1
USE transform_driver_data_mod, ONLY: ztstepmax2, ztstepmin2, ztstepavg2
USE transform_driver_data_mod, ONLY: ztstep, ztstep1, ztstep2
    real(kind=jprd), intent(inout) :: ztloop
  call mpl_allreduce(ztloop,     'sum', ldreprod=.false.)
  call mpl_allreduce(ztstep,     'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepavg,  'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepmax,  'max', ldreprod=.false.)
  call mpl_allreduce(ztstepmin,  'min', ldreprod=.false.)

  call mpl_allreduce(ztstep1,    'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepavg1, 'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepmax1, 'max', ldreprod=.false.)
  call mpl_allreduce(ztstepmin1, 'min', ldreprod=.false.)

  call mpl_allreduce(ztstep2,    'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepavg2, 'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepmax2, 'max', ldreprod=.false.)
  call mpl_allreduce(ztstepmin2, 'min', ldreprod=.false.)
 END SUBROUTINE ectrans_norms_reduce

 SUBROUTINE ectrans_compute_time_stats(nproc,iters)
USE transform_driver_data_mod, ONLY : ztstep, ztstep1, ztstep2
USE transform_driver_data_mod, ONLY : ztstepavg, ztstepavg1, ztstepavg2  
USE transform_driver_data_mod, ONLY : ztstepmed, ztstepmed1, ztstepmed2  
      integer(kind=jpim), intent(in) :: nproc,iters
ztstepavg = (ztstepavg/real(nproc,jprb))/real(iters,jprd)
ztstep(:) = ztstep(:)/real(nproc,jprd)

call sort(ztstep,iters)
ztstepmed = ztstep(iters/2)

ztstepavg1 = (ztstepavg1/real(nproc,jprb))/real(iters,jprd)
ztstep1(:) = ztstep1(:)/real(nproc,jprd)

call sort(ztstep1, iters)
ztstepmed1 = ztstep1(iters/2)

ztstepavg2 = (ztstepavg2/real(nproc,jprb))/real(iters,jprd)
ztstep2(:) = ztstep2(:)/real(nproc,jprd)

call sort(ztstep2,iters)
ztstepmed2 = ztstep2(iters/2)

 END SUBROUTINE ectrans_compute_time_stats

 SUBROUTINE ectrans_print_time_stats(nout,ztloop,nproc)
USE transform_driver_data_mod, ONLY : ztstepavg, ztstepavg1, ztstepavg2  
USE transform_driver_data_mod, ONLY : ztstepmax, ztstepmax1, ztstepmax2  
USE transform_driver_data_mod, ONLY : ztstepmin, ztstepmin1, ztstepmin2  
USE transform_driver_data_mod, ONLY : ztstepmed, ztstepmed1, ztstepmed2  
    real(kind=jprd), intent(inout) :: ztloop
      integer(kind=jpim), intent(in) :: nproc,nout
ztloop = ztloop/real(nproc,jprd)
write(nout,'(a)') '======= Start of time step stats ======='
write(nout,'(" ")')
write(nout,'("Inverse transforms")')
write(nout,'("------------------")')
write(nout,'("avg  (s): ",f8.4)') ztstepavg1
write(nout,'("min  (s): ",f8.4)') ztstepmin1
write(nout,'("max  (s): ",f8.4)') ztstepmax1
write(nout,'("med  (s): ",f8.4)') ztstepmed1
write(nout,'(" ")')
write(nout,'("Direct transforms")')
write(nout,'("-----------------")')
write(nout,'("avg  (s): ",f8.4)') ztstepavg2
write(nout,'("min  (s): ",f8.4)') ztstepmin2
write(nout,'("max  (s): ",f8.4)') ztstepmax2
write(nout,'("med  (s): ",f8.4)') ztstepmed2
write(nout,'(" ")')
write(nout,'("Inverse-direct transforms")')
write(nout,'("-------------------------")')
write(nout,'("avg  (s): ",f8.4)') ztstepavg
write(nout,'("min  (s): ",f8.4)') ztstepmin
write(nout,'("max  (s): ",f8.4)') ztstepmax
write(nout,'("med  (s): ",f8.4)') ztstepmed
write(nout,'("loop (s): ",f8.4)') ztloop
write(nout,'(" ")')
write(nout,'(a)') '======= End of time step stats ======='
write(nout,'(" ")')
END SUBROUTINE ectrans_print_time_stats


subroutine sort(a, n)

  real(kind=jprd), intent(inout) :: a(n)
  integer(kind=jpim), intent(in) :: n

  real(kind=jprd) :: x

  integer :: i, j

  do i = 2, n
    x = a(i)
    j = i - 1
    do while (j >= 1)
      if (a(j) <= x) exit
        a(j + 1) = a(j)
        j = j - 1
      end do
    a(j + 1) = x
  end do

end subroutine sort
!
!===================================================================================================

subroutine initialize_spectral_arrays(nsmax, zsp, sp3d)

  integer,         intent(in)    :: nsmax       ! Spectral truncation
  real(kind=jprb), intent(inout) :: zsp(:,:)    ! Surface pressure
  real(kind=jprb), intent(inout) :: sp3d(:,:,:) ! 3D fields

  integer(kind=jpim) :: nflevl
  integer(kind=jpim) :: nfield

  integer :: i, j

  nflevl = size(sp3d, 1)
  nfield = size(sp3d, 3)

  ! First initialize surface pressure
  call initialize_2d_spectral_field(nsmax, zsp(1,:))

  ! Then initialize all of the 3D fields
  do i = 1, nflevl
    do j = 1, nfield
      call initialize_2d_spectral_field(nsmax, sp3d(i,:,j))
    end do
  end do

end subroutine initialize_spectral_arrays

!===================================================================================================

subroutine initialize_2d_spectral_field(nsmax, field)

  integer,         intent(in)    :: nsmax    ! Spectral truncation
  real(kind=jprb), intent(inout) :: field(:) ! Field to initialize

  integer :: i, index, num_my_zon_wns
  integer, allocatable :: my_zon_wns(:), nasm0(:)

  ! Choose a spherical harmonic to initialize arrays
  integer :: m_num = 4  ! Zonal wavenumber
  integer :: l_num = 19  ! Total wavenumber

  ! First initialise all spectral coefficients to zero
  field(:) = 0.0

  ! Get zonal wavenumbers this rank is responsible for
  call trans_inq(knump=num_my_zon_wns)
  allocate(my_zon_wns(num_my_zon_wns))
  call trans_inq(kmyms=my_zon_wns)

  ! If rank is responsible for the chosen zonal wavenumber...
  if (any(my_zon_wns == m_num) ) then
    ! Get array of spectral array addresses (this maps (m, n=m) to array index)
    allocate(nasm0(0:nsmax))
    call trans_inq(kasm0=nasm0)

    ! Find out local array index of chosen spherical harmonic
    index = nasm0(m_num) + 2 * (l_num - m_num) + 1

    ! Set just that element to a constant value
    field(index) = 1.0
  else
    return
  end if

end subroutine initialize_2d_spectral_field

!===================================================================================================

subroutine dump_gridpoint_field(jstep, myproc, nproma, ngpblks, fld, fldchar, noutdump)

  ! Dump a 2d field to a binary file.

  integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
  integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
  integer(kind=jpim), intent(in) :: nproma ! Size of nproma
  integer(kind=jpim), intent(in) :: ngpblks ! Number of nproma blocks
  real(kind=jprb)   , intent(in) :: fld(nproma,ngpblks) ! 2D field
  character         , intent(in) :: fldchar ! Single character field identifier
  integer(kind=jpim), intent(in) :: noutdump ! Tnit number for output file

  character(len=14) :: filename = "x.xxx.xxxx.dat"

  write(filename(1:1),'(a1)') fldchar
  write(filename(3:5),'(i3.3)') jstep
  write(filename(7:10),'(i4.4)') myproc

  open(noutdump, file=filename, form="unformatted")
  write(noutdump) reshape(fld, (/ nproma*ngpblks /))
  close(noutdump)

end subroutine dump_gridpoint_field


end module transform_driver

!===================================================================================================
