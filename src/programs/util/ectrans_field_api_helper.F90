module ectrans_field_api_helper

use field_module, only:field_1rb, field_2rb, field_3rb, field_4rb
use field_factory_module
use parkind1, only: jpim, jprb, jprd
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"

implicit none

type wrapped_fields

! Set of fields for spectral transform

  class (field_3rb), pointer :: spscalar3      ! spectral scalar fields
  class (field_2rb), pointer :: spscalar2     ! spectral surfacic scalar fields
  class (field_2rb), pointer :: spscalar  ! spectral scalar fields as a single field

  class (field_2rb), pointer :: spvor, spdiv  ! spectral vorticity and divergence

  class (field_3rb), pointer :: vor, div      ! grid-point vorticity and divergence
  class (field_3rb), pointer :: u, v          ! grid-point u and v fields
  class (field_3rb), pointer :: u_ns, v_ns    ! grid-point u and derivatives

  class (field_4rb), pointer :: scalar3      ! grid-point scalar fields
  class (field_4rb), pointer :: scalar3_ew   ! grid-point scalar fields derivatives ew
  class (field_4rb), pointer :: scalar3_ns   ! grid-point scalar fields derivatives ns

  class (field_3rb), pointer :: scalar2      ! grid-point surfacic scalar fields
  class (field_3rb), pointer :: scalar2_ew   ! grid-point surfacic scalar fields derivatives ew
  class (field_3rb), pointer :: scalar2_ns   ! grid-point surfacic scalar fields derivatives ns

  class (field_3rb), pointer :: scalar   ! grid-point scalar fields as a single field
  class (field_3rb), pointer :: scalar_ew! grid-point scalar fields derivatives ew
  class (field_3rb), pointer :: scalar_ns! grid-point scalar fields derivatives ns
end type wrapped_fields

type fields_lists

! List of field lists that will be used as parameter to inv_trans_field_api and dir_trans_field_api

  type (field_basic_ptr), allocatable :: u (:), v (:)                ! grid-point u and v fields
  type (field_basic_ptr), allocatable :: scalar (:)                  ! grid-point scalar fields
  type (field_basic_ptr), allocatable :: spvor (:), spdiv (:)        ! spectral vorticity and divergence
  type (field_basic_ptr), allocatable :: vor (:), div (:)            ! grid-point vorticity and diverence
  type (field_basic_ptr), allocatable :: spscalar (:)                ! spectral scalar fields
  type (field_basic_ptr), allocatable :: u_ns (:), v_ns (:)            ! grid-point u and derivatives we
  type (field_basic_ptr), allocatable :: scalar_ns (:), scalar_ew (:)  ! grid space scalar derivatives ns and ew
  end type fields_lists

contains

function istart(ioffset,inum_fields)
  integer :: istart
  integer, intent(in):: ioffset, inum_fields
  istart = ioffset * inum_fields + 1
end function istart

function iend(ioffset,inum_fields)
  integer :: iend
  integer, intent(in):: ioffset, inum_fields
  iend = (ioffset + 1) * inum_fields
end function iend


subroutine wrap_benchmark_fields_zgp(ywflds, lvordiv, lscders, luvders,&
                                   & zspvor, zspdiv, zspscalar, zgp)

  ! Wrap the arrays given as input in field API objects

    type(wrapped_fields), intent(inout) :: ywflds
    logical, intent(in) :: lvordiv
    logical, intent(in) :: lscders
    logical, intent(in) :: luvders
    real(kind=jprb), intent(in) :: zspvor(:,:)
    real(kind=jprb), intent(in) :: zspdiv(:,:)
    real(kind=jprb), intent(in) :: zspscalar(:,:)
    real(kind=jprb), intent(in) :: zgp (:,:,:)

    integer :: inum_wind_fields,  inum_sc_2d_fields
    integer :: ioffset

    inum_wind_fields = size(zspvor,1)
    inum_sc_2d_fields = size(zspscalar,1)
 ! spectral vector fields
  call field_new(ywflds%spvor,      data=zspvor(:,:))
  call field_new(ywflds%spdiv,      data=zspdiv(:,:))
  ioffset = 0
  if (lvordiv) then
      ! In the benchmark, vorticity is not computed
      call field_new(ywflds%vor, data=zgp(:, ioffset+1:ioffset+inum_wind_fields, :))
      ioffset = ioffset + inum_wind_fields

      call field_new(ywflds%div, data=zgp(:, ioffset+1:ioffset+inum_wind_fields, :))
      ioffset = ioffset + inum_wind_fields
  endif

  call field_new(ywflds%u, data=zgp(:, ioffset+1:ioffset+inum_wind_fields, :))
  ioffset = ioffset + inum_wind_fields

  call field_new(ywflds%v, data=zgp(:, ioffset+1:ioffset+inum_wind_fields, :))
  ioffset = ioffset + inum_wind_fields

  ! grid-point vector derivatives
  if (luvders) then
     call field_new(ywflds%u_ns, data=zgp(:, ioffset+1:ioffset+inum_wind_fields, :))
     ioffset = ioffset + inum_wind_fields
     call field_new(ywflds%v_ns, data=zgp(:, ioffset+1:ioffset+inum_wind_fields, :))
     ioffset = ioffset + inum_wind_fields
  endif

  ! spectral scalar fields
  if (inum_sc_2d_fields > 0) then! spectral surfacic scalar fields
     call field_new(ywflds%spscalar, data=zspscalar(:,:))

     call field_new(ywflds%scalar, data=zgp(:,ioffset+1:ioffset+inum_sc_2d_fields,:))
     ioffset = ioffset + inum_sc_2d_fields
     if (lscders) then

        call field_new(ywflds%scalar_ns, data=zgp(:,ioffset+1:ioffset+inum_sc_2d_fields,:))
        ioffset = ioffset + inum_sc_2d_fields

        call field_new(ywflds%scalar_ew, data=zgp(:,ioffset+1:ioffset+inum_sc_2d_fields,:))
     endif
  endif
end subroutine wrap_benchmark_fields_zgp


subroutine wrap_benchmark_fields(ywflds, lvordiv, lscders, luvders,&
                               & zspvor, zspdiv, zspsc3a, zspsc2, zgpuv, zgp3a, zgp2)

  ! Wrap the arrays given as input in field API objects

    type(wrapped_fields), intent(inout) :: ywflds
    logical, intent(in) :: lvordiv
    logical, intent(in) :: lscders
    logical, intent(in) :: luvders
    real(kind=jprb), intent(in) :: zspvor(:,:)
    real(kind=jprb), intent(in) :: zspdiv(:,:)
    real(kind=jprb), intent(in) :: zspsc3a(:,:,:)
    real(kind=jprb), intent(in) :: zspsc2(:,:)
    real(kind=jprb), intent(in) :: zgpuv (:,:,:,:)
    real(kind=jprb), intent(in) :: zgp3a(:,:,:,:)
    real(kind=jprb), intent(in) :: zgp2(:,:,:)

    integer :: inum_wind_fields, inum_sc_3d_fields, inum_sc_2d_fields
    integer :: ioffset

    inum_wind_fields = 1
    inum_sc_3d_fields = size(zspsc3a,3)
    inum_sc_2d_fields = size(zspsc2,1)
 ! spectral vector fields
  call field_new(ywflds%spvor,      data=zspvor(:,:))
  call field_new(ywflds%spdiv,      data=zspdiv(:,:))

  ioffset = 0
  if (lvordiv) then
      ! In the benchmark, vorticity is not computed
      call field_new(ywflds%vor, data=zgpuv(:,:, istart(ioffset,inum_wind_fields), :))
      ioffset = ioffset + 1
      call field_new(ywflds%div, data=zgpuv(:,:, istart(ioffset,inum_wind_fields), :))
      ioffset = ioffset + 1
  endif

  call field_new(ywflds%u, data=zgpuv(:,:, istart(ioffset,inum_wind_fields), :))
  ioffset = ioffset + 1
  call field_new(ywflds%v, data=zgpuv(:,:, istart(ioffset,inum_wind_fields), :))
  ioffset = ioffset + 1


  ! grid-point vector derivatives
  if (luvders) then
      call field_new(ywflds%u_ns, data=zgpuv(:,:, istart(ioffset,inum_wind_fields), :))
     ioffset = ioffset + 1
     call field_new(ywflds%v_ns, data=zgpuv(:,:, istart(ioffset,inum_wind_fields), :))
  endif


  ! spectral scalar fields
  if (inum_sc_3d_fields > 0) then
     call field_new(ywflds%spscalar3, data=zspsc3a(:,:,:))

     call field_new(ywflds%scalar3,  data=zgp3a(:,:,1:inum_sc_3d_fields,:))
     if (lscders) then
         ioffset = 1
         ! grid-point surfacic scalar derivatives fields
        call field_new(ywflds%scalar3_ns,  data=zgp3a(:,:,istart(ioffset,inum_sc_3d_fields): iend(ioffset,inum_sc_3d_fields),:))
        ioffset = ioffset + 1
        call field_new(ywflds%scalar3_ew,  data=zgp3a(:,:,istart(ioffset,inum_sc_3d_fields): iend(ioffset,inum_sc_3d_fields),:))
    endif
   endif

  if (inum_sc_2d_fields > 0) then! spectral surfacic scalar fields
     call field_new(ywflds%spscalar2, data=zspsc2(:,:))

     call field_new(ywflds%scalar2,   data=zgp2(:,1:inum_sc_2d_fields,:))
!     call field_new(ywflds%scalar2,   data=zgp2(:,1:1,:))
     if (lscders) then
      ioffset = 1
      call field_new(ywflds%scalar2_ns, data=zgp2(:,istart(ioffset,inum_sc_2d_fields): iend(ioffset,inum_sc_2d_fields),:))
      ioffset = ioffset + 1
      call field_new(ywflds%scalar2_ew, data=zgp2(:,istart(ioffset,inum_sc_2d_fields): iend(ioffset,inum_sc_2d_fields),:))
     endif
  endif

end subroutine wrap_benchmark_fields

subroutine create_fields_lists(ywflds,ylf, kvsetuv, kvsetsc,kvsetsc2)

  ! Create field lists in ylf from field API objects in ywflds

  type(wrapped_fields), intent(in) :: ywflds       !input fields api objects
  type(fields_lists), intent(inout), target :: ylf ! output field lists
  integer(kind=jpim), optional, intent(in) :: kvsetuv(:)     ! 'b-set' for vector fields
  integer(kind=jpim), optional, intent(in) :: kvsetsc(:)     ! 'b-set' for scalar fields
  integer(kind=jpim), optional, intent(in) :: kvsetsc2(:)    ! 'b-set' for surfacic fields

  if(associated(ywflds%spvor)) ylf%spvor=[b(ywflds%spvor,'spvor',kvsetuv)]

  if(associated(ywflds%spdiv)) ylf%spdiv= [b(ywflds%spdiv,'spdiv',kvsetuv)]

  if(associated(ywflds%u)) ylf%u = [b(ywflds%u,'u',kvsetuv)]
  if(associated(ywflds%v)) ylf%v = [b(ywflds%v,'v',kvsetuv)]

  if(associated(ywflds%u_ns)) ylf%u_ns=[b(ywflds%u_ns,'u_ns', kvsetuv)]
  if(associated(ywflds%v_ns)) ylf%v_ns=[b(ywflds%v_ns,'v_ns', kvsetuv)]

  if(associated(ywflds%vor))  ylf%vor = [b(ywflds%vor,'vor', kvsetuv)]
  if(associated(ywflds%div))  ylf%div = [b(ywflds%div,'div', kvsetuv)]

  if (associated(ywflds%spscalar)) then
    ylf%spscalar = [b(ywflds%spscalar,'spscalar',kvsetsc)]
  else if (associated(ywflds%spscalar3) .and. associated(ywflds%spscalar2) ) then
    ylf%spscalar = [b(ywflds%spscalar3,'spscalar3',kvsetsc), b(ywflds%spscalar2,'spscalar2',kvsetsc2)]
  else if (associated(ywflds%spscalar3)) then
    ylf%spscalar = [b(ywflds%spscalar3,'spscalar3',kvsetsc)]
  else if (associated(ywflds%spscalar2)) then
    ylf%spscalar = [b(ywflds%spscalar2,'spscalar2',kvsetsc2)]
  endif

  if (associated(ywflds%scalar)) then
    ylf%scalar = [b(ywflds%scalar,'scalar', kvsetsc)]
  else if (associated(ywflds%scalar3) .and. associated(ywflds%scalar2) ) then
    ylf%scalar = [b(ywflds%scalar3,'scalar', kvsetsc), b(ywflds%scalar2,'scalar2', kvsetsc2)]
  else if (associated(ywflds%scalar3)) then
    ylf%scalar = [b(ywflds%scalar3,'scalar', kvsetsc)]
  else if (associated(ywflds%scalar2)) then
    ylf%scalar = [b(ywflds%scalar2,'scalar2', kvsetsc2)]
  endif

  if (associated(ywflds%scalar_ns)) then
    ylf%scalar_ns = [b(ywflds%scalar_ns,'scalar_ns', kvsetsc)]
  else if (associated(ywflds%scalar3_ns) .and. associated(ywflds%scalar2_ns) ) then
    ylf%scalar_ns = [b(ywflds%scalar3_ns,'scalar3_ns', kvsetsc), b(ywflds%scalar2_ns,'scalar2_ns', kvsetsc2)]
  else if (associated(ywflds%scalar3_ns)) then
    ylf%scalar_ns = [b(ywflds%scalar3_ns,'scalar3_ns', kvsetsc)]
  else if (associated(ywflds%scalar2_ns)) then
    ylf%scalar_ns = [b(ywflds%scalar2_ns,'scalar2_ns', kvsetsc2)]
  endif

  if (associated(ywflds%scalar_ew)) then
    ylf%scalar_ew = [b(ywflds%scalar_ew,'scalar_ew', kvsetsc)]
  else if (associated(ywflds%scalar3_ew) .and. associated(ywflds%scalar2_ew) ) then
    ylf%scalar_ew = [b(ywflds%scalar3_ew,'scalar3_ew', kvsetsc), b(ywflds%scalar2_ew,'scalar2_ew', kvsetsc2)]
  else if (associated(ywflds%scalar3_ew)) then
    ylf%scalar_ew = [b(ywflds%scalar3_ew,'scalar3_ew', kvsetsc)]
  else if (associated(ywflds%scalar2_ew)) then
    ylf%scalar_ew = [b(ywflds%scalar2_ew,'scalar2_ew', kvsetsc2)]
  endif
 end subroutine create_fields_lists

 subroutine delete_wrapped_fields(ywflds)

  ! Delete  all fields in ywflds

  type(wrapped_fields), intent(inout) :: ywflds

  if(associated(ywflds%spvor))      call field_delete(ywflds%spvor)
  if(associated(ywflds%spdiv))      call field_delete(ywflds%spdiv)
  if(associated(ywflds%spscalar))   call field_delete(ywflds%spscalar)
  if(associated(ywflds%spscalar3))  call field_delete(ywflds%spscalar3)
  if(associated(ywflds%spscalar2))  call field_delete(ywflds%spscalar2)

  if(associated(ywflds%u))          call field_delete(ywflds%u)
  if(associated(ywflds%v))          call field_delete(ywflds%v)
  if(associated(ywflds%u_ns))       call field_delete(ywflds%u_ns)
  if(associated(ywflds%v_ns))       call field_delete(ywflds%v_ns)

  if(associated(ywflds%vor))        call field_delete(ywflds%vor)
  if(associated(ywflds%div))        call field_delete(ywflds%div)

  if(associated(ywflds%scalar3))     call field_delete(ywflds%scalar3)
  if(associated(ywflds%scalar3_ew))  call field_delete(ywflds%scalar3_ew)
  if(associated(ywflds%scalar3_ns))  call field_delete(ywflds%scalar3_ns)

  if(associated(ywflds%scalar2))    call field_delete(ywflds%scalar2)
  if(associated(ywflds%scalar2_ew)) call field_delete(ywflds%scalar2_ew)
  if(associated(ywflds%scalar2_ns)) call field_delete(ywflds%scalar2_ns )

  if(associated(ywflds%scalar))    call field_delete(ywflds%scalar)
  if(associated(ywflds%scalar_ew)) call field_delete(ywflds%scalar_ew)
  if(associated(ywflds%scalar_ns)) call field_delete(ywflds%scalar_ns )


end subroutine delete_wrapped_fields

subroutine delete_fields_lists(yfl)

  ! Delete  all field lists in yfl

  type(fields_lists), intent(inout) ::yfl
  if (allocated(yfl%u))         deallocate(yfl%u)
  if (allocated(yfl%v))         deallocate(yfl%v)
  if (allocated(yfl%scalar))    deallocate(yfl%scalar)
  if (allocated(yfl%spscalar))  deallocate(yfl%spscalar)
  if (allocated(yfl%spvor))     deallocate(yfl%spvor)
  if (allocated(yfl%spdiv))     deallocate(yfl%spdiv)
  if (allocated(yfl%vor))       deallocate(yfl%vor)
  if (allocated(yfl%div))       deallocate(yfl%div)
  if (allocated(yfl%u_ns))      deallocate(yfl%u_ns)
  if (allocated(yfl%v_ns))      deallocate(yfl%v_ns)
  if (allocated(yfl%scalar_ns)) deallocate(yfl%scalar_ns)
  if (allocated(yfl%scalar_ew)) deallocate(yfl%scalar_ew)
end subroutine delete_fields_lists

subroutine synchost_rdonly_wrapped_fields(ywflds)

  ! Synchronize all field lists on host readonly

  type(wrapped_fields),intent(inout) ::ywflds

  if (associated(ywflds%spvor))      call ywflds%spvor%sync_host_rdonly()
  if (associated(ywflds%spdiv))      call ywflds%spdiv%sync_host_rdonly()
  if (associated(ywflds%spscalar))   call ywflds%spscalar%sync_host_rdonly()
  if (associated(ywflds%spscalar3))  call ywflds%spscalar3%sync_host_rdonly()
  if (associated(ywflds%spscalar2))  call ywflds%spscalar2%sync_host_rdonly()
  if (associated(ywflds%u))          call ywflds%u%sync_host_rdonly()
  if (associated(ywflds%v))          call ywflds%v%sync_host_rdonly()
  if (associated(ywflds%u_ns))       call ywflds%u_ns%sync_host_rdonly()
  if (associated(ywflds%v_ns))       call ywflds%v_ns%sync_host_rdonly()
  if (associated(ywflds%vor))        call ywflds%vor%sync_host_rdonly()
  if (associated(ywflds%div))        call ywflds%div%sync_host_rdonly()
  if (associated(ywflds%scalar3))     call ywflds%scalar3%sync_host_rdonly()
  if (associated(ywflds%scalar3_ew))  call ywflds%scalar3_ew%sync_host_rdonly()
  if (associated(ywflds%scalar3_ns))  call ywflds%scalar3_ns%sync_host_rdonly()
  if (associated(ywflds%scalar2))    call ywflds%scalar2%sync_host_rdonly()
  if (associated(ywflds%scalar2_ew)) call ywflds%scalar2_ew%sync_host_rdonly()
  if (associated(ywflds%scalar2_ns)) call ywflds%scalar2_ns%sync_host_rdonly()
  if (associated(ywflds%scalar))    call ywflds%scalar%sync_host_rdonly()
  if (associated(ywflds%scalar_ew)) call ywflds%scalar_ew%sync_host_rdonly()
  if (associated(ywflds%scalar_ns)) call ywflds%scalar_ns%sync_host_rdonly()


end subroutine synchost_rdonly_wrapped_fields

subroutine synchost_rdwr_wrapped_fields(ywflds)

  ! Synchronize all field lists on host read/write

  type(wrapped_fields),intent(inout) ::ywflds

  if (associated(ywflds%spvor))      call ywflds%spvor%sync_host_rdwr()
  if (associated(ywflds%spdiv))      call ywflds%spdiv%sync_host_rdwr()
  if (associated(ywflds%spscalar))   call ywflds%spscalar%sync_host_rdwr()
  if (associated(ywflds%spscalar3))  call ywflds%spscalar3%sync_host_rdwr()
  if (associated(ywflds%spscalar2))  call ywflds%spscalar2%sync_host_rdwr()
  if (associated(ywflds%u))          call ywflds%u%sync_host_rdwr()
  if (associated(ywflds%v))          call ywflds%v%sync_host_rdwr()
  if (associated(ywflds%u_ns))       call ywflds%u_ns%sync_host_rdwr()
  if (associated(ywflds%v_ns))       call ywflds%v_ns%sync_host_rdwr()
  if (associated(ywflds%vor))        call ywflds%vor%sync_host_rdwr()
  if (associated(ywflds%div))        call ywflds%div%sync_host_rdwr()
  if (associated(ywflds%scalar3))     call ywflds%scalar3%sync_host_rdwr()
  if (associated(ywflds%scalar3_ew))  call ywflds%scalar3_ew%sync_host_rdwr()
  if (associated(ywflds%scalar3_ns))  call ywflds%scalar3_ns%sync_host_rdwr()
  if (associated(ywflds%scalar2))    call ywflds%scalar2%sync_host_rdwr()
  if (associated(ywflds%scalar2_ew)) call ywflds%scalar2_ew%sync_host_rdwr()
  if (associated(ywflds%scalar2_ns)) call ywflds%scalar2_ns%sync_host_rdwr()
  if (associated(ywflds%scalar))    call ywflds%scalar%sync_host_rdwr()
  if (associated(ywflds%scalar_ew)) call ywflds%scalar_ew%sync_host_rdwr()
  if (associated(ywflds%scalar_ns)) call ywflds%scalar_ns%sync_host_rdwr()

end subroutine synchost_rdwr_wrapped_fields

subroutine nullify_wrapped_fields(ywflds)

    ! Nullify all pointers in ywflds

  type(wrapped_fields), intent(inout) :: ywflds

  nullify(ywflds%spvor)
  nullify(ywflds%spdiv)
  nullify(ywflds%spscalar)
  nullify(ywflds%spscalar3)
  nullify(ywflds%spscalar2)

  nullify(ywflds%u)
  nullify(ywflds%v)
  nullify(ywflds%u_ns)
  nullify(ywflds%v_ns)
  nullify(ywflds%vor)
  nullify(ywflds%div)

  nullify(ywflds%scalar3)
  nullify(ywflds%scalar3_ew)
  nullify(ywflds%scalar3_ns)

  nullify(ywflds%scalar2)
  nullify(ywflds%scalar2_ew)
  nullify(ywflds%scalar2_ns)

  nullify(ywflds%scalar)
  nullify(ywflds%scalar_ew)
  nullify(ywflds%scalar_ns)
end subroutine nullify_wrapped_fields


subroutine output_wrapped_fields(nout, ywflds)

  ! output the adress of all data fields in ywflds

  type(wrapped_fields), intent(in) :: ywflds
  integer(kind=jpim), intent(in) :: nout

  write(nout,*) "ywflds%spvor", loc(ywflds%spvor)
  write(nout,*) "ywflds%spdiv", loc(ywflds%spdiv)
  write(nout,*) "ywflds%spscalar", loc(ywflds%spscalar)
  write(nout,*) "ywflds%spscalar3", loc(ywflds%spscalar3)
  write(nout,*) "ywflds%spscalar2", loc(ywflds%spscalar2)

  write(nout,*) "ywflds%u", loc(ywflds%u)
  write(nout,*) "ywflds%v", loc(ywflds%v)
  write(nout,*) "ywflds%u_ns", loc(ywflds%u_ns)
  write(nout,*) "ywflds%v_ns", loc(ywflds%v_ns)
  write(nout,*) "ywflds%vor", loc(ywflds%vor)
  write(nout,*) "ywflds%div", loc(ywflds%div)

  write(nout,*) "ywflds%scalar3", loc(ywflds%scalar3)
  write(nout,*) "ywflds%scalar3_ew", loc(ywflds%scalar3_ew)
  write(nout,*) "ywflds%scalar3_ns", loc(ywflds%scalar3_ns)

  write(nout,*) "ywflds%scalar2", loc(ywflds%scalar2)
  write(nout,*) "ywflds%scalar2_ew", loc(ywflds%scalar2_ew)
  write(nout,*) "ywflds%scalar2_ns", loc(ywflds%scalar2_ns)

  write(nout,*) "ywflds%scalar", loc(ywflds%scalar)
  write(nout,*) "ywflds%scalar_ew", loc(ywflds%scalar_ew)
  write(nout,*) "ywflds%scalar_ns", loc(ywflds%scalar_ns)
end subroutine output_wrapped_fields

subroutine output_fields_lists(nout,yfl)

  ! output the size of all field lists in yfl

  integer(kind=jpim), intent(in) :: nout
  type(fields_lists), intent(in) :: yfl

  if (allocated(yfl%u)) write(nout,*) "yfl%u", size(yfl%u)
  if (allocated(yfl%v)) write(nout,*) "yfl%v", size(yfl%v)
  if (allocated(yfl%scalar)) write(nout,*) "yfl%scalar", size(yfl%scalar)
  if (allocated(yfl%spscalar)) write(nout,*) "yfl%spscalar", size(yfl%spscalar)
  if (allocated(yfl%spvor)) write(nout,*) "yfl%spvor", size(yfl%spvor)
  if (allocated(yfl%spdiv)) write(nout,*) "yfl%spdiv", size(yfl%spdiv)
  if (allocated(yfl%vor)) write(nout,*) "yfl%vor", size(yfl%vor)
  if (allocated(yfl%div)) write(nout,*) "yfl%div", size(yfl%div)
  if (allocated(yfl%u_ns)) write(nout,*) "yfl%u_ns", size(yfl%u_ns)
  if (allocated(yfl%v_ns)) write(nout,*) "yfl%v_ns", size(yfl%v_ns)
  if (allocated(yfl%scalar_ns)) write(nout,*) "yfl%scalar_ns", size(yfl%scalar_ns)
  if (allocated(yfl%scalar_ew)) write(nout,*) "yfl%scalar_ew", size(yfl%scalar_ew)
end subroutine output_fields_lists


end module ectrans_field_api_helper
