module ectrans_field_api_helper

use field_module, only:field_1rb, field_2rb, field_3rb, field_4rb
use field_factory_module
use parkind1, only: jpim, jprb, jprd
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"

implicit none

type wrapped_fields

! Set of fields for spectral transform

  class (field_3rb), pointer :: spscalar      ! spectral scalar fields
  class (field_2rb), pointer :: spscalar2     ! spectral surfacic scalar fields
  class (field_2rb), pointer :: spvor, spdiv  ! spectral vorticity and divergence

  class (field_3rb), pointer :: vor, div      ! grid-point vorticity and divergence
  class (field_3rb), pointer :: u, v          ! grid-point u and v fields
  class (field_4rb), pointer :: u_ns, v_ns    ! grid-point u and derivatives

  class (field_4rb), pointer :: scalar       ! grid-point scalar fields
  class (field_4rb), pointer :: scalar_ew    ! grid-point scalar fields derivatives ew
  class (field_4rb), pointer :: scalar_ns    ! grid-point scalar fields derivatives ns

  class (field_3rb), pointer :: scalar2      ! grid-point surfacic scalar fields
  class (field_3rb), pointer :: scalar2_ew   ! grid-point surfacic scalar fields derivatives ew
  class (field_3rb), pointer :: scalar2_ns   ! grid-point surfacic scalar fields derivatives ns
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

subroutine wrap_benchmark_fields(ywflds, lvordiv, lscders, luvders,&
                               & sp3d, spc2, zgmv, zgmvs, zgp2,    &
                               & jbegin_uv,jend_uv,                &
                               & jbegin_sc,jend_sc,                &
                               & jbegin_scder_ns, jend_scder_ns,   &
                               & jbegin_scder_ew, jend_scder_ew,   &
                               & jbegin_uder_ew, jend_uder_ew,     &
                               & jbegin_vder_ew, jend_vder_ew)

  ! Wrap the arrays given as input in field API objects

    type(wrapped_fields), intent(inout) :: ywflds
    logical, intent(in) :: lvordiv
    logical, intent(in) :: lscders
    logical, intent(in) :: luvders
    real(kind=jprb), intent(in) :: sp3d(:,:,:)
    real(kind=jprb), intent(in) :: spc2(:,:)
    real(kind=jprb), intent(in) :: zgmv(:,:,:,:)
    real(kind=jprb), intent(in) :: zgmvs(:,:,:)
    real(kind=jprb), intent(in) :: zgp2 (:,:,:)
    integer(kind=jpim), intent(in) :: jbegin_uv
    integer(kind=jpim), intent(in) :: jend_uv
    integer(kind=jpim), intent(in) :: jbegin_sc
    integer(kind=jpim), intent(in) :: jend_sc
    integer(kind=jpim), intent(in) :: jbegin_scder_ns
    integer(kind=jpim), intent(in) :: jend_scder_ns
    integer(kind=jpim), intent(in) :: jbegin_scder_ew
    integer(kind=jpim), intent(in) :: jend_scder_ew
    integer(kind=jpim), intent(in) :: jbegin_uder_ew
    integer(kind=jpim), intent(in) :: jend_uder_ew
    integer(kind=jpim), intent(in) :: jbegin_vder_ew
    integer(kind=jpim), intent(in) :: jend_vder_ew

  if (lvordiv) then
      if (jbegin_uv>0 )      call field_new(ywflds%u, data=zgmv(:,:,jbegin_uv,:))
      if (jend_uv>0 )        call field_new(ywflds%v, data=zgmv(:,:,jend_uv,:))

      ! In the benchmark, vorticity is not computed
      !call field_new(ywflds%vor,         data=zgmv(:,:,jbegin_uv,:))
      !call field_new(ywflds%div,         data=zgmv(:,:,jend_uv,:))
  endif

  ! spectral vector fields
  if (size(sp3d,3) >=1 ) call field_new(ywflds%spvor,      data=sp3d(:,:,1))
  if (size(sp3d,3) >=2 ) call field_new(ywflds%spdiv,      data=sp3d(:,:,2))

  ! spectral scalar fields
  if (size(sp3d,3) >=3 ) call field_new(ywflds%spscalar,  data=sp3d(:,:,3:))
  if (size(spc2,2) >=1 ) call field_new(ywflds%spscalar2, data=spc2(:,:))

  ! spectral surfacic scalar fields
  if (size(zgmvs,2)>=1)  call field_new(ywflds%scalar2,   data=zgmvs(:,1:1,:))

  ! grid-point vector derivatives
  if (luvders) then
     if (jend_uder_ew>0 .and. jend_uder_ew>=jbegin_uder_ew ) call field_new(ywflds%u_ns, data=zgmv(:,:,jbegin_uder_ew:jend_uder_ew,:))
     if (jend_vder_ew>0 .and. jend_vder_ew>=jbegin_vder_ew ) call field_new(ywflds%v_ns, data=zgmv(:,:,jbegin_vder_ew:jend_vder_ew,:))
  endif
 
  ! grid-point scalar fields
  if (jend_sc>0 .and. jend_sc>=jbegin_sc ) call field_new(ywflds%scalar,  data=zgmv(:,:,jbegin_sc:jend_sc,:))

  ! grid-point scalar derivatives fields
  if (lscders) then
    if (jend_scder_ew>0 .and. jend_scder_ew>=jbegin_scder_ew ) call field_new(ywflds%scalar_ew,  data=zgmv(:,:,jbegin_scder_ew:jend_scder_ew,:))
    if (jend_scder_ns>0 .and. jend_scder_ns>=jbegin_scder_ns ) call field_new(ywflds%scalar_ns,  data=zgmv(:,:,jbegin_scder_ns:jend_scder_ns,:))

    ! grid-point surfacic scalar derivatives fields
    if (size(zgmvs,2)>=2)     call field_new(ywflds%scalar2_ns, data=zgmvs(:,2:2,:))
    if (size(zgmvs,2)>=3)     call field_new(ywflds%scalar2_ew, data=zgmvs(:,3:3,:))
  endif
end subroutine wrap_benchmark_fields

subroutine create_fields_lists(ywflds,ylf, nbsetlev,nbsetsc2)

  ! Create field lists in ylf from field API objects in ywflds

  type(wrapped_fields), intent(in) :: ywflds       !input fields api objects
  type(fields_lists), intent(inout), target :: ylf ! output field lists
  integer(kind=jpim), intent(in) :: nbsetlev(:)    ! 'b-set' for vector fields
  integer(kind=jpim), intent(in) :: nbsetsc2(:)    ! 'b-set' for surfacic fields

  if(associated(ywflds%spvor)) ylf%spvor=[b(ywflds%spvor,'spvor',nbsetlev)]

  if(associated(ywflds%spdiv)) ylf%spdiv= [b(ywflds%spdiv,'spdiv',nbsetlev)]

  if(associated(ywflds%u)) ylf%u = [b(ywflds%u,'u',nbsetlev)]
  if(associated(ywflds%v)) ylf%v = [b(ywflds%v,'v',nbsetlev)]

  if(associated(ywflds%u_ns)) ylf%u_ns=[b(ywflds%u_ns,'u_ns', nbsetlev)]
  if(associated(ywflds%v_ns)) ylf%v_ns=[b(ywflds%v_ns,'v_ns', nbsetlev)]

  if(associated(ywflds%vor))  ylf%vor = [b(ywflds%vor,'vor', nbsetlev)]
  if(associated(ywflds%div))  ylf%div = [b(ywflds%div,'div', nbsetlev)]

  if (associated(ywflds%spscalar) .and. associated(ywflds%spscalar2) ) then
    ylf%spscalar = [b(ywflds%spscalar,'spscalar',nbsetlev), b(ywflds%spscalar2,'spscalar2',nbsetsc2)]
  else if (associated(ywflds%spscalar)) then
    ylf%spscalar = [b(ywflds%spscalar,'spscalar',nbsetlev)]
  else if (associated(ywflds%spscalar2)) then
    ylf%spscalar = [b(ywflds%spscalar2,'spscalar2',nbsetsc2)]
  endif

  if (associated(ywflds%scalar) .and. associated(ywflds%scalar2) ) then
    ylf%scalar = [b(ywflds%scalar,'scalar', nbsetlev), b(ywflds%scalar2,'scalar2', nbsetsc2)]
  else if (associated(ywflds%scalar)) then
    ylf%scalar = [b(ywflds%scalar,'scalar', nbsetlev)]
  else if (associated(ywflds%scalar2)) then
    ylf%scalar = [b(ywflds%scalar2,'scalar2', nbsetsc2)]
  endif

  if (associated(ywflds%scalar_ns) .and. associated(ywflds%scalar2_ns) ) then
    ylf%scalar_ns = [b(ywflds%scalar_ns,'scalar_ns', nbsetlev), b(ywflds%scalar2_ns,'scalar2_ns', nbsetsc2)]
  else if (associated(ywflds%scalar_ns)) then
    ylf%scalar_ns = [b(ywflds%scalar_ns,'scalar_ns', nbsetlev)]
  else if (associated(ywflds%scalar2_ns)) then
    ylf%scalar_ns = [b(ywflds%scalar2_ns,'scalar2_ns', nbsetsc2)]
  endif

  if (associated(ywflds%scalar_ew) .and. associated(ywflds%scalar2_ew) ) then
    ylf%scalar_ew = [b(ywflds%scalar_ew,'scalar_ew', nbsetlev), b(ywflds%scalar2_ew,'scalar2_ew', nbsetsc2)]
  else if (associated(ywflds%scalar_ew)) then
    ylf%scalar_ew = [b(ywflds%scalar_ew,'scalar_ew', nbsetlev)]
  else if (associated(ywflds%scalar2_ew)) then
    ylf%scalar_ew = [b(ywflds%scalar2_ew,'scalar2_ew', nbsetsc2)]
  endif
 end subroutine create_fields_lists

 subroutine delete_wrapped_fields(ywflds)

  ! Delete  all fields in ywflds

  type(wrapped_fields), intent(inout) :: ywflds

  if(associated(ywflds%spvor))      call field_delete(ywflds%spvor)
  if(associated(ywflds%spdiv))      call field_delete(ywflds%spdiv)
  if(associated(ywflds%spscalar))   call field_delete(ywflds%spscalar)
  if(associated(ywflds%spscalar2))  call field_delete(ywflds%spscalar2)

  if(associated(ywflds%u))          call field_delete(ywflds%u)
  if(associated(ywflds%v))          call field_delete(ywflds%v)
  if(associated(ywflds%u_ns))       call field_delete(ywflds%u_ns)
  if(associated(ywflds%v_ns))       call field_delete(ywflds%v_ns)
  if(associated(ywflds%scalar))     call field_delete(ywflds%scalar)
  if(associated(ywflds%scalar_ew))  call field_delete(ywflds%scalar_ew)
  if(associated(ywflds%scalar_ns))  call field_delete(ywflds%scalar_ns)
  if(associated(ywflds%vor))        call field_delete(ywflds%vor)
  if(associated(ywflds%div))        call field_delete(ywflds%div)

  if(associated(ywflds%scalar2))    call field_delete(ywflds%scalar2)
  if(associated(ywflds%scalar2_ew)) call field_delete(ywflds%scalar2_ew)
  if(associated(ywflds%scalar2_ns)) call field_delete(ywflds%scalar2_ns )

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
  if (associated(ywflds%spscalar2))  call ywflds%spscalar2%sync_host_rdonly()
  if (associated(ywflds%u))          call ywflds%u%sync_host_rdonly()
  if (associated(ywflds%v))          call ywflds%v%sync_host_rdonly()
  if (associated(ywflds%u_ns))       call ywflds%u_ns%sync_host_rdonly()
  if (associated(ywflds%v_ns))       call ywflds%v_ns%sync_host_rdonly()
  if (associated(ywflds%scalar))     call ywflds%scalar%sync_host_rdonly()
  if (associated(ywflds%scalar_ew))  call ywflds%scalar_ew%sync_host_rdonly()
  if (associated(ywflds%scalar_ns))  call ywflds%scalar_ns%sync_host_rdonly()
  if (associated(ywflds%vor))        call ywflds%vor%sync_host_rdonly()
  if (associated(ywflds%div))        call ywflds%div%sync_host_rdonly()
  if (associated(ywflds%scalar2))    call ywflds%scalar2%sync_host_rdonly()
  if (associated(ywflds%scalar2_ew)) call ywflds%scalar2_ew%sync_host_rdonly()
  if (associated(ywflds%scalar2_ns)) call ywflds%scalar2_ns%sync_host_rdonly()
  
end subroutine synchost_rdonly_wrapped_fields

subroutine synchost_rdwr_wrapped_fields(ywflds)

  ! Synchronize all field lists on host read/write

  type(wrapped_fields),intent(inout) ::ywflds
  
  if (associated(ywflds%spvor))      call ywflds%spvor%sync_host_rdwr()
  if (associated(ywflds%spdiv))      call ywflds%spdiv%sync_host_rdwr()
  if (associated(ywflds%spscalar))   call ywflds%spscalar%sync_host_rdwr()
  if (associated(ywflds%spscalar2))  call ywflds%spscalar2%sync_host_rdwr()
  if (associated(ywflds%u))          call ywflds%u%sync_host_rdwr()
  if (associated(ywflds%v))          call ywflds%v%sync_host_rdwr()
  if (associated(ywflds%u_ns))       call ywflds%u_ns%sync_host_rdwr()
  if (associated(ywflds%v_ns))       call ywflds%v_ns%sync_host_rdwr()
  if (associated(ywflds%scalar))     call ywflds%scalar%sync_host_rdwr()
  if (associated(ywflds%scalar_ew))  call ywflds%scalar_ew%sync_host_rdwr()
  if (associated(ywflds%scalar_ns))  call ywflds%scalar_ns%sync_host_rdwr()
  if (associated(ywflds%vor))        call ywflds%vor%sync_host_rdwr()
  if (associated(ywflds%div))        call ywflds%div%sync_host_rdwr()
  if (associated(ywflds%scalar2))    call ywflds%scalar2%sync_host_rdwr()
  if (associated(ywflds%scalar2_ew)) call ywflds%scalar2_ew%sync_host_rdwr()
  if (associated(ywflds%scalar2_ns)) call ywflds%scalar2_ns%sync_host_rdwr()

end subroutine synchost_rdwr_wrapped_fields

subroutine nullify_wrapped_fields(ywflds)

    ! Nullify all pointers in ywflds

  type(wrapped_fields), intent(inout) :: ywflds

  nullify(ywflds%spvor)
  nullify(ywflds%spdiv)
  nullify(ywflds%spscalar)
  nullify(ywflds%spscalar2)

  nullify(ywflds%u)
  nullify(ywflds%v)
  nullify(ywflds%u_ns)
  nullify(ywflds%v_ns)
  nullify(ywflds%scalar)
  nullify(ywflds%scalar_ew)
  nullify(ywflds%scalar_ns)
  nullify(ywflds%vor)
  nullify(ywflds%div)

  nullify(ywflds%scalar2)
  nullify(ywflds%scalar2_ew)
  nullify(ywflds%scalar2_ns)
end subroutine nullify_wrapped_fields


subroutine output_wrapped_fields(nout, ywflds)

  ! output the adress of all data fields in ywflds

  type(wrapped_fields), intent(in) :: ywflds
  integer(kind=jpim), intent(in) :: nout

  write(nout,*) "ywflds%spvor", loc(ywflds%spvor)
  write(nout,*) "ywflds%spdiv", loc(ywflds%spdiv)
  write(nout,*) "ywflds%spscalar", loc(ywflds%spscalar)
  write(nout,*) "ywflds%spscalar2", loc(ywflds%spscalar2)

  write(nout,*) "ywflds%u", loc(ywflds%u)
  write(nout,*) "ywflds%v", loc(ywflds%v)
  write(nout,*) "ywflds%u_ns", loc(ywflds%u_ns)
  write(nout,*) "ywflds%v_ns", loc(ywflds%v_ns)
  write(nout,*) "ywflds%scalar", loc(ywflds%scalar)
  write(nout,*) "ywflds%scalar_ew", loc(ywflds%scalar_ew)
  write(nout,*) "ywflds%scalar_ns", loc(ywflds%scalar_ns)
  write(nout,*) "ywflds%vor", loc(ywflds%vor)
  write(nout,*) "ywflds%div", loc(ywflds%div)

  write(nout,*) "ywflds%scalar2", loc(ywflds%scalar2)
  write(nout,*) "ywflds%scalar2_ew", loc(ywflds%scalar2_ew)
  write(nout,*) "ywflds%scalar2_ns", loc(ywflds%scalar2_ns)
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
