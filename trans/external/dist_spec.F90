subroutine dist_spec(pspecg,kfdistg,kfrom,kvset,pspec)


#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE DIST_SPEC_CONTROL_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN) :: pspecg(:,:)
INTEGER_M          , INTENT(IN) :: kfdistg
INTEGER_M          , INTENT(IN) :: kfrom(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSET(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: pspec(:,:)

!ifndef INTERFACE

INTEGER_M :: IVSET(kfdistg)
INTEGER_M :: IFSEND,IFRECV,J

if(ubound(kfrom,1) < kfdistg) then
 call abor1('dist_spec: kfrom too short!')
endif
 
ifsend = 0
do j=1,kfdistg
  if(kfrom(j) == myproc) ifsend = ifsend+1
enddo

if(ifsend > 0) then
  if(.not.present(pspecg)) then
    call abor1('dist_spec:pspecg missing')
  endif
  if(ubound(pspecg,1) < ifsend) then
    call abor1('dist_spec:first dimension of pspecg too small')
  endif 
 if(ubound(pspecg,2) < R%nspec2_g) then
    call abor1('dist_spec:first dimension of pspecg too small')
  endif
endif

IF(PRESENT(KVSET)) THEN
  if(ubound(kvset,1) < kfdistg) then
    call abor1('dist_spec: kvset too short!')
  endif
  DO J=1,kfdistg
    IF(KVSET(J) > NPRTRV) THEN
      WRITE(NERR,*) 'DIST_SPEC:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABOR1('DIST_SPEC:KVSET CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    if(KVSET(J) == mysetv) then
      ifrecv = ifrecv+1
    endif
  ENDDO
  ivset(:) = kvset(1:kfdistg)
ELSE
  ifrecv = kfdistg
  ivset(:) = mysetv
ENDIF

if(ifrecv > 0 ) then
  if(.not.present(pspec)) then
    call abor1('DIST_SPEC: fields to recieve and pspec not present')
  endif
  if(ubound(pspec,1) < ifrecv) then
    call abor1('DIST_SPEC: first dimension of pspec too small')
  endif
  if(ubound(pspec,2) < d%nspec2 ) then
    call abor1('DIST_SPEC: second dimension of pspec too small')
  endif
endif

call dist_spec_control(pspecg,kfdistg,kfrom,ivset,pspec)

!endif INTERFACE

end subroutine dist_spec
