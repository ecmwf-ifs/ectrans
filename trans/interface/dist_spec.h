subroutine dist_spec(pspecg,kfdistg,kfrom,kvset,pspec)


#include "tsmbkind.h"

REAL_B    ,OPTIONAL, INTENT(IN) :: pspecg(:,:)
INTEGER_M          , INTENT(IN) :: kfdistg
INTEGER_M          , INTENT(IN) :: kfrom(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSET(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: pspec(:,:)
end subroutine dist_spec
