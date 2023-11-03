module ectrans_mod
#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "specnorm.h"
end module ectrans_mod
submodule (ectrans_mod) setup_trans0_mod
#include "setup_trans0.F90"
end submodule setup_trans0_mod
submodule (ectrans_mod) setup_trans_mod
#include "setup_trans.F90"
end submodule setup_trans_mod
submodule (ectrans_mod) inv_trans_mod
#include "inv_trans.F90"
end submodule inv_trans_mod
submodule (ectrans_mod) inv_trans_mod
#include "inv_trans.F90"
end submodule inv_trans_mod
submodule (ectrans_mod) specnorm_mod
#include "specnorm.F90"
end submodule specnorm_mod

