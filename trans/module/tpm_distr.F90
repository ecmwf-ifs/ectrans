MODULE TPM_DISTR
#include "tsmbkind.h"

IMPLICIT NONE

SAVE

INTEGER_M :: NPROC     ! Number of processors (NPRGPNS*NPRGPEW) 
INTEGER_M :: NPRGPNS   ! No. of sets in N-S direction (grid-point space)
INTEGER_M :: NPRGPEW   ! No. of sets in E-W direction (grid-point space)
INTEGER_M :: NPRTRW    ! No. of sets in wave direction (spectral space)
INTEGER_M :: NPRTRV    ! NPROC/NPRTRW
INTEGER_M :: NPRTRNS   ! No. of sets in N-S direction (Fourier space)
INTEGER_M :: MYPROC    ! My processor number
INTEGER_M :: MYSETNS   ! My set number in N-S direction (grid-point space)
INTEGER_M :: MYSETEW   ! My set number in E-W direction (grid-point space)
INTEGER_M :: MYSETW    ! My set number in wave direction (spectral space) 
INTEGER_M :: MYSETV    ! My set number in field direction(S.S and F.S)
INTEGER_M :: NCOMBFLEN ! Size of communication buffer

INTEGER_M :: MTAGLETR
INTEGER_M :: MTAGML
INTEGER_M :: MTAGLG
INTEGER_M :: MTAGGL
INTEGER_M :: MTAGPART
INTEGER_M :: MTAGDISTSP
INTEGER_M :: MTAGLM
INTEGER_M :: MTAGDISTGP
 
INTEGER_M :: MINTET=1
INTEGER_M :: MREALT=2
INTEGER_M :: MLOGIT=3
INTEGER_M :: MCHART=4
INTEGER_M ,ALLOCATABLE :: NPRCIDS(:) ! Array containing the process ids

TYPE DISTR_TYPE
LOGICAL   :: LSPLIT
INTEGER_M :: NAPSETS
!Spectral space
INTEGER_M :: NUMP      ! No. of spectral waves handled by this processor
INTEGER_M :: NSPEC     ! No. of complex spectral coefficients (on this PE)
INTEGER_M :: NSPEC2    ! 2*NSPEC
INTEGER_M :: NSPEC2MX  ! maximun NSPEC2 among all PEs
INTEGER_M :: NTPEC2
INTEGER_M :: NUMTP



INTEGER_M :: NSPOLEGL  ! No. of legendre polynomials on this PE
INTEGER_M :: NLEI3D    ! (NLEI3-1)/NPRTRW+1

INTEGER_M ,POINTER :: MYMS(:)    ! Wave numbers handled by this PE
INTEGER_M ,POINTER :: NUMPP(:)   ! No. of wave numbers each wave set is
                                     ! responsible for
INTEGER_M ,POINTER :: NPOSSP(:) !Not needed in transform?
INTEGER_M ,POINTER :: NPROCM(:) ! Process that does the calc. for certain 
!                                      wavenumber M
INTEGER_M ,POINTER :: NDIM0G(:)

INTEGER_M ,POINTER :: NASM0(:)  ! Address in a spectral array of (m, n=m)
INTEGER_M ,POINTER :: NATM0(:)


INTEGER_M ,POINTER :: NLATLS(:)
INTEGER_M ,POINTER :: NLATLE(:)

INTEGER_M ,POINTER :: NPMT(:)
INTEGER_M ,POINTER :: NPMS(:)
INTEGER_M ,POINTER :: NPMG(:)

INTEGER_M :: NDGL_FS

INTEGER_M ,POINTER :: NULTPP(:) ! No of lats. for each wave_set  (F.S)
INTEGER_M ,POINTER :: NPROCL(:) ! Process responsible for each lat. (F.S)
INTEGER_M ,POINTER :: NPTRLS(:) ! Pointer to first lat. (F.S)
INTEGER_M ,POINTER :: NALLMS(:)
INTEGER_M ,POINTER :: NPTRMS(:)

INTEGER_M ,POINTER :: NSTAGT0B(:)
INTEGER_M ,POINTER :: NSTAGT1B(:)
INTEGER_M ,POINTER :: NPNTGTB0(:,:)
INTEGER_M ,POINTER :: NPNTGTB1(:,:)
INTEGER_M ,POINTER :: NLTSFTB(:)

INTEGER_M ,POINTER :: NLTSGTB(:)
INTEGER_M ,POINTER :: MSTABF(:)

INTEGER_M :: NLENGT0B
INTEGER_M :: NLENGT1B

! GRIDPOINT SPACE

INTEGER_M :: NDGL_GP
INTEGER_M ,POINTER :: NFRSTLAT(:)
INTEGER_M ,POINTER :: NLSTLAT(:)
INTEGER_M :: NFRSTLOFF
INTEGER_M ,POINTER :: NPTRLAT(:)
INTEGER_M ,POINTER :: NPTRFRSTLAT(:)
INTEGER_M ,POINTER :: NPTRLSTLAT(:)
INTEGER_M :: NPTRFLOFF
LOGICAl   ,POINTER :: LSPLITLAT(:)

INTEGER_M ,POINTER :: NSTA(:,:)
INTEGER_M ,POINTER :: NONL(:,:)

INTEGER_M :: NGPTOT
INTEGER_M :: NGPTOTG
INTEGER_M :: NGPTOTMX

INTEGER_M ,POINTER  :: NSTAGTF(:)
INTEGER_M :: NLENGTF
END TYPE DISTR_TYPE

TYPE(DISTR_TYPE),ALLOCATABLE,TARGET :: DISTR_RESOL(:)
TYPE(DISTR_TYPE),POINTER     :: D

end module tpm_distr







