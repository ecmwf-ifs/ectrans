MODULE TPM_DISTR

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!*    Variables describing distributed memory parallelization

INTEGER_M :: NPROC     ! Number of processors (NPRGPNS*NPRGPEW) 
INTEGER_M :: NPRGPNS   ! No. of sets in N-S direction (grid-point space)
INTEGER_M :: NPRGPEW   ! No. of sets in E-W direction (grid-point space)
INTEGER_M :: NPRTRW    ! No. of sets in wave direction (spectral space)
INTEGER_M :: NPRTRV    ! NPROC/NPRTRW
INTEGER_M :: NPRTRNS   ! No. of sets in N-S direction (Fourier space)
                       ! (always equal to NPRTRW)
INTEGER_M :: MYPROC    ! My processor number
INTEGER_M :: MYSETNS   ! My set number in N-S direction (grid-point space)
INTEGER_M :: MYSETEW   ! My set number in E-W direction (grid-point space)
INTEGER_M :: MYSETW    ! My set number in wave direction (spectral space) 
INTEGER_M :: MYSETV    ! My set number in field direction(S.S and F.S)
INTEGER_M :: NCOMBFLEN ! Size of communication buffer

INTEGER_M :: MTAGLETR   ! Tag
INTEGER_M :: MTAGML     ! Tag
INTEGER_M :: MTAGLG     ! Tag
INTEGER_M :: MTAGGL     ! Tag
INTEGER_M :: MTAGPART   ! Tag
INTEGER_M :: MTAGDISTSP ! Tag
INTEGER_M :: MTAGLM     ! Tag
INTEGER_M :: MTAGDISTGP ! Tag
 
INTEGER_M ,ALLOCATABLE :: NPRCIDS(:) ! Array containing the process ids

TYPE DISTR_TYPE
LOGICAL   :: LSPLIT    ! TRUE - latitudes are shared between a-sets
INTEGER_M :: NAPSETS   ! Number of apple sets at the poles. Default is zero.

! SPECTRAL SPACE

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
INTEGER_M ,POINTER :: NPOSSP(:)  ! Not needed in transform?
INTEGER_M ,POINTER :: NPROCM(:)  ! Process that does the calc. for certain 
                                 ! wavenumber M
INTEGER_M ,POINTER :: NDIM0G(:)  ! Defines partitioning of global spectral
                                 ! fields among PEs

INTEGER_M ,POINTER :: NASM0(:)  ! Address in a spectral array of (m, n=m)
INTEGER_M ,POINTER :: NATM0(:)  ! Same as NASM0 but for NTMAX
INTEGER_M ,POINTER :: NALLMS(:) ! Wave numbers for all a-set concatenated
                                ! together to give all wave numbers in a-set
                                ! order. Used when global spectral norms 
                                ! have to be gathered.
INTEGER_M ,POINTER :: NPTRMS(:) ! Pointer to the first wave number of a given
                                ! a-set in nallms array.


! Legendre polynomials

INTEGER_M ,POINTER :: NLATLS(:) ! First latitude for which each a-set calcul.
INTEGER_M ,POINTER :: NLATLE(:) ! Last latitude for which each a-set calcul.

INTEGER_M ,POINTER :: NPMT(:) ! Adress for legendre polynomial for
                              ! given M (NTMAX)
INTEGER_M ,POINTER :: NPMS(:) ! Adress for legendre polynomial for
                              ! given M (NSMAX)
INTEGER_M ,POINTER :: NPMG(:) ! Global version of NPMS

! FOURIER SPACE

INTEGER_M :: NDGL_FS ! Number of rows of latitudes for which this process is
                     ! performing Fourier Space calculations

INTEGER_M ,POINTER  :: NSTAGTF(:) ! Offset for specific latitude in 
                                  ! Fourier/gridpoint buffer
INTEGER_M :: NLENGTF ! Second dimension of Fourier/gridpoint buffer 
                     ! (sum of (NLOEN+3) over local latitudes)

INTEGER_M ,POINTER :: NULTPP(:) ! No of lats. for each wave_set  (F.S)
INTEGER_M ,POINTER :: NPROCL(:) ! Process responsible for each lat. (F.S)
INTEGER_M ,POINTER :: NPTRLS(:) ! Pointer to first lat. (F.S)

INTEGER_M ,POINTER :: NSTAGT0B(:) ! Start adresses for segments within buffer
                                  ! (according to processors to whom data 
                                  ! is going to be sent) 
INTEGER_M ,POINTER :: NSTAGT1B(:) 
INTEGER_M ,POINTER :: NPNTGTB0(:,:)
INTEGER_M ,POINTER :: NPNTGTB1(:,:)
INTEGER_M ,POINTER :: NLTSFTB(:)  

INTEGER_M ,POINTER :: NLTSGTB(:)
INTEGER_M ,POINTER :: MSTABF(:)

INTEGER_M :: NLENGT0B
INTEGER_M :: NLENGT1B

! GRIDPOINT SPACE

INTEGER_M :: NDGL_GP ! D%NLSTLAT(MYSETNS)-D%NFRSTLOFF
INTEGER_M ,POINTER :: NFRSTLAT(:) ! First lat of each a-set 
INTEGER_M ,POINTER :: NLSTLAT(:)  ! Last lat of each a-set 
INTEGER_M :: NFRSTLOFF ! Offset for first lat of own a-set 
                       ! i.e. NFRSTLOFF=NFRSTLAT(MYSETA)-1
INTEGER_M ,POINTER :: NPTRLAT(:) ! Pointer to start of latitude 
INTEGER_M ,POINTER :: NPTRFRSTLAT(:) ! Pointer to the first latitude of each 
                                     ! a-set in NSTA and NONL arrays
INTEGER_M ,POINTER :: NPTRLSTLAT(:) ! Pointer to the last latitude of each
                                    ! a-set in NSTA and NONL arrays
INTEGER_M :: NPTRFLOFF ! Offset for pointer to the first latitude of own a-set
                       ! NSTA and NONL arrays, i.e. NPTRFRSTLAT(MYSETA)-1
LOGICAL   ,POINTER :: LSPLITLAT(:) ! True if latitude is split over 2 a-sets

!  NSTA(R%NDGL+NPRGPNS-1,NPRGPEW) :  Position of first grid column
!             for the latitudes on a processor. The information is
!             available for all processors. The b-sets are distinguished
!             by the last dimension of NSTA(). The latitude band for
!             each a-set is addressed by NPTRFRSTLAT(JASET),
!             NPTRLSTLAT(JASET), and NPTRFLOFF=NPTRFRSTLAT(MYSETA) on
!             this processors a-set. Each split latitude has two entries
!             in NSTA(,:) which necessitates the rather complex
!             addressing of NSTA(,:) and the overdimensioning of NSTA by
!             NPRGPNS.
!  NONL(R%NDGL+NPRGPNS-1,NPRGPEW)  :  Number of grid columns for
!             the latitudes on a processor. Similar to NSTA() in data
!             structure.
INTEGER_M ,POINTER :: NSTA(:,:)
INTEGER_M ,POINTER :: NONL(:,:)

INTEGER_M :: NGPTOT   ! Total number of grid columns on this PE
INTEGER_M :: NGPTOTG  ! Total number of grid columns on the Globe
INTEGER_M :: NGPTOTMX ! Maximum number of grid columns on any of the PEs
INTEGER_M ,POINTER :: NGPTOTL(:,:) ! Number of grid columns on each PE.

END TYPE DISTR_TYPE

TYPE(DISTR_TYPE),ALLOCATABLE,TARGET :: DISTR_RESOL(:)
TYPE(DISTR_TYPE),POINTER     :: D

END MODULE TPM_DISTR







