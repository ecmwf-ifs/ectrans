MODULE TPM_FIELDS

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRB) ,ALLOCATABLE :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRB) ,ALLOCATABLE :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRB) ,ALLOCATABLE :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRB) ,ALLOCATABLE :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRB) ,ALLOCATABLE :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRB) ,ALLOCATABLE :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRB) ,ALLOCATABLE :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRB) ,ALLOCATABLE :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLTN(:) ! R%NTMAX+2-JN

REAL(KIND=JPRB) ,ALLOCATABLE :: RMU2(:)    ! sin(theta) for dual input/output latitudes
REAL(KIND=JPRB) ,ALLOCATABLE :: RACTHE2(:) ! 1./SQRT(R1MU2), 1/(cos(theta)) dual input/output latitudes
END TYPE FIELDS_TYPE

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

END MODULE TPM_FIELDS
