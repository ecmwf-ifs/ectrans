MODULE TPM_FIELDS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL(KIND=JPRB) ,POINTER :: RPNM(:,:) ! Legendre polynomials
REAL(KIND=JPRB) ,POINTER :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL(KIND=JPRB) ,POINTER :: RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRB) ,POINTER :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL(KIND=JPRB) ,POINTER :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL(KIND=JPRB) ,POINTER :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL(KIND=JPRB) ,POINTER :: RN(:)     ! n (to avoid integer to real conversion)
REAL(KIND=JPRB) ,POINTER :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER(KIND=JPIM) ,POINTER :: NLTN(:) ! R%NTMAX+2-JN
END TYPE FIELDS_TYPE

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

END MODULE TPM_FIELDS
