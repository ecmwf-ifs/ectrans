MODULE TPM_FIELDS

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

TYPE FIELDS_TYPE
REAL_B ,POINTER :: RPNM(:,:) ! Legendre polynomials
REAL_B ,POINTER :: RMU(:)    ! sin(theta) for Gaussian latitudes
REAL_B ,POINTER :: RW(:)     ! Weights of the Gaussian quadrature
REAL_B ,POINTER :: R1MU2(:)  ! 1.-MU*MU, cos(theta)**2
REAL_B ,POINTER :: RACTHE(:) ! 1./SQRT(R1MU2), 1/(cos(theta))

REAL_B ,POINTER :: REPSNM(:) ! eps(n,m) used in the Legendre transforms
REAL_B ,POINTER :: RN(:)     ! n (to avoid integer to real conversion)
REAL_B ,POINTER :: RLAPIN(:) ! eigen-values of the inverse Laplace operator
INTEGER_M ,POINTER :: NLTN(:) ! R%NTMAX+2-JN
END TYPE FIELDS_TYPE

TYPE(FIELDS_TYPE),ALLOCATABLE,TARGET :: FIELDS_RESOL(:)
TYPE(FIELDS_TYPE),POINTER     :: F

END MODULE TPM_FIELDS
