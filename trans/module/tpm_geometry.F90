MODULE TPM_GEOMETRY

! Module containing data describing Gaussian grid.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE GEOM_TYPE
INTEGER(KIND=JPIM),ALLOCATABLE :: NLOEN(:) ! NUMBER OF POINTS ON A PARALLEL
INTEGER(KIND=JPIM),ALLOCATABLE :: NMEN(:)  ! ASSOCIATED CUT-OFF WAVE NUMBER
INTEGER(KIND=JPIM),ALLOCATABLE :: NDGLU(:) ! NUMBER OF HEMISPERIC LATITUDES
!                                   FOR A GIVEN WAVE NUMBER M 

LOGICAL :: LAM           ! LAM geometry if T, Global geometry if F
LOGICAL :: LREDUCED_GRID ! Reduced Gaussian grid if T
!                          quadratic Gaussian grid otherwise.
REAL(KIND=JPRB) :: RSTRET ! Stretching factor (for Legendre polynomials
!                           computed on stretched latitudes only)
END TYPE GEOM_TYPE

TYPE(GEOM_TYPE),ALLOCATABLE,TARGET :: GEOM_RESOL(:)
TYPE(GEOM_TYPE),POINTER     :: G

END MODULE TPM_GEOMETRY
