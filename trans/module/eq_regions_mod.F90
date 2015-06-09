MODULE EQ_REGIONS_MOD
!
!     Purpose.
!     --------
!           eq_regions_mod provides the code to perform a high level
!           partitioning of the surface of a sphere into regions of
!           equal area and small diameter.
!           the type.
!
!     Background.
!     -----------
!     This Fortran version of eq_regions is a much cut down version of the
!     "Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox" of the
!     same name developed by Paul Leopardi at the University of New South Wales.
!     This version has been coded specifically for the case of partitioning the
!     surface of a sphere or S^dim (where dim=2) as denoted in the original code.
!     Only a subset of the original eq_regions package has been coded to determine
!     the high level distribution of regions on a sphere, as the detailed
!     distribution of grid points to each region is left to IFS software.
!     This is required to take into account the spatial distribution of grid
!     points in an IFS gaussian grid and provide an optimal (i.e. exact)
!     distribution of grid points over regions.
!
!     The following copyright notice for the eq_regions package is included from
!     the original MatLab release.
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     + Release 1.10 2005-06-26                                                 +
!     +                                                                         +
!     + Copyright (c) 2004, 2005, University of New South Wales                 +
!     +                                                                         +
!     + Permission is hereby granted, free of charge, to any person obtaining   +
!     + a copy of this software and associated documentation files (the         +
!     + "Software"), to deal in the Software without restriction, including     +
!     + without limitation the rights to use, copy, modify, merge, publish,     +
!     + distribute, sublicense, and/or sell copies of the Software, and to      +
!     + permit persons to whom the Software is furnished to do so, subject to   +
!     + the following conditions:                                               +
!     +                                                                         +
!     + The above copyright notice and this permission notice shall be included +
!     + in all copies or substantial portions of the Software.                  +
!     +                                                                         +
!     + THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         +
!     + EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      +
!     + MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  +
!     + IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    +
!     + CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    +
!     + TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       +
!     + SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  +
!     +                                                                         +
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Author.
!     -------
!        George Mozdzynski *ECMWF*
!
!     Modifications.
!     --------------
!        Original : 2006-04-15
!
!--------------------------------------------------------------------------------
!
USE PARKIND1  ,ONLY : JPIM,   JPRB

IMPLICIT NONE

SAVE

PRIVATE

PUBLIC EQ_REGIONS,L_REGIONS_DEBUG,N_REGIONS_NS,N_REGIONS_EW,N_REGIONS,MY_REGION_NS,MY_REGION_EW

REAL(KIND=JPRB) PI
LOGICAL :: L_REGIONS_DEBUG=.FALSE.
INTEGER(KIND=JPIM) :: N_REGIONS_NS
INTEGER(KIND=JPIM) :: N_REGIONS_EW
INTEGER(KIND=JPIM) :: MY_REGION_NS
INTEGER(KIND=JPIM) :: MY_REGION_EW
INTEGER(KIND=JPIM),ALLOCATABLE :: N_REGIONS(:)

CONTAINS

SUBROUTINE EQ_REGIONS(N)
!
! eq_regions uses the zonal equal area sphere partitioning algorithm to partition
! the surface of a sphere into N regions of equal area and small diameter.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N
INTEGER(KIND=JPIM) :: N_COLLARS,J
REAL(KIND=JPRB),ALLOCATABLE :: R_REGIONS(:)
REAL(KIND=JPRB) :: C_POLAR

PI=2.0_JPRB*ASIN(1.0_JPRB)

N_REGIONS(:)=0

IF( N == 1 )THEN

  !
  ! We have only one region, which must be the whole sphere.
  !
  N_REGIONS(1)=1
  N_REGIONS_NS=1

ELSE

  !
  ! Given N, determine c_polar
  ! the colatitude of the North polar spherical cap.
  !
  C_POLAR = POLAR_COLAT(N)
  !
  ! Given N, determine the ideal angle for spherical collars.
  ! Based on N, this ideal angle, and c_polar,
  ! determine n_collars, the number of collars between the polar caps.
  !
  N_COLLARS = NUM_COLLARS(N,C_POLAR,IDEAL_COLLAR_ANGLE(N))
  N_REGIONS_NS=N_COLLARS+2
  !
  ! Given N, c_polar and n_collars, determine r_regions,
  ! a list of the ideal real number of regions in each collar,
  ! plus the polar caps.
  ! The number of elements is n_collars+2.
  ! r_regions[1] is 1.
  ! r_regions[n_collars+2] is 1.
  ! The sum of r_regions is N.
  ALLOCATE(R_REGIONS(N_COLLARS+2))
  CALL IDEAL_REGION_LIST(N,C_POLAR,N_COLLARS,R_REGIONS)
  !
  ! Given N and r_regions, determine n_regions, a list of the natural number
  ! of regions in each collar and the polar caps.
  ! This list is as close as possible to r_regions.
  ! The number of elements is n_collars+2.
  ! n_regions[1] is 1.
  ! n_regions[n_collars+2] is 1.
  ! The sum of n_regions is N.
  !
  CALL ROUND_TO_NATURALS(N,N_COLLARS,R_REGIONS)
  DEALLOCATE(R_REGIONS)
  IF( N /= SUM(N_REGIONS(:)) )THEN
    WRITE(*,'("eq_regions: N=",I10," sum(n_regions(:))=",I10)')N,SUM(N_REGIONS(:))
    CALL ABOR1('eq_regions: N /= sum(n_regions)')
  ENDIF

ENDIF

IF( L_REGIONS_DEBUG )THEN
  WRITE(*,'("eq_regions: N=",I6," n_regions_ns=",I4)') N,N_REGIONS_NS
  DO J=1,N_REGIONS_NS
    WRITE(*,'("eq_regions: n_regions(",I4,")=",I4)') J,N_REGIONS(J)
  ENDDO
ENDIF
N_REGIONS_EW=MAXVAL(N_REGIONS(:))

RETURN
END SUBROUTINE EQ_REGIONS

FUNCTION NUM_COLLARS(N,C_POLAR,A_IDEAL) RESULT(NUM_C)
!
!NUM_COLLARS The number of collars between the polar caps
!
! Given N, an ideal angle, and c_polar,
! determine n_collars, the number of collars between the polar caps.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N
REAL(KIND=JPRB),INTENT(IN) :: A_IDEAL,C_POLAR
INTEGER(KIND=JPIM) :: NUM_C
LOGICAL ENOUGH
ENOUGH = (N > 2) .AND. (A_IDEAL > 0)
IF( ENOUGH )THEN
  NUM_C = MAX(1,NINT((PI-2.*C_POLAR)/A_IDEAL))
ELSE
  NUM_C = 0
ENDIF
RETURN
END FUNCTION NUM_COLLARS

SUBROUTINE IDEAL_REGION_LIST(N,C_POLAR,N_COLLARS,R_REGIONS)
!
!IDEAL_REGION_LIST The ideal real number of regions in each zone
!
! List the ideal real number of regions in each collar, plus the polar caps.
!
! Given N, c_polar and n_collars, determine r_regions, a list of the ideal real
! number of regions in each collar, plus the polar caps.
! The number of elements is n_collars+2.
! r_regions[1] is 1.
! r_regions[n_collars+2] is 1.
! The sum of r_regions is N.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N,N_COLLARS
REAL(KIND=JPRB),INTENT(IN) :: C_POLAR
REAL(KIND=JPRB),INTENT(OUT) :: R_REGIONS(N_COLLARS+2)
INTEGER(KIND=JPIM) :: COLLAR_N
REAL(KIND=JPRB) :: IDEAL_REGION_AREA,IDEAL_COLLAR_AREA
REAL(KIND=JPRB) :: A_FITTING
R_REGIONS(:)=0.0_JPRB
R_REGIONS(1) = 1.0_JPRB
IF( N_COLLARS > 0 )THEN
  !
  ! Based on n_collars and c_polar, determine a_fitting,
  ! the collar angle such that n_collars collars fit between the polar caps.
  !
  A_FITTING = (PI-2.0_JPRB*C_POLAR)/FLOAT(N_COLLARS)
  IDEAL_REGION_AREA = AREA_OF_IDEAL_REGION(N)
  DO COLLAR_N=1,N_COLLARS
    IDEAL_COLLAR_AREA = AREA_OF_COLLAR(C_POLAR+(COLLAR_N-1)*A_FITTING, &
     & C_POLAR+COLLAR_N*A_FITTING)
    R_REGIONS(1+COLLAR_N) = IDEAL_COLLAR_AREA / IDEAL_REGION_AREA
  ENDDO
ENDIF
R_REGIONS(2+N_COLLARS) = 1.
RETURN
END SUBROUTINE IDEAL_REGION_LIST

FUNCTION IDEAL_COLLAR_ANGLE(N) RESULT(IDEAL)
!
! IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
!
! IDEAL_COLLAR_ANGLE(N) sets ANGLE to the ideal angle for the
! spherical collars of an EQ partition of the unit sphere S^2 into N regions.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N
REAL(KIND=JPRB) :: IDEAL
IDEAL = AREA_OF_IDEAL_REGION(N)**(0.5_JPRB)
RETURN
END FUNCTION IDEAL_COLLAR_ANGLE

SUBROUTINE ROUND_TO_NATURALS(N,N_COLLARS,R_REGIONS)
!
! ROUND_TO_NATURALS Round off a given list of numbers of regions
!
! Given N and r_regions, determine n_regions, a list of the natural number
! of regions in each collar and the polar caps.
! This list is as close as possible to r_regions, using rounding.
! The number of elements is n_collars+2.
! n_regions[1] is 1.
! n_regions[n_collars+2] is 1.
! The sum of n_regions is N.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N,N_COLLARS
REAL(KIND=JPRB),INTENT(IN) :: R_REGIONS(N_COLLARS+2)
INTEGER(KIND=JPIM) :: ZONE_N
REAL(KIND=JPRB) :: DISCREPANCY
N_REGIONS(1:N_COLLARS+2) = R_REGIONS(:)
DISCREPANCY = 0.0_JPRB
DO ZONE_N = 1,N_COLLARS+2
    N_REGIONS(ZONE_N) = NINT(R_REGIONS(ZONE_N)+DISCREPANCY);
    DISCREPANCY = DISCREPANCY+R_REGIONS(ZONE_N)-FLOAT(N_REGIONS(ZONE_N));
ENDDO
RETURN
END SUBROUTINE ROUND_TO_NATURALS

FUNCTION POLAR_COLAT(N) RESULT(POLAR_C)
!
! Given N, determine the colatitude of the North polar spherical cap.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N
REAL(KIND=JPRB) :: AREA
REAL(KIND=JPRB) :: POLAR_C
IF( N == 1 ) POLAR_C=PI
IF( N == 2 ) POLAR_C=PI/2.0_JPRB
IF( N > 2 )THEN
  AREA=AREA_OF_IDEAL_REGION(N)
  POLAR_C=SRADIUS_OF_CAP(AREA)
ENDIF
RETURN
END FUNCTION POLAR_COLAT

FUNCTION AREA_OF_IDEAL_REGION(N) RESULT(AREA)
!
! AREA_OF_IDEAL_REGION(N) sets AREA to be the area of one of N equal
! area regions on S^2, that is 1/N times AREA_OF_SPHERE.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: N
REAL(KIND=JPRB) :: AREA_OF_SPHERE
REAL(KIND=JPRB) :: AREA
AREA_OF_SPHERE = (2.0_JPRB*PI**1.5_JPRB/GAMMA(1.5_JPRB))
AREA = AREA_OF_SPHERE/FLOAT(N)
RETURN
END FUNCTION AREA_OF_IDEAL_REGION

FUNCTION SRADIUS_OF_CAP(AREA) RESULT(SRADIUS)
!
! SRADIUS_OF_CAP(AREA) returns the spherical radius of
! an S^2 spherical cap of area AREA.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN) :: AREA
REAL(KIND=JPRB) :: SRADIUS
SRADIUS = 2.0_JPRB*ASIN(SQRT(AREA/PI)/2.0_JPRB)
RETURN
END FUNCTION SRADIUS_OF_CAP

FUNCTION AREA_OF_COLLAR(A_TOP, A_BOT) RESULT(AREA)
!
! AREA_OF_COLLAR Area of spherical collar
!
! AREA_OF_COLLAR(A_TOP, A_BOT) sets AREA to be the area of an S^2 spherical
! collar specified by A_TOP, A_BOT, where A_TOP is top (smaller) spherical radius,
! A_BOT is bottom (larger) spherical radius.
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN) :: A_TOP,A_BOT
REAL(KIND=JPRB) AREA
AREA = AREA_OF_CAP(A_BOT) - AREA_OF_CAP(A_TOP)
RETURN
END FUNCTION AREA_OF_COLLAR

FUNCTION AREA_OF_CAP(S_CAP) RESULT(AREA)
!
! AREA_OF_CAP Area of spherical cap
!
! AREA_OF_CAP(S_CAP) sets AREA to be the area of an S^2 spherical
! cap of spherical radius S_CAP.
!
REAL(KIND=JPRB),INTENT(IN) :: S_CAP
REAL(KIND=JPRB) AREA
AREA = 4.0_JPRB*PI * SIN(S_CAP/2.0_JPRB)**2
RETURN
END FUNCTION AREA_OF_CAP

FUNCTION GAMMA(X) RESULT(GAMMA_RES)
!
USE PARKIND1  ,ONLY : JPIM,   JPRB
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN) :: X
REAL(KIND=JPRB) :: GAMMA_RES
REAL(KIND=JPRB) :: P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13
REAL(KIND=JPRB) :: W,Y
INTEGER(KIND=JPIM) :: K,N
PARAMETER (&
& P0 =   0.999999999999999990E+00_JPRB,&
& P1 =  -0.422784335098466784E+00_JPRB,&
& P2 =  -0.233093736421782878E+00_JPRB,&
& P3 =   0.191091101387638410E+00_JPRB,&
& P4 =  -0.024552490005641278E+00_JPRB,&
& P5 =  -0.017645244547851414E+00_JPRB,&
& P6 =   0.008023273027855346E+00_JPRB)
PARAMETER (&
& P7 =  -0.000804329819255744E+00_JPRB,&
& P8 =  -0.000360837876648255E+00_JPRB,&
& P9 =   0.000145596568617526E+00_JPRB,&
& P10 = -0.000017545539395205E+00_JPRB,&
& P11 = -0.000002591225267689E+00_JPRB,&
& P12 =  0.000001337767384067E+00_JPRB,&
& P13 = -0.000000199542863674E+00_JPRB)
N = NINT(X - 2)
W = X - (N + 2)
Y = ((((((((((((P13 * W + P12) * W + P11) * W + P10) *&
&    W + P9) * W + P8) * W + P7) * W + P6) * W + P5) *&
&    W + P4) * W + P3) * W + P2) * W + P1) * W + P0
IF (N .GT. 0) THEN
  W = X - 1
  DO K = 2, N
    W = W * (X - K)
  END DO
ELSE
  W = 1
  DO K = 0, -N - 1
    Y = Y * (X + K)
  END DO
END IF
GAMMA_RES = W / Y
RETURN
END FUNCTION GAMMA

END MODULE EQ_REGIONS_MOD
