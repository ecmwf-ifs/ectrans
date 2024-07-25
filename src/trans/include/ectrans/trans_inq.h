! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE TRANS_INQ(KRESOL,KSPEC,KSPEC2,KSPEC2G,KSPEC2MX,KNUMP,&
                    &KGPTOT,KGPTOTG,KGPTOTMX,KGPTOTL,&
                    &KMYMS,KASM0,KUMPP,KPOSSP,KPTRMS,KALLMS,KDIM0G,&
                    &KFRSTLAT,KLSTLAT,KFRSTLOFF,KPTRLAT,&
                    &KPTRFRSTLAT,KPTRLSTLAT,KPTRFLOFF,KSTA,KONL,&
                    &KULTPP,KPTRLS,KNMENG,&
                    &KPRTRW,KMYSETW,KMYSETV,KMY_REGION_NS,KMY_REGION_EW,&
                    &LDSPLITLAT,&
                    &KSMAX,PLAPIN,KNVALUE,KDEF_RESOL,LDLAM,&
                    &PMU,PGW,PRPNM,KLEI3,KSPOLEGL,KPMS,KDGLU)

! begin_doc_block
! ## `TRANS_INQ`
!
! ### Signature
!
! ```f90
! SUBROUTINE TRANS_INQ(KRESOL, KSPEC, KSPEC2, KSPEC2G, KSPEC2MX, KNUMP, KGPTOT, KGPTOTG, KGPTOTMX, &
!   &                  KGPTOTL, KMYMS, KASM0, KUMPP, KPOSSP, KPTRMS, KALLMS, KDIM0G, KFRSTLAT, &
!   &                  KLSTLAT, KFRSTLOFF, KPTRLAT, KPTRFRSTLAT, KPTRLSTLAT, KPTRFLOFF, KSTA, &
!   &                  KONL, KULTPP, KPTRLS, KNMENG, KPRTRW, KMYSETW, KMYSETV, KMY_REGION_NS, &
!   &                  KMY_REGION_EW, LDSPLITLAT, KSMAX, PLAPIN, KNVALUE, KDEF_RESOL, LDLAM, PMU, &
!   &                  PGW, PRPNM, KLEI3, KSPOLEGL, KPMS, KDGLU)
! ```
!
! ### Purpose
!
! This subroutine allows one to obtain information on the transform configuration for a given
! `KRESOL` handle (or the default handle if this argument is omitted).
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle for which information is required.
!   *Default*: The first defined resolution handle.
!
! ### `OPTIONAL, INTENT(OUT)` arguments
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC`  
!   Number of complex spectral coefficients on this MPI task.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC2`  
!   Number of real and imaginary spectral coefficients on this MPI task.
!   Equal to `2 * KSPEC`.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC2G`  
!   Number of real and imaginary spectral coefficients across all MPI tasks.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC2MX`  
!   The maximum value of KSPEC2 over all MPI tasks.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KNUMP`  
!   Number of zonal wave numbers handled by this MPI task for the Legendre transform.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KGPTOT`  
!   Number of grid columns on this MPI task.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KGPTOTG`  
!   Total number of grid columns across all MPI tasks.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KGPTOTMX`  
!   Maximum number of grid columns on any of the MPI tasks.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KGPTOTL(:,:)`  
!   TODO: check. Number of grid columns one each PE (dimension N_REGIONS_NS:N_REGIONS_EW)
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KMYMS(:)`  
!   The zonal wavenumbers handled by this MPI task for the Legendre transform.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KASM0(0:)`  
!   Address in a spectral array of (m, n=m).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KUMPP(:)`  
!   TODO: check this. Number of wave numbers each wave set is responsible for
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPOSSP(:)`  
!  Defines partitioning of global spectral fields among PEs
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRMS(:)`  
!   Pointer to the first wave number of a given A-set
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KALLMS(:)`  
!   Wave numbers for all wave-set concatenated together to give all wave numbers in wave-set order
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KDIM0G(0:)`  
!   Defines partitioning of global spectral fields among PEs
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KFRSTLAT(:)`  
!   First latitude of each a-set in grid-point space
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KLSTLAT(:)`  
!   Last latitude of each a-set in grid-point space
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KFRSTLOFF`  
!   Offset for first lat of own a-set in grid-point space
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRLAT(:)`  
!   Pointer to the start of each latitude
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRFRSTLAT(:)`  
!   Pointer to the first latitude of each a-set in NSTA and NONL arrays
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRLSTLAT(:)`  
!   Pointer to the last latitude of each a-set in NSTA and NONL arrays
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRFLOFF`  
!   Offset for pointer to the first latitude of own a-set NSTA and NONL arrays, i.e. nptrfrstlat  
!   (myseta)-1
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSTA(:,:)`  
!   Position of first grid column for the latitudes on a processor. The information is available  
!   for all processors. The b-sets are distinguished by the last dimension of nsta(). The latitude  
!   band for each a-set is addressed by nptrfrstlat(jaset), nptrlstlat(jaset), and  
!   nptrfloff=nptrfrstlat(myseta) on this processors a-set. Each split latitude has two entries in  
!   nsta(,:) which necessitates the rather complex addressing of nsta(,:) and the overdimensioning  
!   of nsta by N_REGIONS_NS.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KONL(:,:)`  
!   Number of grid columns for the latitudes on a processor. Similar to nsta() in data structure.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KULTPP(:)`  
!   number of latitudes for which each a-set is calculating the FFTs.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRLS(:)`  
!   pointer to first global latitude of each a-set for which it performs the Fourier calculations
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KNMENG(:)`  
!   associated (with NLOENG) cut-off zonal wavenumber
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPRTRW`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KMYSETW`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KMYSETV`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KMY_REGION_NS`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KMY_REGION_EW`
! - `LOGICAL, OPTIONAL, INTENT(OUT) :: LDSPLITLAT(:)`  
!   TRUE if latitude is split in grid point space over two a-sets
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSMAX`  
!   spectral truncation
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PLAPIN(-1:)`  
!   Eigenvalues of the inverse Laplace operator
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KNVALUE(:)`  
!   n value for each KSPEC2 spectral coeffient
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KDEF_RESOL`  
!   number or resolutions defined
! - `LOGICAL, OPTIONAL, INTENT(OUT) :: LDLAM`  
!   .T. if the corresponding resolution is LAM, .F. if it is global
! - `REAL(KIND=JPRD), OPTIONAL, INTENT(OUT) :: PMU(:)`  
!   sin(Gaussian latitudes)
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGW(:)`  
!   Gaussian weights
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PRPNM(:,:)`  
!   Legendre polynomials
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KLEI3`  
!   First dimension of Legendre polynomials
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPOLEGL`  
!   Second dimension of Legendre polynomials
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPMS(0:)`  
!   Address for legendre polynomial for given M (NSMAX)
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KDGLU(0:)`  
!   Number of active points in an hemisphere for a given wavenumber "m"
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD


IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC2
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC2G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPEC2MX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNUMP
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOT
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOTG
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOTMX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KGPTOTL(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KFRSTLOFF
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRFLOFF

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYMS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KASM0(0:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KUMPP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPOSSP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRMS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KALLMS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDIM0G(0:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KFRSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KLSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRFRSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRLSTLAT(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSTA(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KONL(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPRTRW
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYSETW
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMYSETV
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMY_REGION_NS
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KMY_REGION_EW
LOGICAL   ,OPTIONAL, INTENT(OUT) :: LDSPLITLAT(:)

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KULTPP(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPTRLS(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNMENG(:)

REAL(KIND=JPRD)    ,OPTIONAL, INTENT(OUT) :: PMU(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PGW(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PRPNM(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KLEI3
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSPOLEGL
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KPMS(0:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDGLU(0:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PLAPIN(-1:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KNVALUE(:)

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(OUT) :: KDEF_RESOL
LOGICAL           ,OPTIONAL,INTENT(OUT)   :: LDLAM

END SUBROUTINE TRANS_INQ






END INTERFACE
