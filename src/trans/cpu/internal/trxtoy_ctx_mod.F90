! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRXTOY_CTX_MOD

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

TYPE TRXTOY_CTX
  INTEGER(KIND=JPIM) :: ISENDCOUNT = -9999
  INTEGER(KIND=JPIM) :: IRECVCOUNT = -9999
  INTEGER(KIND=JPIM) :: INSEND = -9999
  INTEGER(KIND=JPIM) :: INRECV = -9999
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISENDTOT (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IRECVTOT (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISEND    (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IRECV    (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IINDEX(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: INDOFF(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IGPTRSEND(:,:,:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISETAL(:), ISETBL(:), ISETWL(:), ISETVL(:)
END TYPE TRXTOY_CTX

END MODULE TRXTOY_CTX_MOD
