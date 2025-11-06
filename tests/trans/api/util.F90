! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE UTIL

IMPLICIT NONE

PRIVATE
PUBLIC :: DETECT_MPIRUN

CONTAINS

LOGICAL FUNCTION DETECT_MPIRUN()
  USE EC_ENV_MOD, ONLY : EC_PUTENV

  INTEGER :: ILEN
  CHARACTER(LEN=32), DIMENSION(4) :: CMPIRUN_DETECT
  CHARACTER(LEN=4) :: CLENV
  INTEGER :: IVAR

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  CMPIRUN_DETECT(1) = 'OMPI_COMM_WORLD_SIZE'  ! OPENMPI
  CMPIRUN_DETECT(2) = 'ALPS_APP_PE'           ! CRAY PE
  CMPIRUN_DETECT(3) = 'PMI_SIZE'              ! INTEL
  CMPIRUN_DETECT(4) = 'SLURM_NTASKS'          ! SLURM

  DETECT_MPIRUN = .FALSE.
  DO IVAR = 1, 4
    CALL GET_ENVIRONMENT_VARIABLE(NAME=TRIM(CMPIRUN_DETECT(IVAR)), LENGTH=ILEN)
    IF (ILEN > 0) THEN
      DETECT_MPIRUN = .TRUE.
      EXIT ! Break
    ENDIF
  ENDDO

  CALL GET_ENVIRONMENT_VARIABLE(NAME="ECTRANS_USE_MPI", VALUE=CLENV, LENGTH=ILEN )
  IF (ILEN > 0) THEN
    DETECT_MPIRUN = .TRUE.
    IF (TRIM(CLENV) == "0" .OR. TRIM(CLENV) == "OFF" .OR. TRIM(CLENV) == "OFF" &
      & .OR. TRIM(CLENV) == "F") THEN
      DETECT_MPIRUN = .FALSE.
    ENDIF
    CALL EC_PUTENV("DR_HOOK_ASSERT_MPI_INITIALIZED=0", OVERWRITE=.TRUE.)
  ENDIF
END FUNCTION

END MODULE UTIL