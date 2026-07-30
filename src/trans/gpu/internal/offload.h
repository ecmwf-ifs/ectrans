! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http:!www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

#ifdef ACCGPU
#define _DIR !$ACC
#define _DATA DATA
#define _END_DATA END DATA
#define _EXIT_DATA EXIT DATA
#define _PRESENT(...) PRESENT(__VA_ARGS__)
#define _DELETE(...) DELETE(__VA_ARGS__)
#ifdef _CRAYFTN
#define _ASYNC(...)
#else
#define _ASYNC(...) ASYNC(__VA_ARGS__)
#endif
#define _WAIT(...) WAIT(__VA_ARGS__)
#define _PARALLEL PARALLEL LOOP
#define _FIRSTPRIVATE(...) FIRSTPRIVATE(__VA_ARGS__)
#define _SHARED(...)
#endif
#ifdef OMPGPU
#define _DIR !$OMP
#define _DATA TARGET DATA
#define _END_DATA END TARGET DATA
#define _PRESENT(...) MAP(PRESENT,ALLOC: __VA_ARGS__)
#define _EXIT_DATA TARGET EXIT DATA
#define _DELETE(...) MAP(DELETE: __VA_ARGS__)
#define _ASYNC(...)
#define _WAIT(...)
#define _PARALLEL TARGET TEAMS DISTRIBUTE PARALLEL DO
#define _FIRSTPRIVATE(...) MAP(TO: __VA_ARGS__)
#define _SHARED(...) MAP(TO: __VA_ARGS__)
#endif
