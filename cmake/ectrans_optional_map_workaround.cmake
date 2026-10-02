# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# CCE 21's OpenMP offload runtime segfaults in is_contiguous_dv, called from
# cray_acc_new_transfer_list_1, when the transfer list of a TARGET construct
# contains the descriptor of an absent OPTIONAL assumed-shape dummy: it walks the
# null dope vector instead of skipping the item. OpenMP says a non-present
# optional in a map clause is ignored, so this is a runtime defect rather than
# anything the directives ask for.
#
# TRGTOL and TRLTOG trip it on every call that omits one of the PGP* gridpoint
# arrays, which covers all of benchmark call mode 2 and the partial-field
# field_view entry points: 138 of 366 tests segfault without this.
#
# ECTRANS_OPTIONAL_MAP_WORKAROUND makes those two routines map local pointers,
# aimed at the caller's arrays when supplied and at one-element stand-ins when
# not, so every mapped item carries a valid descriptor. See the constructs.
#
# Only CCE is known to be affected and only 21.0.0 has been tested here, so the
# gate stays narrow. Set ECTRANS_HAVE_OPTIONAL_MAP_ISSUE explicitly to override.
if( NOT DEFINED ECTRANS_HAVE_OPTIONAL_MAP_ISSUE )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
    set( ECTRANS_HAVE_OPTIONAL_MAP_ISSUE True )
  else()
    set( ECTRANS_HAVE_OPTIONAL_MAP_ISSUE False )
  endif()
endif()
