// (C) Copyright 2026- ECMWF.
// (C) Copyright 2026- NVIDIA.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#pragma once

extern "C" void growing_allocator_register_free_c(void *,
                                                  void (&)(float *, size_t));
