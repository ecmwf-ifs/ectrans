#!/usr/bin/env bash

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Convert a runner/compiler/build tuple into the stable directory and artifact key
# used for platform-specific checksum references.

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $(basename "$0") RUNNER_OS COMPILER BUILD_TYPE" >&2
  exit 2
fi

key="$1-$2-$3"
key=$(printf '%s' "${key}" | tr '[:upper:]' '[:lower:]')
key=$(printf '%s' "${key}" | sed -E 's/[^a-z0-9]+/-/g; s/^-+//; s/-+$//')

printf '%s\n' "${key}"