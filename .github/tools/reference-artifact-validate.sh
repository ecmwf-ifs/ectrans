#!/usr/bin/env bash

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Validate a packaged checksum reference artifact before it is reused or committed,
# ensuring its metadata matches the expected platform key and PR head SHA.

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $(basename "$0") ARTIFACT_DIR EXPECTED_REFERENCE_KEY EXPECTED_HEAD_SHA" >&2
  exit 2
fi

artifact_dir=$1
expected_reference_key=$2
expected_head_sha=$3
metadata_file="${artifact_dir}/metadata.json"
references_dir="${artifact_dir}/references"

if [[ ! -f "${metadata_file}" ]]; then
  echo "Missing artifact metadata: ${metadata_file}" >&2
  exit 1
fi

if [[ ! -d "${references_dir}" ]]; then
  echo "Missing artifact reference directory: ${references_dir}" >&2
  exit 1
fi

python3 - "${metadata_file}" "${expected_reference_key}" "${expected_head_sha}" <<'PY'
import json
import sys

metadata_file, expected_reference_key, expected_head_sha = sys.argv[1:]
with open(metadata_file, encoding="utf-8") as handle:
    metadata = json.load(handle)

errors = []
if metadata.get("schema") != 1:
    errors.append("metadata schema is not 1")
if metadata.get("reference_key") != expected_reference_key:
    errors.append(f"reference_key is {metadata.get('reference_key')!r}, expected {expected_reference_key!r}")
if metadata.get("head_sha") != expected_head_sha:
    errors.append(f"head_sha is {metadata.get('head_sha')!r}, expected {expected_head_sha!r}")
if metadata.get("status") != "ok":
    errors.append(f"status is {metadata.get('status')!r}, expected 'ok'")
if int(metadata.get("reference_count", 0)) <= 0:
    errors.append("reference_count is not positive")
if not isinstance(metadata.get("build_provenance"), dict):
  errors.append("build_provenance is missing or is not an object")

if errors:
    for error in errors:
        print(error, file=sys.stderr)
    raise SystemExit(1)
PY

if ! find "${references_dir}" -type f -name '*.checksums' | grep -q .; then
  echo "Artifact contains no .checksums files" >&2
  exit 1
fi