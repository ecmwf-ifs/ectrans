#!/usr/bin/env bash

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Package generated checksum reference files into the artifact layout consumed by
# update-references.yml, including metadata that ties the artifact to a PR head SHA
# and platform reference key.

set -euo pipefail

if [[ $# -ne 4 && $# -ne 5 ]]; then
  echo "Usage: $(basename "$0") REFERENCE_KEY HEAD_SHA PR_NUMBER OUTPUT_DIR [BUILD_PROVENANCE]" >&2
  exit 2
fi

reference_key=$1
head_sha=$2
pr_number=$3
output_dir=$4
build_provenance=${5:-}

metadata_file="${output_dir}/metadata.json"
references_dir="${output_dir}/references"

rm -rf "${output_dir}"
mkdir -p "${references_dir}"

status="missing-update-script"
message="No configured ectrans-update-references.sh script was found"
reference_count=0

while IFS= read -r update_script; do
  status="update-failed"
  message="Failed to package generated checksum references"
  rm -rf "${references_dir}"
  mkdir -p "${references_dir}"

  if ! bash "${update_script}" --no-anonymous "${references_dir}"; then
    continue
  fi

  reference_count=$(find "${references_dir}" -type f -name '*.checksums' | wc -l | tr -d ' ')
  if [[ "${reference_count}" != "0" ]]; then
    status="ok"
    message="references packaged"
    break
  fi

  status="empty"
  message="No checksum reference files were packaged"
done < <(find "${GITHUB_WORKSPACE:-$PWD}" -path '*/tests/ectrans-update-references.sh' -type f | sort)

reference_count=$(find "${references_dir}" -type f -name '*.checksums' | wc -l | tr -d ' ')
if [[ -z "${build_provenance}" && -f "${references_dir}/metadata.json" ]]; then
  build_provenance="${references_dir}/metadata.json"
fi

REFERENCE_KEY="${reference_key}" HEAD_SHA="${head_sha}" PR_NUMBER="${pr_number}" \
STATUS="${status}" MESSAGE="${message}" REFERENCE_COUNT="${reference_count}" \
python3 - "${metadata_file}" "${build_provenance}" <<'PY'
import json
import os
import sys

metadata_file, provenance_file = sys.argv[1:]
provenance = {}
if provenance_file:
  with open(provenance_file, encoding="utf-8") as handle:
    provenance = json.load(handle)

metadata = {
  "schema": 1,
  "reference_key": os.environ["REFERENCE_KEY"],
  "head_sha": os.environ["HEAD_SHA"],
  "pr_number": os.environ["PR_NUMBER"],
  "runner_os": os.environ.get("RUNNER_OS", ""),
  "runner_arch": os.environ.get("RUNNER_ARCH", ""),
  "status": os.environ["STATUS"],
  "message": os.environ["MESSAGE"],
  "reference_count": int(os.environ["REFERENCE_COUNT"]),
  "build_provenance": provenance,
}
with open(metadata_file, "w", encoding="utf-8") as handle:
  json.dump(metadata, handle, indent=2, sort_keys=True)
  handle.write("\n")
PY

cat "${metadata_file}"

if [[ "${status}" != "ok" ]]; then
  exit 0
fi