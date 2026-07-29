#!/usr/bin/env bash

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Download a checksum reference candidate artifact for one expected platform key.
# If the artifact is unavailable, print links and GitHub CLI commands for rerunning
# the workflow or matrix job that should produce it.

set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: $(basename "$0") WORKFLOW_FILE BRANCH HEAD_SHA ARTIFACT_NAME MATRIX_NAME BUILD_TYPE" >&2
  exit 2
fi

workflow_file=$1
branch=$2
head_sha=$3
artifact_name=$4
matrix_name=$5
build_type=$6

update_dir="${RUNNER_TEMP}/ectrans-reference-update"
reused_dir="${RUNNER_TEMP}/ectrans-reference-reused"
workflow_branch_url=$(python3 - "${workflow_file}" "${branch}" <<'PY'
import os
import sys
import urllib.parse

workflow_file, branch = sys.argv[1:]
base = f"{os.environ['GITHUB_SERVER_URL']}/{os.environ['GITHUB_REPOSITORY']}/actions/workflows/{workflow_file}"
query = urllib.parse.quote(f"branch:{branch}", safe="")
print(f"{base}?query={query}")
PY
)

rm -rf "${update_dir}" "${reused_dir}"

runs_json=$(gh run list \
  --repo "${GITHUB_REPOSITORY}" \
  --workflow "${workflow_file}" \
  --branch "${branch}" \
  --limit 20 \
  --json databaseId,status,headSha,url)

run_ids=$(RUNS_JSON="${runs_json}" python3 - "${head_sha}" <<'PY'
import json
import os
import sys

head_sha = sys.argv[1]
for run in json.loads(os.environ["RUNS_JSON"]):
    if run.get("status") == "completed" and run.get("headSha") == head_sha:
        print(run["databaseId"])
PY
)

if [[ -z "${run_ids}" ]]; then
  {
    echo "::error title=Reference workflow run not found::No completed ${workflow_file} workflow run was found for PR head ${head_sha}."
    echo "Run or rerun the workflow for branch ${branch}: ${workflow_branch_url}"
    echo "GitHub CLI: gh workflow run ${workflow_file} --repo ${GITHUB_REPOSITORY} --ref ${branch}"
  } >&2
  exit 1
fi

for run_id in ${run_ids}; do
  rm -rf "${reused_dir}"
  mkdir -p "${reused_dir}"
  if ! gh run download "${run_id}" --repo "${GITHUB_REPOSITORY}" --name "${artifact_name}" --dir "${reused_dir}"; then
    continue
  fi

  cp -a "${reused_dir}" "${update_dir}"
  exit 0
done

latest_run_id=$(printf '%s\n' ${run_ids} | head -n 1)
run_url="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${latest_run_id}"
job_hint=$(gh run view "${latest_run_id}" --repo "${GITHUB_REPOSITORY}" --json jobs \
  | jq -r --arg name "${matrix_name}" --arg build_type "${build_type}" '
      .jobs[]
      | select((.name | contains($name)) and (.name | contains($build_type)))
      | "\(.databaseId)\t\(.url)\t\(.name)"
    ' \
  | head -n 1)

{
  echo "::error title=Reference artifact not found::No reusable artifact named ${artifact_name} was found for PR head ${head_sha}."
  echo "Rerun the job that produces ${artifact_name}, then rerun this workflow."
  echo "Workflow for this branch: ${workflow_branch_url}"
  echo "Most recent matching run: ${run_url}"
  echo "GitHub CLI: gh run rerun ${latest_run_id} --repo ${GITHUB_REPOSITORY}"
  if [[ -n "${job_hint}" ]]; then
    IFS=$'\t' read -r job_id job_url job_name <<<"${job_hint}"
    echo "Matching matrix job: ${job_name}"
    echo "Job URL: ${job_url}"
    echo "GitHub CLI, one job: gh run rerun ${latest_run_id} --repo ${GITHUB_REPOSITORY} --job ${job_id}"
  fi
} >&2
exit 1
