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

export GITHUB_REPOSITORY="${GITHUB_REPOSITORY:-ecmwf-ifs/ectrans}"
export GITHUB_SERVER_URL="${GITHUB_SERVER_URL:-https://github.com}"
export RUNNER_TEMP="${RUNNER_TEMP:-${PWD}}"

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

matching_job_hint() {
  local run_id=$1

  gh run view "${run_id}" --repo "${GITHUB_REPOSITORY}" --json jobs \
    | jq -r --arg name "${matrix_name}" --arg build_type "${build_type}" '
        .jobs[]
        | select((.name | contains($name)) and (.name | contains($build_type)))
        | [(.databaseId | tostring), .url, .name, .status, (.conclusion // "")]
        | @tsv
      ' \
    | head -n 1
}

artifact_download_url() {
  local run_id=$1

  gh api "repos/${GITHUB_REPOSITORY}/actions/runs/${run_id}/artifacts" \
    --jq ".artifacts[] | select(.name == \"${artifact_name}\") | .archive_download_url" \
    | head -n 1
}

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
  if run.get("headSha") == head_sha:
    print(run["databaseId"])
PY
)

if [[ -z "${run_ids}" ]]; then
  {
    echo "::error title=Reference workflow run not found::No ${workflow_file} workflow run was found for PR head ${head_sha}."
    echo "Run or rerun the workflow for branch ${branch}: ${workflow_branch_url}"
    echo "GitHub CLI: gh workflow run ${workflow_file} --repo ${GITHUB_REPOSITORY} --ref ${branch}"
  } >&2
  exit 1
fi

for run_id in ${run_ids}; do
  run_url="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${run_id}"
  archive_download_url=$(artifact_download_url "${run_id}" || true)

  echo "Checking reference artifact ${artifact_name} from run ${run_id}: ${run_url}" >&2
  if [[ -n "${archive_download_url}" ]]; then
    echo "Artifact download URL: ${archive_download_url}" >&2
  else
    echo "Artifact download URL: unavailable before gh run download" >&2
  fi
  echo "Download target directory: ${reused_dir}" >&2

  rm -rf "${reused_dir}"
  mkdir -p "${reused_dir}"
  if ! gh run download "${run_id}" --repo "${GITHUB_REPOSITORY}" --name "${artifact_name}" --dir "${reused_dir}"; then
    echo "Artifact ${artifact_name} was not downloaded from run ${run_id}." >&2
    continue
  fi

  cp -a "${reused_dir}" "${update_dir}"
  echo "Downloaded artifact ${artifact_name} to ${reused_dir}" >&2
  echo "Reference files are available at ${update_dir}" >&2
  exit 0
done

# Diagnostic information for the user to rerun the workflow or matrix job that should produce the artifact.

latest_run_id=$(printf '%s\n' ${run_ids} | head -n 1)
run_url="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${latest_run_id}"

job_hint=$(matching_job_hint "${latest_run_id}")

{
  echo "::error title=Reference artifact not found::No reusable artifact named ${artifact_name} was found for PR head ${head_sha}."
  echo "Rerun this workflow after the job that produces ${artifact_name} has uploaded the artifact."
  echo "Workflow for this branch: ${workflow_branch_url}"
  echo "Most recent matching run: ${run_url}"
  if [[ -n "${job_hint}" ]]; then
    IFS=$'\t' read -r job_id job_url job_name job_status job_conclusion <<<"${job_hint}"
    echo "Matching matrix job: ${job_name}"
    echo "Job status: ${job_status}${job_conclusion:+ (${job_conclusion})}"
    echo "Job URL: ${job_url}"
    if [[ "${job_status}" == "completed" ]]; then
      echo "GitHub CLI, one job: gh run rerun ${latest_run_id} --repo ${GITHUB_REPOSITORY} --job ${job_id}"
    fi
  else
    echo "GitHub CLI: gh run rerun ${latest_run_id} --repo ${GITHUB_REPOSITORY}"
  fi
} >&2
exit 1
