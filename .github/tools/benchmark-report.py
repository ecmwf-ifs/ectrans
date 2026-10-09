#!/usr/bin/env python3

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Render a Markdown benchmark comparison from the artifacts of a benchmark run.

Takes the directory that actions/download-artifact wrote into, containing one
subdirectory per platform leg (each with results.csv, runs.csv and
metadata.env as produced by .github/workflows/benchmark.yml), and writes a
pull request comment body to stdout.

The long-format results.csv is pivoted here rather than on disk, so that the
sweep can keep appending one row per run and a killed job still leaves usable
partial results behind.

Standard library only: the report job needs no pip install.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# Lets the report job find its own previous comment and update in place instead
# of posting a new one on every run. Changing this orphans existing comments.
MARKER = "<!-- ectrans-benchmark-report -->"

DEVELOP = "develop"
PR = "pr"


def read_metadata(path):
    """Parse a key=value file, tolerating absence and blank or malformed lines."""
    metadata = {}
    if not path.is_file():
        return metadata
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        metadata[key.strip()] = value.strip()
    return metadata


def read_rows(path):
    """Read a CSV as a list of dicts, returning [] if it is missing."""
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle) if row]


def to_float(value):
    """Floats that may legitimately be 'NA' when a run crashed."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def short_sha(sha):
    return sha[:7] if sha else "unknown"


def format_ms(value):
    return f"{value:.3f}" if value is not None else "n/a"


def format_delta(develop, pr, threshold):
    """Signed percentage change from develop to pr, flagged against a threshold.

    Returns the rendered cell and the signed percentage (None when either side
    is missing, so the caller can leave it out of the summary).
    """
    if develop is None or pr is None or develop == 0.0:
        return "—", None

    pct = (pr - develop) / develop * 100.0

    if pct <= -threshold:
        icon = "🟢"
    elif pct >= threshold:
        icon = "🔴"
    else:
        icon = "⚪"

    sign = "+" if pct > 0 else "−" if pct < 0 else ""
    text = f"{sign}{abs(pct):.1f}%"
    if abs(pct) >= threshold:
        text = f"**{text}**"

    return f"{icon} {text}", pct


def describe_case(runs):
    """One-line description of the benchmark case, taken from the runs index."""
    if not runs:
        return None

    first = runs[0]
    parts = []

    precision = first.get("precision")
    truncations = {row.get("truncation") for row in runs if row.get("truncation")}
    if len(truncations) == 1:
        parts.append(f"T{truncations.pop()}")
    if precision:
        parts.append(precision)

    for label, key in (("nfld", "nfld"), ("nlev", "nlev"), ("niter", "niter")):
        if first.get(key):
            parts.append(f"{label}={first[key]}")

    ranks = first.get("ntasks_per_node")
    threads = first.get("cpus_per_task")
    if ranks and threads:
        parts.append(f"{ranks} ranks/node × {threads} threads")

    return " · ".join(parts) if parts else None


def pivot(results):
    """Long results rows -> {(nodes, truncation): {version: median_ms}}.

    Rows with an unparseable node count are dropped: they cannot be placed in
    the table, and silently mispositioning them would be worse.
    """
    table = {}
    for row in results:
        nodes = to_int(row.get("nodes"))
        if nodes is None:
            continue
        key = (nodes, row.get("truncation") or "")
        table.setdefault(key, {})[row.get("version")] = to_float(row.get("median_ms"))
    return table


def render_leg(directory, threshold):
    """Render one platform leg. Returns (lines, deltas, ok)."""
    results = read_rows(directory / "results.csv")
    runs = read_rows(directory / "runs.csv")
    metadata = read_metadata(directory / "metadata.env")

    name = metadata.get("reference_key") or directory.name
    lines = [f"### `{name}`", ""]

    if not results:
        lines += [
            ":warning: No timings were produced for this platform. "
            "See the workflow run for the build and sweep logs.",
            "",
        ]
        return lines, [], False

    case = describe_case(runs)
    if case:
        lines += [case, ""]

    table = pivot(results)
    truncations = {truncation for _, truncation in table}
    show_truncation = len(truncations) > 1

    if show_truncation:
        lines += [
            "| nodes | truncation | develop (ms) | this PR (ms) | Δ |",
            "|------:|-----------:|-------------:|-------------:|:--|",
        ]
    else:
        lines += [
            "| nodes | develop (ms) | this PR (ms) | Δ |",
            "|------:|-------------:|-------------:|:--|",
        ]

    deltas = []
    for (nodes, truncation), values in sorted(table.items()):
        develop = values.get(DEVELOP)
        pr = values.get(PR)
        cell, pct = format_delta(develop, pr, threshold)
        if pct is not None:
            deltas.append((nodes, pct))

        row = [str(nodes)]
        if show_truncation:
            row.append(f"T{truncation}" if truncation else "—")
        row += [format_ms(develop), format_ms(pr), cell]
        lines.append("| " + " | ".join(row) + " |")

    lines.append("")

    failed = [row for row in runs if to_int(row.get("exit_code")) not in (0, None)]
    if failed:
        detail = ", ".join(
            f"`{row.get('logfile') or '?'}` (exit {row.get('exit_code')})" for row in failed
        )
        lines += [f":warning: {len(failed)} run(s) failed: {detail}", ""]

    return lines, deltas, not failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "artifacts",
        type=Path,
        help="Directory containing one subdirectory per benchmark artifact",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=5.0,
        help="Percentage change flagged as significant (default: 5)",
    )
    parser.add_argument(
        "--run-url",
        default=os.environ.get("BENCHMARK_RUN_URL", ""),
        help="Link back to the workflow run",
    )
    args = parser.parse_args()

    # download-artifact writes one subdirectory per artifact, but accept a bare
    # results.csv too so the script can be pointed at a single unpacked leg.
    if (args.artifacts / "results.csv").is_file():
        directories = [args.artifacts]
    elif args.artifacts.is_dir():
        directories = sorted(d for d in args.artifacts.iterdir() if d.is_dir())
    else:
        directories = []

    out = [MARKER, "", "## ecTrans benchmark", ""]

    if not directories:
        out += [
            ":warning: No benchmark artifacts were produced, so there is nothing "
            "to compare. The build or the sweep most likely failed before writing "
            "any results.",
            "",
        ]
        if args.run_url:
            out += [f"[Workflow run]({args.run_url})", ""]
        sys.stdout.write("\n".join(out))
        return 0

    # Provenance is identical across legs; take it from the first that has any.
    metadata = {}
    for directory in directories:
        metadata = read_metadata(directory / "metadata.env")
        if metadata:
            break

    pr_repo = metadata.get("pr_repo", "")
    pr_sha = short_sha(metadata.get("pr_sha", ""))
    base_ref = metadata.get("base_ref", "base")
    base_sha = short_sha(metadata.get("base_sha", ""))
    if pr_repo:
        out += [f"`{pr_repo}@{pr_sha}` vs `{base_ref}@{base_sha}`", ""]

    all_deltas = []
    for directory in directories:
        lines, deltas, _ = render_leg(directory, args.threshold)
        out += lines
        all_deltas += deltas

    if all_deltas:
        worst_nodes, worst = max(all_deltas, key=lambda item: item[1])
        best_nodes, best = min(all_deltas, key=lambda item: item[1])
        if worst >= args.threshold:
            out += [
                f"**Slowest case is {worst:+.1f}% on {worst_nodes} node(s)** "
                f"(threshold ±{args.threshold:.0f}%).",
                "",
            ]
        elif best <= -args.threshold:
            out += [f"Fastest case is {best:+.1f}% on {best_nodes} node(s).", ""]
        else:
            out += [f"No case moved by more than ±{args.threshold:.0f}%.", ""]

    footer = []
    section = metadata.get("timing_section")
    if section:
        footer.append(f"{section}, median time step")
    if metadata.get("arch"):
        footer.append(f"`{metadata['arch']}`")
    if args.run_url:
        footer.append(f"[workflow run]({args.run_url})")
    if footer:
        out += ["<sub>" + " · ".join(footer) + "</sub>", ""]

    sys.stdout.write("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
