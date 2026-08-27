#!/usr/bin/env python3

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Compare two checksum files with a simple line-by-line diff.

Lines are compared after stripping any trailing comment starting with '#'.
This means differences in comment text alone are ignored.
"""

import argparse
import os
from pathlib import Path

SKIP_MISSING_REFERENCE = 77

class style:
    BOLD = "\033[1m"
    RED_BOLD = "\033[1;31m"
    GREEN_BOLD = "\033[1;32m"
    YELLOW_BOLD = "\033[1;33m"
    END = "\033[0m"


STYLE_MAP = {
    "none": "",
    "bold": style.BOLD,
    "red_bold": style.RED_BOLD,
    "green_bold": style.GREEN_BOLD,
    "yellow_bold": style.YELLOW_BOLD,
}


def strip_comment(line):
    return line.split("#", 1)[0].rstrip()


def split_comment(line):
    head, sep, tail = line.partition("#")
    if sep:
        return head, sep + tail
    return line.rstrip(), ""


def read_display_lines(path):
    with open(path, encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def read_normalized_lines(path):
    with open(path, encoding="utf-8") as handle:
        return [strip_comment(line) for line in handle]


def wrap_text(text, width):
    if width <= 0:
        return [""]

    lines = []
    for start in range(0, len(text), width):
        lines.append(text[start:start + width])

    if not lines:
        lines.append("")

    return lines


def display_path(path):
    try:
        return os.path.relpath(path)
    except ValueError:
        return path


def resolve_overrideable_reference_path(path, overrideable_reference_root):
    if not overrideable_reference_root:
        return path

    path = Path(path)
    configured_root = Path(overrideable_reference_root)
    reference_root = Path(os.environ.get("ECTRANS_TEST_REFERENCE_DIR") or configured_root)

    if not path.is_absolute():
        return str(reference_root / path)

    try:
        relative_path = path.relative_to(configured_root)
    except ValueError:
        return str(path)

    return str(reference_root / relative_path)


def differences(line, other_line):
    chars = []
    max_length = max(len(line), len(other_line))

    for index in range(max_length):
        char = line[index] if index < len(line) else ""
        other_char = other_line[index] if index < len(other_line) else ""
        chars.append((char, char != other_char))

    return chars


def render_differences(differences, style_name):
    parts = []
    in_style = False

    for char, is_changed in differences:
        if is_changed and not in_style:
            prefix = STYLE_MAP[style_name]
            if prefix:
                parts.append(prefix)
                in_style = True
        elif in_style and not is_changed:
            parts.append(style.END)
            in_style = False

        parts.append(char)

    if in_style:
        parts.append(style.END)

    return "".join(parts)


def parse_color_map(value):
    color_map = {"change": "yellow_bold"}
    if not value:
        return color_map

    for entry in value.split(","):
        entry = entry.strip()
        if not entry:
            continue
        if ":" not in entry:
            raise ValueError(f"Invalid color-map entry: {entry}")

        key, style_name = entry.split(":", 1)
        key = key.strip()
        style_name = style_name.strip()

        if key != "change":
            raise ValueError(f"Unsupported color-map key: {key}")
        if style_name not in STYLE_MAP:
            raise ValueError(f"Unsupported color-map style: {style_name}")

        color_map[key] = style_name

    return color_map


def render_display_line(line, other_line, is_match, width=0, change_style="yellow_bold"):
    visible_line, comment = split_comment(line)
    other_visible_line, _ = split_comment(other_line)

    if is_match:
        rendered_line = visible_line
    else:
        rendered_line = render_differences(
            differences(visible_line, other_visible_line),
            change_style,
        )

    rendered_line = rendered_line + comment
    padding = max(0, width - len(visible_line) - len(comment))
    return rendered_line + (" " * padding)


def group_selected_lines(selected_indices):
    if not selected_indices:
        return []

    groups = [[selected_indices[0]]]
    for index in selected_indices[1:]:
        if index == groups[-1][-1] + 1:
            groups[-1].append(index)
        else:
            groups.append([index])

    return groups


def print_selected_lines(first_lines, second_lines, selected_indices, mismatching_indices, left_width,
                         right_width, first_path, second_path, color_map):
    groups = group_selected_lines(selected_indices)
    line_number_width = 10
    left_header_width = max(1, left_width - 10)
    right_header_width = max(1, right_width - 10)

    left_header_lines = wrap_text(first_path, left_header_width)
    right_header_lines = wrap_text(second_path, right_header_width)
    nheader_lines = max(len(left_header_lines), len(right_header_lines))

    for index in range(nheader_lines):
        line_label = "line" if index == 0 else ""
        left_header = left_header_lines[index] if index < len(left_header_lines) else ""
        right_header = right_header_lines[index] if index < len(right_header_lines) else ""
        print(f"{'':<2}{line_label:<{line_number_width}} {left_header:<{left_width}} {right_header}")

    print()

    for group_number, group in enumerate(groups):
        if group_number > 0:
            print(f"{'':<2}{'...':<{line_number_width}}")

        for line_index in group:
            left_line = first_lines[line_index] if line_index < len(first_lines) else ""
            right_line = second_lines[line_index] if line_index < len(second_lines) else ""
            is_match = line_index not in mismatching_indices
            left_rendered = render_display_line(
                left_line,
                right_line,
                is_match=is_match,
                width=left_width,
                change_style=color_map["change"],
            )
            right_rendered = render_display_line(
                right_line,
                left_line,
                is_match=is_match,
                change_style=color_map["change"],
            )
            change_marker = ">" if not is_match else ""
            print(f"{change_marker:<2}{line_index + 1:<{line_number_width}} {left_rendered} {right_rendered}")


def truncate_selected_lines(selected_indices, mismatching_indices, max_output_lines):
    if max_output_lines is None or max_output_lines <= 0:
        return selected_indices, len(mismatching_indices), 0

    truncated_indices = selected_indices[:max_output_lines]
    shown_mismatches = sum(1 for index in truncated_indices if index in mismatching_indices)
    hidden_mismatches = len(mismatching_indices) - shown_mismatches
    return truncated_indices, len(mismatching_indices), hidden_mismatches


def print_truncation_message(total_diffs, hidden_diffs):
    message = (
        f"Output truncated: {total_diffs} total diffs, "
        f"{hidden_diffs} diffs not shown"
    )
    print(f"{style.RED_BOLD}{message}{style.END}")


def select_context_lines(mismatching_indices, total_lines, numlines, whole_file):
    if whole_file:
        return list(range(total_lines))

    selected = set()
    for mismatch_index in mismatching_indices:
        start = max(0, mismatch_index - numlines)
        end = min(total_lines, mismatch_index + numlines + 1)
        selected.update(range(start, end))

    return sorted(selected)


def line_contains_marker(lines, index, marker):
    return index < len(lines) and marker in lines[index]


def find_stop_index(first_display_lines, second_display_lines, mismatching_indices, marker):
    if not marker or not mismatching_indices:
        return None

    first_mismatch = mismatching_indices[0]
    max_lines = max(len(first_display_lines), len(second_display_lines))

    for index in range(first_mismatch + 1, max_lines):
        if line_contains_marker(first_display_lines, index, marker):
            return index
        if line_contains_marker(second_display_lines, index, marker):
            return index

    return None


def find_marker_index(first_display_lines, second_display_lines, marker):
    if not marker:
        return None

    max_lines = max(len(first_display_lines), len(second_display_lines))
    for index in range(max_lines):
        if line_contains_marker(first_display_lines, index, marker):
            return index
        if line_contains_marker(second_display_lines, index, marker):
            return index

    return None


def select_lines_before_marker(selected_indices, first_display_lines, second_display_lines, marker):
    stop_index = find_marker_index(first_display_lines, second_display_lines, marker)
    if stop_index is None:
        return selected_indices

    return [index for index in selected_indices if index < stop_index]


def compare_files(first_path, second_path, numlines, whole_file, show_matching, color_map, stop_marker,
                  max_output_lines):
    first_display_lines = read_display_lines(first_path)
    second_display_lines = read_display_lines(second_path)
    first_lines = read_normalized_lines(first_path)
    second_lines = read_normalized_lines(second_path)
    first_display_path = display_path(first_path)
    second_display_path = display_path(second_path)
    left_width = max((len(line) for line in first_display_lines), default=0) + 30
    right_width = max((len(line) for line in second_display_lines), default=0) + 30

    max_lines = max(len(first_lines), len(second_lines), len(first_display_lines), len(second_display_lines))
    mismatching_indices = []

    for index in range(max_lines):
        first_line = first_lines[index] if index < len(first_lines) else ""
        second_line = second_lines[index] if index < len(second_lines) else ""

        if first_line != second_line:
            mismatching_indices.append(index)

    if not mismatching_indices:
        if show_matching:
            selected_indices = list(range(max_lines))
            selected_indices = select_lines_before_marker(
                selected_indices,
                first_display_lines,
                second_display_lines,
                stop_marker,
            )
            print_selected_lines(
                first_display_lines,
                second_display_lines,
                selected_indices,
                set(),
                left_width,
                right_width,
                first_display_path,
                second_display_path,
                color_map,
            )
        return True

    selected_indices = select_context_lines(mismatching_indices, max_lines, numlines, whole_file)
    stop_index = find_stop_index(first_display_lines, second_display_lines, mismatching_indices, stop_marker)
    if stop_index is not None:
        selected_indices = [index for index in selected_indices if index < stop_index]

    selected_indices, total_diffs, hidden_diffs = truncate_selected_lines(
        selected_indices,
        set(mismatching_indices),
        max_output_lines,
    )

    print_selected_lines(
        first_display_lines,
        second_display_lines,
        selected_indices,
        set(mismatching_indices),
        left_width,
        right_width,
        first_display_path,
        second_display_path,
        color_map,
    )

    if hidden_diffs > 0:
        print_truncation_message(total_diffs, hidden_diffs)

    return False


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("first_file", help="First checksum file (reference)")
    parser.add_argument("second_file", help="Second checksum file")
    parser.add_argument(
        "-U", "--numlines",
        type=int,
        default=5,
        help="Number of matching context lines to show above and below mismatches"
    )
    parser.add_argument(
        "-W", "--whole-file",
        action="store_true",
        help="Show the whole file when mismatches are found"
    )
    parser.add_argument(
        "--show-matching",
        action="store_true",
        help="Show the whole file even when the normalized file contents match"
    )
    parser.add_argument(
        "--color-map",
        default="change:yellow_bold",
        help=(
            "Comma-separated mapping such as 'change:yellow_bold', 'change:bold', "
            "'change:red_bold', 'change:green_bold' or 'change:none'"
        )
    )
    parser.add_argument(
        "-s", "--stop-marker",
        metavar="TEXT",
        help=(
            "Stop printing output before the first line containing TEXT "
            "after the first mismatch. The marker is searched in original lines, "
            "including comments. The full files are still compared."
        )
    )
    parser.add_argument(
        "--overrideable-reference-root",
        help=(
            "Configured reference directory for the first file. When set, "
            "ECTRANS_TEST_REFERENCE_DIR can override this directory at runtime."
        )
    )
    parser.add_argument(
        "--max-output-lines",
        type=int,
        default=0,
        help=(
            "Maximum number of diff output lines to print when mismatches are found. "
            "A non-positive value disables truncation."
        )
    )
    args = parser.parse_args()

    try:
        color_map = parse_color_map(args.color_map)
    except ValueError as exc:
        parser.error(str(exc))

    first_file = resolve_overrideable_reference_path(args.first_file, args.overrideable_reference_root)

    if not Path(first_file).is_file():
        print(f"Skipping checksum comparison: missing reference file: {first_file}")
        raise SystemExit(SKIP_MISSING_REFERENCE)

    if compare_files(
        first_file,
        args.second_file,
        args.numlines,
        args.whole_file,
        args.show_matching,
        color_map,
        args.stop_marker,
        args.max_output_lines,
    ):
        raise SystemExit(0)

    raise SystemExit(1)


if __name__ == "__main__":
    main()