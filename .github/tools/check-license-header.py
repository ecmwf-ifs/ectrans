#!/usr/bin/env python3
"""Check that files contain the ECMWF Apache-2.0 license header as a contiguous block.

The header is stored once in ``license-header.template`` with a ``{comment}``
placeholder standing in for the comment marker. For each file, the placeholder is
replaced with the comment marker appropriate to the file's extension (e.g. ``//``
for C/C++, ``!`` for Fortran, ``#`` for Python/CMake) and the resulting block of
lines is searched for, contiguously and in order, within the file.

Usage:
    check-license-header.py FILE [FILE ...]

Exits non-zero if any file is missing the header block.
"""

from __future__ import annotations

import os
import sys

# Map a file extension (lower-case, with leading dot) to its line-comment marker.
_EXT_COMMENT = {
    # Fortran
    ".f90": "!", ".f": "!", ".f77": "!", ".inc": "!",
    # CMake / Python
    ".cmake": "#", ".py": "#",
    # C / C++ and headers
    ".c": "//", ".cpp": "//", ".cc": "//", ".cxx": "//",
    ".h": "//", ".hpp": "//", ".hh": "//", ".hxx": "//",
}

_TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "license-header.template")


def _comment_marker(path: str) -> str | None:
    if os.path.basename(path) == "CMakeLists.txt":
        return "#"
    # Extensions are matched case-insensitively (.F90 == .f90).
    return _EXT_COMMENT.get(os.path.splitext(path)[1].lower())


def _load_template() -> list[str]:
    with open(_TEMPLATE_PATH, "r", encoding="utf-8") as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


def _expected_block(template: list[str], marker: str) -> list[str]:
    return [line.replace("{comment}", marker) for line in template]


def _has_header_block(path: str, template: list[str]) -> bool:
    marker = _comment_marker(path)
    if marker is None:
        print(f"::error file={path}::No comment marker known for {path}")
        return False

    expected = _expected_block(template, marker)

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            lines = [line.rstrip("\n").rstrip() for line in fh]
    except OSError as exc:
        print(f"::error file={path}::Could not read {path}: {exc}")
        return False

    needed = len(expected)
    for start in range(0, len(lines) - needed + 1):
        if lines[start : start + needed] == expected:
            return True
    return False


def main(argv: list[str]) -> int:
    files = argv[1:]
    template = _load_template()
    missing = [f for f in files if not _has_header_block(f, template)]

    print(f"Checked {len(files)} newly added file(s).")

    if missing:
        print("::error::The following newly added files are missing the required "
              "license header (it must appear as a contiguous block):")
        for f in missing:
            print(f"::error file={f}::Missing license header in {f}")
            print(f"  - {f}")
        return 1

    print("All newly added files contain the required license header.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
