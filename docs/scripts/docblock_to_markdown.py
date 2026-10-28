#!/usr/bin/env python3
#
# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""
Generate Markdown API documentation from the Doxygen-style docblocks in the ecTrans interface
headers (src/trans/include/ectrans/*.h).

A docblock is a block of comment lines starting with "!>" (first line) or "!!" (subsequent lines)
placed immediately before the SUBROUTINE statement. It contains:

  - Free text describing the subroutine (paragraphs, lists, and "@note" paragraphs).
  - An optional "@deprecated" paragraph.
  - One "@param NAME description" entry per dummy argument. The description may continue on
    following (indented) lines until the next "@" command or a blank line.

Within a parameter description, lines beginning with the following are kept on their own line and
formatted specially:

  Default: VALUE ...         -> *Default*: `VALUE` ...
  Dimensions: (...)          -> Dimensions: (...)
  Potentially deprecatable.  -> **Potentially deprecatable.**
  TODO: ...                  -> **TODO:** ...

Lines beginning with "1." or "-" are kept as list items. All other lines are joined into
paragraphs.

The signature and the type, attributes, and shape of each argument are taken from the Fortran
source itself, so they are not repeated in the docblock. Argument names, the names of other API
subroutines, logical literals and identifiers containing "_" or "%" are wrapped in backticks
automatically.

Usage: docblock_to_markdown.py FILE [FILE ...]
"""

import re
import sys
import textwrap
from pathlib import Path

LINE_WIDTH = 100

# Lines within a parameter description which are rendered on their own line
MARKER_RE = re.compile(r"^(Default:|Dimensions:|Potentially deprecatable\.|TODO:)")
LIST_ITEM_RE = re.compile(r"^(\d+\.|-)\s")
SUBROUTINE_RE = re.compile(r"^\s*SUBROUTINE\s+(\w+)\s*\(", re.IGNORECASE)
END_SUBROUTINE_RE = re.compile(r"^\s*END\s*SUBROUTINE", re.IGNORECASE)

# Order in which argument sections are listed: (optional, intent)
SECTION_ORDER = [
    (False, "IN"), (False, "INOUT"), (False, "OUT"),
    (True, "IN"), (True, "INOUT"), (True, "OUT"), (True, "EXTERNAL"), (False, "EXTERNAL"),
]


class DocblockError(Exception):
    pass


# ----------------------------------------------------------------------------------------------
# Fortran parsing
# ----------------------------------------------------------------------------------------------

def split_top_level(text, sep=","):
    """Split text on sep, ignoring separators nested within parentheses."""
    parts, depth, current = [], 0, ""
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == sep and depth == 0:
            parts.append(current.strip())
            current = ""
        else:
            current += char
    parts.append(current.strip())
    return [p for p in parts if p]


def join_statements(lines):
    """Strip comments and join continuation lines into complete Fortran statements."""
    statements, current = [], ""
    for line in lines:
        code = line.split("!", 1)[0].rstrip()
        if not code.strip():
            continue
        code = code.strip()
        if code.startswith("&"):
            code = code[1:]
        if code.endswith("&"):
            current += code[:-1]
        else:
            statements.append(current + code)
            current = ""
    if current:
        statements.append(current)
    return statements


class Argument:
    def __init__(self, name):
        self.name = name
        self.type = None
        self.dims = ""
        self.attrs = []

    @property
    def optional(self):
        return "OPTIONAL" in self.attrs

    @property
    def intent(self):
        if "EXTERNAL" in self.attrs:
            return "EXTERNAL"
        for attr in self.attrs:
            match = re.match(r"INTENT\((\w+)\)", attr)
            if match:
                return match.group(1)
        return ""

    def declaration(self):
        """Reconstruct a normalised declaration, e.g. 'REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: X(:)'."""
        attrs = sorted(self.attrs, key=lambda a: (a != "OPTIONAL", not a.startswith("INTENT")))
        parts = ([self.type] if self.type else []) + attrs
        return f"{', '.join(parts)} :: {self.name}{self.dims}"


def parse_fortran(lines, start):
    """Parse the subroutine starting at line index start. Returns (name, ordered arguments)."""
    end = next(i for i in range(start, len(lines)) if END_SUBROUTINE_RE.match(lines[i]))
    statements = join_statements(lines[start:end])

    header = re.match(r"SUBROUTINE\s+(\w+)\s*\((.*)\)", statements[0], re.IGNORECASE)
    name = header.group(1).upper()
    args = {a.upper(): Argument(a.upper()) for a in split_top_level(header.group(2))}

    for statement in statements[1:]:
        if "::" in statement:
            spec, entities = statement.split("::", 1)
            spec = split_top_level(spec)
            type_spec = re.sub(r"\s+", "", spec[0]).upper()
            attrs = [re.sub(r"\s+", "", a).upper() for a in spec[1:]]
            for entity in split_top_level(entities):
                match = re.match(r"(\w+)\s*(\(.*\))?", entity)
                arg = args.get(match.group(1).upper())
                if arg is None:
                    continue
                arg.type = type_spec
                arg.dims = re.sub(r"\s+", "", match.group(2) or "")
                arg.attrs += [a for a in attrs if a not in arg.attrs]
        else:
            # Attribute statements, e.g. "EXTERNAL FSPGL_PROC" or "OPTIONAL FSPGL_PROC"
            match = re.match(r"(EXTERNAL|OPTIONAL)\s+(.*)", statement, re.IGNORECASE)
            if match:
                for entity in split_top_level(match.group(2)):
                    arg = args.get(entity.upper())
                    if arg and match.group(1).upper() not in arg.attrs:
                        arg.attrs.append(match.group(1).upper())

    return name, list(args.values())


# ----------------------------------------------------------------------------------------------
# Docblock parsing
# ----------------------------------------------------------------------------------------------

class Docblock:
    def __init__(self):
        self.description = []  # Lines of free text
        self.deprecated = []   # Lines of the @deprecated paragraph
        self.params = {}       # Argument name -> list of description lines


def extract_docblock(lines, subroutine_line):
    """Return the comment lines (without the comment markers) preceding the SUBROUTINE statement."""
    start = subroutine_line
    while start > 0 and re.match(r"^\s*!!", lines[start - 1]):
        start -= 1
    if start == 0 or not re.match(r"^\s*!>", lines[start - 1]):
        raise DocblockError("no docblock (starting with '!>') found before SUBROUTINE statement")
    start -= 1
    return [re.sub(r"^\s*!(>|!) ?", "", line).rstrip() for line in lines[start:subroutine_line]]


def parse_docblock(doc_lines):
    doc = Docblock()
    target = doc.description
    for line in doc_lines:
        command = re.match(r"^\s*@(\w+)\s*(.*)$", line)
        if command and command.group(1) == "param":
            match = re.match(r"(\w+)\s*(.*)", command.group(2))
            if not match:
                raise DocblockError("@param without an argument name")
            name = match.group(1).upper()
            if name in doc.params:
                raise DocblockError(f"argument {name} documented more than once")
            target = doc.params[name] = [match.group(2)] if match.group(2) else []
        elif command and command.group(1) == "deprecated":
            target = doc.deprecated
            target.append(command.group(2))
        elif target is not doc.description and not line.strip():
            # A blank line terminates a @param or @deprecated paragraph
            target = doc.description
        elif target is doc.description:
            target.append(line)
        else:
            target.append(line.strip() if not LIST_ITEM_RE.match(line.strip()) else line)
    return doc


# ----------------------------------------------------------------------------------------------
# Markdown generation
# ----------------------------------------------------------------------------------------------

class Formatter:
    def __init__(self, identifiers):
        names = "|".join(sorted(map(re.escape, identifiers), key=len, reverse=True))
        self.code_re = re.compile(
            r"(?<![\w%.`])("
            + (names + "|" if names else "")
            + r"\.TRUE\.|\.FALSE\.|[A-Z]\w*%\w+|[A-Z][A-Z0-9]*_[A-Z0-9_]*[A-Z0-9])(?![\w%`])"
        )

    def code(self, text):
        """Wrap identifiers in backticks, leaving existing code spans untouched."""
        pieces = re.split(r"(`[^`]*`|\]\([^)]*\))", text)
        return "".join(
            p if p.startswith(("`", "](")) else self.code_re.sub(r"`\1`", p) for p in pieces
        )

    def description(self, lines):
        """Format free text, converting Doxygen @note paragraphs to FORD's @note/@endnote."""
        out, in_note = [], False
        for line in lines + [""]:
            note = re.match(r"^\s*@note\s*(.*)$", line)
            if note:
                out.append("@note")
                if note.group(1):
                    out.append(self.code(note.group(1)))
                in_note = True
            elif in_note and (not line.strip() or line.strip() == "@endnote"):
                out += ["@endnote", ""]
                in_note = False
            elif line.strip() == "@endnote":
                continue
            else:
                out.append(self.code(line))
        while out and not out[-1].strip():
            out.pop()
        return out

    def marker(self, line):
        if line.startswith("Default:"):
            value = line[len("Default:"):].strip()
            first, _, rest = value.partition(" ")
            if re.match(r"^([-+]?\d[\w.]*|\.[A-Z]+\.|'[^']*'|[A-Z][A-Z0-9_%]*)$", first):
                value = f"`{first}`" + (f" {self.code(rest)}" if rest else "")
            else:
                value = self.code(value)
            return f"*Default*: {value}"
        if line.startswith("TODO:"):
            return f"**TODO:** {self.code(line[len('TODO:'):].strip())}"
        if line.startswith("Potentially deprecatable."):
            return f"**{line}**"
        return self.code(line)

    def param(self, lines):
        """Format a parameter description as a list of Markdown lines (without indentation)."""
        blocks = []  # (kind, text)
        for line in lines:
            if MARKER_RE.match(line):
                blocks.append(("marker", self.marker(line)))
            elif LIST_ITEM_RE.match(line.strip()):
                blocks.append(("list", self.code(line.strip())))
            elif blocks and blocks[-1][0] == "text":
                blocks[-1] = ("text", blocks[-1][1] + " " + self.code(line))
            elif line:
                blocks.append(("text", self.code(line)))

        out = []
        for i, (kind, text) in enumerate(blocks):
            following = blocks[i + 1][0] if i + 1 < len(blocks) else None
            indent = "    " if kind == "list" else "  "
            wrapped = textwrap.wrap(text, LINE_WIDTH, initial_indent=indent,
                                    subsequent_indent=indent + ("   " if kind == "list" else ""),
                                    break_long_words=False, break_on_hyphens=False)
            # Hard line break unless this is the last line or the next block is a list
            if following not in (None, "list"):
                wrapped[-1] += "  "
            out += wrapped
        return out


def signature(name, args):
    prefix = f"SUBROUTINE {name}("
    indent = "  &" + " " * (len(prefix) - 3)
    lines, current = [], prefix
    for i, arg in enumerate(args):
        token = arg + (", " if i < len(args) - 1 else ")")
        if len((current + token).rstrip()) + 2 > LINE_WIDTH and current.strip() not in ("", "&"):
            lines.append(current.rstrip() + " &")
            current = indent
        current += token
    lines.append(current)
    return lines


def section_heading(optional, intent):
    attrs = (["OPTIONAL"] if optional else []) + (
        ["EXTERNAL"] if intent == "EXTERNAL" else [f"INTENT({intent})"] if intent else []
    )
    return f"### `{', '.join(attrs)}` arguments"


def to_markdown(path, api_names=()):
    lines = Path(path).read_text().splitlines()
    sub_line = next((i for i, l in enumerate(lines) if SUBROUTINE_RE.match(l)), None)
    if sub_line is None:
        raise DocblockError("no SUBROUTINE statement found")
    name, args = parse_fortran(lines, sub_line)
    doc = parse_docblock(extract_docblock(lines, sub_line))

    arg_names = [a.name for a in args]
    for documented in doc.params:
        if documented not in arg_names:
            raise DocblockError(f"@param {documented} is not an argument of {name}")
    if not doc.deprecated:
        for arg in arg_names:
            if arg not in doc.params:
                print(f"{path}: warning: argument {arg} is not documented", file=sys.stderr)

    fmt = Formatter(set(arg_names) | set(api_names) | {name})
    out = [f"## `{name}`", ""]

    if doc.deprecated:
        text = fmt.code(" ".join(line.strip() for line in doc.deprecated))
        out += ["@note", *textwrap.wrap(text, LINE_WIDTH), "@endnote", ""]

    out += ["### Signature", "", "```f90", *signature(name, arg_names), "```", ""]

    description = fmt.description(doc.description)
    while description and not description[0].strip():
        description.pop(0)
    if description:
        out += ["### Purpose", "", *description, ""]

    # Documented arguments first, in docblock order, then any undocumented arguments. Sections
    # appear in the order in which their first argument is documented.
    by_name = {a.name: a for a in args}
    ordered = [by_name[n] for n in doc.params] + [a for a in args if a.name not in doc.params]
    sections = list(dict.fromkeys((a.optional, a.intent) for a in ordered))
    if not doc.params:
        sections.sort(key=lambda s: SECTION_ORDER.index(s) if s in SECTION_ORDER else len(SECTION_ORDER))
    for optional, intent in sections:
        section = [a for a in ordered if (a.optional, a.intent) == (optional, intent)]
        out += [section_heading(optional, intent), ""]
        for arg in section:
            body = fmt.param(doc.params.get(arg.name, []))
            out.append(f"- `{arg.declaration()}`" + ("  " if body else ""))
            out += body
        out.append("")

    while out and not out[-1]:
        out.pop()
    return "\n".join(out) + "\n"


def api_subroutine_names(paths):
    """Names of all subroutines declared in the headers alongside the given files."""
    names = set()
    for directory in {Path(p).resolve().parent for p in paths}:
        for header in directory.glob("*.h"):
            for line in header.read_text().splitlines():
                match = SUBROUTINE_RE.match(line)
                if match:
                    names.add(match.group(1).upper())
    return names


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0 if argv else 1
    api_names = api_subroutine_names(argv)
    status = 0
    for i, path in enumerate(argv):
        try:
            sys.stdout.write(("\n" if i else "") + to_markdown(path, api_names))
        except DocblockError as error:
            print(f"{path}: error: {error}", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
