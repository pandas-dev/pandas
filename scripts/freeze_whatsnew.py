"""
Freeze ipython directives in whatsnew RST files.

Converts ``.. ipython::`` directives to static ``.. code-block:: pycon``
blocks with captured ``>>>`` REPL output, so that old whatsnew notes no
longer execute during doc builds. See GH#6856.

Usage
-----
    python scripts/freeze_whatsnew.py doc/source/whatsnew/v0.18.0.rst
    python scripts/freeze_whatsnew.py doc/source/whatsnew/v1.*.rst
    python scripts/freeze_whatsnew.py --dry-run doc/source/whatsnew/v2.0.0.rst

To produce output with an older pandas version, run the script with
the Python whose environment has that version installed::

    /path/to/old-env/bin/python scripts/freeze_whatsnew.py ...
"""

from __future__ import annotations

import ast
from io import StringIO
from pathlib import Path
import re
import sys
import warnings

# ---------------------------------------------------------------------------
# RST parsing
# ---------------------------------------------------------------------------

_DIRECTIVE_RE = re.compile(r"^(\s*)\.\. ipython::(.*)$")
_OPTION_RE = re.compile(r"^\s+:(\w+):\s*$")


def _parse_blocks(lines: list[str]) -> list[dict]:
    """Return a list of ipython directive blocks found in *lines*."""
    blocks: list[dict] = []
    idx = 0

    while idx < len(lines):
        match = _DIRECTIVE_RE.match(lines[idx])
        if not match:
            idx += 1
            continue

        indent = match.group(1)
        indent_len = len(indent)
        start = idx
        idx += 1

        # Collect options (e.g. :suppress:, :okwarning:, :verbatim:)
        options: set[str] = set()
        while idx < len(lines):
            opt = _OPTION_RE.match(lines[idx])
            if opt:
                options.add(opt.group(1))
                idx += 1
            elif lines[idx].strip() == "":
                idx += 1
                break
            else:
                break

        # Skip additional blank lines before content
        while idx < len(lines) and lines[idx].strip() == "":
            idx += 1

        # Determine content indentation from first content line
        content_indent = ""
        if idx < len(lines) and lines[idx].strip():
            content_indent = re.match(r"^(\s*)", lines[idx]).group(1)

        # Collect body lines — anything indented deeper than the directive
        body_lines: list[str] = []
        while idx < len(lines):
            if lines[idx].strip() == "":
                body_lines.append("")
                idx += 1
            else:
                leading = len(lines[idx]) - len(lines[idx].lstrip())
                if leading > indent_len:
                    body_lines.append(lines[idx])
                    idx += 1
                else:
                    break

        # Strip trailing blank lines
        while body_lines and body_lines[-1] == "":
            body_lines.pop()

        # Dedent the code
        ci_len = len(content_indent)
        code_lines: list[str] = []
        for bline in body_lines:
            if bline == "":
                code_lines.append("")
            elif bline.startswith(content_indent):
                code_lines.append(bline[ci_len:])
            else:
                code_lines.append(bline.lstrip())

        blocks.append(
            {
                "start": start,
                "end": idx,
                "indent": indent,
                "content_indent": content_indent,
                "options": options,
                "code_lines": code_lines,
                "code_text": "\n".join(code_lines),
            }
        )

    return blocks


# ---------------------------------------------------------------------------
# Code splitting — break a block's code into individual cells
# ---------------------------------------------------------------------------


def _split_cells(code_text: str) -> list[tuple[str, bool]]:
    """
    Split *code_text* into individual executable cells.

    The IPython directive executes each top-level statement individually.
    Blank lines separate "groups", and within a group each top-level
    statement is its own cell.  Multi-line constructs (function defs, loops,
    ``with`` blocks, try/except, multi-line expressions via parens/brackets)
    are kept together.

    Returns a list of ``(code, suppressed)`` tuples.  The ``@suppress``
    inline directive causes the *next* statement to be marked as suppressed
    (executed but not shown).
    """
    lines = code_text.split("\n")
    cells: list[tuple[str, bool]] = []
    suppress_next = False

    # First pass: strip @suppress markers and tag following lines
    tagged_lines: list[tuple[str, bool]] = []
    for line in lines:
        if line.strip() == "@suppress":
            suppress_next = True
            continue
        tagged_lines.append((line, suppress_next))
        if suppress_next and line.strip():
            suppress_next = False

    # Second pass: group into statements
    # Blank lines are group separators. Within a group, each top-level
    # statement is a cell. A line that starts with whitespace is a
    # continuation of the previous statement.
    group: list[tuple[str, bool]] = []

    def _flush_group():
        if not group:
            return
        # Split the group into individual statements.
        # A new statement starts when: a line has no leading whitespace AND
        # there are no unclosed brackets or triple-quoted strings from
        # previous lines.
        current_lines: list[str] = []
        current_suppressed = False
        bracket_depth = 0
        in_triple_quote: str | None = None  # None, '"""', or "'''"

        for gline, gsup in group:
            stripped = gline.strip()
            if not stripped:
                continue
            is_indented = gline[0] == " " or gline[0] == "\t"
            # Start a new statement if we're at column 0, brackets are
            # balanced, and we're not inside a triple-quoted string.
            if (
                current_lines
                and not is_indented
                and bracket_depth == 0
                and in_triple_quote is None
            ):
                cells.append(("\n".join(current_lines), current_suppressed))
                current_lines = []
                current_suppressed = False
                bracket_depth = 0
            if not current_lines:
                current_suppressed = gsup
            current_lines.append(gline)
            # Track bracket depth and triple-quote state
            idx = 0
            while idx < len(stripped):
                if in_triple_quote is not None:
                    # Look for the closing triple-quote
                    close_pos = stripped.find(in_triple_quote, idx)
                    if close_pos == -1:
                        break  # rest of line is inside the string
                    in_triple_quote = None
                    idx = close_pos + 3
                elif stripped[idx : idx + 3] in ('"""', "'''"):
                    in_triple_quote = stripped[idx : idx + 3]
                    # Check if closed on the same line (after the opener)
                    close_pos = stripped.find(in_triple_quote, idx + 3)
                    if close_pos != -1:
                        in_triple_quote = None
                        idx = close_pos + 3
                    else:
                        break  # rest of line is inside the string
                elif stripped[idx] in "([{":
                    bracket_depth += 1
                    idx += 1
                elif stripped[idx] in ")]}":
                    bracket_depth = max(0, bracket_depth - 1)
                    idx += 1
                elif stripped[idx] in ('"', "'"):
                    # Regular string — skip to closing quote
                    quote = stripped[idx]
                    idx += 1
                    while idx < len(stripped):
                        if stripped[idx] == "\\":
                            idx += 2
                        elif stripped[idx] == quote:
                            idx += 1
                            break
                        else:
                            idx += 1
                elif stripped[idx] == "#":
                    break  # rest of line is a comment
                else:
                    idx += 1

        if current_lines:
            cells.append(("\n".join(current_lines), current_suppressed))

    for line, sup in tagged_lines:
        if line == "":
            _flush_group()
            group = []
        else:
            group.append((line, sup))

    _flush_group()
    return cells


# ---------------------------------------------------------------------------
# Code execution
# ---------------------------------------------------------------------------


class _Executor:
    """Execute code cells and format output with ``>>>`` REPL prompts."""

    def __init__(self) -> None:
        self.namespace: dict = {"__name__": "__main__", "__builtins__": __builtins__}
        self._silent_exec("import numpy as np")
        self._silent_exec("import pandas as pd")

    # -- helpers -------------------------------------------------------------

    def _silent_exec(self, code: str) -> None:
        exec(compile(code, "<setup>", "exec"), self.namespace)  # noqa: S102

    # -- public API ----------------------------------------------------------

    def run_silent(self, code: str) -> None:
        """Execute *code* without producing any output."""
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                exec(compile(code, "<suppress>", "exec"), self.namespace)  # noqa: S102
        except Exception:
            pass

    def run_cell(self, code: str) -> tuple[str, bool]:
        """
        Execute *code* and return ``(formatted_text, success)``.

        *formatted_text* uses standard ``>>>`` / ``...`` REPL prompts.
        """
        code = code.rstrip()
        cell_lines = code.split("\n")
        if not cell_lines or (len(cell_lines) == 1 and not cell_lines[0].strip()):
            return "", True

        # Format the input prompts
        parts = [f">>> {cell_lines[0]}"]
        parts.extend(f"... {continuation}" for continuation in cell_lines[1:])
        input_text = "\n".join(parts)

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = captured = StringIO()
        sys.stderr = StringIO()

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tree = ast.parse(code)

                if not tree.body:
                    sys.stdout, sys.stderr = old_stdout, old_stderr
                    return input_text, True

                last_is_expr = isinstance(tree.body[-1], ast.Expr)

                if last_is_expr and len(tree.body) > 1:
                    # Execute all-but-last as statements, eval last
                    mod = ast.Module(body=tree.body[:-1], type_ignores=[])
                    exec(compile(mod, "<cell>", "exec"), self.namespace)  # noqa: S102
                    result = eval(
                        compile(ast.Expression(tree.body[-1].value), "<cell>", "eval"),
                        self.namespace,
                    )
                elif last_is_expr:
                    result = eval(
                        compile(ast.Expression(tree.body[0].value), "<cell>", "eval"),
                        self.namespace,
                    )
                else:
                    exec(compile(tree, "<cell>", "exec"), self.namespace)  # noqa: S102
                    result = None

                printed = captured.getvalue()
                sys.stdout, sys.stderr = old_stdout, old_stderr

                output_parts = [input_text]
                if printed:
                    output_parts.append(printed.rstrip())
                if last_is_expr and result is not None:
                    output_parts.append(repr(result))
                return "\n".join(output_parts), True

        except SyntaxError:
            sys.stdout, sys.stderr = old_stdout, old_stderr
            return input_text, False
        except Exception:
            sys.stdout, sys.stderr = old_stdout, old_stderr
            return input_text, False


# ---------------------------------------------------------------------------
# File-level processing
# ---------------------------------------------------------------------------


def freeze_file(
    filepath: str | Path, *, dry_run: bool = False
) -> tuple[int, list[str]]:
    """
    Freeze all ``.. ipython::`` directives in *filepath*.

    Returns ``(blocks_processed, list_of_warnings)``.
    """
    filepath = Path(filepath)
    content = filepath.read_text()
    lines = content.split("\n")
    blocks = _parse_blocks(lines)

    if not blocks:
        return 0, []

    executor = _Executor()
    results: list[dict] = []

    # Forward pass — execute every block in order, collect replacement info
    for block in blocks:
        options = block["options"]
        code_text = block["code_text"]

        if "suppress" in options:
            executor.run_silent(code_text)
            results.append({"action": "remove"})
            continue

        if "verbatim" in options:
            results.append({"action": "verbatim"})
            continue

        cells = _split_cells(code_text)
        output_parts: list[str] = []
        all_ok = True

        for cell_code, suppressed in cells:
            if suppressed:
                executor.run_silent(cell_code)
                continue
            formatted, success = executor.run_cell(cell_code)
            if not success:
                all_ok = False
            output_parts.append(formatted)

        results.append(
            {"action": "replace", "output_parts": output_parts, "all_ok": all_ok}
        )

    # Reverse pass — replace blocks back-to-front so line numbers stay valid
    new_lines = list(lines)
    warn_list: list[str] = []
    processed = 0

    for block, result in reversed(list(zip(blocks, results, strict=True))):
        indent = block["indent"]
        ci = block["content_indent"]

        if result["action"] == "remove":
            end = block["end"]
            # Also consume trailing blank lines that belonged to the block
            while end < len(new_lines) and new_lines[end].strip() == "":
                end += 1
            new_lines[block["start"] : end] = []
            processed += 1
            continue

        if result["action"] == "verbatim":
            replacement = [f"{indent}.. code-block:: pycon", ""]
            for code_line in block["code_lines"]:
                if code_line == "":
                    replacement.append("")
                else:
                    replacement.append(f"{ci}{code_line}")
            replacement.append("")
            new_lines[block["start"] : block["end"]] = replacement
            processed += 1
            continue

        # action == "replace"
        output_parts = result["output_parts"]
        if not output_parts:
            continue

        replacement = [f"{indent}.. code-block:: pycon", ""]
        for part in output_parts:
            for output_line in part.split("\n"):
                replacement.append(f"{ci}{output_line}")
            replacement.append("")

        # Strip trailing blank lines, then add exactly one back so
        # the code-block is separated from the following RST content.
        while replacement and replacement[-1].strip() == "":
            replacement.pop()
        replacement.append("")

        new_lines[block["start"] : block["end"]] = replacement
        processed += 1

        if not result["all_ok"]:
            warn_list.append(
                f"  {filepath}:{block['start'] + 1}: some cells failed to execute"
            )

    if not dry_run:
        filepath.write_text("\n".join(new_lines))

    return processed, warn_list


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Freeze ipython directives in whatsnew RST files (GH#6856)."
    )
    parser.add_argument("files", nargs="+", help="RST files to process")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what would be done without writing files.",
    )
    args = parser.parse_args()

    total = 0
    all_warnings: list[str] = []

    for filepath in args.files:
        processed, file_warnings = freeze_file(filepath, dry_run=args.dry_run)
        total += processed
        all_warnings.extend(file_warnings)
        if processed == 0:
            print(f"  {filepath}: no ipython blocks")
        else:
            tag = "OK" if not file_warnings else f"{len(file_warnings)} warnings"
            print(f"  {filepath}: {processed} blocks frozen — {tag}")

    print(f"\nTotal: {total} blocks frozen across {len(args.files)} files")
    if all_warnings:
        print(f"\nWarnings ({len(all_warnings)}):")
        for warning in all_warnings:
            print(warning)


if __name__ == "__main__":
    main()
