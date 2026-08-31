"""
Check that newly added GitHub references use the full URL form.

pandas has accumulated two dozen spellings of the same thing -- ``GH#1234``,
``GH 1234``, ``GH-1234``, ``GH1234``, ``gh-1234``, ``GH: #1234`` -- and
GH-55461 settled on the full URL::

    # https://github.com/pandas-dev/pandas/issues/1234

Rewriting the ~12k references already in the tree would touch a comment on
every third test and reflow a couple of thousand of them to stay under the
line limit, so they are grandfathered in ``scripts/gh_reference_baseline.txt``
instead.  This hook only rejects a short reference that is *not* in that
baseline, which in practice means a newly added one.

If you move or rename a file, its grandfathered references move with it::

    python scripts/validate_gh_references.py --update-baseline

That command refuses to run when it would grow the baseline: it is there to
relocate existing references, not to silence new ones.  A rebase that picks up
short references from ``main`` is the one case where growth is legitimate, and
``--allow-growth`` overrides the refusal -- the added lines show up in the diff,
so a reviewer can tell the two apart.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run gh-reference-format --all-files
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import (
        Iterator,
        Sequence,
    )

BASELINE_PATH = Path(__file__).parent / "gh_reference_baseline.txt"

# --update-baseline walks the tracked tree, minus the same directory the
# pre-commit hook excludes: this checker and its tests have to spell out the
# forms they reject.
EXCLUDED_ROOT = "scripts"
SUFFIXES = frozenset({".py", ".pyi", ".pyx", ".pxd", ".pxi"})

# The separator class stays narrow on purpose.  ``_`` would match identifiers
# such as ``gh_13141_expected``, and ``,`` would match CSV fixture data.  The
# trailing lookahead keeps us off both of those as well as off longer numbers.
SHORT_REF = re.compile(r"(?<![A-Za-z0-9_])[Gg][Hh][ :#-]{0,4}([0-9]{2,6})(?![0-9_])")

BASELINE_HEADER = """\
# Short-form GitHub references that predate the full-URL convention (GH-55461).
#
# One line per file: the path, then every issue number referenced in it.
# Regenerate with:
#
#     python scripts/validate_gh_references.py --update-baseline
#
# It grows only when a rebase carries short references over from main
# (--allow-growth).  New references use the full URL form:
#
#     # https://github.com/pandas-dev/pandas/issues/1234
#
"""


def find_short_refs(content: str) -> Iterator[tuple[int, int, str, str]]:
    """
    Yield ``(line number, column, matched text, issue number)`` for each
    short-form reference in ``content``.
    """
    for lineno, line in enumerate(content.splitlines(), start=1):
        for match in SHORT_REF.finditer(line):
            yield lineno, match.start() + 1, match.group(0), match.group(1)


def main(content: str, file: str, baseline: dict[str, set[str]]) -> int:
    allowed = baseline.get(Path(file).as_posix(), set())
    ret = 0
    for lineno, col, text, number in find_short_refs(content):
        if number in allowed:
            continue
        print(
            f"{file}:{lineno}:{col} found {text!r}; reference issues as "
            f"https://github.com/pandas-dev/pandas/issues/{number}"
        )
        ret = 1
    return ret


def load_baseline(path: Path = BASELINE_PATH) -> dict[str, set[str]]:
    baseline: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        filename, *numbers = line.split()
        baseline[filename] = set(numbers)
    return baseline


def collect_baseline(filenames: Sequence[str] | None = None) -> dict[str, set[str]]:
    if filenames is None:
        tracked = subprocess.run(
            ["git", "ls-files"], capture_output=True, text=True, check=True
        )
        filenames = tracked.stdout.splitlines()

    found: dict[str, set[str]] = {}
    for filename in sorted(filenames):
        path = Path(filename)
        if path.suffix not in SUFFIXES or path.parts[0] == EXCLUDED_ROOT:
            continue
        content = path.read_text(encoding="utf-8")
        numbers = {number for _, _, _, number in find_short_refs(content)}
        if numbers:
            found[path.as_posix()] = numbers
    return found


def update_baseline(path: Path = BASELINE_PATH, allow_growth: bool = False) -> int:
    # ``None`` only on the initial bootstrap, when there is nothing to grow from
    before = (
        sum(len(numbers) for numbers in load_baseline(path).values())
        if path.exists()
        else None
    )
    found = collect_baseline()
    after = sum(len(numbers) for numbers in found.values())

    if before is not None and after > before and not allow_growth:
        print(
            f"refusing to update: the baseline would grow from {before} to {after} "
            "references. A new reference must use the full URL form, e.g. "
            "https://github.com/pandas-dev/pandas/issues/1234. If the growth came "
            "from rebasing onto main, pass --allow-growth"
        )
        return 1

    lines = [
        f"{filename} {' '.join(sorted(numbers, key=int))}"
        for filename, numbers in sorted(found.items())
    ]
    path.write_text(BASELINE_HEADER + "\n".join(lines) + "\n", encoding="utf-8")
    if before is None:
        print(f"baseline: created with {after} references")
    else:
        print(f"baseline: {before} -> {after} references")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    parser.add_argument(
        "--update-baseline",
        action="store_true",
        help="rewrite the baseline, e.g. after moving a file",
    )
    parser.add_argument(
        "--allow-growth",
        action="store_true",
        help="let --update-baseline absorb references picked up from main",
    )
    args = parser.parse_args()

    if args.update_baseline:
        sys.exit(update_baseline(allow_growth=args.allow_growth))

    loaded = load_baseline()
    ret = 0
    for file in args.paths:
        with open(file, encoding="utf-8") as fd:
            ret |= main(fd.read(), file, loaded)

    sys.exit(ret)
