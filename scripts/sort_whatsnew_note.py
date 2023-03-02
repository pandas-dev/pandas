"""
Sort whatsnew note blocks by issue number.

NOTE: this assumes that each entry is on its own line, and ends with an issue number.
If that's not the case, then an entry might not get sorted. However, virtually all
recent-enough whatsnew entries follow this pattern. So, although not perfect, this
script should be good enough to significantly reduce merge conflicts.

For example:

- Fixed bug in resample (:issue:`321`)
- Fixed bug in groupby (:issue:`123`)

would become

- Fixed bug in groupby (:issue:`123`)
- Fixed bug in resample (:issue:`321`)

The motivation is to reduce merge conflicts by reducing the chances that multiple
contributors will edit the same line of code.

You can run this manually with

    pre-commit run sort-whatsnew-items --all-files
"""
from __future__ import annotations

import argparse
import re
import sys
from typing import Sequence

pattern = re.compile(r"\(:issue:`(\d+)`\)\n$")


def sort_whatsnew_note(content: str) -> int:
    new_lines = []
    block: list[str] = []
    lines = content.splitlines(keepends=True)
    for line in lines:
        if line.startswith("- ") and pattern.search(line) is not None:
            block.append(line)
        else:
            key = lambda x: int(pattern.search(x).group(1))
            block = sorted(block, key=key)
            new_lines.extend(block)
            new_lines.append(line)
            block = []
    if sorted(new_lines) != sorted(lines):  # pragma: no cover
        # Defensive check - this script should only reorder lines, not modify any
        # content.
        raise AssertionError(
            "Script modified content of file. Something is wrong, please don't "
            "trust it."
        )
    return "".join(new_lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)
    ret = 0
    for path in args.paths:
        with open(path) as fd:
            content = fd.read()
        new_content = sort_whatsnew_note(content)
        if content != new_content:
            ret |= 1
            with open(path, "w") as fd:
                fd.write(new_content)
    return ret


if __name__ == "__main__":
    sys.exit(main())
