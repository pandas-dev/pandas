"""
Script to run ``autotyping``, to get around the fact that
pre-commit puts ``args`` before the list of files, whereas
``autotyping`` wants the files to come after, see
https://github.com/pandas-dev/pandas/issues/48808#issuecomment-1259711679.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args(argv)
    if not args.paths:
        sys.exit(0)
    output = subprocess.run(
        [
            "python",
            "-m",
            "libcst.tool",
            "codemod",
            "autotyping.AutotypeCommand",
            *args.paths,
            "--no-format",
            "--safe",
            # all except 'guess-common-names' from 'aggresive'
            "--bool-param",
            "--int-param",
            "--float-param",
            "--str-param",
            "--bytes-param",
            "--annotate-imprecise-magics",
        ],
        check=True,
    )
    sys.exit(output.returncode)


if __name__ == "__main__":
    main()
