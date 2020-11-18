"""
Check that test suite file doesn't use the pandas namespace inconsistently.

We check for cases of ``Series`` and ``pd.Series`` appearing in the same file
(likewise for some other common classes).

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run inconsistent-namespace-usage --all-files
"""

import argparse
from pathlib import Path
import re
from typing import Optional, Sequence

PATTERN = r"""
    (
        (?<!pd\.)(?<!\w)    # check class_name doesn't start with pd. or character
        ([A-Z]\w+)\(        # match DataFrame but not pd.DataFrame or tm.makeDataFrame
        .*                  # match anything
        pd\.\2\(            # only match e.g. pd.DataFrame
    )|
    (
        pd\.([A-Z]\w+)\(    # only match e.g. pd.DataFrame
        .*                  # match anything
        (?<!pd\.)(?<!\w)    # check class_name doesn't start with pd. or character
        \4\(                # match DataFrame but not pd.DataFrame or tm.makeDataFrame
    )
    """
ERROR_MESSAGE = "Found both `pd.{class_name}` and `{class_name}` in {path}"


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args(argv)

    pattern = re.compile(
        PATTERN.encode(),
        flags=re.MULTILINE | re.DOTALL | re.VERBOSE,
    )
    for path in args.paths:
        contents = path.read_bytes()
        match = pattern.search(contents)
        if match is None:
            continue
        if match.group(2) is not None:
            raise AssertionError(
                ERROR_MESSAGE.format(class_name=match.group(2).decode(), path=str(path))
            )
        if match.group(4) is not None:
            raise AssertionError(
                ERROR_MESSAGE.format(class_name=match.group(4).decode(), path=str(path))
            )


if __name__ == "__main__":
    main()
