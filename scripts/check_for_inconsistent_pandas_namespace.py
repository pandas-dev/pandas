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
        (?<!pd\.)(?<!\w)    # check class_name start with pd. or character
        {class_name}\(      # match DataFrame but not pd.DataFrame or tm.makeDataFrame
        .*                  # match anything
        pd\.{class_name}\(  # only match e.g. pd.DataFrame
    )|
    (
        pd\.{class_name}\(  # only match e.g. pd.DataFrame
        .*                  # match anything
        (?<!pd\.)(?<!\w)    # check class_name start with pd. or character
        {class_name}\(      # match DataFrame but not pd.DataFrame or tm.makeDataFrame
    )
    """
CLASS_NAMES = (
    "Series",
    "DataFrame",
    "Index",
    "MultiIndex",
    "Timestamp",
    "Timedelta",
    "TimedeltaIndex",
    "DatetimeIndex",
    "Categorical",
)
ERROR_MESSAGE = "Found both `pd.{class_name}` and `{class_name}` in {path}"


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args(argv)

    for class_name in CLASS_NAMES:
        pattern = re.compile(
            PATTERN.format(class_name=class_name).encode(),
            flags=re.MULTILINE | re.DOTALL | re.VERBOSE,
        )
        for path in args.paths:
            contents = path.read_bytes()
            match = pattern.search(contents)
            assert match is None, ERROR_MESSAGE.format(
                class_name=class_name, path=str(path)
            )


if __name__ == "__main__":
    main()
