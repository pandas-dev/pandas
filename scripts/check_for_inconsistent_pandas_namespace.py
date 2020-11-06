"""
Check that test suite file doesn't contain both of e.g. ``Series`` and ``pd.Series``.

This is meant to be run as a pre-commit hook - run it manually, you can do:

    pre-commit run inconsistent-namespace-usage --all-files
"""

import argparse
import re

PATTERN = r"(?<!pd\.)(?<!\w){class_name}\(.*pd\.{class_name}\("
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()

    for class_name in CLASS_NAMES:
        pattern = re.compile(
            PATTERN.format(class_name=class_name).encode(), re.MULTILINE | re.DOTALL
        )
        for path in args.paths:
            with open(path, "rb") as f:
                contents = f.read()
            match = pattern.search(contents)
            assert (
                match is None
            ), f"Found both `pd.{class_name}` and `{class_name}` in {path}"
