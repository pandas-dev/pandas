#!/usr/bin/env python
"""
GH #23922

Check for the use of bare pytest raise.

For example:

>>> with pytest.raise(ValueError):
...    # Some code that raises ValueError

Instead of:

>>> with pytest.raise(ValueError, match="foo"):
...    # Some code that raises ValueError
"""

import argparse
import os
import sys
import token
import tokenize
from typing import Generator, List, Tuple

FILE_EXTENSIONS_TO_CHECK = ".py"


def main(source_path: str, output_format: str) -> bool:
    """
    Main entry point of the script.

    Parameters
    ----------
    source_path : str
        Source path representing path to a file/directory.
    output_format : str
        Output format of the script.

    Returns
    -------
    bool
        True if found any bare pytest raises.

    Raises
    ------
    ValueError
        If the `source_path` is not pointing to existing file/directory.
    """
    if not os.path.exists(source_path):
        raise ValueError(
            "Please enter a valid path, pointing to a valid file/directory."
        )

    is_failed: bool = False

    msg = "Bare pytests raise have been found."

    if os.path.isfile(source_path):
        for source_path, line_number in bare_pytest_raise(source_path):
            is_failed = True
            print(
                output_format.format(
                    source_path=source_path, line_number=line_number, msg=msg
                )
            )

    for subdir, _, files in os.walk(source_path):
        for file_name in files:
            if any(
                file_name.endswith(extension) for extension in FILE_EXTENSIONS_TO_CHECK
            ):
                for source_path, line_number in bare_pytest_raise(
                    os.path.join(subdir, file_name)
                ):
                    is_failed = True
                    print(
                        output_format.format(
                            source_path=source_path, line_number=line_number, msg=msg
                        )
                    )
    return is_failed


def bare_pytest_raise(source_path: str) -> Generator[Tuple[str, int], None, None]:
    """
    Yielding the files and line numbers of files with bare pytest raise.

    Parameters
    ----------
    source_path : str
        File path pointing to a single file.

    Yields
    ------
    source_path : str
        Source file path.
    line_number : int
        Line number of bare pytests raise.
    """
    with open(source_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for counter, current_token in enumerate(tokens, start=1):
        if current_token[0] == token.NAME and current_token[1] == "raises":
            for next_token in tokens[counter:]:
                if next_token[0] == token.NAME and next_token[1] == "match":
                    break
                if next_token[0] == token.NEWLINE:
                    yield source_path, current_token[2][0]
                    break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate there's no use of bare pytest raise"
    )

    parser.add_argument(
        "path", nargs="?", default=".", help="Source path of file/directory to check."
    )
    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{msg}",
        help="Output format of the error message.",
    )

    args = parser.parse_args()

    sys.exit(main(source_path=args.path, output_format=args.format))
