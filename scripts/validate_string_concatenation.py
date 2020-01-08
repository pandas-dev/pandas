#!/usr/bin/env python
"""
GH #30454

Check where there is a string that needs to be concatenated.

This is necessary after black formating,
where for example black transforms this:

>>> foo = (
...         "bar "
...         "baz"
...     )

into this:

>>> foo = ("bar " "baz")

Black is not considering this as an
issue (see issue https://github.com/psf/black/issues/1051),
so we are checking it here.
"""

import argparse
import os
import sys
import token
import tokenize
from typing import Generator, List, Tuple

FILE_EXTENSIONS_TO_CHECK = (".py", ".pyx", ".pyx.ini", ".pxd")


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
        True if found any strings that needs to be concatenated.

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

    msg = "String unnecessarily split in two by black. Please merge them manually."

    if os.path.isfile(source_path):
        for source_path, line_number in strings_to_concatenate(source_path):
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
                for source_path, line_number in strings_to_concatenate(
                    os.path.join(subdir, file_name)
                ):
                    is_failed = True
                    print(
                        output_format.format(
                            source_path=source_path, line_number=line_number, msg=msg
                        )
                    )
    return is_failed


def strings_to_concatenate(source_path: str) -> Generator[Tuple[str, int], None, None]:
    """
    Yielding the strings that needs to be concatenated in a given file.

    Parameters
    ----------
    source_path : str
        File path pointing to a single file.

    Yields
    ------
    source_path : str
        Source file path.
    line_number : int
        Line number of unconcatenated string.
    """
    with open(source_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for current_token, next_token in zip(tokens, tokens[1:]):
        if current_token[0] == next_token[0] == token.STRING:
            yield source_path, current_token[2][0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate concatenated strings")

    parser.add_argument(
        "path", nargs="?", default=".", help="Source path of file/directory to check."
    )
    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{msg}",
        help="Output format of the unconcatenated strings.",
    )

    args = parser.parse_args()

    sys.exit(main(source_path=args.path, output_format=args.format))
