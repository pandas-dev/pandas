#!/usr/bin/env python
"""
Test case for leading spaces in concated strings.

For example:

# Good
foo = (
        "bar "
        "baz"
)


# Bad
foo = (
        "bar"
        " baz"
)
"""

import argparse
import os
import sys
import token
import tokenize
from typing import Generator, List, Tuple

FILE_EXTENSIONS_TO_CHECK = (".py", ".pyx", ".pyx.ini", ".pxd")

MSG = "String has a space at the beginning instead of the end of the previous string."


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
        True if found any strings that have a leading space in the wrong line.

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

    if os.path.isfile(source_path):
        for source_path, line_number in strings_with_wrong_space(source_path):
            is_failed = True
            print(
                output_format.format(
                    source_path=source_path, line_number=line_number, msg=MSG
                )
            )

    for subdir, _, files in os.walk(source_path):
        for file_name in files:
            if any(
                file_name.endswith(extension) for extension in FILE_EXTENSIONS_TO_CHECK
            ):
                for source_path, line_number in strings_with_wrong_space(
                    os.path.join(subdir, file_name)
                ):
                    is_failed = True
                    print(
                        output_format.format(
                            source_path=source_path, line_number=line_number, msg=MSG
                        )
                    )
    return is_failed


def strings_with_wrong_space(
    source_path: str,
) -> Generator[Tuple[str, int], None, None]:
    """
    Yielding the file path and the line number of the string to fix.

    Parameters
    ----------
    source_path : str
        File path pointing to a single file.

    Yields
    ------
    source_path : str
        Source file path.
    line_number : int
        Line number of the wrong placed space.
    """
    with open(source_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for first_token, second_token, third_token in zip(tokens, tokens[1:], tokens[2:]):
        if (
            first_token[0] == third_token[0] == token.STRING
            and second_token[0] == token.NL
        ):
            # Means we are in a block of concated string

            # Striping the quotes
            first_string = first_token[1][1:-1]
            second_string = third_token[1][1:-1]

            if (not first_string.endswith(" ")) and (second_string.startswith(" ")):
                yield source_path, third_token[2][0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate spaces over concated strings"
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
