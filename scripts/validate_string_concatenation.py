#!/usr/bin/env python
"""
Unwanted patterns test cases.

The reason this file exist despite the fact we already have
`ci/code_checks.sh`,
(see https://github.com/pandas-dev/pandas/blob/master/ci/code_checks.sh)

is that some of the test cases are more complex/imposible to validate via regex.
So this file is somewhat an extensions to `ci/code_checks.sh`
"""

import argparse
import os
import sys
import token
import tokenize
from typing import Callable, Generator, List, Tuple

FILE_EXTENSIONS_TO_CHECK = (".py", ".pyx", ".pyx.ini", ".pxd")
TYPE: int = 0
VALUE: int = 1


def strings_to_concatenate(
    source_path: str,
) -> Generator[Tuple[str, int, str], None, None]:
    """
    This test case is necessary after 'Black' (https://github.com/psf/black),
    is formating strings over multiple lines.

    For example, when this:

    >>> foo = (
    ...     "bar "
    ...     "baz"
    ... )

    Is becoming this:

    >>> foo = ("bar " "baz")

    'Black' is not considering this as an
    issue (see https://github.com/psf/black/issues/1051),
    so we are checking it here instead.

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
    msg : str
        Explenation of the error.

    Notes
    -----
    GH #30454
    """
    with open(source_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for current_token, next_token in zip(tokens, tokens[1:]):
        if current_token[TYPE] == next_token[TYPE] == token.STRING:
            yield (
                source_path,
                current_token[2][0],
                (
                    "String unnecessarily split in two by black. "
                    "Please merge them manually."
                ),
            )


def strings_with_wrong_placed_space(
    source_path: str,
) -> Generator[Tuple[str, int, str], None, None]:
    """
    Test case for leading spaces in concated strings.

    For example:

    >>> rule = (
    ...    "We want the space at the end of the line, "
    ...    "not at the beginning"
    ... )

    Instead of:

    >>> rule = (
    ...    "We want the space at the end of the line,"
    ...    " not at the beginning"
    ... )

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
    msg : str
        Explenation of the error.
    """
    with open(source_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for first_token, second_token, third_token in zip(tokens, tokens[1:], tokens[2:]):
        # Checking if we are in a block of concated string
        if (
            first_token[TYPE] == third_token[TYPE] == token.STRING
            and second_token[0] == token.NL
        ):
            # Striping the quotes
            first_string: str = first_token[1][1:-1]
            second_string: str = third_token[1][1:-1]

            is_first_ends_single: bool = first_string.endswith(" ")
            is_first_ends_newline: bool = first_string.endswith("\n")
            is_second_starts_single: bool = second_string.startswith(" ")
            is_second_starts_double: bool = second_string.startswith("  ")

            if (
                (not is_first_ends_single)
                and is_second_starts_single
                and (not is_second_starts_double)
            ) and (not is_first_ends_newline):
                yield (
                    source_path,
                    third_token[2][0],
                    (
                        "String has a space at the beginning instead "
                        "of the end of the previous string."
                    ),
                )


def bare_pytest_raises(source_path: str) -> Generator[Tuple[str, int, str], None, None]:
    """
    Test Case for bare pytest raises.

    For example, this is wrong:

    >>> with pytest.raise(ValueError):
    ...     # Some code that raises ValueError

    And this is what we want instead:

    >>> with pytest.raise(ValueError, match="foo"):
    ...     # Some code that raises ValueError

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
    msg : str
        Explenation of the error.

    Notes
    -----
    GH #23922
    """
    with open(source_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for counter, current_token in enumerate(tokens, start=1):
        if not (current_token[TYPE] == token.NAME and current_token[VALUE] == "raises"):
            continue
        for next_token in tokens[counter:]:
            if next_token[TYPE] == token.NAME and next_token[VALUE] == "match":
                break
            # token.NEWLINE refers to end of a logical line
            # unlike token.NL or "\n" which represents a newline
            if next_token[TYPE] == token.NEWLINE:
                yield (
                    source_path,
                    current_token[2][0],
                    "Bare pytests raise have been found.",
                )
                break


def main(
    function: Callable[[str], Generator[Tuple[str, int, str], None, None]],
    source_path: str,
    output_format: str,
) -> bool:
    """
    Main entry point of the script.

    Parameters
    ----------
    function : Callable
        Function to execute for the test case.
    source_path : str
        Source path representing path to a file/directory.
    output_format : str
        Output format of the error message.

    Returns
    -------
    bool
        True if found any patterns are found related to the given function.

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
        for source_path, line_number, msg in function(source_path):
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
                for source_path, line_number, msg in function(
                    os.path.join(subdir, file_name)
                ):
                    is_failed = True
                    print(
                        output_format.format(
                            source_path=source_path, line_number=line_number, msg=msg
                        )
                    )
    return is_failed


if __name__ == "__main__":
    available_validation_types: List[str] = [
        "strings_to_concatenate",
        "strings_with_wrong_placed_space",
        "bare_pytest_raises",
    ]

    parser = argparse.ArgumentParser(description="Unwanted patterns checker.")

    parser.add_argument(
        "path", nargs="?", default=".", help="Source path of file/directory to check."
    )
    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{msg}.",
        help="Output format of the error message.",
    )
    parser.add_argument(
        "--validation-type",
        "-vt",
        choices=available_validation_types,
        required=True,
        help="Validation test case to check.",
    )

    args = parser.parse_args()

    sys.exit(
        main(
            function=globals().get(args.validation_type),
            source_path=args.path,
            output_format=args.format,
        )
    )
