#!/usr/bin/env python3
"""
Unwanted patterns test cases.

The reason this file exist despite the fact we already have
`ci/code_checks.sh`,
(see https://github.com/pandas-dev/pandas/blob/master/ci/code_checks.sh)

is that some of the test cases are more complex/imposible to validate via regex.
So this file is somewhat an extensions to `ci/code_checks.sh`
"""

import argparse
import ast
import os
import sys
import token
import tokenize
from typing import IO, Callable, FrozenSet, Iterable, List, Tuple

PATHS_TO_IGNORE: Tuple[str, ...] = ("asv_bench/env",)


def _get_literal_string_prefix_len(token_string: str) -> int:
    """
    Getting the length of the literal string prefix.

    Parameters
    ----------
    token_string : str
        String to check.

    Returns
    -------
    int
        Length of the literal string prefix.

    Examples
    --------
    >>> example_string = "'Hello world'"
    >>> _get_literal_string_prefix_len(example_string)
    0
    >>> example_string = "r'Hello world'"
    >>> _get_literal_string_prefix_len(example_string)
    1
    """
    try:
        return min(
            token_string.find(quote)
            for quote in (r"'", r'"')
            if token_string.find(quote) >= 0
        )
    except ValueError:
        return 0


def bare_pytest_raises(file_obj: IO[str]) -> Iterable[Tuple[int, str]]:
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
    file_obj : IO
        File-like object containing the Python code to validate.

    Yields
    ------
    line_number : int
        Line number of unconcatenated string.
    msg : str
        Explenation of the error.

    Notes
    -----
    GH #23922
    """
    contents = file_obj.read()
    tree = ast.parse(contents)

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        try:
            if not (node.func.value.id == "pytest" and node.func.attr == "raises"):
                continue
        except AttributeError:
            continue

        if not node.keywords:
            yield (
                node.lineno,
                "Bare pytests raise have been found. "
                "Please pass in the argument 'match' as well the exception.",
            )
        else:
            # Means that there are arguments that are being passed in,
            # now we validate that `match` is one of the passed in arguments
            if not any(keyword.arg == "match" for keyword in node.keywords):
                yield (
                    node.lineno,
                    "Bare pytests raise have been found. "
                    "Please pass in the argument 'match' as well the exception.",
                )


def strings_to_concatenate(file_obj: IO[str]) -> Iterable[Tuple[int, str]]:
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
    file_obj : IO
        File-like object containing the Python code to validate.

    Yields
    ------
    line_number : int
        Line number of unconcatenated string.
    msg : str
        Explenation of the error.

    Notes
    -----
    GH #30454
    """
    tokens: List = list(tokenize.generate_tokens(file_obj.readline))

    for current_token, next_token in zip(tokens, tokens[1:]):
        if current_token.type == next_token.type == token.STRING:
            yield (
                current_token.start[0],
                (
                    "String unnecessarily split in two by black. "
                    "Please merge them manually."
                ),
            )


def strings_with_wrong_placed_whitespace(
    file_obj: IO[str],
) -> Iterable[Tuple[int, str]]:
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
    file_obj : IO
        File-like object containing the Python code to validate.

    Yields
    ------
    line_number : int
        Line number of unconcatenated string.
    msg : str
        Explenation of the error.
    """

    def has_wrong_whitespace(first_line: str, second_line: str) -> bool:
        """
        Checking if the two lines are mattching the unwanted pattern.

        Parameters
        ----------
        first_line : str
            First line to check.
        second_line : str
            Second line to check.

        Returns
        -------
        bool
            True if the two recived string match, an unwanted pattern.

        Notes
        -----
        The unwanted pattern that we are trying to catch is if the spaces in
        a string that is concatenated over multiple lines are placed at the
        end of each string, unless this string is ending with a
        newline character (\n).

        For example, this is bad:

        >>> rule = (
        ...    "We want the space at the end of the line,"
        ...    " not at the beginning"
        ... )

        And what we want is:

        >>> rule = (
        ...    "We want the space at the end of the line, "
        ...    "not at the beginning"
        ... )

        And if the string is ending with a new line character (\n) we
        do not want any trailing whitespaces after it.

        For example, this is bad:

        >>> rule = (
        ...    "We want the space at the begging of "
        ...    "the line if the previous line is ending with a \n "
        ...    "not at the end, like always"
        ... )

        And what we do want is:

        >>> rule = (
        ...    "We want the space at the begging of "
        ...    "the line if the previous line is ending with a \n"
        ...    " not at the end, like always"
        ... )
        """
        if first_line.endswith(r"\n"):
            return False
        elif first_line.startswith("  ") or second_line.startswith("  "):
            return False
        elif first_line.endswith("  ") or second_line.endswith("  "):
            return False
        elif (not first_line.endswith(" ")) and second_line.startswith(" "):
            return True
        return False

    tokens: List = list(tokenize.generate_tokens(file_obj.readline))

    for first_token, second_token, third_token in zip(tokens, tokens[1:], tokens[2:]):
        # Checking if we are in a block of concated string
        if (
            first_token.type == third_token.type == token.STRING
            and second_token.type == token.NL
        ):
            # Striping the quotes, with the string litteral prefix
            first_string: str = first_token.string[
                _get_literal_string_prefix_len(first_token.string) + 1 : -1
            ]
            second_string: str = third_token.string[
                _get_literal_string_prefix_len(third_token.string) + 1 : -1
            ]

            if has_wrong_whitespace(first_string, second_string):
                yield (
                    third_token.start[0],
                    (
                        "String has a space at the beginning instead "
                        "of the end of the previous string."
                    ),
                )


def main(
    function: Callable[[IO[str]], Iterable[Tuple[int, str]]],
    source_path: str,
    output_format: str,
    file_extensions_to_check: str,
) -> bool:
    """
    Main entry point of the script.

    Parameters
    ----------
    function : Callable
        Function to execute for the specified validation type.
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
        raise ValueError("Please enter a valid path, pointing to a file/directory.")

    is_failed: bool = False
    file_path: str = ""

    FILE_EXTENSIONS_TO_CHECK: FrozenSet[str] = frozenset(
        file_extensions_to_check.split(",")
    )

    if os.path.isfile(source_path):
        file_path = source_path
        with open(file_path, "r") as file_obj:
            for line_number, msg in function(file_obj):
                is_failed = True
                print(
                    output_format.format(
                        source_path=file_path, line_number=line_number, msg=msg
                    )
                )

    for subdir, _, files in os.walk(source_path):
        if any(path in subdir for path in PATHS_TO_IGNORE):
            continue
        for file_name in files:
            if not any(
                file_name.endswith(extension) for extension in FILE_EXTENSIONS_TO_CHECK
            ):
                continue

            file_path = os.path.join(subdir, file_name)
            with open(file_path, "r") as file_obj:
                for line_number, msg in function(file_obj):
                    is_failed = True
                    print(
                        output_format.format(
                            source_path=file_path, line_number=line_number, msg=msg
                        )
                    )

    return is_failed


if __name__ == "__main__":
    available_validation_types: List[str] = [
        "bare_pytest_raises",
        "strings_to_concatenate",
        "strings_with_wrong_placed_whitespace",
    ]

    parser = argparse.ArgumentParser(description="Unwanted patterns checker.")

    parser.add_argument(
        "path", nargs="?", default=".", help="Source path of file/directory to check."
    )
    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{msg}",
        help="Output format of the error message.",
    )
    parser.add_argument(
        "--validation-type",
        "-vt",
        choices=available_validation_types,
        required=True,
        help="Validation test case to check.",
    )
    parser.add_argument(
        "--included-file-extensions",
        default="py,pyx,pxd,pxi",
        help="Coma seperated file extensions to check.",
    )

    args = parser.parse_args()

    sys.exit(
        main(
            function=globals().get(args.validation_type),  # type: ignore
            source_path=args.path,
            output_format=args.format,
            file_extensions_to_check=args.included_file_extensions,
        )
    )
