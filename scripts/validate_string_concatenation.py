#!/usr/bin/env python3
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

import os
import sys
import token
import tokenize
from typing import FrozenSet, Generator, List

FILE_EXTENSIONS_TO_CHECK: FrozenSet[str] = frozenset(
    (".pxd", ".py", ".pyx", ".pyx.ini")
)


def strings_to_concatenate(file_path: str) -> Generator[str, None, None]:
    """
    Yielding the strings that needs to be concatenated in a given file.

    Parameters
    ----------
    file_path : str
        File path pointing to a single file.

    Yields
    ------
    str
        Message containing info about the string that needs to be concatenated.
    """
    with open(file_path, "r") as file_name:
        tokens: List = list(tokenize.generate_tokens(file_name.readline))

    for current_token, next_token in zip(tokens, tokens[1:]):
        if current_token[0] == next_token[0] == token.STRING:
            line_number = current_token[2][0]
            start = current_token[1]
            end = next_token[1]
            yield f"{file_path}:{line_number}:\t between {start} and {end}\n"


if __name__ == "__main__":
    path: str = sys.argv[1]

    if not os.path.exists(path):
        raise ValueError("Please enter a valid path, to a file/directory.")

    failed: bool = False

    if os.path.isfile(path):
        for msg in strings_to_concatenate(path):
            if msg:
                failed = True
                print(msg)

    for subdir, _, files in os.walk(path):
        for file_name in files:
            if any(
                file_name.endswith(extension) for extension in FILE_EXTENSIONS_TO_CHECK
            ):
                file_extension = os.path.join(subdir, file_name)

                for msg in strings_to_concatenate(os.path.join(subdir, file_name)):
                    if msg:
                        failed = True
                        print(msg)
    sys.exit(failed)
