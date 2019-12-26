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
issue (see https://github.com/psf/black/issues/1051), so we are checking
it here.
"""

import os
import sys
import token
import tokenize

# Can be annotated as typing.FrozenSet[str]
FILE_EXTENSIONS_TO_CHECK = frozenset((".pxd", ".py", ".pyx", ".pyx.ini"))


def is_concatenated(file_path):
    """
    Checking if the file containing strings that needs to be concatenated.

    Parameters
    ----------
    file_path : str
        File path pointing to a single file.

    Returns
    -------
    int
        Status code representing if the file needs a fix.
        0 - All good.
        1 - Needs to be fixed.
    """
    need_fix = False
    with open(file_path, "r") as file_name:
        tokens = list(tokenize.generate_tokens(file_name.readline))
        for current_token, next_token in zip(tokens, tokens[1:]):
            if current_token[0] == next_token[0] == token.STRING:
                need_fix = True
                print(
                    "{file_path}:{line_number}:\t{start} and {end}".format(
                        file_path=file_path,
                        line_number=current_token[2][0],
                        start=current_token[1],
                        end=next_token[1],
                    )
                )

    return int(need_fix)


if __name__ == "__main__":
    path = sys.argv[1]

    if not os.path.exists(path):
        raise ValueError("Please enter a valid path, to a file/directory.")

    if os.path.isfile(path):
        # Means that the given path is of a single file.
        sys.exit(is_concatenated(path))

    failures = 0
    # Means that the given path is of a directory.
    for subdir, _, files in os.walk(path):
        for file_name in files:
            if any(
                file_name.endswith(extension) for extension in FILE_EXTENSIONS_TO_CHECK
            ):
                file_extension = os.path.join(subdir, file_name)
                failures += is_concatenated(os.path.join(subdir, file_name))

    exit_code = 1 if failures >= 1 else 0
    sys.exit(exit_code)
