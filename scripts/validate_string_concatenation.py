#!/usr/bin/env python
"""
Check where there is a string that needs to be concatenated.
"""

import os
import sys
import token
import tokenize

FILE_EXTENTIONS_TO_CHECK = [".py", ".pyx"]


def main():
    path = sys.argv[1]

    if not os.path.exists(path):
        raise ValueError("Please enter a valid path, to a file/directory.")

    if os.path.isfile(path):
        # Means that the given path is of a single file.
        sys.exit(is_concatenated(path))

    status_codes = set()
    # Means that the given path is of a directory.
    for subdir, _, files in os.walk(path):
        for file_name in files:
            ext = os.path.splitext(os.path.join(subdir, file_name))[1]
            if ext in FILE_EXTENTIONS_TO_CHECK:
                status_codes.add(is_concatenated(os.path.join(subdir, file_name)))

    if 1 in status_codes:
        sys.exit(1)

    sys.exit(0)


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
        toks = list(tokenize.generate_tokens(file_name.readline))
        for i in range(len(toks) - 1):
            tok = toks[i]
            tok2 = toks[i + 1]
            if tok[0] == token.STRING and tok[0] == tok2[0]:
                need_fix = True
                print(
                    "{file_path}:{line_number}:\t{start} and {end}".format(
                        file_path=file_path,
                        line_number=tok[2][0],
                        start=tok[1],
                        end=tok2[1],
                    )
                )

    return int(need_fix)


if __name__ == "__main__":
    main()
