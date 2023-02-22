"""
Python3.9 introduces removesuffix and remove prefix.

They're reimplemented here for use in Python3.8.

NOTE: when pyupgrade --py39-plus removes nearly everything in this file,
this file and the associated tests should be removed.
"""
from __future__ import annotations

import sys

if sys.version_info < (3, 9):

    def removesuffix(string: str, suffix: str) -> str:
        if string.endswith(suffix):
            return string[: -len(suffix)]
        return string

    def removeprefix(string: str, prefix: str) -> str:
        if string.startswith(prefix):
            return string[len(prefix) :]
        return string

else:
    # NOTE: remove this file when pyupgrade --py39-plus removes
    # the above block
    pass
