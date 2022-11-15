"""
Python3.9 introduces removesuffix and remove prefix.

They're reimplemented here for use in Python3.8.
"""
from __future__ import annotations


def removesuffix(string: str, suffix: str) -> str:
    if string.endswith(suffix):
        return string[: -len(suffix)]
    return string


def removeprefix(string: str, prefix: str) -> str:
    if string.startswith(prefix):
        return string[len(prefix) :]
    return string
