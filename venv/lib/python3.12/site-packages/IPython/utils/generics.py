# encoding: utf-8
"""Generic functions for extending IPython."""

from __future__ import annotations

from typing import Any

from IPython.core.error import TryNext
from functools import singledispatch


@singledispatch
def inspect_object(obj: Any) -> None:
    """Called when you do obj?"""
    raise TryNext


@singledispatch
def complete_object(obj: Any, prev_completions: list[str]) -> list[str]:
    """Custom completer dispatching for python objects.

    Parameters
    ----------
    obj : object
        The object to complete.
    prev_completions : list
        List of attributes discovered so far.
    This should return the list of attributes in obj. If you only wish to
    add to the attributes already discovered normally, return
    own_attrs + prev_completions.
    """
    raise TryNext
