# encoding: utf-8
"""Generic functions for extending IPython."""

from __future__ import annotations

import warnings
from typing import Any

from IPython.core.error import TryNext
from functools import singledispatch


@singledispatch
def _inspect_object(obj: Any) -> None:
    """Called when you do obj?

    .. deprecated:: 9.15
        `inspect_object` is deprecated and will be removed in a future
        version. It is no longer used within IPython, so registering
        handlers on it has no effect.
    """
    raise TryNext


def __getattr__(name: str) -> Any:
    if name == "inspect_object":
        warnings.warn(
            "inspect_object is deprecated since IPython 9.15 and will be "
            "removed in a future version. It is no longer used within "
            "IPython, so registering handlers on it has no effect.",
            DeprecationWarning,
            stacklevel=2,
        )
        return _inspect_object
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


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
