# encoding: utf-8
"""Utilities for working with data structures like lists, dicts and tuples.
"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

import warnings
from collections.abc import Iterable, Sequence
from typing import TypeVar


T = TypeVar("T")


def uniq_stable(elems: Iterable[T]) -> list[T]:
    """uniq_stable(elems) -> list

    .. deprecated:: 9.8
        This function is deprecated and will be removed in a future version.
        It is not used within IPython and was never part of the public API.

    Return from an iterable, a list of all the unique elements in the input,
    but maintaining the order in which they first appear.

    Note: All elements in the input must be hashable for this routine
    to work, as it internally uses a set for efficiency reasons.
    """
    warnings.warn(
        "uniq_stable is deprecated since IPython 9.8 and will be removed in a future version. "
        "It was never part of the public API.",
        DeprecationWarning,
        stacklevel=2,
    )
    seen: set[T] = set()
    result: list[T] = []
    for x in elems:
        if x not in seen:
            seen.add(x)
            result.append(x)
    return result


def chop(seq: Sequence[T], size: int) -> list[Sequence[T]]:
    """Chop a sequence into chunks of the given size."""
    return [seq[i : i + size] for i in range(0, len(seq), size)]
