"""JupyterLab command handler."""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import annotations

import itertools
from collections.abc import Callable
from typing import Any

from .jupyterlab_semver import Range, gt, gte, lt, lte

_CmpFn = Callable[[str | Any, str | Any, bool], bool]


def _test_overlap(
    spec1: str,
    spec2: str,
    drop_prerelease1: bool = False,
    drop_prerelease2: bool = False,
) -> bool | None:
    """Test whether two version specs overlap.

    Returns `None` if we cannot determine compatibility,
    otherwise whether there is an overlap
    """
    cmp = _compare_ranges(
        spec1,
        spec2,
        drop_prerelease1=drop_prerelease1,
        drop_prerelease2=drop_prerelease2,
    )
    if cmp is None:
        return None
    return cmp == 0


def _compare_ranges(  # noqa: C901, PLR0912
    spec1: str,
    spec2: str,
    drop_prerelease1: bool = False,
    drop_prerelease2: bool = False,
) -> int | None:
    """Test whether two version specs overlap.

    Returns `None` if we cannot determine compatibility,
    otherwise return 0 if there is an overlap, 1 if
    spec1 is lower/older than spec2, and -1 if spec1
    is higher/newer than spec2.
    """
    # Test for overlapping semver ranges.
    r1 = Range(spec1, True)  # type: ignore[no-untyped-call]
    r2 = Range(spec2, True)  # type: ignore[no-untyped-call]

    # If either range is empty, we cannot verify.
    if not r1.range or not r2.range:
        return None

    # Set return_value to a sentinel value
    return_value: int | bool | None = False

    # r1.set may be a list of ranges if the range involved an ||,
    # so we need to test for overlaps between each pair.
    for r1set, r2set in itertools.product(r1.set, r2.set):
        x1 = r1set[0].semver
        x2 = r1set[-1].semver
        y1 = r2set[0].semver
        y2 = r2set[-1].semver

        if x1.prerelease and drop_prerelease1:
            x1 = x1.inc("patch")

        if y1.prerelease and drop_prerelease2:
            y1 = y1.inc("patch")

        o1 = r1set[0].operator
        o2 = r2set[0].operator

        # We do not handle (<) specifiers.
        if o1.startswith("<") or o2.startswith("<"):
            continue

        # Handle single value specifiers.
        lx: _CmpFn = lte if x1 == x2 else lt
        ly: _CmpFn = lte if y1 == y2 else lt
        gx: _CmpFn = gte if x1 == x2 else gt
        gy: _CmpFn = gte if x1 == x2 else gt

        # Handle unbounded (>) specifiers.
        def noop(_x: str | Any, _y: str | Any, _z: bool) -> bool:
            """Return True unconditionally, representing an unbounded range."""
            return True

        if x1 == x2 and o1.startswith(">"):
            lx = noop
        if y1 == y2 and o2.startswith(">"):
            ly = noop

        # Check for overlap.
        if (
            (gte(x1, y1, True) and ly(x1, y2, True))  # type: ignore[no-untyped-call]
            or (gy(x2, y1, True) and ly(x2, y2, True))
            or (gte(y1, x1, True) and lx(y1, x2, True))  # type: ignore[no-untyped-call]
            or (gx(y2, x1, True) and lx(y2, x2, True))
        ):
            # if we ever find an overlap, we can return immediately
            return 0

        if gte(y1, x2, True):  # type: ignore[no-untyped-call]
            if return_value is False:
                # We can possibly return 1
                return_value = 1
            elif return_value == -1:
                # conflicting information, so we must return None
                return_value = None
            continue

        if gte(x1, y2, True):  # type: ignore[no-untyped-call]
            if return_value is False:
                return_value = -1
            elif return_value == 1:
                # conflicting information, so we must return None
                return_value = None
            continue

        msg = "Unexpected case comparing version ranges"
        raise AssertionError(msg)

    if return_value is False:
        return_value = None
    return return_value
