from __future__ import annotations

import datetime

CombinedOutput = type("CombinedOutput", (tuple,), {})


def overlap_chunk(func, before, after, *args, **kwargs):
    dfs = [df for df in args if isinstance(df, CombinedOutput)]
    combined, prev_part_length, next_part_length = dfs[0]

    args = [arg[0] if isinstance(arg, CombinedOutput) else arg for arg in args]

    out = func(*args, **kwargs)

    if prev_part_length is None:
        before = None
    if isinstance(before, datetime.timedelta):
        before = prev_part_length

    expansion = None
    if combined.shape[0] != 0:
        expansion = out.shape[0] // combined.shape[0]
    if before and expansion:
        before *= expansion
    if next_part_length is None:
        return out.iloc[before:]
    if isinstance(after, datetime.timedelta):
        after = next_part_length
    if after and expansion:
        after *= expansion
    return out.iloc[before:-after]


def _head_timedelta(current, next_, after):
    """Return rows of ``next_`` whose index is before the last
    observation in ``current`` + ``after``.

    Parameters
    ----------
    current : DataFrame
    next_ : DataFrame
    after : timedelta

    Returns
    -------
    overlapped : DataFrame
    """
    return next_[next_.index < (current.index.max() + after)]
