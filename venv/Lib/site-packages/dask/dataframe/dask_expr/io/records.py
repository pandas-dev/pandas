from __future__ import annotations

from dask.utils import M


def to_records(df):
    """Create Dask Array from a Dask Dataframe

    Warning: This creates a dask.array without precise shape information.
    Operations that depend on shape information, like slicing or reshaping,
    will not work.

    Examples
    --------
    >>> df.to_records()  # doctest: +SKIP

    See Also
    --------
    dask.dataframe.DataFrame.values
    dask.dataframe.from_dask_array
    """
    return df.map_partitions(M.to_records)
