from __future__ import annotations

from dask.dataframe.dask_expr import FrameBase
from dask.dataframe.io.io import _df_to_bag
from dask.tokenize import tokenize


def to_bag(df, index=False, format="tuple"):
    """Create Dask Bag from a Dask DataFrame

    Parameters
    ----------
    index : bool, optional
        If True, the elements are tuples of ``(index, value)``, otherwise
        they're just the ``value``.  Default is False.
    format : {"tuple", "dict", "frame"}, optional
        Whether to return a bag of tuples, dictionaries, or
        dataframe-like objects. Default is "tuple". If "frame",
        the original partitions of ``df`` will not be transformed
        in any way.


    Examples
    --------
    >>> bag = df.to_bag()  # doctest: +SKIP
    """
    from dask.bag.core import Bag

    df = df.optimize()

    if not isinstance(df, FrameBase):
        raise TypeError("df must be either DataFrame or Series")
    name = "to_bag-" + tokenize(df._name, index, format)
    if format == "frame":
        dsk = df.dask
        name = df._name
    else:
        dsk = {
            (name, i): (_df_to_bag, block, index, format)
            for (i, block) in enumerate(df.__dask_keys__())
        }
        dsk.update(df.__dask_graph__())
    return Bag(dsk, name, df.npartitions)
