from __future__ import annotations

import contextlib
import logging
import tempfile

import numpy as np

from dask import config
from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.dispatch import group_split_dispatch, hash_object_dispatch

logger = logging.getLogger(__name__)


class maybe_buffered_partd:
    """
    If serialized, will return non-buffered partd. Otherwise returns a buffered partd
    """

    def __init__(self, encode_cls=None, buffer=True, tempdir=None):
        self.tempdir = tempdir or config.get("temporary_directory", None)
        self.buffer = buffer
        self.compression = config.get("dataframe.shuffle.compression", None)
        self.encode_cls = encode_cls
        if encode_cls is None:
            import partd

            self.encode_cls = partd.PandasBlocks

    def __reduce__(self):
        if self.tempdir:
            return (maybe_buffered_partd, (self.encode_cls, False, self.tempdir))
        else:
            return (maybe_buffered_partd, (self.encode_cls, False))

    def __call__(self, *args, **kwargs):
        import partd

        path = tempfile.mkdtemp(suffix=".partd", dir=self.tempdir)

        try:
            partd_compression = (
                getattr(partd.compressed, self.compression)
                if self.compression
                else None
            )
        except AttributeError as e:
            raise ImportError(
                "Not able to import and load {} as compression algorithm."
                "Please check if the library is installed and supported by Partd.".format(
                    self.compression
                )
            ) from e
        file = partd.File(path)
        partd.file.cleanup_files.append(path)
        # Envelope partd file with compression, if set and available
        if partd_compression:
            file = partd_compression(file)
        if self.buffer:
            return self.encode_cls(partd.Buffer(partd.Dict(), file))
        else:
            return self.encode_cls(file)


########################################################
# Various convenience functions to be run by the above #
########################################################


def partitioning_index(df, npartitions, cast_dtype=None):
    """
    Computes a deterministic index mapping each record to a partition.

    Identical rows are mapped to the same partition.

    Parameters
    ----------
    df : DataFrame/Series/Index
    npartitions : int
        The number of partitions to group into.
    cast_dtype : dtype, optional
        The dtype to cast to to avoid nullability issues

    Returns
    -------
    partitions : ndarray
        An array of int64 values mapping each record to a partition.
    """
    if cast_dtype is not None:
        # Fixme: astype raises with strings in numeric columns, but raising
        # here might be very noisy
        df = df.astype(cast_dtype, errors="ignore")
    res = hash_object_dispatch(df, index=False) % int(npartitions)
    # Note: Use a signed integer since pandas is more efficient at handling
    # this since there is not always a fastpath for uints
    return res.astype(np.min_scalar_type(-(npartitions - 1)))


def barrier(args):
    list(args)
    return 0


def collect(p, part, meta, barrier_token):
    """Collect partitions from partd, yield dataframes"""
    with ensure_cleanup_on_exception(p):
        res = p.get(part)
        return res if len(res) > 0 else meta


def set_partitions_pre(s, divisions, ascending=True, na_position="last"):
    try:
        if ascending:
            partitions = divisions.searchsorted(s, side="right") - 1
        else:
            partitions = len(divisions) - divisions.searchsorted(s, side="right") - 1
    except (TypeError, ValueError):
        # `searchsorted` fails if either `divisions` or `s` contains nulls and strings
        partitions = np.empty(len(s), dtype="int32")
        not_null = s.notna()
        divisions_notna = divisions[divisions.notna()]
        if ascending:
            partitions[not_null] = (
                divisions_notna.searchsorted(s[not_null], side="right") - 1
            )
        else:
            partitions[not_null] = (
                len(divisions)
                - divisions_notna.searchsorted(s[not_null], side="right")
                - 1
            )
    partitions[(partitions < 0) | (partitions >= len(divisions) - 1)] = (
        len(divisions) - 2 if ascending else 0
    )
    nas = s.isna()
    # We could be a ndarray already (datetime dtype)
    nas = getattr(nas, "values", nas)
    partitions[nas] = len(divisions) - 2 if na_position == "last" else 0
    return partitions


def shuffle_group_2(df, cols, ignore_index, nparts):
    if not len(df):
        return {}, df

    if isinstance(cols, str):
        cols = [cols]

    if cols and cols[0] == "_partitions":
        ind = df[cols[0]].astype(np.int32)
    else:
        ind = (
            hash_object_dispatch(df[cols] if cols else df, index=False) % int(nparts)
        ).astype(np.int32)

    n = ind.max() + 1
    result2 = group_split_dispatch(df, ind, n, ignore_index=ignore_index)
    return result2, df.iloc[:0]


def shuffle_group_get(g_head, i):
    g, head = g_head
    if i in g:
        return g[i]
    else:
        return head


def shuffle_group(df, cols, stage, k, npartitions, ignore_index, nfinal):
    """Splits dataframe into groups

    The group is determined by their final partition, and which stage we are in
    in the shuffle

    Parameters
    ----------
    df: DataFrame
    cols: str or list
        Column name(s) on which to split the dataframe. If ``cols`` is not
        "_partitions", hashing will be used to determine target partition
    stage: int
        We shuffle dataframes with many partitions we in a few stages to avoid
        a quadratic number of tasks.  This number corresponds to which stage
        we're in, starting from zero up to some small integer
    k: int
        Desired number of splits from this dataframe
    npartition: int
        Total number of output partitions for the full dataframe
    nfinal: int
        Total number of output partitions after repartitioning

    Returns
    -------
    out: Dict[int, DataFrame]
        A dictionary mapping integers in {0..k} to dataframes such that the
        hash values of ``df[col]`` are well partitioned.
    """
    if isinstance(cols, str):
        cols = [cols]

    if cols and cols[0] == "_partitions":
        ind = df[cols[0]]
    else:
        ind = hash_object_dispatch(df[cols] if cols else df, index=False)
        if nfinal and nfinal != npartitions:
            ind = ind % int(nfinal)

    typ = np.min_scalar_type(npartitions * 2)
    # Here we convert the final output index `ind` into the output index
    # for the current stage.
    kwargs = {} if PANDAS_GE_300 else {"copy": False}
    ind = (ind % npartitions).astype(typ, **kwargs) // k**stage % k
    return group_split_dispatch(df, ind, k, ignore_index=ignore_index)


@contextlib.contextmanager
def ensure_cleanup_on_exception(p):
    """Ensure a partd.File is cleaned up.

    We have several tasks referring to a `partd.File` instance. We want to
    ensure that the file is cleaned up if and only if there's an exception
    in the tasks using the `partd.File`.
    """
    try:
        yield
    except Exception:
        # the function (e.g. shuffle_group_3) had an internal exception.
        # We'll cleanup our temporary files and re-raise.
        try:
            p.drop()
        except Exception:
            logger.exception("ignoring exception in ensure_cleanup_on_exception")
        raise


def drop_overlap(df, index):
    return df.drop(index) if index in df.index else df


def get_overlap(df, index):
    return df.loc[[index]] if index in df.index else df._constructor()
