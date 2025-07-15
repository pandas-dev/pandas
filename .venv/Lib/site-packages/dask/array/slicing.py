from __future__ import annotations

import bisect
import functools
import math
import warnings
from itertools import product
from numbers import Integral, Number
from operator import itemgetter

import numpy as np
from tlz import concat, memoize, merge, pluck

from dask import core
from dask._task_spec import Alias, DataNode, Task, TaskRef
from dask.array.chunk import getitem
from dask.base import is_dask_collection, tokenize
from dask.highlevelgraph import HighLevelGraph
from dask.utils import _deprecated, cached_cumsum, disable_gc, is_arraylike

colon = slice(None, None, None)


def _sanitize_index_element(ind):
    """Sanitize a one-element index."""
    if isinstance(ind, Number):
        ind2 = int(ind)
        if ind2 != ind:
            raise IndexError("Bad index.  Must be integer-like: %s" % ind)
        else:
            return ind2
    elif ind is None:
        return None
    elif is_dask_collection(ind):
        if ind.dtype.kind != "i" or ind.size != 1:
            raise IndexError(f"Bad index. Must be integer-like: {ind}")
        return ind
    else:
        raise TypeError("Invalid index type", type(ind), ind)


def sanitize_index(ind):
    """Sanitize the elements for indexing along one axis

    >>> sanitize_index([2, 3, 5])
    array([2, 3, 5])
    >>> sanitize_index([True, False, True, False])
    array([0, 2])
    >>> sanitize_index(np.array([1, 2, 3]))
    array([1, 2, 3])
    >>> sanitize_index(np.array([False, True, True]))
    array([1, 2])
    >>> type(sanitize_index(np.int32(0)))
    <class 'int'>
    >>> sanitize_index(1.0)
    1
    >>> sanitize_index(0.5)
    Traceback (most recent call last):
    ...
    IndexError: Bad index.  Must be integer-like: 0.5
    """
    from dask.array.utils import asanyarray_safe

    if ind is None:
        return None
    elif isinstance(ind, slice):
        return slice(
            _sanitize_index_element(ind.start),
            _sanitize_index_element(ind.stop),
            _sanitize_index_element(ind.step),
        )
    elif isinstance(ind, Number):
        return _sanitize_index_element(ind)
    elif is_dask_collection(ind):
        return ind
    index_array = asanyarray_safe(ind, like=ind)
    if index_array.dtype == bool:
        nonzero = np.nonzero(index_array)
        if len(nonzero) == 1:
            # If a 1-element tuple, unwrap the element
            nonzero = nonzero[0]
        if is_arraylike(nonzero):
            return nonzero
        else:
            return np.asanyarray(nonzero)
    elif np.issubdtype(index_array.dtype, np.integer):
        return index_array
    elif np.issubdtype(index_array.dtype, np.floating):
        int_index = index_array.astype(np.intp)
        if np.allclose(index_array, int_index):
            return int_index
        else:
            check_int = np.isclose(index_array, int_index)
            first_err = index_array.ravel()[np.flatnonzero(~check_int)[0]]
            raise IndexError("Bad index.  Must be integer-like: %s" % first_err)
    else:
        raise TypeError("Invalid index type", type(ind), ind)


@disable_gc()
def slice_array(out_name, in_name, blockdims, index):
    """
    Main function for array slicing

    This function makes a new dask that slices blocks along every
    dimension and aggregates (via cartesian product) each dimension's
    slices so that the resulting block slices give the same results
    as the original slice on the original structure

    Index must be a tuple.  It may contain the following types

        int, slice, list (at most one list), None

    Parameters
    ----------
    in_name - string
      This is the dask variable name that will be used as input
    out_name - string
      This is the dask variable output name
    blockshape - iterable of integers
    index - iterable of integers, slices, lists, or None

    Returns
    -------
    Dict where the keys are tuples of

        (out_name, dim_index[, dim_index[, ...]])

    and the values are

        (function, (in_name, dim_index, dim_index, ...),
                   (slice(...), [slice()[,...]])

    Also new blockdims with shapes of each block

        ((10, 10, 10, 10), (20, 20))

    Examples
    --------
    >>> from pprint import pprint
    >>> dsk, blockdims = slice_array('y', 'x', [(20, 20, 20, 20, 20)],
    ...                              (slice(10, 35),))
    >>> pprint(dsk)  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    {('y', 0): <Task ('y', 0) getitem(TaskRef(('x', 0)), (slice(10, 20, 1),))>,
     ('y', 1): <Task ('y', 1) getitem(TaskRef(('x', 1)), (slice(0, 15, 1),))>}
    >>> blockdims
    ((10, 15),)

    See Also
    --------
    This function works by successively unwrapping cases and passing down
    through a sequence of functions.

    slice_with_newaxis : handle None/newaxis case
    slice_wrap_lists : handle fancy indexing with lists
    slice_slices_and_integers : handle everything else
    """
    blockdims = tuple(map(tuple, blockdims))

    # x[:, :, :] - Punt and return old value
    if all(
        isinstance(index, slice) and index == slice(None, None, None) for index in index
    ):
        raise SlicingNoop()

    # Add in missing colons at the end as needed.  x[5] -> x[5, :, :]
    not_none_count = sum(i is not None for i in index)
    missing = len(blockdims) - not_none_count
    index += (slice(None, None, None),) * missing

    # Pass down to next function
    dsk_out, bd_out = slice_with_newaxes(out_name, in_name, blockdims, index)

    bd_out = tuple(map(tuple, bd_out))
    return dsk_out, bd_out


def slice_with_newaxes(out_name, in_name, blockdims, index):
    """
    Handle indexing with Nones

    Strips out Nones then hands off to slice_wrap_lists
    """
    # Strip Nones from index
    index2 = tuple(ind for ind in index if ind is not None)
    where_none = [i for i, ind in enumerate(index) if ind is None]
    where_none_orig = list(where_none)
    for i, x in enumerate(where_none):
        n = sum(isinstance(ind, Integral) for ind in index[:x])
        if n:
            where_none[i] -= n

    # Pass down and do work
    dsk, blockdims2 = slice_wrap_lists(
        out_name, in_name, blockdims, index2, not where_none
    )

    if where_none:
        expand = expander(where_none)
        expand_orig = expander(where_none_orig)

        # Insert ",0" into the key:  ('x', 2, 3) -> ('x', 0, 2, 0, 3)
        dsk2 = {}
        for k, v in dsk.items():
            if k[0] == out_name:
                k2 = (out_name,) + expand(k[1:], 0)
                if isinstance(v.args[1], TaskRef):
                    # positional indexing with newaxis
                    indexer = expand_orig(dsk[v.args[1].key].value[1], None)
                    tok = "shuffle-taker-" + tokenize(indexer)
                    dsk2[tok] = DataNode(tok, (1, indexer))
                    arg = TaskRef(tok)
                else:
                    arg = expand_orig(v.args[1], None)
                # raise NotImplementedError
                dsk2[k2] = Task(k2, v.func, v.args[0], arg)
            else:
                dsk2[k] = v

        # Insert (1,) into blockdims:  ((2, 2), (3, 3)) -> ((2, 2), (1,), (3, 3))
        blockdims3 = expand(blockdims2, (1,))
        return dsk2, blockdims3

    else:
        return dsk, blockdims2


def slice_wrap_lists(out_name, in_name, blockdims, index, allow_getitem_optimization):
    """
    Fancy indexing along blocked array dasks

    Handles index of type list.  Calls slice_slices_and_integers for the rest

    See Also
    --------

    take : handle slicing with lists ("fancy" indexing)
    slice_slices_and_integers : handle slicing with slices and integers
    """
    assert all(isinstance(i, (slice, list, Integral)) or is_arraylike(i) for i in index)
    if not len(blockdims) == len(index):
        raise IndexError("Too many indices for array")

    # Do we have more than one list in the index?
    where_list = [
        i for i, ind in enumerate(index) if is_arraylike(ind) and ind.ndim > 0
    ]
    if len(where_list) > 1:
        raise NotImplementedError("Don't yet support nd fancy indexing")
    # Is the single list an empty list? In this case just treat it as a zero
    # length slice
    if where_list and not index[where_list[0]].size:
        index = list(index)
        index[where_list.pop()] = slice(0, 0, 1)
        index = tuple(index)

    # No lists, hooray! just use slice_slices_and_integers
    if not where_list:
        return slice_slices_and_integers(
            out_name, in_name, blockdims, index, allow_getitem_optimization
        )

    # Replace all lists with full slices  [3, 1, 0] -> slice(None, None, None)
    index_without_list = tuple(
        slice(None, None, None) if is_arraylike(i) else i for i in index
    )

    # lists and full slices.  Just use take
    if all(is_arraylike(i) or i == slice(None, None, None) for i in index):
        axis = where_list[0]
        blockdims2, dsk3 = take(
            out_name, in_name, blockdims, index[where_list[0]], axis=axis
        )
    # Mixed case. Both slices/integers and lists. slice/integer then take
    else:
        # Do first pass without lists
        tmp = "slice-" + tokenize((out_name, in_name, blockdims, index))
        dsk, blockdims2 = slice_slices_and_integers(
            tmp,
            in_name,
            blockdims,
            index_without_list,
            allow_getitem_optimization=False,
        )

        # After collapsing some axes due to int indices, adjust axis parameter
        axis = where_list[0]
        axis2 = axis - sum(
            1 for i, ind in enumerate(index) if i < axis and isinstance(ind, Integral)
        )

        # Do work
        blockdims2, dsk2 = take(out_name, tmp, blockdims2, index[axis], axis=axis2)
        dsk3 = merge(dsk, dsk2)

    return dsk3, blockdims2


def slice_slices_and_integers(
    out_name, in_name, blockdims, index, allow_getitem_optimization=False
):
    """
    Dask array indexing with slices and integers

    See Also
    --------

    _slice_1d
    """
    from dask.array.core import unknown_chunk_message

    shape = tuple(cached_cumsum(dim, initial_zero=True)[-1] for dim in blockdims)

    for dim, ind in zip(shape, index):
        if np.isnan(dim) and ind != slice(None, None, None):
            raise ValueError(
                f"Arrays chunk sizes are unknown: {shape}{unknown_chunk_message}"
            )

    assert all(isinstance(ind, (slice, Integral)) for ind in index)
    assert len(index) == len(blockdims)

    # Get a list (for each dimension) of dicts{blocknum: slice()}
    block_slices = list(map(_slice_1d, shape, blockdims, index))
    sorted_block_slices = [sorted(i.items()) for i in block_slices]

    # (in_name, 1, 1, 2), (in_name, 1, 1, 4), (in_name, 2, 1, 2), ...
    in_names = product([in_name], *[pluck(0, s) for s in sorted_block_slices])

    # (out_name, 0, 0, 0), (out_name, 0, 0, 1), (out_name, 0, 1, 0), ...
    out_names = product(
        [out_name],
        *[
            range(len(d))[::-1] if i.step and i.step < 0 else range(len(d))
            for d, i in zip(block_slices, index)
            if not isinstance(i, Integral)
        ],
    )

    all_slices = product(*(pluck(1, s) for s in sorted_block_slices))
    if not allow_getitem_optimization:
        dsk_out = {
            out_name: (Task(out_name, getitem, TaskRef(in_name), slices))
            for out_name, in_name, slices in zip(out_names, in_names, all_slices)
        }
    else:
        empty_slice = slice(None, None, None)

        all_empty_slices = map(
            all,
            product(
                *(
                    (sp == empty_slice for sp in pluck(1, s))
                    for s in sorted_block_slices
                )
            ),
        )
        dsk_out = {
            out_name: (
                Task(out_name, getitem, TaskRef(in_name), slices)
                if not all_empty
                else Alias(out_name, in_name)
            )
            for out_name, in_name, all_empty, slices in zip(
                out_names, in_names, all_empty_slices, all_slices
            )
        }

    new_blockdims = [
        new_blockdim(d, db, i)
        for d, i, db in zip(shape, index, blockdims)
        if not isinstance(i, Integral)
    ]

    return dsk_out, new_blockdims


def _slice_1d(dim_shape, lengths, index):
    """Returns a dict of {blocknum: slice}

    This function figures out where each slice should start in each
    block for a single dimension. If the slice won't return any elements
    in the block, that block will not be in the output.

    Parameters
    ----------

    dim_shape - the number of elements in this dimension.
      This should be a positive, non-zero integer
    blocksize - the number of elements per block in this dimension
      This should be a positive, non-zero integer
    index - a description of the elements in this dimension that we want
      This might be an integer, a slice(), or an Ellipsis

    Returns
    -------

    dictionary where the keys are the integer index of the blocks that
      should be sliced and the values are the slices

    Examples
    --------

    Trivial slicing

    >>> _slice_1d(100, [60, 40], slice(None, None, None))
    {0: slice(None, None, None), 1: slice(None, None, None)}

    100 length array cut into length 20 pieces, slice 0:35

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(0, 35))
    {0: slice(None, None, None), 1: slice(0, 15, 1)}

    Support irregular blocks and various slices

    >>> _slice_1d(100, [20, 10, 10, 10, 25, 25], slice(10, 35))
    {0: slice(10, 20, 1), 1: slice(None, None, None), 2: slice(0, 5, 1)}

    Support step sizes

    >>> _slice_1d(100, [15, 14, 13], slice(10, 41, 3))
    {0: slice(10, 15, 3), 1: slice(1, 14, 3), 2: slice(2, 12, 3)}

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(0, 100, 40))  # step > blocksize
    {0: slice(0, 20, 40), 2: slice(0, 20, 40), 4: slice(0, 20, 40)}

    Also support indexing single elements

    >>> _slice_1d(100, [20, 20, 20, 20, 20], 25)
    {1: 5}

    And negative slicing

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(100, 0, -3)) # doctest: +NORMALIZE_WHITESPACE
    {4: slice(-1, -21, -3),
     3: slice(-2, -21, -3),
     2: slice(-3, -21, -3),
     1: slice(-1, -21, -3),
     0: slice(-2, -20, -3)}

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(100, 12, -3)) # doctest: +NORMALIZE_WHITESPACE
    {4: slice(-1, -21, -3),
     3: slice(-2, -21, -3),
     2: slice(-3, -21, -3),
     1: slice(-1, -21, -3),
     0: slice(-2, -8, -3)}

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(100, -12, -3))
    {4: slice(-1, -12, -3)}
    """
    chunk_boundaries = cached_cumsum(lengths)

    if isinstance(index, Integral):
        # use right-side search to be consistent with previous result
        i = bisect.bisect_right(chunk_boundaries, index)
        if i > 0:
            # the very first chunk has no relative shift
            ind = index - chunk_boundaries[i - 1]
        else:
            ind = index
        return {int(i): int(ind)}

    assert isinstance(index, slice)

    if index == colon:
        return dict.fromkeys(range(len(lengths)), colon)

    step = index.step or 1
    if step > 0:
        start = index.start or 0
        stop = index.stop if index.stop is not None else dim_shape
    else:
        start = index.start if index.start is not None else dim_shape - 1
        start = dim_shape - 1 if start >= dim_shape else start
        stop = -(dim_shape + 1) if index.stop is None else index.stop

    # posify start and stop
    if start < 0:
        start += dim_shape
    if stop < 0:
        stop += dim_shape

    d = dict()
    if step > 0:
        istart = bisect.bisect_right(chunk_boundaries, start)
        istop = bisect.bisect_left(chunk_boundaries, stop)

        # the bound is not exactly tight; make it tighter?
        istop = min(istop + 1, len(lengths))

        # jump directly to istart
        if istart > 0:
            start = start - chunk_boundaries[istart - 1]
            stop = stop - chunk_boundaries[istart - 1]

        for i in range(istart, istop):
            length = lengths[i]
            if start < length and stop > 0:
                d[i] = slice(start, min(stop, length), step)
                start = (start - length) % step
            else:
                start = start - length
            stop -= length
    else:
        rstart = start  # running start

        istart = bisect.bisect_left(chunk_boundaries, start)
        istop = bisect.bisect_right(chunk_boundaries, stop)

        # the bound is not exactly tight; make it tighter?
        istart = min(istart + 1, len(chunk_boundaries) - 1)
        istop = max(istop - 1, -1)

        for i in range(istart, istop, -1):
            chunk_stop = chunk_boundaries[i]
            # create a chunk start and stop
            if i == 0:
                chunk_start = 0
            else:
                chunk_start = chunk_boundaries[i - 1]

            # if our slice is in this chunk
            if (chunk_start <= rstart < chunk_stop) and (rstart > stop):
                d[i] = slice(
                    rstart - chunk_stop,
                    max(chunk_start - chunk_stop - 1, stop - chunk_stop),
                    step,
                )

                # compute the next running start point,
                offset = (rstart - (chunk_start - 1)) % step
                rstart = chunk_start + offset - 1

    # replace 0:20:1 with : if appropriate
    for k, v in d.items():
        if v == slice(0, lengths[k], 1):
            d[k] = slice(None, None, None)

    if not d:  # special case x[:0]
        d[0] = slice(0, 0, 1)

    return d


def partition_by_size(sizes, seq):
    """

    >>> partition_by_size([10, 20, 10], [1, 5, 9, 12, 29, 35])
    [array([1, 5, 9]), array([ 2, 19]), array([5])]
    """
    if not is_arraylike(seq):
        seq = np.asanyarray(seq)
    left = np.empty(len(sizes) + 1, dtype=int)
    left[0] = 0

    right = np.cumsum(sizes, out=left[1:])
    locations = np.empty(len(sizes) + 1, dtype=int)
    locations[0] = 0
    locations[1:] = np.searchsorted(seq, right)
    return [(seq[j:k] - l) for j, k, l in zip(locations[:-1], locations[1:], left)]


def issorted(seq):
    """Is sequence sorted?

    >>> issorted([1, 2, 3])
    np.True_
    >>> issorted([3, 1, 2])
    np.False_
    """
    if len(seq) == 0:
        return True
    return np.all(seq[:-1] <= seq[1:])


class SlicingNoop(Exception):
    """This indicates that a slicing operation is a no-op. The caller has to handle this"""

    pass


def take(outname, inname, chunks, index, axis=0):
    """Index array with an iterable of index

    Handles a single index by a single list

    Mimics ``np.take``

    >>> from pprint import pprint
    >>> chunks, dsk = take('y', 'x', [(20, 20, 20, 20)], [5, 1, 47, 3], axis=0)
    >>> chunks
    ((4,),)

    When list is sorted we still try to preserve properly sized chunks.

    >>> chunks, dsk = take('y', 'x', [(20, 20, 20, 20)], [1, 3, 5, 47], axis=0)
    >>> chunks
    ((4,),)

    """

    if not np.isnan(chunks[axis]).any():
        from dask.array._shuffle import _shuffle
        from dask.array.utils import arange_safe, asarray_safe

        arange = arange_safe(np.sum(chunks[axis]), like=index)
        if len(index) == len(arange) and np.abs(index - arange).sum() == 0:
            # TODO: This should be a real no-op, but the call stack is
            # too deep to do this efficiently for now
            chunk_tuples = product(*(range(len(c)) for i, c in enumerate(chunks)))
            graph = {
                (outname,) + c: Alias((outname,) + c, (inname,) + c)
                for c in chunk_tuples
            }
            return tuple(chunks), graph

        average_chunk_size = int(sum(chunks[axis]) / len(chunks[axis]))

        indexer = []
        index = asarray_safe(index, like=index)
        for i in range(0, len(index), average_chunk_size):
            indexer.append(index[i : i + average_chunk_size].tolist())

        token = (
            outname.split("-")[-1]
            if "-" in outname
            else tokenize(outname, chunks, index, axis)
        )
        chunks, graph = _shuffle(chunks, indexer, axis, inname, outname, token)
        return chunks, graph
    elif len(chunks[axis]) == 1:
        slices = [slice(None)] * len(chunks)
        slices[axis] = list(index)
        slices = tuple(slices)
        chunk_tuples = product(*(range(len(c)) for i, c in enumerate(chunks)))
        dsk = {
            (outname,)
            + ct: Task((outname,) + ct, getitem, TaskRef((inname,) + ct), slices)
            for ct in chunk_tuples
        }
        return chunks, dsk
    else:
        from dask.array.core import unknown_chunk_message

        raise ValueError(
            f"Array chunk size or shape is unknown. {unknown_chunk_message}"
        )


def posify_index(shape, ind):
    """Flip negative indices around to positive ones

    >>> posify_index(10, 3)
    3
    >>> posify_index(10, -3)
    7
    >>> posify_index(10, [3, -3])
    array([3, 7])

    >>> posify_index((10, 20), (3, -3))
    (3, 17)
    >>> posify_index((10, 20), (3, [3, 4, -3]))  # doctest: +NORMALIZE_WHITESPACE
    (3, array([ 3,  4, 17]))
    """
    if isinstance(ind, tuple):
        return tuple(map(posify_index, shape, ind))
    if isinstance(ind, Integral):
        if ind < 0 and not math.isnan(shape):
            return ind + shape
        else:
            return ind
    if isinstance(ind, (np.ndarray, list)) and not math.isnan(shape):
        ind = np.asanyarray(ind)
        return np.where(ind < 0, ind + shape, ind)
    return ind


@memoize
def _expander(where):
    if not where:

        def expand(seq, val):
            return seq

        return expand
    else:
        decl = """def expand(seq, val):
            return ({left}) + tuple({right})
        """
        left = []
        j = 0
        for i in range(max(where) + 1):
            if i in where:
                left.append("val, ")
            else:
                left.append("seq[%d], " % j)
                j += 1
        right = "seq[%d:]" % j
        left = "".join(left)
        decl = decl.format(**locals())
        ns = {}
        exec(compile(decl, "<dynamic>", "exec"), ns, ns)
        return ns["expand"]


def expander(where):
    """Create a function to insert value at many locations in sequence.

    >>> expander([0, 2])(['a', 'b', 'c'], 'z')
    ('z', 'a', 'z', 'b', 'c')
    """
    return _expander(tuple(where))


def new_blockdim(dim_shape, lengths, index):
    """

    >>> new_blockdim(100, [20, 10, 20, 10, 40], slice(0, 90, 2))
    [10, 5, 10, 5, 15]

    >>> new_blockdim(100, [20, 10, 20, 10, 40], [5, 1, 30, 22])
    [4]

    >>> new_blockdim(100, [20, 10, 20, 10, 40], slice(90, 10, -2))
    [16, 5, 10, 5, 4]
    """
    if index == slice(None, None, None):
        return lengths
    if isinstance(index, list):
        return [len(index)]
    assert not isinstance(index, Integral)
    pairs = sorted(_slice_1d(dim_shape, lengths, index).items(), key=itemgetter(0))
    slices = [
        slice(0, lengths[i], 1) if slc == slice(None, None, None) else slc
        for i, slc in pairs
    ]
    if isinstance(index, slice) and index.step and index.step < 0:
        slices = slices[::-1]
    return [int(math.ceil((1.0 * slc.stop - slc.start) / slc.step)) for slc in slices]


def replace_ellipsis(n, index):
    """Replace ... with slices, :, : ,:

    >>> replace_ellipsis(4, (3, Ellipsis, 2))
    (3, slice(None, None, None), slice(None, None, None), 2)

    >>> replace_ellipsis(2, (Ellipsis, None))
    (slice(None, None, None), slice(None, None, None), None)
    """
    # Careful about using in or index because index may contain arrays
    isellipsis = [i for i, ind in enumerate(index) if ind is Ellipsis]
    if not isellipsis:
        return index
    else:
        loc = isellipsis[0]
    extra_dimensions = n - (len(index) - sum(i is None for i in index) - 1)
    return (
        index[:loc] + (slice(None, None, None),) * extra_dimensions + index[loc + 1 :]
    )


def normalize_slice(idx, dim):
    """Normalize slices to canonical form

    Parameters
    ----------
    idx: slice or other index
    dim: dimension length

    Examples
    --------
    >>> normalize_slice(slice(0, 10, 1), 10)
    slice(None, None, None)
    """

    if isinstance(idx, slice):
        if math.isnan(dim):
            return idx
        start, stop, step = idx.indices(dim)
        if step > 0:
            if start == 0:
                start = None
            if stop >= dim:
                stop = None
            if step == 1:
                step = None
            if stop is not None and start is not None and stop < start:
                stop = start
        elif step < 0:
            if start >= dim - 1:
                start = None
            if stop < 0:
                stop = None
        return slice(start, stop, step)
    return idx


def normalize_index(idx, shape):
    """Normalize slicing indexes

    1.  Replaces ellipses with many full slices
    2.  Adds full slices to end of index
    3.  Checks bounding conditions
    4.  Replace multidimensional numpy arrays with dask arrays
    5.  Replaces numpy arrays with lists
    6.  Posify's integers and lists
    7.  Normalizes slices to canonical form

    Examples
    --------
    >>> normalize_index(1, (10,))
    (1,)
    >>> normalize_index(-1, (10,))
    (9,)
    >>> normalize_index([-1], (10,))
    (array([9]),)
    >>> normalize_index(slice(-3, 10, 1), (10,))
    (slice(7, None, None),)
    >>> normalize_index((Ellipsis, None), (10,))
    (slice(None, None, None), None)
    >>> normalize_index(np.array([[True, False], [False, True], [True, True]]), (3, 2))
    (dask.array<array, shape=(3, 2), dtype=bool, chunksize=(3, 2), chunktype=numpy.ndarray>,)
    """
    from dask.array import Array, from_array

    if not isinstance(idx, tuple):
        idx = (idx,)

    # if a > 1D numpy.array is provided, cast it to a dask array
    if len(idx) > 0 and len(shape) > 1:
        i = idx[0]
        if is_arraylike(i) and not isinstance(i, Array) and i.shape == shape:
            idx = (from_array(i), *idx[1:])

    idx = replace_ellipsis(len(shape), idx)
    n_sliced_dims = 0
    for i in idx:
        if hasattr(i, "ndim") and i.ndim >= 1:
            n_sliced_dims += i.ndim
        elif i is None:
            continue
        else:
            n_sliced_dims += 1

    idx = idx + (slice(None),) * (len(shape) - n_sliced_dims)
    if len([i for i in idx if i is not None]) > len(shape):
        raise IndexError("Too many indices for array")

    none_shape = []
    i = 0
    for ind in idx:
        if ind is not None:
            none_shape.append(shape[i])
            i += 1
        else:
            none_shape.append(None)

    for axis, (i, d) in enumerate(zip(idx, none_shape)):
        if d is not None:
            check_index(axis, i, d)
    idx = tuple(map(sanitize_index, idx))
    idx = tuple(map(normalize_slice, idx, none_shape))
    idx = posify_index(none_shape, idx)
    return idx


def check_index(axis, ind, dimension):
    """Check validity of index for a given dimension

    Examples
    --------
    >>> check_index(0, 3, 5)
    >>> check_index(0, 5, 5)
    Traceback (most recent call last):
    ...
    IndexError: Index 5 is out of bounds for axis 0 with size 5

    >>> check_index(1, 6, 5)
    Traceback (most recent call last):
    ...
    IndexError: Index 6 is out of bounds for axis 1 with size 5

    >>> check_index(1, -1, 5)
    >>> check_index(1, -6, 5)
    Traceback (most recent call last):
    ...
    IndexError: Index -6 is out of bounds for axis 1 with size 5

    >>> check_index(0, [1, 2], 5)
    >>> check_index(0, [6, 3], 5)
    Traceback (most recent call last):
    ...
    IndexError: Index is out of bounds for axis 0 with size 5

    >>> check_index(1, slice(0, 3), 5)

    >>> check_index(0, [True], 1)
    >>> check_index(0, [True, True], 3)
    Traceback (most recent call last):
    ...
    IndexError: Boolean array with size 2 is not long enough for axis 0 with size 3
    >>> check_index(0, [True, True, True], 1)
    Traceback (most recent call last):
    ...
    IndexError: Boolean array with size 3 is not long enough for axis 0 with size 1
    """
    if isinstance(ind, list):
        ind = np.asanyarray(ind)

    # unknown dimension, assumed to be in bounds
    if np.isnan(dimension):
        return
    elif is_dask_collection(ind):
        return
    elif is_arraylike(ind):
        if ind.dtype == bool:
            if ind.size != dimension:
                raise IndexError(
                    f"Boolean array with size {ind.size} is not long enough "
                    f"for axis {axis} with size {dimension}"
                )
        elif (ind >= dimension).any() or (ind < -dimension).any():
            raise IndexError(
                f"Index is out of bounds for axis {axis} with size {dimension}"
            )
    elif isinstance(ind, slice):
        return
    elif ind is None:
        return

    elif ind >= dimension or ind < -dimension:
        raise IndexError(
            f"Index {ind} is out of bounds for axis {axis} with size {dimension}"
        )


def slice_with_int_dask_array(x, index):
    """Slice x with at most one 1D dask arrays of ints.

    This is a helper function of :meth:`Array.__getitem__`.

    Parameters
    ----------
    x: Array
    index: tuple with as many elements as x.ndim, among which there are
           one or more Array's with dtype=int

    Returns
    -------
    tuple of (sliced x, new index)

    where the new index is the same as the input, but with slice(None)
    replaced to the original slicer where a 1D filter has been applied and
    one less element where a zero-dimensional filter has been applied.
    """
    from dask.array.core import Array

    assert len(index) == x.ndim
    fancy_indexes = [
        isinstance(idx, (tuple, list))
        or (isinstance(idx, (np.ndarray, Array)) and idx.ndim > 0)
        for idx in index
    ]
    if sum(fancy_indexes) > 1:
        raise NotImplementedError("Don't yet support nd fancy indexing")

    out_index = []
    dropped_axis_cnt = 0
    for in_axis, idx in enumerate(index):
        out_axis = in_axis - dropped_axis_cnt
        if isinstance(idx, Array) and idx.dtype.kind in "iu":
            if idx.ndim == 0:
                idx = idx[np.newaxis]
                x = slice_with_int_dask_array_on_axis(x, idx, out_axis)
                x = x[tuple(0 if i == out_axis else slice(None) for i in range(x.ndim))]
                dropped_axis_cnt += 1
            elif idx.ndim == 1:
                x = slice_with_int_dask_array_on_axis(x, idx, out_axis)
                out_index.append(slice(None))
            else:
                raise NotImplementedError(
                    "Slicing with dask.array of ints only permitted when "
                    "the indexer has zero or one dimensions"
                )
        else:
            out_index.append(idx)
    return x, tuple(out_index)


def slice_with_int_dask_array_on_axis(x, idx, axis):
    """Slice a ND dask array with a 1D dask arrays of ints along the given
    axis.

    This is a helper function of :func:`slice_with_int_dask_array`.
    """
    from dask.array import chunk
    from dask.array.core import Array, blockwise, from_array
    from dask.array.utils import asarray_safe

    assert 0 <= axis < x.ndim

    if np.isnan(x.chunks[axis]).any():
        raise NotImplementedError(
            "Slicing an array with unknown chunks with "
            "a dask.array of ints is not supported"
        )

    # Calculate the offset at which each chunk starts along axis
    # e.g. chunks=(..., (5, 3, 4), ...) -> offset=[0, 5, 8]
    offset = np.roll(np.cumsum(asarray_safe(x.chunks[axis], like=x._meta)), 1)
    offset[0] = 0
    offset = from_array(offset, chunks=1)
    # Tamper with the declared chunks of offset to make blockwise align it with
    # x[axis]
    offset = Array(
        offset.dask, offset.name, (x.chunks[axis],), offset.dtype, meta=x._meta
    )

    # Define axis labels for blockwise
    x_axes = tuple(range(x.ndim))
    idx_axes = (x.ndim,)  # arbitrary index not already in x_axes
    offset_axes = (axis,)
    p_axes = x_axes[: axis + 1] + idx_axes + x_axes[axis + 1 :]
    y_axes = x_axes[:axis] + idx_axes + x_axes[axis + 1 :]

    # Calculate the cartesian product of every chunk of x vs every chunk of idx
    p = blockwise(
        chunk.slice_with_int_dask_array,
        p_axes,
        x,
        x_axes,
        idx,
        idx_axes,
        offset,
        offset_axes,
        x_size=x.shape[axis],
        axis=axis,
        dtype=x.dtype,
        meta=x._meta,
    )

    # Aggregate on the chunks of x along axis
    y = blockwise(
        chunk.slice_with_int_dask_array_aggregate,
        y_axes,
        idx,
        idx_axes,
        p,
        p_axes,
        concatenate=True,
        x_chunks=x.chunks[axis],
        axis=axis,
        dtype=x.dtype,
        meta=x._meta,
    )
    return y


def slice_with_bool_dask_array(x, index):
    """Slice x with one or more dask arrays of bools

    This is a helper function of `Array.__getitem__`.

    Parameters
    ----------
    x: Array
    index: tuple with as many elements as x.ndim, among which there are
           one or more Array's with dtype=bool

    Returns
    -------
    tuple of (sliced x, new index)

    where the new index is the same as the input, but with slice(None)
    replaced to the original slicer when a filter has been applied.

    Note: The sliced x will have nan chunks on the sliced axes.
    """
    from dask.array.core import Array, blockwise, elemwise

    out_index = [
        slice(None) if isinstance(ind, Array) and ind.dtype == bool else ind
        for ind in index
    ]

    if len(index) == 1 and index[0].ndim == x.ndim:
        if not np.isnan(x.shape).any() and not np.isnan(index[0].shape).any():
            x = x.ravel()
            index = tuple(i.ravel() for i in index)
        elif x.ndim > 1:
            warnings.warn(
                "When slicing a Dask array of unknown chunks with a boolean mask "
                "Dask array, the output array may have a different ordering "
                "compared to the equivalent NumPy operation. This will raise an "
                "error in a future release of Dask.",
                stacklevel=3,
            )
        y = elemwise(getitem, x, *index, dtype=x.dtype)
        name = "getitem-" + tokenize(x, index)
        dsk = {(name, i): k for i, k in enumerate(core.flatten(y.__dask_keys__()))}
        chunks = ((np.nan,) * y.npartitions,)
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[y])
        return Array(graph, name, chunks, x.dtype), out_index

    if any(
        isinstance(ind, Array) and ind.dtype == bool and ind.ndim != 1 for ind in index
    ):
        raise NotImplementedError(
            "Slicing with dask.array of bools only permitted when "
            "the indexer has only one dimension or when "
            "it has the same dimension as the sliced "
            "array"
        )
    indexes = [
        ind if isinstance(ind, Array) and ind.dtype == bool else slice(None)
        for ind in index
    ]

    arginds = []
    i = 0
    for ind in indexes:
        if isinstance(ind, Array) and ind.dtype == bool:
            new = (ind, tuple(range(i, i + ind.ndim)))
            i += x.ndim
        else:
            new = (slice(None), None)
            i += 1
        arginds.append(new)

    arginds = list(concat(arginds))

    out = blockwise(
        getitem_variadic,
        tuple(range(x.ndim)),
        x,
        tuple(range(x.ndim)),
        *arginds,
        dtype=x.dtype,
    )

    chunks = []
    for ind, chunk in zip(index, out.chunks):
        if isinstance(ind, Array) and ind.dtype == bool:
            chunks.append((np.nan,) * len(chunk))
        else:
            chunks.append(chunk)
    out._chunks = tuple(chunks)
    return out, tuple(out_index)


def getitem_variadic(x, *index):
    return x[index]


def make_block_sorted_slices(index, chunks):
    """Generate blockwise-sorted index pairs for shuffling an array.

    Parameters
    ----------
    index : ndarray
        An array of index positions.
    chunks : tuple
        Chunks from the original dask array

    Returns
    -------
    index2 : ndarray
        Same values as `index`, but each block has been sorted
    index3 : ndarray
        The location of the values of `index` in `index2`

    Examples
    --------
    >>> index = np.array([6, 0, 4, 2, 7, 1, 5, 3])
    >>> chunks = ((4, 4),)
    >>> a, b = make_block_sorted_slices(index, chunks)

    Notice that the first set of 4 items are sorted, and the
    second set of 4 items are sorted.

    >>> a
    array([0, 2, 4, 6, 1, 3, 5, 7])
    >>> b
    array([3, 0, 2, 1, 7, 4, 6, 5])
    """
    from dask.array.core import slices_from_chunks

    slices = slices_from_chunks(chunks)

    if len(slices[0]) > 1:
        slices = [slice_[0] for slice_ in slices]

    offsets = np.roll(np.cumsum(chunks[0]), 1)
    offsets[0] = 0

    index2 = np.empty_like(index)
    index3 = np.empty_like(index)

    for slice_, offset in zip(slices, offsets):
        a = index[slice_]
        b = np.sort(a)
        c = offset + np.argsort(b.take(np.argsort(a)))
        index2[slice_] = b
        index3[slice_] = c

    return index2, index3


def shuffle_slice(x, index):
    """A relatively efficient way to shuffle `x` according to `index`.

    Parameters
    ----------
    x : Array
    index : ndarray
        This should be an ndarray the same length as `x` containing
        each index position in ``range(0, len(x))``.

    Returns
    -------
    Array
    """
    from dask.array.core import PerformanceWarning

    chunks1 = chunks2 = x.chunks
    if x.ndim > 1:
        chunks1 = (chunks1[0],)
    index2, index3 = make_block_sorted_slices(index, chunks1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PerformanceWarning)
        return x[index2].rechunk(chunks2)[index3]


def parse_assignment_indices(indices, shape):
    """Reformat the indices for assignment.

    The aim of this is to convert the indices to a standardised form
    so that it is easier to ascertain which chunks are touched by the
    indices.

    This function is intended to be called by `setitem_array`.

    A slice object that is decreasing (i.e. with a negative step), is
    recast as an increasing slice (i.e. with a positive step. For
    example ``slice(7,3,-1)`` would be cast as ``slice(4,8,1)``. This
    is to facilitate finding which blocks are touched by the
    index. The dimensions for which this has occurred are returned by
    the function.

    Parameters
    ----------
    indices : numpy-style indices
        Indices to array defining the elements to be assigned.
    shape : sequence of `int`
        The shape of the array.

    Returns
    -------
    parsed_indices : `list`
        The reformatted indices that are equivalent to the input
        indices.
    implied_shape : `list`
        The shape implied by the parsed indices. For instance, indices
        of ``(slice(0,2), 5, [4,1,-1])`` will have implied shape
        ``[2,3]``.
    reverse : `list`
        The positions of the dimensions whose indices in the
        parsed_indices output are reversed slices.
    implied_shape_positions: `list`
        The positions of the dimensions whose indices contribute to
        the implied_shape. For instance, indices of ``(slice(0,2), 5,
        [4,1,-1])`` will have implied_shape ``[2,3]`` and
        implied_shape_positions ``[0,2]``.

    Examples
    --------
    >>> parse_assignment_indices((slice(1, -1),), (8,))
    ([slice(1, 7, 1)], [6], [], [0])

    >>> parse_assignment_indices(([1, 2, 6, 5],), (8,))
    ([array([1, 2, 6, 5])], [4], [], [0])

    >>> parse_assignment_indices((3, slice(-1, 2, -1)), (7, 8))
    ([3, slice(3, 8, 1)], [5], [1], [1])

    >>> parse_assignment_indices((slice(-1, 2, -1), 3, [1, 2]), (7, 8, 9))
    ([slice(3, 7, 1), 3, array([1, 2])], [4, 2], [0], [0, 2])

    >>> parse_assignment_indices((slice(0, 5), slice(3, None, 2)), (5, 4))
    ([slice(0, 5, 1), slice(3, 4, 2)], [5, 1], [], [0, 1])

    >>> parse_assignment_indices((slice(0, 5), slice(3, 3, 2)), (5, 4))
    ([slice(0, 5, 1), slice(3, 3, 2)], [5, 0], [], [0])

    """
    if not isinstance(indices, tuple):
        indices = (indices,)

    # Disallow scalar boolean indexing, and also indexing by scalar
    # numpy or dask array.
    #
    # numpy allows these, but Array.__getitem__ does not yet implement
    # them properly, so disallow it for now in __setitem__
    for index in indices:
        if index is True or index is False:
            raise NotImplementedError(
                "dask does not yet implement assignment to a scalar "
                f"boolean index: {index!r}"
            )

        if (is_arraylike(index) or is_dask_collection(index)) and not index.ndim:
            raise NotImplementedError(
                "dask does not yet implement assignment to a scalar "
                f"numpy or dask array index: {index!r}"
            )

    # Initialize output variables
    implied_shape = []
    implied_shape_positions = []
    reverse = []
    parsed_indices = list(normalize_index(indices, shape))

    n_lists = 0

    for i, (index, size) in enumerate(zip(parsed_indices, shape)):
        is_slice = isinstance(index, slice)
        if is_slice:
            # Index is a slice
            start, stop, step = index.indices(size)
            if step < 0 and stop == -1:
                stop = None

            index = slice(start, stop, step)

            if step < 0:
                # When the slice step is negative, transform the
                # original slice to a new slice with a positive step
                # such that the result of the new slice is the reverse
                # of the result of the original slice.
                #
                # For example, if the original slice is slice(6,0,-2)
                # then the new slice will be slice(2,7,2).
                #
                # >>> a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
                # >>> a[slice(6, 0, -2)]
                # [6, 4, 2]
                # >>> a[slice(2, 7, 2)]
                # [2, 4, 6]
                # >>> a[slice(6, 0, -2)] == list(reversed(a[slice(2, 7, 2)]))
                # True
                start, stop, step = index.indices(size)
                step *= -1
                div, mod = divmod(start - stop - 1, step)
                div_step = div * step
                start -= div_step
                stop = start + div_step + 1

                index = slice(start, stop, step)
                reverse.append(i)

            start, stop, step = index.indices(size)

            # Note: We now have stop >= start and step >= 0

            div, mod = divmod(stop - start, step)
            if not div and not mod:
                # stop equals start => zero-sized slice for this
                # dimension
                implied_shape.append(0)
            else:
                if mod != 0:
                    div += 1

                implied_shape.append(div)
                implied_shape_positions.append(i)

        elif isinstance(index, (int, np.integer)):
            # Index is an integer
            index = int(index)

        elif is_arraylike(index) or is_dask_collection(index):
            # Index is 1-d array
            n_lists += 1
            if n_lists > 1:
                raise NotImplementedError(
                    "dask is currently limited to at most one "
                    "dimension's assignment index being a "
                    "1-d array of integers or booleans. "
                    f"Got: {indices}"
                )

            if index.ndim != 1:
                raise IndexError(
                    f"Incorrect shape ({index.shape}) of integer "
                    f"indices for dimension with size {size}"
                )

            index_size = index.size
            if (
                index.dtype == bool
                and not math.isnan(index_size)
                and index_size != size
            ):
                raise IndexError(
                    "boolean index did not match indexed array along "
                    f"dimension {i}; dimension is {size} but "
                    f"corresponding boolean dimension is {index_size}"
                )

            # Posify an integer dask array (integer numpy arrays were
            # posified in `normalize_index`)
            if is_dask_collection(index):
                if index.dtype == bool:
                    index_size = np.nan
                else:
                    index = np.where(index < 0, index + size, index)

            implied_shape.append(index_size)
            implied_shape_positions.append(i)

        parsed_indices[i] = index

    return parsed_indices, implied_shape, reverse, implied_shape_positions


def concatenate_array_chunks(x):
    """Concatenate the multidimensional chunks of an array.

    Can be used on chunks with unknown sizes.

    Parameters
    ----------
    x : dask array

    Returns
    -------
    dask array
        The concatenated dask array with one chunk.

    """
    from dask.array.core import Array, concatenate_shaped

    if x.npartitions == 1:
        return x

    name = "concatenate-shaped-" + tokenize(x)
    d = {
        (name, 0): (
            concatenate_shaped,
            list(core.flatten(x.__dask_keys__())),
            x.numblocks,
        )
    }
    graph = HighLevelGraph.from_collections(name, d, dependencies=[x])
    chunks = x.shape
    if not chunks:
        chunks = (1,)

    return Array(graph, name, chunks=(chunks,), dtype=x.dtype)


def setitem_array(out_name, array, indices, value):
    """Master function for array assignment.

    This function, that is intended to be called by
    `Array.__setitem__`, creates a new dask that assigns values to
    each block that is touched by the indices, leaving other blocks
    unchanged.

    Each block that overlaps the indices is assigned from the
    appropriate part of the assignment value. The dasks of these value
    parts are included in the output dask dictionary, as are the dasks
    of any 1-d dask array indices. This ensures that the dask array
    assignment value and any dask array indices are not computed until
    the `Array.__setitem__` operation is computed.

    The part of the assignment value applies to block is created as a
    "getitem" slice of the full assignment value.

    Parameters
    ----------
    out_name : `str`
        The dask variable output name.
    array : dask array
        The dask array that is being assigned to.
    indices : numpy-style indices
        Indices to array defining the elements to be assigned.
    value : dask array
        The assignment value, i.e. the values which will be assigned
        to elements of array.

    Returns
    -------
    dsk : `dict`
        A dictionary where the keys are new unique tokens for each
        block of the form

            (out_name, dim_index[, dim_index[, ...]])

       and the values are either

            (key,)

        or

            (setitem, key, v_key, block_indices)

        where key is an existing top-level dask key of array.

        The first case occurs when the block represented by key does
        not overlap the indices.

        The second case occurs when the block represented by key does
        overlap the indices. setitem is the chunk assignment function;
        v_key is the dask key of the the part of the assignment value
        that corresponds to the block; and block_indices are the
        assignment indices that apply to the block.

        The dictionary also includes any additional key/value pairs
        needed to define v_key, as well as any any additional
        key/value pairs needed to define dask keys contained in the
        block_indices list as references to dask array indices.

    """

    @functools.lru_cache
    def block_index_from_1d_index(dim, loc0, loc1, is_bool):
        """The positions of index elements in the range values loc0 and loc1.

        The index is the input assignment index that is defined in the
        namespace of the caller. It is assumed that negative elements
        of an integer array have already been posified.

        The non-hashable dsk is the output dask dictionary that is
        defined in the namespace of the caller.

        Parameters
        ----------
        dim : `int`
           The dimension position of the index that is used as a proxy
           for the non-hashable index to define the LRU cache key.
        loc0 : `int`
            The start index of the block along the dimension.
        loc1 : `int`
            The stop index of the block along the dimension.
        is_bool : `bool`
            Whether or not the index is of boolean data type.

        Returns
        -------
        numpy array or `str`
            If index is a numpy array then a numpy array is
            returned.

            If index is a dask array then the dask of the block index
            is inserted into the output dask dictionary, and its
            unique top-layer key is returned.

        """
        if is_bool:
            # Boolean array (dask or numpy)
            i = index[loc0:loc1]
        elif is_dask_collection(index):
            # Integer dask array
            #
            # Check for values in [loc0,loc1).
            #
            # Use the 3-argument "where" to insert place-holder
            # elements that will be searched for and removed in the
            # `setitem` function at compute time. The place-holder
            # value must be the size of the block, i.e. loc1-loc0. We
            # can't use a 1-argument "where" here because that won't
            # work if index has unknown chunk sizes.
            i = np.where((loc0 <= index) & (index < loc1), index, loc1)
            i -= loc0
        else:
            # Integer numpy array
            #
            # Check for positive values in [loc0,loc1).
            i = np.where((loc0 <= index) & (index < loc1))[0]
            i = index[i] - loc0

        if is_dask_collection(i):
            # Return dask key instead of dask array
            i = concatenate_array_chunks(i)
            dsk.update(dict(i.dask))
            i = next(flatten(i.__dask_keys__()))

        return i

    @functools.lru_cache
    def block_index_shape_from_1d_bool_index(dim, loc0, loc1):
        """Number of True index elements between positions loc0 and loc1.

        The index is the input assignment index that is defined in the
        namespace of the caller.

        Parameters
        ----------
        dim : `int`
           The dimension position of the index that is used as a proxy
           for the non-hashable index to define the LRU cache key.
        loc0 : `int`
            The start index of the block along the dimension.
        loc1 : `int`
            The stop index of the block along the dimension.

        Returns
        -------
        numpy array or dask array
            If index is a numpy array then a numpy array is
            returned.

            If index is dask array then a dask array is returned.

        """
        return np.sum(index[loc0:loc1])

    @functools.lru_cache
    def n_preceding_from_1d_bool_index(dim, loc0):
        """Number of True index elements preceding position loc0.

        The index is the input assignment index that is defined in the
        namespace of the caller.

        Parameters
        ----------
        dim : `int`
           The dimension position of the index that is used as a proxy
           for the non-hashable index to define the LRU cache key.
        loc0 : `int`
            The start index of the block along the dimension.

        Returns
        -------
        numpy array or dask array
            If index is a numpy array then a numpy array is
            returned.

            If index is dask array then a dask array is returned.

        """
        return np.sum(index[:loc0])

    @_deprecated(message=("Please use `n_preceding_from_1d_bool_index` instead."))
    def n_preceeding_from_1d_bool_index(dim, loc0):
        return n_preceding_from_1d_bool_index(dim, loc0)

    @functools.lru_cache
    def value_indices_from_1d_int_index(dim, vsize, loc0, loc1):
        """Value indices for index elements between loc0 and loc1.

        The index is the input assignment index that is defined in the
        namespace of the caller. It is assumed that negative elements
        have already been posified.

        Parameters
        ----------
        dim : `int`
           The dimension position of the index that is used as a proxy
           for the non-hashable index to define the LRU cache key.
        vsize : `int`
            The full size of the dimension of the assignment value.
        loc0 : `int`
            The start index of the block along the dimension.
        loc1 : `int`
            The stop index of the block along the dimension.

        Returns
        -------
        numpy array or dask array
            If index is a numpy array then a numpy array is
            returned.

            If index is dask array then a dask array is returned.

        """
        # Check for values in [loc0,loc1)
        if is_dask_collection(index):
            if np.isnan(index.size):
                # Integer dask array with unknown size.
                #
                # The 1-argument "where" won't work, so use the
                # 3-argument "where" and convert to a boolean
                # array. We concatenate the resulting boolean index
                # and set the chunk size (which must be the full size
                # of the dimension of the assignment value) which
                # allows the returned array to be used as a
                # __getitem__ index to the assignment value.
                i = np.where((loc0 <= index) & (index < loc1), True, False)
                i = concatenate_array_chunks(i)
                i._chunks = ((vsize,),)
            else:
                # Integer dask array with known size
                i = np.where((loc0 <= index) & (index < loc1))[0]
                i = concatenate_array_chunks(i)
        else:
            # Integer numpy array.
            i = np.where((loc0 <= index) & (index < loc1))[0]

        return i

    from dask.core import flatten

    array_shape = array.shape
    value_shape = value.shape
    value_ndim = len(value_shape)

    # Reformat input indices
    indices, implied_shape, reverse, implied_shape_positions = parse_assignment_indices(
        indices, array_shape
    )

    # Empty slices can only be assigned size 1 values
    if 0 in implied_shape and value_shape and max(value_shape) > 1:
        raise ValueError(
            f"shape mismatch: value array of shape {value_shape} "
            "could not be broadcast to indexing result "
            f"of shape {tuple(implied_shape)}"
        )

    # Set variables needed when creating the part of the assignment
    # value that applies to each block.
    #
    #  offset: The additive offset to the assignment value dimension
    #          positions that results in the positions of the
    #          corresponding dimensions in the array. offset is a
    #          non-negative integer, and a positive value means that
    #          the array has more dimensions than the assignment
    #          value.
    #
    #  value_offset: The additive offset to the array dimension
    #                positions that results in the positions of the
    #                corresponding dimensions in the assignment
    #                value. value_offset is a non-negative integer,
    #                and a positive value means that the assignment
    #                value has more dimensions than the array.
    #
    #          For example:
    #
    #          array.shape   value.shape   offset  value_offset
    #          ------------  ------------  ------  ------------
    #          (3, 4)        (3, 4)        0       0
    #          (1, 1, 3, 4)  (3, 4)        2       0
    #          (3, 4)        (1, 1, 3, 4)  0       2
    #          ------------  ------------  ------  ------------
    #
    #  array_common_shape: The shape of those dimensions of array
    #                      which correspond to dimensions of the
    #                      assignment value.
    #
    #  value_common_shape: The shape of those dimensions of the
    #                      assignment value which correspond to
    #                      dimensions of the array.
    #
    #  base_value_indices: The indices used for initialising the
    #                      selection of the part of the assignment
    #                      value that applies to each block of
    #                      array. An element of `None` will end up
    #                      being replaced by an appropriate slice on a
    #                      block-by-block basis.
    #
    # non_broadcast_dimensions: The integer positions of
    #                           array_common_shape which do not
    #                           correspond to broadcast dimensions in
    #                           the assignment value.
    #
    # Note that array_common_shape and value_common_shape may be
    # different if there are any size 1 dimensions being brodacast.
    offset = len(implied_shape) - value_ndim
    if offset >= 0:
        # The array has the same number or more dimensions than the
        # assignment value
        array_common_shape = implied_shape[offset:]
        value_common_shape = value_shape
        value_offset = 0
        reverse = [i - offset for i in reverse if i >= offset]
    else:
        # The assigmment value has more dimensions than the array
        value_offset = -offset
        array_common_shape = implied_shape
        value_common_shape = value_shape[value_offset:]
        offset = 0

        # All of the extra leading dimensions must have size 1
        if value_shape[:value_offset] != (1,) * value_offset:
            raise ValueError(
                "could not broadcast input array from shape"
                f"{value_shape} into shape {tuple(implied_shape)}"
            )

    base_value_indices = []
    non_broadcast_dimensions = []

    for i, (a, b, j) in enumerate(
        zip(array_common_shape, value_common_shape, implied_shape_positions)
    ):
        index = indices[j]
        if is_dask_collection(index) and index.dtype == bool:
            if math.isnan(b) or b <= index.size:
                base_value_indices.append(None)
                non_broadcast_dimensions.append(i)
            else:
                raise ValueError(
                    f"shape mismatch: value array dimension size of {b} is "
                    "greater then corresponding boolean index size of "
                    f"{index.size}"
                )

            continue

        if b == 1:
            base_value_indices.append(slice(None))
        elif a == b:
            base_value_indices.append(None)
            non_broadcast_dimensions.append(i)
        elif math.isnan(a):
            base_value_indices.append(None)
            non_broadcast_dimensions.append(i)
        else:
            raise ValueError(
                f"shape mismatch: value array of shape {value_shape} "
                "could not be broadcast to indexing result of shape "
                f"{tuple(implied_shape)}"
            )

    # Translate chunks tuple to a set of array locations in product
    # order
    chunks = array.chunks
    cumdims = [cached_cumsum(bds, initial_zero=True) for bds in chunks]
    array_locations = [
        [(s, s + dim) for s, dim in zip(starts, shapes)]
        for starts, shapes in zip(cumdims, chunks)
    ]
    array_locations = product(*array_locations)

    # Get the dask keys of the most recent layer in the same order as
    # the array locations.
    in_keys = list(flatten(array.__dask_keys__()))

    # Create a new "setitem" dask entry for each block in the array
    dsk = {}
    out_name = (out_name,)
    for in_key, locations in zip(in_keys, array_locations):
        # Now loop round each block dimension.
        #
        # If the block overlaps the indices then set the following
        # (which will be used to define a new dask entry):
        #
        # block_indices: The indices that will be used to assign to
        #                this block.
        #
        # block_indices_shape: The shape implied by block_indices.
        #
        # block_preceding_sizes: How many assigned elements precede
        #                        this block along each dimension that
        #                        doesn't have an integer. It is
        #                        assumed that a slice will have a
        #                        positive step, as will be the case
        #                        for reformatted indices. `None` is
        #                        used for dimensions with 1-d integer
        #                        arrays.
        block_indices = []
        block_indices_shape = []
        block_preceding_sizes = []

        local_offset = offset

        # Assume, until demonstrated otherwise, that this block
        # overlaps the assignment indices.
        overlaps = True

        # Note which dimension, if any, has 1-d integer array index
        dim_1d_int_index = None

        for dim, (index, (loc0, loc1)) in enumerate(zip(indices, locations)):
            integer_index = isinstance(index, int)
            if isinstance(index, slice):
                # Index is a slice
                stop = loc1 - loc0
                if index.stop < loc1:
                    stop -= loc1 - index.stop

                start = index.start - loc0
                if start < 0:
                    # Make start positive
                    start %= index.step

                if start >= stop:
                    # This block does not overlap the slice index
                    overlaps = False
                    break

                step = index.step
                block_index = slice(start, stop, step)
                block_index_size, rem = divmod(stop - start, step)
                if rem:
                    block_index_size += 1

                pre = index.indices(loc0)
                n_preceding, rem = divmod(pre[1] - pre[0], step)
                if rem:
                    n_preceding += 1

            elif integer_index:
                # Index is an integer
                local_offset += 1
                if not loc0 <= index < loc1:
                    # This block does not overlap the integer index
                    overlaps = False
                    break

                block_index = index - loc0

            else:
                # Index is a 1-d array
                is_bool = index.dtype == bool
                block_index = block_index_from_1d_index(dim, loc0, loc1, is_bool)
                if is_bool:
                    block_index_size = block_index_shape_from_1d_bool_index(
                        dim, loc0, loc1
                    )
                    n_preceding = n_preceding_from_1d_bool_index(dim, loc0)
                else:
                    block_index_size = None
                    n_preceding = None
                    dim_1d_int_index = dim
                    loc0_loc1 = loc0, loc1

                if not is_dask_collection(index) and not block_index.size:
                    # This block does not overlap the 1-d numpy array
                    # index
                    overlaps = False
                    break

                # Note: When the 1-d array index is a dask array then
                #       we can't tell if this block overlaps it, so we
                #       assume that it does. If it in fact doesn't
                #       overlap then the part of the assignment value
                #       that corresponds to this block will have zero
                #       size which, at compute time, will indicate to
                #       the `setitem` function to pass the block
                #       through unchanged.

            # Still here? This block overlaps the index for this
            # dimension.
            block_indices.append(block_index)
            if not integer_index:
                block_indices_shape.append(block_index_size)
                block_preceding_sizes.append(n_preceding)

        # The new dask key
        out_key = out_name + in_key[1:]

        if not overlaps:
            # This block does not overlap the indices for all
            # dimensions, so pass the block through unchanged.
            dsk[out_key] = in_key
            continue

        # Still here? Then this block overlaps the indices for all
        # dimensions and so needs to have some of its elements
        # assigned.

        # Initialise the indices of the assignment value that define
        # the parts of it which are to be assigned to this block
        value_indices = base_value_indices[:]
        for i in non_broadcast_dimensions:
            j = i + offset
            if j == dim_1d_int_index:
                # Index is a 1-d integer array
                #
                # Define index in the current namespace for use in
                # `value_indices_from_1d_int_index`
                index = indices[j]

                value_indices[i] = value_indices_from_1d_int_index(
                    dim_1d_int_index, value_shape[i + value_offset], *loc0_loc1
                )
            else:
                # Index is a slice or 1-d boolean array
                start = block_preceding_sizes[j]
                value_indices[i] = slice(start, start + block_indices_shape[j])

        # If required as a consequence of reformatting any slice
        # objects of the original indices to have a positive steps,
        # reverse the indices to assignment value.
        for i in reverse:
            size = value_common_shape[i]
            start, stop, step = value_indices[i].indices(size)
            size -= 1
            start = size - start
            stop = size - stop
            if stop < 0:
                stop = None

            value_indices[i] = slice(start, stop, -1)

        if value_ndim > len(indices):
            # The assignment value has more dimensions than array, so
            # add a leading Ellipsis to the indices of value.
            value_indices.insert(0, Ellipsis)

        # Create the part of the full assignment value that is to be
        # assigned to elements of this block and make sure that it has
        # just one chunk (so we can represent it with a single key in
        # the argument list of setitem).
        v = value[tuple(value_indices)]
        v = concatenate_array_chunks(v)
        v_key = next(flatten(v.__dask_keys__()))

        # Insert into the output dask dictionary the dask of the part
        # of assignment value for this block (not minding when we
        # overwrite any existing keys as the values will be the same).
        dsk = merge(dict(v.dask), dsk)

        # Define the assignment function for this block.
        dsk[out_key] = (setitem, in_key, v_key, block_indices)

    block_index_from_1d_index.cache_clear()
    block_index_shape_from_1d_bool_index.cache_clear()
    n_preceding_from_1d_bool_index.cache_clear()
    value_indices_from_1d_int_index.cache_clear()

    return dsk


def setitem(x, v, indices):
    """Chunk function of `setitem_array`.

    Assign v to indices of x.

    Parameters
    ----------
    x : numpy/cupy/etc. array
        The array to be assigned to.
    v : numpy/cupy/etc. array
        The values which will be assigned.
    indices : list of `slice`, `int`, or numpy array
        The indices describing the elements of x to be assigned from
        v. One index per axis.

        Note that an individual index can not be a `list`, use a 1-d
        numpy array instead.

        If a 1-d numpy array index contains the non-valid value of the
        size of the corresponding dimension of x, then those index
        elements will be removed prior to the assignment (see
        `block_index_from_1d_index` function).

    Returns
    -------
    numpy/cupy/etc. array
        A new independent array with assigned elements, unless v is
        empty (i.e. has zero size) in which case then the input array
        is returned and the indices are ignored.

    Examples
    --------
    >>> x = np.arange(8).reshape(2, 4)
    >>> setitem(x, np.array(-99), [np.array([False, True])])
    array([[  0,   1,   2,   3],
           [-99, -99, -99, -99]])
    >>> x
    array([[0, 1, 2, 3],
           [4, 5, 6, 7]])
    >>> setitem(x, np.array([-88, -99]), [slice(None), np.array([1, 3])])
    array([[  0, -88,   2, -99],
           [  4, -88,   6, -99]])
    >>> setitem(x, -x, [slice(None)])
    array([[ 0, -1, -2, -3],
           [-4, -5, -6, -7]])
    >>> x
    array([[0, 1, 2, 3],
           [4, 5, 6, 7]])
    >>> setitem(x, np.array([-88, -99]), [slice(None), np.array([4, 4, 3, 4, 1, 4])])
    array([[  0, -99,   2, -88],
           [  4, -99,   6, -88]])
    >>> value = np.where(x < 0)[0]
    >>> value.size
    0
    >>> y = setitem(x, value, [Ellipsis])
    >>> y is x
    True
    """
    if not v.size:
        return x

    # Normalize integer array indices
    for i, (index, block_size) in enumerate(zip(indices, x.shape)):
        if isinstance(index, np.ndarray) and index.dtype != bool:
            # Strip out any non-valid place-holder values
            index = index[np.where(index < block_size)[0]]
            indices[i] = index

    # If x is not masked but v is, then turn the x into a masked
    # array.
    if not np.ma.isMA(x) and np.ma.isMA(v):
        x = x.view(np.ma.MaskedArray)

    # Copy the array to guarantee no other objects are corrupted.
    # When x is the output of a scalar __getitem__ call, it is a
    # np.generic, which is read-only. Convert it to a (writeable)
    # 0-d array. x could also be a cupy array etc.
    x = np.asarray(x) if isinstance(x, np.generic) else x.copy()

    # Do the assignment
    try:
        x[tuple(indices)] = v
    except ValueError as e:
        raise ValueError(
            "shape mismatch: value array could not be broadcast to indexing result"
        ) from e

    return x
