from __future__ import annotations

import operator
import warnings
from functools import partial
from numbers import Number

import numpy as np
import tlz as toolz

from dask.array.core import Array, concatenate, dotmany, from_delayed
from dask.array.creation import eye
from dask.array.random import RandomState, default_rng
from dask.array.utils import (
    array_safe,
    meta_from_array,
    solve_triangular_safe,
    svd_flip,
)
from dask.base import tokenize, wait
from dask.blockwise import blockwise
from dask.delayed import delayed
from dask.highlevelgraph import HighLevelGraph
from dask.utils import apply, derived_from


def _cumsum_blocks(it):
    total = 0
    for x in it:
        total_previous = total
        total += x
        yield (total_previous, total)


def _cumsum_part(last, new):
    return (last[1], last[1] + new)


def _nanmin(m, n):
    k_0 = min([m, n])
    k_1 = m if np.isnan(n) else n
    return k_1 if np.isnan(k_0) else k_0


def _wrapped_qr(a):
    """
    A wrapper for np.linalg.qr that handles arrays with 0 rows

    Notes: Created for tsqr so as to manage cases with uncertain
    array dimensions. In particular, the case where arrays have
    (uncertain) chunks with 0 rows.
    """
    # workaround may be removed when numpy stops rejecting edge cases
    if a.shape[0] == 0:
        return np.zeros_like(a, shape=(0, 0)), np.zeros_like(a, shape=(0, a.shape[1]))
    else:
        return np.linalg.qr(a)


def tsqr(data, compute_svd=False, _max_vchunk_size=None):
    """Direct Tall-and-Skinny QR algorithm

    As presented in:

        A. Benson, D. Gleich, and J. Demmel.
        Direct QR factorizations for tall-and-skinny matrices in
        MapReduce architectures.
        IEEE International Conference on Big Data, 2013.
        https://arxiv.org/abs/1301.1071

    This algorithm is used to compute both the QR decomposition and the
    Singular Value Decomposition.  It requires that the input array have a
    single column of blocks, each of which fit in memory.

    Parameters
    ----------
    data: Array
    compute_svd: bool
        Whether to compute the SVD rather than the QR decomposition
    _max_vchunk_size: Integer
        Used internally in recursion to set the maximum row dimension
        of chunks in subsequent recursive calls.

    Notes
    -----
    With ``k`` blocks of size ``(m, n)``, this algorithm has memory use that
    scales as ``k * n * n``.

    The implementation here is the recursive variant due to the ultimate
    need for one "single core" QR decomposition. In the non-recursive version
    of the algorithm, given ``k`` blocks, after ``k`` ``m * n`` QR
    decompositions, there will be a "single core" QR decomposition that will
    have to work with a ``(k * n, n)`` matrix.

    Here, recursion is applied as necessary to ensure that ``k * n`` is not
    larger than ``m`` (if ``m / n >= 2``). In particular, this is done
    to ensure that single core computations do not have to work on blocks
    larger than ``(m, n)``.

    Where blocks are irregular, the above logic is applied with the "height" of
    the "tallest" block used in place of ``m``.

    Consider use of the ``rechunk`` method to control this behavior.
    Taller blocks will reduce overall memory use (assuming that many of them
    still fit in memory at once).

    See Also
    --------
    dask.array.linalg.qr
        Powered by this algorithm
    dask.array.linalg.svd
        Powered by this algorithm
    dask.array.linalg.sfqr
        Variant for short-and-fat arrays
    """
    nr, nc = len(data.chunks[0]), len(data.chunks[1])
    cr_max, cc = max(data.chunks[0]), data.chunks[1][0]

    if not (data.ndim == 2 and nc == 1):  # Is a matrix  # Only one column block
        raise ValueError(
            "Input must have the following properties:\n"
            "  1. Have two dimensions\n"
            "  2. Have only one column of blocks\n\n"
            "Note: This function (tsqr) supports QR decomposition in the case of\n"
            "tall-and-skinny matrices (single column chunk/block; see qr)\n"
            "Current shape: {},\nCurrent chunksize: {}".format(
                data.shape, data.chunksize
            )
        )

    token = "-" + tokenize(data, compute_svd)

    m, n = data.shape
    numblocks = (nr, 1)

    qq, rr = np.linalg.qr(np.ones(shape=(1, 1), dtype=data.dtype))

    layers = data.__dask_graph__().layers.copy()
    dependencies = data.__dask_graph__().dependencies.copy()

    # Block qr
    name_qr_st1 = "qr" + token
    dsk_qr_st1 = blockwise(
        _wrapped_qr,
        name_qr_st1,
        "ij",
        data.name,
        "ij",
        numblocks={data.name: numblocks},
    )
    layers[name_qr_st1] = dsk_qr_st1
    dependencies[name_qr_st1] = set(data.__dask_layers__())

    # Block qr[0]
    name_q_st1 = "getitem" + token + "-q1"
    dsk_q_st1 = {
        (name_q_st1, i, 0): (operator.getitem, (name_qr_st1, i, 0), 0)
        for i in range(numblocks[0])
    }
    layers[name_q_st1] = dsk_q_st1
    dependencies[name_q_st1] = {name_qr_st1}

    # Block qr[1]
    name_r_st1 = "getitem" + token + "-r1"
    dsk_r_st1 = {
        (name_r_st1, i, 0): (operator.getitem, (name_qr_st1, i, 0), 1)
        for i in range(numblocks[0])
    }
    layers[name_r_st1] = dsk_r_st1
    dependencies[name_r_st1] = {name_qr_st1}

    # Next step is to obtain a QR decomposition for the stacked R factors, so either:
    # - gather R factors into a single core and do a QR decomposition
    # - recurse with tsqr (if single core computation too large and a-priori "meaningful
    #   reduction" possible, meaning that chunks have to be well defined)

    single_core_compute_m = nr * cc
    chunks_well_defined = not any(np.isnan(c) for cs in data.chunks for c in cs)
    prospective_blocks = np.ceil(single_core_compute_m / cr_max)
    meaningful_reduction_possible = (
        cr_max if _max_vchunk_size is None else _max_vchunk_size
    ) >= 2 * cc
    can_distribute = chunks_well_defined and int(prospective_blocks) > 1

    if chunks_well_defined and meaningful_reduction_possible and can_distribute:
        # stack chunks into blocks and recurse using tsqr

        # Prepare to stack chunks into blocks (from block qr[1])
        all_blocks = []
        curr_block = []
        curr_block_sz = 0
        for idx, a_m in enumerate(data.chunks[0]):
            m_q = a_m
            n_q = min(m_q, cc)
            m_r = n_q
            # n_r = cc
            if curr_block_sz + m_r > cr_max:
                all_blocks.append(curr_block)
                curr_block = []
                curr_block_sz = 0
            curr_block.append((idx, m_r))
            curr_block_sz += m_r
        if len(curr_block) > 0:
            all_blocks.append(curr_block)

        # R_stacked
        name_r_stacked = "stack" + token + "-r1"
        dsk_r_stacked = {
            (name_r_stacked, i, 0): (
                np.vstack,
                (tuple, [(name_r_st1, idx, 0) for idx, _ in sub_block_info]),
            )
            for i, sub_block_info in enumerate(all_blocks)
        }
        layers[name_r_stacked] = dsk_r_stacked
        dependencies[name_r_stacked] = {name_r_st1}

        # retrieve R_stacked for recursion with tsqr
        vchunks_rstacked = tuple(
            sum(map(lambda x: x[1], sub_block_info)) for sub_block_info in all_blocks
        )
        graph = HighLevelGraph(layers, dependencies)
        # dsk.dependencies[name_r_stacked] = {data.name}
        r_stacked_meta = meta_from_array(
            data, len((sum(vchunks_rstacked), n)), dtype=rr.dtype
        )
        r_stacked = Array(
            graph,
            name_r_stacked,
            shape=(sum(vchunks_rstacked), n),
            chunks=(vchunks_rstacked, n),
            meta=r_stacked_meta,
        )

        # recurse
        q_inner, r_inner = tsqr(r_stacked, _max_vchunk_size=cr_max)
        layers = toolz.merge(q_inner.dask.layers, r_inner.dask.layers)
        dependencies = toolz.merge(q_inner.dask.dependencies, r_inner.dask.dependencies)

        # Q_inner: "unstack"
        name_q_st2 = "getitem" + token + "-q2"
        dsk_q_st2 = {
            (name_q_st2, j, 0): (
                operator.getitem,
                (q_inner.name, i, 0),
                ((slice(e[0], e[1])), (slice(0, n))),
            )
            for i, sub_block_info in enumerate(all_blocks)
            for j, e in zip(
                [x[0] for x in sub_block_info],
                _cumsum_blocks([x[1] for x in sub_block_info]),
            )
        }
        layers[name_q_st2] = dsk_q_st2
        dependencies[name_q_st2] = set(q_inner.__dask_layers__())

        # R: R_inner
        name_r_st2 = "r-inner" + token
        dsk_r_st2 = {(name_r_st2, 0, 0): (r_inner.name, 0, 0)}
        layers[name_r_st2] = dsk_r_st2
        dependencies[name_r_st2] = set(r_inner.__dask_layers__())

        # Q: Block qr[0] (*) Q_inner
        name_q_st3 = "dot" + token + "-q3"
        dsk_q_st3 = blockwise(
            np.dot,
            name_q_st3,
            "ij",
            name_q_st1,
            "ij",
            name_q_st2,
            "ij",
            numblocks={name_q_st1: numblocks, name_q_st2: numblocks},
        )
        layers[name_q_st3] = dsk_q_st3
        dependencies[name_q_st3] = {name_q_st1, name_q_st2}
    else:
        # Do single core computation

        # Stacking for in-core QR computation
        to_stack = [(name_r_st1, i, 0) for i in range(numblocks[0])]
        name_r_st1_stacked = "stack" + token + "-r1"
        dsk_r_st1_stacked = {(name_r_st1_stacked, 0, 0): (np.vstack, (tuple, to_stack))}
        layers[name_r_st1_stacked] = dsk_r_st1_stacked
        dependencies[name_r_st1_stacked] = {name_r_st1}

        # In-core QR computation
        name_qr_st2 = "qr" + token + "-qr2"
        dsk_qr_st2 = blockwise(
            np.linalg.qr,
            name_qr_st2,
            "ij",
            name_r_st1_stacked,
            "ij",
            numblocks={name_r_st1_stacked: (1, 1)},
        )
        layers[name_qr_st2] = dsk_qr_st2
        dependencies[name_qr_st2] = {name_r_st1_stacked}

        # In-core qr[0]
        name_q_st2_aux = "getitem" + token + "-q2-aux"
        dsk_q_st2_aux = {
            (name_q_st2_aux, 0, 0): (operator.getitem, (name_qr_st2, 0, 0), 0)
        }
        layers[name_q_st2_aux] = dsk_q_st2_aux
        dependencies[name_q_st2_aux] = {name_qr_st2}

        chucks_are_all_known = not any(np.isnan(c) for cs in data.chunks for c in cs)
        if chucks_are_all_known:
            # when chunks are all known...
            # obtain slices on q from in-core compute (e.g.: (slice(10, 20), slice(0, 5)))
            q2_block_sizes = [min(e, n) for e in data.chunks[0]]
            block_slices = [
                (slice(e[0], e[1]), slice(0, n)) for e in _cumsum_blocks(q2_block_sizes)
            ]
            dsk_q_blockslices = {}
            deps = set()
        else:
            # when chunks are not already known...

            # request shape information: vertical chunk sizes & column dimension (n)
            name_q2bs = "shape" + token + "-q2"
            dsk_q2_shapes = {
                (name_q2bs, i): (min, (getattr, (data.name, i, 0), "shape"))
                for i in range(numblocks[0])
            }
            name_n = "getitem" + token + "-n"
            dsk_n = {
                name_n: (operator.getitem, (getattr, (data.name, 0, 0), "shape"), 1)
            }

            # cumulative sums (start, end)
            name_q2cs = "cumsum" + token + "-q2"
            dsk_q2_cumsum = {(name_q2cs, 0): [0, (name_q2bs, 0)]}

            for i in range(1, numblocks[0]):
                dsk_q2_cumsum[(name_q2cs, i)] = (
                    _cumsum_part,
                    (name_q2cs, i - 1),
                    (name_q2bs, i),
                )

            # obtain slices on q from in-core compute (e.g.: (slice(10, 20), slice(0, 5)))
            name_blockslice = "slice" + token + "-q"
            dsk_block_slices = {
                (name_blockslice, i): (
                    tuple,
                    [(apply, slice, (name_q2cs, i)), (slice, 0, name_n)],
                )
                for i in range(numblocks[0])
            }

            dsk_q_blockslices = toolz.merge(
                dsk_n, dsk_q2_shapes, dsk_q2_cumsum, dsk_block_slices
            )

            deps = {data.name}
            block_slices = [(name_blockslice, i) for i in range(numblocks[0])]

        layers["q-blocksizes" + token] = dsk_q_blockslices
        dependencies["q-blocksizes" + token] = deps

        # In-core qr[0] unstacking
        name_q_st2 = "getitem" + token + "-q2"
        dsk_q_st2 = {
            (name_q_st2, i, 0): (operator.getitem, (name_q_st2_aux, 0, 0), b)
            for i, b in enumerate(block_slices)
        }
        layers[name_q_st2] = dsk_q_st2
        if chucks_are_all_known:
            dependencies[name_q_st2] = {name_q_st2_aux}
        else:
            dependencies[name_q_st2] = {name_q_st2_aux, "q-blocksizes" + token}

        # Q: Block qr[0] (*) In-core qr[0]
        name_q_st3 = "dot" + token + "-q3"
        dsk_q_st3 = blockwise(
            np.dot,
            name_q_st3,
            "ij",
            name_q_st1,
            "ij",
            name_q_st2,
            "ij",
            numblocks={name_q_st1: numblocks, name_q_st2: numblocks},
        )
        layers[name_q_st3] = dsk_q_st3
        dependencies[name_q_st3] = {name_q_st1, name_q_st2}

        # R: In-core qr[1]
        name_r_st2 = "getitem" + token + "-r2"
        dsk_r_st2 = {(name_r_st2, 0, 0): (operator.getitem, (name_qr_st2, 0, 0), 1)}
        layers[name_r_st2] = dsk_r_st2
        dependencies[name_r_st2] = {name_qr_st2}

    if not compute_svd:
        is_unknown_m = np.isnan(data.shape[0]) or any(
            np.isnan(c) for c in data.chunks[0]
        )
        is_unknown_n = np.isnan(data.shape[1]) or any(
            np.isnan(c) for c in data.chunks[1]
        )

        if is_unknown_m and is_unknown_n:
            # assumption: m >= n
            q_shape = data.shape
            q_chunks = (data.chunks[0], (np.nan,))
            r_shape = (np.nan, np.nan)
            r_chunks = ((np.nan,), (np.nan,))
        elif is_unknown_m and not is_unknown_n:
            # assumption: m >= n
            q_shape = data.shape
            q_chunks = (data.chunks[0], (n,))
            r_shape = (n, n)
            r_chunks = (n, n)
        elif not is_unknown_m and is_unknown_n:
            # assumption: m >= n
            q_shape = data.shape
            q_chunks = (data.chunks[0], (np.nan,))
            r_shape = (np.nan, np.nan)
            r_chunks = ((np.nan,), (np.nan,))
        else:
            q_shape = (
                data.shape
                if data.shape[0] >= data.shape[1]
                else (data.shape[0], data.shape[0])
            )
            q_chunks = (
                data.chunks
                if data.shape[0] >= data.shape[1]
                else (data.chunks[0], data.chunks[0])
            )
            r_shape = (n, n) if data.shape[0] >= data.shape[1] else data.shape
            r_chunks = r_shape

        # dsk.dependencies[name_q_st3] = {data.name}
        # dsk.dependencies[name_r_st2] = {data.name}
        graph = HighLevelGraph(layers, dependencies)
        q_meta = meta_from_array(data, len(q_shape), qq.dtype)
        r_meta = meta_from_array(data, len(r_shape), rr.dtype)
        q = Array(graph, name_q_st3, shape=q_shape, chunks=q_chunks, meta=q_meta)
        r = Array(graph, name_r_st2, shape=r_shape, chunks=r_chunks, meta=r_meta)
        return q, r
    else:
        # In-core SVD computation
        name_svd_st2 = "svd" + token + "-2"
        dsk_svd_st2 = blockwise(
            np.linalg.svd,
            name_svd_st2,
            "ij",
            name_r_st2,
            "ij",
            numblocks={name_r_st2: (1, 1)},
        )
        # svd[0]
        name_u_st2 = "getitem" + token + "-u2"
        dsk_u_st2 = {(name_u_st2, 0, 0): (operator.getitem, (name_svd_st2, 0, 0), 0)}
        # svd[1]
        name_s_st2 = "getitem" + token + "-s2"
        dsk_s_st2 = {(name_s_st2, 0): (operator.getitem, (name_svd_st2, 0, 0), 1)}
        # svd[2]
        name_v_st2 = "getitem" + token + "-v2"
        dsk_v_st2 = {(name_v_st2, 0, 0): (operator.getitem, (name_svd_st2, 0, 0), 2)}
        # Q * U
        name_u_st4 = "getitem" + token + "-u4"
        dsk_u_st4 = blockwise(
            dotmany,
            name_u_st4,
            "ij",
            name_q_st3,
            "ik",
            name_u_st2,
            "kj",
            numblocks={name_q_st3: numblocks, name_u_st2: (1, 1)},
        )

        layers[name_svd_st2] = dsk_svd_st2
        dependencies[name_svd_st2] = {name_r_st2}
        layers[name_u_st2] = dsk_u_st2
        dependencies[name_u_st2] = {name_svd_st2}
        layers[name_u_st4] = dsk_u_st4
        dependencies[name_u_st4] = {name_q_st3, name_u_st2}
        layers[name_s_st2] = dsk_s_st2
        dependencies[name_s_st2] = {name_svd_st2}
        layers[name_v_st2] = dsk_v_st2
        dependencies[name_v_st2] = {name_svd_st2}

        uu, ss, vvh = np.linalg.svd(np.ones(shape=(1, 1), dtype=data.dtype))

        k = _nanmin(m, n)  # avoid RuntimeWarning with np.nanmin([m, n])

        m_u = m
        n_u = int(k) if not np.isnan(k) else k
        n_s = n_u
        m_vh = n_u
        n_vh = n
        d_vh = max(m_vh, n_vh)  # full matrix returned: but basically n
        graph = HighLevelGraph(layers, dependencies)
        u_meta = meta_from_array(data, len((m_u, n_u)), uu.dtype)
        s_meta = meta_from_array(data, len((n_s,)), ss.dtype)
        vh_meta = meta_from_array(data, len((d_vh, d_vh)), vvh.dtype)
        u = Array(
            graph,
            name_u_st4,
            shape=(m_u, n_u),
            chunks=(data.chunks[0], (n_u,)),
            meta=u_meta,
        )
        s = Array(graph, name_s_st2, shape=(n_s,), chunks=((n_s,),), meta=s_meta)
        vh = Array(
            graph, name_v_st2, shape=(d_vh, d_vh), chunks=((n,), (n,)), meta=vh_meta
        )
        return u, s, vh


def sfqr(data, name=None):
    """Direct Short-and-Fat QR

    Currently, this is a quick hack for non-tall-and-skinny matrices which
    are one chunk tall and (unless they are one chunk wide) have chunks
    that are wider than they are tall

    Q [R_1 R_2 ...] = [A_1 A_2 ...]

    it computes the factorization Q R_1 = A_1, then computes the other
    R_k's in parallel.

    Parameters
    ----------
    data: Array

    See Also
    --------
    dask.array.linalg.qr
        Main user API that uses this function
    dask.array.linalg.tsqr
        Variant for tall-and-skinny case
    """
    nr, nc = len(data.chunks[0]), len(data.chunks[1])
    cr, cc = data.chunks[0][0], data.chunks[1][0]

    if not (
        (data.ndim == 2)
        and (nr == 1)  # Is a matrix
        and (  # Has exactly one block row
            (cr <= cc)
            or (nc == 1)  # Chunking dimension on rows is at least that on cols or...
        )
    ):  # ... only one block col
        raise ValueError(
            "Input must have the following properties:\n"
            "  1. Have two dimensions\n"
            "  2. Have only one row of blocks\n"
            "  3. Either one column of blocks or (first) chunk size on cols\n"
            "     is at most that on rows (e.g.: for a 5x20 matrix,\n"
            "     chunks=((5), (8,4,8)) is fine, but chunks=((5), (4,8,8)) is not;\n"
            "     still, prefer something simple like chunks=(5,10) or chunks=5)\n\n"
            "Note: This function (sfqr) supports QR decomposition in the case\n"
            "of short-and-fat matrices (single row chunk/block; see qr)"
        )

    prefix = name or "sfqr-" + tokenize(data)
    prefix += "_"

    m, n = data.shape

    qq, rr = np.linalg.qr(np.ones(shape=(1, 1), dtype=data.dtype))

    layers = data.__dask_graph__().layers.copy()
    dependencies = data.__dask_graph__().dependencies.copy()

    # data = A = [A_1 A_rest]
    name_A_1 = prefix + "A_1"
    name_A_rest = prefix + "A_rest"
    layers[name_A_1] = {(name_A_1, 0, 0): (data.name, 0, 0)}
    dependencies[name_A_1] = set(data.__dask_layers__())
    layers[name_A_rest] = {
        (name_A_rest, 0, idx): (data.name, 0, 1 + idx) for idx in range(nc - 1)
    }
    if len(layers[name_A_rest]) > 0:
        dependencies[name_A_rest] = set(data.__dask_layers__())
    else:
        dependencies[name_A_rest] = set()

    # Q R_1 = A_1
    name_Q_R1 = prefix + "Q_R_1"
    name_Q = prefix + "Q"
    name_R_1 = prefix + "R_1"
    layers[name_Q_R1] = {(name_Q_R1, 0, 0): (np.linalg.qr, (name_A_1, 0, 0))}
    dependencies[name_Q_R1] = {name_A_1}

    layers[name_Q] = {(name_Q, 0, 0): (operator.getitem, (name_Q_R1, 0, 0), 0)}
    dependencies[name_Q] = {name_Q_R1}

    layers[name_R_1] = {(name_R_1, 0, 0): (operator.getitem, (name_Q_R1, 0, 0), 1)}
    dependencies[name_R_1] = {name_Q_R1}

    graph = HighLevelGraph(layers, dependencies)

    Q_meta = meta_from_array(data, len((m, min(m, n))), dtype=qq.dtype)
    R_1_meta = meta_from_array(data, len((min(m, n), cc)), dtype=rr.dtype)
    Q = Array(graph, name_Q, shape=(m, min(m, n)), chunks=(m, min(m, n)), meta=Q_meta)
    R_1 = Array(graph, name_R_1, shape=(min(m, n), cc), chunks=(cr, cc), meta=R_1_meta)

    # R = [R_1 Q'A_rest]
    Rs = [R_1]

    if nc > 1:
        A_rest_meta = meta_from_array(data, len((min(m, n), n - cc)), dtype=rr.dtype)
        A_rest = Array(
            graph,
            name_A_rest,
            shape=(min(m, n), n - cc),
            chunks=(cr, data.chunks[1][1:]),
            meta=A_rest_meta,
        )
        Rs.append(Q.T.dot(A_rest))

    R = concatenate(Rs, axis=1)

    return Q, R


def compression_level(n, q, n_oversamples=10, min_subspace_size=20):
    """Compression level to use in svd_compressed

    Given the size ``n`` of a space, compress that that to one of size
    ``q`` plus n_oversamples.

    The oversampling allows for greater flexibility in finding an
    appropriate subspace, a low value is often enough (10 is already a
    very conservative choice, it can be further reduced).
    ``q + oversampling`` should not be larger than ``n``.  In this
    specific implementation, ``q + n_oversamples`` is at least
    ``min_subspace_size``.

    Parameters
    ----------
    n: int
        Column/row dimension of original matrix
    q: int
        Size of the desired subspace (the actual size will be bigger,
        because of oversampling, see ``da.linalg.compression_level``)
    n_oversamples: int, default=10
        Number of oversamples used for generating the sampling matrix.
    min_subspace_size : int, default=20
        Minimum subspace size.

     Examples
    --------
    >>> compression_level(100, 10)
    20
    """
    return min(max(min_subspace_size, q + n_oversamples), n)


def compression_matrix(
    data,
    q,
    iterator="power",
    n_power_iter=0,
    n_oversamples=10,
    seed=None,
    compute=False,
):
    """Randomly sample matrix to find most active subspace

    This compression matrix returned by this algorithm can be used to
    compute both the QR decomposition and the Singular Value
    Decomposition.

    Parameters
    ----------
    data: Array
    q: int
        Size of the desired subspace (the actual size will be bigger,
        because of oversampling, see ``da.linalg.compression_level``)
    iterator: {'power', 'QR'}, default='power'
        Define the technique used for iterations to cope with flat
        singular spectra or when the input matrix is very large.
    n_power_iter: int
        Number of power iterations, useful when the singular values
        decay slowly. Error decreases exponentially as `n_power_iter`
        increases. In practice, set `n_power_iter` <= 4.
    n_oversamples: int, default=10
        Number of oversamples used for generating the sampling matrix.
        This value increases the size of the subspace computed, which is more
        accurate at the cost of efficiency.  Results are rarely sensitive to this choice
        though and in practice a value of 10 is very commonly high enough.
    compute : bool
        Whether or not to compute data at each use.
        Recomputing the input while performing several passes reduces memory
        pressure, but means that we have to compute the input multiple times.
        This is a good choice if the data is larger than memory and cheap to
        recreate.

    References
    ----------
    N. Halko, P. G. Martinsson, and J. A. Tropp.
    Finding structure with randomness: Probabilistic algorithms for
    constructing approximate matrix decompositions.
    SIAM Rev., Survey and Review section, Vol. 53, num. 2,
    pp. 217-288, June 2011
    https://arxiv.org/abs/0909.4061
    """
    if iterator not in ["power", "QR"]:
        raise ValueError(
            f"Iterator '{iterator}' not valid, must one one of ['power', 'QR']"
        )
    m, n = data.shape
    comp_level = compression_level(min(m, n), q, n_oversamples=n_oversamples)
    if isinstance(seed, RandomState):
        state = seed
    else:
        state = default_rng(seed)
    datatype = np.float64
    if (data.dtype).type in {np.float32, np.complex64}:
        datatype = np.float32
    omega = state.standard_normal(
        size=(n, comp_level), chunks=(data.chunks[1], (comp_level,))
    ).astype(datatype, copy=False)
    mat_h = data.dot(omega)
    if iterator == "power":
        for _ in range(n_power_iter):
            if compute:
                mat_h = mat_h.persist()
                wait(mat_h)
            tmp = data.T.dot(mat_h)
            if compute:
                tmp = tmp.persist()
                wait(tmp)
            mat_h = data.dot(tmp)
        q, _ = tsqr(mat_h)
    else:
        q, _ = tsqr(mat_h)
        for _ in range(n_power_iter):
            if compute:
                q = q.persist()
                wait(q)
            q, _ = tsqr(data.T.dot(q))
            if compute:
                q = q.persist()
                wait(q)
            q, _ = tsqr(data.dot(q))
    return q.T


def svd_compressed(
    a,
    k,
    iterator="power",
    n_power_iter=0,
    n_oversamples=10,
    seed=None,
    compute=False,
    coerce_signs=True,
):
    """Randomly compressed rank-k thin Singular Value Decomposition.

    This computes the approximate singular value decomposition of a large
    array.  This algorithm is generally faster than the normal algorithm
    but does not provide exact results.  One can balance between
    performance and accuracy with input parameters (see below).

    Parameters
    ----------
    a: Array
        Input array
    k: int
        Rank of the desired thin SVD decomposition.
    iterator: {'power', 'QR'}, default='power'
        Define the technique used for iterations to cope with flat
        singular spectra or when the input matrix is very large.
    n_power_iter: int, default=0
        Number of power iterations, useful when the singular values
        decay slowly. Error decreases exponentially as `n_power_iter`
        increases. In practice, set `n_power_iter` <= 4.
    n_oversamples: int, default=10
        Number of oversamples used for generating the sampling matrix.
        This value increases the size of the subspace computed, which is more
        accurate at the cost of efficiency.  Results are rarely sensitive to this choice
        though and in practice a value of 10 is very commonly high enough.
    compute : bool
        Whether or not to compute data at each use.
        Recomputing the input while performing several passes reduces memory
        pressure, but means that we have to compute the input multiple times.
        This is a good choice if the data is larger than memory and cheap to
        recreate.
    coerce_signs : bool
        Whether or not to apply sign coercion to singular vectors in
        order to maintain deterministic results, by default True.


    Examples
    --------
    >>> u, s, v = svd_compressed(x, 20)  # doctest: +SKIP

    Returns
    -------
    u:  Array, unitary / orthogonal
    s:  Array, singular values in decreasing order (largest first)
    v:  Array, unitary / orthogonal

    References
    ----------
    N. Halko, P. G. Martinsson, and J. A. Tropp.
    Finding structure with randomness: Probabilistic algorithms for
    constructing approximate matrix decompositions.
    SIAM Rev., Survey and Review section, Vol. 53, num. 2,
    pp. 217-288, June 2011
    https://arxiv.org/abs/0909.4061
    """
    comp = compression_matrix(
        a,
        k,
        iterator=iterator,
        n_power_iter=n_power_iter,
        n_oversamples=n_oversamples,
        seed=seed,
        compute=compute,
    )
    if compute:
        comp = comp.persist()
        wait(comp)
    a_compressed = comp.dot(a)
    v, s, u = tsqr(a_compressed.T, compute_svd=True)
    u = comp.T.dot(u.T)
    v = v.T
    u = u[:, :k]
    s = s[:k]
    v = v[:k, :]
    if coerce_signs:
        u, v = svd_flip(u, v)
    return u, s, v


def qr(a):
    """
    Compute the qr factorization of a matrix.

    Parameters
    ----------
    a : Array

    Returns
    -------
    q:  Array, orthonormal
    r:  Array, upper-triangular

    Examples
    --------
    >>> q, r = da.linalg.qr(x)  # doctest: +SKIP

    See Also
    --------
    numpy.linalg.qr: Equivalent NumPy Operation
    dask.array.linalg.tsqr: Implementation for tall-and-skinny arrays
    dask.array.linalg.sfqr: Implementation for short-and-fat arrays
    """

    if len(a.chunks[1]) == 1 and len(a.chunks[0]) > 1:
        return tsqr(a)
    elif len(a.chunks[0]) == 1:
        return sfqr(a)
    else:
        raise NotImplementedError(
            "qr currently supports only tall-and-skinny (single column chunk/block; see tsqr)\n"
            "and short-and-fat (single row chunk/block; see sfqr) matrices\n\n"
            "Consider use of the rechunk method. For example,\n\n"
            "x.rechunk({0: -1, 1: 'auto'}) or x.rechunk({0: 'auto', 1: -1})\n\n"
            "which rechunk one shorter axis to a single chunk, while allowing\n"
            "the other axis to automatically grow/shrink appropriately."
        )


def svd(a, coerce_signs=True):
    """
    Compute the singular value decomposition of a matrix.

    Parameters
    ----------
    a : (M, N) Array
    coerce_signs : bool
        Whether or not to apply sign coercion to singular vectors in
        order to maintain deterministic results, by default True.

    Examples
    --------

    >>> u, s, v = da.linalg.svd(x)  # doctest: +SKIP

    Returns
    -------

    u : (M, K) Array, unitary / orthogonal
        Left-singular vectors of `a` (in columns) with shape (M, K)
        where K = min(M, N).
    s : (K,) Array, singular values in decreasing order (largest first)
        Singular values of `a`.
    v : (K, N) Array, unitary / orthogonal
        Right-singular vectors of `a` (in rows) with shape (K, N)
        where K = min(M, N).

    Warnings
    --------

    SVD is only supported for arrays with chunking in one dimension.
    This requires that all inputs either contain a single column
    of chunks (tall-and-skinny) or a single row of chunks (short-and-fat).
    For arrays with chunking in both dimensions, see da.linalg.svd_compressed.

    See Also
    --------

    np.linalg.svd : Equivalent NumPy Operation
    da.linalg.svd_compressed : Randomized SVD for fully chunked arrays
    dask.array.linalg.tsqr : QR factorization for tall-and-skinny arrays
    dask.array.utils.svd_flip : Sign normalization for singular vectors
    """
    nb = a.numblocks
    if a.ndim != 2:
        raise ValueError(
            "Array must be 2D.\n"
            "Input shape: {}\n"
            "Input ndim: {}\n".format(a.shape, a.ndim)
        )
    if nb[0] > 1 and nb[1] > 1:
        raise NotImplementedError(
            "Array must be chunked in one dimension only. "
            "This function (svd) only supports tall-and-skinny or short-and-fat "
            "matrices (see da.linalg.svd_compressed for SVD on fully chunked arrays).\n"
            "Input shape: {}\n"
            "Input numblocks: {}\n".format(a.shape, nb)
        )

    # Single-chunk case
    if nb[0] == nb[1] == 1:
        m, n = a.shape
        k = min(a.shape)
        mu, ms, mv = np.linalg.svd(
            np.ones_like(a._meta, shape=(1, 1), dtype=a._meta.dtype)
        )
        u, s, v = delayed(np.linalg.svd, nout=3)(a, full_matrices=False)
        u = from_delayed(u, shape=(m, k), meta=mu)
        s = from_delayed(s, shape=(k,), meta=ms)
        v = from_delayed(v, shape=(k, n), meta=mv)
    # Multi-chunk cases
    else:
        # Tall-and-skinny case
        if nb[0] > nb[1]:
            u, s, v = tsqr(a, compute_svd=True)
            truncate = a.shape[0] < a.shape[1]
        # Short-and-fat case
        else:
            vt, s, ut = tsqr(a.T, compute_svd=True)
            u, v = ut.T, vt.T
            truncate = a.shape[0] > a.shape[1]
        # Only when necessary, remove extra singular vectors if array
        # has shape that contradicts chunking, e.g. the array is a
        # column of chunks but still has more columns than rows overall
        if truncate:
            k = min(a.shape)
            u, v = u[:, :k], v[:k, :]
    if coerce_signs:
        u, v = svd_flip(u, v)
    return u, s, v


def _solve_triangular_lower(a, b):
    return solve_triangular_safe(a, b, lower=True)


def lu(a):
    """
    Compute the lu decomposition of a matrix.

    Examples
    --------

    >>> p, l, u = da.linalg.lu(x)  # doctest: +SKIP

    Returns
    -------

    p:  Array, permutation matrix
    l:  Array, lower triangular matrix with unit diagonal.
    u:  Array, upper triangular matrix
    """
    import scipy.linalg

    if a.ndim != 2:
        raise ValueError("Dimension must be 2 to perform lu decomposition")

    xdim, ydim = a.shape
    if xdim != ydim:
        raise ValueError("Input must be a square matrix to perform lu decomposition")
    if len(set(a.chunks[0] + a.chunks[1])) != 1:
        msg = (
            "All chunks must be a square matrix to perform lu decomposition. "
            "Use .rechunk method to change the size of chunks."
        )
        raise ValueError(msg)

    vdim = len(a.chunks[0])
    hdim = len(a.chunks[1])

    token = tokenize(a)
    name_lu = "lu-lu-" + token

    name_p = "lu-p-" + token
    name_l = "lu-l-" + token
    name_u = "lu-u-" + token

    # for internal calculation
    name_p_inv = "lu-p-inv-" + token
    name_l_permuted = "lu-l-permute-" + token
    name_u_transposed = "lu-u-transpose-" + token
    name_plu_dot = "lu-plu-dot-" + token
    name_lu_dot = "lu-lu-dot-" + token

    dsk = {}
    for i in range(min(vdim, hdim)):
        target = (a.name, i, i)
        if i > 0:
            prevs = []
            for p in range(i):
                prev = name_plu_dot, i, p, p, i
                dsk[prev] = (np.dot, (name_l_permuted, i, p), (name_u, p, i))
                prevs.append(prev)
            target = (operator.sub, target, (sum, prevs))
        # diagonal block
        dsk[name_lu, i, i] = (scipy.linalg.lu, target)

        # sweep to horizontal
        for j in range(i + 1, hdim):
            target = (np.dot, (name_p_inv, i, i), (a.name, i, j))
            if i > 0:
                prevs = []
                for p in range(i):
                    prev = name_lu_dot, i, p, p, j
                    dsk[prev] = (np.dot, (name_l, i, p), (name_u, p, j))
                    prevs.append(prev)
                target = (operator.sub, target, (sum, prevs))
            dsk[name_lu, i, j] = (_solve_triangular_lower, (name_l, i, i), target)

        # sweep to vertical
        for k in range(i + 1, vdim):
            target = (a.name, k, i)
            if i > 0:
                prevs = []
                for p in range(i):
                    prev = name_plu_dot, k, p, p, i
                    dsk[prev] = (np.dot, (name_l_permuted, k, p), (name_u, p, i))
                    prevs.append(prev)
                target = (operator.sub, target, (sum, prevs))
            # solving x.dot(u) = target is equal to u.T.dot(x.T) = target.T
            dsk[name_lu, k, i] = (
                np.transpose,
                (
                    _solve_triangular_lower,
                    (name_u_transposed, i, i),
                    (np.transpose, target),
                ),
            )

    for i in range(min(vdim, hdim)):
        for j in range(min(vdim, hdim)):
            if i == j:
                dsk[name_p, i, j] = (operator.getitem, (name_lu, i, j), 0)
                dsk[name_l, i, j] = (operator.getitem, (name_lu, i, j), 1)
                dsk[name_u, i, j] = (operator.getitem, (name_lu, i, j), 2)

                # permuted l is required to be propagated to i > j blocks
                dsk[name_l_permuted, i, j] = (np.dot, (name_p, i, j), (name_l, i, j))
                dsk[name_u_transposed, i, j] = (np.transpose, (name_u, i, j))
                # transposed permutation matrix is equal to its inverse
                dsk[name_p_inv, i, j] = (np.transpose, (name_p, i, j))
            elif i > j:
                dsk[name_p, i, j] = (np.zeros, (a.chunks[0][i], a.chunks[1][j]))
                # calculations are performed using permuted l,
                # thus the result should be reverted by inverted (=transposed) p
                # to have the same row order as diagonal blocks
                dsk[name_l, i, j] = (np.dot, (name_p_inv, i, i), (name_lu, i, j))
                dsk[name_u, i, j] = (np.zeros, (a.chunks[0][i], a.chunks[1][j]))
                dsk[name_l_permuted, i, j] = (name_lu, i, j)
            else:
                dsk[name_p, i, j] = (np.zeros, (a.chunks[0][i], a.chunks[1][j]))
                dsk[name_l, i, j] = (np.zeros, (a.chunks[0][i], a.chunks[1][j]))
                dsk[name_u, i, j] = (name_lu, i, j)
                # l_permuted is not referred in upper triangulars

    pp, ll, uu = scipy.linalg.lu(np.ones(shape=(1, 1), dtype=a.dtype))
    pp_meta = meta_from_array(a, dtype=pp.dtype)
    ll_meta = meta_from_array(a, dtype=ll.dtype)
    uu_meta = meta_from_array(a, dtype=uu.dtype)

    graph = HighLevelGraph.from_collections(name_p, dsk, dependencies=[a])
    p = Array(graph, name_p, shape=a.shape, chunks=a.chunks, meta=pp_meta)

    graph = HighLevelGraph.from_collections(name_l, dsk, dependencies=[a])
    l = Array(graph, name_l, shape=a.shape, chunks=a.chunks, meta=ll_meta)

    graph = HighLevelGraph.from_collections(name_u, dsk, dependencies=[a])
    u = Array(graph, name_u, shape=a.shape, chunks=a.chunks, meta=uu_meta)

    return p, l, u


def solve_triangular(a, b, lower=False):
    """
    Solve the equation `a x = b` for `x`, assuming a is a triangular matrix.

    Parameters
    ----------
    a : (M, M) array_like
        A triangular matrix
    b : (M,) or (M, N) array_like
        Right-hand side matrix in `a x = b`
    lower : bool, optional
        Use only data contained in the lower triangle of `a`.
        Default is to use upper triangle.

    Returns
    -------
    x : (M,) or (M, N) array
        Solution to the system `a x = b`. Shape of return matches `b`.
    """

    if a.ndim != 2:
        raise ValueError("a must be 2 dimensional")
    if b.ndim <= 2:
        if a.shape[1] != b.shape[0]:
            raise ValueError("a.shape[1] and b.shape[0] must be equal")
        if a.chunks[1] != b.chunks[0]:
            msg = (
                "a.chunks[1] and b.chunks[0] must be equal. "
                "Use .rechunk method to change the size of chunks."
            )
            raise ValueError(msg)
    else:
        raise ValueError("b must be 1 or 2 dimensional")

    vchunks = len(a.chunks[1])
    hchunks = 1 if b.ndim == 1 else len(b.chunks[1])
    token = tokenize(a, b, lower)
    name = "solve-triangular-" + token

    # for internal calculation
    # (name, i, j, k, l) corresponds to a_ij.dot(b_kl)
    name_mdot = "solve-tri-dot-" + token

    def _b_init(i, j):
        if b.ndim == 1:
            return b.name, i
        else:
            return b.name, i, j

    def _key(i, j):
        if b.ndim == 1:
            return name, i
        else:
            return name, i, j

    dsk = {}
    if lower:
        for i in range(vchunks):
            for j in range(hchunks):
                target = _b_init(i, j)
                if i > 0:
                    prevs = []
                    for k in range(i):
                        prev = name_mdot, i, k, k, j
                        dsk[prev] = (np.dot, (a.name, i, k), _key(k, j))
                        prevs.append(prev)
                    target = (operator.sub, target, (sum, prevs))
                dsk[_key(i, j)] = (_solve_triangular_lower, (a.name, i, i), target)
    else:
        for i in range(vchunks):
            for j in range(hchunks):
                target = _b_init(i, j)
                if i < vchunks - 1:
                    prevs = []
                    for k in range(i + 1, vchunks):
                        prev = name_mdot, i, k, k, j
                        dsk[prev] = (np.dot, (a.name, i, k), _key(k, j))
                        prevs.append(prev)
                    target = (operator.sub, target, (sum, prevs))
                dsk[_key(i, j)] = (
                    solve_triangular_safe,
                    (a.name, i, i),
                    target,
                )

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[a, b])

    a_meta = meta_from_array(a)
    b_meta = meta_from_array(b)
    res = _solve_triangular_lower(
        array_safe([[1, 0], [1, 2]], dtype=a.dtype, like=a_meta),
        array_safe([0, 1], dtype=b.dtype, like=b_meta),
    )
    meta = meta_from_array(a, b.ndim, dtype=res.dtype)
    return Array(graph, name, shape=b.shape, chunks=b.chunks, meta=meta)


def solve(a, b, sym_pos=None, assume_a="gen"):
    """
    Solve the equation ``a x = b`` for ``x``. By default, use LU
    decomposition and forward / backward substitutions. When ``assume_a = "pos"``
    use Cholesky decomposition.

    Parameters
    ----------
    a : (M, M) array_like
        A square matrix.
    b : (M,) or (M, N) array_like
        Right-hand side matrix in ``a x = b``.
    sym_pos : bool, optional
        Assume a is symmetric and positive definite. If ``True``, use Cholesky
        decomposition.

        .. note::
            ``sym_pos`` is deprecated and will be removed in a future version.
            Use ``assume_a = 'pos'`` instead.

    assume_a : {'gen', 'pos'}, optional
        Type of data matrix. It is used to choose the dedicated solver.
        Note that Dask does not support 'her' and 'sym' types.

        .. versionchanged:: 2022.8.0
            ``assume_a = 'pos'`` was previously defined as ``sym_pos = True``.

    Returns
    -------
    x : (M,) or (M, N) Array
        Solution to the system ``a x = b``.  Shape of the return matches the
        shape of `b`.

    See Also
    --------
    scipy.linalg.solve
    """
    if sym_pos is not None:
        warnings.warn(
            "The sym_pos keyword is deprecated and should be replaced by using ``assume_a = 'pos'``."
            "``sym_pos`` will be removed in a future version.",
            category=FutureWarning,
        )
        if sym_pos:
            assume_a = "pos"

    if assume_a == "pos":
        l, u = _cholesky(a)
    elif assume_a == "gen":
        p, l, u = lu(a)
        b = p.T.dot(b)
    else:
        raise ValueError(
            f"{assume_a = } is not a recognized matrix structure, "  # noqa: E251
            "valid structures in Dask are 'pos' and 'gen'."
        )

    uy = solve_triangular(l, b, lower=True)
    return solve_triangular(u, uy)


def inv(a):
    """
    Compute the inverse of a matrix with LU decomposition and
    forward / backward substitutions.

    Parameters
    ----------
    a : array_like
        Square matrix to be inverted.

    Returns
    -------
    ainv : Array
        Inverse of the matrix `a`.
    """
    return solve(a, eye(a.shape[0], chunks=a.chunks[0][0]))


def _cholesky_lower(a):
    return np.linalg.cholesky(a)


def cholesky(a, lower=False):
    """
    Returns the Cholesky decomposition, :math:`A = L L^*` or
    :math:`A = U^* U` of a Hermitian positive-definite matrix A.

    Parameters
    ----------
    a : (M, M) array_like
        Matrix to be decomposed
    lower : bool, optional
        Whether to compute the upper or lower triangular Cholesky
        factorization.  Default is upper-triangular.

    Returns
    -------
    c : (M, M) Array
        Upper- or lower-triangular Cholesky factor of `a`.
    """

    l, u = _cholesky(a)
    if lower:
        return l
    else:
        return u


def _cholesky(a):
    """
    Private function to perform Cholesky decomposition, which returns both
    lower and upper triangulars.
    """

    if a.ndim != 2:
        raise ValueError("Dimension must be 2 to perform cholesky decomposition")

    xdim, ydim = a.shape
    if xdim != ydim:
        raise ValueError(
            "Input must be a square matrix to perform cholesky decomposition"
        )
    if len(set(a.chunks[0] + a.chunks[1])) != 1:
        msg = (
            "All chunks must be a square matrix to perform cholesky decomposition. "
            "Use .rechunk method to change the size of chunks."
        )
        raise ValueError(msg)

    vdim = len(a.chunks[0])
    hdim = len(a.chunks[1])

    token = tokenize(a)
    name = "cholesky-" + token

    # (name_lt_dot, i, j, k, l) corresponds to l_ij.dot(l_kl.T)
    name_lt_dot = "cholesky-lt-dot-" + token
    # because transposed results are needed for calculation,
    # we can build graph for upper triangular simultaneously
    name_upper = "cholesky-upper-" + token

    # calculates lower triangulars because subscriptions get simpler
    dsk = {}
    for i in range(vdim):
        for j in range(hdim):
            if i < j:
                dsk[name, i, j] = (
                    partial(np.zeros_like, shape=(a.chunks[0][i], a.chunks[1][j])),
                    meta_from_array(a),
                )
                dsk[name_upper, j, i] = (name, i, j)
            elif i == j:
                target = (a.name, i, j)
                if i > 0:
                    prevs = []
                    for p in range(i):
                        prev = name_lt_dot, i, p, i, p
                        dsk[prev] = (np.dot, (name, i, p), (name_upper, p, i))
                        prevs.append(prev)
                    target = (operator.sub, target, (sum, prevs))
                dsk[name, i, i] = (_cholesky_lower, target)
                dsk[name_upper, i, i] = (np.transpose, (name, i, i))
            else:
                # solving x.dot(L11.T) = (A21 - L20.dot(L10.T)) is equal to
                # L11.dot(x.T) = A21.T - L10.dot(L20.T)
                # L11.dot(x.T) = A12 - L10.dot(L02)
                target = (a.name, j, i)
                if j > 0:
                    prevs = []
                    for p in range(j):
                        prev = name_lt_dot, j, p, i, p
                        dsk[prev] = (np.dot, (name, j, p), (name_upper, p, i))
                        prevs.append(prev)
                    target = (operator.sub, target, (sum, prevs))
                dsk[name_upper, j, i] = (_solve_triangular_lower, (name, j, j), target)
                dsk[name, i, j] = (np.transpose, (name_upper, j, i))

    graph_upper = HighLevelGraph.from_collections(name_upper, dsk, dependencies=[a])
    graph_lower = HighLevelGraph.from_collections(name, dsk, dependencies=[a])
    a_meta = meta_from_array(a)
    cho = np.linalg.cholesky(array_safe([[1, 2], [2, 5]], dtype=a.dtype, like=a_meta))
    meta = meta_from_array(a, dtype=cho.dtype)

    lower = Array(graph_lower, name, shape=a.shape, chunks=a.chunks, meta=meta)
    # do not use .T, because part of transposed blocks are already calculated
    upper = Array(graph_upper, name_upper, shape=a.shape, chunks=a.chunks, meta=meta)
    return lower, upper


def _reverse(x):
    return x[::-1]


def lstsq(a, b):
    """
    Return the least-squares solution to a linear matrix equation using
    QR decomposition.

    Solves the equation `a x = b` by computing a vector `x` that
    minimizes the Euclidean 2-norm `|| b - a x ||^2`.  The equation may
    be under-, well-, or over- determined (i.e., the number of
    linearly independent rows of `a` can be less than, equal to, or
    greater than its number of linearly independent columns).  If `a`
    is square and of full rank, then `x` (but for round-off error) is
    the "exact" solution of the equation.

    Parameters
    ----------
    a : (M, N) array_like
        "Coefficient" matrix.
    b : {(M,), (M, K)} array_like
        Ordinate or "dependent variable" values. If `b` is two-dimensional,
        the least-squares solution is calculated for each of the `K` columns
        of `b`.

    Returns
    -------
    x : {(N,), (N, K)} Array
        Least-squares solution. If `b` is two-dimensional,
        the solutions are in the `K` columns of `x`.
    residuals : {(1,), (K,)} Array
        Sums of residuals; squared Euclidean 2-norm for each column in
        ``b - a*x``.
        If `b` is 1-dimensional, this is a (1,) shape array.
        Otherwise the shape is (K,).
    rank : Array
        Rank of matrix `a`.
    s : (min(M, N),) Array
        Singular values of `a`.
    """
    q, r = qr(a)
    x = solve_triangular(r, q.T.conj().dot(b))
    residuals = b - a.dot(x)
    residuals = abs(residuals**2).sum(axis=0, keepdims=b.ndim == 1)

    token = tokenize(a, b)

    # r must be a triangular with single block

    # rank
    rname = "lstsq-rank-" + token
    rdsk = {(rname,): (np.linalg.matrix_rank, (r.name, 0, 0))}
    graph = HighLevelGraph.from_collections(rname, rdsk, dependencies=[r])
    # rank must be an integer
    rank = Array(graph, rname, shape=(), chunks=(), dtype=int)

    # singular
    sname = "lstsq-singular-" + token
    rt = r.T.conj()
    sdsk = {
        (sname, 0): (
            _reverse,
            (np.sqrt, (np.linalg.eigvalsh, (np.dot, (rt.name, 0, 0), (r.name, 0, 0)))),
        )
    }
    graph = HighLevelGraph.from_collections(sname, sdsk, dependencies=[rt, r])
    meta = meta_from_array(residuals, 1)
    s = Array(graph, sname, shape=(r.shape[0],), chunks=r.shape[0], meta=meta)

    return x, residuals, rank, s


@derived_from(np.linalg)
def norm(x, ord=None, axis=None, keepdims=False):
    if axis is None:
        axis = tuple(range(x.ndim))
    elif isinstance(axis, Number):
        axis = (int(axis),)
    else:
        axis = tuple(axis)

    if len(axis) > 2:
        raise ValueError("Improper number of dimensions to norm.")

    if ord == "fro":
        ord = None
        if len(axis) == 1:
            raise ValueError("Invalid norm order for vectors.")

    r = abs(x)

    if ord is None:
        r = (r**2).sum(axis=axis, keepdims=keepdims) ** 0.5
    elif ord == "nuc":
        if len(axis) == 1:
            raise ValueError("Invalid norm order for vectors.")
        if x.ndim > 2:
            raise NotImplementedError("SVD based norm not implemented for ndim > 2")

        r = svd(x)[1][None].sum(keepdims=keepdims)
    elif ord == np.inf:
        if len(axis) == 1:
            r = r.max(axis=axis, keepdims=keepdims)
        else:
            r = r.sum(axis=axis[1], keepdims=True).max(axis=axis[0], keepdims=True)
            if keepdims is False:
                r = r.squeeze(axis=axis)
    elif ord == -np.inf:
        if len(axis) == 1:
            r = r.min(axis=axis, keepdims=keepdims)
        else:
            r = r.sum(axis=axis[1], keepdims=True).min(axis=axis[0], keepdims=True)
            if keepdims is False:
                r = r.squeeze(axis=axis)
    elif ord == 0:
        if len(axis) == 2:
            raise ValueError("Invalid norm order for matrices.")

        r = (r != 0).astype(r.dtype).sum(axis=axis, keepdims=keepdims)
    elif ord == 1:
        if len(axis) == 1:
            r = r.sum(axis=axis, keepdims=keepdims)
        else:
            r = r.sum(axis=axis[0], keepdims=True).max(axis=axis[1], keepdims=True)
            if keepdims is False:
                r = r.squeeze(axis=axis)
    elif len(axis) == 2 and ord == -1:
        r = r.sum(axis=axis[0], keepdims=True).min(axis=axis[1], keepdims=True)
        if keepdims is False:
            r = r.squeeze(axis=axis)
    elif len(axis) == 2 and ord == 2:
        if x.ndim > 2:
            raise NotImplementedError("SVD based norm not implemented for ndim > 2")
        r = svd(x)[1][None].max(keepdims=keepdims)
    elif len(axis) == 2 and ord == -2:
        if x.ndim > 2:
            raise NotImplementedError("SVD based norm not implemented for ndim > 2")
        r = svd(x)[1][None].min(keepdims=keepdims)
    else:
        if len(axis) == 2:
            raise ValueError("Invalid norm order for matrices.")

        r = (r**ord).sum(axis=axis, keepdims=keepdims) ** (1.0 / ord)

    return r
