"""
Deferred gathers of string and binary Arrow data, for ``ArrowExtensionArray``.

A ``take`` or a boolean mask of such an array stores the parent array and the
selection (see ``PendingTake``) and gathers on first use. Everything here
works on ``pyarrow.ChunkedArray`` and NumPy arrays only. ``ArrowExtensionArray``
holds a ``PendingTake`` in ``_pending_take`` instead of ``_pa_array`` until
the gather runs, and hooks into it at these points:

* creating a deferred array: ``take`` and ``__getitem__`` with a mask call
  ``defer``, or ``PendingTake.select`` if the array is already deferred;
* gathering: ``__getattr__`` runs ``PendingTake.materialize`` on the first
  access to ``_pa_array``, so every other method works unchanged;
* pickling: ``__getstate__`` gathers, and never pickles the parent.

and, to avoid gathering where it is not needed:

* ``__len__`` returns the length of the selection;
* ``view`` shares the ``PendingTake``, so views gather once between them;
* ``__getitem__`` with a slice calls ``PendingTake.select_slice``, and with an
  integer ``PendingTake.value_at``;
* ``copy`` calls ``PendingTake.gather_copy``.

A frame applies a selection to each of its columns separately, but after the
first selection the deferred columns all hold the same selection array (see
``BlockManager.reindex_indexer``). Inside ``shared_selections``, the columns
compose it with the next selection's positions once between them, and keep
one array for the result, so a chain of selections costs about the same
whatever the number of deferred columns.
"""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
import threading
from typing import (
    TYPE_CHECKING,
    TypeVar,
    cast,
)

import numpy as np

from pandas.compat import HAS_PYARROW

if HAS_PYARROW:
    import pyarrow as pa

if TYPE_CHECKING:
    from collections.abc import Iterator

    from pandas._typing import npt

    _Selection = npt.NDArray[np.intp] | npt.NDArray[np.bool_]
    # the two arrays composed, held so their ids cannot be reused while the
    # entry exists, and the composed selection
    _SharedEntry = tuple[object, object, _Selection]

_ScalarT = TypeVar("_ScalarT", bound=np.generic)

# the selections composed in the innermost ``shared_selections`` block, by
# the ids of the two arrays composed; None outside one
_shared: ContextVar[dict[tuple[int, int], _SharedEntry] | None] = ContextVar(
    "_shared_selections", default=None
)


@contextmanager
def shared_selections() -> Iterator[None]:
    """
    Within this block, deferred takes that compose the same selection array
    with the same positions array share the result, if the positions can
    never change (see ``freeze``).

    For a frame operation that selects the same rows from each column. The
    composed selections are kept only until the block exits, along with the
    arrays they were composed from, so an id cannot be reused meanwhile.
    """
    token = _shared.set({})
    try:
        yield
    finally:
        _shared.reset(token)


def _take_by_chunk(
    parent: pa.ChunkedArray, indices: npt.NDArray[np.intp]
) -> pa.ChunkedArray:
    """
    ``parent.take(indices)`` for sorted ``indices``, one chunk at a time.

    The result keeps the chunk boundaries of ``parent``, as ``filter`` does.
    """
    lengths = np.fromiter((len(chunk) for chunk in parent.chunks), dtype=np.intp)
    ends = np.cumsum(lengths)
    pieces = np.split(indices, np.searchsorted(indices, ends[:-1]))
    chunks = [
        chunk.take(pa.array(piece - (end - length)))
        for chunk, piece, end, length in zip(
            parent.chunks, pieces, ends, lengths, strict=True
        )
        if len(piece)
    ]
    return pa.chunked_array(chunks, type=parent.type)


def _filter(parent: pa.ChunkedArray, mask: npt.NDArray[np.bool_]) -> pa.ChunkedArray:
    """``parent.filter(mask)``, applied chunk by chunk."""
    # pa.array(mask) converts one value at a time; packing the bits with
    # NumPy is over 10x faster
    bits = pa.py_buffer(np.packbits(mask, bitorder="little"))
    return parent.filter(pa.Array.from_buffers(pa.bool_(), len(mask), [None, bits]))


def _gather(parent: pa.ChunkedArray, indices: npt.NDArray[np.intp]) -> pa.ChunkedArray:
    """
    ``parent.take(indices)``, without concatenating ``parent`` where possible.

    ``ChunkedArray.take`` concatenates the chunks first. That copies all of
    ``parent``, however few values are taken, and overflows the 32-bit
    offsets of ``string`` and ``binary`` once ``parent`` holds more than
    2 GiB. Positions that only increase, which every boolean mask of a
    frame produces, select what a mask does, so they are gathered with
    ``filter``, which works chunk by chunk. Other positions are taken as
    they are, and only if that overflows, gathered chunk by chunk, which
    loops over the chunks in Python. The overflow is found before any
    values are copied, so trying costs little.
    """
    if parent.num_chunks <= 1:
        return parent.take(pa.array(indices))
    if bool((indices[1:] > indices[:-1]).all()):
        mask = np.zeros(len(parent), dtype=bool)
        mask[indices] = True
        return _filter(parent, mask)
    try:
        return parent.take(pa.array(indices))
    except pa.ArrowInvalid:
        if pa.types.is_large_string(parent.type) or pa.types.is_large_binary(
            parent.type
        ):
            raise
    if bool((indices[1:] >= indices[:-1]).all()):
        return _take_by_chunk(parent, indices)
    # gather in sorted order, then restore the requested order; this only
    # fails if the result itself is over 2 GiB
    order = np.argsort(indices, kind="stable")
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    return _take_by_chunk(parent, indices[order]).take(pa.array(inverse))


def _select(
    parent: pa.ChunkedArray,
    selection: npt.NDArray[np.intp] | npt.NDArray[np.bool_],
) -> pa.ChunkedArray:
    """
    The rows of ``parent`` at ``selection``: positions, or a boolean mask.

    A mask is applied with ``filter`` and positions with ``take``, as an
    eager selection would: ``filter`` copies runs of selected values at
    once, which is much faster than ``take`` for a mask that keeps most rows.
    """
    if selection.dtype == np.bool_:
        return _filter(parent, cast("npt.NDArray[np.bool_]", selection))
    # a step slice of a deferred take leaves a strided view, which pa.array
    # must copy, and copies more slowly than NumPy (for 1.7M positions,
    # 2.9 ms against 1.7 ms); a contiguous array is used without a copy
    return _gather(parent, np.ascontiguousarray(selection, dtype=np.intp))


def _drop_unused_validity(values: pa.ChunkedArray) -> pa.ChunkedArray:
    """
    ``values`` of a string or binary type, without the validity bitmap of
    any chunk that has no nulls.

    ``take`` and ``filter`` allocate a bitmap even when nothing is null,
    which a copy of an eager array drops, since ``pa.concat_arrays`` does.
    Dropping it, which copies nothing, gives a copy of a deferred take the
    same buffers, and so the same ``nbytes``.
    """
    chunks = [
        pa.Array.from_buffers(
            chunk.type,
            len(chunk),
            [None, *chunk.buffers()[1:]],
            null_count=0,
            offset=chunk.offset,
        )
        if chunk.null_count == 0 and chunk.buffers()[0] is not None
        else chunk
        for chunk in values.chunks
    ]
    return pa.chunked_array(chunks, type=values.type)


class PendingTake:
    """
    A ``take`` or boolean mask on string or binary data that has not run yet.

    Gathering string and binary values copies every byte of every selected
    value, which dominates merges and boolean-mask filtering of wide frames,
    while many of the gathered columns are never used afterwards. An
    ``ArrowExtensionArray`` built by ``take`` or a boolean mask therefore
    stores the parent array and the selection (the positions, or the mask)
    instead of ``_pa_array``, and gathers on first access to ``_pa_array``
    (see ``ArrowExtensionArray.__getattr__``). Further takes, slices and
    masks compose with the selection, so a chain of them runs as a single
    gather. Views share one instance, so the gather runs once for all of them
    and they share its buffers.

    A slice of a deferred take that has not run gathers only its own rows,
    so displaying a large result stays cheap, but it keeps a link to the
    take it was cut from (``origin``). Overlapping slices would otherwise
    gather the rows they share once each. So once the rows gathered through
    slices of a take would exceed its length, the take is gathered whole,
    and that slice and later ones become views of its result, as slices of
    an eager result are. A slice of such a slice links to the same take,
    so the links are never more than one deep and every slice counts
    against the same rows. Gathering costs about one copy per row, so this
    never gathers more than twice the rows an eager take would.

    Like a slice, which is a view of its parent's buffers, a deferred take
    keeps its parent alive until it is gathered. ``copy`` returns an array
    gathered into new buffers, which does not keep the parent alive, but
    leaves the array it was called on deferred.
    """

    __slots__ = (
        "length",
        "lock",
        "origin",
        "origin_start",
        "parent",
        "positions",
        "result",
        "selection",
        "sliced_rows",
    )

    def __init__(
        self,
        parent: pa.ChunkedArray,
        selection: npt.NDArray[np.intp] | npt.NDArray[np.bool_],
    ) -> None:
        self.parent: pa.ChunkedArray | None = parent
        self.selection: npt.NDArray[np.intp] | npt.NDArray[np.bool_] | None = selection
        # the positions of a mask's selected rows, computed when needed
        self.positions: npt.NDArray[np.intp] | None = None
        self.length = _selection_length(selection)
        self.result: pa.ChunkedArray | None = None
        self.lock = threading.Lock()
        # the deferred take this is a contiguous slice of, from origin_start
        self.origin: PendingTake | None = None
        self.origin_start = 0
        # rows gathered so far by slices of this take, on their own
        self.sliced_rows = 0

    def __len__(self) -> int:
        return self.length

    def slice_origin(self) -> tuple[PendingTake, int]:
        """
        The take a contiguous slice of this one is cut from, and where this
        one starts in it: its own origin if it is a slice of a take that has
        not been gathered, else itself.
        """
        with self.lock:
            if self.origin is not None:
                return self.origin, self.origin_start
            return self, 0

    def _slice_gathers_whole(self, length: int) -> bool:
        """
        Whether a slice of ``length`` rows should be cut from this take's
        result, gathering it first if needed, rather than gather its rows on
        its own. The slice's rows are counted if it gathers on its own.
        """
        with self.lock:
            if self.result is not None:
                return True
            if self.sliced_rows + length > self.length:
                return True
            self.sliced_rows += length
            return False

    def source(
        self, *, positions: bool = False
    ) -> tuple[pa.ChunkedArray, npt.NDArray[np.intp] | npt.NDArray[np.bool_] | None]:
        """
        The parent and the selection if the gather has not run, else the
        result and None. With ``positions``, a mask is returned as the
        positions it selects.
        """
        with self.lock:
            if self.result is not None:
                return self.result, None
            assert self.parent is not None
            assert self.selection is not None
            if positions and self.selection.dtype == np.bool_:
                if self.positions is None:
                    self.positions = np.flatnonzero(self.selection)
                return self.parent, self.positions
            return self.parent, self.selection

    def materialize(self) -> pa.ChunkedArray:
        with self.lock:
            if self.result is None:
                assert self.parent is not None
                assert self.selection is not None
                origin = self.origin
                # a slice only ever locks the take it was cut from, never
                # the other way round, so the locks cannot deadlock
                if origin is not None and origin._slice_gathers_whole(self.length):
                    self.result = origin.materialize().slice(
                        self.origin_start, self.length
                    )
                else:
                    self.result = _select(self.parent, self.selection)
                # release the parent, the selection and the origin
                self.parent = self.selection = self.positions = None
                self.origin = None
            return self.result

    def select(
        self, positions: npt.NDArray[np.integer] | npt.NDArray[np.bool_] | slice
    ) -> PendingTake | pa.ChunkedArray | None:
        """
        The rows of this take at ``positions``, composed with its selection
        so they are gathered straight from the parent: a ``PendingTake``, or
        the gathered values if there are none. None if this take has already
        been gathered, so the caller should select from its result.

        A mask is kept as a mask, to be gathered with ``filter``, while only
        masks are applied to it; a take or a slice of it composes with the
        positions it selects.
        """
        is_mask = isinstance(positions, np.ndarray) and positions.dtype == np.bool_
        parent, base = self.source(positions=not is_mask)
        if base is None:
            return None
        shared = None
        if isinstance(positions, np.ndarray) and is_immutable(positions):
            # positions that can never change: equal arrays compose equally
            shared = _shared.get()
        key = (id(base), id(positions))
        if shared is not None:
            entry = shared.get(key)
            if entry is not None:
                # another column composed this pair already: share its result
                return _pending_or_empty(parent, entry[2])
        composed: npt.NDArray[np.intp] | npt.NDArray[np.bool_]
        if base.dtype == np.bool_:
            # a mask of a mask: select among the rows ``base`` keeps
            assert isinstance(positions, np.ndarray)
            composed = base.copy()
            composed[base] = positions
        else:
            # indexing a NumPy array with an array or a slice: a new array
            # or a view of ``base``, which is never modified
            composed = base[positions]
        if shared is not None:
            # columns share it from now on. Nothing writes to a selection,
            # and only the deferred takes hold this new array, so the flag
            # keeps it unchanged without copying it into ``bytes``.
            composed.flags.writeable = False
            shared[key] = (base, positions, composed)
        return _pending_or_empty(parent, composed)

    def select_slice(self, item: slice) -> PendingTake | pa.ChunkedArray | None:
        """
        ``select(item)`` for a slice: this take itself for a full slice,
        which is a view; a contiguous slice keeps a link to the take it is
        cut from, so overlapping slices share one gather.
        """
        start, stop, step = item.indices(self.length)
        if (start, stop, step) == (0, self.length, 1):
            return self
        result = self.select(item)
        if isinstance(result, PendingTake) and step == 1:
            result.origin, offset = self.slice_origin()
            result.origin_start = offset + start
        return result

    def value_at(self, i: int | np.integer) -> pa.Scalar | None:
        """
        The value at position ``i``, read from the parent without gathering;
        None if this take has been gathered or ``i`` is out of bounds.
        """
        parent, base = self.source(positions=True)
        if base is None or not -len(base) <= i < len(base):
            return None
        return parent[int(base[i])]

    def gather_copy(self) -> pa.ChunkedArray | None:
        """
        The values gathered into new buffers, which do not keep the parent
        alive, leaving this take deferred; None if it has been gathered.
        """
        parent, base = self.source()
        if base is None:
            return None
        return _drop_unused_validity(_select(parent, base))


def freeze(values: npt.NDArray[_ScalarT]) -> npt.NDArray[_ScalarT]:
    """
    A copy of ``values`` in memory that can never change, so that deferred
    takes can keep it as it is (see ``is_immutable``).
    """
    return np.frombuffer(values.tobytes(), dtype=values.dtype)


def is_immutable(values: np.ndarray) -> bool:
    """
    Whether the memory of ``values`` can never change: it belongs to a
    ``bytes`` object, as that of an array made by ``freeze`` does. A
    read-only flag is not enough, since the owner of the memory can make it
    writeable again.
    """
    base: object = values
    while isinstance(base, np.ndarray):
        base = base.base
    return isinstance(base, bytes)


def _selection_length(
    selection: npt.NDArray[np.intp] | npt.NDArray[np.bool_],
) -> int:
    if selection.dtype == np.bool_:
        return int(np.count_nonzero(selection))
    return len(selection)


def defer(
    parent: pa.ChunkedArray,
    positions: npt.NDArray[np.integer] | npt.NDArray[np.bool_],
    *,
    copy: bool = True,
) -> PendingTake | pa.ChunkedArray:
    """
    The rows of ``parent`` at ``positions``, as a ``PendingTake``, or the
    gathered values if there are none.

    Integer ``positions`` must be in bounds and non-negative, and a boolean
    mask must have the length of ``parent``. Unless ``copy`` is False,
    either is copied before being stored, since the caller may modify its
    indexer afterwards.
    """
    # without ``copy``, copy only to convert the dtype
    copy_if = True if copy else None
    composed: npt.NDArray[np.intp] | npt.NDArray[np.bool_]
    if positions.dtype == np.bool_:
        composed = np.array(positions, dtype=np.bool_, copy=copy_if)
    else:
        composed = np.array(positions, dtype=np.intp, copy=copy_if)
    return _pending_or_empty(parent, composed)


def _pending_or_empty(
    parent: pa.ChunkedArray,
    selection: npt.NDArray[np.intp] | npt.NDArray[np.bool_],
) -> PendingTake | pa.ChunkedArray:
    if not _selection_length(selection):
        # nothing to defer, and no reason to keep the parent alive
        return _select(parent, selection)
    return PendingTake(parent, selection)
