"""Read-only browsing helpers for b2view."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import Any

import numpy as np

import blosc2

# Above this uncompressed size, plot_series does not read the whole series at
# once for an exact min/max envelope.  Local objects are instead streamed in
# bounded spans (still exact); only remote c2arrays fall back to a strided
# sample to avoid many network round-trips.  ~1 GB ≈ 125M float64.
_PLOT_FULL_READ_MAX_BYTES = 1_000_000_000

# Target size of a single streamed read in the exact-but-bounded envelope path.
_PLOT_STREAM_BUFFER_BYTES = 64_000_000  # ~64 MB


def _minmax_buckets(
    vmin: np.ndarray, vmax: np.ndarray, positions: np.ndarray, n: int, max_points: int
) -> dict[str, np.ndarray]:
    """Reduce per-element (or per-block) min/max into <= *max_points* buckets.

    *vmin*/*vmax* are the per-source-unit minima/maxima, *positions* the global
    row index of each unit's start, *n* the total row count.  NaN units (already
    NaN in *vmin*/*vmax*) are ignored within a bucket; an all-NaN bucket stays
    NaN.  Returns ``{"x", "ymin", "ymax"}`` with bucket-center x positions.
    """
    nunits = vmin.shape[0]
    if nunits == 0:
        empty = np.empty(0)
        return {"x": empty, "ymin": empty, "ymax": empty}
    if nunits <= max_points:
        starts = np.arange(nunits)
    else:
        group = -(-nunits // max_points)  # ceil
        starts = np.arange(0, nunits, group)
    # NaN-aware reduceat: +inf/-inf neutralizes NaN units, then mapped back.
    lo = np.where(np.isnan(vmin), np.inf, vmin)
    hi = np.where(np.isnan(vmax), -np.inf, vmax)
    ymin = np.minimum.reduceat(lo, starts)
    ymax = np.maximum.reduceat(hi, starts)
    ymin = np.where(np.isinf(ymin), np.nan, ymin)
    ymax = np.where(np.isinf(ymax), np.nan, ymax)
    x = np.minimum(positions[starts], max(0, n - 1))
    return {"x": x, "ymin": ymin, "ymax": ymax}


def _reduce_envelope(vals: np.ndarray, n: int, max_points: int) -> dict[str, np.ndarray]:
    """Per-bucket min/max envelope of an in-memory 1-D series."""
    vals = np.asarray(vals)
    if vals.shape[0] == 0:
        empty = np.empty(0)
        return {"x": empty, "ymin": empty, "ymax": empty}
    return _minmax_buckets(vals, vals, np.arange(vals.shape[0]), n, max_points)


def _bucket_geometry(n: int, max_points: int) -> tuple[int, int]:
    """Return ``(group, nbuckets)`` for an *n*-row series cut into ``<= max_points``
    contiguous buckets, matching the grouping policy of :func:`_minmax_buckets`
    so the streamed and full-read envelopes bucket rows identically."""
    if n <= 0:
        return 1, 0
    group = 1 if n <= max_points else -(-n // max_points)  # ceil
    nbuckets = -(-n // group)  # ceil
    return group, nbuckets


def _minmax_buckets_streaming(read_chunk, n: int, max_points: int, *, span: int) -> dict[str, np.ndarray]:
    """Exact per-bucket min/max envelope, read in row spans of at most *span*.

    Equivalent to reading the whole series and calling :func:`_reduce_envelope`,
    but never holds more than *span* rows in memory.  *read_chunk* is a callable
    ``(start, stop) -> 1-D array`` for that row range.  Because min/max are
    associative, arbitrary span boundaries (including ones that fall inside a
    bucket) yield a result identical to the single-read path.
    """
    group, nbuckets = _bucket_geometry(n, max_points)
    if nbuckets == 0:
        empty = np.empty(0)
        return {"x": empty, "ymin": empty, "ymax": empty}
    ymin = np.full(nbuckets, np.inf)
    ymax = np.full(nbuckets, -np.inf)
    span = max(1, int(span))
    for s in range(0, n, span):
        e = min(s + span, n)
        vals = np.asarray(read_chunk(s, e), dtype=float).ravel()[: e - s]
        bidx = np.arange(s, e) // group  # global bucket per row, non-decreasing
        seg_starts = np.concatenate(([0], np.flatnonzero(np.diff(bidx)) + 1))
        buckets = bidx[seg_starts]  # unique within this span (contiguous runs)
        lo = np.where(np.isnan(vals), np.inf, vals)
        hi = np.where(np.isnan(vals), -np.inf, vals)
        ymin[buckets] = np.minimum(ymin[buckets], np.minimum.reduceat(lo, seg_starts))
        ymax[buckets] = np.maximum(ymax[buckets], np.maximum.reduceat(hi, seg_starts))
    ymin = np.where(np.isinf(ymin), np.nan, ymin)
    ymax = np.where(np.isinf(ymax), np.nan, ymax)
    x = np.minimum(np.arange(nbuckets) * group, max(0, n - 1))
    return {"x": x, "ymin": ymin, "ymax": ymax}


def _stream_span(n: int, itemsize: int, chunklen: int | None) -> int:
    """Rows per streamed read: ~``_PLOT_STREAM_BUFFER_BYTES`` worth, aligned to
    whole native chunks when possible (chunks are the decompression unit)."""
    budget = max(1, _PLOT_STREAM_BUFFER_BYTES // max(1, itemsize))
    if chunklen and chunklen > 0:
        return chunklen if chunklen >= budget else (budget // chunklen) * chunklen
    return budget


@dataclass(frozen=True)
class NodeInfo:
    """Lightweight description of one TreeStore child."""

    path: str
    name: str
    kind: str
    has_children: bool


@dataclass(frozen=True)
class ObjectInfo:
    """Metadata for a TreeStore object or group."""

    path: str
    kind: str
    metadata: dict[str, Any]
    user_attrs: dict[str, Any] | None = None


@dataclass
class DataSliceLayout:
    """Describes the fixed/navigable state for slicing an N-D array into a 2-D table view.

    At most 2 dimensions can be navigable (shown as table rows/columns).
    All other dimensions must be fixed at a specific index value.
    """

    shape: tuple[int, ...]
    fixed_values: dict[int, int]  # dim_index → fixed index value
    navigable_dims: list[int]  # sorted list of up to 2 navigable dim indices

    # Current scroll positions for navigable dims
    # (index 0 → rows, index 1 → cols if present)
    row_start: int = 0
    row_stop: int = 0
    col_start: int = 0
    col_stop: int = 0

    # Optional locked window (absolute [start, stop)) on the navigable row dim.
    # When set, the grid sees a row dimension of length ``stop - start`` whose
    # logical row 0 maps to absolute row ``start`` (see ``preview_array_from_layout``).
    row_window: tuple[int, int] | None = None

    @classmethod
    def from_shape(cls, shape: tuple[int, ...]) -> DataSliceLayout:
        """Create a default layout: leading dims fixed at 0, last up-to-2 dims navigable."""
        ndim = len(shape)
        if ndim <= 2:
            navigable = list(range(ndim))
            fixed: dict[int, int] = {}
        else:
            navigable = list(range(ndim - 2, ndim))
            fixed = dict.fromkeys(range(ndim - 2), 0)
        return cls(
            shape=shape,
            fixed_values=fixed,
            navigable_dims=navigable,
        )

    def make_slices(self, max_rows: int = 20, max_cols: int = 10) -> tuple[int | slice, ...]:
        """Build the tuple of index expressions for slicing into the array.

        Uses *max_rows* and *max_cols* to size the navigable dimensions when
        ``row_stop <= row_start`` (i.e. no explicit stop was set).
        """
        slices: list[int | slice] = []
        for i in range(len(self.shape)):
            if i in self.fixed_values:
                slices.append(self.fixed_values[i])
            elif self.navigable_dims and i == self.navigable_dims[0]:
                start = max(0, min(self.row_start, self.shape[i]))
                if self.row_stop > self.row_start:
                    stop = min(self.row_stop, self.shape[i])
                else:
                    stop = min(start + max_rows, self.shape[i])
                slices.append(slice(start, stop))
            elif len(self.navigable_dims) > 1 and i == self.navigable_dims[1]:
                start = max(0, min(self.col_start, self.shape[i]))
                if self.col_stop > self.col_start:
                    stop = min(self.col_stop, self.shape[i])
                else:
                    stop = min(start + max_cols, self.shape[i])
                slices.append(slice(start, stop))
            else:
                slices.append(slice(0, self.shape[i]))
        return tuple(slices)

    def copy_with(
        self,
        *,
        fixed_values: dict[int, int] | None = None,
        navigable_dims: list[int] | None = None,
        row_start: int | None = None,
        row_stop: int | None = None,
        col_start: int | None = None,
        col_stop: int | None = None,
    ) -> DataSliceLayout:
        """Return a new layout with specified fields overridden."""
        return DataSliceLayout(
            shape=self.shape,
            fixed_values=self.fixed_values if fixed_values is None else fixed_values,
            navigable_dims=list(self.navigable_dims) if navigable_dims is None else navigable_dims,
            row_start=self.row_start if row_start is None else row_start,
            row_stop=self.row_stop if row_stop is None else row_stop,
            col_start=self.col_start if col_start is None else col_start,
            col_stop=self.col_stop if col_stop is None else col_stop,
            row_window=self.row_window,
        )

    def row_window_bounds(self, row_dim: int | None) -> tuple[int, int]:
        """Return the absolute [start, stop) extent of the navigable row dim.

        Narrowed to ``row_window`` when one is set; otherwise the full dim.
        """
        full = self.shape[row_dim] if row_dim is not None else 1
        if row_dim is None or self.row_window is None:
            return 0, full
        w0, w1 = self.row_window
        w0 = max(0, min(w0, full))
        w1 = max(w0, min(w1, full))
        return w0, w1

    def total_for_dim(self, dim: int) -> int:
        """Return the total size of *dim*."""
        if 0 <= dim < len(self.shape):
            return self.shape[dim]
        return 0


class StoreBrowser:
    """Small, read-only adapter used by the b2view UI.

    The adapter intentionally exposes a narrow API so the TUI does not depend
    on TreeStore internals.  It accepts either a TreeStore hierarchy or a
    single top-level Blosc2 object (for example a standalone CTable).  It
    performs bounded previews only; callers must explicitly request pages or
    slices.
    """

    def __init__(self, urlpath: str):
        self.urlpath = urlpath
        self.store = blosc2.open(urlpath, mode="r")
        self.is_tree = isinstance(self.store, blosc2.TreeStore)
        # Per-path row filters for CTable nodes (path -> expr / where() view)
        self._filters: dict[str, str] = {}
        self._filter_views: dict[str, Any] = {}
        # Per-path locked row windows for CTable nodes (path -> slice() view)
        self._window_views: dict[str, Any] = {}
        # Per-path sort order for CTable nodes (path -> (column, reverse) / view)
        self._sorts: dict[str, tuple[str, bool]] = {}
        self._sort_views: dict[str, Any] = {}
        # Per-path group-by for CTable nodes (path -> (key, op, value_col|None) /
        # materialized result CTable).  Exclusive with filter/window/column-sel.
        self._groups: dict[str, tuple[str, str, str | None]] = {}
        self._group_views: dict[str, Any] = {}
        # Optional sort applied on top of a grouped result (path -> (col, reverse)
        # / sorted view of the small result).  Lives only while grouped.
        self._group_sorts: dict[str, tuple[str, bool]] = {}
        self._group_sort_views: dict[str, Any] = {}
        # Memoized grouped results, keyed by (path, key, op, value_col).  Survives
        # ungrouping so re-running a prior group-by is instant (store is read-only).
        self._group_result_cache: dict[tuple[str, str, str, str | None], Any] = {}
        # Per-path column filters (path -> substring pattern / matched names)
        self._column_filters: dict[str, str] = {}
        self._column_selections: dict[str, list[str]] = {}

    def close(self) -> None:
        close = getattr(self.store, "close", None)
        if close is not None:
            close()

    def __enter__(self) -> StoreBrowser:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    @staticmethod
    def normalize_path(path: str) -> str:
        """Return an absolute TreeStore path."""
        if not path:
            return "/"
        if not path.startswith("/"):
            path = "/" + path
        normalized = str(PurePosixPath(path))
        return "/" if normalized == "." else normalized

    def list_children(self, path: str = "/") -> list[NodeInfo]:
        """Return direct children for *path*."""
        path = self.normalize_path(path)
        if not self.is_tree:
            self._check_root_path(path)
            return []

        children = []
        for child_path in self.store.get_children(path):
            descendants = self.store.get_descendants(child_path)
            has_children = bool(descendants)
            kind = "group" if has_children else self.kind(child_path)
            children.append(
                NodeInfo(
                    path=child_path,
                    name=child_path.rsplit("/", 1)[-1] or "/",
                    kind=kind,
                    has_children=has_children,
                )
            )
        return children

    def kind(self, path: str) -> str:
        """Classify a browser path."""
        path = self.normalize_path(path)
        if not self.is_tree:
            self._check_root_path(path)
            return object_kind(self.store)
        if path == "/" or self.store.get_descendants(path):
            return "group"
        obj = self.store[path]
        return object_kind(obj)

    def get_info(self, path: str) -> ObjectInfo:
        """Return metadata for *path*."""
        path = self.normalize_path(path)
        kind = self.kind(path)
        if kind == "group":
            metadata: dict[str, Any] = {
                "type": "TreeStore group",
                "children": len(self.store.get_children(path)),
                "descendants": len(self.store.get_descendants(path)),
            }
            user_attrs = self._vlmeta_dict(self.store.vlmeta)
            return ObjectInfo(path=path, kind=kind, metadata=metadata, user_attrs=user_attrs)

        obj = self._get_object(path)
        metadata = object_metadata(obj)
        metadata.setdefault("type", type(obj).__name__)
        user_attrs = self._vlmeta_dict(getattr(obj, "vlmeta", None))
        if user_attrs is None and self.is_tree:
            user_attrs = self._vlmeta_dict(self.store.vlmeta)
        return ObjectInfo(path=path, kind=kind, metadata=metadata, user_attrs=user_attrs)

    def preview(
        self,
        path: str,
        *,
        start: int = 0,
        stop: int | None = None,
        columns: list[str] | None = None,
        slices: tuple[Any, ...] | None = None,
        max_rows: int = 20,
        max_cols: int = 10,
        col_start: int = 0,
        slice_indices: list[int] | None = None,
        layout: DataSliceLayout | None = None,
    ) -> Any:
        """Return a bounded data preview for *path*.

        For N-D arrays (N >= 3) a *layout* may be provided instead of the
        legacy *slice_indices*, *start*/*stop*, *col_start* parameters.
        """
        path = self.normalize_path(path)
        obj = self._get_object(path)
        kind = object_kind(obj)
        if kind in {"ndarray", "c2array"}:
            shape = tuple(getattr(obj, "shape", ()) or ())
            if slices is None:
                if layout is not None:
                    return preview_array_from_layout(
                        obj, layout=layout, max_rows=max_rows, max_cols=max_cols
                    )
                if len(shape) >= 3:
                    return preview_array_nd_slice(
                        obj,
                        slice_indices=slice_indices,
                        start=start,
                        stop=stop,
                        col_start=col_start,
                        max_cols=max_cols,
                    )
                if len(shape) == 2:
                    stop = min(start + max_rows, shape[0]) if stop is None else stop
                    return preview_array_2d(
                        obj, start=start, stop=stop, col_start=col_start, max_cols=max_cols
                    )
                if len(shape) == 1:
                    stop = min(start + max_rows, shape[0]) if stop is None else stop
                    return preview_array_1d(obj, start=start, stop=stop)
            return preview_array(obj, slices=slices, max_rows=max_rows, max_cols=max_cols)
        if kind == "ctable":
            # Read precedence: a group-by (smaller aggregated result) wins, then
            # a locked row window (set by 'v'), then a row filter, then a sort
            # view; all already fold in their predecessors.
            obj = self._ordered_object(path, obj)
            if columns is None:
                columns = self._column_selections.get(path)
            stop = min(start + max_rows, len(obj)) if stop is None else stop
            return preview_ctable(
                obj, start=start, stop=stop, columns=columns, max_cols=max_cols, col_start=col_start
            )
        if kind == "schunk":
            stop = start + max_rows if stop is None else stop
            return preview_schunk(obj, start=start, stop=stop)
        return {"message": f"Preview is not supported for {kind!r} objects."}

    def plot_series(
        self,
        path: str,
        *,
        column: str | int | None = None,
        layout: DataSliceLayout | None = None,
        max_points: int = 2000,
        row_start: int = 0,
        row_stop: int | None = None,
    ) -> dict[str, Any]:
        """Return a peak-preserving overview of one series for plotting.

        The result is ``{"x", "ymin", "ymax", "n", "row_start", "row_stop",
        "method"}`` with at most *max_points* buckets; ``ymin``/``ymax`` are the
        per-bucket extremes so a plotted envelope never hides a peak or trough,
        ``n`` is the *total* series length and ``row_start``/``row_stop`` the
        plotted range.  Three tiers, cheapest first:

        - ``"summary"``: read precomputed per-block min/max from the column's
          SUMMARY index — no data decompression (whole series only; CTable
          columns, no active filter, numeric).
        - ``"reduce"``: read the (sub)series and reduce per bucket — exact,
          O(range), streamed in bounded spans above ``_PLOT_FULL_READ_MAX_BYTES``
          for local objects.
        - ``"sample"``: strided sample for a remote series too large to read
          fully; this may miss extremes, so callers should label it.

        Pass *row_start*/*row_stop* to zoom into a sub-range (always read
        exactly; ``x`` stays in absolute row coordinates).  The series is a
        CTable column (*column* is its name; a locked row window takes
        precedence, otherwise an active row filter is honored) or an array
        (*column* is the global index along the column dimension of *layout*,
        or None for 1-D arrays).
        """
        path = self.normalize_path(path)
        obj = self._get_object(path)
        kind = object_kind(obj)

        if kind == "ctable":
            # Window/filter/sort precedence mirrors preview()/read_cell(): a plot
            # shows exactly the rows (and order) the grid is showing.  The SUMMARY
            # fast-path spans the whole column in original order, so it is only
            # valid when nothing narrows *or reorders* the series.
            view = self._ordered_object(path, obj)
            narrowed = view is not obj
            n = len(view)
            start, stop = self._clamp_range(row_start, row_stop, n)
            if start == 0 and stop == n and not narrowed:
                env = self._column_summary_envelope(obj, column, n, max_points)
                if env is not None:
                    return {**env, "n": n, "row_start": start, "row_stop": stop, "method": "summary"}
            # Range ordered purely by the very column we're plotting: the sort
            # view is monotonic over any contiguous range, so we can read just
            # the bucket boundaries instead of gathering every value (this also
            # accelerates zooming, where start/stop is a sub-range).
            sort = self._sorts.get(path)
            sorted_only = path not in self._window_views and path not in self._filter_views
            if sorted_only and sort is not None and sort[0] == column:
                env = self._sorted_column_envelope(view[column], start, stop, n, max_points)
                return {**env, "n": n, "row_start": start, "row_stop": stop, "method": "sorted"}
            col = view[column]
            chunks = getattr(col, "chunks", None)
            return self._range_envelope(
                lambda s, e, st=1: safe_asarray(col[s:e:st]),
                start,
                stop,
                n,
                np.dtype(col.dtype).itemsize,
                chunks[0] if chunks else None,
                remote=False,
                max_points=max_points,
            )

        if kind in {"ndarray", "c2array"}:
            shape = tuple(getattr(obj, "shape", ()) or ())
            ndim = len(shape)
            if ndim == 0:
                raise ValueError("Cannot plot a scalar")
            row_dim = layout.navigable_dims[0] if layout is not None and layout.navigable_dims else 0
            n = shape[row_dim]
            start, stop = self._clamp_range(row_start, row_stop, n)

            def _row_index(row_slice):
                idx: list[int | slice] = []
                for i in range(ndim):
                    if i == row_dim:
                        idx.append(row_slice)
                    elif layout is not None and i in layout.fixed_values:
                        idx.append(layout.fixed_values[i])
                    elif (
                        layout is not None
                        and len(layout.navigable_dims) > 1
                        and i == layout.navigable_dims[1]
                    ):
                        idx.append(int(column))
                    else:
                        idx.append(0)
                return tuple(idx)

            chunks = getattr(obj, "chunks", None)
            return self._range_envelope(
                lambda s, e, st=1: np.asarray(obj[_row_index(slice(s, e, st))]),
                start,
                stop,
                n,
                np.dtype(obj.dtype).itemsize,
                chunks[row_dim] if chunks else None,
                remote=(kind == "c2array"),
                max_points=max_points,
            )

        raise ValueError(f"Cannot plot {kind!r} objects")

    def read_series(
        self,
        path: str,
        *,
        column: str | int | None = None,
        layout: DataSliceLayout | None = None,
        row_start: int = 0,
        row_stop: int | None = None,
        max_points: int | None = None,
    ) -> dict[str, Any]:
        """Return the *raw* values of one series over ``[row_start, row_stop)``.

        Same series selection as :meth:`plot_series` (CTable column honoring a
        locked row window then an active filter, or an array column via
        *layout*) but with no bucketing —
        every value is read exactly, for the high-res ``h``/``r`` view.  The
        result is ``{"x", "y", "n", "row_start", "row_stop", "stride", "shown",
        "sampled"}`` with ``x`` in absolute row coordinates.

        When *max_points* is given and the range is wider, the read is
        strided-sampled (``stride = ceil(width / max_points)``, like
        :meth:`read_xy`) so a wide raw range stays bounded; otherwise it is read
        exactly (``stride=1``, ``sampled=False``).
        """
        path = self.normalize_path(path)
        obj = self._get_object(path)
        kind = object_kind(obj)

        if kind == "ctable":
            # Window > filter > sort, matching preview()/read_cell() so the
            # hi-res view tracks the visible grid (rows and order).
            view = self._ordered_object(path, obj)
            n = len(view)
            start, stop = self._clamp_range(row_start, row_stop, n)
            stride = self._series_stride(stop - start, max_points)
            y = safe_asarray(view[column][start:stop:stride])
        elif kind in {"ndarray", "c2array"}:
            shape = tuple(getattr(obj, "shape", ()) or ())
            ndim = len(shape)
            if ndim == 0:
                raise ValueError("Cannot plot a scalar")
            row_dim = layout.navigable_dims[0] if layout is not None and layout.navigable_dims else 0
            n = shape[row_dim]
            start, stop = self._clamp_range(row_start, row_stop, n)
            stride = self._series_stride(stop - start, max_points)
            # Same column/fixed-dim selection as plot_series' array branch.
            idx: list[int | slice] = []
            for i in range(ndim):
                if i == row_dim:
                    idx.append(slice(start, stop, stride))
                elif layout is not None and i in layout.fixed_values:
                    idx.append(layout.fixed_values[i])
                elif layout is not None and len(layout.navigable_dims) > 1 and i == layout.navigable_dims[1]:
                    idx.append(int(column))
                else:
                    idx.append(0)
            y = np.asarray(obj[tuple(idx)])
        else:
            raise ValueError(f"Cannot plot {kind!r} objects")

        return {
            "x": np.arange(start, stop, stride),
            "y": y,
            "n": n,
            "row_start": start,
            "row_stop": stop,
            "stride": stride,
            "shown": len(y),
            "sampled": stride > 1,
        }

    @staticmethod
    def _series_stride(width: int, max_points: int | None) -> int:
        """Stride to keep a raw read within *max_points* (1 == exact)."""
        if max_points is None or width <= max_points:
            return 1
        return max(1, -(-width // max_points))

    def read_xy(
        self,
        path: str,
        *,
        xcol: str,
        ycol: str,
        layout: DataSliceLayout | None = None,
        row_start: int = 0,
        row_stop: int | None = None,
        max_points: int = 50_000,
    ) -> dict[str, Any]:
        """Return two row-aligned CTable columns over ``[row_start, row_stop)``.

        For the col-vs-col scatter (``s`` in the plot panel): *xcol* and *ycol*
        are read over the **same** live-row range, using the same window→filter
        precedence as :meth:`read_series`, so the points are row-aligned for
        free.  When the range is wider than *max_points*, both columns are
        strided-sampled (same stride) and ``sampled`` is set.  Both columns must
        be numeric; otherwise a ``ValueError`` is raised.  Result keys:
        ``{"x", "y", "n", "row_start", "row_stop", "stride", "shown",
        "sampled"}``.
        """
        path = self.normalize_path(path)
        obj = self._get_object(path)
        kind = object_kind(obj)
        if kind != "ctable":
            raise ValueError("Scatter requires a CTable source")

        # Window > filter > sort, matching read_series() so the scatter tracks
        # exactly the visible rows (and order).
        view = self._ordered_object(path, obj)
        n = len(view)
        start, stop = self._clamp_range(row_start, row_stop, n)
        width = stop - start
        stride = max(1, -(-width // max_points)) if width > max_points else 1

        x = safe_asarray(view[xcol][start:stop:stride])
        y = safe_asarray(view[ycol][start:stop:stride])
        for nm, arr in ((xcol, x), (ycol, y)):
            if arr.dtype.kind not in "iufb":
                raise ValueError(f"Column {nm!r} is not numeric")

        return {
            "x": x,
            "y": y,
            "n": n,
            "row_start": start,
            "row_stop": stop,
            "stride": stride,
            "shown": len(x),
            "sampled": stride > 1,
        }

    def read_cell(self, path: str, column: str, row: int) -> Any:
        """Decode a single CTable cell — the on-demand path for expensive columns.

        *row* is in the same live-row space as :meth:`preview` (it mirrors the
        window/filter view precedence), so the row the grid shows is the cell
        that gets decoded.  Returns the native Python value (list/dict/array/…),
        not a NumPy-wrapped one, so callers can pretty-print its structure.
        """
        path = self.normalize_path(path)
        obj = self._get_object(path)
        if object_kind(obj) == "ctable":
            # Same precedence as preview() (window > filter > sort) so the
            # visible row index resolves the same cell.
            obj = self._ordered_object(path, obj)
        values = obj[column][row : row + 1]
        if len(values) == 0:
            raise IndexError(f"row {row} is out of range")
        return values[0]

    @staticmethod
    def _clamp_range(row_start: int, row_stop: int | None, n: int) -> tuple[int, int]:
        start = 0 if row_start is None else max(0, min(int(row_start), n))
        stop = n if row_stop is None else max(0, min(int(row_stop), n))
        return (stop, start) if stop < start else (start, stop)

    def _range_envelope(
        self,
        read,
        start: int,
        stop: int,
        n_total: int,
        itemsize: int,
        chunklen: int | None,
        *,
        remote: bool,
        max_points: int,
    ) -> dict[str, Any]:
        """Envelope of rows ``[start, stop)`` via *read(s, e, step=1)``, with ``x``
        in absolute row coordinates.  Reads the range exactly (reduce/stream),
        falling back to a strided sample only for large *remote* ranges."""
        rng = stop - start
        base = {"n": n_total, "row_start": start, "row_stop": stop}
        if rng <= 0:
            empty = np.empty(0)
            return {"x": empty, "ymin": empty, "ymax": empty, **base, "method": "reduce"}
        if rng * itemsize > _PLOT_FULL_READ_MAX_BYTES:
            if remote:
                step = max(1, -(-rng // max_points))
                y = np.asarray(read(start, stop, step))
                x = np.arange(start, stop, step)
                m = min(len(x), len(y))
                return {"x": x[:m], "ymin": y[:m], "ymax": y[:m], **base, "method": "sample"}
            span = _stream_span(rng, itemsize, chunklen)
            env = _minmax_buckets_streaming(
                lambda s, e: read(start + s, start + e), rng, max_points, span=span
            )
        else:
            env = _reduce_envelope(np.asarray(read(start, stop)), rng, max_points)
        env["x"] = np.asarray(env["x"]) + start
        return {**env, **base, "method": "reduce"}

    def _sorted_column_envelope(
        self, col: Any, start: int, stop: int, n: int, max_points: int
    ) -> dict[str, np.ndarray]:
        """Exact min/max envelope of rows ``[start, stop)`` of a column read
        through its own sort view, with ``x`` in absolute row coordinates.

        When the plotted column *is* the column the view is sorted by, the values
        are monotonic over any contiguous range (NaNs land in a contiguous block
        at the very end), so each bucket's min/max are just its two endpoints —
        no need to gather every value.  We read only the ~2*nbuckets boundary
        values plus, for the lone finite/NaN transition bucket, that one bucket
        in full.  Bit-identical to the full-read reduce path, but ~50x cheaper.

        ponytail: relies on the sort view being monotonic; only call it when the
        plotted column equals the sort column (see :meth:`plot_series`).
        """
        rng = stop - start
        group, nbuckets = _bucket_geometry(rng, max_points)
        if nbuckets == 0:
            empty = np.empty(0)
            return {"x": empty, "ymin": empty, "ymax": empty}
        offs = np.arange(nbuckets) * group  # bucket starts, relative to *start*
        ends = start + np.minimum(offs + group - 1, rng - 1)
        vstart = np.asarray(col[start:stop:group], dtype=float)[:nbuckets]
        vend = np.asarray(col[ends], dtype=float)
        ymin = np.fmin(vstart, vend)
        ymax = np.fmax(vstart, vend)
        # A bucket straddling the finite/NaN boundary has one NaN endpoint; its
        # interior extreme is hidden, so read that one bucket exactly.
        for i in np.flatnonzero(np.isnan(vstart) != np.isnan(vend)):
            s = start + int(offs[i])
            seg = np.asarray(col[s : min(s + group, stop)], dtype=float)
            ymin[i] = np.nanmin(seg)
            ymax[i] = np.nanmax(seg)
        x = start + np.minimum(offs, max(0, rng - 1))
        return {"x": x, "ymin": ymin, "ymax": ymax}

    def _column_summary_envelope(
        self, table: Any, column: str | int | None, n: int, max_points: int
    ) -> dict[str, np.ndarray] | None:
        """Build a min/max envelope from a column's index summaries, or None.

        Reads precomputed per-block ``(min, max)`` from the index — no data
        decompression.  Every index kind (SUMMARY, FULL, PARTIAL, BUCKET, OPSI)
        persists the same block-level ``(min, max, flags)`` sidecars in its
        ``levels`` descriptor, so any indexed numeric column plots instantly
        without a dedicated summary index.  Returns None when there is no usable
        summary (non-string column, no index, non-numeric, or no block level).
        """
        if not isinstance(column, str):
            return None
        try:
            idx = table.index(column)
        except Exception:
            return None
        try:
            desc = idx.descriptor
            levels = desc.get("levels") or {}
            # Prefer the finest whole-column level available (block), else any.
            level = "block" if "block" in levels else next(iter(levels), None)
            if level is None or np.dtype(desc["dtype"]).kind not in "iuf":
                return None
            path = levels[level].get("path")
            if path is None:
                return None  # in-memory sidecar: nothing to fast-read
            from blosc2.indexing import _INDEX_MMAP_MODE, FLAG_ALL_NAN, _open_sidecar_file

            # Drop the handle after reading: the cached _open_level_summary_handle
            # would hold a file descriptor open for the whole session (one per
            # plotted column), exhausting the FD limit on large test runs.
            handle = _open_sidecar_file(path, _INDEX_MMAP_MODE)
            bmin = np.asarray(handle["min"][:])
            bmax = np.asarray(handle["max"][:])
            flags = np.asarray(handle["flags"][:])
            del handle
        except Exception:
            return None
        if bmin.shape[0] == 0:
            return None
        all_nan = (flags & FLAG_ALL_NAN) != 0
        if all_nan.any():
            bmin = np.where(all_nan, np.nan, bmin)
            bmax = np.where(all_nan, np.nan, bmax)
        block = int(desc["blocks"][0])
        positions = np.arange(bmin.shape[0]) * block
        return _minmax_buckets(bmin, bmax, positions, n, max_points)

    def column_names(self, path: str) -> list[str] | None:
        """Return the column names for a CTable path, or None for other kinds.

        When a column filter is active, only the matching names are returned
        (navigation operates on the filtered universe).
        """
        path = self.normalize_path(path)
        if path in self._group_views:
            return list(getattr(self._group_views[path], "col_names", []) or []) or None
        selection = self._column_selections.get(path)
        if selection is not None:
            return list(selection)
        names = list(getattr(self._get_object(path), "col_names", []) or [])
        return names or None

    def set_filter(self, path: str, expr: str | None) -> int:
        """Set or clear the row filter of a CTable path; return its row count.

        An empty (or None) *expr* clears the filter.  Errors from ``where()``
        propagate to the caller and leave any previous filter untouched.
        """
        path = self.normalize_path(path)
        expr = (expr or "").strip()
        if not expr:
            self._filters.pop(path, None)
            self._filter_views.pop(path, None)
            return len(self._get_object(path))
        # A filter redefines the row set, so any existing sort over the old rows
        # is dropped; the user re-sorts the filtered rows on top if they want.
        self.clear_sort(path)
        view = self._get_object(path).where(expr)
        self._filters[path] = expr
        self._filter_views[path] = view
        return len(view)

    def get_filter(self, path: str) -> str | None:
        """Return the active filter expression for *path*, if any."""
        return self._filters.get(self.normalize_path(path))

    def full_index_columns(self, path: str) -> list[str]:
        """Names of CTable columns at *path* that carry a FULL index (sortable)."""
        obj = self._get_object(path)
        return [ix.col_name for ix in getattr(obj, "indexes", []) if getattr(ix, "kind", None) == "full"]

    def _row_source(self, path: str) -> Any:
        """Base rows that sort/group build on: the active row filter view if any,
        else the underlying table.  This is what makes a filter compose — sort
        and group operate on the filtered rows, not the whole table."""
        return self._filter_views.get(path, self._get_object(path))

    def set_sort(self, path: str, column: str, reverse: bool) -> None:
        """Order the CTable at *path* by *column* as a zero-copy sorted view.

        Composes over any active row filter (sorts the filtered rows); replaces
        any locked window.  With no filter and a FULL index on *column*, the view
        streams from the index, so the full table is never materialised.
        """
        path = self.normalize_path(path)
        self._window_views.pop(path, None)
        view = self._row_source(path).sort_by(column, ascending=not reverse, view=True)
        self._sorts[path] = (column, reverse)
        self._sort_views[path] = view

    def clear_sort(self, path: str) -> None:
        """Drop any sort order from *path*, restoring the original row order."""
        path = self.normalize_path(path)
        self._sorts.pop(path, None)
        self._sort_views.pop(path, None)

    def get_sort(self, path: str) -> tuple[str, bool] | None:
        """Return the active ``(column, reverse)`` sort for *path*, if any."""
        return self._sorts.get(self.normalize_path(path))

    def _ordered_object(self, path: str, obj: Any) -> Any:
        """CTable read precedence for *path*: group > window > sort > filter > base.

        Sort and group build on the filtered rows (see :meth:`_row_source`), so a
        sort view already incorporates the filter and ranks above the bare filter
        view; the bare filter view is read only when no sort is active.
        """
        if path in self._group_sort_views:
            return self._group_sort_views[path]
        if path in self._group_views:
            return self._group_views[path]
        if path in self._window_views:
            return self._window_views[path]
        if path in self._sort_views:
            return self._sort_views[path]
        return self._filter_views.get(path, obj)

    def set_row_window(self, path: str, start: int, stop: int) -> int:
        """Lock the CTable at *path* to live rows ``[start:stop]``; return its length.

        The window is a zero-copy :meth:`CTable.slice` view of whatever is
        currently visible (so it composes over any active row filter).  Paging
        then cannot leave the range because the view reports only its own rows.
        """
        path = self.normalize_path(path)
        # The sort view already incorporates any filter, so prefer it; fall back
        # to the bare filter view, then the base table.
        if path in self._sort_views:
            base = self._sort_views[path]
        elif path in self._filter_views:
            base = self._filter_views[path]
        else:
            base = self._get_object(path)
        view = base.slice(start, stop, copy=False)
        self._window_views[path] = view
        return len(view)

    def clear_row_window(self, path: str) -> None:
        """Remove any locked row window from *path*."""
        self._window_views.pop(self.normalize_path(path), None)

    def get_row_window(self, path: str) -> bool:
        """Return whether *path* currently has a locked row window."""
        return self.normalize_path(path) in self._window_views

    def group_key_columns(self, path: str) -> list[str]:
        """CTable columns at *path* usable as group-by keys (dictionary or numeric).

        Includes floats: grouping by a numeric column (e.g. trip duration/distance
        or a timestamp) buckets rows by exact value — handy for spotting rush
        hours or best-profit windows.  Cardinality is the user's call; we don't
        bin.  Non-dictionary text and nested/ndarray columns are excluded.
        """
        obj = self._get_object(path)
        out = []
        for name in getattr(obj, "col_names", []) or []:
            col = obj[name]
            dt = getattr(col, "dtype", None)
            if getattr(col, "is_dictionary", False) or (dt is not None and dt.kind in "iuf"):
                out.append(name)
        return out

    def group_value_columns(self, path: str) -> list[str]:
        """CTable columns at *path* usable as aggregation value columns (numeric scalar)."""
        obj = self._get_object(path)
        out = []
        for name in getattr(obj, "col_names", []) or []:
            col = obj[name]
            if getattr(col, "is_dictionary", False) or getattr(col, "is_ndarray", False):
                continue
            dt = getattr(col, "dtype", None)
            if dt is not None and dt.kind in "iuf":
                out.append(name)
        return out

    def set_group(self, path: str, key: str, op: str, value_col: str | None) -> int:
        """Group the CTable at *path* by *key*, aggregating with *op*; return group count.

        Builds a materialized result CTable (one row per group, columns = key +
        aggregate).  Composes over any active row filter (groups the filtered
        rows); drops any window, sort, and column selection.  Errors from
        ``group_by`` propagate.
        """
        path = self.normalize_path(path)
        # Memoize the materialized result: the store is read-only, so a given
        # (path, filter, key, op, value_col) always aggregates to the same tiny
        # CTable.  The active filter expr is part of the key so a filtered group
        # never collides with the unfiltered one.
        # ponytail: unbounded dict, but results are one-row-per-group and a
        # session tries only a handful of configs — add an LRU cap if that ever
        # stops being true.
        cache_key = (path, self._filters.get(path), key, op, value_col)
        result = self._group_result_cache.get(cache_key)
        if result is None:
            gb = self._row_source(path).group_by(key)
            result = gb.size() if op == "size" else gb.agg({value_col: op})
            self._group_result_cache[cache_key] = result
        self._window_views.pop(path, None)
        self._sorts.pop(path, None)
        self._sort_views.pop(path, None)
        self._column_filters.pop(path, None)
        self._column_selections.pop(path, None)
        self._group_sorts.pop(path, None)  # a fresh result invalidates any group sort
        self._group_sort_views.pop(path, None)
        self._groups[path] = (key, op, value_col)
        self._group_views[path] = result
        return len(result)

    def get_group(self, path: str) -> tuple[str, str, str | None] | None:
        """Return the active ``(key, op, value_col)`` group-by for *path*, if any."""
        return self._groups.get(self.normalize_path(path))

    def clear_group(self, path: str) -> None:
        """Drop any group-by (and its sort) from *path*, restoring the base table."""
        path = self.normalize_path(path)
        self._groups.pop(path, None)
        self._group_views.pop(path, None)
        self._group_sorts.pop(path, None)
        self._group_sort_views.pop(path, None)

    def set_group_sort(self, path: str, column: str, reverse: bool) -> None:
        """Sort the (already grouped) result at *path* by one of its columns.

        The grouped result is tiny and carries no index, so this lexsorts it as a
        zero-copy view.  No-op if *path* is not grouped.
        """
        path = self.normalize_path(path)
        result = self._group_views.get(path)
        if result is None:
            return
        self._group_sorts[path] = (column, reverse)
        self._group_sort_views[path] = result.sort_by(column, ascending=not reverse, view=True)

    def get_group_sort(self, path: str) -> tuple[str, bool] | None:
        """Return the active ``(column, reverse)`` sort on the grouped result, if any."""
        return self._group_sorts.get(self.normalize_path(path))

    def clear_group_sort(self, path: str) -> None:
        """Drop any sort on the grouped result, keeping the group itself."""
        path = self.normalize_path(path)
        self._group_sorts.pop(path, None)
        self._group_sort_views.pop(path, None)

    def group_agg_column(self, path: str) -> str | None:
        """Name of the aggregate column in the grouped result for *path*, if grouped."""
        group = self._groups.get(self.normalize_path(path))
        if group is None:
            return None
        key, op, value_col = group
        return "size" if op == "size" else f"{value_col}_{op}"

    def group_bars(self, path: str, top_n: int = 50) -> dict[str, Any]:
        """Grouped result as plot data, chosen by key dtype.

        Follows the grid's active sort on the grouped result (so the plot matches
        what the user sees — Pareto when sorted by the aggregate, key-order when
        sorted by the key); with no sort, ranks by the aggregate descending.

        A **categorical** key (dictionary/string) yields a *bar* chart: the first
        ``top_n`` groups as ``labels`` + ``values``.  A **numeric** key yields a
        *line* curve over **all** groups (``x`` + ``values``): ``x`` is the key
        value when sorted by the key (a spacing-honest distribution), else the
        rank index (a Pareto curve).  ``numeric`` says which.
        """
        path = self.normalize_path(path)
        result = self._group_views.get(path)
        group = self._groups.get(path)
        empty = {
            "numeric": False,
            "labels": [],
            "x": [],
            "values": [],
            "total": 0,
            "key": "",
            "agg": "",
            "xlabel": "",
        }
        if result is None or group is None:
            return empty
        key = group[0]
        agg = self.group_agg_column(path)
        total = len(result)
        sorted_result = self._group_sort_views.get(path)
        sort = self._group_sorts.get(path)
        source = sorted_result if sorted_result is not None else result

        if self._is_numeric_key(result, key):
            # Full curve, no top-N cap.  Honour the active sort; default to
            # aggregate-descending (a Pareto curve) when none is set.
            keys = np.asarray(source[key][:], dtype=float)
            vals = np.asarray(source[agg][:], dtype=float)
            if sorted_result is None:
                order = np.argsort(vals)[::-1]
                keys, vals = keys[order], vals[order]
            if sort is not None and sort[0] == key:  # sorted by the key → key on X
                x, xlabel = keys.tolist(), key
            else:  # ranked by the aggregate → rank on X (Pareto)
                x, xlabel = list(range(len(vals))), f"rank (by {agg})"
            return {
                "numeric": True,
                "labels": [],
                "x": [float(v) for v in x],
                "values": [float(v) for v in vals],
                "total": total,
                "key": key,
                "agg": agg,
                "xlabel": xlabel,
            }

        # Categorical key → top-N bars.
        if sorted_result is not None:
            keys = safe_asarray(sorted_result[key][:top_n])
            vals = np.asarray(sorted_result[agg][:top_n], dtype=float)
            order = range(len(vals))
        else:
            keys = safe_asarray(result[key][:])
            vals = np.asarray(result[agg][:], dtype=float)
            order = np.argsort(vals)[::-1][:top_n]
        return {
            "numeric": False,
            "x": [],
            "labels": [str(keys[i]) for i in order],
            "values": [float(vals[i]) for i in order],
            "total": total,
            "key": key,
            "agg": agg,
            "xlabel": key,
        }

    @staticmethod
    def _is_numeric_key(result: Any, key: str) -> bool:
        """True when the grouped key column is numeric (not a decoded dictionary)."""
        col = result[key]
        if getattr(col, "is_dictionary", False):
            return False
        dt = getattr(col, "dtype", None)
        return dt is not None and dt.kind in "iuf"

    def base_nrows(self, path: str) -> int:
        """Return the unfiltered row count of the CTable at *path*."""
        return len(self._get_object(path))

    def set_column_selection(self, path: str, names: list[str] | None) -> int:
        """Restrict a CTable path to an explicit ordered set of *names*.

        Keeps an arbitrary chosen subset, in table order.  Unknown names are
        dropped; an empty (or None) selection — or one naming every column —
        clears the filter (all columns visible).  Returns the number of columns
        now visible.  A status-chip descriptor is stored in ``_column_filters``
        so :meth:`get_column_filter` stays truthy while a selection is active.
        """
        path = self.normalize_path(path)
        all_names = list(getattr(self._get_object(path), "col_names", []) or [])
        chosen = set(names or [])
        selection = [name for name in all_names if name in chosen]
        if not selection or len(selection) == len(all_names):
            # Nothing chosen, or everything chosen -> no narrowing; clear.
            self._column_filters.pop(path, None)
            self._column_selections.pop(path, None)
            return len(all_names)
        self._column_filters[path] = f"{len(selection)} of {len(all_names)}"
        self._column_selections[path] = selection
        return len(selection)

    def get_column_filter(self, path: str) -> str | None:
        """Return the active column filter descriptor for *path*, if any."""
        return self._column_filters.get(self.normalize_path(path))

    def base_ncols(self, path: str) -> int:
        """Return the unfiltered column count of the CTable at *path*."""
        return len(list(getattr(self._get_object(path), "col_names", []) or []))

    def base_column_names(self, path: str) -> list[str]:
        """Return all column names of the CTable at *path*, ignoring any filter."""
        return list(getattr(self._get_object(path), "col_names", []) or [])

    def _get_object(self, path: str) -> Any:
        """Return the object represented by *path*."""
        path = self.normalize_path(path)
        if self.is_tree:
            return self.store[path]
        self._check_root_path(path)
        return self.store

    @staticmethod
    def _check_root_path(path: str) -> None:
        if path != "/":
            raise KeyError(f"Standalone objects only expose the root path '/', got {path!r}")

    _INTERNAL_VLMETA_KEYS = frozenset(
        {
            "kind",
            "version",
            "schema",
            "n_rows",
            "value_epoch",
            "computed_columns",
            "materialized_columns",
        }
    )

    @staticmethod
    def _vlmeta_dict(vlmeta) -> dict[str, Any] | None:
        if vlmeta is None:
            return None
        try:
            data = vlmeta[:]
        except Exception:
            try:
                data = {name: vlmeta[name] for name in vlmeta}
            except Exception:
                return None
        if data is None:
            return None
        # Filter out internal blosc2 metadata keys (schema, version, etc.)
        return {k: v for k, v in data.items() if k not in StoreBrowser._INTERNAL_VLMETA_KEYS}


def object_kind(obj: Any) -> str:
    """Return a stable b2view kind string for *obj*."""
    if isinstance(obj, blosc2.TreeStore):
        return "group"
    if isinstance(obj, blosc2.NDArray):
        return "ndarray"
    if isinstance(obj, blosc2.CTable):
        return "ctable"
    if hasattr(blosc2, "C2Array") and isinstance(obj, blosc2.C2Array):
        return "c2array"
    if isinstance(obj, blosc2.SChunk):
        return "schunk"
    return "unknown"


def object_metadata(obj: Any) -> dict[str, Any]:
    """Extract lightweight metadata from a supported object."""
    kind = object_kind(obj)
    if kind in {"ndarray", "c2array"}:
        return {
            "shape": getattr(obj, "shape", None),
            "ndim": len(getattr(obj, "shape", ()) or ()),
            "dtype": str(getattr(obj, "dtype", None)),
            "chunks": getattr(obj, "chunks", None),
            "blocks": getattr(obj, "blocks", None),
            "nbytes": getattr(obj, "nbytes", None),
            "cbytes": getattr(obj, "cbytes", None),
        }
    if kind == "ctable":
        try:
            return dict(obj.info_items)
        except Exception:
            return {
                "nrows": getattr(obj, "nrows", len(obj)),
                "ncols": getattr(obj, "ncols", len(getattr(obj, "col_names", []))),
                "columns": {
                    name: str(getattr(obj[name], "dtype", None)) for name in getattr(obj, "col_names", [])
                },
            }
    if kind == "schunk":
        return {
            "chunks": getattr(obj, "nchunks", None),
            "nbytes": getattr(obj, "nbytes", None),
            "cbytes": getattr(obj, "cbytes", None),
        }
    return {"repr": repr(obj)}


def preview_array_from_layout(
    obj: Any,
    *,
    layout: DataSliceLayout,
    max_rows: int = 20,
    max_cols: int = 10,
) -> dict[str, Any]:
    """Return a bounded preview for an N-D array using a *layout*.

    The layout describes which dimensions are fixed (slider) vs navigable
    (table rows/columns).  At most 2 navigable dimensions are allowed.
    """
    shape = tuple(getattr(obj, "shape", ()) or ())
    if len(shape) != len(layout.shape):
        raise ValueError(f"Layout shape {layout.shape} does not match object shape {shape}")
    ndim = len(shape)
    navigable = layout.navigable_dims

    # Determine row and col navigable dims
    row_dim = navigable[0] if len(navigable) >= 1 else None
    col_dim = navigable[1] if len(navigable) >= 2 else None

    # Page sizes.  A locked row window narrows the navigable row dim to
    # [win_lo, win_hi): the grid sees only ``nrows`` rows (so paging cannot
    # leave it) and every read is offset by ``win_lo``.
    win_lo, win_hi = layout.row_window_bounds(row_dim)
    nrows = (win_hi - win_lo) if row_dim is not None else 1
    ncols = shape[col_dim] if col_dim is not None else 1

    # Clamp fixed values
    fixed_values = {}
    for d, val in layout.fixed_values.items():
        total = shape[d]
        fixed_values[d] = max(0, min(val, total - 1)) if total > 0 else 0

    # Ensure every non-navigable dim is fixed at 0 (safety catch)
    for i in range(ndim):
        if i not in fixed_values and (row_dim is None or i != row_dim) and (col_dim is None or i != col_dim):
            fixed_values[i] = 0

    # Build slicing tuple
    idx: list[int | slice] = []
    for i in range(ndim):
        if i in fixed_values:
            idx.append(fixed_values[i])
        elif row_dim is not None and i == row_dim:
            # ``layout.row_start`` is window-relative; offset into the array.
            start = max(0, min(layout.row_start, nrows))
            stop = min(start + max_rows, nrows)
            idx.append(slice(win_lo + start, win_lo + stop))
        elif col_dim is not None and i == col_dim:
            col_start = max(0, min(layout.col_start, ncols))
            col_stop = min(col_start + max_cols, ncols)
            idx.append(slice(col_start, col_stop))
        else:
            # Shouldn't happen: non-navigable dims are caught above
            idx.append(slice(0, shape[i]))

    values = np.asarray(obj[tuple(idx)])

    # Build column labels — match data keys below
    if col_dim is not None:
        col_start = max(0, min(layout.col_start, ncols))
        col_stop = min(col_start + max_cols, ncols)
        columns = [str(i) for i in range(col_start, col_stop)]
    elif row_dim is not None:
        columns = ["value"]
    else:
        columns = ["value"]

    # Extract 2-D data from result
    data: dict[str, Any] = {}
    if row_dim is not None and col_dim is not None:
        # 2-D navigable → 2-D table
        col_start = max(0, min(layout.col_start, ncols))
        col_stop = min(col_start + max_cols, ncols)
        for i, c in enumerate(range(col_start, col_stop)):
            data[str(c)] = values[:, i]
    elif row_dim is not None:
        # Only rows navigable → 1-D view
        data["value"] = values
    else:
        # 0 navigable → scalar
        data["value"] = np.asarray([values.item()]) if np.ndim(values) == 0 else np.asarray([values])

    row_start_val = max(0, min(layout.row_start, nrows)) if row_dim is not None else 0
    row_stop_val = min(row_start_val + max_rows, nrows) if row_dim is not None else 1
    col_start_val = max(0, min(layout.col_start, ncols)) if col_dim is not None else 0
    col_stop_val = min(col_start_val + max_cols, ncols) if col_dim is not None else 1

    result: dict[str, Any] = {
        "start": row_start_val,
        "stop": row_stop_val,
        "nrows": nrows,
        "columns": columns,
        "hidden_columns": max(0, ncols - (col_stop_val - col_start_val)),
        "data": data,
        "source_kind": "ndarray_slice",
        "shape": shape,
        "col_start": col_start_val,
        "col_stop": col_stop_val,
        "ncols": ncols,
        "layout": layout,
        "slice_indices": [fixed_values.get(i, 0) for i in range(min(ndim - 2, ndim))],
        "n_slices_per_dim": [shape[i] for i in range(ndim) if i in fixed_values],
    }
    # Keep legacy fields for backward compat
    result["slice_indices"] = [fixed_values.get(i, 0) for i in range(ndim) if i in fixed_values]
    result["n_slices_per_dim"] = [shape[i] for i in range(ndim) if i in fixed_values]
    return result


def preview_array_nd_slice(
    obj: Any,
    *,
    slice_indices: list[int] | None = None,
    start: int = 0,
    stop: int = 20,
    col_start: int = 0,
    max_cols: int = 10,
) -> dict[str, Any]:
    """Return a bounded 2-D slice preview for N-D arrays (N >= 3)."""
    shape = tuple(getattr(obj, "shape", ()) or ())
    ndim = len(shape)
    if ndim < 3:
        raise ValueError(f"Expected an N-D array with N >= 3, got shape {shape!r}")
    n_leading = ndim - 2
    n_slices_per_dim = list(shape[:n_leading])
    if slice_indices is None or len(slice_indices) != n_leading:
        slice_indices = [0] * n_leading
    # Clamp
    slice_indices = [
        min(max(0, idx), n_slices_per_dim[i] - 1) if n_slices_per_dim[i] > 0 else 0
        for i, idx in enumerate(slice_indices)
    ]
    nrows, ncols = shape[-2], shape[-1]
    if stop is None:
        stop = min(start + 20, nrows)
    start = max(0, min(start, nrows))
    stop = min(max(start, stop), nrows)
    col_start = max(0, min(col_start, ncols))
    col_stop = min(col_start + max_cols, ncols)
    columns = [str(i) for i in range(col_start, col_stop)]
    idx = tuple(slice_indices) + (slice(start, stop), slice(col_start, col_stop))
    values = np.asarray(obj[idx])
    data = {str(col): values[:, i] for i, col in enumerate(range(col_start, col_stop))}
    return {
        "start": start,
        "stop": stop,
        "nrows": nrows,
        "columns": columns,
        "hidden_columns": max(0, ncols - (col_stop - col_start)),
        "data": data,
        "source_kind": "ndarray_slice",
        "shape": shape,
        "col_start": col_start,
        "col_stop": col_stop,
        "ncols": ncols,
        "slice_indices": slice_indices,
        "n_slices_per_dim": n_slices_per_dim,
    }


def preview_array_2d(
    obj: Any, *, start: int = 0, stop: int = 20, col_start: int = 0, max_cols: int = 10
) -> dict[str, Any]:
    """Return a bounded row/column preview for a 2-D array."""
    shape = tuple(getattr(obj, "shape", ()) or ())
    if len(shape) != 2:
        raise ValueError(f"Expected a 2-D array, got shape {shape!r}")
    nrows, ncols = shape
    start = max(0, min(start, nrows))
    stop = min(max(start, stop), nrows)
    col_start = max(0, min(col_start, ncols))
    col_stop = min(col_start + max_cols, ncols)
    columns = [str(i) for i in range(col_start, col_stop)]
    values = np.asarray(obj[(slice(start, stop), slice(col_start, col_stop))])
    data = {str(col): values[:, i] for i, col in enumerate(range(col_start, col_stop))}
    return {
        "start": start,
        "stop": stop,
        "nrows": nrows,
        "columns": columns,
        "hidden_columns": max(0, ncols - (col_stop - col_start)),
        "data": data,
        "source_kind": "ndarray2d",
        "shape": shape,
        "col_start": col_start,
        "col_stop": col_stop,
        "ncols": ncols,
    }


def preview_array_1d(obj: Any, *, start: int = 0, stop: int = 20, **kwargs) -> dict[str, Any]:
    """Return a bounded row preview for a 1-D array."""
    shape = tuple(getattr(obj, "shape", ()) or ())
    if len(shape) != 1:
        raise ValueError(f"Expected a 1-D array, got shape {shape!r}")
    nrows = shape[0]
    start = max(0, min(start, nrows))
    stop = min(max(start, stop), nrows)
    data = {
        "value": np.asarray(obj[start:stop]),
    }
    return {
        "start": start,
        "stop": stop,
        "nrows": nrows,
        "columns": ["value"],
        "hidden_columns": 0,
        "data": data,
        "source_kind": "ndarray1d",
        "shape": shape,
    }


def preview_array(
    obj: Any, *, slices: tuple[Any, ...] | None = None, max_rows: int = 20, max_cols: int = 10
):
    """Return a small NumPy preview from an NDArray/C2Array-like object."""
    shape = tuple(getattr(obj, "shape", ()) or ())
    if slices is None:
        if len(shape) == 0:
            slices = ()
        elif len(shape) == 1:
            slices = (slice(0, min(shape[0], max_rows)),)
        elif len(shape) == 2:
            slices = (slice(0, min(shape[0], max_rows)), slice(0, min(shape[1], max_cols)))
        else:
            leading = tuple(0 for _ in shape[:-2])
            slices = leading + (
                slice(0, min(shape[-2], max_rows)),
                slice(0, min(shape[-1], max_cols)),
            )
    return np.asarray(obj[slices])


def preview_ctable(
    obj: Any,
    *,
    start: int = 0,
    stop: int = 20,
    columns: list[str] | None = None,
    max_cols: int = 10,
    col_start: int = 0,
    include_expensive: bool = False,
) -> dict[str, Any]:
    """Return a bounded column-oriented preview from a CTable.

    *col_start* selects the first visible column, so wide tables can be
    paged horizontally just like 2-D arrays.

    Complex nested/list/object columns may require one variable-length block
    read per row.  By default, keep table navigation responsive by showing a
    placeholder for those columns instead of decoding them eagerly.
    """
    all_columns = list(getattr(obj, "col_names", []))
    selectable = all_columns if columns is None else [name for name in columns if name in all_columns]
    ncols = len(selectable)
    col_start = max(0, min(col_start, max(0, ncols - 1)))
    col_stop = min(col_start + max_cols, ncols)
    visible_columns = selectable[col_start:col_stop]
    hidden_columns = max(0, ncols - len(visible_columns))
    start = max(0, start)
    stop = min(max(start, stop), len(obj))
    data = {}
    skipped_columns = {}
    nrows = stop - start
    for name in visible_columns:
        if not include_expensive and is_expensive_ctable_column(obj, name):
            label = ctable_column_label(obj, name)
            placeholder = f"<{label}; skipped>"
            data[name] = np.full(nrows, placeholder, dtype=object)
            skipped_columns[name] = label
        else:
            data[name] = safe_asarray(obj[name][start:stop])
    return {
        "start": start,
        "stop": stop,
        "nrows": len(obj),
        "columns": visible_columns,
        "hidden_columns": hidden_columns,
        "skipped_columns": skipped_columns,
        "data": data,
        "source_kind": "ctable",
        "col_start": col_start,
        "col_stop": col_stop,
        "ncols": ncols,
    }


def schunk_row_geometry(typesize: int) -> tuple[int, int]:
    """Return ``(items_per_row, bytes_per_row)`` for the hex dump.

    Bytes are grouped into ``typesize``-wide items (so a 4-byte typesize shows
    32-bit words); ``items_per_row`` is chosen so a row is ~16 bytes wide, and
    never below one whole item.
    """
    typesize = max(1, int(typesize or 1))
    items_per_row = max(1, 16 // typesize)
    return items_per_row, items_per_row * typesize


def preview_schunk(obj: Any, *, start: int = 0, stop: int = 20) -> dict[str, Any]:
    """Return a bounded ``xxd``-style hex dump of an SChunk's raw bytes.

    Each grid row is one ``bytes_per_row`` span (a multiple of ``typesize``);
    *start*/*stop* are in those row units, so the existing row-paging machinery
    applies unchanged.  Only the visible byte span is read (``obj[a:b]``), so a
    multi-GB SChunk previews instantly.  The byte offset is the row label.
    """
    nbytes = int(getattr(obj, "nbytes", 0) or 0)
    typesize = max(1, int(getattr(obj, "typesize", 1) or 1))
    items_per_row, bytes_per_row = schunk_row_geometry(typesize)
    total_rows = (nbytes + bytes_per_row - 1) // bytes_per_row
    start = max(0, start)
    stop = min(max(start, stop), total_rows)
    byte_start = start * bytes_per_row
    byte_stop = min(stop * bytes_per_row, nbytes)
    raw = bytes(obj[byte_start:byte_stop]) if byte_stop > byte_start else b""
    hex_width = items_per_row * typesize * 2 + (items_per_row - 1)
    hex_col: list[str] = []
    ascii_col: list[str] = []
    labels: list[str] = []
    for r in range(stop - start):
        chunk = raw[r * bytes_per_row : (r + 1) * bytes_per_row]
        items = [chunk[k : k + typesize].hex() for k in range(0, len(chunk), typesize)]
        hex_col.append(" ".join(items).ljust(hex_width))
        ascii_col.append("".join(chr(b) if 0x20 <= b <= 0x7E else "." for b in chunk))
        labels.append(format(byte_start + r * bytes_per_row, "08x"))
    return {
        "start": start,
        "stop": stop,
        "nrows": total_rows,
        "columns": ["hex", "ascii"],
        "hidden_columns": 0,
        "row_labels": labels,
        "data": {
            "hex": np.array(hex_col, dtype=object),
            "ascii": np.array(ascii_col, dtype=object),
        },
        "source_kind": "schunk",
        "typesize": typesize,
        "nbytes": nbytes,
    }


def is_expensive_ctable_column(obj: Any, name: str) -> bool:
    """Return whether previewing a CTable column is likely row-by-row expensive."""
    try:
        schema = obj.schema_dict()
    except Exception:
        return False
    for column in schema.get("columns", []):
        if column.get("name") != name:
            continue
        return column.get("kind") in {"list", "struct", "object", "ndarray"}
    return False


def ctable_column_label(obj: Any, name: str) -> str:
    """Return a compact schema label for *name*."""
    try:
        columns = dict(obj.info_items).get("columns", {})
        label = columns.get(name)
        if label is not None:
            # Strip the trailing size annotation, e.g. "list[struct] (cbytes: ...)".
            return str(label).split(" (", 1)[0]
    except Exception:
        pass
    try:
        for column in obj.schema_dict().get("columns", []):
            if column.get("name") == name:
                return str(column.get("kind", "complex"))
    except Exception:
        pass
    return "complex"


def safe_asarray(values: Any) -> np.ndarray:
    """Convert preview values to an array, preserving ragged/nested values.

    NumPy 2 raises for ragged nested sequences unless ``dtype=object`` is
    requested explicitly.  CTable columns can legitimately contain list/struct
    values, so previews must keep those as object cells instead of failing.
    """
    try:
        return np.asarray(values)
    except ValueError:
        return np.asarray(values, dtype=object)
