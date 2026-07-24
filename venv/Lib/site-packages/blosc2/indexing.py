#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import ast
import contextlib
import enum
import hashlib
import math
import os
import re
import sys
import tempfile
import warnings
import weakref
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import blosc2

from . import indexing_ext

INDEXES_VLMETA_KEY = "blosc2_indexes"
INDEX_FORMAT_VERSION = 1
SELF_TARGET_NAME = "__self__"

# On Windows, mmap holds file locks that prevent later writes (vlmeta updates,
# sidecar recreation during rebuild_index, etc.).  Disable mmap for all index
# I/O on that platform.
_INDEX_MMAP_MODE = None if sys.platform == "win32" else "r"

FLAG_ALL_NAN = np.uint8(1 << 0)
FLAG_HAS_NAN = np.uint8(1 << 1)

SEGMENT_LEVELS_BY_KIND = {
    # SUMMARY stores per-segment min/max.  Block granularity prunes far more
    # finely than whole-chunk min/max (chunks can be ~10-100x larger than
    # blocks), so it is the default; callers may override per index via the
    # ``granularity`` option (see ``SUMMARY_GRANULARITIES``).
    "summary": ("block",),
    "bucket": ("chunk", "block"),
    "partial": ("chunk", "block", "subblock"),
    "full": ("chunk", "block", "subblock"),
    "opsi": ("chunk", "block", "subblock"),
}

# Valid ``granularity`` values for SUMMARY indexes, coarsest to finest.
SUMMARY_GRANULARITIES = ("chunk", "block")


def _resolve_summary_levels(granularity: str | None) -> tuple[str, ...] | None:
    """Map a SUMMARY ``granularity`` to the level tuple, or ``None`` for default."""
    if granularity is None:
        return None
    if granularity not in SUMMARY_GRANULARITIES:
        raise ValueError(f"granularity must be one of {SUMMARY_GRANULARITIES}, got {granularity!r}")
    return (granularity,)


_IN_MEMORY_INDEXES: dict[int, dict] = {}
_IN_MEMORY_INDEX_FINALIZERS: dict[int, weakref.finalize] = {}
_PERSISTENT_INDEXES: dict[tuple[str, str | int], dict] = {}
# Records which array object owns the descriptor cached at a given
# ``(array_key, token)``.  In compact stores (.b2z/.b2d) every CTable column
# shares one urlpath, so ``_array_key`` collides across columns and a column's
# index store can be returned for a *different* column's predicate.  This map
# lets ``_descriptor_for_target`` reject a descriptor when the querying array is
# not the one the descriptor was registered for.  It is in-memory only and is
# never persisted to disk.
_DESCRIPTOR_OWNERS: dict[tuple, int] = {}
_DATA_CACHE: dict[tuple[int, str | None, str, str], np.ndarray] = {}
_SIDECAR_HANDLE_CACHE: dict[tuple[int, str | None, str, str], object] = {}

# ---------------------------------------------------------------------------
# Query-result cache constants and global state
# ---------------------------------------------------------------------------
QUERY_CACHE_VLMETA_KEY = "_blosc2_query_cache"
QUERY_CACHE_FORMAT_VERSION = 1
QUERY_CACHE_MAX_ENTRY_NBYTES = 65_536  # 64 KB of logical int64 positions per persistent entry
QUERY_CACHE_MAX_MEM_NBYTES = 131_072  # 128 KB for the in-process hot cache
QUERY_CACHE_MAX_PERSISTENT_NBYTES = 4 * 1024 * 1024  # 4 MB of logical int64 positions in the payload store


@dataclass(frozen=True, slots=True)
class _CompressedHotCoords:
    dtype: str
    nrows: int
    compressed: bool
    data: bytes

    @property
    def nbytes(self) -> int:
        return len(self.data)


# In-process hot cache: (array-scope, digest) -> compressed coordinate payload.
_HOT_CACHE: dict[tuple[tuple[str, str | int], str], _CompressedHotCoords] = {}
# Insertion-order list for LRU eviction.
_HOT_CACHE_ORDER: list[tuple[tuple[str, str | int], str]] = []
# Total bytes of arrays currently in the hot cache.
_HOT_CACHE_BYTES: int = 0
# Legacy query-cache sidecar handles: resolved urlpath -> open ObjectArray object.
# Query caches are hot-cache-only now, but we keep this state so invalidation can
# still drop stale artifacts produced by older versions.
_QUERY_CACHE_STORE_HANDLES: dict[str, object] = {}
# Cached mmap handles for data arrays used in full-query gather: urlpath -> NDArray.
_GATHER_MMAP_HANDLES: dict[str, object] = {}
# Last-seen on-disk fingerprint (mtime_ns, size) for persistent paths whose query
# handles/coordinates are cached above. Lets an in-place overwrite (same path, new
# contents) be detected and invalidated, not just an outright deletion.
_PERSISTENT_FINGERPRINTS: dict[str, tuple[int, int]] = {}
# Registry for index sidecar files stored inside a .b2z bundle.
# Maps absolute sidecar path -> (b2z_path, data_offset_in_zip).
# Populated by the storage layer so indexing code can open sidecars without
# extracting them to a temporary directory first.
_SIDECAR_ZIP_REGISTRY: dict[str, tuple[str, int]] = {}
_HOT_CACHE_GLOBAL_SCOPE = ("global", 0)

FULL_OOC_RUN_ITEMS = 10_000_000
FULL_OOC_MERGE_BUFFER_ITEMS = 500_000
FULL_SELECTIVE_OOC_MAX_SPANS = 128
FULL_RUN_BOUNDED_FALLBACK_RUNS = 8
FULL_RUN_BOUNDED_FALLBACK_ITEMS = 1_000_000
INDEX_QUERY_MIN_CHUNKS_PER_THREAD = 8


def _python_executor_threads(requested_threads: int) -> int:
    # wasm32 builds do not support spawning Python worker threads reliably.
    if blosc2.IS_WASM:
        return 1
    return max(1, int(requested_threads))


def _sanitize_token(token: str) -> str:
    return re.sub(r"[^0-9A-Za-z_.-]+", "_", token)


def _cleanup_in_memory_store(key: int) -> None:
    _IN_MEMORY_INDEXES.pop(key, None)
    _IN_MEMORY_INDEX_FINALIZERS.pop(key, None)
    scope = ("memory", key)
    stale_data = [cache_key for cache_key in _DATA_CACHE if cache_key[0] == scope]
    for cache_key in stale_data:
        _DATA_CACHE.pop(cache_key, None)
    stale_handles = [cache_key for cache_key in _SIDECAR_HANDLE_CACHE if cache_key[0] == scope]
    for cache_key in stale_handles:
        _SIDECAR_HANDLE_CACHE.pop(cache_key, None)
    _hot_cache_clear(scope=("memory", key))


def _persistent_cache_path_exists(path: str | int) -> bool:
    if not isinstance(path, str):
        return False
    path_obj = Path(path)
    return path_obj.exists() or path_obj.parent.exists()


def _purge_stale_persistent_caches() -> None:
    stale_scopes = {
        key
        for key in tuple(_PERSISTENT_INDEXES)
        if key[0] == "persistent" and not _persistent_cache_path_exists(key[1])
    }
    for key in stale_scopes:
        _PERSISTENT_INDEXES.pop(key, None)

    stale_data_keys = {
        key
        for key in tuple(_DATA_CACHE)
        if key[0][0] == "persistent" and not _persistent_cache_path_exists(key[0][1])
    }
    stale_scopes.update(key[0] for key in stale_data_keys)
    for key in stale_data_keys:
        _DATA_CACHE.pop(key, None)

    stale_handle_keys = {
        key
        for key in tuple(_SIDECAR_HANDLE_CACHE)
        if key[0][0] == "persistent" and not _persistent_cache_path_exists(key[0][1])
    }
    stale_scopes.update(key[0] for key in stale_handle_keys)
    for key in stale_handle_keys:
        _SIDECAR_HANDLE_CACHE.pop(key, None)

    stale_query_paths = [path for path in tuple(_QUERY_CACHE_STORE_HANDLES) if not Path(path).exists()]
    for path in stale_query_paths:
        _QUERY_CACHE_STORE_HANDLES.pop(path, None)

    stale_gather_paths = [path for path in tuple(_GATHER_MMAP_HANDLES) if not Path(path).exists()]
    for path in stale_gather_paths:
        _GATHER_MMAP_HANDLES.pop(path, None)

    stale_fingerprints = [path for path in tuple(_PERSISTENT_FINGERPRINTS) if not Path(path).exists()]
    for path in stale_fingerprints:
        _PERSISTENT_FINGERPRINTS.pop(path, None)

    for scope in stale_scopes:
        _hot_cache_clear(scope=scope)


def _file_fingerprint(path: str) -> tuple[int, int] | None:
    """Return a cheap change-detection fingerprint ``(mtime_ns, size)`` for *path*."""
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_mtime_ns, st.st_size)


def _drop_persistent_path_caches(raw_path: str, resolved: str) -> None:
    """Drop every process-global cache entry keyed by the persistent file *path*."""
    scope = ("persistent", resolved)
    _PERSISTENT_INDEXES.pop(scope, None)
    for cache in (_DATA_CACHE, _SIDECAR_HANDLE_CACHE):
        for key in [k for k in tuple(cache) if k[0] == scope]:
            cache.pop(key, None)
    _hot_cache_clear(scope=scope)
    for handles in (_QUERY_CACHE_STORE_HANDLES, _GATHER_MMAP_HANDLES):
        handles.pop(raw_path, None)
        handles.pop(resolved, None)


def _refresh_persistent_caches(path: str | None) -> None:
    """Invalidate cached handles/coordinates for *path* if its file changed on disk.

    Reuse of the process-global query caches is keyed by urlpath, and the entries are
    otherwise only dropped once their file is *deleted* (see
    :func:`_purge_stale_persistent_caches`).  A file that is overwritten in place (for
    instance re-uploaded with new contents at the same path) would therefore still be
    served from the stale mmap handle and coordinate cache of the previous file.

    Comparing the current ``(mtime_ns, size)`` against the fingerprint recorded when the
    caches were populated lets us detect that and drop the affected entries; they simply
    repopulate on the next query.
    """
    if not path:
        return
    raw_path = str(path)
    try:
        resolved = str(Path(raw_path).resolve())
    except OSError:
        return
    fingerprint = _file_fingerprint(resolved)
    if fingerprint is None:
        return
    previous = _PERSISTENT_FINGERPRINTS.get(resolved)
    if previous is not None and previous != fingerprint:
        _drop_persistent_path_caches(raw_path, resolved)
    _PERSISTENT_FINGERPRINTS[resolved] = fingerprint


def evict_cached_index_handles(root: str | None) -> None:
    """Drop cached sidecar handles/data for the persistent store at *root*.

    Index reads cache file-backed handles in process-global dicts for query
    reuse; they are normally only purged once their files are deleted, so a
    table closed but left on disk keeps its descriptors open — one per table,
    which exhausts the file-descriptor limit over a large session.  Closing a
    table calls this to pop (and thereby release) the handles it owns; the
    caches simply repopulate on the next query.
    """
    if not root:
        return
    try:
        resolved = str(Path(root).resolve())
    except Exception:
        return
    prefix = resolved + os.sep

    def _owned_scope(scope) -> bool:
        # scope is an _array_key: ("persistent", path) or ("memory", id).
        return (
            isinstance(scope, tuple)
            and len(scope) == 2
            and scope[0] == "persistent"
            and isinstance(scope[1], str)
            and (scope[1] == resolved or scope[1].startswith(prefix))
        )

    def _owned_path(path) -> bool:
        return isinstance(path, str) and (path == resolved or path.startswith(prefix))

    for cache in (_SIDECAR_HANDLE_CACHE, _DATA_CACHE, _HOT_CACHE):
        for key in [k for k in tuple(cache) if _owned_scope(k[0])]:
            cache.pop(key, None)
    for handles in (_QUERY_CACHE_STORE_HANDLES, _GATHER_MMAP_HANDLES):
        for path in [p for p in tuple(handles) if _owned_path(p)]:
            handles.pop(path, None)


def _open_sidecar_file(path: str, mmap_mode=None) -> blosc2.NDArray:
    """Open an index sidecar file, using zip-offset access when registered."""
    reg = _SIDECAR_ZIP_REGISTRY.get(path)
    if reg is not None:
        b2z_path, offset = reg
        from blosc2.schunk import process_opened_object

        return process_opened_object(
            blosc2.blosc2_ext.open(b2z_path, mode="r", offset=offset, mmap_mode=mmap_mode)
        )
    return blosc2.open(path, mode="r", mmap_mode=mmap_mode)


@dataclass(slots=True)
class IndexPlan:
    usable: bool
    reason: str
    descriptor: dict | None = None
    base: blosc2.NDArray | None = None
    target: dict | None = None
    field: str | None = None
    level: str | None = None
    segment_len: int | None = None
    candidate_units: np.ndarray | None = None
    total_units: int = 0
    selected_units: int = 0
    exact_positions: np.ndarray | None = None
    bucket_masks: np.ndarray | None = None
    bucket_len: int | None = None
    chunk_len: int | None = None
    block_len: int | None = None
    lower: object | None = None
    lower_inclusive: bool = True
    upper: object | None = None
    upper_inclusive: bool = True
    candidate_chunks: int = 0
    candidate_nav_segments: int = 0
    candidate_base_spans: int = 0
    lookup_path: str | None = None
    # Cross-column refinement: exact positions from one indexed column that
    # still need refinement against other predicates (different columns).
    partial_exact_positions: np.ndarray | None = None


@dataclass(slots=True)
class SegmentPredicatePlan:
    base: blosc2.NDArray
    candidate_units: np.ndarray
    descriptor: dict
    target: dict
    field: str | None
    level: str
    segment_len: int


@dataclass(slots=True)
class ExactPredicatePlan:
    base: blosc2.NDArray
    descriptor: dict
    target: dict
    field: str | None
    lower: object | None = None
    lower_inclusive: bool = True
    upper: object | None = None
    upper_inclusive: bool = True


@dataclass(slots=True)
class SortedRun:
    values_path: Path
    positions_path: Path
    length: int


@dataclass(slots=True)
class TempRunTracker:
    current_disk_bytes: int = 0
    peak_disk_bytes: int = 0
    total_written_bytes: int = 0


@dataclass(slots=True)
class OrderedIndexPlan:
    usable: bool
    reason: str
    descriptor: dict | None = None
    base: blosc2.NDArray | None = None
    field: str | None = None
    order_fields: list[str | None] | None = None
    total_rows: int = 0
    selected_rows: int = 0
    secondary_refinement: bool = False


@dataclass(frozen=True, slots=True)
class IndexComponent:
    label: str
    category: str
    name: str
    path: str | None


def _default_index_store() -> dict:
    return {"version": INDEX_FORMAT_VERSION, "indexes": {}}


def _array_key(array: blosc2.NDArray) -> tuple[str, str | int]:
    if _is_persistent_array(array):
        return ("persistent", str(Path(array.urlpath).resolve()))
    return ("memory", id(array))


def _register_descriptor_owner(array: blosc2.NDArray, token: str) -> None:
    """Record that *array* owns the descriptor cached at ``(_array_key, token)``.

    Called when a CTable injects a column's index descriptor into the shared
    (urlpath-keyed) persistent store, so that :func:`_descriptor_for_target` can
    later tell that descriptor apart from a sibling column whose ``_array_key``
    collides.
    """
    _DESCRIPTOR_OWNERS[(_array_key(array), token)] = id(array)


def _field_token(field: str | None) -> str:
    return "__self__" if field is None else field


def _target_token(target: dict) -> str:
    source = target.get("source")
    if source == "field":
        return _field_token(target.get("field"))
    if source == "expression":
        digest = hashlib.sha1(target["expression_key"].encode("utf-8")).hexdigest()[:12]
        return f"__expr__{digest}"
    raise ValueError(f"unsupported index target source {source!r}")


def _copy_nested_dict(value: dict | None) -> dict | None:
    if value is None:
        return None
    copied = value.copy()
    for key, item in list(copied.items()):
        if isinstance(item, dict):
            copied[key] = item.copy()
    return copied


def _copy_descriptor(descriptor: dict) -> dict:
    copied = descriptor.copy()
    if descriptor.get("cparams") is not None:
        copied["cparams"] = descriptor["cparams"].copy()
    copied["levels"] = _copy_nested_dict(descriptor.get("levels"))
    if descriptor.get("target") is not None:
        copied["target"] = descriptor["target"].copy()
    if descriptor.get("bucket") is not None:
        copied["bucket"] = descriptor["bucket"].copy()
    if descriptor.get("partial") is not None:
        copied["partial"] = descriptor["partial"].copy()
    if descriptor.get("full") is not None:
        copied["full"] = descriptor["full"].copy()
        if "runs" in copied["full"]:
            copied["full"]["runs"] = [run.copy() for run in copied["full"]["runs"]]
    if descriptor.get("opsi") is not None:
        copied["opsi"] = descriptor["opsi"].copy()
    return copied


def _descriptor_for_token(array: blosc2.NDArray, token: str) -> dict:
    descriptor = _load_store(array)["indexes"].get(token)
    if descriptor is None:
        raise KeyError("index not found")
    return descriptor


def _copy_descriptor_for_token(array: blosc2.NDArray, token: str) -> dict:
    return _copy_descriptor(_descriptor_for_token(array, token))


def _is_persistent_array(array: blosc2.NDArray) -> bool:
    return getattr(array, "urlpath", None) is not None


def _tmpdir_for_array(array: blosc2.NDArray) -> str | None:
    """Return a directory on the same filesystem as *array* for temp files.

    When the array is backed by a file on disk we place temporaries next
    to it so that they share the same filesystem.  This avoids hitting
    size limits on ``/tmp`` (commonly a tmpfs with only a few GB).
    For in-memory arrays we fall back to the system default (``None``).
    """
    urlpath = getattr(array, "urlpath", None)
    if urlpath is not None:
        return str(Path(urlpath).resolve().parent)
    return None


def _resolve_full_index_tmpdir(array: blosc2.NDArray, tmpdir: str | None) -> str | None:
    """Resolve the workspace for temporary files used by OOC sorted index builds.

    An explicit ``tmpdir`` wins. Otherwise, persistent arrays default to the
    directory that stores the array so temporary sidecars stay on the same
    filesystem. In-memory arrays fall back to the system temporary directory.
    """
    if tmpdir is not None:
        return tmpdir
    return _tmpdir_for_array(array)


def _load_store(array: blosc2.NDArray) -> dict:
    _purge_stale_persistent_caches()
    if _is_persistent_array(array):
        key = _array_key(array)
        cached = _PERSISTENT_INDEXES.get(key)
        if cached is not None:
            return cached
        try:
            store = array.schunk.vlmeta[INDEXES_VLMETA_KEY]
        except KeyError:
            store = _default_index_store()
        if not isinstance(store, dict):
            store = _default_index_store()
        store.setdefault("version", INDEX_FORMAT_VERSION)
        store.setdefault("indexes", {})
        _PERSISTENT_INDEXES[key] = store
        return store

    key = id(array)
    cached = _IN_MEMORY_INDEXES.get(key)
    if cached is not None:
        return cached
    store = _default_index_store()
    _IN_MEMORY_INDEXES[key] = store
    _IN_MEMORY_INDEX_FINALIZERS[key] = weakref.finalize(array, _cleanup_in_memory_store, key)
    return store


def _save_store(array: blosc2.NDArray, store: dict) -> None:
    store.setdefault("version", INDEX_FORMAT_VERSION)
    store.setdefault("indexes", {})
    if _is_persistent_array(array):
        _PERSISTENT_INDEXES[_array_key(array)] = store
        array.schunk.vlmeta[INDEXES_VLMETA_KEY] = store
    else:
        key = id(array)
        _IN_MEMORY_INDEXES[key] = store
        _IN_MEMORY_INDEX_FINALIZERS.setdefault(key, weakref.finalize(array, _cleanup_in_memory_store, key))


# ---------------------------------------------------------------------------
# Stage 1 – Query cache: metadata helpers and container plumbing
# ---------------------------------------------------------------------------


def _query_cache_payload_path(array: blosc2.NDArray) -> str:
    """Return the path for the persistent query-cache ObjectArray payload store."""
    path, root = _sanitize_sidecar_root(array.urlpath)
    return str(path.with_name(f"{root}.__query_cache__.b2frame"))


def _query_cache_owner(array: blosc2.NDArray) -> blosc2.NDArray:
    owner = getattr(array, "ndarr", None)
    return owner if owner is not None else array


def _ensure_in_memory_array_finalizer(array: blosc2.NDArray) -> None:
    if _is_persistent_array(array):
        return
    key = id(array)
    _IN_MEMORY_INDEX_FINALIZERS.setdefault(key, weakref.finalize(array, _cleanup_in_memory_store, key))


def _query_cache_scope(array: blosc2.NDArray) -> tuple[str, str | int]:
    owner = _query_cache_owner(array)
    _ensure_in_memory_array_finalizer(owner)
    return _array_key(owner)


def _default_query_cache_catalog(payload_path: str) -> dict:
    return {
        "version": QUERY_CACHE_FORMAT_VERSION,
        "payload_ref": {"kind": "urlpath", "version": 1, "urlpath": payload_path},
        "max_entry_nbytes": QUERY_CACHE_MAX_ENTRY_NBYTES,
        "max_mem_nbytes": QUERY_CACHE_MAX_MEM_NBYTES,
        "max_persistent_nbytes": QUERY_CACHE_MAX_PERSISTENT_NBYTES,
        "persistent_nbytes": 0,
        "next_slot": 0,
        "entries": {},
    }


def _normalize_query_cache_catalog(catalog: dict) -> dict:
    """Ensure the prototype query-cache catalog has the current nbytes schema."""
    if not isinstance(catalog, dict):
        return _default_query_cache_catalog("")
    catalog.setdefault("version", QUERY_CACHE_FORMAT_VERSION)
    catalog.setdefault("payload_ref", {"kind": "urlpath", "version": 1, "urlpath": ""})
    catalog.setdefault("max_entry_nbytes", QUERY_CACHE_MAX_ENTRY_NBYTES)
    catalog.setdefault("max_mem_nbytes", QUERY_CACHE_MAX_MEM_NBYTES)
    catalog.setdefault("max_persistent_nbytes", QUERY_CACHE_MAX_PERSISTENT_NBYTES)
    catalog.setdefault("persistent_nbytes", 0)
    catalog.setdefault("next_slot", 0)
    catalog.setdefault("entries", {})
    return catalog


def _load_query_cache_catalog(array: blosc2.NDArray) -> dict | None:
    """Return ``None`` because query caches are intentionally not persisted."""
    return None


def _save_query_cache_catalog(array: blosc2.NDArray, catalog: dict) -> None:
    """No-op: query caches are intentionally not persisted."""
    return


def _open_query_cache_store(array: blosc2.NDArray, *, create: bool = False):
    """Return ``None`` because query caches are intentionally not persisted."""
    return


def _close_query_cache_store(path: str) -> None:
    """Drop a cached ObjectArray handle for *path*."""
    _QUERY_CACHE_STORE_HANDLES.pop(path, None)


# ---------------------------------------------------------------------------
# Stage 2 – Cache key normalization
# ---------------------------------------------------------------------------


def _normalize_query_descriptor(
    expression: str,
    tokens: list[str],
    order: list[str] | None,
) -> dict:
    """Build a canonical, order-stable query descriptor for cache keying."""
    try:
        normalized_expr = ast.unparse(ast.parse(expression, mode="eval"))
    except Exception:
        normalized_expr = expression
    return {
        "version": QUERY_CACHE_FORMAT_VERSION,
        "kind": "indices",
        "tokens": sorted(tokens),
        "expr": normalized_expr,
        "order": list(order) if order is not None else None,
    }


def _query_cache_digest(descriptor: dict) -> str:
    """Return a 32-character hex digest for *descriptor*."""
    import json

    canonical = json.dumps(descriptor, sort_keys=True, separators=(",", ":"))
    return hashlib.blake2b(canonical.encode(), digest_size=16).hexdigest()


# ---------------------------------------------------------------------------
# Stage 3 – Payload encode/decode and hot/persistent cache helpers
# ---------------------------------------------------------------------------


def _encode_coords_payload(coords: np.ndarray) -> dict:
    """Encode a coordinate array as a compact msgpack-safe mapping."""
    if coords.size == 0:
        dtype = np.dtype("<u4")
    else:
        dtype = np.dtype("<u4") if coords.max() <= np.iinfo(np.uint32).max else np.dtype("<u8")
    return {
        "version": QUERY_CACHE_FORMAT_VERSION,
        "dtype": dtype.str,
        "nrows": len(coords),
        "data": coords.astype(dtype).tobytes(),
    }


def _decode_coords_payload(payload: dict) -> np.ndarray:
    """Reconstruct a coordinate array from a cached payload mapping."""
    return np.frombuffer(payload["data"], dtype=np.dtype(payload["dtype"])).copy()


def _hot_cache_key(
    digest: str, scope: tuple[str, str | int] | None = None
) -> tuple[tuple[str, str | int], str]:
    return (_HOT_CACHE_GLOBAL_SCOPE if scope is None else scope, digest)


def _compress_hot_coords(coords: np.ndarray) -> _CompressedHotCoords:
    payload = _encode_coords_payload(np.asarray(coords))
    raw = payload["data"]
    compressed = False
    data = raw
    if len(raw) != 0:
        dtype = np.dtype(payload["dtype"])
        candidate = blosc2.compress2(raw, typesize=dtype.itemsize, codec=blosc2.Codec.LZ4, clevel=5)
        if len(candidate) < len(raw):
            data = candidate
            compressed = True
    return _CompressedHotCoords(
        dtype=payload["dtype"], nrows=int(payload["nrows"]), compressed=compressed, data=data
    )


def _decompress_hot_coords(entry: _CompressedHotCoords) -> np.ndarray:
    dtype = np.dtype(entry.dtype)
    if entry.nrows == 0:
        return np.empty((0,), dtype=dtype)
    raw = blosc2.decompress2(entry.data) if entry.compressed else entry.data
    return np.frombuffer(raw, dtype=dtype, count=entry.nrows).copy()


def _hot_cache_get(digest: str, scope: tuple[str, str | int] | None = None) -> np.ndarray | None:
    """Return the cached coordinate array for *digest*, or ``None``."""
    key = _hot_cache_key(digest, scope)
    entry = _HOT_CACHE.get(key)
    if entry is None:
        return None
    # Move to most-recently-used position.
    with contextlib.suppress(ValueError):
        _HOT_CACHE_ORDER.remove(key)
    _HOT_CACHE_ORDER.append(key)
    return _decompress_hot_coords(entry)


def _hot_cache_put(digest: str, coords: np.ndarray, scope: tuple[str, str | int] | None = None) -> None:
    """Insert *coords* into the hot cache, evicting LRU entries if needed."""
    global _HOT_CACHE_BYTES
    key = _hot_cache_key(digest, scope)
    entry = _compress_hot_coords(coords)
    entry_bytes = entry.nbytes
    if entry_bytes > QUERY_CACHE_MAX_MEM_NBYTES:
        # Single entry too large; skip.
        return
    # If already present, remove old accounting first.
    if key in _HOT_CACHE:
        _HOT_CACHE_BYTES -= _HOT_CACHE[key].nbytes
        with contextlib.suppress(ValueError):
            _HOT_CACHE_ORDER.remove(key)
    # Evict LRU entries until there is room.
    while _HOT_CACHE_ORDER and _HOT_CACHE_BYTES + entry_bytes > QUERY_CACHE_MAX_MEM_NBYTES:
        oldest = _HOT_CACHE_ORDER.pop(0)
        evicted = _HOT_CACHE.pop(oldest, None)
        if evicted is not None:
            _HOT_CACHE_BYTES -= evicted.nbytes
    _HOT_CACHE[key] = entry
    _HOT_CACHE_ORDER.append(key)
    _HOT_CACHE_BYTES += entry_bytes


def _hot_cache_clear(scope: tuple[str, str | int] | None = None) -> None:
    """Clear all in-process hot cache entries for *scope* (or all scopes)."""
    global _HOT_CACHE_BYTES
    if scope is not None:
        keys = [key for key in _HOT_CACHE if key[0] == scope]
        for key in keys:
            _HOT_CACHE_BYTES -= _HOT_CACHE.pop(key).nbytes
        _HOT_CACHE_ORDER[:] = [key for key in _HOT_CACHE_ORDER if key[0] != scope]
        return
    _HOT_CACHE.clear()
    _HOT_CACHE_ORDER.clear()
    _HOT_CACHE_BYTES = 0


def _persistent_cache_lookup(array: blosc2.NDArray, digest: str) -> np.ndarray | None:
    """Return ``None`` because query caches are intentionally not persisted."""
    return None


def _query_cache_entry_nbytes(coords: np.ndarray) -> int:
    """Return the logical int64 position bytes used for persistent budget accounting."""
    return int(np.asarray(coords).size) * np.dtype(np.int64).itemsize


def _reset_persistent_query_cache_catalog(array: blosc2.NDArray, catalog: dict | None = None) -> dict:
    """Drop persistent cache storage and return a fresh empty catalog preserving limits."""
    payload_path = _query_cache_payload_path(array)
    _close_query_cache_store(payload_path)
    blosc2.remove_urlpath(payload_path)

    fresh = _default_query_cache_catalog(payload_path)
    if catalog is not None:
        fresh["max_entry_nbytes"] = int(catalog.get("max_entry_nbytes", QUERY_CACHE_MAX_ENTRY_NBYTES))
        fresh["max_mem_nbytes"] = int(catalog.get("max_mem_nbytes", QUERY_CACHE_MAX_MEM_NBYTES))
        fresh["max_persistent_nbytes"] = int(
            catalog.get("max_persistent_nbytes", QUERY_CACHE_MAX_PERSISTENT_NBYTES)
        )
    _save_query_cache_catalog(array, fresh)
    return fresh


def _persistent_cache_insert(
    array: blosc2.NDArray,
    digest: str,
    coords: np.ndarray,
    query_descriptor: dict,
) -> bool:
    """Return ``False`` because query caches are intentionally not persisted."""
    return False


# ---------------------------------------------------------------------------
# Stage 5 – Query cache invalidation
# ---------------------------------------------------------------------------


def _invalidate_query_cache(array: blosc2.NDArray) -> None:
    """Drop the entire query cache for *array* (persistent file + hot cache)."""
    scope = _query_cache_scope(array)
    if not _is_persistent_array(array):
        _hot_cache_clear(scope=scope)
        return
    payload_path = _query_cache_payload_path(array)
    _close_query_cache_store(payload_path)
    blosc2.remove_urlpath(payload_path)
    with contextlib.suppress(KeyError, Exception):
        del array.schunk.vlmeta[QUERY_CACHE_VLMETA_KEY]
    _hot_cache_clear(scope=scope)
    # Drop any cached mmap handle for this array's data file so a re-opened or
    # extended array is not served from a stale mapping.
    urlpath = getattr(array, "urlpath", None)
    if urlpath is not None:
        _GATHER_MMAP_HANDLES.pop(str(urlpath), None)


# ---------------------------------------------------------------------------
# Public helper: cached coordinate lookup (used by lazyexpr.py integration)
# ---------------------------------------------------------------------------


def get_cached_coords(
    array: blosc2.NDArray,
    expression: str,
    tokens: list[str],
    order: list[str] | None,
) -> np.ndarray | None:
    """Return cached coordinates for *expression*/*tokens*/*order*, or ``None``."""
    owner = _query_cache_owner(array)
    scope = _query_cache_scope(owner)
    # Drop stale coordinates if the underlying file was overwritten in place.
    if _is_persistent_array(owner):
        _refresh_persistent_caches(getattr(owner, "urlpath", None))
    descriptor = _normalize_query_descriptor(expression, tokens, order)
    digest = _query_cache_digest(descriptor)
    return _hot_cache_get(digest, scope=scope)


def store_cached_coords(
    array: blosc2.NDArray,
    expression: str,
    tokens: list[str],
    order: list[str] | None,
    coords: np.ndarray,
) -> None:
    """Store *coords* in the in-process hot cache only."""
    owner = _query_cache_owner(array)
    scope = _query_cache_scope(owner)
    descriptor = _normalize_query_descriptor(expression, tokens, order)
    digest = _query_cache_digest(descriptor)
    _hot_cache_put(digest, coords, scope=scope)


def _supported_index_dtype(dtype: np.dtype) -> bool:
    return np.dtype(dtype).kind in {"b", "i", "u", "f", "m", "M", "S", "U"}


def _field_target_descriptor(field: str | None) -> dict:
    return {"source": "field", "field": field}


def _expression_target_descriptor(expression: str, expression_key: str, dependencies: list[str]) -> dict:
    return {
        "source": "expression",
        "expression": expression,
        "expression_key": expression_key,
        "dependencies": list(dependencies),
    }


def _target_field(target: dict) -> str | None:
    return target.get("field") if target.get("source") == "field" else None


def _field_dtype(array: blosc2.NDArray, field: str | None) -> np.dtype:
    if field is None:
        return np.dtype(array.dtype)
    if array.dtype.fields is None:
        raise TypeError("field indexes require a structured dtype")
    if field not in array.dtype.fields:
        raise ValueError(f"field {field!r} is not present in the dtype")
    return np.dtype(array.dtype.fields[field][0])


def _validate_index_target(array: blosc2.NDArray, field: str | None) -> np.dtype:
    if not isinstance(array, blosc2.NDArray):
        raise TypeError("indexes are only supported on NDArray")
    if array.ndim != 1:
        raise ValueError("indexes are only supported on 1-D NDArray objects")
    dtype = _field_dtype(array, field)
    if not _supported_index_dtype(dtype):
        raise TypeError(f"dtype {dtype} is not supported by the current index engine")
    return dtype


def _normalize_index_kind(kind: blosc2.IndexKind | str) -> str:
    if isinstance(kind, enum.Enum):
        kind = kind.value
    if kind not in SEGMENT_LEVELS_BY_KIND:
        raise NotImplementedError(f"unsupported index kind {kind!r}")
    return kind


def _normalize_build_mode(build: str) -> str:
    if build not in {"auto", "memory", "ooc"}:
        raise ValueError("build must be one of 'auto', 'memory', or 'ooc'")
    return build


def _normalize_full_build_method(method: str | None) -> str:
    if method is None:
        return "global-sort"
    if method != "global-sort":
        raise ValueError("full index method must be 'global-sort'; use kind=IndexKind.OPSI for OPSI indexes")
    return method


def _field_name_exists(array: blosc2.NDArray, field: str | None) -> bool:
    return field is not None and array.dtype.fields is not None and field in array.dtype.fields


def _normalize_create_index_target(
    array: blosc2.NDArray,
    field: str | None,
    expression: str | None,
    operands: dict | None,
) -> tuple[dict, np.dtype]:
    if field is not None and expression is not None:
        raise ValueError("field and expression are mutually exclusive")

    if expression is not None:
        if operands is None:
            operands = array.fields if array.dtype.fields is not None else {"value": array}
        base, target, dtype = _normalize_expression_target(expression, operands)
        if base is not array:
            raise ValueError(
                "expression index operands must resolve to the same array passed to create_index()"
            )
        return target, dtype

    if operands is not None:
        raise ValueError("operands can only be provided together with expression")

    if _field_name_exists(array, field):
        return _field_target_descriptor(field), _validate_index_target(array, field)

    if field is not None:
        base_operands = array.fields if array.dtype.fields is not None else {"value": array}
        base, target, dtype = _normalize_expression_target(field, base_operands)
        if base is not array:
            raise ValueError(
                "expression index operands must resolve to the same array passed to create_index()"
            )
        return target, dtype

    return _field_target_descriptor(None), _validate_index_target(array, None)


class _OperandCanonicalizer(ast.NodeTransformer):
    def __init__(self, operands: dict):
        self.operands = operands
        self.base: blosc2.NDArray | None = None
        self.dependencies: list[str] = []
        self.valid = True

    def visit_Name(self, node: ast.Name) -> ast.AST:
        operand = self.operands.get(node.id)
        if operand is None:
            return node
        target = _operand_target(operand)
        if target is None:
            self.valid = False
            return node
        base, field = target
        if self.base is None:
            self.base = base
        elif self.base is not base:
            self.valid = False
            return node
        canonical = SELF_TARGET_NAME if field is None else field
        self.dependencies.append(canonical)
        return ast.copy_location(ast.Name(id=canonical, ctx=node.ctx), node)


def _normalize_expression_node(
    node: ast.AST, operands: dict
) -> tuple[blosc2.NDArray, str, list[str]] | None:
    canonicalizer = _OperandCanonicalizer(operands)
    normalized = canonicalizer.visit(
        ast.fix_missing_locations(ast.parse(ast.unparse(node), mode="eval")).body
    )
    if not canonicalizer.valid or canonicalizer.base is None or not canonicalizer.dependencies:
        return None
    dependencies = list(dict.fromkeys(canonicalizer.dependencies))
    return canonicalizer.base, ast.unparse(normalized), dependencies


def _normalize_expression_target(expression: str, operands: dict) -> tuple[blosc2.NDArray, dict, np.dtype]:
    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError as exc:
        raise ValueError("expression is not valid Python syntax") from exc

    normalized = _normalize_expression_node(tree.body, operands)
    if normalized is None:
        raise ValueError("expression indexes require operands from a single 1-D NDArray target")
    base, expression_key, dependencies = normalized
    if base.ndim != 1:
        raise ValueError("expression indexes are only supported on 1-D NDArray objects")
    target = _expression_target_descriptor(expression, expression_key, dependencies)
    sample_stop = min(int(base.shape[0]), max(1, int(base.blocks[0]) if base.blocks else 1))
    sample = _slice_values_for_target(base, target, 0, sample_stop)
    dtype = np.dtype(sample.dtype)
    if sample.ndim != 1:
        raise ValueError("expression indexes require expressions returning a 1-D scalar stream")
    if not _supported_index_dtype(dtype):
        raise TypeError(f"dtype {dtype} is not supported by the current index engine")
    return base, target, dtype


def _sanitize_sidecar_root(urlpath: str | Path) -> tuple[Path, str]:
    path = Path(urlpath)
    suffix = "".join(path.suffixes)
    root = path.name[: -len(suffix)] if suffix else path.name
    return path, root


def _sidecar_path(array: blosc2.NDArray, token: str, kind: str, name: str) -> str:
    path, root = _sanitize_sidecar_root(array.urlpath)
    # Drop redundant kind prefix when category == kind (e.g. "summary.block" under kind="summary").
    if name.startswith(f"{kind}."):
        name = name[len(kind) + 1 :]
    # Omit __self__ token: it is the common case (self-indexed column) and adds no information.
    token_part = f".{_sanitize_token(token)}" if token != SELF_TARGET_NAME else ""
    # CTable index anchors use "_anchor" as a placeholder; the _indexes/{col}/ directory
    # already identifies the column, so no prefix is needed in the filename.
    # For standalone NDArray indexes the root.__index__ prefix avoids collisions with sibling files.
    if root == "_anchor":
        return str(path.parent / f"{kind}{token_part}.{name}.b2nd")
    return str(path.parent / f"{root}.__index__{token_part}.{kind}.{name}.b2nd")


def _segment_len(array: blosc2.NDArray, level: str) -> int:
    if level == "chunk":
        return int(array.chunks[0])
    if level == "block":
        return int(array.blocks[0])
    if level == "subblock":
        return max(1, int(array.blocks[0]) // 8)
    raise ValueError(f"unknown level {level!r}")


def _data_cache_key(array: blosc2.NDArray, token: str, category: str, name: str):
    return (_array_key(array), token, category, name)


def _clear_cached_data(array: blosc2.NDArray, token: str) -> None:
    prefix = (_array_key(array), token)
    keys = [key for key in _DATA_CACHE if key[:2] == prefix]
    for key in keys:
        _DATA_CACHE.pop(key, None)
    handle_keys = [key for key in _SIDECAR_HANDLE_CACHE if key[:2] == prefix]
    for key in handle_keys:
        _SIDECAR_HANDLE_CACHE.pop(key, None)


def _sidecar_handle_cache_key(
    array: blosc2.NDArray, token: str, category: str, name: str, path: str | None = None
):
    return (_array_key(array), token, category, name, path)


def _sidecar_storage_category(category: str) -> str:
    return category.removesuffix("_handle")


def _invalidate_sidecar_cache_entries(array: blosc2.NDArray, token: str, category: str, name: str) -> None:
    """Drop cached data/handles for a sidecar before replacing its storage."""
    storage_category = _sidecar_storage_category(category)
    categories = {storage_category, f"{storage_category}_handle"}
    for cache_category in categories:
        _DATA_CACHE.pop(_data_cache_key(array, token, cache_category, name), None)
        prefix = _sidecar_handle_cache_key(array, token, cache_category, name)
        for key in [k for k in _SIDECAR_HANDLE_CACHE if k[:4] == prefix[:4]]:
            _SIDECAR_HANDLE_CACHE.pop(key, None)


def _open_sidecar_handle(array: blosc2.NDArray, token: str, category: str, name: str, path: str | None):
    _purge_stale_persistent_caches()
    cache_key = _sidecar_handle_cache_key(array, token, category, name, path)
    cached = _SIDECAR_HANDLE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    if path is None:
        storage_category = _sidecar_storage_category(category)
        legacy = _SIDECAR_HANDLE_CACHE.get(_sidecar_handle_cache_key(array, token, storage_category, name))
        if legacy is not None:
            _SIDECAR_HANDLE_CACHE[cache_key] = legacy
            return legacy
        legacy = _DATA_CACHE.get(_data_cache_key(array, token, storage_category, name))
        if legacy is None:
            raise RuntimeError("sidecar handle path is not available")
        handle = legacy if isinstance(legacy, blosc2.NDArray) else blosc2.asarray(np.asarray(legacy))
    else:
        handle = _open_sidecar_file(path, _INDEX_MMAP_MODE)
    _SIDECAR_HANDLE_CACHE[cache_key] = handle
    return handle


def _operands_for_dependencies(values: np.ndarray, dependencies: list[str]) -> dict[str, np.ndarray]:
    operands = {}
    for dependency in dependencies:
        if dependency == SELF_TARGET_NAME:
            operands[dependency] = values
        else:
            operands[dependency] = values[dependency]
    return operands


def _values_from_numpy_target(values: np.ndarray, target: dict) -> np.ndarray:
    if target["source"] == "field":
        field = target.get("field")
        return values if field is None else values[field]
    if target["source"] == "expression":
        from .lazyexpr import ne_evaluate

        result = ne_evaluate(
            target["expression_key"], _operands_for_dependencies(values, target["dependencies"])
        )
        return np.asarray(result)
    raise ValueError(f"unsupported index target source {target['source']!r}")


def _values_for_target(array: blosc2.NDArray, target: dict) -> np.ndarray:
    return _slice_values_for_target(array, target, 0, int(array.shape[0]))


def _slice_values_for_target(array: blosc2.NDArray, target: dict, start: int, stop: int) -> np.ndarray:
    return _values_from_numpy_target(array[start:stop], target)


def _summary_dtype(dtype: np.dtype) -> np.dtype:
    return np.dtype([("min", dtype), ("max", dtype), ("flags", np.uint8)])


def _boundary_dtype(dtype: np.dtype) -> np.dtype:
    return np.dtype([("start", dtype), ("end", dtype)])


def _segment_summary(segment: np.ndarray, dtype: np.dtype):
    flags = np.uint8(0)
    if dtype.kind == "f":
        valid = ~np.isnan(segment)
        if not np.all(valid):
            flags |= FLAG_HAS_NAN
        if not np.any(valid):
            flags |= FLAG_ALL_NAN
            zero = np.zeros((), dtype=dtype)[()]
            return zero, zero, flags
        segment = segment[valid]
    if dtype.kind in "US":
        # String dtypes: ufunc 'minimum'/'maximum' lack a loop.
        mn = segment[0]
        mx = segment[0]
        for v in segment[1:]:
            if v < mn:
                mn = v
            if v > mx:
                mx = v
        return mn, mx, flags
    return segment.min(), segment.max(), flags


def _compute_segment_summaries(values: np.ndarray, dtype: np.dtype, segment_len: int) -> np.ndarray:
    nsegments = math.ceil(values.shape[0] / segment_len)
    summary_dtype = _summary_dtype(dtype)
    summaries = np.empty(nsegments, dtype=summary_dtype)

    for idx in range(nsegments):
        start = idx * segment_len
        stop = min(start + segment_len, values.shape[0])
        segment = values[start:stop]
        summaries[idx] = _segment_summary(segment, dtype)
    return summaries


def _fill_summaries_from_2d(
    data_2d: np.ndarray,
    summaries_arr: np.ndarray,
    offset: int,
    dtype: np.dtype,
) -> None:
    """Fill summaries_arr[offset:offset+n] from data_2d (shape n×segment_len) with vectorized ops."""
    n = data_2d.shape[0]
    if n == 0:
        return
    if dtype.kind == "f":
        # All-NaN blocks make np.nanmin/nanmax emit "All-NaN slice encountered";
        # their results are immediately overwritten with zero below, so silence
        # the (purely cosmetic) RuntimeWarning.
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"All-NaN slice encountered", RuntimeWarning)
            has_nan = np.any(np.isnan(data_2d), axis=1)
            all_nan = np.all(np.isnan(data_2d), axis=1)
            mins = np.nanmin(data_2d, axis=1)
            maxs = np.nanmax(data_2d, axis=1)
        flags = np.where(has_nan, FLAG_HAS_NAN, np.uint8(0)).astype(np.uint8)
        flags = np.where(all_nan, np.uint8(FLAG_ALL_NAN | FLAG_HAS_NAN), flags)
        zero = dtype.type(0)
        mins = np.where(all_nan, zero, mins).astype(dtype)
        maxs = np.where(all_nan, zero, maxs).astype(dtype)
    else:
        if dtype.kind in "US":
            # String dtypes: numpy ufunc 'minimum'/'maximum' lack a loop for <U/S.
            # Use manual per-row comparison (cheap for small segment_len).
            mins = np.empty(n, dtype=dtype)
            maxs = np.empty(n, dtype=dtype)
            for i in range(n):
                row = data_2d[i]
                mn = row[0]
                mx = row[0]
                for v in row[1:]:
                    if v < mn:
                        mn = v
                    if v > mx:
                        mx = v
                mins[i] = mn
                maxs[i] = mx
        else:
            mins = data_2d.min(axis=1)
            maxs = data_2d.max(axis=1)
        flags = np.zeros(n, dtype=np.uint8)
    summaries_arr["min"][offset : offset + n] = mins
    summaries_arr["max"][offset : offset + n] = maxs
    summaries_arr["flags"][offset : offset + n] = flags


def _sorted_segment_boundary(segment: np.ndarray, dtype: np.dtype):
    if dtype.kind == "f":
        valid = ~np.isnan(segment)
        if not np.any(valid):
            nan = np.asarray(np.nan, dtype=dtype)[()]
            return nan, nan
        segment = segment[valid]
    return segment[0], segment[-1]


def _compute_sorted_boundaries(values: np.ndarray, dtype: np.dtype, segment_len: int) -> np.ndarray:
    nsegments = math.ceil(values.shape[0] / segment_len)
    boundaries = np.empty(nsegments, dtype=_boundary_dtype(dtype))

    for idx in range(nsegments):
        start = idx * segment_len
        stop = min(start + segment_len, values.shape[0])
        segment = values[start:stop]
        boundaries[idx] = _sorted_segment_boundary(segment, dtype)
    return boundaries


def _compute_sorted_boundaries_from_sidecar(
    path: str, dtype: np.dtype, length: int, segment_len: int
) -> np.ndarray:
    nsegments = math.ceil(length / segment_len)
    boundaries = np.empty(nsegments, dtype=_boundary_dtype(dtype))
    sidecar = _open_sidecar_file(path, _INDEX_MMAP_MODE)
    start_value = np.empty(1, dtype=dtype)
    end_value = np.empty(1, dtype=dtype)
    for idx in range(nsegments):
        start = idx * segment_len
        stop = min(start + segment_len, length)
        _read_ndarray_linear_span(sidecar, start, start_value)
        _read_ndarray_linear_span(sidecar, stop - 1, end_value)
        boundaries[idx] = (start_value[0], end_value[0])
    return boundaries


def _compute_sorted_boundaries_from_handle(
    handle, dtype: np.dtype, length: int, segment_len: int
) -> np.ndarray:
    nsegments = math.ceil(length / segment_len)
    boundaries = np.empty(nsegments, dtype=_boundary_dtype(dtype))
    start_value = np.empty(1, dtype=dtype)
    end_value = np.empty(1, dtype=dtype)
    for idx in range(nsegments):
        start = idx * segment_len
        stop = min(start + segment_len, length)
        _read_ndarray_linear_span(handle, start, start_value)
        _read_ndarray_linear_span(handle, stop - 1, end_value)
        boundaries[idx] = (start_value[0], end_value[0])
    return boundaries


def _store_array_sidecar(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    category: str,
    name: str,
    data: np.ndarray,
    persistent: bool,
    *,
    chunks: tuple[int, ...] | None = None,
    blocks: tuple[int, ...] | None = None,
    cparams: dict | None = None,
) -> dict:
    cache_key = _data_cache_key(array, token, category, name)
    handle_cache_key = _sidecar_handle_cache_key(array, token, category, name)
    _invalidate_sidecar_cache_entries(array, token, category, name)
    if persistent:
        path = _sidecar_path(array, token, kind, f"{category}.{name}")
        blosc2.remove_urlpath(path)
        kwargs = {"urlpath": path, "mode": "w"}
        if chunks is not None:
            kwargs["chunks"] = chunks
        if blocks is not None:
            kwargs["blocks"] = blocks
        if cparams is not None:
            kwargs["cparams"] = cparams
        # Do not retain writable persistent handles in the process-wide cache.
        # They keep native resources alive after index construction and can
        # accumulate badly across tests on macOS/Python 3.14.
        handle = blosc2.asarray(data, **kwargs)
        del handle
        _DATA_CACHE.pop(cache_key, None)
    else:
        path = None
        kwargs = {}
        if chunks is not None:
            kwargs["chunks"] = chunks
        if blocks is not None:
            kwargs["blocks"] = blocks
        if cparams is not None:
            kwargs["cparams"] = cparams
        handle_data = np.array(data, copy=True) if isinstance(data, np.memmap) else data
        handle = blosc2.asarray(handle_data, **kwargs)
        _SIDECAR_HANDLE_CACHE[handle_cache_key] = handle
        _DATA_CACHE.pop(cache_key, None)
    return {"path": path, "dtype": data.dtype.descr if data.dtype.fields else data.dtype.str}


def _create_persistent_sidecar_handle(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    category: str,
    name: str,
    length: int,
    dtype: np.dtype,
    *,
    chunks: tuple[int, ...] | None = None,
    blocks: tuple[int, ...] | None = None,
    cparams: dict | None = None,
) -> tuple[blosc2.NDArray | None, dict]:
    _invalidate_sidecar_cache_entries(array, token, category, name)
    path = _sidecar_path(array, token, kind, f"{category}.{name}")
    blosc2.remove_urlpath(path)
    kwargs = {"urlpath": path, "mode": "w"}
    if chunks is not None:
        kwargs["chunks"] = chunks
    if blocks is not None:
        kwargs["blocks"] = blocks
    if cparams is not None:
        kwargs["cparams"] = cparams
    if length == 0:
        handle = blosc2.asarray(np.empty(0, dtype=dtype), **kwargs)
        del handle
        return None, {"path": path, "dtype": dtype.descr if dtype.fields else dtype.str}
    handle = blosc2.empty((length,), dtype=dtype, **kwargs)
    return handle, {"path": path, "dtype": dtype.descr if dtype.fields else dtype.str}


def _normalize_index_cparams(cparams) -> blosc2.CParams | None:
    if cparams is None:
        return None
    if isinstance(cparams, blosc2.CParams):
        return cparams
    return blosc2.CParams(**cparams)


def _plain_index_cparams(cparams: dict | blosc2.CParams | None) -> dict | None:
    if cparams is None:
        return None

    def _plain_value(value):
        if isinstance(value, enum.Enum):
            return value.value
        if isinstance(value, dict):
            return {key: _plain_value(item) for key, item in value.items()}
        if isinstance(value, list | tuple):
            return type(value)(_plain_value(item) for item in value)
        return value

    if isinstance(cparams, blosc2.CParams):
        cparams = asdict(cparams)
    else:
        cparams = cparams.copy()
    return {key: _plain_value(value) for key, value in cparams.items()}


def _load_array_sidecar(
    array: blosc2.NDArray, token: str, category: str, name: str, path: str | None = None
) -> np.ndarray:
    cache_key = _data_cache_key(array, token, category, name)
    cached = _DATA_CACHE.get(cache_key)
    if cached is not None:
        return cached
    handle = _SIDECAR_HANDLE_CACHE.get(_sidecar_handle_cache_key(array, token, category, name, path))
    if handle is None:
        raise RuntimeError("in-memory index metadata is missing from the current process")
    data = _read_sidecar_span(handle, 0, int(handle.shape[0]))
    _DATA_CACHE[cache_key] = data
    return data


def _build_levels_descriptor(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    values: np.ndarray,
    persistent: bool,
    cparams: dict | None = None,
    summary_levels: tuple[str, ...] | None = None,
) -> dict:
    levels = {}
    levels_to_build = summary_levels if summary_levels is not None else SEGMENT_LEVELS_BY_KIND[kind]
    for level in levels_to_build:
        segment_len = _segment_len(array, level)
        summaries = _compute_segment_summaries(values, dtype, segment_len)
        sidecar = _store_array_sidecar(
            array, token, kind, "summary", level, summaries, persistent, cparams=cparams
        )
        levels[level] = {
            "segment_len": segment_len,
            "nsegments": len(summaries),
            "path": sidecar["path"],
            "dtype": sidecar["dtype"],
        }
    return levels


def _build_levels_descriptor_ooc(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    persistent: bool,
    cparams: dict | None = None,
    summary_levels: tuple[str, ...] | None = None,
    precomputed_summaries: dict[str, np.ndarray] | None = None,
) -> dict:
    size = int(array.shape[0])
    summary_dtype = _summary_dtype(dtype)
    chunk_len = int(array.chunks[0])
    levels_to_build = summary_levels if summary_levels is not None else SEGMENT_LEVELS_BY_KIND[kind]
    segment_lens = {level: _segment_len(array, level) for level in levels_to_build}
    nsegments_total = {level: math.ceil(size / slen) for level, slen in segment_lens.items()}
    all_summaries = {level: np.empty(n, dtype=summary_dtype) for level, n in nsegments_total.items()}

    # Incremental fast path: summaries already accumulated during the write phase
    # (e.g. CTable import folds per-block min/max as data is appended, when it is
    # still uncompressed in memory).  Using them avoids decompressing the whole
    # column back just to recompute min/max.  Only trusted when every requested
    # level is present with the exact expected segment count and dtype; otherwise
    # fall through to the decompression pass below.
    use_precomputed = precomputed_summaries is not None and all(
        level in precomputed_summaries
        and precomputed_summaries[level].dtype == summary_dtype
        and len(precomputed_summaries[level]) == nsegments_total[level]
        for level in levels_to_build
    )
    if use_precomputed:
        for level in levels_to_build:
            all_summaries[level] = np.ascontiguousarray(precomputed_summaries[level])

    # Fast path: all segment sizes are ≤ chunk_len and divide it evenly, so no segment
    # spans a chunk boundary.  A single decompression pass over the data suffices.
    can_fast = all(slen <= chunk_len and chunk_len % slen == 0 for slen in segment_lens.values())

    if use_precomputed:
        pass
    elif can_fast:
        seg_offsets = dict.fromkeys(levels_to_build, 0)
        nchunks = math.ceil(size / chunk_len)
        for chunk_id in range(nchunks):
            chunk_start = chunk_id * chunk_len
            chunk_stop = min(chunk_start + chunk_len, size)
            chunk_values = _slice_values_for_target(array, target, chunk_start, chunk_stop)
            chunk_size = chunk_stop - chunk_start
            for level in levels_to_build:
                slen = segment_lens[level]
                summaries_arr = all_summaries[level]
                offset = seg_offsets[level]
                n_complete = chunk_size // slen
                remainder = chunk_size % slen
                if n_complete > 0:
                    data_2d = chunk_values[: n_complete * slen].reshape(n_complete, slen)
                    _fill_summaries_from_2d(data_2d, summaries_arr, offset, dtype)
                if remainder > 0:
                    summaries_arr[offset + n_complete] = _segment_summary(
                        chunk_values[n_complete * slen :], dtype
                    )
                    seg_offsets[level] = offset + n_complete + 1
                else:
                    seg_offsets[level] = offset + n_complete
    else:
        # Fallback: original segment-by-segment approach
        for level in levels_to_build:
            slen = segment_lens[level]
            for idx in range(nsegments_total[level]):
                start = idx * slen
                stop = min(start + slen, size)
                all_summaries[level][idx] = _segment_summary(
                    _slice_values_for_target(array, target, start, stop), dtype
                )

    levels = {}
    for level in levels_to_build:
        sidecar = _store_array_sidecar(
            array, token, kind, "summary", level, all_summaries[level], persistent, cparams=cparams
        )
        levels[level] = {
            "segment_len": segment_lens[level],
            "nsegments": nsegments_total[level],
            "path": sidecar["path"],
            "dtype": sidecar["dtype"],
        }
    return levels


def _sidecar_storage_geometry(
    path: str | None, fallback_chunk_len: int, fallback_block_len: int
) -> tuple[int, int]:
    if path is None:
        return fallback_chunk_len, fallback_block_len
    sidecar = _open_sidecar_file(path, _INDEX_MMAP_MODE)
    return int(sidecar.chunks[0]), int(sidecar.blocks[0])


def _rebuild_full_navigation_sidecars(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    full: dict,
    sorted_values: np.ndarray,
    persistent: bool,
    cparams: dict | None = None,
) -> None:
    chunk_len, block_len = _sidecar_storage_geometry(
        full.get("values_path"),
        int(full.get("sidecar_chunk_len", array.chunks[0])),
        int(full.get("sidecar_block_len", array.blocks[0])),
    )
    l1 = _compute_sorted_boundaries(sorted_values, np.dtype(sorted_values.dtype), chunk_len)
    l2 = _compute_sorted_boundaries(sorted_values, np.dtype(sorted_values.dtype), block_len)
    l1_sidecar = _store_array_sidecar(array, token, kind, "full_nav", "l1", l1, persistent, cparams=cparams)
    l2_sidecar = _store_array_sidecar(array, token, kind, "full_nav", "l2", l2, persistent, cparams=cparams)
    full["l1_path"] = l1_sidecar["path"]
    full["l2_path"] = l2_sidecar["path"]
    full["sidecar_chunk_len"] = int(chunk_len)
    full["sidecar_block_len"] = int(block_len)


def _rebuild_full_navigation_sidecars_from_path(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    full: dict,
    values_path: str,
    dtype: np.dtype,
    length: int,
    persistent: bool,
    cparams: dict | None = None,
) -> None:
    chunk_len, block_len = _sidecar_storage_geometry(values_path, int(array.chunks[0]), int(array.blocks[0]))
    l1 = _compute_sorted_boundaries_from_sidecar(values_path, dtype, length, chunk_len)
    l2 = _compute_sorted_boundaries_from_sidecar(values_path, dtype, length, block_len)
    l1_sidecar = _store_array_sidecar(array, token, kind, "full_nav", "l1", l1, persistent, cparams=cparams)
    l2_sidecar = _store_array_sidecar(array, token, kind, "full_nav", "l2", l2, persistent, cparams=cparams)
    full["l1_path"] = l1_sidecar["path"]
    full["l2_path"] = l2_sidecar["path"]
    full["sidecar_chunk_len"] = int(chunk_len)
    full["sidecar_block_len"] = int(block_len)
    full["l1_dtype"] = l1_sidecar["dtype"]
    full["l2_dtype"] = l2_sidecar["dtype"]


def _rebuild_full_navigation_sidecars_from_handle(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    full: dict,
    values_handle,
    dtype: np.dtype,
    length: int,
    persistent: bool,
    cparams: dict | None = None,
) -> None:
    chunk_len = int(values_handle.chunks[0]) if hasattr(values_handle, "chunks") else int(array.chunks[0])
    block_len = int(values_handle.blocks[0]) if hasattr(values_handle, "blocks") else int(array.blocks[0])
    l1 = _compute_sorted_boundaries_from_handle(values_handle, dtype, length, chunk_len)
    l2 = _compute_sorted_boundaries_from_handle(values_handle, dtype, length, block_len)
    l1_sidecar = _store_array_sidecar(array, token, kind, "full_nav", "l1", l1, persistent, cparams=cparams)
    l2_sidecar = _store_array_sidecar(array, token, kind, "full_nav", "l2", l2, persistent, cparams=cparams)
    full["l1_path"] = l1_sidecar["path"]
    full["l2_path"] = l2_sidecar["path"]
    full["sidecar_chunk_len"] = int(chunk_len)
    full["sidecar_block_len"] = int(block_len)
    full["l1_dtype"] = l1_sidecar["dtype"]
    full["l2_dtype"] = l2_sidecar["dtype"]


def _stream_copy_sidecar_array(
    source_path: Path | str,
    dest_path: Path | str,
    length: int,
    dtype: np.dtype,
    chunks: tuple[int, ...],
    blocks: tuple[int, ...],
    cparams: dict | None = None,
) -> None:
    source = blosc2.open(str(source_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)
    blosc2.remove_urlpath(str(dest_path))
    kwargs = {"chunks": chunks, "blocks": blocks, "urlpath": str(dest_path), "mode": "w"}
    if cparams is not None:
        kwargs["cparams"] = cparams
    dest = blosc2.empty((length,), dtype=dtype, **kwargs)
    chunk_len = int(dest.chunks[0])
    for start in range(0, length, chunk_len):
        stop = min(start + chunk_len, length)
        span = np.empty(stop - start, dtype=dtype)
        _read_ndarray_linear_span(source, start, span)
        dest[start:stop] = span
    del source, dest


def _full_sidecar_geometry(array: blosc2.NDArray, optlevel: int) -> tuple[tuple[int], tuple[int], int]:
    chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
    block_len = int(array.blocks[0])
    chunk_len = _opsi_storage_chunk_len(int(array.chunks[0]), block_len, chunk_multiplier)
    return (chunk_len,), (block_len,), chunk_multiplier


def _stream_copy_temp_run_to_full_sidecars(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    full: dict,
    run: SortedRun,
    dtype: np.dtype,
    persistent: bool,
    tracker: TempRunTracker | None = None,
    cparams: dict | None = None,
    optlevel: int = 5,
) -> None:
    if not persistent:
        raise ValueError("temp-run streaming only supports persistent runs")

    values_path = _sidecar_path(array, token, kind, "full.values")
    positions_path = _sidecar_path(array, token, kind, "full.positions")
    sidecar_chunks, sidecar_blocks, chunk_multiplier = _full_sidecar_geometry(array, optlevel)
    _remove_sidecar_path(values_path)
    _remove_sidecar_path(positions_path)
    _stream_copy_sidecar_array(
        run.values_path,
        values_path,
        run.length,
        dtype,
        sidecar_chunks,
        sidecar_blocks,
        cparams,
    )
    _stream_copy_sidecar_array(
        run.positions_path,
        positions_path,
        run.length,
        np.dtype(np.int64),
        sidecar_chunks,
        sidecar_blocks,
        cparams,
    )
    _tracker_register_delete(tracker, run.values_path, run.positions_path)
    run.values_path.unlink(missing_ok=True)
    run.positions_path.unlink(missing_ok=True)
    full["values_path"] = values_path
    full["positions_path"] = positions_path
    full["chunk_multiplier"] = chunk_multiplier
    full["runs"] = []
    full["next_run_id"] = 0
    _rebuild_full_navigation_sidecars_from_path(
        array, token, kind, full, values_path, dtype, run.length, persistent, cparams
    )


def _keysort_values_with_positions(values: np.ndarray, start: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Return values sorted with original/global positions carried in lockstep."""
    sorted_values = np.array(values, copy=True)
    positions = np.arange(start, start + sorted_values.size, dtype=np.int64)
    try:
        indexing_ext.keysort_values_positions(sorted_values, positions)
        return sorted_values, positions
    except (AttributeError, TypeError):
        order = np.argsort(values, kind="stable")
        return values[order], (order + start).astype(np.int64, copy=False)


def _build_full_descriptor(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    values: np.ndarray,
    persistent: bool,
    cparams: dict | None = None,
    optlevel: int = 5,
) -> dict:
    sorted_values, positions = _keysort_values_with_positions(values)
    sidecar_chunks, sidecar_blocks, chunk_multiplier = _full_sidecar_geometry(array, optlevel)
    values_sidecar = _store_array_sidecar(
        array,
        token,
        kind,
        "full",
        "values",
        sorted_values,
        persistent,
        chunks=sidecar_chunks,
        blocks=sidecar_blocks,
        cparams=cparams,
    )
    positions_sidecar = _store_array_sidecar(
        array,
        token,
        kind,
        "full",
        "positions",
        positions,
        persistent,
        chunks=sidecar_chunks,
        blocks=sidecar_blocks,
        cparams=cparams,
    )
    full = {
        "values_path": values_sidecar["path"],
        "positions_path": positions_sidecar["path"],
        "chunk_multiplier": chunk_multiplier,
        "sidecar_chunk_len": sidecar_chunks[0],
        "sidecar_block_len": sidecar_blocks[0],
        "runs": [],
        "next_run_id": 0,
    }
    _rebuild_full_navigation_sidecars(array, token, kind, full, sorted_values, persistent, cparams)
    return full


def _position_dtype(max_value: int) -> np.dtype:
    if max_value <= np.iinfo(np.uint8).max:
        return np.dtype(np.uint8)
    if max_value <= np.iinfo(np.uint16).max:
        return np.dtype(np.uint16)
    if max_value <= np.iinfo(np.uint32).max:
        return np.dtype(np.uint32)
    return np.dtype(np.uint64)


def _resolve_ooc_mode(kind: str, build: str) -> bool:
    if kind not in {"summary", "bucket", "partial", "full", "opsi"}:
        return False
    build = _normalize_build_mode(build)
    return build != "memory"


def _build_block_sorted_payload(
    values: np.ndarray, block_len: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.dtype]:
    nblocks = math.ceil(values.shape[0] / block_len)
    position_dtype = _position_dtype(block_len - 1)
    offsets = np.empty(nblocks + 1, dtype=np.int64)
    offsets[0] = 0
    sorted_values = np.empty_like(values)
    positions = np.empty(values.shape[0], dtype=position_dtype)
    cursor = 0

    for block_id in range(nblocks):
        start = block_id * block_len
        stop = min(start + block_len, values.shape[0])
        block = values[start:stop]
        order = np.argsort(block, kind="stable")
        block_size = stop - start
        next_cursor = cursor + block_size
        sorted_values[cursor:next_cursor] = block[order]
        positions[cursor:next_cursor] = order.astype(position_dtype, copy=False)
        cursor = next_cursor
        offsets[block_id + 1] = cursor

    return sorted_values, positions, offsets, position_dtype


def _build_partial_descriptor(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    values: np.ndarray,
    optlevel: int,
    persistent: bool,
    cparams: dict | None = None,
) -> dict:
    chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
    chunk_len = int(array.chunks[0]) * chunk_multiplier
    nav_segment_len, nav_segment_divisor = _partial_nav_segment_len(
        int(array.blocks[0]), chunk_len, optlevel
    )
    sorted_values, positions, offsets, l2, _ = _build_chunk_sorted_payload(
        values, chunk_len, nav_segment_len, cparams
    )
    l1 = _compute_sorted_boundaries(sorted_values, np.dtype(values.dtype), chunk_len)
    partial = _chunk_index_payload_storage(
        array,
        token,
        kind,
        "partial",
        "values",
        sorted_values,
        "positions",
        positions,
        offsets,
        l1,
        l2,
        persistent,
        chunk_len,
        nav_segment_len,
        cparams,
    )
    partial["position_dtype"] = positions.dtype.str
    partial["nav_segment_divisor"] = nav_segment_divisor
    partial["chunk_multiplier"] = chunk_multiplier
    return partial


def _segment_row_count(chunk_len: int, nav_segment_len: int) -> int:
    return max(1, math.ceil(chunk_len / nav_segment_len))


def _chunk_offsets(size: int, chunk_len: int) -> np.ndarray:
    nchunks = math.ceil(size / chunk_len)
    offsets = np.empty(nchunks + 1, dtype=np.int64)
    offsets[0] = 0
    if nchunks == 0:
        return offsets
    offsets[1:] = np.minimum(np.arange(1, nchunks + 1, dtype=np.int64) * chunk_len, size)
    return offsets


def _index_build_threads(cparams: dict | blosc2.CParams | None = None) -> int:
    if blosc2.IS_WASM:
        return 1
    forced = os.getenv("BLOSC2_INDEX_BUILD_THREADS")
    if forced is not None:
        try:
            forced_threads = int(forced)
        except ValueError:
            forced_threads = 1
        return _python_executor_threads(forced_threads)
    if cparams is not None:
        nthreads = cparams.nthreads if isinstance(cparams, blosc2.CParams) else cparams.get("nthreads")
    else:
        nthreads = None
    if nthreads is not None:
        try:
            cparams_threads = int(nthreads)
        except (TypeError, ValueError):
            cparams_threads = 1
        return _python_executor_threads(cparams_threads)
    return _python_executor_threads(int(getattr(blosc2, "nthreads", 1) or 1))


def _partial_nav_segment_divisor(optlevel: int) -> int:
    if optlevel <= 1:
        return 1
    if optlevel <= 3:
        return 2
    if optlevel <= 6:
        return 4
    if optlevel == 9:
        return 16
    return 8


def _partial_nav_segment_len(block_len: int, chunk_len: int, optlevel: int) -> tuple[int, int]:
    divisor = min(block_len, _partial_nav_segment_divisor(int(optlevel)))
    max_segments_per_chunk = 2048
    chunk_floor = max(1, math.ceil(int(chunk_len) / max_segments_per_chunk))
    return max(1, block_len // divisor, chunk_floor), divisor


def _build_chunk_sorted_payload(
    values: np.ndarray,
    chunk_len: int,
    nav_segment_len: int,
    cparams: dict | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.dtype]:
    size = values.shape[0]
    nchunks = math.ceil(size / chunk_len)
    position_dtype = _position_dtype(chunk_len - 1)
    offsets = np.empty(nchunks + 1, dtype=np.int64)
    offsets[0] = 0
    sorted_values = np.empty_like(values)
    positions = np.empty(size, dtype=position_dtype)
    l1 = np.empty(nchunks, dtype=_boundary_dtype(values.dtype))
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(values.dtype))

    cursor = 0
    thread_count = _index_build_threads(cparams)
    for chunk_id in range(nchunks):
        start = chunk_id * chunk_len
        stop = min(start + chunk_len, size)
        chunk = values[start:stop]
        chunk_size = stop - start
        next_cursor = cursor + chunk_size
        chunk_sorted, chunk_positions = _sort_chunk_intra_chunk(
            chunk, position_dtype, thread_count=thread_count
        )
        sorted_values[cursor:next_cursor] = chunk_sorted
        positions[cursor:next_cursor] = chunk_positions
        offsets[chunk_id + 1] = next_cursor
        l1[chunk_id] = _sorted_segment_boundary(chunk_sorted, np.dtype(values.dtype))

        row_start = chunk_id * nsegments_per_chunk
        segment_count = _segment_row_count(chunk_size, nav_segment_len)
        for segment_id in range(segment_count):
            seg_start = cursor + segment_id * nav_segment_len
            seg_stop = min(seg_start + nav_segment_len, next_cursor)
            l2[row_start + segment_id] = _sorted_segment_boundary(
                sorted_values[seg_start:seg_stop], np.dtype(values.dtype)
            )
        for segment_id in range(segment_count, nsegments_per_chunk):
            l2[row_start + segment_id] = l2[row_start + segment_count - 1]
        cursor = next_cursor

    return sorted_values, positions, offsets, l2, position_dtype


def _build_chunk_sorted_payload_direct(
    array: blosc2.NDArray,
    target: dict,
    dtype: np.dtype,
    chunk_len: int,
    nav_segment_len: int,
    *,
    payload_dtype: np.dtype | None = None,
    aux_dtype: np.dtype | None = None,
    value_transform=None,
    aux_transform=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    size = int(array.shape[0])
    nchunks = math.ceil(size / chunk_len)
    payload_dtype = np.dtype(dtype if payload_dtype is None else payload_dtype)
    aux_dtype = np.dtype(_position_dtype(chunk_len - 1) if aux_dtype is None else aux_dtype)
    offsets = np.empty(nchunks + 1, dtype=np.int64)
    offsets[0] = 0
    payload = np.empty(size, dtype=payload_dtype)
    aux = np.empty(size, dtype=aux_dtype)
    l1 = np.empty(nchunks, dtype=_boundary_dtype(payload_dtype))
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(payload_dtype))

    cursor = 0
    for chunk_id in range(nchunks):
        start = chunk_id * chunk_len
        stop = min(start + chunk_len, size)
        chunk = _slice_values_for_target(array, target, start, stop)
        order = np.argsort(chunk, kind="stable")
        chunk_size = stop - start
        next_cursor = cursor + chunk_size
        chunk_payload = chunk[order]
        if value_transform is not None:
            chunk_payload = value_transform(chunk_payload)
        chunk_aux = order.astype(_position_dtype(chunk_len - 1), copy=False)
        if aux_transform is not None:
            chunk_aux = aux_transform(chunk_aux)
        payload[cursor:next_cursor] = chunk_payload
        aux[cursor:next_cursor] = chunk_aux
        offsets[chunk_id + 1] = next_cursor
        if chunk_size > 0:
            l1[chunk_id] = _sorted_segment_boundary(chunk_payload, payload_dtype)
            row_start = chunk_id * nsegments_per_chunk
            segment_count = _segment_row_count(chunk_size, nav_segment_len)
            for segment_id in range(segment_count):
                seg_start = segment_id * nav_segment_len
                seg_stop = min(seg_start + nav_segment_len, chunk_size)
                l2[row_start + segment_id] = _sorted_segment_boundary(
                    chunk_payload[seg_start:seg_stop], payload_dtype
                )
            for segment_id in range(segment_count, nsegments_per_chunk):
                l2[row_start + segment_id] = l2[row_start + segment_count - 1]
        cursor = next_cursor

    return payload, aux, offsets, l1, l2


def _intra_chunk_run_ranges(chunk_size: int, thread_count: int) -> list[tuple[int, int]]:
    if chunk_size <= 0:
        return []
    run_count = max(1, min(thread_count, chunk_size))
    boundaries = np.linspace(0, chunk_size, run_count + 1, dtype=np.int64)
    return [(int(boundaries[idx]), int(boundaries[idx + 1])) for idx in range(run_count)]


def _sort_chunk_run(
    chunk: np.ndarray, run_start: int, run_stop: int, position_dtype: np.dtype
) -> tuple[np.ndarray, np.ndarray]:
    run = chunk[run_start:run_stop]
    try:
        return indexing_ext.intra_chunk_sort_run(run, run_start, position_dtype)
    except TypeError:
        order = np.argsort(run, kind="stable")
        return run[order], (order + run_start).astype(position_dtype, copy=False)


def _merge_sorted_run_pair(
    left_values: np.ndarray,
    left_positions: np.ndarray,
    right_values: np.ndarray,
    right_positions: np.ndarray,
    dtype: np.dtype,
    position_dtype: np.dtype,
) -> tuple[np.ndarray, np.ndarray]:
    try:
        merged_values, merged_positions = indexing_ext.intra_chunk_merge_sorted_slices(
            left_values, left_positions, right_values, right_positions, position_dtype
        )
    except TypeError:
        merged_values, merged_positions = _merge_sorted_slices(
            left_values, left_positions, right_values, right_positions, dtype
        )
    return merged_values, merged_positions.astype(position_dtype, copy=False)


def _sort_chunk_intra_chunk(
    chunk: np.ndarray,
    position_dtype: np.dtype,
    *,
    thread_count: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    chunk_size = chunk.shape[0]
    if chunk_size == 0:
        return np.empty(0, dtype=chunk.dtype), np.empty(0, dtype=position_dtype)
    if thread_count is None:
        thread_count = _index_build_threads()
    thread_count = max(1, min(int(thread_count), chunk_size))
    if thread_count <= 1:
        order = np.argsort(chunk, kind="stable")
        return chunk[order], order.astype(position_dtype, copy=False)

    def sort_run(run_range: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
        return _sort_chunk_run(chunk, run_range[0], run_range[1], position_dtype)

    run_ranges = _intra_chunk_run_ranges(chunk_size, thread_count)
    with ThreadPoolExecutor(max_workers=thread_count) as executor:
        runs = list(executor.map(sort_run, run_ranges))

    while len(runs) > 1:
        pair_specs = [(runs[idx], runs[idx + 1]) for idx in range(0, len(runs) - 1, 2)]

        def merge_pair(
            pair_spec: tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]],
        ) -> tuple[np.ndarray, np.ndarray]:
            (left_values, left_positions), (right_values, right_positions) = pair_spec
            return _merge_sorted_run_pair(
                left_values, left_positions, right_values, right_positions, chunk.dtype, position_dtype
            )

        if pair_specs:
            merge_workers = min(thread_count, len(pair_specs))
            if merge_workers <= 1:
                merged_runs = [merge_pair(pair_spec) for pair_spec in pair_specs]
            else:
                with ThreadPoolExecutor(max_workers=merge_workers) as executor:
                    merged_runs = list(executor.map(merge_pair, pair_specs))
        else:
            merged_runs = []
        if len(runs) % 2 == 1:
            merged_runs.append(runs[-1])
        runs = merged_runs

    return runs[0]


def _build_partial_chunk_payloads_intra_chunk(
    array: blosc2.NDArray,
    target: dict,
    dtype: np.dtype,
    chunk_len: int,
    nav_segment_len: int,
    cparams: dict | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    size = int(array.shape[0])
    nchunks = math.ceil(size / chunk_len)
    position_dtype = _position_dtype(chunk_len - 1)
    sorted_values = np.empty(size, dtype=dtype)
    positions = np.empty(size, dtype=position_dtype)
    offsets = np.empty(nchunks + 1, dtype=np.int64)
    offsets[0] = 0
    l1 = np.empty(nchunks, dtype=_boundary_dtype(dtype))
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(dtype))
    cursor = 0
    thread_count = _index_build_threads(cparams)

    for chunk_id in range(nchunks):
        start = chunk_id * chunk_len
        stop = min(start + chunk_len, size)
        chunk_sorted, local_positions = _sort_chunk_intra_chunk(
            _slice_values_for_target(array, target, start, stop), position_dtype, thread_count=thread_count
        )
        chunk_size = stop - start
        next_cursor = cursor + chunk_size
        sorted_values[cursor:next_cursor] = chunk_sorted
        positions[cursor:next_cursor] = local_positions
        offsets[chunk_id + 1] = next_cursor
        if chunk_size > 0:
            l1[chunk_id] = _sorted_segment_boundary(chunk_sorted, dtype)
            row_start = chunk_id * nsegments_per_chunk
            segment_count = _segment_row_count(chunk_size, nav_segment_len)
            for segment_id in range(segment_count):
                seg_start = segment_id * nav_segment_len
                seg_stop = min(seg_start + nav_segment_len, chunk_size)
                l2[row_start + segment_id] = _sorted_segment_boundary(
                    chunk_sorted[seg_start:seg_stop], dtype
                )
            for segment_id in range(segment_count, nsegments_per_chunk):
                l2[row_start + segment_id] = l2[row_start + segment_count - 1]
        cursor = next_cursor

    return sorted_values, positions, offsets, l1, l2


def _chunk_index_payload_storage(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    category: str,
    payload_name: str,
    payload: np.ndarray,
    aux_name: str,
    aux_payload: np.ndarray,
    offsets: np.ndarray,
    l1: np.ndarray,
    l2: np.ndarray,
    persistent: bool,
    chunk_len: int,
    nav_segment_len: int,
    cparams: dict | None = None,
) -> dict:
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    payload_sidecar = _store_array_sidecar(
        array,
        token,
        kind,
        category,
        payload_name,
        payload,
        persistent,
        chunks=(chunk_len,),
        blocks=(nav_segment_len,),
        cparams=cparams,
    )
    aux_sidecar = _store_array_sidecar(
        array,
        token,
        kind,
        category,
        aux_name,
        aux_payload,
        persistent,
        chunks=(chunk_len,),
        blocks=(nav_segment_len,),
        cparams=cparams,
    )
    offsets_sidecar = _store_array_sidecar(
        array, token, kind, category, "offsets", offsets, persistent, cparams=cparams
    )
    l1_sidecar = _store_array_sidecar(
        array, token, kind, f"{category}_nav", "l1", l1, persistent, cparams=cparams
    )
    l2_sidecar = _store_array_sidecar(
        array,
        token,
        kind,
        f"{category}_nav",
        "l2",
        l2,
        persistent,
        chunks=(nsegments_per_chunk,),
        blocks=(min(nsegments_per_chunk, max(1, nsegments_per_chunk)),),
        cparams=cparams,
    )
    return {
        "layout": "chunk-local-v1",
        "chunk_len": chunk_len,
        "nav_segment_len": nav_segment_len,
        "nsegments_per_chunk": nsegments_per_chunk,
        "values_path": payload_sidecar["path"],
        f"{aux_name}_path": aux_sidecar["path"],
        "offsets_path": offsets_sidecar["path"],
        "l1_path": l1_sidecar["path"],
        "l2_path": l2_sidecar["path"],
    }


def _prepare_chunk_index_payload_sidecars(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    category: str,
    payload_name: str,
    payload_dtype: np.dtype,
    aux_name: str,
    aux_dtype: np.dtype,
    size: int,
    chunk_len: int,
    nav_segment_len: int,
    cparams: dict | None = None,
) -> tuple[blosc2.NDArray | None, dict, blosc2.NDArray | None, dict]:
    payload_handle, payload_sidecar = _create_persistent_sidecar_handle(
        array,
        token,
        kind,
        category,
        payload_name,
        size,
        payload_dtype,
        chunks=(chunk_len,),
        blocks=(nav_segment_len,),
        cparams=cparams,
    )
    aux_handle, aux_sidecar = _create_persistent_sidecar_handle(
        array,
        token,
        kind,
        category,
        aux_name,
        size,
        aux_dtype,
        chunks=(chunk_len,),
        blocks=(nav_segment_len,),
        cparams=cparams,
    )
    return payload_handle, payload_sidecar, aux_handle, aux_sidecar


def _finalize_chunk_index_payload_storage(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    category: str,
    aux_name: str,
    offsets: np.ndarray,
    l1: np.ndarray,
    l2: np.ndarray,
    payload_sidecar: dict,
    aux_sidecar: dict,
    chunk_len: int,
    nav_segment_len: int,
    cparams: dict | None = None,
) -> dict:
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    offsets_sidecar = _store_array_sidecar(
        array, token, kind, category, "offsets", offsets, True, cparams=cparams
    )
    l1_sidecar = _store_array_sidecar(array, token, kind, f"{category}_nav", "l1", l1, True, cparams=cparams)
    l2_sidecar = _store_array_sidecar(
        array,
        token,
        kind,
        f"{category}_nav",
        "l2",
        l2,
        True,
        chunks=(nsegments_per_chunk,),
        blocks=(min(nsegments_per_chunk, max(1, nsegments_per_chunk)),),
        cparams=cparams,
    )
    return {
        "layout": "chunk-local-v1",
        "chunk_len": chunk_len,
        "nav_segment_len": nav_segment_len,
        "nsegments_per_chunk": nsegments_per_chunk,
        "values_path": payload_sidecar["path"],
        f"{aux_name}_path": aux_sidecar["path"],
        "offsets_path": offsets_sidecar["path"],
        "l1_path": l1_sidecar["path"],
        "l2_path": l2_sidecar["path"],
    }


def _build_partial_descriptor_ooc(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    optlevel: int,
    persistent: bool,
    cparams: dict | None = None,
) -> dict:
    if persistent:
        size = int(array.shape[0])
        chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
        chunk_len = int(array.chunks[0]) * chunk_multiplier
        nav_segment_len, nav_segment_divisor = _partial_nav_segment_len(
            int(array.blocks[0]), chunk_len, optlevel
        )
        nchunks = math.ceil(size / chunk_len)
        position_dtype = _position_dtype(chunk_len - 1)
        offsets = np.empty(nchunks + 1, dtype=np.int64)
        offsets[0] = 0
        l1 = np.empty(nchunks, dtype=_boundary_dtype(dtype))
        nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
        l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(dtype))
        values_handle = positions_handle = None
        values_sidecar = positions_sidecar = None
        try:
            values_handle, values_sidecar, positions_handle, positions_sidecar = (
                _prepare_chunk_index_payload_sidecars(
                    array,
                    token,
                    kind,
                    "partial",
                    "values",
                    dtype,
                    "positions",
                    position_dtype,
                    size,
                    chunk_len,
                    nav_segment_len,
                    cparams,
                )
            )
            cursor = 0
            for chunk_id in range(nchunks):
                start = chunk_id * chunk_len
                stop = min(start + chunk_len, size)
                chunk_size = stop - start
                next_cursor = cursor + chunk_size
                chunk_sorted, local_positions = _sort_chunk_intra_chunk(
                    _slice_values_for_target(array, target, start, stop), position_dtype
                )
                if values_handle is not None:
                    values_handle[cursor:next_cursor] = chunk_sorted
                if positions_handle is not None:
                    positions_handle[cursor:next_cursor] = local_positions
                offsets[chunk_id + 1] = next_cursor
                if chunk_size > 0:
                    l1[chunk_id] = _sorted_segment_boundary(chunk_sorted, dtype)
                    row_start = chunk_id * nsegments_per_chunk
                    segment_count = _segment_row_count(chunk_size, nav_segment_len)
                    for segment_id in range(segment_count):
                        seg_start = segment_id * nav_segment_len
                        seg_stop = min(seg_start + nav_segment_len, chunk_size)
                        l2[row_start + segment_id] = _sorted_segment_boundary(
                            chunk_sorted[seg_start:seg_stop], dtype
                        )
                    for segment_id in range(segment_count, nsegments_per_chunk):
                        l2[row_start + segment_id] = l2[row_start + segment_count - 1]
                cursor = next_cursor
            del values_handle, positions_handle
            partial = _finalize_chunk_index_payload_storage(
                array,
                token,
                kind,
                "partial",
                "positions",
                offsets,
                l1,
                l2,
                values_sidecar,
                positions_sidecar,
                chunk_len,
                nav_segment_len,
                cparams,
            )
        except Exception:
            if values_sidecar is not None:
                _remove_sidecar_path(values_sidecar["path"])
            if positions_sidecar is not None:
                _remove_sidecar_path(positions_sidecar["path"])
            raise
        partial["position_dtype"] = position_dtype.str
        partial["nav_segment_divisor"] = nav_segment_divisor
        partial["chunk_multiplier"] = chunk_multiplier
        return partial

    chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
    chunk_len = int(array.chunks[0]) * chunk_multiplier
    nav_segment_len, nav_segment_divisor = _partial_nav_segment_len(
        int(array.blocks[0]), chunk_len, optlevel
    )
    sorted_values, positions, offsets, l1, l2 = _build_partial_chunk_payloads_intra_chunk(
        array, target, dtype, chunk_len, nav_segment_len, cparams
    )
    partial = _chunk_index_payload_storage(
        array,
        token,
        kind,
        "partial",
        "values",
        sorted_values,
        "positions",
        positions,
        offsets,
        l1,
        l2,
        persistent,
        chunk_len,
        nav_segment_len,
        cparams,
    )
    partial["position_dtype"] = positions.dtype.str
    partial["nav_segment_divisor"] = nav_segment_divisor
    partial["chunk_multiplier"] = chunk_multiplier
    return partial


def _bucket_bucket_count(block_len: int) -> int:
    return max(1, min(64, block_len))


def _pack_bucket_mask(bucket_ids: np.ndarray) -> np.uint64:
    mask = np.uint64(0)
    for bucket_id in np.unique(bucket_ids):
        mask |= np.uint64(1) << np.uint64(int(bucket_id))
    return mask


def _bucket_value_lossy_bits(dtype: np.dtype, optlevel: int) -> int:
    dtype = np.dtype(dtype)
    if dtype.kind in {"i", "u"} or dtype == np.dtype(np.float32) or dtype == np.dtype(np.float64):
        max_bits = dtype.itemsize
    else:
        return 0
    return min(max(0, 9 - int(optlevel)), max_bits)


def _quantize_integer_array(values: np.ndarray, bits: int) -> np.ndarray:
    if bits <= 0:
        return values
    dtype = np.dtype(values.dtype)
    base_mask = np.iinfo(dtype).max if dtype.kind == "u" else -1
    mask = np.asarray(base_mask ^ ((1 << bits) - 1), dtype=dtype)[()]
    quantized = values.copy()
    np.bitwise_and(quantized, mask, out=quantized)
    return quantized


def _quantize_integer_scalar(value, dtype: np.dtype, bits: int):
    scalar = np.asarray(value, dtype=dtype)[()]
    if bits <= 0:
        return scalar
    base_mask = np.iinfo(dtype).max if dtype.kind == "u" else -1
    mask = np.asarray(base_mask ^ ((1 << bits) - 1), dtype=dtype)[()]
    return np.bitwise_and(scalar, mask, dtype=dtype)


def _float_order_uint_dtype(dtype: np.dtype) -> np.dtype:
    if dtype == np.dtype(np.float32):
        return np.dtype(np.uint32)
    if dtype == np.dtype(np.float64):
        return np.dtype(np.uint64)
    raise TypeError(f"unsupported float dtype {dtype}")


def _ordered_uint_from_float(values: np.ndarray) -> np.ndarray:
    dtype = np.dtype(values.dtype)
    uint_dtype = _float_order_uint_dtype(dtype)
    bits = values.view(uint_dtype).copy()
    sign_mask = np.asarray(1 << (dtype.itemsize * 8 - 1), dtype=uint_dtype)[()]
    negative = (bits & sign_mask) != 0
    bits[negative] = ~bits[negative]
    bits[~negative] ^= sign_mask
    return bits


def _float_from_ordered_uint(ordered: np.ndarray, dtype: np.dtype) -> np.ndarray:
    uint_dtype = _float_order_uint_dtype(dtype)
    bits = ordered.astype(uint_dtype, copy=True)
    sign_mask = np.asarray(1 << (dtype.itemsize * 8 - 1), dtype=uint_dtype)[()]
    positive = (bits & sign_mask) != 0
    bits[positive] ^= sign_mask
    bits[~positive] = ~bits[~positive]
    return bits.view(dtype)


def _quantize_float_array(values: np.ndarray, bits: int) -> np.ndarray:
    if bits <= 0:
        return values
    quantized = values.copy()
    finite = np.isfinite(quantized)
    if not np.any(finite):
        return quantized
    ordered = _ordered_uint_from_float(quantized[finite])
    uint_dtype = ordered.dtype
    mask = np.asarray(np.iinfo(uint_dtype).max ^ ((1 << bits) - 1), dtype=uint_dtype)[()]
    np.bitwise_and(ordered, mask, out=ordered)
    quantized[finite] = _float_from_ordered_uint(ordered, quantized.dtype)
    return quantized


def _quantize_float_scalar(value, dtype: np.dtype, bits: int):
    scalar = np.asarray(value, dtype=dtype)[()]
    if bits <= 0 or not np.isfinite(scalar):
        return scalar
    ordered = _ordered_uint_from_float(np.asarray([scalar], dtype=dtype))
    uint_dtype = ordered.dtype
    mask = np.asarray(np.iinfo(uint_dtype).max ^ ((1 << bits) - 1), dtype=uint_dtype)[()]
    np.bitwise_and(ordered, mask, out=ordered)
    return _float_from_ordered_uint(ordered, dtype)[0]


def _quantize_bucket_values_array(values: np.ndarray, bits: int) -> np.ndarray:
    dtype = np.dtype(values.dtype)
    if bits <= 0:
        return values
    if dtype.kind in {"i", "u"}:
        return _quantize_integer_array(values, bits)
    if dtype == np.dtype(np.float32) or dtype == np.dtype(np.float64):
        return _quantize_float_array(values, bits)
    return values


def _quantize_bucket_value_scalar(value, dtype: np.dtype, bits: int):
    dtype = np.dtype(dtype)
    if bits <= 0:
        return np.asarray(value, dtype=dtype)[()]
    if dtype.kind in {"i", "u"}:
        return _quantize_integer_scalar(value, dtype, bits)
    if dtype == np.dtype(np.float32) or dtype == np.dtype(np.float64):
        return _quantize_float_scalar(value, dtype, bits)
    return np.asarray(value, dtype=dtype)[()]


def _build_bucket_chunk_payloads(
    array: blosc2.NDArray,
    target: dict,
    dtype: np.dtype,
    chunk_len: int,
    nav_segment_len: int,
    value_lossy_bits: int,
    bucket_len: int,
    bucket_dtype: np.dtype,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    size = int(array.shape[0])
    nchunks = math.ceil(size / chunk_len)
    offsets = np.empty(nchunks + 1, dtype=np.int64)
    offsets[0] = 0
    sorted_values = np.empty(size, dtype=dtype)
    bucket_positions = np.empty(size, dtype=bucket_dtype)
    l1 = np.empty(nchunks, dtype=_boundary_dtype(dtype))
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(dtype))
    position_dtype = _position_dtype(chunk_len - 1)
    cursor = 0

    for chunk_id in range(nchunks):
        start = chunk_id * chunk_len
        stop = min(start + chunk_len, size)
        chunk = _slice_values_for_target(array, target, start, stop)
        order = np.argsort(chunk, kind="stable")
        chunk_size = stop - start
        next_cursor = cursor + chunk_size
        chunk_sorted = chunk[order]
        stored_chunk_sorted = chunk_sorted
        if value_lossy_bits > 0:
            stored_chunk_sorted = _quantize_bucket_values_array(chunk_sorted, value_lossy_bits)
        local_positions = order.astype(position_dtype, copy=False)
        sorted_values[cursor:next_cursor] = stored_chunk_sorted
        bucket_positions[cursor:next_cursor] = (local_positions // bucket_len).astype(
            bucket_dtype, copy=False
        )
        offsets[chunk_id + 1] = next_cursor
        if chunk_size > 0:
            l1[chunk_id] = _sorted_segment_boundary(stored_chunk_sorted, dtype)
            row_start = chunk_id * nsegments_per_chunk
            segment_count = _segment_row_count(chunk_size, nav_segment_len)
            for segment_id in range(segment_count):
                seg_start = segment_id * nav_segment_len
                seg_stop = min(seg_start + nav_segment_len, chunk_size)
                l2[row_start + segment_id] = _sorted_segment_boundary(
                    chunk_sorted[seg_start:seg_stop], dtype
                )
            for segment_id in range(segment_count, nsegments_per_chunk):
                l2[row_start + segment_id] = l2[row_start + segment_count - 1]
        cursor = next_cursor

    return sorted_values, bucket_positions, offsets, l1, l2


def _build_bucket_chunk_payloads_intra_chunk(
    array: blosc2.NDArray,
    target: dict,
    dtype: np.dtype,
    chunk_len: int,
    nav_segment_len: int,
    value_lossy_bits: int,
    bucket_len: int,
    bucket_dtype: np.dtype,
    cparams: dict | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    size = int(array.shape[0])
    nchunks = math.ceil(size / chunk_len)
    offsets = np.empty(nchunks + 1, dtype=np.int64)
    offsets[0] = 0
    sorted_values = np.empty(size, dtype=dtype)
    bucket_positions = np.empty(size, dtype=bucket_dtype)
    l1 = np.empty(nchunks, dtype=_boundary_dtype(dtype))
    nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
    l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(dtype))
    position_dtype = _position_dtype(chunk_len - 1)
    cursor = 0
    thread_count = _index_build_threads(cparams)

    for chunk_id in range(nchunks):
        start = chunk_id * chunk_len
        stop = min(start + chunk_len, size)
        chunk_sorted, local_positions = _sort_chunk_intra_chunk(
            _slice_values_for_target(array, target, start, stop), position_dtype, thread_count=thread_count
        )
        chunk_size = stop - start
        next_cursor = cursor + chunk_size
        stored_chunk_sorted = chunk_sorted
        if value_lossy_bits > 0:
            stored_chunk_sorted = _quantize_bucket_values_array(chunk_sorted, value_lossy_bits)
        sorted_values[cursor:next_cursor] = stored_chunk_sorted
        bucket_positions[cursor:next_cursor] = (local_positions // bucket_len).astype(
            bucket_dtype, copy=False
        )
        offsets[chunk_id + 1] = next_cursor
        if chunk_size > 0:
            l1[chunk_id] = _sorted_segment_boundary(stored_chunk_sorted, dtype)
            row_start = chunk_id * nsegments_per_chunk
            segment_count = _segment_row_count(chunk_size, nav_segment_len)
            for segment_id in range(segment_count):
                seg_start = segment_id * nav_segment_len
                seg_stop = min(seg_start + nav_segment_len, chunk_size)
                l2[row_start + segment_id] = _sorted_segment_boundary(
                    chunk_sorted[seg_start:seg_stop], dtype
                )
            for segment_id in range(segment_count, nsegments_per_chunk):
                l2[row_start + segment_id] = l2[row_start + segment_count - 1]
        cursor = next_cursor

    return sorted_values, bucket_positions, offsets, l1, l2


def _build_bucket_descriptor(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    values: np.ndarray,
    optlevel: int,
    persistent: bool,
    cparams: dict | None = None,
) -> dict:
    chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
    chunk_len = int(array.chunks[0]) * chunk_multiplier
    nav_segment_len = int(array.blocks[0])
    bucket_len = max(1, math.ceil(nav_segment_len / 64))
    bucket_count = math.ceil(chunk_len / bucket_len)
    value_lossy_bits = _bucket_value_lossy_bits(values.dtype, optlevel)
    sorted_values, positions, offsets, l2, _ = _build_chunk_sorted_payload(
        values, chunk_len, nav_segment_len, cparams
    )
    if value_lossy_bits > 0:
        sorted_values = _quantize_bucket_values_array(sorted_values, value_lossy_bits)
    bucket_dtype = _position_dtype(bucket_count - 1)
    bucket_positions = (positions // bucket_len).astype(bucket_dtype, copy=False)
    l1 = _compute_sorted_boundaries(sorted_values, np.dtype(sorted_values.dtype), chunk_len)
    bucket = _chunk_index_payload_storage(
        array,
        token,
        kind,
        "bucket",
        "values",
        sorted_values,
        "bucket_positions",
        bucket_positions,
        offsets,
        l1,
        l2,
        persistent,
        chunk_len,
        nav_segment_len,
        cparams,
    )
    bucket["bucket_count"] = bucket_count
    bucket["bucket_len"] = bucket_len
    bucket["chunk_multiplier"] = chunk_multiplier
    bucket["value_lossy_bits"] = value_lossy_bits
    bucket["bucket_dtype"] = bucket_positions.dtype.str
    return bucket


def _build_bucket_descriptor_ooc(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    optlevel: int,
    persistent: bool,
    cparams: dict | None = None,
) -> dict:
    chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
    chunk_len = int(array.chunks[0]) * chunk_multiplier
    nav_segment_len = int(array.blocks[0])
    bucket_len = max(1, math.ceil(nav_segment_len / 64))
    bucket_count = math.ceil(chunk_len / bucket_len)
    value_lossy_bits = _bucket_value_lossy_bits(dtype, optlevel)
    bucket_dtype = _position_dtype(bucket_count - 1)

    if persistent:
        # Streaming path: write directly into sidecars chunk by chunk so that
        # we never allocate full-row-count payload arrays in RAM.
        size = int(array.shape[0])
        nchunks = math.ceil(size / chunk_len)
        nsegments_per_chunk = _segment_row_count(chunk_len, nav_segment_len)
        position_dtype = _position_dtype(chunk_len - 1)
        thread_count = _index_build_threads(cparams)

        offsets = np.empty(nchunks + 1, dtype=np.int64)
        offsets[0] = 0
        l1 = np.empty(nchunks, dtype=_boundary_dtype(dtype))
        l2 = np.empty(nchunks * nsegments_per_chunk, dtype=_boundary_dtype(dtype))

        values_handle = bucket_handle = None
        values_sidecar = bucket_sidecar = None
        try:
            values_handle, values_sidecar, bucket_handle, bucket_sidecar = (
                _prepare_chunk_index_payload_sidecars(
                    array,
                    token,
                    kind,
                    "bucket",
                    "values",
                    dtype,
                    "bucket_positions",
                    bucket_dtype,
                    size,
                    chunk_len,
                    nav_segment_len,
                    cparams,
                )
            )
            cursor = 0
            for chunk_id in range(nchunks):
                start = chunk_id * chunk_len
                stop = min(start + chunk_len, size)
                chunk_size = stop - start
                next_cursor = cursor + chunk_size
                chunk_sorted, local_positions = _sort_chunk_intra_chunk(
                    _slice_values_for_target(array, target, start, stop),
                    position_dtype,
                    thread_count=thread_count,
                )
                stored_chunk_sorted = chunk_sorted
                if value_lossy_bits > 0:
                    stored_chunk_sorted = _quantize_bucket_values_array(chunk_sorted, value_lossy_bits)
                if values_handle is not None:
                    values_handle[cursor:next_cursor] = stored_chunk_sorted
                if bucket_handle is not None:
                    bucket_handle[cursor:next_cursor] = (local_positions // bucket_len).astype(
                        bucket_dtype, copy=False
                    )
                offsets[chunk_id + 1] = next_cursor
                if chunk_size > 0:
                    l1[chunk_id] = _sorted_segment_boundary(stored_chunk_sorted, dtype)
                    row_start = chunk_id * nsegments_per_chunk
                    segment_count = _segment_row_count(chunk_size, nav_segment_len)
                    for segment_id in range(segment_count):
                        seg_start = segment_id * nav_segment_len
                        seg_stop = min(seg_start + nav_segment_len, chunk_size)
                        l2[row_start + segment_id] = _sorted_segment_boundary(
                            chunk_sorted[seg_start:seg_stop], dtype
                        )
                    for segment_id in range(segment_count, nsegments_per_chunk):
                        l2[row_start + segment_id] = l2[row_start + segment_count - 1]
                cursor = next_cursor
            del values_handle, bucket_handle
            bucket = _finalize_chunk_index_payload_storage(
                array,
                token,
                kind,
                "bucket",
                "bucket_positions",
                offsets,
                l1,
                l2,
                values_sidecar,
                bucket_sidecar,
                chunk_len,
                nav_segment_len,
                cparams,
            )
        except Exception:
            if values_sidecar is not None:
                _remove_sidecar_path(values_sidecar["path"])
            if bucket_sidecar is not None:
                _remove_sidecar_path(bucket_sidecar["path"])
            raise
        bucket["bucket_count"] = bucket_count
        bucket["bucket_len"] = bucket_len
        bucket["chunk_multiplier"] = chunk_multiplier
        bucket["value_lossy_bits"] = value_lossy_bits
        bucket["bucket_dtype"] = bucket_dtype.str
        return bucket

    # Non-persistent path: full staging in RAM is acceptable.
    sorted_values, bucket_positions, offsets, l1, l2 = _build_bucket_chunk_payloads_intra_chunk(
        array,
        target,
        dtype,
        chunk_len,
        nav_segment_len,
        value_lossy_bits,
        bucket_len,
        bucket_dtype,
        cparams,
    )
    bucket = _chunk_index_payload_storage(
        array,
        token,
        kind,
        "bucket",
        "values",
        sorted_values,
        "bucket_positions",
        bucket_positions,
        offsets,
        l1,
        l2,
        persistent,
        chunk_len,
        nav_segment_len,
        cparams,
    )
    bucket["bucket_count"] = bucket_count
    bucket["bucket_len"] = bucket_len
    bucket["chunk_multiplier"] = chunk_multiplier
    bucket["value_lossy_bits"] = value_lossy_bits
    bucket["bucket_dtype"] = bucket_positions.dtype.str
    return bucket


def _scalar_compare(left, right, dtype: np.dtype) -> int:
    dtype = np.dtype(dtype)
    if dtype.kind == "f":
        left_nan = np.isnan(left)
        right_nan = np.isnan(right)
        if left_nan and right_nan:
            return 0
        if left_nan:
            return 1
        if right_nan:
            return -1
    if left < right:
        return -1
    if left > right:
        return 1
    return 0


def _pair_le(left_value, left_position: int, right_value, right_position: int, dtype: np.dtype) -> bool:
    cmp = _scalar_compare(left_value, right_value, dtype)
    if cmp < 0:
        return True
    if cmp > 0:
        return False
    return int(left_position) <= int(right_position)


def _pair_record_dtype(dtype: np.dtype) -> np.dtype:
    return np.dtype([("value", dtype), ("position", np.int64)])


def _pair_records(values: np.ndarray, positions: np.ndarray, dtype: np.dtype) -> np.ndarray:
    records = np.empty(values.shape[0], dtype=_pair_record_dtype(dtype))
    records["value"] = values
    records["position"] = positions
    return records


def _merge_sorted_slices(
    left_values: np.ndarray,
    left_positions: np.ndarray,
    right_values: np.ndarray,
    right_positions: np.ndarray,
    dtype: np.dtype,
) -> tuple[np.ndarray, np.ndarray]:
    if left_values.size == 0:
        return right_values, right_positions
    if right_values.size == 0:
        return left_values, left_positions
    values = np.concatenate((left_values, right_values))
    positions = np.concatenate((left_positions, right_positions))
    order = np.lexsort((positions, values))
    return values[order], positions[order]


def _pair_searchsorted_right(values: np.ndarray, positions: np.ndarray, value, position: int) -> int:
    records = _pair_records(values, positions, values.dtype)
    needle = np.asarray((value, position), dtype=records.dtype)[()]
    return int(np.searchsorted(records, needle, side="right"))


def _temp_run_storage_geometry(
    length: int, dtype: np.dtype, buffer_items: int
) -> tuple[tuple[int], tuple[int]]:
    chunk_items = max(1, min(length, buffer_items))
    target_block_bytes = 256 * 1024
    block_items = max(1, min(chunk_items, target_block_bytes // max(1, dtype.itemsize)))
    return (chunk_items,), (block_items,)


def _path_disk_bytes(path: Path | str) -> int:
    path = Path(path)
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size
    return sum(entry.stat().st_size for entry in path.rglob("*") if entry.is_file())


def _tracker_register_create(tracker: TempRunTracker | None, *paths: Path) -> None:
    if tracker is None:
        return
    delta = sum(_path_disk_bytes(path) for path in paths)
    tracker.current_disk_bytes += delta
    tracker.total_written_bytes += delta
    tracker.peak_disk_bytes = max(tracker.peak_disk_bytes, tracker.current_disk_bytes)


def _tracker_register_delete(tracker: TempRunTracker | None, *paths: Path) -> None:
    if tracker is None:
        return
    delta = sum(_path_disk_bytes(path) for path in paths)
    tracker.current_disk_bytes = max(0, tracker.current_disk_bytes - delta)


def _create_blosc2_temp_array(
    path: Path, length: int, dtype: np.dtype, buffer_items: int, cparams: dict | None = None
):
    chunks, blocks = _temp_run_storage_geometry(length, dtype, buffer_items)
    if cparams is None:
        cparams = blosc2.CParams(codec=blosc2.Codec.ZSTD, clevel=1)
    return blosc2.empty(
        (length,),
        dtype=dtype,
        chunks=chunks,
        blocks=blocks,
        urlpath=str(path),
        mode="w",
        cparams=cparams,
    )


def _read_ndarray_linear_span(array: blosc2.NDArray | np.ndarray, start: int, out: np.ndarray) -> None:
    if len(out) == 0:
        return
    if isinstance(array, np.ndarray):
        out[...] = array[start : start + len(out)]
        return
    chunk_len = int(array.chunks[0])
    cursor = int(start)
    out_cursor = 0
    while out_cursor < len(out):
        chunk_id = cursor // chunk_len
        local_start = cursor % chunk_len
        take = min(len(out) - out_cursor, chunk_len - local_start)
        target = out[out_cursor : out_cursor + take]
        try:
            array.get_1d_span_numpy(target, int(chunk_id), int(local_start), int(take))
        except RuntimeError:
            # Newer c-blosc2 builds can reject getitem on some tiny persistent
            # sidecars even though normal ndarray slicing still succeeds.
            target[...] = array[cursor : cursor + take]
        cursor += take
        out_cursor += take


def _write_ndarray_linear_span(array: blosc2.NDArray | np.ndarray, start: int, values: np.ndarray) -> None:
    if len(values) == 0:
        return
    stop = int(start) + len(values)
    if stop > int(array.shape[0]):
        raise RuntimeError(
            f"attempted to write past the end of temporary array: stop={stop}, length={int(array.shape[0])}"
        )
    if isinstance(array, np.ndarray):
        array[start:stop] = values
        return
    chunk_len = int(array.chunks[0])
    cursor = int(start)
    in_cursor = 0
    while in_cursor < len(values):
        chunk_id = cursor // chunk_len
        local_start = cursor % chunk_len
        take = min(len(values) - in_cursor, chunk_len - local_start)
        try:
            array[cursor : cursor + take] = values[in_cursor : in_cursor + take]
        except Exception as exc:
            raise RuntimeError(
                "failed temporary sidecar span write: "
                f"array_len={int(array.shape[0])}, chunk_len={chunk_len}, "
                f"write_start={cursor}, write_stop={cursor + take}, "
                f"write_items={take}, local_start={local_start}, chunk_id={chunk_id}, "
                f"input_offset={in_cursor}, input_len={len(values)}, dtype={np.dtype(array.dtype)}"
            ) from exc
        cursor += take
        in_cursor += take


def _read_sidecar_span(handle, start: int, stop: int) -> np.ndarray:
    if stop <= start:
        return np.empty(0, dtype=np.dtype(handle.dtype))
    out = np.empty(stop - start, dtype=np.dtype(handle.dtype))
    _read_ndarray_linear_span(handle, start, out)
    return out


def _materialize_sorted_run(
    values: np.ndarray,
    positions: np.ndarray,
    length: int,
    value_dtype: np.dtype,
    workdir: Path,
    prefix: str,
    tracker: TempRunTracker | None = None,
    cparams: dict | None = None,
) -> SortedRun:
    values_path = workdir / f"{prefix}.values.b2nd"
    positions_path = workdir / f"{prefix}.positions.b2nd"
    run_values = _create_blosc2_temp_array(
        values_path, length, value_dtype, FULL_OOC_MERGE_BUFFER_ITEMS, cparams
    )
    run_positions = _create_blosc2_temp_array(
        positions_path, length, np.dtype(np.int64), FULL_OOC_MERGE_BUFFER_ITEMS, cparams
    )
    _write_ndarray_linear_span(run_values, 0, values)
    _write_ndarray_linear_span(run_positions, 0, positions)
    del run_values, run_positions
    _tracker_register_create(tracker, values_path, positions_path)
    return SortedRun(values_path, positions_path, length)


def _copy_sidecar_to_temp_run(
    path: str,
    length: int,
    dtype: np.dtype,
    workdir: Path,
    prefix: str,
    tracker: TempRunTracker | None = None,
    cparams: dict | None = None,
) -> Path:
    out_path = workdir / f"{prefix}.b2nd"
    sidecar = _open_sidecar_file(path, _INDEX_MMAP_MODE)
    output = _create_blosc2_temp_array(out_path, length, dtype, FULL_OOC_MERGE_BUFFER_ITEMS, cparams)
    chunk_len = int(sidecar.chunks[0])
    for chunk_id, start in enumerate(range(0, length, chunk_len)):
        stop = min(start + chunk_len, length)
        span = np.empty(stop - start, dtype=dtype)
        sidecar.get_1d_span_numpy(span, chunk_id, 0, stop - start)
        output[start:stop] = span
    del output
    _tracker_register_create(tracker, out_path)
    return out_path


def _copy_sidecar_handle_to_temp_run(
    handle,
    length: int,
    dtype: np.dtype,
    workdir: Path,
    prefix: str,
    tracker: TempRunTracker | None = None,
    cparams: dict | None = None,
) -> Path:
    out_path = workdir / f"{prefix}.b2nd"
    output = _create_blosc2_temp_array(out_path, length, dtype, FULL_OOC_MERGE_BUFFER_ITEMS, cparams)
    chunk_len = int(handle.chunks[0]) if hasattr(handle, "chunks") else length
    for start in range(0, length, chunk_len):
        stop = min(start + chunk_len, length)
        output[start:stop] = _read_sidecar_span(handle, start, stop)
    del output
    _tracker_register_create(tracker, out_path)
    return out_path


def _refill_run_buffer(
    values_src, positions_src, cursor: int, buffer_items: int
) -> tuple[np.ndarray, np.ndarray, int]:
    if cursor >= len(values_src):
        values_dtype = values_src.dtype if hasattr(values_src, "dtype") else np.float64
        positions_dtype = positions_src.dtype if hasattr(positions_src, "dtype") else np.int64
        return np.empty(0, dtype=values_dtype), np.empty(0, dtype=positions_dtype), cursor
    stop = min(cursor + buffer_items, len(values_src))
    if isinstance(values_src, np.ndarray):
        return np.asarray(values_src[cursor:stop]), np.asarray(positions_src[cursor:stop]), stop
    values = np.empty(stop - cursor, dtype=np.dtype(values_src.dtype))
    positions = np.empty(stop - cursor, dtype=np.dtype(positions_src.dtype))
    _read_ndarray_linear_span(values_src, cursor, values)
    _read_ndarray_linear_span(positions_src, cursor, positions)
    return values, positions, stop


def _merge_run_pair(
    left: SortedRun,
    right: SortedRun,
    workdir: Path,
    dtype: np.dtype,
    merge_id: int,
    buffer_items: int,
    tracker: TempRunTracker | None = None,
    cparams: dict | None = None,
) -> SortedRun:
    left_values_mm = blosc2.open(str(left.values_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)
    left_positions_mm = blosc2.open(str(left.positions_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)
    right_values_mm = blosc2.open(str(right.values_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)
    right_positions_mm = blosc2.open(str(right.positions_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)

    out_values_path = workdir / f"full_merge_values_{merge_id}.b2nd"
    out_positions_path = workdir / f"full_merge_positions_{merge_id}.b2nd"
    out_values = _create_blosc2_temp_array(
        out_values_path, left.length + right.length, dtype, buffer_items, cparams
    )
    out_positions = _create_blosc2_temp_array(
        out_positions_path, left.length + right.length, np.dtype(np.int64), buffer_items, cparams
    )
    out_total = left.length + right.length

    left_cursor = 0
    right_cursor = 0
    out_cursor = 0
    left_values = np.empty(0, dtype=dtype)
    left_positions = np.empty(0, dtype=np.int64)
    right_values = np.empty(0, dtype=dtype)
    right_positions = np.empty(0, dtype=np.int64)
    while True:
        if left_values.size == 0:
            left_values, left_positions, left_cursor = _refill_run_buffer(
                left_values_mm, left_positions_mm, left_cursor, buffer_items
            )
        if right_values.size == 0:
            right_values, right_positions, right_cursor = _refill_run_buffer(
                right_values_mm, right_positions_mm, right_cursor, buffer_items
            )

        if left_values.size == 0 and right_values.size == 0:
            break
        if left_values.size == 0:
            take = right_values.size
            try:
                _write_ndarray_linear_span(out_values, out_cursor, right_values)
                _write_ndarray_linear_span(out_positions, out_cursor, right_positions)
            except Exception as exc:
                raise RuntimeError(
                    "sorted index OOC merge write failed while flushing right run remainder: "
                    f"merge_id={merge_id}, left_len={left.length}, right_len={right.length}, "
                    f"out_total={out_total}, out_cursor={out_cursor}, take={take}, "
                    f"left_cursor={left_cursor}, right_cursor={right_cursor}, buffer_items={buffer_items}"
                ) from exc
            out_cursor += take
            right_values = np.empty(0, dtype=dtype)
            right_positions = np.empty(0, dtype=np.int64)
            continue
        if right_values.size == 0:
            take = left_values.size
            try:
                _write_ndarray_linear_span(out_values, out_cursor, left_values)
                _write_ndarray_linear_span(out_positions, out_cursor, left_positions)
            except Exception as exc:
                raise RuntimeError(
                    "sorted index OOC merge write failed while flushing left run remainder: "
                    f"merge_id={merge_id}, left_len={left.length}, right_len={right.length}, "
                    f"out_total={out_total}, out_cursor={out_cursor}, take={take}, "
                    f"left_cursor={left_cursor}, right_cursor={right_cursor}, buffer_items={buffer_items}"
                ) from exc
            out_cursor += take
            left_values = np.empty(0, dtype=dtype)
            left_positions = np.empty(0, dtype=np.int64)
            continue

        if _pair_le(left_values[-1], left_positions[-1], right_values[-1], right_positions[-1], dtype):
            left_cut = left_values.size
            right_cut = _pair_searchsorted_right(
                right_values, right_positions, left_values[-1], int(left_positions[-1])
            )
        else:
            left_cut = _pair_searchsorted_right(
                left_values, left_positions, right_values[-1], int(right_positions[-1])
            )
            right_cut = right_values.size

        try:
            merged_values, merged_positions = indexing_ext.intra_chunk_merge_sorted_slices(
                left_values[:left_cut],
                left_positions[:left_cut],
                right_values[:right_cut],
                right_positions[:right_cut],
                np.int64,
            )
        except (TypeError, AttributeError):
            # ponytail: fallback for non-numeric dtypes (strings, etc.) that the
            # Cython merge doesn't support.  _merge_sorted_slices uses np.lexsort.
            merged_values, merged_positions = _merge_sorted_slices(
                left_values[:left_cut],
                left_positions[:left_cut],
                right_values[:right_cut],
                right_positions[:right_cut],
                dtype,
            )
        take = merged_values.size
        try:
            _write_ndarray_linear_span(out_values, out_cursor, merged_values)
            _write_ndarray_linear_span(out_positions, out_cursor, merged_positions)
        except Exception as exc:
            raise RuntimeError(
                "sorted index OOC merge write failed for merged batch: "
                f"merge_id={merge_id}, left_len={left.length}, right_len={right.length}, "
                f"out_total={out_total}, out_cursor={out_cursor}, take={take}, "
                f"left_cursor={left_cursor}, right_cursor={right_cursor}, "
                f"left_buffer={left_values.size}, right_buffer={right_values.size}, "
                f"left_cut={left_cut}, right_cut={right_cut}, buffer_items={buffer_items}"
            ) from exc
        out_cursor += take
        left_values = left_values[left_cut:]
        left_positions = left_positions[left_cut:]
        right_values = right_values[right_cut:]
        right_positions = right_positions[right_cut:]

    if out_cursor != out_total:
        raise RuntimeError(
            "sorted index OOC merge produced an unexpected output length: "
            f"merge_id={merge_id}, left_len={left.length}, right_len={right.length}, "
            f"expected={out_total}, written={out_cursor}"
        )

    del out_values, out_positions
    _tracker_register_create(tracker, out_values_path, out_positions_path)
    del left_values_mm, left_positions_mm, right_values_mm, right_positions_mm
    _tracker_register_delete(
        tracker, left.values_path, left.positions_path, right.values_path, right.positions_path
    )
    left.values_path.unlink(missing_ok=True)
    left.positions_path.unlink(missing_ok=True)
    right.values_path.unlink(missing_ok=True)
    right.positions_path.unlink(missing_ok=True)
    return SortedRun(out_values_path, out_positions_path, out_cursor)


@dataclass
class OpsiStageSidecars:
    values: blosc2.NDArray | np.ndarray | None
    positions: blosc2.NDArray | np.ndarray | None
    mins: blosc2.NDArray | np.ndarray | None
    medians: blosc2.NDArray | np.ndarray | None
    maxs: blosc2.NDArray | np.ndarray | None
    values_path: str | None
    positions_path: str | None
    mins_path: str | None
    medians_path: str | None
    maxs_path: str | None
    chunk_len: int
    block_len: int
    nblocks: int


def _opsi_stage_create(
    array: blosc2.NDArray,
    token: str,
    kind: str,
    cycle: int,
    size: int,
    nblocks: int,
    dtype: np.dtype,
    persistent: bool,
    chunk_len: int,
    block_len: int,
    cparams: dict | blosc2.CParams | None = None,
) -> OpsiStageSidecars:
    category = f"opsi_cycle{cycle}"
    if persistent:
        values, values_info = _create_persistent_sidecar_handle(
            array,
            token,
            kind,
            category,
            "values",
            size,
            dtype,
            chunks=(chunk_len,),
            blocks=(block_len,),
            cparams=cparams,
        )
        positions, positions_info = _create_persistent_sidecar_handle(
            array,
            token,
            kind,
            category,
            "positions",
            size,
            np.dtype(np.int64),
            chunks=(chunk_len,),
            blocks=(block_len,),
            cparams=cparams,
        )
        minmax_chunk = max(1, min(nblocks, max(1, chunk_len // block_len)))
        mins, mins_info = _create_persistent_sidecar_handle(
            array,
            token,
            kind,
            category,
            "mins",
            nblocks,
            dtype,
            chunks=(minmax_chunk,),
            blocks=(minmax_chunk,),
            cparams=cparams,
        )
        medians, medians_info = _create_persistent_sidecar_handle(
            array,
            token,
            kind,
            category,
            "medians",
            nblocks,
            dtype,
            chunks=(minmax_chunk,),
            blocks=(minmax_chunk,),
            cparams=cparams,
        )
        maxs, maxs_info = _create_persistent_sidecar_handle(
            array,
            token,
            kind,
            category,
            "maxs",
            nblocks,
            dtype,
            chunks=(minmax_chunk,),
            blocks=(minmax_chunk,),
            cparams=cparams,
        )
        return OpsiStageSidecars(
            values,
            positions,
            mins,
            medians,
            maxs,
            values_info["path"],
            positions_info["path"],
            mins_info["path"],
            medians_info["path"],
            maxs_info["path"],
            chunk_len,
            block_len,
            nblocks,
        )

    return OpsiStageSidecars(
        np.empty(size, dtype=dtype),
        np.empty(size, dtype=np.int64),
        np.empty(nblocks, dtype=dtype),
        np.empty(nblocks, dtype=dtype),
        np.empty(nblocks, dtype=dtype),
        None,
        None,
        None,
        None,
        None,
        chunk_len,
        block_len,
        nblocks,
    )


def _opsi_drop_stage(stage: OpsiStageSidecars | None, *, keep_payload: bool = False) -> None:
    if stage is None:
        return
    values_path = stage.values_path
    positions_path = stage.positions_path
    mins_path = stage.mins_path
    medians_path = stage.medians_path
    maxs_path = stage.maxs_path
    # Release writable B2ND handles before removing their paths.  Otherwise,
    # when BLOSC_TRACE is enabled, c-blosc2 may report noisy short-write errors
    # from handle finalizers trying to flush files that have already been
    # removed.
    stage.values = None
    stage.positions = None
    stage.mins = None
    stage.medians = None
    stage.maxs = None
    if not keep_payload:
        _remove_sidecar_path(values_path)
        _remove_sidecar_path(positions_path)
    _remove_sidecar_path(mins_path)
    _remove_sidecar_path(medians_path)
    _remove_sidecar_path(maxs_path)


def _opsi_sort_values_positions(values: np.ndarray, positions: np.ndarray) -> None:
    try:
        indexing_ext.keysort_values_positions(values, positions)
    except (AttributeError, TypeError):
        order = np.lexsort((positions, values))
        values[...] = values[order]
        positions[...] = positions[order]


def _opsi_write_block_boundaries(
    stage: OpsiStageSidecars, chunk_start: int, chunk_values: np.ndarray, size: int
) -> None:
    if chunk_values.size == 0:
        return
    block_len = stage.block_len
    first_block = chunk_start // block_len
    nblocks = math.ceil(chunk_values.size / block_len)
    mins = np.empty(nblocks, dtype=chunk_values.dtype)
    medians = np.empty(nblocks, dtype=chunk_values.dtype)
    maxs = np.empty(nblocks, dtype=chunk_values.dtype)
    for idx in range(nblocks):
        local_start = idx * block_len
        local_stop = min(local_start + block_len, chunk_values.size)
        mins[idx] = chunk_values[local_start]
        medians[idx] = chunk_values[local_start + (local_stop - local_start - 1) // 2]
        maxs[idx] = chunk_values[local_stop - 1]
    _write_ndarray_linear_span(stage.mins, first_block, mins)
    _write_ndarray_linear_span(stage.medians, first_block, medians)
    _write_ndarray_linear_span(stage.maxs, first_block, maxs)


def _index_chunk_multiplier_for_optlevel(optlevel: int) -> int:
    """Return the chunk-local index chunk multiplier for *optlevel*."""
    if optlevel <= 1:
        return 1
    if optlevel <= 6:
        return 2
    if optlevel == 9:
        return 4
    return 3


def _opsi_max_cycles_for_optlevel(optlevel: int) -> int:
    """Return the OPSI max cycles for *optlevel*."""
    if optlevel <= 1:
        return 1
    if optlevel <= 6:
        return 2
    return 3


def _opsi_storage_chunk_len(chunk_len: int, block_len: int, multiplier: int = 1) -> int:
    """Return an OPSI sidecar chunk length aligned to whole OPSI blocks."""
    chunk_len = max(1, int(chunk_len))
    block_len = max(1, int(block_len))
    multiplier = max(1, int(multiplier))
    if chunk_len % block_len == 0:
        aligned_chunk_len = chunk_len
    elif chunk_len >= block_len:
        aligned_chunk_len = (chunk_len // block_len) * block_len
    else:
        aligned_chunk_len = block_len
    return max(block_len, aligned_chunk_len * multiplier)


def _opsi_build_stage1(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    persistent: bool,
    cycle: int,
    previous: OpsiStageSidecars | None,
    block_order: np.ndarray | None,
    cparams: dict | blosc2.CParams | None = None,
    chunk_multiplier: int = 1,
) -> OpsiStageSidecars:
    size = int(array.shape[0])
    source_chunk_len = int(array.chunks[0])
    block_len = int(array.blocks[0])
    chunk_len = _opsi_storage_chunk_len(source_chunk_len, block_len, chunk_multiplier)
    nblocks = math.ceil(size / block_len)
    stage = _opsi_stage_create(
        array, token, kind, cycle, size, nblocks, dtype, persistent, chunk_len, block_len, cparams
    )

    for chunk_start in range(0, size, chunk_len):
        chunk_stop = min(chunk_start + chunk_len, size)
        chunk_size = chunk_stop - chunk_start
        if previous is None:
            chunk_values = _slice_values_for_target(array, target, chunk_start, chunk_stop).astype(
                dtype, copy=True
            )
            chunk_positions = np.arange(chunk_start, chunk_stop, dtype=np.int64)
        else:
            chunk_values = np.empty(chunk_size, dtype=dtype)
            chunk_positions = np.empty(chunk_size, dtype=np.int64)
            local = 0
            while local < chunk_size:
                dst_global = chunk_start + local
                dst_block_id = dst_global // block_len
                src_block_id = int(block_order[dst_block_id])
                src_start = src_block_id * block_len
                src_stop = min(src_start + block_len, size)
                take = min(src_stop - src_start, chunk_size - local)
                _read_ndarray_linear_span(previous.values, src_start, chunk_values[local : local + take])
                _read_ndarray_linear_span(
                    previous.positions, src_start, chunk_positions[local : local + take]
                )
                local += take
        _opsi_sort_values_positions(chunk_values, chunk_positions)
        _write_ndarray_linear_span(stage.values, chunk_start, chunk_values)
        _write_ndarray_linear_span(stage.positions, chunk_start, chunk_positions)
        _opsi_write_block_boundaries(stage, chunk_start, chunk_values, size)
    return stage


def _opsi_ordered_le_array(left: np.ndarray, right: np.ndarray, dtype: np.dtype) -> np.ndarray:
    dtype = np.dtype(dtype)
    if dtype.kind != "f":
        return left <= right
    left_nan = np.isnan(left)
    right_nan = np.isnan(right)
    return (~left_nan & right_nan) | (left_nan & right_nan) | (~left_nan & ~right_nan & (left <= right))


def _opsi_is_csi(stage: OpsiStageSidecars, dtype: np.dtype) -> bool:
    if stage.nblocks <= 1:
        return True
    maxs = _read_sidecar_span(stage.maxs, 0, stage.nblocks - 1)
    mins = _read_sidecar_span(stage.mins, 1, stage.nblocks)
    return bool(np.all(_opsi_ordered_le_array(maxs, mins, dtype)))


def _opsi_block_pruning_score(stage: OpsiStageSidecars, dtype: np.dtype) -> tuple[float, int, float, float]:
    """Estimate point-query block lookup cost; lower is better.

    A candidate block count alone is not enough for out-of-core lookups: many
    isolated candidate blocks can be slower than a slightly larger but more
    contiguous set because each run becomes a separate sidecar read.  Score both
    total candidate blocks and run fragmentation, using a small read-run penalty
    so higher optlevels do not select layouts that prune more blocks but scatter
    the remaining candidates across the file.
    """
    if stage.nblocks <= 1:
        return (1.0, 1, 1.0, 1.0)
    mins = _read_sidecar_span(stage.mins, 0, stage.nblocks)
    medians = _read_sidecar_span(stage.medians, 0, stage.nblocks)
    maxs = _read_sidecar_span(stage.maxs, 0, stage.nblocks)
    sample_count = min(129, stage.nblocks)
    try:
        samples = np.sort(medians, kind="stable")
    except TypeError:
        samples = medians
    if sample_count < stage.nblocks:
        sample_ids = np.linspace(0, stage.nblocks - 1, sample_count, dtype=np.int64)
        samples = samples[sample_ids]
    costs = []
    counts = []
    runs = []
    read_run_penalty = 8.0
    for sample in samples:
        candidates = (mins <= sample) & (maxs >= sample)
        true_ids = np.flatnonzero(candidates)
        count = int(true_ids.size)
        if count == 0:
            run_count = 0
        else:
            run_count = int(np.count_nonzero(np.diff(true_ids) != 1) + 1)
        counts.append(count)
        runs.append(run_count)
        costs.append(count + read_run_penalty * run_count)
    nsamples = max(1, len(costs))
    return (
        float(sum(costs)) / nsamples,
        max(costs, default=0),
        float(sum(counts)) / nsamples,
        float(sum(runs)) / nsamples,
    )


def _opsi_compute_block_order(stage: OpsiStageSidecars, key_kind: str, size: int) -> np.ndarray:
    # Whole-block relocation assumes fixed-size destination block slots.  Keep
    # a final partial block, when present, in the final slot for now.
    sortable_blocks = stage.nblocks - 1 if size % stage.block_len else stage.nblocks
    if key_kind == "min":
        key_sidecar = stage.mins
    elif key_kind == "median":
        key_sidecar = stage.medians
    elif key_kind == "max":
        key_sidecar = stage.maxs
    else:
        raise ValueError(f"unsupported OPSI block-order key kind {key_kind!r}")
    keys = _read_sidecar_span(key_sidecar, 0, sortable_blocks)
    block_ids = np.arange(sortable_blocks, dtype=np.int64)
    try:
        indexing_ext.keysort_keys_indices(keys, block_ids)
    except (AttributeError, TypeError):
        order = np.lexsort((block_ids, keys))
        block_ids = block_ids[order]
    if sortable_blocks != stage.nblocks:
        block_ids = np.concatenate((block_ids, np.asarray([stage.nblocks - 1], dtype=np.int64)))
    return block_ids


def _build_opsi_descriptor(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    persistent: bool,
    cparams: dict | blosc2.CParams | None = None,
    max_cycles: int = 3,
    optlevel: int = 5,
) -> dict:
    max_cycles = max(0, int(max_cycles))
    chunk_multiplier = _index_chunk_multiplier_for_optlevel(optlevel)
    size = int(array.shape[0])
    if size == 0:
        values_sidecar = _store_array_sidecar(
            array, token, kind, "opsi", "values", np.empty(0, dtype=dtype), persistent, cparams=cparams
        )
        positions_sidecar = _store_array_sidecar(
            array,
            token,
            kind,
            "opsi",
            "positions",
            np.empty(0, dtype=np.int64),
            persistent,
            cparams=cparams,
        )
        mins_sidecar = _store_array_sidecar(
            array, token, kind, "opsi_nav", "mins", np.empty(0, dtype=dtype), persistent, cparams=cparams
        )
        maxs_sidecar = _store_array_sidecar(
            array, token, kind, "opsi_nav", "maxs", np.empty(0, dtype=dtype), persistent, cparams=cparams
        )
        return {
            "values_path": values_sidecar["path"],
            "positions_path": positions_sidecar["path"],
            "mins_path": mins_sidecar["path"],
            "maxs_path": maxs_sidecar["path"],
            "chunk_len": _opsi_storage_chunk_len(
                int(array.chunks[0]), int(array.blocks[0]), chunk_multiplier
            ),
            "block_len": int(array.blocks[0]),
            "chunk_multiplier": chunk_multiplier,
            "nblocks": 0,
            "cycles": 0,
            "max_cycles": max_cycles,
            "is_csi": True,
        }

    previous = None
    block_order = None
    current = None
    best = None
    best_cycle = 0
    best_is_csi = False
    best_score = (math.inf, math.inf, math.inf, math.inf)
    no_improve_patience = 3
    try:
        for cycle in range(max_cycles):
            old_previous = previous
            current = _opsi_build_stage1(
                array,
                target,
                token,
                kind,
                dtype,
                persistent,
                cycle,
                previous,
                block_order,
                cparams,
                chunk_multiplier,
            )
            is_csi = _opsi_is_csi(current, dtype)
            score = (1.0, 1, 1.0, 1.0) if is_csi else _opsi_block_pruning_score(current, dtype)
            if score < best_score:
                if best is not None and best is not old_previous and best is not current:
                    _opsi_drop_stage(best)
                best = current
                best_cycle = cycle + 1
                best_is_csi = is_csi
                best_score = score

            exhausted_budget = cycle == max_cycles - 1
            stalled = (cycle + 1 - best_cycle) >= no_improve_patience
            if is_csi or exhausted_budget or stalled:
                final = best if best is not None else current
                mins = _read_sidecar_span(final.mins, 0, final.nblocks)
                maxs = _read_sidecar_span(final.maxs, 0, final.nblocks)
                nav_chunk = max(1, min(final.nblocks, final.chunk_len // final.block_len))
                mins_sidecar = _store_array_sidecar(
                    array,
                    token,
                    kind,
                    "opsi_nav",
                    "mins",
                    mins,
                    persistent,
                    chunks=(nav_chunk,),
                    blocks=(nav_chunk,),
                    cparams=cparams,
                )
                maxs_sidecar = _store_array_sidecar(
                    array,
                    token,
                    kind,
                    "opsi_nav",
                    "maxs",
                    maxs,
                    persistent,
                    chunks=(nav_chunk,),
                    blocks=(nav_chunk,),
                    cparams=cparams,
                )
                opsi = {
                    "values_path": final.values_path,
                    "positions_path": final.positions_path,
                    "mins_path": mins_sidecar["path"],
                    "maxs_path": maxs_sidecar["path"],
                    "chunk_len": final.chunk_len,
                    "block_len": final.block_len,
                    "chunk_multiplier": chunk_multiplier,
                    "nblocks": final.nblocks,
                    "cycles": best_cycle,
                    "attempted_cycles": cycle + 1,
                    "max_cycles": max_cycles,
                    "is_csi": best_is_csi,
                    "block_pruning_score": best_score[0],
                    "block_pruning_blocks": best_score[2],
                    "block_pruning_runs": best_score[3],
                }
                if not persistent:
                    values = np.asarray(final.values)
                    positions = np.asarray(final.positions)
                    values_sidecar = _store_array_sidecar(
                        array,
                        token,
                        kind,
                        "opsi",
                        "values",
                        values,
                        persistent,
                        chunks=(final.chunk_len,),
                        blocks=(final.block_len,),
                        cparams=cparams,
                    )
                    positions_sidecar = _store_array_sidecar(
                        array,
                        token,
                        kind,
                        "opsi",
                        "positions",
                        positions,
                        persistent,
                        chunks=(final.chunk_len,),
                        blocks=(final.block_len,),
                        cparams=cparams,
                    )
                    opsi["values_path"] = values_sidecar["path"]
                    opsi["positions_path"] = positions_sidecar["path"]
                if old_previous is not final:
                    _opsi_drop_stage(old_previous)
                if current is not final:
                    _opsi_drop_stage(current)
                _opsi_drop_stage(final, keep_payload=persistent)
                return opsi
            key_kind = ("min", "median", "max")[cycle % 3]
            block_order = _opsi_compute_block_order(current, key_kind=key_kind, size=size)
            if old_previous is not best:
                _opsi_drop_stage(old_previous)
            previous = current
            current = None
        raise RuntimeError("OPSI max_cycles must be greater than zero")
    except Exception:
        _opsi_drop_stage(previous)
        _opsi_drop_stage(current)
        if best is not previous and best is not current:
            _opsi_drop_stage(best)
        raise


def _full_ooc_run_item_budget_for_optlevel(optlevel: int) -> int:
    return FULL_OOC_RUN_ITEMS * _index_chunk_multiplier_for_optlevel(optlevel)


def _build_full_descriptor_ooc(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    dtype: np.dtype,
    persistent: bool,
    workdir: Path,
    cparams: dict | None = None,
    optlevel: int = 5,
) -> dict:
    size = int(array.shape[0])
    tracker = TempRunTracker()
    if size == 0:
        sorted_values = np.empty(0, dtype=dtype)
        positions = np.empty(0, dtype=np.int64)
        sidecar_chunks, sidecar_blocks, chunk_multiplier = _full_sidecar_geometry(array, optlevel)
        values_sidecar = _store_array_sidecar(
            array,
            token,
            kind,
            "full",
            "values",
            sorted_values,
            persistent,
            chunks=sidecar_chunks,
            blocks=sidecar_blocks,
            cparams=cparams,
        )
        positions_sidecar = _store_array_sidecar(
            array,
            token,
            kind,
            "full",
            "positions",
            positions,
            persistent,
            chunks=sidecar_chunks,
            blocks=sidecar_blocks,
            cparams=cparams,
        )
        full = {
            "values_path": values_sidecar["path"],
            "positions_path": positions_sidecar["path"],
            "chunk_multiplier": chunk_multiplier,
            "sidecar_chunk_len": sidecar_chunks[0],
            "sidecar_block_len": sidecar_blocks[0],
            "runs": [],
            "next_run_id": 0,
            "ooc_run_items": 0,
            "ooc_run_item_budget": 0,
            "ooc_run_item_budget_source": "empty",
        }
        _rebuild_full_navigation_sidecars(array, token, kind, full, sorted_values, persistent, cparams)
        return full
    run_items_env = os.getenv("BLOSC2_FULL_OOC_RUN_ITEMS")
    run_item_budget_source = "optlevel"
    if run_items_env is not None:
        try:
            run_item_budget = max(1, int(run_items_env))
            run_item_budget_source = "env"
        except ValueError:
            run_item_budget = _full_ooc_run_item_budget_for_optlevel(optlevel)
    else:
        run_item_budget = _full_ooc_run_item_budget_for_optlevel(optlevel)
    run_items = max(int(array.chunks[0]), min(size, run_item_budget))
    runs = []
    for run_id, start in enumerate(range(0, size, run_items)):
        stop = min(start + run_items, size)
        values = _slice_values_for_target(array, target, start, stop)
        sorted_values, sorted_positions = _keysort_values_with_positions(values, start)
        runs.append(
            _materialize_sorted_run(
                sorted_values,
                sorted_positions,
                stop - start,
                dtype,
                workdir,
                f"full_run_{run_id}",
                tracker,
                cparams,
            )
        )

    merge_buffer_items = max(int(array.chunks[0]), FULL_OOC_MERGE_BUFFER_ITEMS)
    merge_id = 0
    while len(runs) > 1:
        next_runs = []
        for idx in range(0, len(runs), 2):
            if idx + 1 >= len(runs):
                next_runs.append(runs[idx])
                continue
            next_runs.append(
                _merge_run_pair(
                    runs[idx],
                    runs[idx + 1],
                    workdir,
                    dtype,
                    merge_id,
                    merge_buffer_items,
                    tracker,
                    cparams,
                )
            )
            merge_id += 1
        runs = next_runs

    final_run = runs[0]
    full = {
        "values_path": None,
        "positions_path": None,
        "runs": [],
        "next_run_id": 0,
        "temp_backend": "blosc2",
        "temp_peak_disk_bytes": tracker.peak_disk_bytes,
        "temp_total_written_bytes": tracker.total_written_bytes,
        "ooc_run_items": run_items,
        "ooc_run_item_budget": run_item_budget,
        "ooc_run_item_budget_source": run_item_budget_source,
    }
    if persistent:
        _stream_copy_temp_run_to_full_sidecars(
            array, token, kind, full, final_run, dtype, persistent, tracker, cparams, optlevel
        )
    else:
        sorted_values = blosc2.open(str(final_run.values_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)[:]
        positions = blosc2.open(str(final_run.positions_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)[:]
        sidecar_chunks, sidecar_blocks, chunk_multiplier = _full_sidecar_geometry(array, optlevel)
        values_sidecar = _store_array_sidecar(
            array,
            token,
            kind,
            "full",
            "values",
            sorted_values,
            persistent,
            chunks=sidecar_chunks,
            blocks=sidecar_blocks,
            cparams=cparams,
        )
        positions_sidecar = _store_array_sidecar(
            array,
            token,
            kind,
            "full",
            "positions",
            positions,
            persistent,
            chunks=sidecar_chunks,
            blocks=sidecar_blocks,
            cparams=cparams,
        )
        full["values_path"] = values_sidecar["path"]
        full["positions_path"] = positions_sidecar["path"]
        full["chunk_multiplier"] = chunk_multiplier
        full["sidecar_chunk_len"] = sidecar_chunks[0]
        full["sidecar_block_len"] = sidecar_blocks[0]
        _rebuild_full_navigation_sidecars(array, token, kind, full, sorted_values, persistent, cparams)
        del sorted_values, positions
        _tracker_register_delete(tracker, final_run.values_path, final_run.positions_path)
        final_run.values_path.unlink(missing_ok=True)
        final_run.positions_path.unlink(missing_ok=True)
    return full


def _build_descriptor(
    array: blosc2.NDArray,
    target: dict,
    token: str,
    kind: str,
    optlevel: int,
    persistent: bool,
    ooc: bool,
    name: str | None,
    dtype: np.dtype,
    levels: dict,
    bucket: dict | None,
    partial: dict | None,
    full: dict | None,
    cparams: dict | None = None,
    opsi: dict | None = None,
) -> dict:
    return {
        "name": name
        or (target["expression"] if target["source"] == "expression" else _field_token(target.get("field"))),
        "token": token,
        "target": target.copy(),
        "field": _target_field(target),
        "kind": kind,
        "version": INDEX_FORMAT_VERSION,
        "optlevel": optlevel,
        "persistent": persistent,
        "ooc": ooc,
        "stale": False,
        "dtype": np.dtype(dtype).str,
        "shape": tuple(array.shape),
        "chunks": tuple(array.chunks),
        "blocks": tuple(array.blocks),
        "levels": levels,
        "bucket": bucket,
        "partial": partial,
        "full": full,
        "opsi": opsi,
        "cparams": _plain_index_cparams(cparams),
    }


def create_index(
    array: blosc2.NDArray,
    field: str | None = None,
    *,
    expression: str | None = None,
    operands: dict | None = None,
    kind: blosc2.IndexKind = blosc2.IndexKind.BUCKET,
    optlevel: int = 5,
    persistent: bool | None = None,
    build: str = "auto",
    name: str | None = None,
    tmpdir: str | None = None,
    **kwargs,
) -> dict:
    """Create an index descriptor for a 1-D array, field, or expression.

    Parameters are equivalent to :meth:`blosc2.NDArray.create_index`.
    ``tmpdir`` controls where temporary files for out-of-core ``kind=FULL``
    builds are created. If ``None``, persistent arrays default to their own
    directory and in-memory arrays use the system temporary directory.
    """
    cparams = _normalize_index_cparams(kwargs.pop("cparams", None))
    method = kwargs.pop("method", None)
    opsi_max_cycles_arg = kwargs.pop("opsi_max_cycles", None)
    summary_levels = kwargs.pop("summary_levels", None)
    precomputed_summaries = kwargs.pop("precomputed_summaries", None)
    if kwargs:
        unexpected = ", ".join(sorted(kwargs))
        raise TypeError(f"unexpected keyword argument(s): {unexpected}")
    if not isinstance(kind, blosc2.IndexKind):
        raise TypeError("kind must be a blosc2.IndexKind")
    kind = _normalize_index_kind(kind)
    build = _normalize_build_mode(build)
    if opsi_max_cycles_arg is None:
        opsi_max_cycles = _opsi_max_cycles_for_optlevel(optlevel)
    else:
        opsi_max_cycles = max(1, int(opsi_max_cycles_arg))
    if kind == "full":
        _normalize_full_build_method(method)
    if method is not None and kind != "full":
        raise ValueError("method is only supported for kind=IndexKind.FULL")
    target, dtype = _normalize_create_index_target(array, field, expression, operands)
    token = _target_token(target)
    if persistent is None:
        persistent = _is_persistent_array(array)
    use_ooc = _resolve_ooc_mode(kind, build)

    if use_ooc:
        levels = _build_levels_descriptor_ooc(
            array,
            target,
            token,
            kind,
            dtype,
            persistent,
            cparams,
            summary_levels=summary_levels,
            precomputed_summaries=precomputed_summaries,
        )
        bucket = (
            _build_bucket_descriptor_ooc(array, target, token, kind, dtype, optlevel, persistent, cparams)
            if kind == "bucket"
            else None
        )
        partial = (
            _build_partial_descriptor_ooc(array, target, token, kind, dtype, optlevel, persistent, cparams)
            if kind == "partial"
            else None
        )
        full = None
        opsi = None
        if kind == "full":
            with tempfile.TemporaryDirectory(
                prefix="blosc2-index-ooc-", dir=_resolve_full_index_tmpdir(array, tmpdir)
            ) as tmpdir:
                full = _build_full_descriptor_ooc(
                    array, target, token, kind, dtype, persistent, Path(tmpdir), cparams, optlevel
                )
                full["build_method"] = "global-sort"
        if kind == "opsi":
            opsi = _build_opsi_descriptor(
                array, target, token, kind, dtype, persistent, cparams, opsi_max_cycles, optlevel
            )
        descriptor = _build_descriptor(
            array,
            target,
            token,
            kind,
            optlevel,
            persistent,
            True,
            name,
            dtype,
            levels,
            bucket,
            partial,
            full,
            cparams,
            opsi,
        )
    else:
        values = _values_for_target(array, target)
        levels = _build_levels_descriptor(
            array, target, token, kind, dtype, values, persistent, cparams, summary_levels=summary_levels
        )
        bucket = (
            _build_bucket_descriptor(array, token, kind, values, optlevel, persistent, cparams)
            if kind == "bucket"
            else None
        )
        partial = (
            _build_partial_descriptor(array, token, kind, values, optlevel, persistent, cparams)
            if kind == "partial"
            else None
        )
        full = None
        opsi = None
        if kind == "full":
            full = _build_full_descriptor(array, token, kind, values, persistent, cparams, optlevel)
            full["build_method"] = "global-sort"
        if kind == "opsi":
            opsi = _build_opsi_descriptor(
                array, target, token, kind, dtype, persistent, cparams, opsi_max_cycles, optlevel
            )
        descriptor = _build_descriptor(
            array,
            target,
            token,
            kind,
            optlevel,
            persistent,
            False,
            name,
            dtype,
            levels,
            bucket,
            partial,
            full,
            cparams,
            opsi,
        )

    store = _load_store(array)
    store["indexes"][token] = descriptor
    _save_store(array, store)
    return _copy_descriptor(descriptor)


def _resolve_index_token(store: dict, field: str | None, name: str | None) -> str:
    token = None
    if field is not None:
        token = _field_token(field)
    elif name is None and len(store["indexes"]) == 1:
        token = next(iter(store["indexes"]))
    if token is None:
        for key, descriptor in store["indexes"].items():
            if descriptor.get("name") == name:
                token = key
                break
    if token is None or token not in store["indexes"]:
        raise KeyError("index not found")
    return token


def iter_index_components(array: blosc2.NDArray, descriptor: dict):
    for level in descriptor["levels"]:
        level_info = descriptor["levels"][level]
        yield IndexComponent(f"summary.{level}", "summary", level, level_info.get("path"))

    bucket = descriptor.get("bucket")
    if bucket is not None:
        yield IndexComponent("bucket.values", "bucket", "values", bucket.get("values_path"))
        yield IndexComponent(
            "bucket.bucket_positions", "bucket", "bucket_positions", bucket.get("bucket_positions_path")
        )
        yield IndexComponent("bucket.offsets", "bucket", "offsets", bucket.get("offsets_path"))
        yield IndexComponent("bucket_nav.l1", "bucket_nav", "l1", bucket.get("l1_path"))
        yield IndexComponent("bucket_nav.l2", "bucket_nav", "l2", bucket.get("l2_path"))

    partial = descriptor.get("partial")
    if partial is not None:
        yield IndexComponent("partial.values", "partial", "values", partial.get("values_path"))
        yield IndexComponent("partial.positions", "partial", "positions", partial.get("positions_path"))
        yield IndexComponent("partial.offsets", "partial", "offsets", partial.get("offsets_path"))
        yield IndexComponent("partial_nav.l1", "partial_nav", "l1", partial.get("l1_path"))
        yield IndexComponent("partial_nav.l2", "partial_nav", "l2", partial.get("l2_path"))

    full = descriptor.get("full")
    if full is not None:
        yield IndexComponent("full.values", "full", "values", full.get("values_path"))
        yield IndexComponent("full.positions", "full", "positions", full.get("positions_path"))
        yield IndexComponent("full_nav.l1", "full_nav", "l1", full.get("l1_path"))
        yield IndexComponent("full_nav.l2", "full_nav", "l2", full.get("l2_path"))
        for run in full.get("runs", ()):
            run_id = int(run["id"])
            yield IndexComponent(
                f"full_run.{run_id}.values",
                "full_run",
                f"{run_id}.values",
                run.get("values_path"),
            )
            yield IndexComponent(
                f"full_run.{run_id}.positions",
                "full_run",
                f"{run_id}.positions",
                run.get("positions_path"),
            )

    opsi = descriptor.get("opsi")
    if opsi is not None:
        yield IndexComponent("opsi.values", "opsi", "values", opsi.get("values_path"))
        yield IndexComponent("opsi.positions", "opsi", "positions", opsi.get("positions_path"))
        yield IndexComponent("opsi_nav.mins", "opsi_nav", "mins", opsi.get("mins_path"))
        yield IndexComponent("opsi_nav.maxs", "opsi_nav", "maxs", opsi.get("maxs_path"))


def _component_nbytes(array: blosc2.NDArray, descriptor: dict, component: IndexComponent) -> int:
    if component.path is not None:
        return int(_open_sidecar_file(component.path, _INDEX_MMAP_MODE).nbytes)
    token = descriptor["token"]
    return int(_load_array_sidecar(array, token, component.category, component.name, component.path).nbytes)


def _component_cbytes(array: blosc2.NDArray, descriptor: dict, component: IndexComponent) -> int:
    if component.path is not None:
        return int(_open_sidecar_file(component.path, _INDEX_MMAP_MODE).cbytes)
    token = descriptor["token"]
    sidecar = _load_array_sidecar(array, token, component.category, component.name, component.path)
    kwargs = {}
    cparams = descriptor.get("cparams")
    if cparams is not None:
        kwargs["cparams"] = cparams
    return int(blosc2.asarray(sidecar, **kwargs).cbytes)


class Index(Mapping):
    """Handle for an index attached to an :class:`blosc2.NDArray`.

    ``Index`` objects are returned by NDArray indexing helpers such as
    :meth:`blosc2.NDArray.create_index`, :meth:`blosc2.NDArray.index`, and
    :attr:`blosc2.NDArray.indexes`, and by the equivalent :class:`blosc2.CTable`
    helpers.  They expose descriptor metadata and convenience methods for
    dropping, rebuilding, or compacting the underlying index.  Users should not
    instantiate this class directly.
    """

    def __init__(
        self,
        array: blosc2.NDArray | None,
        token: str,
        *,
        table=None,
        descriptor: dict | None = None,
    ):
        self._array = array
        self._token = token
        self._table = table
        self._descriptor_snapshot = _copy_descriptor(descriptor) if descriptor is not None else None

    @classmethod
    def _from_table(cls, table, lookup_key: str, descriptor: dict) -> Index:
        """Create an index handle for a CTable catalog entry."""
        return cls(None, lookup_key, table=table, descriptor=descriptor)

    def _is_table_index(self) -> bool:
        return self._table is not None

    def _descriptor(self) -> dict:
        if self._table is None:
            return _descriptor_for_token(self._array, self._token)
        _, descriptor = self._table._root_table._resolve_index_catalog_entry(self._token)
        self._descriptor_snapshot = _copy_descriptor(descriptor)
        return descriptor

    def _target_array(self) -> blosc2.NDArray:
        if self._table is None:
            return self._array
        return self._table._root_table._index_target_array(self._token, self._descriptor())

    @property
    def descriptor(self) -> dict:
        """Copy of the index descriptor dictionary."""
        if self._table is None:
            return _copy_descriptor_for_token(self._array, self._token)
        return _copy_descriptor(self._descriptor())

    @property
    def kind(self) -> blosc2.IndexKind | str:
        """Kind of index.

        NDArray indexes return :class:`blosc2.IndexKind`; table indexes keep the
        catalog kind string used by CTable descriptors.
        """
        kind = self._descriptor()["kind"]
        if self._table is not None:
            return kind
        return blosc2.IndexKind(kind)

    @property
    def col_name(self) -> str | None:
        """CTable lookup key for this index target, if this is a table index."""
        return self._token if self._table is not None else None

    @property
    def field(self) -> str | None:
        """Structured-array field indexed by this handle, if any."""
        return self._descriptor().get("field")

    @property
    def name(self) -> str | None:
        """Optional human-readable name assigned at creation time."""
        return self._descriptor()["name"]

    @property
    def target(self) -> dict:
        """Descriptor of the indexed target."""
        return self.descriptor["target"]

    @property
    def persistent(self) -> bool:
        """True when index sidecars are stored persistently on disk."""
        return bool(self._descriptor()["persistent"])

    @property
    def stale(self) -> bool:
        """True if the index needs rebuilding before it can be reused."""
        return bool(self._descriptor()["stale"])

    @property
    def nbytes(self) -> int:
        """Total uncompressed size in bytes for this index payload."""
        descriptor = self._descriptor()
        array = self._target_array()
        return sum(
            _component_nbytes(array, descriptor, component)
            for component in iter_index_components(array, descriptor)
        )

    @property
    def cbytes(self) -> int:
        """Total compressed size in bytes for this index payload."""
        descriptor = self._descriptor()
        array = self._target_array()
        return sum(
            _component_cbytes(array, descriptor, component)
            for component in iter_index_components(array, descriptor)
        )

    @property
    def cratio(self) -> float:
        """Compression ratio for this index payload."""
        cbytes = self.cbytes
        if cbytes == 0:
            return math.inf
        return self.nbytes / cbytes

    def storage_stats(self) -> tuple[int, int, float] | None:
        """Return ``(nbytes, cbytes, cratio)`` when sidecars are directly measurable."""
        try:
            nbytes = self.nbytes
            cbytes = self.cbytes
        except (OSError, RuntimeError, KeyError, ValueError):
            if self._table is None:
                return None
            root = self._table._root_table
            storage = getattr(root, "_storage", None)
            if storage.__class__.__name__ != "FileTableStorage":
                return None

            descriptor = self._descriptor_snapshot or self.descriptor
            target_arr = root._index_target_array(self._token, descriptor)
            store = storage._open_store()
            nbytes = 0
            cbytes = 0
            try:
                for component in iter_index_components(target_arr, descriptor):
                    if component.path is None:
                        return None
                    key = self._component_store_key(component.path)
                    obj = store[key]
                    nbytes += int(obj.nbytes)
                    cbytes += int(obj.cbytes)
            except (OSError, RuntimeError, KeyError, ValueError):
                return None
        cratio = math.inf if cbytes == 0 else nbytes / cbytes
        return nbytes, cbytes, cratio

    @staticmethod
    def _component_store_key(path: str) -> str:
        """Return the logical TreeStore key for an index component path."""
        normalized = path.replace("\\", "/")
        marker = "_indexes/"
        idx = normalized.find(marker)
        if idx < 0:
            raise KeyError(f"Cannot resolve index component path {path!r} inside table store.")
        relpath = normalized[idx:]
        for suffix in (".b2nd", ".b2f"):
            if relpath.endswith(suffix):
                relpath = relpath[: -len(suffix)]
                break
        return "/" + relpath.lstrip("/")

    def drop(self) -> None:
        """Drop this index and delete its sidecar payloads."""
        if self._table is None:
            drop_index(self._array, field=self.field, name=self.name)
            return
        target = self._descriptor().get("target") or {}
        if target.get("source") == "expression":
            self._table.drop_index(expression=target.get("expression"), name=self.name)
        else:
            self._table.drop_index(self._token)

    def rebuild(self) -> Index:
        """Rebuild this index and return the updated handle."""
        if self._table is None:
            rebuild_index(self._array, field=self.field, name=self.name)
            return self
        target = self._descriptor().get("target") or {}
        if target.get("source") == "expression":
            return self._table.rebuild_index(expression=target.get("expression"), name=self.name)
        return self._table.rebuild_index(self._token)

    def compact(self) -> Index:
        """Compact this index, merging incremental runs, and return the updated handle."""
        if self._table is None:
            compact_index(self._array, field=self.field, name=self.name)
            return self
        target = self._descriptor().get("target") or {}
        if target.get("source") == "expression":
            return self._table.compact_index(expression=target.get("expression"), name=self.name)
        return self._table.compact_index(self._token)

    def __getitem__(self, key):
        """Return a descriptor value by key."""
        return self.descriptor[key]

    def __iter__(self):
        """Iterate over descriptor keys."""
        return iter(self.descriptor)

    def __len__(self) -> int:
        """Number of keys in the descriptor mapping."""
        return len(self.descriptor)

    def __repr__(self) -> str:
        """Return a concise representation of this index handle."""
        try:
            descriptor = self._descriptor()
        except KeyError:
            return "Index(<dropped>)"
        if self._table is not None:
            target = descriptor.get("target") or {}
            target_label = self._token if target.get("source") != "expression" else target.get("expression")
            return (
                f"Index(kind={descriptor['kind']!r}, col_name={target_label!r}, "
                f"name={descriptor.get('name')!r}, stale={descriptor.get('stale')!r})"
            )
        return (
            f"Index(kind={descriptor['kind']!r}, field={descriptor.get('field')!r}, "
            f"name={descriptor.get('name')!r}, stale={descriptor.get('stale')!r})"
        )


def _remove_sidecar_path(path: str | None) -> None:
    if path:
        blosc2.remove_urlpath(path)


def _drop_descriptor_sidecars(descriptor: dict) -> None:
    for level_info in descriptor["levels"].values():
        _remove_sidecar_path(level_info["path"])
    if descriptor.get("bucket") is not None:
        _remove_sidecar_path(descriptor["bucket"]["values_path"])
        _remove_sidecar_path(descriptor["bucket"]["bucket_positions_path"])
        _remove_sidecar_path(descriptor["bucket"]["offsets_path"])
        _remove_sidecar_path(descriptor["bucket"].get("l1_path"))
        _remove_sidecar_path(descriptor["bucket"].get("l2_path"))
    if descriptor.get("partial") is not None:
        _remove_sidecar_path(descriptor["partial"]["values_path"])
        _remove_sidecar_path(descriptor["partial"]["positions_path"])
        _remove_sidecar_path(descriptor["partial"]["offsets_path"])
        _remove_sidecar_path(descriptor["partial"].get("l1_path"))
        _remove_sidecar_path(descriptor["partial"].get("l2_path"))
    if descriptor.get("full") is not None:
        _remove_sidecar_path(descriptor["full"]["values_path"])
        _remove_sidecar_path(descriptor["full"]["positions_path"])
        _remove_sidecar_path(descriptor["full"].get("l1_path"))
        _remove_sidecar_path(descriptor["full"].get("l2_path"))
        for run in descriptor["full"].get("runs", ()):
            _remove_sidecar_path(run.get("values_path"))
            _remove_sidecar_path(run.get("positions_path"))
    if descriptor.get("opsi") is not None:
        _remove_sidecar_path(descriptor["opsi"]["values_path"])
        _remove_sidecar_path(descriptor["opsi"]["positions_path"])
        _remove_sidecar_path(descriptor["opsi"].get("mins_path"))
        _remove_sidecar_path(descriptor["opsi"].get("maxs_path"))


def _replace_levels_descriptor_tail(
    array: blosc2.NDArray, descriptor: dict, kind: str, old_size: int, persistent: bool
) -> None:
    target = descriptor["target"]
    token = descriptor["token"]
    dtype = np.dtype(descriptor["dtype"])
    new_size = int(array.shape[0])
    cparams = _normalize_index_cparams(descriptor.get("cparams"))
    for level, level_info in descriptor["levels"].items():
        segment_len = int(level_info["segment_len"])
        start_segment = old_size // segment_len
        prefix = _open_level_summary_handle(array, descriptor, level)[:start_segment]
        tail_start = start_segment * segment_len
        tail_values = _slice_values_for_target(array, target, tail_start, new_size)
        tail_summaries = _compute_segment_summaries(tail_values, dtype, segment_len)
        summaries = np.concatenate((prefix, tail_summaries)) if len(prefix) else tail_summaries
        sidecar = _store_array_sidecar(
            array, token, kind, "summary", level, summaries, persistent, cparams=cparams
        )
        level_info["path"] = sidecar["path"]
        level_info["dtype"] = sidecar["dtype"]
        level_info["nsegments"] = len(summaries)


def _replace_partial_descriptor_tail(
    array: blosc2.NDArray, descriptor: dict, old_size: int, persistent: bool
) -> None:
    del old_size
    target = descriptor["target"]
    partial = descriptor["partial"]
    cparams = _normalize_index_cparams(descriptor.get("cparams"))
    for key in ("values_path", "positions_path", "offsets_path", "l1_path", "l2_path"):
        _remove_sidecar_path(partial.get(key))
    if descriptor.get("ooc", False):
        rebuilt = _build_partial_descriptor_ooc(
            array,
            target,
            descriptor["token"],
            descriptor["kind"],
            np.dtype(descriptor["dtype"]),
            descriptor["optlevel"],
            persistent,
            cparams,
        )
    else:
        rebuilt = _build_partial_descriptor(
            array,
            descriptor["token"],
            descriptor["kind"],
            _values_for_target(array, target),
            descriptor["optlevel"],
            persistent,
            cparams,
        )
    descriptor["partial"] = rebuilt


def _replace_bucket_descriptor_tail(
    array: blosc2.NDArray, descriptor: dict, old_size: int, persistent: bool
) -> None:
    del old_size
    target = descriptor["target"]
    bucket = descriptor["bucket"]
    cparams = _normalize_index_cparams(descriptor.get("cparams"))
    for key in ("values_path", "bucket_positions_path", "offsets_path", "l1_path", "l2_path"):
        _remove_sidecar_path(bucket.get(key))
    if descriptor.get("ooc", False):
        rebuilt = _build_bucket_descriptor_ooc(
            array,
            target,
            descriptor["token"],
            descriptor["kind"],
            np.dtype(descriptor["dtype"]),
            descriptor["optlevel"],
            persistent,
            cparams,
        )
    else:
        rebuilt = _build_bucket_descriptor(
            array,
            descriptor["token"],
            descriptor["kind"],
            _values_for_target(array, target),
            descriptor["optlevel"],
            persistent,
            cparams,
        )
    descriptor["bucket"] = rebuilt


def _replace_full_descriptor(
    array: blosc2.NDArray,
    descriptor: dict,
    sorted_values: np.ndarray,
    positions: np.ndarray,
    persistent: bool,
) -> None:
    kind = descriptor["kind"]
    token = descriptor["token"]
    full = descriptor["full"]
    cparams = _normalize_index_cparams(descriptor.get("cparams"))
    for run in full.get("runs", ()):
        _remove_sidecar_path(run.get("values_path"))
        _remove_sidecar_path(run.get("positions_path"))
    _remove_sidecar_path(full.get("l1_path"))
    _remove_sidecar_path(full.get("l2_path"))
    _clear_cached_data(array, token)
    values_sidecar = _store_array_sidecar(
        array, token, kind, "full", "values", sorted_values, persistent, cparams=cparams
    )
    positions_sidecar = _store_array_sidecar(
        array, token, kind, "full", "positions", positions, persistent, cparams=cparams
    )
    full["values_path"] = values_sidecar["path"]
    full["positions_path"] = positions_sidecar["path"]
    full["runs"] = []
    full["next_run_id"] = 0
    _rebuild_full_navigation_sidecars(array, token, kind, full, sorted_values, persistent, cparams)


def _replace_full_descriptor_from_paths(
    array: blosc2.NDArray,
    descriptor: dict,
    values_path: Path,
    positions_path: Path,
    length: int,
) -> None:
    kind = descriptor["kind"]
    token = descriptor["token"]
    full = descriptor["full"]
    persistent = descriptor["persistent"]
    cparams = _normalize_index_cparams(descriptor.get("cparams"))
    if not persistent:
        raise ValueError("path-based full replacement requires persistent indexes")
    for run in full.get("runs", ()):
        _remove_sidecar_path(run.get("values_path"))
        _remove_sidecar_path(run.get("positions_path"))
    _remove_sidecar_path(full.get("l1_path"))
    _remove_sidecar_path(full.get("l2_path"))
    _clear_cached_data(array, token)
    final_values_path = _sidecar_path(array, token, kind, "full.values")
    final_positions_path = _sidecar_path(array, token, kind, "full.positions")
    _remove_sidecar_path(final_values_path)
    _remove_sidecar_path(final_positions_path)
    _stream_copy_sidecar_array(
        values_path,
        final_values_path,
        length,
        np.dtype(descriptor["dtype"]),
        (int(array.chunks[0]),),
        (int(array.blocks[0]),),
        cparams,
    )
    _stream_copy_sidecar_array(
        positions_path,
        final_positions_path,
        length,
        np.dtype(np.int64),
        (int(array.chunks[0]),),
        (int(array.blocks[0]),),
        cparams,
    )
    values_path.unlink(missing_ok=True)
    positions_path.unlink(missing_ok=True)
    full["values_path"] = final_values_path
    full["positions_path"] = final_positions_path
    full["runs"] = []
    full["next_run_id"] = 0
    _rebuild_full_navigation_sidecars_from_path(
        array,
        token,
        kind,
        full,
        final_values_path,
        np.dtype(descriptor["dtype"]),
        length,
        persistent,
        cparams,
    )


def _store_full_run_descriptor(
    array: blosc2.NDArray,
    descriptor: dict,
    run_id: int,
    sorted_values: np.ndarray,
    positions: np.ndarray,
) -> dict:
    kind = descriptor["kind"]
    token = descriptor["token"]
    persistent = descriptor["persistent"]
    cparams = _normalize_index_cparams(descriptor.get("cparams"))
    values_sidecar = _store_array_sidecar(
        array, token, kind, "full_run", f"{run_id}.values", sorted_values, persistent, cparams=cparams
    )
    positions_sidecar = _store_array_sidecar(
        array, token, kind, "full_run", f"{run_id}.positions", positions, persistent, cparams=cparams
    )
    return {
        "id": run_id,
        "length": len(sorted_values),
        "values_path": values_sidecar["path"],
        "positions_path": positions_sidecar["path"],
    }


def _append_full_descriptor(
    array: blosc2.NDArray, descriptor: dict, old_size: int, appended_values: np.ndarray
) -> None:
    full = descriptor.get("full")
    if full is None:
        raise RuntimeError("full-index metadata is not available")
    appended_positions = np.arange(old_size, old_size + len(appended_values), dtype=np.int64)
    order = np.lexsort((appended_positions, appended_values))
    run_id = int(full.get("next_run_id", 0))
    run = _store_full_run_descriptor(
        array,
        descriptor,
        run_id,
        appended_values[order],
        appended_positions[order],
    )
    runs = list(full.get("runs", ()))
    runs.append(run)
    full["runs"] = runs
    full["next_run_id"] = run_id + 1
    _clear_full_merge_cache(array, descriptor["token"])


def append_to_indexes(array: blosc2.NDArray, old_size: int, appended_values: np.ndarray) -> None:
    store = _load_store(array)
    if not store["indexes"]:
        return

    for descriptor in store["indexes"].values():
        kind = descriptor["kind"]
        persistent = descriptor["persistent"]
        target = descriptor["target"]
        target_values = _values_from_numpy_target(appended_values, target)
        if descriptor.get("stale", False):
            continue
        if kind == "full":
            _append_full_descriptor(array, descriptor, old_size, target_values)
        elif kind == "partial":
            _replace_partial_descriptor_tail(array, descriptor, old_size, persistent)
        elif kind == "bucket":
            _replace_bucket_descriptor_tail(array, descriptor, old_size, persistent)
        _replace_levels_descriptor_tail(array, descriptor, kind, old_size, persistent)
        descriptor["shape"] = tuple(array.shape)
        descriptor["chunks"] = tuple(array.chunks)
        descriptor["blocks"] = tuple(array.blocks)
        descriptor["stale"] = False
    _save_store(array, store)
    _invalidate_query_cache(array)


def drop_index(array: blosc2.NDArray, field: str | None = None, name: str | None = None) -> None:
    store = _load_store(array)
    token = _resolve_index_token(store, field, name)
    descriptor = store["indexes"][token]
    _clear_cached_data(array, descriptor["token"])
    descriptor = store["indexes"].pop(token)
    _save_store(array, store)
    _drop_descriptor_sidecars(descriptor)
    _invalidate_query_cache(array)


def rebuild_index(array: blosc2.NDArray, field: str | None = None, name: str | None = None) -> dict:
    store = _load_store(array)
    token = _resolve_index_token(store, field, name)
    descriptor = store["indexes"][token]
    rebuild_kwargs = {}
    if descriptor["kind"] == "full":
        rebuild_kwargs["method"] = descriptor.get("full", {}).get("build_method")
    if descriptor["kind"] == "opsi":
        rebuild_kwargs["opsi_max_cycles"] = descriptor.get("opsi", {}).get("max_cycles")
    drop_index(array, field=descriptor["field"], name=descriptor["name"])
    if descriptor["target"]["source"] == "expression":
        operands = array.fields if array.dtype.fields is not None else {SELF_TARGET_NAME: array}
        return create_index(
            array,
            expression=descriptor["target"]["expression_key"],
            operands=operands,
            kind=blosc2.IndexKind(descriptor["kind"]),
            optlevel=descriptor["optlevel"],
            persistent=descriptor["persistent"],
            build="ooc" if descriptor.get("ooc", False) else "memory",
            name=descriptor["name"],
            **rebuild_kwargs,
        )
    return create_index(
        array,
        field=descriptor["field"],
        kind=blosc2.IndexKind(descriptor["kind"]),
        optlevel=descriptor["optlevel"],
        persistent=descriptor["persistent"],
        build="ooc" if descriptor.get("ooc", False) else "memory",
        name=descriptor["name"],
        **rebuild_kwargs,
    )


def _full_compaction_runs(array: blosc2.NDArray, descriptor: dict, workdir: Path) -> list[SortedRun]:
    full = descriptor["full"]
    dtype = np.dtype(descriptor["dtype"])
    runs = []
    base_length = int(array.shape[0]) - sum(int(run["length"]) for run in full.get("runs", ()))
    if full["values_path"] is not None and full["positions_path"] is not None:
        base_values_path = _copy_sidecar_to_temp_run(
            full["values_path"], base_length, dtype, workdir, "compact_base_values"
        )
        base_positions_path = _copy_sidecar_to_temp_run(
            full["positions_path"], base_length, np.dtype(np.int64), workdir, "compact_base_positions"
        )
        runs.append(SortedRun(base_values_path, base_positions_path, base_length))
    else:
        values_handle, positions_handle = _load_full_sidecar_handles(array, descriptor)
        base_values_path = _copy_sidecar_handle_to_temp_run(
            values_handle, base_length, dtype, workdir, "compact_base_values"
        )
        base_positions_path = _copy_sidecar_handle_to_temp_run(
            positions_handle, base_length, np.dtype(np.int64), workdir, "compact_base_positions"
        )
        runs.append(SortedRun(base_values_path, base_positions_path, base_length))

    for run in full.get("runs", ()):
        run_length = int(run["length"])
        run_id = int(run["id"])
        if run["values_path"] is not None and run["positions_path"] is not None:
            run_values_path = _copy_sidecar_to_temp_run(
                run["values_path"], run_length, dtype, workdir, f"run_{run_id}_values"
            )
            run_positions_path = _copy_sidecar_to_temp_run(
                run["positions_path"], run_length, np.dtype(np.int64), workdir, f"run_{run_id}_positions"
            )
            runs.append(SortedRun(run_values_path, run_positions_path, run_length))
            continue
        run_values_handle, run_positions_handle = _load_full_run_sidecar_handles(array, descriptor, run)
        run_values_path = _copy_sidecar_handle_to_temp_run(
            run_values_handle, run_length, dtype, workdir, f"run_{run_id}_values"
        )
        run_positions_path = _copy_sidecar_handle_to_temp_run(
            run_positions_handle, run_length, np.dtype(np.int64), workdir, f"run_{run_id}_positions"
        )
        runs.append(SortedRun(run_values_path, run_positions_path, run_length))
    return runs


def compact_index(array: blosc2.NDArray, field: str | None = None, name: str | None = None) -> dict:
    store = _load_store(array)
    token = _resolve_index_token(store, field, name)
    descriptor = store["indexes"][token]
    if descriptor["kind"] != "full":
        raise NotImplementedError("compact_index() is currently only implemented for full indexes")
    if descriptor.get("stale", False):
        raise RuntimeError("cannot compact a stale index; rebuild it first")

    full = descriptor["full"]
    if not full.get("runs"):
        if full.get("l1_path") is None or full.get("l2_path") is None:
            cparams = _normalize_index_cparams(descriptor.get("cparams"))
            dtype = np.dtype(descriptor["dtype"])
            _remove_sidecar_path(full.get("l1_path"))
            _remove_sidecar_path(full.get("l2_path"))
            if descriptor["persistent"] and full.get("values_path") is not None:
                _rebuild_full_navigation_sidecars_from_path(
                    array,
                    descriptor["token"],
                    descriptor["kind"],
                    full,
                    full["values_path"],
                    dtype,
                    int(array.shape[0]),
                    descriptor["persistent"],
                    cparams,
                )
            else:
                values_handle, _ = _load_full_sidecar_handles(array, descriptor)
                _rebuild_full_navigation_sidecars_from_handle(
                    array,
                    descriptor["token"],
                    descriptor["kind"],
                    full,
                    values_handle,
                    dtype,
                    int(array.shape[0]),
                    descriptor["persistent"],
                    cparams,
                )
        _clear_full_merge_cache(array, descriptor["token"])
        _save_store(array, store)
        _invalidate_query_cache(array)
        return _copy_descriptor(descriptor)

    dtype = np.dtype(descriptor["dtype"])
    with tempfile.TemporaryDirectory(prefix="blosc2-index-compact-", dir=_tmpdir_for_array(array)) as tmpdir:
        workdir = Path(tmpdir)
        runs = _full_compaction_runs(array, descriptor, workdir)
        merge_buffer_items = max(int(array.chunks[0]), FULL_OOC_MERGE_BUFFER_ITEMS)
        merge_id = 0
        while len(runs) > 1:
            next_runs = []
            for idx in range(0, len(runs), 2):
                if idx + 1 >= len(runs):
                    next_runs.append(runs[idx])
                    continue
                next_runs.append(
                    _merge_run_pair(runs[idx], runs[idx + 1], workdir, dtype, merge_id, merge_buffer_items)
                )
                merge_id += 1
            runs = next_runs
        final_run = runs[0]
        if descriptor["persistent"]:
            _replace_full_descriptor_from_paths(
                array, descriptor, final_run.values_path, final_run.positions_path, final_run.length
            )
        else:
            sorted_values = blosc2.open(str(final_run.values_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)[:]
            positions = blosc2.open(str(final_run.positions_path), mode="r", mmap_mode=_INDEX_MMAP_MODE)[:]
            _replace_full_descriptor(array, descriptor, sorted_values, positions, descriptor["persistent"])
            del sorted_values, positions
            final_run.values_path.unlink(missing_ok=True)
            final_run.positions_path.unlink(missing_ok=True)

    _clear_full_merge_cache(array, descriptor["token"])
    _save_store(array, store)
    _invalidate_query_cache(array)
    return _copy_descriptor(descriptor)


def get_indexes(array: blosc2.NDArray) -> list[Index]:
    store = _load_store(array)
    return [Index(array, key) for key in sorted(store["indexes"])]


def get_index(array: blosc2.NDArray, field: str | None = None, name: str | None = None) -> Index:
    store = _load_store(array)
    token = _resolve_index_token(store, field, name)
    return Index(array, token)


def mark_indexes_stale(array: blosc2.NDArray) -> None:
    store = _load_store(array)
    if not store["indexes"]:
        return
    changed = False
    for descriptor in store["indexes"].values():
        if not descriptor.get("stale", False):
            descriptor["stale"] = True
            changed = True
    if changed:
        _save_store(array, store)
        _invalidate_query_cache(array)


def _descriptor_for(array: blosc2.NDArray, field: str | None) -> dict | None:
    return _descriptor_for_target(array, _field_target_descriptor(field))


def _descriptor_for_target(
    array: blosc2.NDArray, target: dict, array_to_col: dict | None = None
) -> dict | None:
    token = _target_token(target)
    if array_to_col is not None and token == SELF_TARGET_NAME:
        col_name = array_to_col.get(id(array))
        if col_name is not None:
            token = col_name
    descriptor = _load_store(array)["indexes"].get(token)
    if descriptor is None or descriptor.get("stale", False):
        return None
    # In compact stores every column shares one urlpath, so the urlpath-keyed
    # index store can hand back a *sibling* column's descriptor.  When an owner
    # was registered for this (array_key, token), require *array* to be it;
    # otherwise this predicate is on a different column that merely shares the
    # urlpath and (post column-alignment) the same shape/chunks.
    owner_id = _DESCRIPTOR_OWNERS.get((_array_key(array), token))
    if owner_id is not None and owner_id != id(array):
        return None
    if descriptor.get("version") != INDEX_FORMAT_VERSION:
        return None
    if descriptor.get("kind") == "bucket":
        bucket = descriptor.get("bucket", {})
        if bucket.get("layout") != "chunk-local-v1" or "values_path" not in bucket:
            return None
    if descriptor.get("kind") == "partial":
        partial = descriptor.get("partial", {})
        if partial.get("layout") != "chunk-local-v1" or "values_path" not in partial:
            return None
    if tuple(descriptor.get("shape", ())) != tuple(array.shape):
        return None
    if tuple(descriptor.get("chunks", ())) != tuple(array.chunks):
        return None
    return descriptor


def _open_level_summary_handle(array: blosc2.NDArray, descriptor: dict, level: str):
    level_info = descriptor["levels"][level]
    return _open_sidecar_handle(array, descriptor["token"], "summary_handle", level, level_info["path"])


def _candidate_units_from_summary_handle(summary_handle, op: str, value, dtype: np.dtype) -> np.ndarray:
    length = int(summary_handle.shape[0])
    if length == 0:
        return np.zeros(0, dtype=bool)
    chunk_len = int(summary_handle.chunks[0]) if hasattr(summary_handle, "chunks") else length
    candidate = np.empty(length, dtype=bool)
    for start in range(0, length, chunk_len):
        stop = min(start + chunk_len, length)
        candidate[start:stop] = _candidate_units_from_summary(
            _read_sidecar_span(summary_handle, start, stop), op, value, dtype
        )
    return candidate


def _candidate_units_from_exact_plan_handle(
    summary_handle, dtype: np.dtype, plan: ExactPredicatePlan
) -> np.ndarray:
    length = int(summary_handle.shape[0])
    candidate_units = np.ones(length, dtype=bool)
    if plan.lower is not None:
        lower_op = ">=" if plan.lower_inclusive else ">"
        candidate_units &= _candidate_units_from_summary_handle(summary_handle, lower_op, plan.lower, dtype)
    if plan.upper is not None:
        upper_op = "<=" if plan.upper_inclusive else "<"
        candidate_units &= _candidate_units_from_summary_handle(summary_handle, upper_op, plan.upper, dtype)
    return candidate_units


def _candidate_units_from_boundaries_handle(boundaries_handle, plan: ExactPredicatePlan) -> np.ndarray:
    length = int(boundaries_handle.shape[0])
    if length == 0:
        return np.zeros(0, dtype=bool)
    chunk_len = int(boundaries_handle.chunks[0]) if hasattr(boundaries_handle, "chunks") else length
    candidate = np.empty(length, dtype=bool)
    for start in range(0, length, chunk_len):
        stop = min(start + chunk_len, length)
        candidate[start:stop] = _candidate_units_from_boundaries(
            _read_sidecar_span(boundaries_handle, start, stop), plan
        )
    return candidate


def _read_offset_pair(offsets_handle, index: int) -> tuple[int, int]:
    pair = _read_sidecar_span(offsets_handle, index, index + 2)
    return int(pair[0]), int(pair[1])


def _full_merge_cache_key(array: blosc2.NDArray, token: str, name: str):
    return _data_cache_key(array, token, "full_merged", name)


def _clear_full_merge_cache(array: blosc2.NDArray, token: str) -> None:
    _DATA_CACHE.pop(_full_merge_cache_key(array, token, "values"), None)
    _DATA_CACHE.pop(_full_merge_cache_key(array, token, "positions"), None)


def _load_full_navigation_handles(array: blosc2.NDArray, descriptor: dict):
    full = descriptor.get("full")
    if full is None:
        raise RuntimeError("full-index metadata is not available")
    token = descriptor["token"]
    l1 = _open_sidecar_handle(array, token, "full_nav_handle", "l1", full.get("l1_path"))
    l2 = _open_sidecar_handle(array, token, "full_nav_handle", "l2", full.get("l2_path"))
    return l1, l2


def _load_full_sidecar_handles(array: blosc2.NDArray, descriptor: dict):
    full = descriptor.get("full")
    if full is None:
        raise RuntimeError("full-index metadata is not available")
    token = descriptor["token"]
    values_sidecar = _open_sidecar_handle(array, token, "full_handle", "values", full["values_path"])
    positions_sidecar = _open_sidecar_handle(
        array, token, "full_handle", "positions", full["positions_path"]
    )
    return values_sidecar, positions_sidecar


def _load_opsi_sidecar_handles(array: blosc2.NDArray, descriptor: dict):
    opsi = descriptor.get("opsi")
    if opsi is None:
        raise RuntimeError("OPSI-index metadata is not available")
    token = descriptor["token"]
    values_sidecar = _open_sidecar_handle(array, token, "opsi_handle", "values", opsi["values_path"])
    positions_sidecar = _open_sidecar_handle(
        array, token, "opsi_handle", "positions", opsi["positions_path"]
    )
    return values_sidecar, positions_sidecar


def _load_opsi_navigation_handles(array: blosc2.NDArray, descriptor: dict):
    opsi = descriptor.get("opsi")
    if opsi is None:
        raise RuntimeError("OPSI-index metadata is not available")
    token = descriptor["token"]
    mins_sidecar = _open_sidecar_handle(array, token, "opsi_nav_handle", "mins", opsi.get("mins_path"))
    maxs_sidecar = _open_sidecar_handle(array, token, "opsi_nav_handle", "maxs", opsi.get("maxs_path"))
    return mins_sidecar, maxs_sidecar


def _load_full_run_sidecar_handles(array: blosc2.NDArray, descriptor: dict, run: dict):
    run_id = int(run["id"])
    token = descriptor["token"]
    values_sidecar = _open_sidecar_handle(
        array, token, "full_run_handle", f"{run_id}.values", run["values_path"]
    )
    positions_sidecar = _open_sidecar_handle(
        array, token, "full_run_handle", f"{run_id}.positions", run["positions_path"]
    )
    return values_sidecar, positions_sidecar


def _load_partial_l1_handle(array: blosc2.NDArray, descriptor: dict):
    partial = descriptor.get("partial")
    if partial is None:
        raise RuntimeError("partial-index metadata is not available")
    token = descriptor["token"]
    return _open_sidecar_handle(array, token, "partial_nav_handle", "l1", partial["l1_path"])


def _load_partial_sidecar_handles(array: blosc2.NDArray, descriptor: dict):
    partial = descriptor.get("partial")
    if partial is None:
        raise RuntimeError("partial-index metadata is not available")
    token = descriptor["token"]
    values_sidecar = _open_sidecar_handle(array, token, "partial_handle", "values", partial["values_path"])
    positions_sidecar = _open_sidecar_handle(
        array, token, "partial_handle", "positions", partial["positions_path"]
    )
    l2_sidecar = _open_sidecar_handle(array, token, "partial_nav_handle", "l2", partial["l2_path"])
    return values_sidecar, positions_sidecar, l2_sidecar


def _load_partial_offsets_handle(array: blosc2.NDArray, descriptor: dict):
    partial = descriptor.get("partial")
    if partial is None:
        raise RuntimeError("partial-index metadata is not available")
    token = descriptor["token"]
    return _open_sidecar_handle(array, token, "partial_handle", "offsets", partial["offsets_path"])


def _load_bucket_l1_handle(array: blosc2.NDArray, descriptor: dict):
    bucket = descriptor.get("bucket")
    if bucket is None:
        raise RuntimeError("bucket index metadata is not available")
    token = descriptor["token"]
    return _open_sidecar_handle(array, token, "bucket_nav_handle", "l1", bucket["l1_path"])


def _load_bucket_sidecar_handles(array: blosc2.NDArray, descriptor: dict):
    bucket = descriptor.get("bucket")
    if bucket is None:
        raise RuntimeError("bucket index metadata is not available")
    token = descriptor["token"]
    values_sidecar = _open_sidecar_handle(array, token, "bucket_handle", "values", bucket["values_path"])
    bucket_sidecar = _open_sidecar_handle(
        array, token, "bucket_handle", "bucket_positions", bucket["bucket_positions_path"]
    )
    l2_sidecar = _open_sidecar_handle(array, token, "bucket_nav_handle", "l2", bucket["l2_path"])
    return values_sidecar, bucket_sidecar, l2_sidecar


def _load_bucket_offsets_handle(array: blosc2.NDArray, descriptor: dict):
    bucket = descriptor.get("bucket")
    if bucket is None:
        raise RuntimeError("bucket index metadata is not available")
    token = descriptor["token"]
    return _open_sidecar_handle(array, token, "bucket_handle", "offsets", bucket["offsets_path"])


def _normalize_scalar(value, dtype: np.dtype):
    if isinstance(value, np.generic):
        return value.item()
    if dtype.kind == "f" and isinstance(value, float) and np.isnan(value):
        raise ValueError("NaN comparisons are not indexable")
    return np.asarray(value, dtype=dtype)[()]


def _candidate_units_from_summary(summaries: np.ndarray, op: str, value, dtype: np.dtype) -> np.ndarray:
    mins = summaries["min"]
    maxs = summaries["max"]
    flags = summaries["flags"]
    valid = (flags & FLAG_ALL_NAN) == 0
    value = _normalize_scalar(value, dtype)
    if op == "==":
        return valid & (mins <= value) & (value <= maxs)
    if op == "<":
        return valid & (mins < value)
    if op == "<=":
        return valid & (mins <= value)
    if op == ">":
        return valid & (maxs > value)
    if op == ">=":
        return valid & (maxs >= value)
    raise ValueError(f"unsupported comparison operator {op!r}")


def _intervals_from_sorted(values: np.ndarray, op: str, value, dtype: np.dtype) -> list[tuple[int, int]]:
    value = _normalize_scalar(value, dtype)
    if op == "==":
        lo = np.searchsorted(values, value, side="left")
        hi = np.searchsorted(values, value, side="right")
    elif op == "<":
        lo = 0
        hi = np.searchsorted(values, value, side="left")
    elif op == "<=":
        lo = 0
        hi = np.searchsorted(values, value, side="right")
    elif op == ">":
        lo = np.searchsorted(values, value, side="right")
        hi = len(values)
    elif op == ">=":
        lo = np.searchsorted(values, value, side="left")
        hi = len(values)
    else:
        raise ValueError(f"unsupported comparison operator {op!r}")
    return [] if lo >= hi else [(int(lo), int(hi))]


def _operand_target(operand) -> tuple[blosc2.NDArray, str | None] | None:
    if isinstance(operand, blosc2.NDField):
        return operand.ndarr, operand.field
    if isinstance(operand, blosc2.NDArray):
        return operand, None
    return None


def _literal_value(node: ast.AST):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        value = _literal_value(node.operand)
        if isinstance(value, bool):
            raise ValueError("boolean negation is not a scalar literal here")
        return -value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
        return _literal_value(node.operand)
    raise ValueError("node is not a supported scalar literal")


def _flip_operator(op: str) -> str:
    return {"<": ">", "<=": ">=", ">": "<", ">=": "<=", "==": "=="}[op]


def _compare_operator(node: ast.AST) -> str | None:
    if isinstance(node, ast.Eq):
        return "=="
    if isinstance(node, ast.Lt):
        return "<"
    if isinstance(node, ast.LtE):
        return "<="
    if isinstance(node, ast.Gt):
        return ">"
    if isinstance(node, ast.GtE):
        return ">="
    return None


def _compare_target_from_node(node: ast.AST, operands: dict) -> tuple[blosc2.NDArray, dict] | None:
    if isinstance(node, ast.Name):
        operand = operands.get(node.id)
        target = _operand_target(operand) if operand is not None else None
        if target is None:
            return None
        base, field = target
        if base.ndim != 1:
            return None
        return base, _field_target_descriptor(field)

    normalized = _normalize_expression_node(node, operands)
    if normalized is None:
        return None
    base, expression_key, dependencies = normalized
    return base, _expression_target_descriptor(ast.unparse(node), expression_key, dependencies)


def _target_from_compare(
    node: ast.Compare, operands: dict
) -> tuple[blosc2.NDArray, dict, str, object] | None:
    if len(node.ops) != 1 or len(node.comparators) != 1:
        return None
    op = _compare_operator(node.ops[0])
    if op is None:
        return None

    try:
        left_target = _compare_target_from_node(node.left, operands)
        right_target = _compare_target_from_node(node.comparators[0], operands)
        if left_target is not None:
            value = _literal_value(node.comparators[0])
        elif right_target is not None:
            value = _literal_value(node.left)
            op = _flip_operator(op)
        else:
            return None
    except ValueError:
        return None

    base, target = left_target if left_target is not None else right_target
    return base, target, op, value


def _finest_level(descriptor: dict) -> str:
    level_names = tuple(descriptor["levels"])
    return level_names[-1]


def _plan_segment_compare(
    node: ast.Compare, operands: dict, array_to_col: dict | None = None
) -> SegmentPredicatePlan | None:
    target = _target_from_compare(node, operands)
    if target is None:
        return None
    base, target_info, op, value = target
    descriptor = _descriptor_for_target(base, target_info, array_to_col)
    if descriptor is None:
        return None
    level = _finest_level(descriptor)
    level_info = descriptor["levels"][level]
    dtype = np.dtype(descriptor["dtype"])
    try:
        summaries = _open_level_summary_handle(base, descriptor, level)
        candidate_units = _candidate_units_from_summary_handle(summaries, op, value, dtype)
    except (RuntimeError, ValueError, TypeError):
        return None
    return SegmentPredicatePlan(
        base=base,
        candidate_units=candidate_units,
        descriptor=descriptor,
        target=target_info,
        field=_target_field(target_info),
        level=level,
        segment_len=level_info["segment_len"],
    )


def _grid_compatible_segment_plans(left: SegmentPredicatePlan, right: SegmentPredicatePlan) -> bool:
    """Whether two segment plans share an identical row→segment grid.

    Segments are row-aligned (segment ``i`` always covers rows
    ``[i*segment_len, (i+1)*segment_len)``), so two plans with the same level,
    ``segment_len``, candidate-mask shape, and row count map every segment index
    to the same rows — even when they belong to *different* columns.  This is
    what lets cross-column AND/OR intersect their per-segment candidate masks.
    Compact-store columns that join the auto-aligned fast-eval grid share
    ``segment_len``; columns excluded from alignment may not, and fall back.
    """
    return (
        left.level == right.level
        and left.segment_len == right.segment_len
        and left.candidate_units.shape == right.candidate_units.shape
        and int(left.base.shape[0]) == int(right.base.shape[0])
    )


def _segment_plan_scan_rows(plan: SegmentPredicatePlan) -> int:
    """Estimated rows the downstream scan must visit for *plan* (selectivity)."""
    return int(np.count_nonzero(plan.candidate_units)) * int(plan.segment_len)


def _merge_segment_plans(
    left: SegmentPredicatePlan, right: SegmentPredicatePlan, op: str
) -> SegmentPredicatePlan | None:
    if not _grid_compatible_segment_plans(left, right):
        # Different columns can carry indexes on incompatible grids (e.g. a
        # column excluded from the aligned fast-eval set).  We cannot combine
        # their per-segment masks directly.
        if op == "and":
            # AND pruning by a *superset* of the true candidates is always
            # correct, because the downstream scan re-evaluates the full
            # predicate per row.  Keep whichever single side prunes more so we
            # never do worse than single-column pruning.
            return left if _segment_plan_scan_rows(left) <= _segment_plan_scan_rows(right) else right
        # OR cannot be pruned to one side's segments — the other column may
        # match rows in segments this side discarded.  Fall back to full scan.
        return None
    if op == "and":
        candidate_units = left.candidate_units & right.candidate_units
    else:
        candidate_units = left.candidate_units | right.candidate_units
    return SegmentPredicatePlan(
        base=left.base,
        candidate_units=candidate_units,
        descriptor=left.descriptor,
        target=left.target,
        field=left.field,
        level=left.level,
        segment_len=left.segment_len,
    )


def _plan_segment_boolop(
    node: ast.BoolOp, operands: dict, array_to_col: dict | None = None
) -> SegmentPredicatePlan | None:
    op = "and" if isinstance(node.op, ast.And) else "or" if isinstance(node.op, ast.Or) else None
    if op is None:
        return None
    plans = [_plan_segment_node(value, operands, array_to_col) for value in node.values]
    if op == "and":
        plans = [plan for plan in plans if plan is not None]
        if not plans:
            return None
    elif any(plan is None for plan in plans):
        return None

    plan = plans[0]
    for other in plans[1:]:
        merged = _merge_segment_plans(plan, other, op)
        if merged is None:
            return None
        plan = merged
    return plan


def _plan_segment_bitop(
    node: ast.BinOp, operands: dict, array_to_col: dict | None = None
) -> SegmentPredicatePlan | None:
    if isinstance(node.op, ast.BitAnd):
        op = "and"
    elif isinstance(node.op, ast.BitOr):
        op = "or"
    else:
        return None

    left = _plan_segment_node(node.left, operands, array_to_col)
    right = _plan_segment_node(node.right, operands, array_to_col)
    if op == "and":
        if left is None:
            return right
        if right is None:
            return left
        return _merge_segment_plans(left, right, op)
    if left is None or right is None:
        return None
    return _merge_segment_plans(left, right, op)


def _plan_segment_node(
    node: ast.AST, operands: dict, array_to_col: dict | None = None
) -> SegmentPredicatePlan | None:
    if isinstance(node, ast.Compare):
        return _plan_segment_compare(node, operands, array_to_col)
    if isinstance(node, ast.BoolOp):
        return _plan_segment_boolop(node, operands, array_to_col)
    if isinstance(node, ast.BinOp):
        return _plan_segment_bitop(node, operands, array_to_col)
    return None


def _plan_exact_compare(
    node: ast.Compare, operands: dict, array_to_col: dict | None = None
) -> ExactPredicatePlan | None:
    target = _target_from_compare(node, operands)
    if target is None:
        return None
    base, target_info, op, value = target
    descriptor = _descriptor_for_target(base, target_info, array_to_col)
    if descriptor is None or descriptor.get("kind") not in {"bucket", "partial", "full", "opsi"}:
        return None
    try:
        value = _normalize_scalar(value, np.dtype(descriptor["dtype"]))
    except (RuntimeError, ValueError, TypeError):
        return None
    if op == "==":
        return ExactPredicatePlan(
            base=base,
            descriptor=descriptor,
            target=target_info,
            field=_target_field(target_info),
            lower=value,
            lower_inclusive=True,
            upper=value,
            upper_inclusive=True,
        )
    if op == ">":
        return ExactPredicatePlan(
            base=base,
            descriptor=descriptor,
            target=target_info,
            field=_target_field(target_info),
            lower=value,
            lower_inclusive=False,
        )
    if op == ">=":
        return ExactPredicatePlan(
            base=base,
            descriptor=descriptor,
            target=target_info,
            field=_target_field(target_info),
            lower=value,
            lower_inclusive=True,
        )
    if op == "<":
        return ExactPredicatePlan(
            base=base,
            descriptor=descriptor,
            target=target_info,
            field=_target_field(target_info),
            upper=value,
            upper_inclusive=False,
        )
    if op == "<=":
        return ExactPredicatePlan(
            base=base,
            descriptor=descriptor,
            target=target_info,
            field=_target_field(target_info),
            upper=value,
            upper_inclusive=True,
        )
    return None


def _same_base(left: ExactPredicatePlan, right: ExactPredicatePlan) -> bool:
    return left.base is right.base and left.descriptor["token"] == right.descriptor["token"]


def _merge_lower_bound(
    left: object | None, left_inclusive: bool, right: object | None, right_inclusive: bool
) -> tuple[object | None, bool]:
    if left is None:
        return right, right_inclusive
    if right is None:
        return left, left_inclusive
    if left < right:
        return right, right_inclusive
    if left > right:
        return left, left_inclusive
    return left, left_inclusive and right_inclusive


def _merge_upper_bound(
    left: object | None, left_inclusive: bool, right: object | None, right_inclusive: bool
) -> tuple[object | None, bool]:
    if left is None:
        return right, right_inclusive
    if right is None:
        return left, left_inclusive
    if left < right:
        return left, left_inclusive
    if left > right:
        return right, right_inclusive
    return left, left_inclusive and right_inclusive


def _merge_exact_plans(
    left: ExactPredicatePlan, right: ExactPredicatePlan, op: str
) -> ExactPredicatePlan | None:
    if op != "and" or not _same_base(left, right):
        return None
    lower, lower_inclusive = _merge_lower_bound(
        left.lower, left.lower_inclusive, right.lower, right.lower_inclusive
    )
    upper, upper_inclusive = _merge_upper_bound(
        left.upper, left.upper_inclusive, right.upper, right.upper_inclusive
    )
    return ExactPredicatePlan(
        base=left.base,
        descriptor=left.descriptor,
        target=left.target,
        field=left.field,
        lower=lower,
        lower_inclusive=lower_inclusive,
        upper=upper,
        upper_inclusive=upper_inclusive,
    )


def _plan_exact_conjunction(
    node: ast.AST, operands: dict, array_to_col: dict | None = None
) -> list[ExactPredicatePlan] | None:
    if isinstance(node, ast.Compare):
        plan = _plan_exact_compare(node, operands, array_to_col)
        return None if plan is None else [plan]
    if isinstance(node, ast.BoolOp):
        if not isinstance(node.op, ast.And):
            return None
        plans = []
        for value in node.values:
            subplans = _plan_exact_conjunction(value, operands, array_to_col)
            if subplans is None:
                return None
            plans.extend(subplans)
        return plans
    if isinstance(node, ast.BinOp):
        if not isinstance(node.op, ast.BitAnd):
            return None
        left = _plan_exact_conjunction(node.left, operands, array_to_col)
        right = _plan_exact_conjunction(node.right, operands, array_to_col)
        if left is None and right is None:
            return None
        if left is None:
            return right
        if right is None:
            return left
        return left + right
    return None


def _plan_exact_boolop(
    node: ast.BoolOp, operands: dict, array_to_col: dict | None = None
) -> ExactPredicatePlan | None:
    if not isinstance(node.op, ast.And):
        return None
    plans = [_plan_exact_node(value, operands, array_to_col) for value in node.values]
    if any(plan is None for plan in plans):
        return None
    plan = plans[0]
    for other in plans[1:]:
        merged = _merge_exact_plans(plan, other, "and")
        if merged is None:
            return None
        plan = merged
    return plan


def _plan_exact_bitop(
    node: ast.BinOp, operands: dict, array_to_col: dict | None = None
) -> ExactPredicatePlan | None:
    if not isinstance(node.op, ast.BitAnd):
        return None
    left = _plan_exact_node(node.left, operands, array_to_col)
    right = _plan_exact_node(node.right, operands, array_to_col)
    if left is None or right is None:
        return None
    return _merge_exact_plans(left, right, "and")


def _plan_exact_node(
    node: ast.AST, operands: dict, array_to_col: dict | None = None
) -> ExactPredicatePlan | None:
    if isinstance(node, ast.Compare):
        return _plan_exact_compare(node, operands, array_to_col)
    if isinstance(node, ast.BoolOp):
        return _plan_exact_boolop(node, operands, array_to_col)
    if isinstance(node, ast.BinOp):
        return _plan_exact_bitop(node, operands, array_to_col)
    return None


def _range_is_empty(plan: ExactPredicatePlan) -> bool:
    if plan.lower is None or plan.upper is None:
        return False
    if plan.lower < plan.upper:
        return False
    if plan.lower > plan.upper:
        return True
    return not (plan.lower_inclusive and plan.upper_inclusive)


def _candidate_units_from_exact_plan(
    summaries: np.ndarray, dtype: np.dtype, plan: ExactPredicatePlan
) -> np.ndarray:
    candidate_units = np.ones(len(summaries), dtype=bool)
    if plan.lower is not None:
        lower_op = ">=" if plan.lower_inclusive else ">"
        candidate_units &= _candidate_units_from_summary(summaries, lower_op, plan.lower, dtype)
    if plan.upper is not None:
        upper_op = "<=" if plan.upper_inclusive else "<"
        candidate_units &= _candidate_units_from_summary(summaries, upper_op, plan.upper, dtype)
    return candidate_units


def _search_bounds(values: np.ndarray, plan: ExactPredicatePlan) -> tuple[int, int]:
    try:
        return indexing_ext.index_search_bounds(
            values, plan.lower, plan.lower_inclusive, plan.upper, plan.upper_inclusive
        )
    except TypeError:
        lo = 0
        hi = len(values)
        if plan.lower is not None:
            side = "left" if plan.lower_inclusive else "right"
            lo = int(np.searchsorted(values, plan.lower, side=side))
        if plan.upper is not None:
            side = "right" if plan.upper_inclusive else "left"
            hi = int(np.searchsorted(values, plan.upper, side=side))
        return lo, hi


def _candidate_units_from_boundaries(boundaries: np.ndarray, plan: ExactPredicatePlan) -> np.ndarray:
    if len(boundaries) == 0:
        return np.zeros(0, dtype=bool)
    starts = boundaries["start"]
    ends = boundaries["end"]
    if ends.dtype.kind == "f":
        nan_mask = np.isnan(ends)
        if np.any(nan_mask):
            ends = ends.copy()
            ends[nan_mask] = np.inf
        nan_mask_start = np.isnan(starts)
        if np.any(nan_mask_start):
            starts = starts.copy()
            starts[nan_mask_start] = -np.inf
    candidate = np.ones(len(boundaries), dtype=bool)
    if plan.lower is not None:
        candidate &= ends >= plan.lower if plan.lower_inclusive else ends > plan.lower
    if plan.upper is not None:
        candidate &= starts <= plan.upper if plan.upper_inclusive else starts < plan.upper
    return candidate


def _full_runs_need_bounded_fallback(descriptor: dict) -> bool:
    full = descriptor.get("full")
    if full is None:
        return False
    runs = tuple(full.get("runs", ()))
    if not runs:
        return False
    if len(runs) >= FULL_RUN_BOUNDED_FALLBACK_RUNS:
        return True
    return sum(int(run["length"]) for run in runs) >= FULL_RUN_BOUNDED_FALLBACK_ITEMS


def _full_query_mode_override() -> str:
    mode = os.getenv("BLOSC2_FULL_EXACT_QUERY_MODE", "auto").strip().lower()
    if mode not in {"auto", "selective-ooc", "whole-load"}:
        return "auto"
    return mode


def _contiguous_true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    true_ids = np.flatnonzero(mask)
    if len(true_ids) == 0:
        return []
    breaks = np.nonzero(np.diff(true_ids) != 1)[0] + 1
    runs = []
    start = 0
    for stop in (*breaks, len(true_ids)):
        part = true_ids[start:stop]
        runs.append((int(part[0]), int(part[-1]) + 1))
        start = stop
    return runs


def _sorted_chunk_boundaries_from_handle(
    array: blosc2.NDArray,
    token: str,
    category: str,
    name: str,
    values_sidecar,
    dtype: np.dtype,
) -> np.ndarray:
    cache_key = _data_cache_key(array, token, category, name)
    cached = _DATA_CACHE.get(cache_key)
    if cached is not None:
        return cached

    size = int(values_sidecar.shape[0])
    chunk_len = int(values_sidecar.chunks[0])
    nchunks = math.ceil(size / chunk_len)
    boundaries = np.empty(nchunks, dtype=_boundary_dtype(dtype))
    start_value = np.empty(1, dtype=dtype)
    end_value = np.empty(1, dtype=dtype)
    for chunk_id in range(nchunks):
        chunk_start = chunk_id * chunk_len
        chunk_stop = min(chunk_start + chunk_len, size)
        values_sidecar.get_1d_span_numpy(start_value, chunk_id, 0, 1)
        values_sidecar.get_1d_span_numpy(end_value, chunk_id, chunk_stop - chunk_start - 1, 1)
        boundaries[chunk_id] = (start_value[0], end_value[0])
    _DATA_CACHE[cache_key] = boundaries
    return boundaries


def _full_supports_selective_ooc_lookup(array: blosc2.NDArray, descriptor: dict) -> bool:
    full = descriptor.get("full")
    if full is None or full.get("runs"):
        return False
    try:
        values_sidecar, positions_sidecar = _load_full_sidecar_handles(array, descriptor)
        l1_sidecar, l2_sidecar = _load_full_navigation_handles(array, descriptor)
    except Exception:
        return False
    if int(values_sidecar.chunks[0]) != int(full.get("sidecar_chunk_len", values_sidecar.chunks[0])):
        return False
    return (
        _supports_block_reads(values_sidecar)
        and _supports_block_reads(positions_sidecar)
        and _supports_block_reads(l1_sidecar)
        and _supports_block_reads(l2_sidecar)
    )


def _exact_positions_from_sorted_chunks(
    values_sidecar,
    positions_sidecar,
    boundaries: np.ndarray,
    plan: ExactPredicatePlan,
    chunk_len: int,
    dtype: np.dtype,
) -> np.ndarray:
    candidate_chunks = _candidate_units_from_boundaries(boundaries, plan)
    if not np.any(candidate_chunks):
        return np.empty(0, dtype=np.int64)

    parts = []
    size = int(values_sidecar.shape[0])
    for chunk_id in np.flatnonzero(candidate_chunks):
        chunk_start = int(chunk_id) * chunk_len
        chunk_stop = min(chunk_start + chunk_len, size)
        span_items = chunk_stop - chunk_start
        span_values = np.empty(span_items, dtype=dtype)
        values_sidecar.get_1d_span_numpy(span_values, int(chunk_id), 0, span_items)
        lo, hi = _search_bounds(span_values, plan)
        if lo >= hi:
            continue
        matched = np.empty(hi - lo, dtype=np.int64)
        _read_ndarray_linear_span(positions_sidecar, chunk_start + lo, matched)
        parts.append(matched)

    if not parts:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(parts) if len(parts) > 1 else parts[0]


def _exact_positions_from_compact_full_base(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> np.ndarray:
    full = descriptor["full"]
    l1_sidecar, l2_sidecar = _load_full_navigation_handles(array, descriptor)
    candidate_chunks = _candidate_units_from_boundaries_handle(l1_sidecar, plan)
    if not np.any(candidate_chunks):
        return np.empty(0, dtype=np.int64)

    candidate_blocks = _candidate_units_from_boundaries_handle(l2_sidecar, plan)
    if not np.any(candidate_blocks):
        return np.empty(0, dtype=np.int64)

    values_sidecar, positions_sidecar = _load_full_sidecar_handles(array, descriptor)
    dtype = np.dtype(descriptor["dtype"])
    chunk_len = int(full["sidecar_chunk_len"])
    block_len = int(full["sidecar_block_len"])
    size = int(values_sidecar.shape[0])
    parts = []
    span_count = 0

    for chunk_id in np.flatnonzero(candidate_chunks):
        chunk_start = int(chunk_id) * chunk_len
        chunk_stop = min(chunk_start + chunk_len, size)
        first_block = chunk_start // block_len
        nblocks = math.ceil((chunk_stop - chunk_start) / block_len)
        block_mask = np.asarray(candidate_blocks[first_block : first_block + nblocks], dtype=bool)
        if not np.any(block_mask):
            continue
        span_runs = _contiguous_true_runs(block_mask)
        span_count += len(span_runs)
        if span_count > FULL_SELECTIVE_OOC_MAX_SPANS:
            raise RuntimeError("too many candidate spans for selective full lookup")

        for block_start_idx, block_stop_idx in span_runs:
            span_start = chunk_start + block_start_idx * block_len
            span_stop = min(chunk_start + block_stop_idx * block_len, chunk_stop)
            local_start = span_start - chunk_start
            span_items = span_stop - span_start
            span_values = np.empty(span_items, dtype=dtype)
            values_sidecar.get_1d_span_numpy(span_values, int(chunk_id), local_start, span_items)
            lo, hi = _search_bounds(span_values, plan)
            if lo >= hi:
                continue
            matched = np.empty(hi - lo, dtype=np.int64)
            _read_ndarray_linear_span(positions_sidecar, span_start + lo, matched)
            parts.append(matched)

    if not parts:
        return np.empty(0, dtype=np.int64)
    positions = np.concatenate(parts) if len(parts) > 1 else parts[0]
    return np.sort(positions.astype(np.int64, copy=False), kind="stable")


def _exact_positions_from_full_runs_bounded(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> np.ndarray:
    full = descriptor["full"]
    dtype = np.dtype(descriptor["dtype"])
    parts = []

    base_descriptor = descriptor.copy()
    base_full = full.copy()
    base_full["runs"] = []
    base_descriptor["full"] = base_full
    if _full_supports_selective_ooc_lookup(array, base_descriptor):
        base_positions = _exact_positions_from_compact_full_base(array, base_descriptor, plan)
        if len(base_positions):
            parts.append(base_positions)
    else:
        base_values_sidecar, base_positions_sidecar = _load_full_sidecar_handles(array, base_descriptor)
        base_chunk_boundaries = _sorted_chunk_boundaries_from_handle(
            array,
            descriptor["token"],
            "full_bounds",
            "chunks",
            base_values_sidecar,
            dtype,
        )
        base_positions = _exact_positions_from_sorted_chunks(
            base_values_sidecar,
            base_positions_sidecar,
            base_chunk_boundaries,
            plan,
            int(base_values_sidecar.chunks[0]),
            dtype,
        )
        if len(base_positions):
            parts.append(np.sort(base_positions.astype(np.int64, copy=False), kind="stable"))

    for run in full.get("runs", ()):
        run_values_sidecar, run_positions_sidecar = _load_full_run_sidecar_handles(array, descriptor, run)
        chunk_boundaries = _sorted_chunk_boundaries_from_handle(
            array,
            descriptor["token"],
            "full_run_bounds",
            f"{int(run['id'])}.chunks",
            run_values_sidecar,
            dtype,
        )
        run_positions = _exact_positions_from_sorted_chunks(
            run_values_sidecar,
            run_positions_sidecar,
            chunk_boundaries,
            plan,
            int(run_values_sidecar.chunks[0]),
            dtype,
        )
        if len(run_positions):
            parts.append(run_positions)

    if not parts:
        return np.empty(0, dtype=np.int64)
    positions = np.concatenate(parts) if len(parts) > 1 else parts[0]
    return np.sort(positions.astype(np.int64, copy=False), kind="stable")


def _exact_positions_from_full_selective_ooc(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> np.ndarray:
    return _exact_positions_from_compact_full_base(array, descriptor, plan)


def _opsi_block_boundaries_from_handles(
    array: blosc2.NDArray, descriptor: dict, mins_sidecar, maxs_sidecar, dtype: np.dtype
) -> np.ndarray:
    cache_key = _data_cache_key(array, descriptor["token"], "opsi_bounds", "blocks")
    cached = _DATA_CACHE.get(cache_key)
    if cached is not None:
        return cached
    nblocks = int(mins_sidecar.shape[0])
    boundaries = np.empty(nblocks, dtype=_boundary_dtype(dtype))
    if nblocks:
        boundaries["start"] = _read_sidecar_span(mins_sidecar, 0, nblocks)
        boundaries["end"] = _read_sidecar_span(maxs_sidecar, 0, nblocks)
    _DATA_CACHE[cache_key] = boundaries
    return boundaries


def _exact_positions_from_opsi_block_nav(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> np.ndarray:
    dtype = np.dtype(descriptor["dtype"])
    opsi = descriptor["opsi"]
    values_sidecar, positions_sidecar = _load_opsi_sidecar_handles(array, descriptor)
    mins_sidecar, maxs_sidecar = _load_opsi_navigation_handles(array, descriptor)
    boundaries = _opsi_block_boundaries_from_handles(array, descriptor, mins_sidecar, maxs_sidecar, dtype)
    candidate_blocks = _candidate_units_from_boundaries(boundaries, plan)
    if not np.any(candidate_blocks):
        return np.empty(0, dtype=np.int64)

    chunk_len = int(opsi.get("chunk_len", values_sidecar.chunks[0]))
    block_len = int(opsi.get("block_len", values_sidecar.blocks[0]))
    blocks_per_chunk = max(1, chunk_len // block_len)
    size = int(values_sidecar.shape[0])
    parts = []
    span_count = 0

    first_chunk = int(np.flatnonzero(candidate_blocks)[0]) // blocks_per_chunk
    last_chunk = int(np.flatnonzero(candidate_blocks)[-1]) // blocks_per_chunk
    for chunk_id in range(first_chunk, last_chunk + 1):
        first_block = chunk_id * blocks_per_chunk
        last_block = min(first_block + blocks_per_chunk, len(candidate_blocks))
        block_mask = np.asarray(candidate_blocks[first_block:last_block], dtype=bool)
        if not np.any(block_mask):
            continue
        for block_start_idx, block_stop_idx in _contiguous_true_runs(block_mask):
            span_count += 1
            span_start = (first_block + block_start_idx) * block_len
            span_stop = min((first_block + block_stop_idx) * block_len, size)
            if span_start >= span_stop:
                continue
            local_start = span_start - chunk_id * chunk_len
            span_items = span_stop - span_start
            span_values = np.empty(span_items, dtype=dtype)
            values_sidecar.get_1d_span_numpy(span_values, chunk_id, local_start, span_items)
            lo, hi = _search_bounds(span_values, plan)
            if lo >= hi:
                continue
            matched = np.empty(hi - lo, dtype=np.int64)
            _read_ndarray_linear_span(positions_sidecar, span_start + lo, matched)
            parts.append(matched)

    if not parts:
        return np.empty(0, dtype=np.int64)
    positions = np.concatenate(parts) if len(parts) > 1 else parts[0]
    return np.sort(positions.astype(np.int64, copy=False), kind="stable")


def _exact_positions_from_opsi(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> np.ndarray:
    if _range_is_empty(plan):
        return np.empty(0, dtype=np.int64)
    try:
        return _exact_positions_from_opsi_block_nav(array, descriptor, plan)
    except Exception:
        pass
    dtype = np.dtype(descriptor["dtype"])
    values_sidecar, positions_sidecar = _load_opsi_sidecar_handles(array, descriptor)
    chunk_boundaries = _sorted_chunk_boundaries_from_handle(
        array,
        descriptor["token"],
        "opsi_bounds",
        "chunks",
        values_sidecar,
        dtype,
    )
    positions = _exact_positions_from_sorted_chunks(
        values_sidecar,
        positions_sidecar,
        chunk_boundaries,
        plan,
        int(descriptor["opsi"].get("chunk_len", values_sidecar.chunks[0])),
        dtype,
    )
    if len(positions) == 0:
        return np.empty(0, dtype=np.int64)
    return np.sort(positions.astype(np.int64, copy=False), kind="stable")


def _exact_positions_from_full(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> np.ndarray:
    if _range_is_empty(plan):
        return np.empty(0, dtype=np.int64)
    if _full_run_count(descriptor):
        return _exact_positions_from_full_runs_bounded(array, descriptor, plan)
    if _full_supports_selective_ooc_lookup(array, descriptor):
        try:
            return _exact_positions_from_full_selective_ooc(array, descriptor, plan)
        except RuntimeError:
            pass
    dtype = np.dtype(descriptor["dtype"])
    values_sidecar, positions_sidecar = _load_full_sidecar_handles(array, descriptor)
    chunk_boundaries = _sorted_chunk_boundaries_from_handle(
        array,
        descriptor["token"],
        "full_bounds",
        "chunks",
        values_sidecar,
        dtype,
    )
    positions = _exact_positions_from_sorted_chunks(
        values_sidecar,
        positions_sidecar,
        chunk_boundaries,
        plan,
        int(values_sidecar.chunks[0]),
        dtype,
    )
    if len(positions) == 0:
        return np.empty(0, dtype=np.int64)
    return np.sort(positions.astype(np.int64, copy=False), kind="stable")


def _chunk_nav_supports_selective_ooc_lookup(array: blosc2.NDArray, descriptor: dict, kind: str) -> bool:
    if descriptor.get("kind") != kind or not descriptor.get("persistent", False):
        return False
    meta = descriptor.get("bucket" if kind == "bucket" else "partial")
    if meta is None or meta.get("layout") != "chunk-local-v1":
        return False
    required_paths = ("values_path", "l1_path", "l2_path")
    if any(meta.get(name) is None for name in required_paths):
        return False
    if kind == "bucket":
        if meta.get("bucket_positions_path") is None:
            return False
        try:
            values_sidecar, bucket_sidecar, l2_sidecar = _load_bucket_sidecar_handles(array, descriptor)
        except Exception:
            return False
        return (
            _supports_block_reads(array)
            and _supports_block_reads(values_sidecar)
            and _supports_block_reads(bucket_sidecar)
            and _supports_block_reads(l2_sidecar)
        )
    if meta.get("positions_path") is None:
        return False
    try:
        values_sidecar, positions_sidecar, l2_sidecar = _load_partial_sidecar_handles(array, descriptor)
    except Exception:
        return False
    return (
        _supports_block_reads(array)
        and _supports_block_reads(values_sidecar)
        and _supports_block_reads(positions_sidecar)
        and _supports_block_reads(l2_sidecar)
    )


def _chunk_nav_candidate_runs(
    l2_row: np.ndarray, segment_count: int, plan: ExactPredicatePlan
) -> tuple[list[tuple[int, int]], int]:
    segment_lo, segment_hi = _sorted_boundary_search_bounds(l2_row[:segment_count], plan)
    if segment_lo >= segment_hi:
        return [], 0
    return [(segment_lo, segment_hi)], segment_hi - segment_lo


def _index_query_thread_count(task_count: int) -> int:
    if blosc2.IS_WASM:
        return 1
    if task_count < INDEX_QUERY_MIN_CHUNKS_PER_THREAD:
        return 1
    configured_threads = int(getattr(blosc2, "nthreads", 1) or 1)
    return _python_executor_threads(min(configured_threads, task_count // INDEX_QUERY_MIN_CHUNKS_PER_THREAD))


def _chunk_batches(chunk_ids: np.ndarray, thread_count: int) -> list[np.ndarray]:
    if thread_count <= 1 or len(chunk_ids) == 0:
        return [chunk_ids]
    batch_size = max(1, math.ceil(len(chunk_ids) / thread_count))
    return [chunk_ids[start : start + batch_size] for start in range(0, len(chunk_ids), batch_size)]


def _downstream_query_thread_count(task_count: int, plan: IndexPlan) -> int:
    if plan.lookup_path == "chunk-nav-ooc":
        return 1
    return _index_query_thread_count(task_count)


def _merge_position_batches(position_batches: list[np.ndarray]) -> np.ndarray:
    if not position_batches:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(position_batches) if len(position_batches) > 1 else position_batches[0]


def _run_position_batches(chunk_ids: np.ndarray, thread_count: int, process_batch) -> tuple[np.ndarray, int]:
    if thread_count <= 1:
        return process_batch(chunk_ids)
    batches = _chunk_batches(chunk_ids, thread_count)
    position_batches = []
    total_candidate_segments = 0
    with ThreadPoolExecutor(max_workers=thread_count) as executor:
        for positions_part, batch_candidate_segments in executor.map(process_batch, batches):
            total_candidate_segments += batch_candidate_segments
            if len(positions_part) > 0:
                position_batches.append(positions_part)
    return _merge_position_batches(position_batches), total_candidate_segments


def _bucket_batch_result_dtype(where_x) -> np.dtype:
    return _where_output_dtype(where_x)


def _bucket_worker_source(where_x):
    if _supports_block_reads(where_x) and getattr(where_x, "urlpath", None) is not None:
        urlpath = str(where_x.urlpath)
        # Arrays opened from a b2z TreeStore/CTable are offset-backed leaves whose
        # urlpath points at the outer bundle, not at a standalone .b2nd file.
        # Reopening that path would materialize the whole TreeStore/CTable.
        if not urlpath.endswith(".b2z"):
            return blosc2.open(urlpath, mode="r", mmap_mode=_INDEX_MMAP_MODE)
    return where_x


def _gather_mmap_source(where_x):
    """Return a cached mmap handle for *where_x* for use in repeated gather operations.

    On Windows mmap is disabled (see ``_INDEX_MMAP_MODE``), so the original handle
    is returned unchanged.
    """
    if _INDEX_MMAP_MODE is None:
        return where_x
    urlpath = getattr(where_x, "urlpath", None)
    if not _supports_block_reads(where_x) or urlpath is None:
        return where_x
    _purge_stale_persistent_caches()
    urlpath = str(urlpath)
    # Drop the cached mapping if the file was overwritten in place since we mapped it.
    _refresh_persistent_caches(urlpath)
    handle = _GATHER_MMAP_HANDLES.get(urlpath)
    if handle is None:
        handle = blosc2.open(urlpath, mode="r", mmap_mode=_INDEX_MMAP_MODE)
        _GATHER_MMAP_HANDLES[urlpath] = handle
    return handle


def _bucket_match_from_span(span: np.ndarray, plan: IndexPlan) -> np.ndarray:
    if plan.target is not None and plan.target.get("source") == "expression":
        field_values = _values_from_numpy_target(span, plan.target)
    else:
        field_values = span if plan.field is None else span[plan.field]
    match = np.ones(len(field_values), dtype=bool)
    if plan.lower is not None:
        match &= field_values >= plan.lower if plan.lower_inclusive else field_values > plan.lower
    if plan.upper is not None:
        match &= field_values <= plan.upper if plan.upper_inclusive else field_values < plan.upper
    return match


def _process_bucket_chunk_batch(
    chunk_ids: np.ndarray,
    where_x,
    plan: IndexPlan,
    total_len: int,
    return_positions: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    value_parts = []
    position_parts = []
    local_where_x = _bucket_worker_source(where_x)
    for chunk_id in chunk_ids:
        bucket_mask = plan.bucket_masks[int(chunk_id)]
        chunk_start = int(chunk_id) * plan.chunk_len
        chunk_stop = min(chunk_start + plan.chunk_len, total_len)
        for run_start, run_stop in _contiguous_true_runs(np.asarray(bucket_mask, dtype=bool)):
            start = chunk_start + run_start * plan.bucket_len
            stop = min(chunk_start + run_stop * plan.bucket_len, chunk_stop)
            if start >= stop:
                continue
            if _supports_block_reads(local_where_x):
                span = np.empty(stop - start, dtype=local_where_x.dtype)
                _read_ndarray_linear_span(local_where_x, start, span)
            else:
                span = local_where_x[start:stop]
            match = _bucket_match_from_span(span, plan)
            if np.any(match):
                value_parts.append(np.require(span[match], requirements="C"))
                if return_positions:
                    position_parts.append(np.flatnonzero(match).astype(np.int64, copy=False) + start)
    if return_positions:
        return _merge_value_position_batches(
            value_parts, position_parts, _bucket_batch_result_dtype(where_x)
        )
    if not value_parts:
        return np.empty(0, dtype=_bucket_batch_result_dtype(where_x))
    return np.concatenate(value_parts) if len(value_parts) > 1 else value_parts[0]


def _merge_result_batches(parts: list[np.ndarray], dtype: np.dtype) -> np.ndarray:
    parts = [part for part in parts if len(part) > 0]
    if not parts:
        return np.empty(0, dtype=dtype)
    return np.concatenate(parts) if len(parts) > 1 else parts[0]


def _merge_value_position_batches(
    value_batches: list[np.ndarray], position_batches: list[np.ndarray], dtype: np.dtype
) -> tuple[np.ndarray, np.ndarray]:
    return _merge_result_batches(value_batches, dtype), _merge_position_batches(position_batches)


def _merge_segment_query_batches(
    parts: list[np.ndarray] | list[tuple[np.ndarray, np.ndarray]],
    dtype: np.dtype,
    *,
    return_positions: bool,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    if return_positions:
        value_batches = []
        position_batches = []
        for values, positions in parts:
            if len(values) > 0:
                value_batches.append(values)
            if len(positions) > 0:
                position_batches.append(positions)
        return _merge_value_position_batches(value_batches, position_batches, dtype)

    value_batches = [part for part in parts if len(part) > 0]
    if value_batches:
        return np.concatenate(value_batches) if len(value_batches) > 1 else value_batches[0]
    return np.empty(0, dtype=dtype)


def _process_segment_query_batch(
    units: np.ndarray,
    expression: str,
    operands: dict,
    ne_args: dict,
    where: dict,
    plan: IndexPlan,
    result_dtype: np.dtype,
    return_positions: bool,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    from .lazyexpr import _get_result, ne_evaluate
    from .utils import get_chunk_operands

    chunk_operands = {}
    value_parts = []
    position_parts = []
    for unit in units:
        start = int(unit) * plan.segment_len
        stop = min(start + plan.segment_len, plan.base.shape[0])
        cslice = (slice(start, stop, 1),)
        get_chunk_operands(operands, cslice, chunk_operands, plan.base.shape)
        if return_positions:
            match = ne_evaluate(expression, chunk_operands, **ne_args)
            if np.any(match):
                value_parts.append(np.require(chunk_operands["_where_x"][match], requirements="C"))
                absolute = np.arange(start, stop, dtype=np.int64)
                position_parts.append(absolute[match])
        else:
            result, _ = _get_result(expression, chunk_operands, ne_args, where)
            if len(result) > 0:
                value_parts.append(np.require(result, requirements="C"))
    if return_positions:
        return _merge_value_position_batches(value_parts, position_parts, result_dtype)
    return _merge_result_batches(value_parts, result_dtype)


def _reduced_positions_from_cython_batches(
    candidate_chunk_ids: np.ndarray, thread_count: int, process_batch
) -> tuple[np.ndarray, int]:
    return _run_position_batches(candidate_chunk_ids, thread_count, process_batch)


def _reduced_positions_from_python_batches(
    candidate_chunk_ids: np.ndarray, thread_count: int, process_batch
) -> tuple[list[np.ndarray], int]:
    if thread_count <= 1:
        return process_batch(candidate_chunk_ids)
    parts = []
    total_candidate_segments = 0
    batches = _chunk_batches(candidate_chunk_ids, thread_count)
    with ThreadPoolExecutor(max_workers=thread_count) as executor:
        for batch_parts, batch_candidate_segments in executor.map(process_batch, batches):
            total_candidate_segments += batch_candidate_segments
            parts.extend(batch_parts)
    return parts, total_candidate_segments


def _sorted_boundary_search_bounds(boundaries: np.ndarray, plan: ExactPredicatePlan) -> tuple[int, int]:
    if len(boundaries) == 0:
        return 0, 0
    starts = boundaries["start"]
    ends = boundaries["end"]
    try:
        lo, hi = indexing_ext.index_search_boundary_bounds(
            starts, ends, plan.lower, plan.lower_inclusive, plan.upper, plan.upper_inclusive
        )
    except TypeError:
        lo = 0
        hi = len(boundaries)
        if plan.lower is not None:
            lo = int(np.searchsorted(ends, plan.lower, side="left" if plan.lower_inclusive else "right"))
        if plan.upper is not None:
            hi = int(np.searchsorted(starts, plan.upper, side="right" if plan.upper_inclusive else "left"))
    if lo < 0:
        lo = 0
    if hi > len(boundaries):
        hi = len(boundaries)
    return lo, hi


def _bucket_search_plan(
    plan: ExactPredicatePlan, dtype: np.dtype, value_lossy_bits: int
) -> ExactPredicatePlan:
    if value_lossy_bits <= 0 or plan.lower is None:
        return plan
    if dtype.kind in {"i", "u"}:
        next_lower = plan.lower if plan.lower_inclusive else min(int(plan.lower) + 1, np.iinfo(dtype).max)
    else:
        next_lower = (
            plan.lower
            if plan.lower_inclusive
            else np.nextafter(np.asarray(plan.lower, dtype=dtype)[()], np.inf)
        )
    return ExactPredicatePlan(
        base=plan.base,
        descriptor=plan.descriptor,
        target=plan.target,
        field=plan.field,
        lower=_quantize_bucket_value_scalar(next_lower, dtype, value_lossy_bits),
        lower_inclusive=True,
        upper=plan.upper,
        upper_inclusive=plan.upper_inclusive,
    )


def _bucket_masks_from_bucket_chunk_nav_ooc(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> tuple[np.ndarray, int, int]:
    bucket = descriptor["bucket"]
    dtype = np.dtype(descriptor["dtype"])
    value_lossy_bits = int(bucket.get("value_lossy_bits", 0))
    search_plan = _bucket_search_plan(plan, dtype, value_lossy_bits)
    offsets_handle = _load_bucket_offsets_handle(array, descriptor)
    l1_handle = _load_bucket_l1_handle(array, descriptor)
    candidate_chunks = _candidate_units_from_boundaries_handle(l1_handle, search_plan)
    bucket_masks = np.zeros((int(l1_handle.shape[0]), int(bucket["bucket_count"])), dtype=bool)
    if not np.any(candidate_chunks):
        return bucket_masks, 0, 0

    values_sidecar, bucket_sidecar, l2_sidecar = _load_bucket_sidecar_handles(array, descriptor)
    chunk_len = int(bucket["chunk_len"])
    nav_segment_len = int(bucket["nav_segment_len"])
    nsegments_per_chunk = int(bucket["nsegments_per_chunk"])
    bucket_dtype = np.dtype(bucket.get("bucket_dtype", np.uint16))
    total_candidate_segments = 0
    candidate_chunk_ids = np.flatnonzero(candidate_chunks).astype(np.intp, copy=False)

    def process_batch(chunk_ids: np.ndarray) -> tuple[list[tuple[int, np.ndarray]], int]:
        if len(chunk_ids) == 0:
            return [], 0
        batch_values = (
            values_sidecar
            if bucket.get("values_path") is None
            else _open_sidecar_file(bucket["values_path"], _INDEX_MMAP_MODE)
        )
        batch_buckets = (
            bucket_sidecar
            if bucket.get("bucket_positions_path") is None
            else _open_sidecar_file(bucket["bucket_positions_path"], _INDEX_MMAP_MODE)
        )
        batch_l2 = (
            l2_sidecar
            if bucket.get("l2_path") is None
            else _open_sidecar_file(bucket["l2_path"], _INDEX_MMAP_MODE)
        )
        batch_results = []
        batch_candidate_segments = 0
        l2_row = np.empty(nsegments_per_chunk, dtype=_boundary_dtype(dtype))
        span_values = np.empty(chunk_len, dtype=dtype)
        bucket_ids = np.empty(chunk_len, dtype=bucket_dtype)
        for chunk_id in chunk_ids:
            offset_start, offset_stop = _read_offset_pair(offsets_handle, int(chunk_id))
            chunk_items = offset_stop - offset_start
            segment_count = _segment_row_count(chunk_items, nav_segment_len)
            batch_l2.get_1d_span_numpy(l2_row, int(chunk_id), 0, nsegments_per_chunk)
            segment_runs, candidate_segments = _chunk_nav_candidate_runs(l2_row, segment_count, plan)
            batch_candidate_segments += candidate_segments
            if not segment_runs:
                continue
            matched_buckets = np.zeros(int(bucket["bucket_count"]), dtype=bool)
            for seg_start_idx, seg_stop_idx in segment_runs:
                local_start = seg_start_idx * nav_segment_len
                local_stop = min(seg_stop_idx * nav_segment_len, chunk_items)
                span_items = local_stop - local_start
                values_view = span_values[:span_items]
                batch_values.get_1d_span_numpy(values_view, int(chunk_id), local_start, span_items)
                lo, hi = _search_bounds(values_view, search_plan)
                if lo >= hi:
                    continue
                bucket_view = bucket_ids[: hi - lo]
                batch_buckets.get_1d_span_numpy(bucket_view, int(chunk_id), local_start + lo, hi - lo)
                matched_buckets[bucket_view.astype(np.intp, copy=False)] = True
            if np.any(matched_buckets):
                batch_results.append((int(chunk_id), matched_buckets))
        return batch_results, batch_candidate_segments

    thread_count = _index_query_thread_count(len(candidate_chunk_ids))
    if thread_count <= 1:
        batch_results, total_candidate_segments = process_batch(candidate_chunk_ids)
        for chunk_id, matched_buckets in batch_results:
            bucket_masks[chunk_id] = matched_buckets
    else:
        batches = _chunk_batches(candidate_chunk_ids, thread_count)
        with ThreadPoolExecutor(max_workers=thread_count) as executor:
            for batch_results, batch_candidate_segments in executor.map(process_batch, batches):
                total_candidate_segments += batch_candidate_segments
                for chunk_id, matched_buckets in batch_results:
                    bucket_masks[chunk_id] = matched_buckets

    return bucket_masks, int(np.count_nonzero(candidate_chunks)), total_candidate_segments


def _exact_positions_from_partial_chunk_nav_ooc(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> tuple[np.ndarray, int, int]:
    partial = descriptor["partial"]
    offsets_handle = _load_partial_offsets_handle(array, descriptor)
    l1_handle = _load_partial_l1_handle(array, descriptor)
    candidate_chunks = _candidate_units_from_boundaries_handle(l1_handle, plan)
    if not np.any(candidate_chunks):
        return np.empty(0, dtype=np.int64), 0, 0

    dtype = np.dtype(descriptor["dtype"])
    chunk_len = int(partial["chunk_len"])
    nav_segment_len = int(partial["nav_segment_len"])
    nsegments_per_chunk = int(partial["nsegments_per_chunk"])
    local_position_dtype = np.dtype(partial.get("position_dtype", np.uint32))
    candidate_chunk_ids = np.flatnonzero(candidate_chunks).astype(np.intp, copy=False)
    l2_boundary_dtype = _boundary_dtype(dtype)
    values_sidecar, positions_sidecar, l2_sidecar = _load_partial_sidecar_handles(array, descriptor)
    thread_count = _index_query_thread_count(len(candidate_chunk_ids))

    try:
        positions, total_candidate_segments = _partial_chunk_nav_positions_cython(
            partial,
            offsets_handle,
            candidate_chunk_ids,
            thread_count,
            dtype,
            chunk_len,
            nav_segment_len,
            nsegments_per_chunk,
            local_position_dtype,
            l2_boundary_dtype,
            plan,
        )
        if len(positions) == 0:
            return np.empty(0, dtype=np.int64), int(candidate_chunk_ids.size), total_candidate_segments
        return np.sort(positions, kind="stable"), int(candidate_chunk_ids.size), total_candidate_segments
    except TypeError:
        pass

    parts, total_candidate_segments = _partial_chunk_nav_positions_python(
        partial,
        offsets_handle,
        candidate_chunk_ids,
        thread_count,
        dtype,
        chunk_len,
        nav_segment_len,
        nsegments_per_chunk,
        local_position_dtype,
        l2_boundary_dtype,
        values_sidecar,
        positions_sidecar,
        l2_sidecar,
        plan,
    )

    if not parts:
        return np.empty(0, dtype=np.int64), int(candidate_chunk_ids.size), total_candidate_segments
    positions = np.concatenate(parts) if len(parts) > 1 else parts[0]
    return (
        np.sort(positions, kind="stable"),
        int(candidate_chunk_ids.size),
        total_candidate_segments,
    )


def _partial_chunk_nav_positions_cython(
    partial: dict,
    offsets_handle,
    candidate_chunk_ids: np.ndarray,
    thread_count: int,
    dtype: np.dtype,
    chunk_len: int,
    nav_segment_len: int,
    nsegments_per_chunk: int,
    local_position_dtype: np.dtype,
    l2_boundary_dtype: np.dtype,
    plan: ExactPredicatePlan,
) -> tuple[np.ndarray, int]:
    if partial.get("values_path") is None:
        raise TypeError("cython chunk-nav path requires reopenable sidecars")

    offsets = _read_sidecar_span(offsets_handle, 0, int(offsets_handle.shape[0]))

    def process_cython_batch(chunk_ids: np.ndarray) -> tuple[np.ndarray, int]:
        if len(chunk_ids) == 0:
            return np.empty(0, dtype=np.int64), 0
        batch_values = blosc2.open(partial["values_path"], mode="r", mmap_mode=_INDEX_MMAP_MODE)
        batch_positions = blosc2.open(partial["positions_path"], mode="r", mmap_mode=_INDEX_MMAP_MODE)
        batch_l2 = blosc2.open(partial["l2_path"], mode="r", mmap_mode=_INDEX_MMAP_MODE)
        batch_l2_row = np.empty(nsegments_per_chunk, dtype=l2_boundary_dtype)
        batch_span_values = np.empty(chunk_len, dtype=dtype)
        batch_local_positions = np.empty(chunk_len, dtype=local_position_dtype)
        return indexing_ext.index_collect_reduced_chunk_nav_positions(
            offsets,
            chunk_ids,
            batch_values,
            batch_positions,
            batch_l2,
            batch_l2_row,
            batch_span_values,
            batch_local_positions,
            chunk_len,
            nav_segment_len,
            nsegments_per_chunk,
            plan.lower,
            plan.lower_inclusive,
            plan.upper,
            plan.upper_inclusive,
        )

    return _reduced_positions_from_cython_batches(candidate_chunk_ids, thread_count, process_cython_batch)


def _partial_chunk_nav_positions_python(
    partial: dict,
    offsets_handle,
    candidate_chunk_ids: np.ndarray,
    thread_count: int,
    dtype: np.dtype,
    chunk_len: int,
    nav_segment_len: int,
    nsegments_per_chunk: int,
    local_position_dtype: np.dtype,
    l2_boundary_dtype: np.dtype,
    values_sidecar,
    positions_sidecar,
    l2_sidecar,
    plan: ExactPredicatePlan,
) -> tuple[list[np.ndarray], int]:
    def process_batch(chunk_ids: np.ndarray) -> tuple[list[np.ndarray], int]:
        if len(chunk_ids) == 0:
            return [], 0
        batch_values = (
            values_sidecar
            if partial.get("values_path") is None
            else blosc2.open(partial["values_path"], mode="r", mmap_mode=_INDEX_MMAP_MODE)
        )
        batch_positions = (
            positions_sidecar
            if partial.get("positions_path") is None
            else blosc2.open(partial["positions_path"], mode="r", mmap_mode=_INDEX_MMAP_MODE)
        )
        batch_l2 = (
            l2_sidecar
            if partial.get("l2_path") is None
            else blosc2.open(partial["l2_path"], mode="r", mmap_mode=_INDEX_MMAP_MODE)
        )
        batch_parts = []
        batch_candidate_segments = 0
        l2_row = np.empty(nsegments_per_chunk, dtype=l2_boundary_dtype)
        span_values = np.empty(chunk_len, dtype=dtype)
        local_positions = np.empty(chunk_len, dtype=local_position_dtype)
        for chunk_id in chunk_ids:
            offset_start, offset_stop = _read_offset_pair(offsets_handle, int(chunk_id))
            chunk_items = offset_stop - offset_start
            segment_count = _segment_row_count(chunk_items, nav_segment_len)
            batch_l2.get_1d_span_numpy(l2_row, int(chunk_id), 0, nsegments_per_chunk)
            segment_runs, candidate_segments = _chunk_nav_candidate_runs(l2_row, segment_count, plan)
            batch_candidate_segments += candidate_segments
            if not segment_runs:
                continue
            for seg_start_idx, seg_stop_idx in segment_runs:
                local_start = seg_start_idx * nav_segment_len
                local_stop = min(seg_stop_idx * nav_segment_len, chunk_items)
                span_items = local_stop - local_start
                values_view = span_values[:span_items]
                batch_values.get_1d_span_numpy(values_view, int(chunk_id), local_start, span_items)
                lo, hi = _search_bounds(values_view, plan)
                if lo >= hi:
                    continue
                positions_view = local_positions[: hi - lo]
                batch_positions.get_1d_span_numpy(positions_view, int(chunk_id), local_start + lo, hi - lo)
                batch_parts.append(chunk_id * chunk_len + positions_view.astype(np.int64, copy=False))
        return batch_parts, batch_candidate_segments

    return _reduced_positions_from_python_batches(candidate_chunk_ids, thread_count, process_batch)


def _bit_count_sum(masks: np.ndarray) -> int:
    if masks.dtype == bool:
        return int(np.count_nonzero(masks))
    return sum(int(mask).bit_count() for mask in masks.tolist())


def _bucket_masks_from_bucket(
    array: blosc2.NDArray, descriptor: dict, plan: ExactPredicatePlan
) -> tuple[np.ndarray, int, int]:
    if _range_is_empty(plan):
        return np.empty((0, 0), dtype=bool), 0, 0
    return _bucket_masks_from_bucket_chunk_nav_ooc(array, descriptor, plan)


def _exact_positions_from_partial(
    array: blosc2.NDArray, descriptor: dict, dtype: np.dtype, plan: ExactPredicatePlan
) -> tuple[np.ndarray, int, int]:
    if _range_is_empty(plan):
        return np.empty(0, dtype=np.int64), 0, 0
    return _exact_positions_from_partial_chunk_nav_ooc(array, descriptor, plan)


def _exact_positions_from_plan(plan: ExactPredicatePlan) -> np.ndarray | None:
    kind = plan.descriptor["kind"]
    if kind == "full":
        return _exact_positions_from_full(plan.base, plan.descriptor, plan)
    if kind == "opsi":
        return _exact_positions_from_opsi(plan.base, plan.descriptor, plan)
    if kind == "partial":
        return _exact_positions_from_partial(
            plan.base, plan.descriptor, np.dtype(plan.descriptor["dtype"]), plan
        )[0]
    return None


def _multi_exact_positions(plans: list[ExactPredicatePlan]) -> tuple[blosc2.NDArray, np.ndarray] | None:
    if not plans:
        return None
    base = plans[0].base
    merged_by_target: dict[str, ExactPredicatePlan] = {}
    for plan in plans:
        if plan.base is not base:
            return None
        key = plan.descriptor["token"]
        current = merged_by_target.get(key)
        if current is None:
            merged_by_target[key] = plan
            continue
        merged = _merge_exact_plans(current, plan, "and")
        if merged is None:
            return None
        merged_by_target[key] = merged

    exact_arrays = []
    for plan in merged_by_target.values():
        positions = _exact_positions_from_plan(plan)
        if positions is None:
            return None
        exact_arrays.append(np.asarray(positions, dtype=np.int64))

    result = exact_arrays[0]
    for other in exact_arrays[1:]:
        result = np.intersect1d(result, other, assume_unique=False)
    return base, result


def _plan_cross_column_exact(plans: list[ExactPredicatePlan]) -> IndexPlan | None:
    """Try cross-column refinement when plans span different base NDArrays.

    When a conjunction has predicates on different columns (e.g. ``tips > 100
    AND km > 0 AND sec > 0``) and one of those columns has a FULL/PARTIAL
    index that can return compact exact positions, use those positions as a
    pre-filter.  The remaining comparison predicates are evaluated only on
    the candidate positions, skipping a full table scan.

    Returns an ``IndexPlan`` with ``partial_exact_positions`` and
    ``refine_plans`` when a compact candidate set is found, or ``None``.
    """
    threshold = _cross_column_threshold(plans)
    for plan in plans:
        positions = _exact_positions_from_plan(plan)
        if positions is None:
            continue
        positions = np.asarray(positions, dtype=np.int64)
        if len(positions) == 0 or len(positions) > threshold:
            continue
        # Compact exact positions found for this column's index.
        # The caller (``_try_index_where``) will evaluate the full
        # expression on these positions to refine against other
        # predicates.
        descriptor = _copy_descriptor(plan.descriptor)
        return IndexPlan(
            True,
            "cross-column exact refinement",
            descriptor=descriptor,
            base=plan.base,
            target=plan.descriptor.get("target"),
            field=plan.field,
            level=plan.descriptor["kind"],
            total_units=int(plan.base.shape[0]),
            selected_units=len(positions),
            partial_exact_positions=positions,
        )
    return None


def _cross_column_threshold(plans: list[ExactPredicatePlan]) -> int:
    """Return the maximum number of exact positions worth refining against.

    If the candidate set is too large the refinement cost (reading extra
    columns at every candidate position) outweighs the benefit over a
    full scan.  Cap at the smallest column size to stay conservative.
    """
    if not plans:
        return 100_000
    min_nrows = min(int(p.base.shape[0]) for p in plans)
    return min(100_000, max(1, min_nrows // 10))


def _plan_multi_exact_query(plans: list[ExactPredicatePlan]) -> IndexPlan | None:
    multi_exact = _multi_exact_positions(plans)
    if multi_exact is not None:
        base, exact_positions = multi_exact
        if len(exact_positions) >= int(base.shape[0]):
            return None
        descriptor = _copy_descriptor(plans[0].descriptor)
        lookup_path = None
        if descriptor["kind"] == "partial":
            lookup_path = (
                "chunk-nav-ooc"
                if _chunk_nav_supports_selective_ooc_lookup(base, descriptor, "partial")
                else "chunk-nav"
            )
        return IndexPlan(
            True,
            "multi-field positional indexes selected",
            descriptor=descriptor,
            base=base,
            target=plans[0].descriptor.get("target"),
            field=None,
            level="partial",
            total_units=int(base.shape[0]),
            selected_units=len(exact_positions),
            exact_positions=exact_positions,
            lookup_path=lookup_path,
        )
    # Cross-column fallback: try each plan's index individually.  If one
    # column's index produces compact exact positions we can use them as a
    # pre-filter and refine the remaining predicates cheaply.
    cross_col = _plan_cross_column_exact(plans)
    if cross_col is not None:
        return cross_col
    return None


def _plan_single_exact_query(exact_plan: ExactPredicatePlan) -> IndexPlan:
    kind = exact_plan.descriptor["kind"]
    if kind in {"full", "opsi"}:
        exact_positions = (
            _exact_positions_from_full(exact_plan.base, exact_plan.descriptor, exact_plan)
            if kind == "full"
            else _exact_positions_from_opsi(exact_plan.base, exact_plan.descriptor, exact_plan)
        )
        return IndexPlan(
            True,
            f"{kind} index selected",
            descriptor=_copy_descriptor(exact_plan.descriptor),
            base=exact_plan.base,
            target=exact_plan.descriptor.get("target"),
            field=exact_plan.field,
            level=kind,
            total_units=exact_plan.base.shape[0],
            selected_units=len(exact_positions),
            exact_positions=exact_positions,
        )
    if kind == "partial":
        dtype = np.dtype(exact_plan.descriptor["dtype"])
        exact_positions, candidate_chunks, candidate_nav_segments = _exact_positions_from_partial(
            exact_plan.base, exact_plan.descriptor, dtype, exact_plan
        )
        return IndexPlan(
            True,
            f"{kind} index selected",
            descriptor=_copy_descriptor(exact_plan.descriptor),
            base=exact_plan.base,
            target=exact_plan.descriptor.get("target"),
            field=exact_plan.field,
            level=kind,
            total_units=exact_plan.base.shape[0],
            selected_units=len(exact_positions),
            exact_positions=exact_positions,
            chunk_len=int(exact_plan.descriptor["partial"]["chunk_len"]),
            candidate_chunks=candidate_chunks,
            candidate_nav_segments=candidate_nav_segments,
            lookup_path="chunk-nav-ooc"
            if _chunk_nav_supports_selective_ooc_lookup(exact_plan.base, exact_plan.descriptor, "partial")
            else "chunk-nav",
        )
    bucket_masks, candidate_chunks, candidate_nav_segments = _bucket_masks_from_bucket(
        exact_plan.base, exact_plan.descriptor, exact_plan
    )
    bucket = exact_plan.descriptor["bucket"]
    total_units = bucket_masks.size
    selected_units = _bit_count_sum(bucket_masks)
    if selected_units < total_units:
        return IndexPlan(
            True,
            "bucket approximate-order index selected",
            descriptor=_copy_descriptor(exact_plan.descriptor),
            base=exact_plan.base,
            target=exact_plan.descriptor.get("target"),
            field=exact_plan.field,
            level=kind,
            total_units=total_units,
            selected_units=selected_units,
            bucket_masks=bucket_masks,
            bucket_len=int(bucket["bucket_len"]),
            chunk_len=int(bucket["chunk_len"]),
            lower=exact_plan.lower,
            lower_inclusive=exact_plan.lower_inclusive,
            upper=exact_plan.upper,
            upper_inclusive=exact_plan.upper_inclusive,
            candidate_chunks=candidate_chunks,
            candidate_nav_segments=candidate_nav_segments,
            lookup_path="chunk-nav-ooc"
            if _chunk_nav_supports_selective_ooc_lookup(exact_plan.base, exact_plan.descriptor, "bucket")
            else "chunk-nav",
        )
    return IndexPlan(False, "available positional index does not prune any units for this predicate")


def plan_query(
    expression: str,
    operands: dict,
    where: dict | None,
    *,
    use_index: bool = True,
    array_to_col: dict | None = None,
) -> IndexPlan:
    if not use_index:
        return IndexPlan(False, "index usage disabled for this query")
    if where is None or len(where) != 1:
        return IndexPlan(False, "indexing is only available for where(x) style filtering")

    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError:
        return IndexPlan(False, "expression is not valid Python syntax for planning")

    exact_terms = _plan_exact_conjunction(tree.body, operands, array_to_col)
    if exact_terms is not None and len(exact_terms) > 1:
        multi_exact_plan = _plan_multi_exact_query(exact_terms)
        if multi_exact_plan is not None:
            return multi_exact_plan

    exact_plan = _plan_exact_node(tree.body, operands, array_to_col)
    if exact_plan is not None:
        exact_query_plan = _plan_single_exact_query(exact_plan)
        if exact_query_plan.usable:
            return exact_query_plan

    # Cross-column refinement: the full expression tree has no single-plan
    # exact equivalent (different columns, only some indexed), but a
    # partial conjunction may give us compact exact positions from one
    # indexed column that we can refine against the other predicates.
    if exact_terms is not None and len(exact_terms) == 1:
        cross_col = _plan_cross_column_exact(exact_terms)
        if cross_col is not None:
            return cross_col

    segment_plan = _plan_segment_node(tree.body, operands, array_to_col)
    if segment_plan is None:
        return IndexPlan(False, "no usable index was found for this predicate")

    total_units = len(segment_plan.candidate_units)
    selected_units = int(np.count_nonzero(segment_plan.candidate_units))
    if selected_units == total_units:
        return IndexPlan(
            False,
            "available index does not prune any units for this predicate",
            descriptor=_copy_descriptor(segment_plan.descriptor),
            base=segment_plan.base,
            target=segment_plan.descriptor.get("target"),
            field=segment_plan.field,
            level=segment_plan.level,
            segment_len=segment_plan.segment_len,
            candidate_units=segment_plan.candidate_units,
            total_units=total_units,
            selected_units=selected_units,
        )

    return IndexPlan(
        True,
        f"{segment_plan.level} summaries selected",
        descriptor=_copy_descriptor(segment_plan.descriptor),
        base=segment_plan.base,
        target=segment_plan.descriptor.get("target"),
        field=segment_plan.field,
        level=segment_plan.level,
        segment_len=segment_plan.segment_len,
        candidate_units=segment_plan.candidate_units,
        total_units=total_units,
        selected_units=selected_units,
    )


def _where_output_dtype(where_x) -> np.dtype:
    return where_x.dtype if hasattr(where_x, "dtype") else np.asarray(where_x).dtype


def evaluate_segment_query(
    expression: str,
    operands: dict,
    ne_args: dict,
    where: dict,
    plan: IndexPlan,
    *,
    return_positions: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    if plan.base is None or plan.candidate_units is None or plan.segment_len is None:
        raise ValueError("segment evaluation requires a segment-based plan")

    candidate_units = np.flatnonzero(plan.candidate_units).astype(np.intp, copy=False)
    result_dtype = _where_output_dtype(where["_where_x"])

    thread_count = _downstream_query_thread_count(len(candidate_units), plan)
    if thread_count <= 1:
        parts = [
            _process_segment_query_batch(
                candidate_units,
                expression,
                operands,
                ne_args,
                where,
                plan,
                result_dtype,
                return_positions=return_positions,
            )
        ]
    else:
        batches = _chunk_batches(candidate_units, thread_count)
        with ThreadPoolExecutor(max_workers=thread_count) as executor:
            parts = list(
                executor.map(
                    _process_segment_query_batch,
                    batches,
                    [expression] * len(batches),
                    [operands] * len(batches),
                    [ne_args] * len(batches),
                    [where] * len(batches),
                    [plan] * len(batches),
                    [result_dtype] * len(batches),
                    [return_positions] * len(batches),
                )
            )

    return _merge_segment_query_batches(parts, result_dtype, return_positions=return_positions)


def evaluate_bucket_query(
    expression: str,
    operands: dict,
    ne_args: dict,
    where: dict,
    plan: IndexPlan,
    *,
    return_positions: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    del expression, operands, ne_args

    if plan.base is None or plan.bucket_masks is None or plan.chunk_len is None or plan.bucket_len is None:
        raise ValueError("bucket evaluation requires bucket masks and chunk geometry")

    total_len = int(plan.base.shape[0])
    where_x = where["_where_x"]
    candidate_chunk_ids = np.flatnonzero(np.any(plan.bucket_masks, axis=1)).astype(np.intp, copy=False)
    result_dtype = _where_output_dtype(where["_where_x"])

    thread_count = _downstream_query_thread_count(len(candidate_chunk_ids), plan)
    if thread_count <= 1:
        parts = [
            _process_bucket_chunk_batch(candidate_chunk_ids, where_x, plan, total_len, return_positions)
        ]
    else:
        batches = _chunk_batches(candidate_chunk_ids, thread_count)
        with ThreadPoolExecutor(max_workers=thread_count) as executor:
            parts = list(
                executor.map(
                    _process_bucket_chunk_batch,
                    batches,
                    [where_x] * len(batches),
                    [plan] * len(batches),
                    [total_len] * len(batches),
                    [return_positions] * len(batches),
                )
            )

    if return_positions:
        value_batches = []
        position_batches = []
        for values, positions in parts:
            if len(values) > 0:
                value_batches.append(values)
            if len(positions) > 0:
                position_batches.append(positions)
        return _merge_value_position_batches(value_batches, position_batches, result_dtype)

    return _merge_result_batches(parts, result_dtype)


def _gather_positions(where_x, positions: np.ndarray) -> np.ndarray:
    if len(positions) == 0:
        return np.empty(0, dtype=_where_output_dtype(where_x))

    positions = np.asarray(positions, dtype=np.int64)
    breaks = np.nonzero(np.diff(positions) != 1)[0] + 1
    runs = np.split(positions, breaks)
    parts = []
    for run in runs:
        start = int(run[0])
        stop = int(run[-1]) + 1
        parts.append(where_x[start:stop])
    return np.concatenate(parts) if len(parts) > 1 else parts[0]


def _gather_positions_by_chunk(where_x, positions: np.ndarray, chunk_len: int) -> np.ndarray:
    if len(positions) == 0:
        return np.empty(0, dtype=_where_output_dtype(where_x))

    positions = np.asarray(positions, dtype=np.int64)
    output = np.empty(len(positions), dtype=_where_output_dtype(where_x))
    chunk_ids = positions // chunk_len
    breaks = np.nonzero(np.diff(chunk_ids) != 0)[0] + 1
    start_idx = 0
    for stop_idx in (*breaks, len(positions)):
        chunk_positions = positions[start_idx:stop_idx]
        chunk_id = int(chunk_ids[start_idx])
        chunk_start = chunk_id * chunk_len
        chunk_stop = chunk_start + chunk_len
        chunk_values = where_x[chunk_start:chunk_stop]
        local_positions = chunk_positions - chunk_start
        output[start_idx:stop_idx] = chunk_values[local_positions]
        start_idx = stop_idx
    return output


def _supports_block_reads(where_x) -> bool:
    return isinstance(where_x, blosc2.NDArray) and hasattr(where_x, "get_1d_span_numpy")


def _gather_positions_by_block(
    where_x, positions: np.ndarray, chunk_len: int, block_len: int, total_len: int
) -> np.ndarray:
    if len(positions) == 0:
        return np.empty(0, dtype=_where_output_dtype(where_x))
    if not _supports_block_reads(where_x):
        return _gather_positions_by_chunk(where_x, positions, chunk_len)

    positions = np.asarray(positions, dtype=np.int64)
    output = np.empty(len(positions), dtype=_where_output_dtype(where_x))
    chunk_ids = positions // chunk_len
    chunk_breaks = np.nonzero(np.diff(chunk_ids) != 0)[0] + 1
    chunk_start_idx = 0
    for chunk_stop_idx in (*chunk_breaks, len(positions)):
        chunk_positions = positions[chunk_start_idx:chunk_stop_idx]
        chunk_id = int(chunk_ids[chunk_start_idx])
        chunk_origin = chunk_id * chunk_len
        local_positions = chunk_positions - chunk_origin
        if np.any(np.diff(local_positions) < 0):
            order = np.argsort(local_positions, kind="stable")
            sorted_local_positions = local_positions[order]
        else:
            order = None
            sorted_local_positions = local_positions

        sorted_output = (
            output[chunk_start_idx:chunk_stop_idx]
            if order is None
            else np.empty(len(chunk_positions), dtype=output.dtype)
        )
        block_ids = sorted_local_positions // block_len
        block_breaks = np.nonzero(np.diff(block_ids) != 0)[0] + 1
        block_start_idx = 0
        for block_stop_idx in (*block_breaks, len(sorted_local_positions)):
            block_positions = sorted_local_positions[block_start_idx:block_stop_idx]
            span_start = int(block_positions[0])
            span_stop = int(block_positions[-1]) + 1
            span_items = span_stop - span_start
            span_values = np.empty(span_items, dtype=output.dtype)
            where_x.get_1d_span_numpy(span_values, chunk_id, span_start, span_items)
            sorted_output[block_start_idx:block_stop_idx] = span_values[block_positions - span_start]
            block_start_idx = block_stop_idx

        if order is None:
            output[chunk_start_idx:chunk_stop_idx] = sorted_output
        else:
            inverse = np.empty(len(order), dtype=np.intp)
            inverse[order] = np.arange(len(order), dtype=np.intp)
            output[chunk_start_idx:chunk_stop_idx] = sorted_output[inverse]
        chunk_start_idx = chunk_stop_idx
    return output


def evaluate_full_query(where: dict, plan: IndexPlan) -> np.ndarray:
    if plan.exact_positions is None:
        raise ValueError("full evaluation requires positional matches")
    if plan.base is not None:
        # Use a cached mmap handle when available so blosc2_schunk_get_lazychunk can return
        # a zero-copy pointer into the mapped region instead of malloc+pread per block.
        gather_source = _gather_mmap_source(where["_where_x"])
        block_gather_threshold = int(plan.base.blocks[0])
        if len(plan.exact_positions) <= block_gather_threshold:
            return _gather_positions_by_block(
                gather_source,
                plan.exact_positions,
                int(plan.base.chunks[0]),
                int(plan.base.blocks[0]),
                int(plan.base.shape[0]),
            )
        return _gather_positions_by_chunk(gather_source, plan.exact_positions, int(plan.base.chunks[0]))
    return _gather_positions(where["_where_x"], plan.exact_positions)


def _normalize_primary_order_target(array: blosc2.NDArray, order: str | None) -> tuple[dict, str | None]:
    if order is None:
        return _field_target_descriptor(None), None
    if array.dtype.fields is not None and order in array.dtype.fields:
        return _field_target_descriptor(order), order
    operands = array.fields if array.dtype.fields is not None else {SELF_TARGET_NAME: array}
    base, target, _ = _normalize_expression_target(order, operands)
    if base is not array:
        raise ValueError("ordered expressions must resolve to the target array")
    return target, None


def _full_run_count(descriptor: dict | None) -> int:
    if descriptor is None or descriptor.get("full") is None:
        return 0
    return len(descriptor["full"].get("runs", ()))


def _full_lookup_path(descriptor: dict | None, *, ordered: bool) -> str | None:
    if descriptor is None or descriptor.get("kind") != "full":
        return None
    if _full_run_count(descriptor):
        return "ordered-stream-merge" if ordered else "run-bounded-ooc"
    if ordered:
        return "ordered-stream"
    if (
        descriptor.get("persistent")
        and descriptor["full"].get("l1_path")
        and descriptor["full"].get("l2_path")
    ):
        return "compact-selective-ooc"
    return "sidecar-stream"


def _normalize_order_fields(
    array: blosc2.NDArray, order: str | list[str] | None
) -> tuple[dict, list[str | None]]:
    if order is None:
        if array.dtype.fields is None:
            return _field_target_descriptor(None), [None]
        return _field_target_descriptor(array.dtype.names[0]), list(array.dtype.names)
    if isinstance(order, list):
        fields = list(order)
    else:
        fields = [order]
    primary_target, primary_field = _normalize_primary_order_target(array, fields[0])
    normalized_order = [primary_field if primary_field is not None else fields[0]]
    if len(fields) > 1:
        if array.dtype.fields is None:
            raise ValueError("secondary order keys are only supported for structured arrays")
        for field in fields[1:]:
            if field not in array.dtype.fields:
                raise ValueError(f"field {field!r} is not present in the dtype")
        normalized_order.extend(fields[1:])
    return primary_target, normalized_order


def is_expression_order(array: blosc2.NDArray, order: str | list[str] | None) -> bool:
    if order is None:
        return False
    primary = order[0] if isinstance(order, list) else order
    try:
        target, _ = _normalize_primary_order_target(array, primary)
    except (TypeError, ValueError):
        return False
    return target["source"] == "expression"


def plan_array_order(
    array: blosc2.NDArray, order: str | list[str] | None = None, *, require_full: bool = False
) -> OrderedIndexPlan:
    try:
        primary_target, order_fields = _normalize_order_fields(array, order)
    except (TypeError, ValueError) as exc:
        return OrderedIndexPlan(False, str(exc))
    primary_field = _target_field(primary_target)
    descriptor = _full_descriptor_for_order(array, primary_target)
    if descriptor is None:
        if require_full:
            label = primary_field if primary_field is not None else primary_target.get("expression")
            return OrderedIndexPlan(False, f"order target {label!r} must have an associated full index")
        return OrderedIndexPlan(False, "no matching full index was found for ordered access")
    return OrderedIndexPlan(
        True,
        "ordered access will reuse a full index",
        descriptor=_copy_descriptor(descriptor),
        base=array,
        field=primary_field,
        order_fields=order_fields,
        total_rows=int(array.shape[0]),
        selected_rows=int(array.shape[0]),
        secondary_refinement=len(order_fields) > 1,
    )


def _positions_in_input_order(
    positions: np.ndarray, start: int | None, stop: int | None, step: int | None
) -> np.ndarray:
    if step is None:
        step = 1
    if step == 0:
        raise ValueError("step cannot be zero")
    return positions[slice(start, stop, step)]


def _full_descriptor_for_order(array: blosc2.NDArray, target: dict) -> dict | None:
    descriptor = _descriptor_for_target(array, target)
    if descriptor is None or descriptor.get("kind") != "full":
        return None
    return descriptor


def _equal_primary_values(left, right, dtype: np.dtype) -> bool:
    return _scalar_compare(left, right, dtype) == 0


def _refine_secondary_order(
    array: blosc2.NDArray,
    positions: np.ndarray,
    primary_values: np.ndarray,
    primary_dtype: np.dtype,
    secondary_fields: list[str],
) -> np.ndarray:
    if not secondary_fields or len(positions) <= 1:
        return positions

    refined = positions.copy()
    start = 0
    while start < len(refined):
        stop = start + 1
        while stop < len(refined) and _equal_primary_values(
            primary_values[start], primary_values[stop], primary_dtype
        ):
            stop += 1
        if stop - start > 1:
            tied_positions = refined[start:stop]
            tied_rows = array[tied_positions]
            tie_order = np.argsort(tied_rows, order=secondary_fields, kind="stable")
            refined[start:stop] = tied_positions[tie_order]
        start = stop
    return refined


def _concat_order_parts(parts: list[np.ndarray], dtype: np.dtype) -> np.ndarray:
    if not parts:
        return np.empty(0, dtype=dtype)
    return np.concatenate(parts) if len(parts) > 1 else np.asarray(parts[0], dtype=dtype)


def _ordered_selection_filter(
    exact_positions: np.ndarray, total_rows: int
) -> tuple[np.ndarray | None, set[int] | None]:
    normalized = np.asarray(exact_positions, dtype=np.int64)
    if len(normalized) == 0:
        return np.empty(0, dtype=np.int64), set()
    unique_positions = np.unique(normalized)
    if len(unique_positions) == total_rows:
        return None, None
    return unique_positions, set(unique_positions.tolist())


def _ordered_positions_from_single_sorted_handle(
    values_handle,
    positions_handle,
    dtype: np.dtype,
    exact_positions: np.ndarray,
    total_rows: int,
    need_values: bool,
) -> tuple[np.ndarray, np.ndarray]:
    selected_positions, _ = _ordered_selection_filter(exact_positions, total_rows)
    length = int(positions_handle.shape[0])
    chunk_len = int(positions_handle.chunks[0]) if hasattr(positions_handle, "chunks") else length
    position_parts = []
    value_parts = [] if need_values else None

    for start in range(0, length, chunk_len):
        stop = min(start + chunk_len, length)
        chunk_positions = _read_sidecar_span(positions_handle, start, stop).astype(np.int64, copy=False)
        if selected_positions is None:
            mask = None
            kept_positions = chunk_positions
        else:
            mask = np.isin(chunk_positions, selected_positions, assume_unique=True)
            if not np.any(mask):
                continue
            kept_positions = chunk_positions[mask]
        position_parts.append(kept_positions)
        if need_values:
            chunk_values = _read_sidecar_span(values_handle, start, stop)
            value_parts.append(chunk_values if mask is None else chunk_values[mask])

    positions = _concat_order_parts(position_parts, np.dtype(np.int64))
    if not need_values:
        return positions, np.empty(0, dtype=dtype)
    return positions, _concat_order_parts(value_parts, dtype)


class _SortedSidecarCursor:
    def __init__(self, values_handle, positions_handle, dtype: np.dtype):
        self.values_handle = values_handle
        self.positions_handle = positions_handle
        self.dtype = np.dtype(dtype)
        self.length = int(values_handle.shape[0])
        self.chunk_len = int(values_handle.chunks[0]) if hasattr(values_handle, "chunks") else self.length
        self.offset = 0
        self.local_index = 0
        self.values_chunk = np.empty(0, dtype=self.dtype)
        self.positions_chunk = np.empty(0, dtype=np.int64)
        self._fill_chunk()

    def _fill_chunk(self) -> None:
        if self.offset >= self.length:
            self.values_chunk = np.empty(0, dtype=self.dtype)
            self.positions_chunk = np.empty(0, dtype=np.int64)
            self.local_index = 0
            return
        stop = min(self.offset + self.chunk_len, self.length)
        self.values_chunk = _read_sidecar_span(self.values_handle, self.offset, stop)
        self.positions_chunk = _read_sidecar_span(self.positions_handle, self.offset, stop).astype(
            np.int64, copy=False
        )
        self.offset = stop
        self.local_index = 0

    @property
    def exhausted(self) -> bool:
        return len(self.values_chunk) == 0

    def current_value(self):
        return self.values_chunk[self.local_index]

    def current_position(self) -> int:
        return int(self.positions_chunk[self.local_index])

    def advance(self) -> None:
        self.local_index += 1
        if self.local_index >= len(self.values_chunk):
            self._fill_chunk()


def _ordered_positions_from_stream_merge(
    array: blosc2.NDArray,
    descriptor: dict,
    exact_positions: np.ndarray,
    need_values: bool,
) -> tuple[np.ndarray, np.ndarray]:
    dtype = np.dtype(descriptor["dtype"])
    total_rows = int(array.shape[0])
    selected_positions, remaining = _ordered_selection_filter(exact_positions, total_rows)
    full = descriptor["full"]
    streams = [_load_full_sidecar_handles(array, descriptor)]
    for run in full.get("runs", ()):
        streams.append(_load_full_run_sidecar_handles(array, descriptor, run))

    if len(streams) == 1:
        return _ordered_positions_from_single_sorted_handle(
            streams[0][0],
            streams[0][1],
            dtype,
            exact_positions,
            total_rows,
            need_values,
        )

    cursors = [
        _SortedSidecarCursor(values_handle, positions_handle, dtype)
        for values_handle, positions_handle in streams
    ]
    cursors = [cursor for cursor in cursors if not cursor.exhausted]
    if not cursors:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=dtype)

    positions = []
    values = [] if need_values else None

    while cursors:
        best_idx = 0
        best_cursor = cursors[0]
        best_value = best_cursor.current_value()
        best_position = best_cursor.current_position()
        for idx in range(1, len(cursors)):
            candidate = cursors[idx]
            candidate_value = candidate.current_value()
            candidate_position = candidate.current_position()
            if _pair_le(candidate_value, candidate_position, best_value, best_position, dtype):
                best_idx = idx
                best_cursor = candidate
                best_value = candidate_value
                best_position = candidate_position

        if selected_positions is None or best_position in remaining:
            positions.append(best_position)
            if need_values:
                values.append(best_value)
            if remaining is not None:
                remaining.remove(best_position)
                if not remaining:
                    break

        best_cursor.advance()
        if best_cursor.exhausted:
            cursors.pop(best_idx)

    ordered_positions = np.asarray(positions, dtype=np.int64)
    if not need_values:
        return ordered_positions, np.empty(0, dtype=dtype)
    return ordered_positions, np.asarray(values, dtype=dtype)


def _ordered_positions_from_exact_positions(
    array: blosc2.NDArray, descriptor: dict, exact_positions: np.ndarray, order_fields: list[str | None]
) -> np.ndarray:
    secondary_fields = [field for field in order_fields[1:] if field is not None]
    selected_positions, selected_values = _ordered_positions_from_stream_merge(
        array, descriptor, exact_positions, need_values=bool(secondary_fields)
    )
    if secondary_fields:
        selected_positions = _refine_secondary_order(
            array, selected_positions, selected_values, np.dtype(descriptor["dtype"]), secondary_fields
        )
    return selected_positions


def _shortcut_ordered_indices(
    array: blosc2.NDArray,
    descriptor: dict,
    total_rows: int,
    start: int | None,
    stop: int | None,
    step: int | None,
    order_fields: list[str | None],
) -> np.ndarray | None:
    """Fast path: read head or tail of a FULL index directly.

    Returns positions array, or ``None`` if the shortcut does not apply
    (multi-field order, expression index, strided or mid-array slice).
    """
    if len(order_fields) > 1:
        return None  # secondary refinement needed, fall back
    full = descriptor.get("full")
    if full is None:
        return None

    step = 1 if step is None else step
    if abs(step) != 1:
        return None  # strided slice, fall back

    # Use Python range slicing (O(1) len + indexing) to determine the exact
    # set of indices the slice covers, without materializing them.
    actual = range(total_rows)[start:stop:step]
    k = len(actual)
    if k == 0:
        return np.empty(0, dtype=np.int64)
    first = actual[0]
    last = actual[-1]

    if step > 0:
        sidecar_start, sidecar_stop = first, last + 1
    else:
        sidecar_start, sidecar_stop = last, first + 1

    return _read_sidecar_range(array, descriptor, sidecar_start, sidecar_stop, step)


def _read_sidecar_range(
    array: blosc2.NDArray, descriptor: dict, sidecar_start: int, sidecar_stop: int, step: int
) -> np.ndarray:
    """Read positions [*sidecar_start*:*sidecar_stop*] from the FULL index."""
    full = descriptor["full"]
    streams = [_load_full_sidecar_handles(array, descriptor)]
    for run in full.get("runs", ()):
        streams.append(_load_full_run_sidecar_handles(array, descriptor, run))

    if len(streams) != 1:
        # Multi-run index (after append + rebuild): not worth the complexity.
        # Fall back to the normal path, which is correct just slower.
        return None

    # Single run: any contiguous slice is a direct read.
    _, pos_handle = streams[0]
    chunk = pos_handle[sidecar_start:sidecar_stop][:]
    if step < 0:
        chunk = chunk[::-1].copy()
    return np.asarray(chunk, dtype=np.int64)


def ordered_indices(
    array: blosc2.NDArray,
    order: str | list[str] | None = None,
    *,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    require_full: bool = False,
) -> np.ndarray | None:
    ordered_plan = plan_array_order(array, order=order, require_full=require_full)
    if not ordered_plan.usable:
        if require_full:
            raise ValueError(ordered_plan.reason)
        return None
    order_fields = ordered_plan.order_fields
    descriptor = ordered_plan.descriptor
    total_rows = int(array.shape[0])

    shortcut = _shortcut_ordered_indices(array, descriptor, total_rows, start, stop, step, order_fields)
    if shortcut is not None:
        return shortcut

    positions = _ordered_positions_from_exact_positions(
        array, descriptor, np.arange(total_rows, dtype=np.int64), order_fields
    )
    return _positions_in_input_order(positions, start, stop, step)


def plan_ordered_query(
    expression: str, operands: dict, where: dict, order: str | list[str]
) -> OrderedIndexPlan:
    if len(where) != 1:
        return OrderedIndexPlan(False, "ordered index reuse is only available for where(x) style filtering")
    base = where["_where_x"]
    if not isinstance(base, blosc2.NDArray) or base.ndim != 1:
        return OrderedIndexPlan(False, "ordered index reuse requires a 1-D NDArray target")

    base_order_plan = plan_array_order(base, order=order, require_full=False)
    if not base_order_plan.usable:
        return base_order_plan

    filter_plan = plan_query(expression, operands, where, use_index=True)
    if not filter_plan.usable:
        return OrderedIndexPlan(
            False, f"ordered access cannot reuse an index because filtering does not: {filter_plan.reason}"
        )
    if filter_plan.base is not base or filter_plan.exact_positions is None:
        return OrderedIndexPlan(
            False, "ordered access currently requires positional row matches from filtering"
        )

    return OrderedIndexPlan(
        True,
        "ordered access will reuse a full index after positional filtering",
        descriptor=base_order_plan.descriptor,
        base=base,
        field=base_order_plan.field,
        order_fields=base_order_plan.order_fields,
        total_rows=int(base.shape[0]),
        selected_rows=len(filter_plan.exact_positions),
        secondary_refinement=base_order_plan.secondary_refinement,
    )


def ordered_query_indices(
    expression: str,
    operands: dict,
    where: dict,
    order: str | list[str],
    *,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
) -> np.ndarray | None:
    ordered_plan = plan_ordered_query(expression, operands, where, order)
    if not ordered_plan.usable:
        return None
    base = ordered_plan.base
    order_fields = ordered_plan.order_fields
    descriptor = ordered_plan.descriptor

    plan = plan_query(expression, operands, where, use_index=True)

    positions = _ordered_positions_from_exact_positions(base, descriptor, plan.exact_positions, order_fields)
    return _positions_in_input_order(positions, start, stop, step)


def read_sorted(
    array: blosc2.NDArray,
    order: str | list[str] | None = None,
    *,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    require_full: bool = False,
) -> np.ndarray | None:
    positions = ordered_indices(
        array, order=order, start=start, stop=stop, step=step, require_full=require_full
    )
    if positions is None:
        return None
    return _gather_positions_by_block(
        array, positions, int(array.chunks[0]), int(array.blocks[0]), int(array.shape[0])
    )


def iter_sorted(
    array: blosc2.NDArray,
    order: str | list[str] | None = None,
    *,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    batch_size: int | None = None,
) -> np.ndarray:
    positions = ordered_indices(array, order=order, start=start, stop=stop, step=step, require_full=True)
    if batch_size is None:
        batch_size = max(1, int(array.blocks[0]))
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    for idx in range(0, len(positions), batch_size):
        batch = _gather_positions_by_block(
            array,
            positions[idx : idx + batch_size],
            int(array.chunks[0]),
            int(array.blocks[0]),
            int(array.shape[0]),
        )
        yield from batch


def will_use_index(expr) -> bool:
    where = getattr(expr, "_where_args", None)
    order = getattr(expr, "_order", None)
    if order is not None:
        return plan_ordered_query(expr.expression, expr.operands, where, order).usable
    return plan_query(expr.expression, expr.operands, where).usable


def explain_query(expr) -> dict:
    """Return planning details for a lazy query.

    This is an internal helper behind :meth:`blosc2.LazyExpr.explain`. The
    returned mapping summarizes whether indexing can be used, which index kind
    was selected, and additional diagnostics such as candidate counts and the
    lookup path chosen for ``full`` indexes.
    """
    where = getattr(expr, "_where_args", None)
    order = getattr(expr, "_order", None)
    if order is not None:
        ordered_plan = plan_ordered_query(expr.expression, expr.operands, where, order)
        filter_plan = plan_query(expr.expression, expr.operands, where)
        return {
            "will_use_index": ordered_plan.usable,
            "reason": ordered_plan.reason,
            "target": None if ordered_plan.descriptor is None else ordered_plan.descriptor.get("target"),
            "field": ordered_plan.field,
            "kind": None if ordered_plan.descriptor is None else ordered_plan.descriptor["kind"],
            "level": "full" if ordered_plan.usable else None,
            "ordered_access": True,
            "order": ordered_plan.order_fields,
            "secondary_refinement": ordered_plan.secondary_refinement,
            "candidate_units": ordered_plan.selected_rows,
            "total_units": ordered_plan.total_rows,
            "candidate_chunks": ordered_plan.selected_rows,
            "total_chunks": ordered_plan.total_rows,
            "exact_rows": ordered_plan.selected_rows if ordered_plan.usable else None,
            "filter_reason": filter_plan.reason,
            "filter_level": filter_plan.level,
            "full_runs": _full_run_count(ordered_plan.descriptor),
            "lookup_path": _full_lookup_path(ordered_plan.descriptor, ordered=True),
            "descriptor": ordered_plan.descriptor,
        }

    plan = plan_query(expr.expression, expr.operands, where)
    return {
        "will_use_index": plan.usable,
        "reason": plan.reason,
        "target": None if plan.descriptor is None else plan.descriptor.get("target"),
        "field": plan.field,
        "kind": None if plan.descriptor is None else plan.descriptor["kind"],
        "level": plan.level,
        "ordered_access": False,
        "order": None,
        "secondary_refinement": False,
        "candidate_units": plan.selected_units,
        "total_units": plan.total_units,
        "candidate_chunks": plan.candidate_chunks if plan.candidate_chunks else plan.selected_units,
        "total_chunks": plan.total_units,
        "candidate_nav_segments": plan.candidate_nav_segments or None,
        "candidate_base_spans": plan.candidate_base_spans or None,
        "exact_rows": None if plan.exact_positions is None else len(plan.exact_positions),
        "full_runs": _full_run_count(plan.descriptor),
        "lookup_path": plan.lookup_path or _full_lookup_path(plan.descriptor, ordered=False),
        "descriptor": plan.descriptor,
    }
