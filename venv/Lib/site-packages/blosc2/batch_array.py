#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import copy
import pathlib
import statistics
from collections.abc import Iterator, Sequence
from dataclasses import asdict
from functools import lru_cache
from typing import Any

import numpy as np

import blosc2
from blosc2.info import InfoReporter, format_nbytes_info
from blosc2.msgpack_utils import msgpack_packb, msgpack_unpackb

_BATCHARRAY_META = {"version": 1, "serializer": "msgpack", "items_per_block": None, "arrow_schema": None}
_SUPPORTED_SERIALIZERS = {"msgpack", "arrow"}
_BATCHARRAY_VLMETA_KEY = "_batch_array_metadata"


def _check_serialized_size(buffer: bytes) -> None:
    if len(buffer) > blosc2.MAX_BUFFERSIZE:
        raise ValueError(f"Serialized objects cannot be larger than {blosc2.MAX_BUFFERSIZE} bytes")


class Batch(Sequence[Any]):
    """A lazy sequence representing one batch in a :class:`BatchArray`.

    ``Batch`` provides sequence-style access to the items stored in a single
    batch. Integer indexing can use block-local reads when possible, while
    slicing materializes the full batch into Python items.

    Batch instances are normally obtained via :class:`BatchArray` indexing or
    iteration rather than constructed directly.
    """

    def __init__(self, parent: BatchArray, nbatch: int, lazybatch: bytes) -> None:
        self._parent = parent
        self._nbatch = nbatch
        self._lazybatch = lazybatch
        self._items: list[Any] | None = None
        self._cached_block_index: int | None = None
        self._cached_block: list[Any] | None = None
        self._cached_block_column_index: int | None = None
        self._cached_block_column = None
        self._nbytes, self._cbytes, self._nblocks = blosc2.get_cbuffer_sizes(lazybatch)

    def _normalize_index(self, index: int) -> int:
        if not isinstance(index, int):
            raise TypeError("Batch indices must be integers")
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError("Batch index out of range")
        return index

    def _decode_items(self) -> list[Any]:
        if self._items is None:
            blocks = self._parent._decode_blocks(self._nbatch)
            self._items = [item for block in blocks for item in block]
        return self._items

    def _get_block(self, block_index: int) -> list[Any]:
        if self._cached_block_index == block_index and self._cached_block is not None:
            return self._cached_block
        block = self._parent._deserialize_block(self._parent.schunk.get_vlblock(self._nbatch, block_index))
        self._cached_block_index = block_index
        self._cached_block = block
        return block

    def _get_block_item(self, block_index: int, item_index: int) -> Any:
        if self._cached_block_index == block_index and self._cached_block is not None:
            return self._cached_block[item_index]
        if self._parent._serializer != "arrow":
            return self._get_block(block_index)[item_index]
        if self._cached_block_column_index != block_index or self._cached_block_column is None:
            payload = self._parent.schunk.get_vlblock(self._nbatch, block_index)
            self._cached_block_column = self._parent._deserialize_arrow_block_column(payload)
            self._cached_block_column_index = block_index
        return self._cached_block_column[item_index].as_py()

    def __getitem__(self, index: int | slice) -> Any | list[Any]:
        if isinstance(index, slice):
            items = self._decode_items()
            return items[index]
        if index < 0:
            items = self._decode_items()
            index = self._normalize_index(index)
            return items[index]
        items_per_block = self._parent.items_per_block
        if items_per_block is not None:
            block_index, item_index = divmod(index, items_per_block)
            if block_index >= self._nblocks:
                raise IndexError("Batch index out of range")
            try:
                return self._get_block_item(block_index, item_index)
            except IndexError as exc:
                raise IndexError("Batch index out of range") from exc
        items = self._decode_items()
        index = self._normalize_index(index)
        return items[index]

    def __len__(self) -> int:
        batch_length = self._parent._batch_length(self._nbatch)
        if batch_length is not None:
            return batch_length
        return len(self._decode_items())

    def __iter__(self) -> Iterator[Any]:
        for i in range(len(self)):
            yield self[i]

    @property
    def lazybatch(self) -> bytes:
        return self._lazybatch

    @property
    def nbytes(self) -> int:
        return self._nbytes

    @property
    def cbytes(self) -> int:
        return self._cbytes

    @property
    def cratio(self) -> float:
        return self._nbytes / self._cbytes

    def __repr__(self) -> str:
        return f"Batch(len={len(self)}, nbytes={self.nbytes}, cbytes={self.cbytes})"


class BatchArrayItems(Sequence[Any]):
    """A read-only flat view over the items stored in a :class:`BatchArray`."""

    def __init__(self, parent: BatchArray) -> None:
        self._parent = parent

    def __getitem__(self, index: int | slice) -> Any | list[Any]:
        return self._parent._get_flat_item(index)

    def __len__(self) -> int:
        return self._parent._get_total_item_count()


class BatchArray:
    """A batched container for variable-length Python items.

    BatchArray stores data as a sequence of *batches*, where each batch contains
    one or more Python items. Each batch is stored in one compressed chunk, and
    each chunk is internally split into one or more variable-length blocks for
    efficient item access.

    The main abstraction is batch-oriented:

    - indexing the store returns batches
    - iterating the store yields batches
    - :meth:`iter_items` provides flat item-wise traversal

    BatchArray is a good fit when:

    - data arrives naturally in batches
    - batch-level append/update operations are important
    - occasional item-level reads are needed inside a batch

    Parameters
    ----------
    items_per_block : int, optional
        Maximum number of items stored in each internal variable-length block.
        The last block in a batch may contain fewer items than this cap. If not
        provided, a value is inferred from the first batch using the serialized
        item sizes and the compression level. The heuristic uses fixed byte
        budgets so that the layout (and hence the compression ratio) does not
        depend on the CPU it was created on: 1 MiB for ``clevel`` 1 through 3,
        8 MiB for ``clevel`` 4 through 6, and 16 MiB for ``clevel`` 7 and 8.
        At ``clevel`` 9, the whole batch is kept as one block. Smaller blocks
        generally improve random access, while larger blocks generally improve
        compression ratio.
    serializer : {"msgpack", "arrow"}, optional
        Serializer used for batch payloads. ``"msgpack"`` is the default and is
        the general-purpose choice for Python items, including nested Blosc2
        containers such as :class:`blosc2.NDArray`, :class:`blosc2.SChunk`,
        :class:`blosc2.ObjectArray`, :class:`blosc2.BatchArray`, and
        :class:`blosc2.EmbedStore`, which are serialized transparently via
        :meth:`to_cframe` / :func:`blosc2.from_cframe`. Msgpack also supports
        structured Blosc2 reference objects, currently
        :class:`blosc2.C2Array`, :class:`blosc2.LazyExpr`, and
        :class:`blosc2.LazyUDF` backed by :func:`blosc2.dsl_kernel`. These lazy
        objects preserve reference semantics, so only persistent local
        operands, :class:`blosc2.C2Array` operands, and
        :class:`blosc2.DictStore` members are supported; purely in-memory
        operands are rejected. Plain Python :class:`blosc2.LazyUDF` callables
        are not serialized by msgpack. ``"arrow"`` is optional and requires
        ``pyarrow``.
    _from_schunk : blosc2.SChunk, optional
        Internal hook used when reopening an already-tagged BatchArray.
    **kwargs
        Storage, compression, and decompression arguments accepted by the
        constructor.
    """

    @staticmethod
    def _set_typesize_one(cparams: blosc2.CParams | dict | None) -> blosc2.CParams | dict:
        if cparams is None:
            cparams = blosc2.CParams()
        elif isinstance(cparams, blosc2.CParams):
            cparams = copy.deepcopy(cparams)
        else:
            cparams = dict(cparams)

        if isinstance(cparams, blosc2.CParams):
            cparams.typesize = 1
        else:
            cparams["typesize"] = 1
        return cparams

    @staticmethod
    def _coerce_storage(storage: blosc2.Storage | dict | None, kwargs: dict[str, Any]) -> blosc2.Storage:
        if storage is not None:
            storage_keys = set(blosc2.Storage.__annotations__)
            storage_kwargs = storage_keys.intersection(kwargs)
            if storage_kwargs:
                unexpected = ", ".join(sorted(storage_kwargs))
                raise AttributeError(
                    f"Cannot pass both `storage` and other kwargs already included in Storage: {unexpected}"
                )
            if isinstance(storage, blosc2.Storage):
                return copy.deepcopy(storage)
            return blosc2.Storage(**storage)

        storage_kwargs = {
            name: kwargs.pop(name) for name in list(blosc2.Storage.__annotations__) if name in kwargs
        }
        return blosc2.Storage(**storage_kwargs)

    @staticmethod
    def _validate_storage(storage: blosc2.Storage) -> None:
        if storage.mmap_mode not in (None, "r"):
            raise ValueError("For BatchArray containers, mmap_mode must be None or 'r'")
        if storage.mmap_mode == "r" and storage.mode != "r":
            raise ValueError("For BatchArray containers, mmap_mode='r' requires mode='r'")

    def _attach_schunk(self, schunk: blosc2.SChunk) -> None:
        self.schunk = schunk
        self.mode = schunk.mode
        self.mmap_mode = getattr(schunk, "mmap_mode", None)
        try:
            batcharray_meta = self.schunk.meta["batcharray"]
        except KeyError:
            batcharray_meta = {}
        self._serializer = batcharray_meta.get("serializer", self._serializer)
        self._items_per_block = batcharray_meta.get("items_per_block", self._items_per_block)
        self._arrow_schema = batcharray_meta.get("arrow_schema", self._arrow_schema)
        self._arrow_schema_obj = None
        self._batch_lengths = self._load_batch_lengths()
        self._items = BatchArrayItems(self)
        self._item_prefix_sums: np.ndarray | None = None
        self._validate_tag()

    def _maybe_open_existing(self, storage: blosc2.Storage) -> bool:
        urlpath = storage.urlpath
        if urlpath is None or storage.mode not in ("r", "a") or not pathlib.Path(urlpath).exists():
            return False

        schunk = blosc2.blosc2_ext.open(urlpath, mode=storage.mode, offset=0, mmap_mode=storage.mmap_mode)
        self._attach_schunk(schunk)
        return True

    def _make_storage(self) -> blosc2.Storage:
        meta = {name: self.meta[name] for name in self.meta}
        return blosc2.Storage(
            contiguous=self.schunk.contiguous,
            urlpath=self.urlpath,
            mode=self.mode,
            mmap_mode=self.mmap_mode,
            meta=meta,
        )

    def __init__(
        self,
        items_per_block: int | None = None,
        serializer: str = "msgpack",
        _from_schunk: blosc2.SChunk | None = None,
        **kwargs: Any,
    ) -> None:
        """Create a new BatchArray or reopen an existing one.

        When a persistent ``urlpath`` points to an existing BatchArray and the
        mode is ``"r"`` or ``"a"``, the container is reopened automatically.
        Otherwise a new empty store is created.
        """
        if items_per_block is not None and items_per_block <= 0:
            raise ValueError("items_per_block must be a positive integer")
        if serializer not in _SUPPORTED_SERIALIZERS:
            raise ValueError(f"Unsupported BatchArray serializer: {serializer!r}")
        self._items_per_block: int | None = items_per_block
        self._serializer = serializer
        self._arrow_schema: bytes | None = None
        self._arrow_schema_obj = None
        self._batch_lengths: list[int] | None = None
        if _from_schunk is not None:
            if kwargs:
                unexpected = ", ".join(sorted(kwargs))
                raise ValueError(f"Cannot pass {unexpected} together with `_from_schunk`")
            self._attach_schunk(_from_schunk)
            return
        cparams = kwargs.pop("cparams", None)
        dparams = kwargs.pop("dparams", None)
        storage = kwargs.pop("storage", None)
        storage = self._coerce_storage(storage, kwargs)

        if kwargs:
            unexpected = ", ".join(sorted(kwargs))
            raise ValueError(f"Unsupported BatchArray keyword argument(s): {unexpected}")

        self._validate_storage(storage)
        cparams = self._set_typesize_one(cparams)

        if dparams is None:
            dparams = blosc2.DParams()

        if self._maybe_open_existing(storage):
            return

        fixed_meta = dict(storage.meta or {})
        fixed_meta["batcharray"] = {
            **_BATCHARRAY_META,
            "serializer": self._serializer,
            "items_per_block": self._items_per_block,
            "arrow_schema": self._arrow_schema,
        }
        storage.meta = fixed_meta
        schunk = blosc2.SChunk(chunksize=-1, data=None, cparams=cparams, dparams=dparams, storage=storage)
        self._attach_schunk(schunk)

    def _validate_tag(self) -> None:
        if "batcharray" not in self.schunk.meta:
            raise ValueError("The supplied SChunk is not tagged as a BatchArray")
        if self._serializer not in _SUPPORTED_SERIALIZERS:
            raise ValueError(f"Unsupported BatchArray serializer in metadata: {self._serializer!r}")
        if self._serializer == "arrow":
            self._require_pyarrow()

    @staticmethod
    @lru_cache(maxsize=1)
    def _require_pyarrow():
        try:
            import pyarrow as pa
            import pyarrow.ipc as pa_ipc
        except ImportError as exc:
            raise ImportError("BatchArray serializer='arrow' requires pyarrow") from exc
        return pa, pa_ipc

    def _check_writable(self) -> None:
        if self.mode == "r":
            raise ValueError("Cannot modify a BatchArray opened in read-only mode")

    def _normalize_index(self, index: int) -> int:
        if not isinstance(index, int):
            raise TypeError("BatchArray indices must be integers")
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError("BatchArray index out of range")
        return index

    def _normalize_insert_index(self, index: int) -> int:
        if not isinstance(index, int):
            raise TypeError("BatchArray indices must be integers")
        if index < 0:
            index += len(self)
            if index < 0:
                return 0
        if index > len(self):
            return len(self)
        return index

    def _slice_indices(self, index: slice) -> list[int]:
        return list(range(*index.indices(len(self))))

    def _copy_meta(self) -> dict[str, Any]:
        return {name: self.meta[name] for name in self.meta}

    def _load_batch_lengths(self) -> list[int] | None:
        try:
            metadata = self.schunk.vlmeta[_BATCHARRAY_VLMETA_KEY]
        except KeyError:
            return None
        batch_lengths = metadata.get("batch_lengths")
        if not isinstance(batch_lengths, list):
            return None
        return [int(length) for length in batch_lengths]

    def _persist_batch_lengths(self) -> None:
        if self._batch_lengths is None:
            return
        if len(self._batch_lengths) == 0:
            if _BATCHARRAY_VLMETA_KEY in self.vlmeta:
                del self.vlmeta[_BATCHARRAY_VLMETA_KEY]
            return
        self.schunk.vlmeta[_BATCHARRAY_VLMETA_KEY] = {"batch_lengths": list(self._batch_lengths)}

    def _get_batch_lengths(self) -> list[int] | None:
        return self._batch_lengths

    def _ensure_batch_lengths(self) -> list[int]:
        if self._batch_lengths is None:
            self._batch_lengths = []
        return self._batch_lengths

    def _load_or_compute_batch_lengths(self) -> list[int]:
        if self._batch_lengths is None:
            self._batch_lengths = [len(self._get_batch(i)) for i in range(len(self))]
            if self.mode != "r":
                self._persist_batch_lengths()
        return self._batch_lengths

    def _batch_length(self, index: int) -> int | None:
        if self._batch_lengths is None:
            return None
        return self._batch_lengths[index]

    def _invalidate_item_cache(self) -> None:
        self._item_prefix_sums = None

    def _get_item_prefix_sums(self) -> np.ndarray:
        if self._item_prefix_sums is None:
            batch_lengths = np.asarray(self._load_or_compute_batch_lengths(), dtype=np.int64)
            prefix_sums = np.empty(len(batch_lengths) + 1, dtype=np.int64)
            prefix_sums[0] = 0
            prefix_sums[1:] = np.cumsum(batch_lengths, dtype=np.int64)
            self._item_prefix_sums = prefix_sums
        return self._item_prefix_sums

    def _get_total_item_count(self) -> int:
        return int(self._get_item_prefix_sums()[-1])

    def _get_flat_item(self, index: int | slice) -> Any | list[Any]:
        if isinstance(index, slice):
            return [self._get_flat_item(i) for i in range(*index.indices(self._get_total_item_count()))]
        if not isinstance(index, int):
            raise TypeError("BatchArray item indices must be integers")
        nitems = self._get_total_item_count()
        if index < 0:
            index += nitems
        if index < 0 or index >= nitems:
            raise IndexError("BatchArray item index out of range")

        prefix_sums = self._get_item_prefix_sums()
        batch_index = int(np.searchsorted(prefix_sums, index, side="right") - 1)
        item_index = int(index - prefix_sums[batch_index])
        return self[batch_index][item_index]

    def _block_sizes_from_batch_length(self, batch_length: int, nblocks: int) -> list[int]:
        if self._items_per_block is None or nblocks <= 0:
            return []
        full_blocks, remainder = divmod(batch_length, self._items_per_block)
        block_sizes = [self._items_per_block] * full_blocks
        if remainder:
            block_sizes.append(remainder)
        if not block_sizes and batch_length > 0:
            block_sizes.append(batch_length)
        if len(block_sizes) != nblocks:
            return []
        return block_sizes

    def _get_block_sizes(self, batch_sizes: list[int]) -> list[int] | None:
        if self._items_per_block is None:
            return None
        block_sizes: list[int] = []
        for index, batch_length in enumerate(batch_sizes):
            lazychunk = self.schunk.get_lazychunk(index)
            _, _, nblocks = blosc2.get_cbuffer_sizes(lazychunk)
            sizes = self._block_sizes_from_batch_length(batch_length, nblocks)
            if not sizes:
                return None
            block_sizes.extend(sizes)
        return block_sizes

    def _total_nblocks(self) -> int:
        total = 0
        for index in range(len(self)):
            lazychunk = self.schunk.get_lazychunk(index)
            _, _, nblocks = blosc2.get_cbuffer_sizes(lazychunk)
            total += nblocks
        return total

    def _user_vlmeta_items(self) -> dict[str, Any]:
        return {key: value for key, value in self.vlmeta.getall().items() if key != _BATCHARRAY_VLMETA_KEY}

    def _normalize_msgpack_batch(self, value: object) -> list[Any]:
        if isinstance(value, (str, bytes, bytearray, memoryview)):
            raise TypeError("BatchArray entries must be sequences of Python objects")
        if not isinstance(value, Sequence):
            raise TypeError("BatchArray entries must be sequences of Python objects")
        values = list(value)
        if len(values) == 0:
            raise ValueError("BatchArray entries cannot be empty")
        return values

    def _normalize_arrow_batch(self, value: object):
        pa, _ = self._require_pyarrow()
        if isinstance(value, pa.ChunkedArray):
            value = value.combine_chunks()
        elif isinstance(value, pa.RecordBatch):
            if value.num_columns != 1:
                raise TypeError("Arrow RecordBatch inputs for BatchArray must have exactly one column")
            value = value.column(0)
        elif not isinstance(value, pa.Array):
            if isinstance(value, (str, bytes, bytearray, memoryview)):
                raise TypeError("BatchArray entries must be Arrow arrays or sequences of Python objects")
            if not isinstance(value, Sequence):
                raise TypeError("BatchArray entries must be Arrow arrays or sequences of Python objects")
            value = pa.array(list(value))
        if len(value) == 0:
            raise ValueError("BatchArray entries cannot be empty")
        self._ensure_arrow_schema(value)
        return value

    def _ensure_arrow_schema(self, batch) -> None:
        if self._serializer != "arrow":
            return
        pa, _ = self._require_pyarrow()
        schema = pa.schema([pa.field("values", batch.type)])
        if self._arrow_schema is None:
            self._arrow_schema = schema.serialize().to_pybytes()
            self._arrow_schema_obj = schema
            return
        existing_schema = self._get_arrow_schema()
        if not existing_schema.equals(schema):
            raise TypeError("All Arrow batches in a BatchArray must share the same schema")

    def _get_arrow_schema(self):
        if self._serializer != "arrow":
            return None
        if self._arrow_schema is None:
            raise RuntimeError("Arrow schema is not initialized")
        if self._arrow_schema_obj is None:
            pa, pa_ipc = self._require_pyarrow()
            self._arrow_schema_obj = pa_ipc.read_schema(pa.BufferReader(self._arrow_schema))
        return self._arrow_schema_obj

    def _normalize_batch(self, value: object) -> Any:
        if self._serializer == "arrow":
            return self._normalize_arrow_batch(value)
        return self._normalize_msgpack_batch(value)

    def _batch_len(self, batch: Any) -> int:
        return len(batch)

    def _payload_sizes_for_batch(self, batch: Any) -> list[int]:
        if self._serializer == "arrow":
            total_size = batch.get_total_buffer_size()
            avg_size = max(1, total_size // max(1, len(batch)))
            return [avg_size] * len(batch)
        return [len(msgpack_packb(item)) for item in batch]

    def _ensure_layout_for_batch(self, batch: Any) -> None:
        layout_changed = False
        if self._items_per_block is None:
            payload_sizes = self._payload_sizes_for_batch(batch)
            self._items_per_block = self._guess_blocksize(payload_sizes)
            layout_changed = True
        if self._serializer == "arrow" and self._arrow_schema is not None:
            layout_changed = layout_changed or len(self) == 0
        if layout_changed:
            self._persist_layout_metadata()

    def _persist_layout_metadata(self) -> None:
        if len(self) > 0:
            return
        batch_lengths = None if self._batch_lengths is None else list(self._batch_lengths)
        user_vlmeta = self._user_vlmeta_items() if len(self.vlmeta) > 0 else {}
        storage = self._make_storage()
        fixed_meta = dict(storage.meta or {})
        fixed_meta["batcharray"] = {
            **dict(fixed_meta.get("batcharray", {})),
            "items_per_block": self._items_per_block,
            "serializer": self._serializer,
            "arrow_schema": self._arrow_schema,
        }
        storage.meta = fixed_meta
        schunk = blosc2.SChunk(
            chunksize=-1,
            data=None,
            cparams=copy.deepcopy(self.cparams),
            dparams=copy.deepcopy(self.dparams),
            storage=storage,
        )
        self._attach_schunk(schunk)
        for key, value in user_vlmeta.items():
            self.vlmeta[key] = value
        if batch_lengths is not None and self._batch_lengths is None:
            self._batch_lengths = batch_lengths

    def _guess_blocksize(self, payload_sizes: list[int]) -> int:
        if not payload_sizes:
            raise ValueError("BatchArray entries cannot be empty")
        clevel = self.cparams.clevel
        # For serialized batch payloads, especially Arrow IPC, L1-sized blocks are often
        # too small for codecs like Zstd to exploit cross-row redundancy.  Use larger
        # cache-budget tiers as clevel increases, while avoiding full L2 blocks at the
        # default clevel to keep random access reasonably granular.
        # UPDATE: to avoid cratio differences between CPUs, better use fixed budgets instead
        # of CPU cache sizes.
        if clevel == 9:
            return len(payload_sizes)
        if 0 < clevel <= 3:
            # budget = blosc2.cpu_info.get("l1_data_cache_size")
            budget = 2**20  #  1 MB
        elif 3 < clevel <= 6:
            # budget = blosc2.cpu_info.get("l2_cache_size") // 2
            budget = 2**23  #  8 MB
        elif 6 < clevel < 9:
            # budget = blosc2.cpu_info.get("l2_cache_size")
            budget = 2**24  # 16 MB
        else:
            return len(payload_sizes)
        if not isinstance(budget, int) or budget <= 0:
            return len(payload_sizes)
        total = 0
        count = 0
        for payload_size in payload_sizes:
            if count > 0 and total + payload_size > budget:
                break
            total += payload_size
            count += 1
        if count == 0:
            count = 1
        return min(count, len(payload_sizes))

    def _serialize_batch(self, value: object) -> Any:
        batch = self._normalize_batch(value)
        self._ensure_layout_for_batch(batch)
        return batch

    def _serialize_msgpack_block(self, items: list[Any]) -> bytes:
        payload = msgpack_packb(items)
        _check_serialized_size(payload)
        return payload

    def _serialize_arrow_block(self, items) -> bytes:
        pa, pa_ipc = self._require_pyarrow()
        batch = pa.record_batch([items], schema=self._get_arrow_schema())
        sink = pa.BufferOutputStream()
        with pa_ipc.new_stream(sink, batch.schema) as writer:
            writer.write_batch(batch)
        payload = sink.getvalue().to_pybytes()
        _check_serialized_size(payload)
        return payload

    def _serialize_block(self, items: Any) -> bytes:
        if self._serializer == "arrow":
            return self._serialize_arrow_block(items)
        return self._serialize_msgpack_block(items)

    def _deserialize_msgpack_block(self, payload: bytes) -> list[Any]:
        return msgpack_unpackb(payload)

    def _deserialize_arrow_block_column(self, payload: bytes):
        pa, pa_ipc = self._require_pyarrow()
        try:
            reader = pa_ipc.open_stream(pa.BufferReader(payload))
            batch = reader.read_next_batch()
        except (pa.ArrowInvalid, OSError):
            # Backward compatibility for older arrow-serializer blocks written
            # as bare serialized RecordBatch payloads.  Those cannot represent
            # dictionary batches reliably, so new blocks use IPC streams.
            batch = pa_ipc.read_record_batch(pa.BufferReader(payload), self._get_arrow_schema())
        return batch.column(0)

    def _deserialize_arrow_block(self, payload: bytes) -> list[Any]:
        return self._deserialize_arrow_block_column(payload).to_pylist()

    def _deserialize_block(self, payload: bytes) -> list[Any]:
        if self._serializer == "arrow":
            return self._deserialize_arrow_block(payload)
        return self._deserialize_msgpack_block(payload)

    def _deserialize_arrow_block_item(self, payload: bytes, item_index: int) -> Any:
        return self._deserialize_arrow_block_column(payload)[item_index].as_py()

    def _deserialize_block_item(self, payload: bytes, item_index: int) -> Any:
        if self._serializer == "arrow":
            return self._deserialize_arrow_block_item(payload, item_index)
        return self._deserialize_msgpack_block(payload)[item_index]

    def _vl_cparams_kwargs(self) -> dict[str, Any]:
        return asdict(self.schunk.cparams)

    def _vl_dparams_kwargs(self) -> dict[str, Any]:
        return asdict(self.schunk.dparams)

    def _compress_batch(self, batch: Any) -> bytes:
        if self._items_per_block is None:
            raise RuntimeError("BatchArray items_per_block is not initialized")
        blocks = [
            self._serialize_block(batch[i : i + self._items_per_block])
            for i in range(0, self._batch_len(batch), self._items_per_block)
        ]
        return blosc2.blosc2_ext.vlcompress(blocks, **self._vl_cparams_kwargs())

    def _decode_blocks(self, nbatch: int) -> list[list[Any]]:
        block_payloads = blosc2.blosc2_ext.vldecompress(
            self.schunk.get_chunk(nbatch), **self._vl_dparams_kwargs()
        )
        return [self._deserialize_block(payload) for payload in block_payloads]

    def _get_batch(self, index: int) -> Batch:
        return Batch(self, index, self.schunk.get_lazychunk(index))

    def append(self, value: object) -> int:
        """Append one batch and return the new number of batches."""
        self._check_writable()
        batch = self._serialize_batch(value)
        batch_payload = self._compress_batch(batch)
        length = self._batch_len(batch)
        new_len = self.schunk.append_chunk(batch_payload)
        self._ensure_batch_lengths().append(length)
        self._persist_batch_lengths()
        self._invalidate_item_cache()
        return new_len

    def insert(self, index: int, value: object) -> int:
        """Insert one batch at ``index`` and return the new number of batches."""
        self._check_writable()
        index = self._normalize_insert_index(index)
        batch = self._serialize_batch(value)
        batch_payload = self._compress_batch(batch)
        length = self._batch_len(batch)
        new_len = self.schunk.insert_chunk(index, batch_payload)
        self._ensure_batch_lengths().insert(index, length)
        self._persist_batch_lengths()
        self._invalidate_item_cache()
        return new_len

    def delete(self, index: int | slice) -> int:
        """Delete the batch at ``index`` and return the new number of batches."""
        self._check_writable()
        if isinstance(index, slice):
            # Delete in descending order so earlier deletions don't shift
            # the indices of chunks yet to be deleted (negative-step slices
            # produce ascending indices when merely reversed).
            for idx in sorted(self._slice_indices(index), reverse=True):
                self.schunk.delete_chunk(idx)
                if self._batch_lengths is not None:
                    del self._batch_lengths[idx]
            self._persist_batch_lengths()
            self._invalidate_item_cache()
            return len(self)
        index = self._normalize_index(index)
        new_len = self.schunk.delete_chunk(index)
        if self._batch_lengths is not None:
            del self._batch_lengths[index]
            self._persist_batch_lengths()
        self._invalidate_item_cache()
        return new_len

    def pop(self, index: int = -1) -> list[Any]:
        """Remove and return the batch at ``index`` as a Python list."""
        self._check_writable()
        if isinstance(index, slice):
            raise NotImplementedError("Slicing is not supported for BatchArray")
        index = self._normalize_index(index)
        value = self[index][:]
        self.delete(index)
        return value

    def extend(self, values: object) -> None:
        """Append all batches from an iterable of batches."""
        self._check_writable()
        for value in values:
            batch = self._serialize_batch(value)
            batch_payload = self._compress_batch(batch)
            self.schunk.append_chunk(batch_payload)
            self._ensure_batch_lengths().append(self._batch_len(batch))
        self._persist_batch_lengths()
        self._invalidate_item_cache()

    def clear(self) -> None:
        """Remove all entries from the container."""
        self._check_writable()
        storage = self._make_storage()
        if storage.urlpath is not None:
            blosc2.remove_urlpath(storage.urlpath)
        schunk = blosc2.SChunk(
            chunksize=-1,
            data=None,
            cparams=copy.deepcopy(self.cparams),
            dparams=copy.deepcopy(self.dparams),
            storage=storage,
        )
        self._attach_schunk(schunk)
        self._batch_lengths = []
        self._persist_batch_lengths()
        self._invalidate_item_cache()

    def __getitem__(self, index: int | slice) -> Batch | list[Batch]:
        """Return one batch or a list of batches."""
        if isinstance(index, slice):
            return [self[i] for i in self._slice_indices(index)]
        index = self._normalize_index(index)
        return self._get_batch(index)

    def __setitem__(self, index: int | slice, value: object) -> None:
        if isinstance(index, slice):
            self._check_writable()
            indices = self._slice_indices(index)
            values = list(value)
            step = 1 if index.step is None else index.step
            if step == 1:
                start = self._normalize_insert_index(0 if index.start is None else index.start)
                for idx in reversed(indices):
                    self.schunk.delete_chunk(idx)
                    if self._batch_lengths is not None:
                        del self._batch_lengths[idx]
                for offset, item in enumerate(values):
                    batch = self._serialize_batch(item)
                    batch_payload = self._compress_batch(batch)
                    self.schunk.insert_chunk(start + offset, batch_payload)
                    self._ensure_batch_lengths().insert(start + offset, self._batch_len(batch))
                self._persist_batch_lengths()
                self._invalidate_item_cache()
                return
            if len(values) != len(indices):
                raise ValueError(
                    f"attempt to assign sequence of size {len(values)} to extended slice of size {len(indices)}"
                )
            for idx, item in zip(indices, values, strict=True):
                batch = self._serialize_batch(item)
                batch_payload = self._compress_batch(batch)
                self.schunk.update_chunk(idx, batch_payload)
                if self._batch_lengths is not None:
                    self._batch_lengths[idx] = self._batch_len(batch)
            self._persist_batch_lengths()
            self._invalidate_item_cache()
            return
        self._check_writable()
        index = self._normalize_index(index)
        batch = self._serialize_batch(value)
        batch_payload = self._compress_batch(batch)
        self.schunk.update_chunk(index, batch_payload)
        if self._batch_lengths is not None:
            self._batch_lengths[index] = self._batch_len(batch)
            self._persist_batch_lengths()
        self._invalidate_item_cache()

    def __delitem__(self, index: int | slice) -> None:
        self.delete(index)

    def __len__(self) -> int:
        """Return the number of batches stored in the container."""
        return self.schunk.nchunks

    def iter_items(self) -> Iterator[Any]:
        """Iterate over all items across all batches in order."""
        for batch in self:
            yield from batch

    def __iter__(self) -> Iterator[Batch]:
        for i in range(len(self)):
            yield self[i]

    @property
    def meta(self):
        return self.schunk.meta

    @property
    def vlmeta(self):
        return self.schunk.vlmeta

    @property
    def cparams(self):
        return self.schunk.cparams

    @property
    def dparams(self):
        return self.schunk.dparams

    @property
    def items_per_block(self) -> int | None:
        """Maximum number of items per internal block.

        The last block in a batch may contain fewer items.
        """
        return self._items_per_block

    @property
    def items(self) -> BatchArrayItems:
        return self._items

    @property
    def typesize(self) -> int:
        return self.schunk.typesize

    @property
    def nbytes(self) -> int:
        return self.schunk.nbytes

    @property
    def cbytes(self) -> int:
        return self.schunk.cbytes

    @property
    def cratio(self) -> float:
        return self.schunk.cratio

    @property
    def urlpath(self) -> str | None:
        return self.schunk.urlpath

    @property
    def contiguous(self) -> bool:
        return self.schunk.contiguous

    @property
    def info(self) -> InfoReporter:
        """Return an info reporter with a compact summary of the store."""
        return InfoReporter(self)

    @property
    def info_items(self) -> list:
        """Return summary information as ``(name, value)`` pairs."""
        batch_sizes = self._get_batch_lengths()
        if batch_sizes is None:
            batch_sizes = [len(batch) for batch in self]
        block_sizes = self._get_block_sizes(batch_sizes)
        if batch_sizes:
            batch_stats = (
                f"mean={statistics.fmean(batch_sizes):.2f}, max={max(batch_sizes)}, min={min(batch_sizes)}"
            )
            nbatches_value = f"{len(self)} (items per batch: {batch_stats})"
        else:
            nbatches_value = f"{len(self)} (items per batch: n/a)"
        if block_sizes:
            block_stats = (
                f"mean={statistics.fmean(block_sizes):.2f}, max={max(block_sizes)}, min={min(block_sizes)}"
            )
            nblocks_value = f"{self._total_nblocks()} (items per block: {block_stats})"
        else:
            nblocks_value = f"{self._total_nblocks()} (items per block: n/a)"
        return [
            ("type", f"{self.__class__.__name__}"),
            ("serializer", self.serializer),
            ("items_per_block", self.items_per_block),
            ("nbatches", nbatches_value),
            ("nblocks", nblocks_value),
            ("nitems", sum(batch_sizes)),
            ("nbytes", format_nbytes_info(self.nbytes)),
            ("cbytes", format_nbytes_info(self.cbytes)),
            ("cratio", f"{self.cratio:.2f}x"),
            ("cparams", self.cparams),
            ("dparams", self.dparams),
        ]

    def to_cframe(self) -> bytes:
        """Serialize the full store to a Blosc2 cframe buffer."""
        return self.schunk.to_cframe()

    def chunk_copy(self, **kwargs: Any) -> BatchArray:
        """Create a copy by transferring compressed chunks directly at the C level.

        This is significantly faster than :meth:`copy` because it bypasses all
        Python-level serialisation/deserialisation: each blosc2 chunk is read
        from the source SChunk and appended to the destination SChunk as raw
        compressed bytes.

        The destination is created with the **same** ``cparams`` as the source
        so existing chunks are accepted without recompression.  Passing a
        ``cparams`` override is not allowed (raises :class:`ValueError`); use
        :meth:`copy` instead if you need to change the compression settings.

        Parameters
        ----------
        **kwargs:
            Forwarded to the :class:`BatchArray` constructor.  Typical use
            cases: ``urlpath`` / ``mode`` (persistent copy), ``contiguous``,
            ``dparams``.  Do **not** pass ``cparams``, ``meta``, ``serializer``
            or ``items_per_block``.

        Returns
        -------
        BatchArray
            A new standalone copy with identical data and storage metadata.

        Raises
        ------
        ValueError
            If ``cparams`` is in *kwargs* (recompression is not supported by
            this method).

        See Also
        --------
        copy : Element-wise copy that supports cparams overrides.
        """
        if "cparams" in kwargs:
            raise ValueError(
                "chunk_copy() does not support a cparams override because it transfers "
                "pre-compressed chunks as-is.  Use copy() if you need to change cparams."
            )
        if "meta" in kwargs:
            raise ValueError("meta should not be passed to chunk_copy")
        kwargs["cparams"] = copy.deepcopy(self.cparams)
        kwargs.setdefault("dparams", copy.deepcopy(self.dparams))
        kwargs.setdefault("items_per_block", self.items_per_block)
        kwargs.setdefault("serializer", self.serializer)
        kwargs.setdefault("contiguous", self.schunk.contiguous)
        if "urlpath" in kwargs and "mode" not in kwargs:
            kwargs["mode"] = "w"
        kwargs["meta"] = self._copy_meta()

        out = BatchArray(**kwargs)

        src_sc = self.schunk
        dst_sc = out.schunk
        for i in range(src_sc.nchunks):
            dst_sc.append_chunk(src_sc.get_chunk(i))

        # Persist batch_lengths so reopening the file skips the recompute scan.
        out._batch_lengths = list(self._load_or_compute_batch_lengths()) if src_sc.nchunks > 0 else []
        out._persist_batch_lengths()
        out._invalidate_item_cache()

        # Preserve any user-defined vlmeta items (batch-lengths key is internal).
        for key, value in self._user_vlmeta_items().items():
            out.vlmeta[key] = value

        return out

    def copy(self, **kwargs: Any) -> BatchArray:
        """Create a copy of the store with optional constructor overrides."""
        if "meta" in kwargs:
            raise ValueError("meta should not be passed to copy")
        kwargs["cparams"] = kwargs.get("cparams", copy.deepcopy(self.cparams))
        kwargs["dparams"] = kwargs.get("dparams", copy.deepcopy(self.dparams))
        kwargs["items_per_block"] = kwargs.get("items_per_block", self.items_per_block)
        kwargs["serializer"] = kwargs.get("serializer", self.serializer)
        user_vlmeta = self._user_vlmeta_items() if len(self.vlmeta) > 0 else {}

        if "storage" in kwargs:
            storage = self._coerce_storage(kwargs["storage"], {})
            fixed_meta = self._copy_meta()
            if storage.meta is not None:
                fixed_meta.update(storage.meta)
            storage.meta = fixed_meta
            kwargs["storage"] = storage
        else:
            kwargs["meta"] = self._copy_meta()
            kwargs["contiguous"] = kwargs.get("contiguous", self.schunk.contiguous)
            if "urlpath" in kwargs and "mode" not in kwargs:
                kwargs["mode"] = "w"

        out = BatchArray(**kwargs)
        for key, value in user_vlmeta.items():
            out.vlmeta[key] = value
        out.extend(self)
        return out

    def __enter__(self) -> BatchArray:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        return False

    def __repr__(self) -> str:
        return f"BatchArray(len={len(self)}, urlpath={self.urlpath!r})"

    @property
    def serializer(self) -> str:
        """Serializer name used for batch payloads."""
        return self._serializer
