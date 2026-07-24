#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import copy
import pathlib
from typing import TYPE_CHECKING, Any

import blosc2
from blosc2.info import InfoReporter, format_nbytes_info
from blosc2.msgpack_utils import msgpack_packb, msgpack_unpackb

if TYPE_CHECKING:
    from collections.abc import Iterator

    from blosc2.schunk import SChunk

# On-disk metadata tag is kept as "vlarray" for backward compatibility with
# existing stored files.  The public class name is ObjectArray.
_VLARRAY_META = {"version": 1, "serializer": "msgpack"}


def _check_serialized_size(buffer: bytes) -> None:
    if len(buffer) > blosc2.MAX_BUFFERSIZE:
        raise ValueError(f"Serialized objects cannot be larger than {blosc2.MAX_BUFFERSIZE} bytes")


class ObjectArray:
    """A variable-length array backed by an :class:`blosc2.SChunk`.

    Entries are serialized with msgpack before compression. Standard Python
    objects are supported, and Blosc2 containers such as
    :class:`blosc2.NDArray`, :class:`blosc2.SChunk`, :class:`blosc2.ObjectArray`,
    :class:`blosc2.BatchArray`, and :class:`blosc2.EmbedStore` are serialized
    transparently via :meth:`to_cframe` / :func:`blosc2.from_cframe`.

    Msgpack also supports structured Blosc2 reference objects. Currently this
    includes :class:`blosc2.C2Array`, :class:`blosc2.LazyExpr`, and
    :class:`blosc2.LazyUDF` backed by :func:`blosc2.dsl_kernel`. Lazy
    expressions and supported lazy UDFs are serialized as recipes plus durable
    operand references, so only persistent local operands,
    :class:`blosc2.C2Array` operands, and :class:`blosc2.DictStore` members are
    supported. Purely in-memory operands are intentionally rejected. Plain
    Python :class:`blosc2.LazyUDF` callables are not serialized by msgpack.
    """

    @staticmethod
    def _set_typesize_one(cparams: blosc2.CParams | dict | None) -> blosc2.CParams | dict:
        auto_use_dict = cparams is None
        if cparams is None:
            cparams = blosc2.CParams()
        elif isinstance(cparams, blosc2.CParams):
            cparams = copy.deepcopy(cparams)
        else:
            cparams = dict(cparams)
            auto_use_dict = "use_dict" not in cparams

        if isinstance(cparams, blosc2.CParams):
            cparams.typesize = 1
            if auto_use_dict and cparams.codec == blosc2.Codec.ZSTD and cparams.clevel > 0:
                # ObjectArray stores many small serialized payloads; Zstd dicts help materially.
                cparams.use_dict = True
        else:
            cparams["typesize"] = 1
            codec = cparams.get("codec", blosc2.Codec.ZSTD)
            clevel = cparams.get("clevel", 5)
            if auto_use_dict and codec == blosc2.Codec.ZSTD and clevel > 0:
                # ObjectArray stores many small serialized payloads; Zstd dicts help materially.
                cparams["use_dict"] = True
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
            raise ValueError("For ObjectArray containers, mmap_mode must be None or 'r'")
        if storage.mmap_mode == "r" and storage.mode != "r":
            raise ValueError("For ObjectArray containers, mmap_mode='r' requires mode='r'")

    def _attach_schunk(self, schunk: SChunk) -> None:
        self.schunk = schunk
        self.mode = schunk.mode
        self.mmap_mode = getattr(schunk, "mmap_mode", None)
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
        chunksize: int | None = None,
        _from_schunk: SChunk | None = None,
        **kwargs: Any,
    ) -> None:
        if _from_schunk is not None:
            if chunksize is not None:
                raise ValueError("Cannot pass `chunksize` together with `_from_schunk`")
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
            raise ValueError(f"Unsupported ObjectArray keyword argument(s): {unexpected}")

        self._validate_storage(storage)
        cparams = self._set_typesize_one(cparams)

        if dparams is None:
            dparams = blosc2.DParams()

        if self._maybe_open_existing(storage):
            return

        fixed_meta = dict(storage.meta or {})
        fixed_meta["vlarray"] = dict(_VLARRAY_META)
        storage.meta = fixed_meta
        if chunksize is None:
            chunksize = -1
        schunk = blosc2.SChunk(
            chunksize=chunksize, data=None, cparams=cparams, dparams=dparams, storage=storage
        )
        self._attach_schunk(schunk)

    def _validate_tag(self) -> None:
        if "vlarray" not in self.schunk.meta:
            raise ValueError("The supplied SChunk is not tagged as an ObjectArray")

    def _check_writable(self) -> None:
        if self.mode == "r":
            raise ValueError("Cannot modify an ObjectArray opened in read-only mode")

    def _normalize_index(self, index: int) -> int:
        if not isinstance(index, int):
            raise TypeError("ObjectArray indices must be integers")
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError("ObjectArray index out of range")
        return index

    def _normalize_insert_index(self, index: int) -> int:
        if not isinstance(index, int):
            raise TypeError("ObjectArray indices must be integers")
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

    def _item_size_stats(self) -> tuple[list[int], list[int]]:
        item_nbytes = []
        chunk_cbytes = []
        for i in range(len(self)):
            nbytes, cbytes, _ = blosc2.get_cbuffer_sizes(self.schunk.get_lazychunk(i))
            item_nbytes.append(nbytes)
            chunk_cbytes.append(cbytes)
        return item_nbytes, chunk_cbytes

    def _serialize(self, value: Any) -> bytes:
        payload = msgpack_packb(value)
        _check_serialized_size(payload)
        return payload

    def _compress(self, payload: bytes) -> bytes:
        return blosc2.compress2(payload, cparams=self.schunk.cparams)

    def append(self, value: Any) -> int:
        """Append one value and return the new number of entries."""
        self._check_writable()
        chunk = self._compress(self._serialize(value))
        return self.schunk.append_chunk(chunk)

    def insert(self, index: int, value: Any) -> int:
        """Insert one value at ``index`` and return the new number of entries."""
        self._check_writable()
        index = self._normalize_insert_index(index)
        chunk = self._compress(self._serialize(value))
        return self.schunk.insert_chunk(index, chunk)

    def delete(self, index: int) -> int:
        """Delete the value at ``index`` and return the new number of entries."""
        self._check_writable()
        if isinstance(index, slice):
            # Delete in descending order so earlier deletions don't shift
            # the indices of chunks yet to be deleted (negative-step slices
            # produce ascending indices when merely reversed).
            for idx in sorted(self._slice_indices(index), reverse=True):
                self.schunk.delete_chunk(idx)
            return len(self)
        index = self._normalize_index(index)
        return self.schunk.delete_chunk(index)

    def pop(self, index: int = -1) -> Any:
        """Remove and return the value at ``index``."""
        self._check_writable()
        if isinstance(index, slice):
            raise NotImplementedError("Slicing is not supported for ObjectArray.pop()")
        index = self._normalize_index(index)
        value = self[index]
        self.schunk.delete_chunk(index)
        return value

    def extend(self, values: object) -> None:
        """Append all values from an iterable."""
        self._check_writable()
        for value in values:
            chunk = self._compress(self._serialize(value))
            self.schunk.append_chunk(chunk)

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

    def __getitem__(self, index: int) -> Any:
        if isinstance(index, slice):
            return [self[i] for i in self._slice_indices(index)]
        index = self._normalize_index(index)
        payload = self.schunk.decompress_chunk(index)
        return msgpack_unpackb(payload)

    def __setitem__(self, index: int, value: Any) -> None:
        if isinstance(index, slice):
            self._check_writable()
            indices = self._slice_indices(index)
            values = list(value)
            step = 1 if index.step is None else index.step
            if step == 1:
                start = self._normalize_insert_index(0 if index.start is None else index.start)
                for idx in reversed(indices):
                    self.schunk.delete_chunk(idx)
                for offset, item in enumerate(values):
                    chunk = self._compress(self._serialize(item))
                    self.schunk.insert_chunk(start + offset, chunk)
                return
            if len(values) != len(indices):
                raise ValueError(
                    f"attempt to assign sequence of size {len(values)} to extended slice of size {len(indices)}"
                )
            for idx, item in zip(indices, values, strict=True):
                chunk = self._compress(self._serialize(item))
                self.schunk.update_chunk(idx, chunk)
            return
        self._check_writable()
        index = self._normalize_index(index)
        chunk = self._compress(self._serialize(value))
        self.schunk.update_chunk(index, chunk)

    def __delitem__(self, index: int) -> None:
        self.delete(index)

    def __len__(self) -> int:
        return self.schunk.nchunks

    def __iter__(self) -> Iterator[Any]:
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
    def chunksize(self) -> int:
        return self.schunk.chunksize

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
        """Print information about this ObjectArray."""
        return InfoReporter(self)

    @property
    def info_items(self) -> list:
        """A list of tuples with summary information about this ObjectArray."""
        item_nbytes, chunk_cbytes = self._item_size_stats()
        avg_item_nbytes = sum(item_nbytes) / len(item_nbytes) if item_nbytes else 0.0
        avg_chunk_cbytes = sum(chunk_cbytes) / len(chunk_cbytes) if chunk_cbytes else 0.0
        return [
            ("type", f"{self.__class__.__name__}"),
            ("entries", len(self)),
            ("item_nbytes_min", min(item_nbytes) if item_nbytes else 0),
            ("item_nbytes_max", max(item_nbytes) if item_nbytes else 0),
            ("item_nbytes_avg", f"{avg_item_nbytes:.2f}"),
            ("chunk_cbytes_min", min(chunk_cbytes) if chunk_cbytes else 0),
            ("chunk_cbytes_max", max(chunk_cbytes) if chunk_cbytes else 0),
            ("chunk_cbytes_avg", f"{avg_chunk_cbytes:.2f}"),
            ("nbytes", format_nbytes_info(self.nbytes)),
            ("cbytes", format_nbytes_info(self.cbytes)),
            ("cratio", f"{self.cratio:.2f}x"),
            ("cparams", self.cparams),
            ("dparams", self.dparams),
        ]

    def to_cframe(self) -> bytes:
        return self.schunk.to_cframe()

    def copy(self, **kwargs: Any) -> ObjectArray:
        """Create a copy of the container with optional constructor overrides."""
        if "meta" in kwargs:
            raise ValueError("meta should not be passed to copy")

        kwargs["cparams"] = kwargs.get("cparams", copy.deepcopy(self.cparams))
        kwargs["dparams"] = kwargs.get("dparams", copy.deepcopy(self.dparams))
        kwargs["chunksize"] = kwargs.get("chunksize", -1)

        if "storage" not in kwargs:
            kwargs["meta"] = self._copy_meta()
            kwargs["contiguous"] = kwargs.get("contiguous", self.schunk.contiguous)
            if "urlpath" in kwargs and "mode" not in kwargs:
                kwargs["mode"] = "w"

        out = ObjectArray(**kwargs)
        out.extend(self)
        return out

    def __enter__(self) -> ObjectArray:
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        return False

    def __repr__(self) -> str:
        return f"ObjectArray(len={len(self)}, urlpath={self.urlpath!r})"


def objectarray_from_cframe(cframe: bytes, copy: bool = True) -> ObjectArray:
    """Deserialize a CFrame buffer into an :class:`ObjectArray`."""

    schunk = blosc2.schunk_from_cframe(cframe, copy=copy)
    return ObjectArray(_from_schunk=schunk)
