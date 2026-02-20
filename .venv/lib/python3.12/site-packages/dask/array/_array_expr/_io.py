from __future__ import annotations

import functools
from itertools import product

import numpy as np

from dask import istask
from dask.array._array_expr._expr import ArrayExpr
from dask.array.core import (
    getter,
    getter_nofancy,
    graph_from_arraylike,
    normalize_chunks,
    slices_from_chunks,
)
from dask.array.utils import meta_from_array
from dask.utils import SerializableLock


class IO(ArrayExpr):
    pass


class FromGraph(IO):
    _parameters = ["layer", "_meta", "chunks", "keys", "name_prefix"]

    @functools.cached_property
    def _meta(self):
        return self.operand("_meta")

    @functools.cached_property
    def chunks(self):
        return self.operand("chunks")

    @functools.cached_property
    def _name(self):
        return self.operand("name_prefix") + "-" + self.deterministic_token

    def _layer(self):
        dsk = dict(self.operand("layer"))
        # The name may not actually match the layers name therefore rewrite this
        # using an alias
        for k in self.operand("keys"):
            if not isinstance(k, tuple):
                raise TypeError(f"Expected tuple, got {type(k)}")
            orig = dsk[k]
            if not istask(orig):
                del dsk[k]
                dsk[(self._name, *k[1:])] = orig
            else:
                dsk[(self._name, *k[1:])] = k
        return dsk


class FromArray(IO):
    _parameters = [
        "array",
        "chunks",
        "lock",
        "getitem",
        "inline_array",
        "meta",
        "asarray",
        "fancy",
    ]
    _defaults = {
        "getitem": None,
        "inline_array": False,
        "meta": None,
        "asarray": None,
        "fancy": True,
        "lock": False,
    }

    @property
    def chunks(self):
        return normalize_chunks(
            self.operand("chunks"), self.array.shape, dtype=self.array.dtype
        )

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not None:
            return meta_from_array(self.operand("meta"), dtype=self.array.dtype)
        return meta_from_array(self.array, dtype=getattr(self.array, "dtype", None))

    @functools.cached_property
    def asarray_arg(self):
        if self.operand("asarray") is None:
            return not hasattr(self.array, "__array_function__")
        else:
            return self.operand("asarray")

    def _layer(self):
        lock = self.operand("lock")
        if lock is True:
            lock = SerializableLock()

        is_ndarray = type(self.array) in (np.ndarray, np.ma.core.MaskedArray)
        is_single_block = all(len(c) == 1 for c in self.chunks)
        # Always use the getter for h5py etc. Not using isinstance(x, np.ndarray)
        # because np.matrix is a subclass of np.ndarray.
        if is_ndarray and not is_single_block and not lock:
            # eagerly slice numpy arrays to prevent memory blowup
            # GH5367, GH5601
            slices = slices_from_chunks(self.chunks)
            keys = product([self._name], *(range(len(bds)) for bds in self.chunks))
            values = [self.array[slc] for slc in slices]
            dsk = dict(zip(keys, values))
        elif is_ndarray and is_single_block:
            # No slicing needed
            dsk = {(self._name,) + (0,) * self.array.ndim: self.array}
        else:
            getitem = self.operand("getitem")
            if getitem is None:
                if self.operand("fancy"):
                    getitem = getter
                else:
                    getitem = getter_nofancy

            dsk = graph_from_arraylike(
                self.array,
                chunks=self.chunks,
                shape=self.array.shape,
                name=self._name,
                lock=lock,
                getitem=getitem,
                asarray=self.asarray_arg,
                inline_array=self.inline_array,
                dtype=self.array.dtype,
            )
        return dict(dsk)  # this comes as a legacy HLG for now

    def __str__(self):
        return "FromArray(...)"
