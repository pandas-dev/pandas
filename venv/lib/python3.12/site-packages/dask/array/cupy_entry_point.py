from __future__ import annotations

import dask.array as da
from dask import config
from dask.array.backends import ArrayBackendEntrypoint, register_cupy
from dask.array.core import Array
from dask.array.dispatch import to_cupy_dispatch


def _cupy(strict=True):
    try:
        import cupy
    except ImportError:
        if strict:
            raise ImportError("Please install `cupy` to use `CupyBackendEntrypoint`")
        return None
    return cupy


def _da_with_cupy_meta(attr, *args, meta=None, **kwargs):
    # Call the dask.array api with cupy-based meta
    meta = _cupy().empty(()) if meta is None else meta
    with config.set({"array.backend": "numpy"}):
        return getattr(da, attr)(*args, meta=meta, **kwargs)


class CupyBackendEntrypoint(ArrayBackendEntrypoint):
    def __init__(self):
        """Register data-directed dispatch functions"""
        if _cupy(strict=False):
            register_cupy()

    @classmethod
    def to_backend_dispatch(cls):
        return to_cupy_dispatch

    @classmethod
    def to_backend(cls, data: Array, **kwargs):
        if isinstance(data._meta, _cupy().ndarray):
            # Already a cupy-backed collection
            return data
        return data.map_blocks(cls.to_backend_dispatch(), **kwargs)

    @property
    def RandomState(self):
        return _cupy().random.RandomState

    @property
    def default_bit_generator(self):
        return _cupy().random.XORWOW

    @staticmethod
    def ones(*args, **kwargs):
        return _da_with_cupy_meta("ones", *args, **kwargs)

    @staticmethod
    def zeros(*args, **kwargs):
        return _da_with_cupy_meta("zeros", *args, **kwargs)

    @staticmethod
    def empty(*args, **kwargs):
        return _da_with_cupy_meta("empty", *args, **kwargs)

    @staticmethod
    def full(*args, **kwargs):
        return _da_with_cupy_meta("full", *args, **kwargs)

    @staticmethod
    def arange(*args, like=None, **kwargs):
        like = _cupy().empty(()) if like is None else like
        with config.set({"array.backend": "numpy"}):
            return da.arange(*args, like=like, **kwargs)
