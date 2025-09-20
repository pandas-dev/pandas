from __future__ import annotations

from importlib import import_module
from types import ModuleType
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from packaging.version import Version

from xarray.core.utils import is_scalar
from xarray.namedarray.utils import is_duck_array, is_duck_dask_array

integer_types = (int, np.integer)

if TYPE_CHECKING:
    ModType = Literal["dask", "pint", "cupy", "sparse", "cubed", "numbagg"]
    DuckArrayTypes = tuple[type[Any], ...]  # TODO: improve this? maybe Generic
    from xarray.namedarray._typing import _DType, _ShapeType, duckarray


class DuckArrayModule:
    """
    Solely for internal isinstance and version checks.

    Motivated by having to only import pint when required (as pint currently imports xarray)
    https://github.com/pydata/xarray/pull/5561#discussion_r664815718
    """

    module: ModuleType | None
    version: Version
    type: DuckArrayTypes
    available: bool

    def __init__(self, mod: ModType) -> None:
        duck_array_module: ModuleType | None
        duck_array_version: Version
        duck_array_type: DuckArrayTypes
        try:
            duck_array_module = import_module(mod)
            duck_array_version = Version(duck_array_module.__version__)

            if mod == "dask":
                duck_array_type = (import_module("dask.array").Array,)
            elif mod == "pint":
                duck_array_type = (duck_array_module.Quantity,)
            elif mod == "cupy":
                duck_array_type = (duck_array_module.ndarray,)
            elif mod == "sparse":
                duck_array_type = (duck_array_module.SparseArray,)
            elif mod == "cubed":
                duck_array_type = (duck_array_module.Array,)
            # Not a duck array module, but using this system regardless, to get lazy imports
            elif mod == "numbagg":
                duck_array_type = ()
            else:
                raise NotImplementedError

        except (ImportError, AttributeError):  # pragma: no cover
            duck_array_module = None
            duck_array_version = Version("0.0.0")
            duck_array_type = ()

        self.module = duck_array_module
        self.version = duck_array_version
        self.type = duck_array_type
        self.available = duck_array_module is not None


_cached_duck_array_modules: dict[ModType, DuckArrayModule] = {}


def _get_cached_duck_array_module(mod: ModType) -> DuckArrayModule:
    if mod not in _cached_duck_array_modules:
        duckmod = DuckArrayModule(mod)
        _cached_duck_array_modules[mod] = duckmod
        return duckmod
    else:
        return _cached_duck_array_modules[mod]


def array_type(mod: ModType) -> DuckArrayTypes:
    """Quick wrapper to get the array class of the module."""
    return _get_cached_duck_array_module(mod).type


def mod_version(mod: ModType) -> Version:
    """Quick wrapper to get the version of the module."""
    return _get_cached_duck_array_module(mod).version


def is_chunked_array(x: duckarray[Any, Any]) -> bool:
    return is_duck_dask_array(x) or (is_duck_array(x) and hasattr(x, "chunks"))


def is_0d_dask_array(x: duckarray[Any, Any]) -> bool:
    return is_duck_dask_array(x) and is_scalar(x)


def to_numpy(
    data: duckarray[Any, Any], **kwargs: dict[str, Any]
) -> np.ndarray[Any, np.dtype[Any]]:
    from xarray.core.indexing import ExplicitlyIndexed
    from xarray.namedarray.parallelcompat import get_chunked_array_type

    try:
        # for tests only at the moment
        return data.to_numpy()  # type: ignore[no-any-return,union-attr]
    except AttributeError:
        pass

    if isinstance(data, ExplicitlyIndexed):
        data = data.get_duck_array()  # type: ignore[no-untyped-call]

    # TODO first attempt to call .to_numpy() once some libraries implement it
    if is_chunked_array(data):
        chunkmanager = get_chunked_array_type(data)
        data, *_ = chunkmanager.compute(data, **kwargs)
    if isinstance(data, array_type("cupy")):
        data = data.get()
    # pint has to be imported dynamically as pint imports xarray
    if isinstance(data, array_type("pint")):
        data = data.magnitude
    if isinstance(data, array_type("sparse")):
        data = data.todense()
    data = np.asarray(data)

    return data


def to_duck_array(data: Any, **kwargs: dict[str, Any]) -> duckarray[_ShapeType, _DType]:
    from xarray.core.indexing import (
        ExplicitlyIndexed,
        ImplicitToExplicitIndexingAdapter,
    )
    from xarray.namedarray.parallelcompat import get_chunked_array_type

    if is_chunked_array(data):
        chunkmanager = get_chunked_array_type(data)
        loaded_data, *_ = chunkmanager.compute(data, **kwargs)  # type: ignore[var-annotated]
        return loaded_data

    if isinstance(data, ExplicitlyIndexed | ImplicitToExplicitIndexingAdapter):
        return data.get_duck_array()  # type: ignore[no-untyped-call, no-any-return]
    elif is_duck_array(data):
        return data
    else:
        return np.asarray(data)  # type: ignore[return-value]


async def async_to_duck_array(
    data: Any, **kwargs: dict[str, Any]
) -> duckarray[_ShapeType, _DType]:
    from xarray.core.indexing import (
        ExplicitlyIndexed,
        ImplicitToExplicitIndexingAdapter,
    )

    if isinstance(data, ExplicitlyIndexed | ImplicitToExplicitIndexingAdapter):
        return await data.async_get_duck_array()  # type: ignore[union-attr, no-any-return]
    else:
        return to_duck_array(data, **kwargs)
