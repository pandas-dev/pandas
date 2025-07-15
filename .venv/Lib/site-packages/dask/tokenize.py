from __future__ import annotations

import copyreg
import dataclasses
import datetime
import decimal
import hashlib
import inspect
import pathlib
import pickle
import threading
import types
import uuid
from collections import OrderedDict
from collections.abc import Iterable
from contextvars import ContextVar
from functools import partial

import cloudpickle
from tlz import curry, identity
from tlz.functoolz import Compose

from dask import config
from dask.core import literal
from dask.hashing import hash_buffer_hex
from dask.utils import Dispatch


class TokenizationError(RuntimeError):
    pass


def _tokenize(*args: object, **kwargs: object) -> str:
    token: object = _normalize_seq_func(args)
    if kwargs:
        token = token, _normalize_seq_func(sorted(kwargs.items()))

    # Pass `usedforsecurity=False` to support FIPS builds of Python
    return hashlib.md5(str(token).encode(), usedforsecurity=False).hexdigest()


tokenize_lock = threading.RLock()
_SEEN: dict[int, tuple[int, object]] = {}
_ENSURE_DETERMINISTIC: ContextVar[bool | None] = ContextVar("_ENSURE_DETERMINISTIC")


def tokenize(
    *args: object, ensure_deterministic: bool | None = None, **kwargs: object
) -> str:
    """Deterministic token

    >>> tokenize([1, 2, '3'])  # doctest: +SKIP
    '06961e8de572e73c2e74b51348177918'

    >>> tokenize('Hello') == tokenize('Hello')
    True

    Parameters
    ----------
    args, kwargs:
        objects to tokenize
    ensure_deterministic: bool, optional
        If True, raise TokenizationError if the objects cannot be deterministically
        tokenized, e.g. two identical objects will return different tokens.
        Defaults to the `tokenize.ensure-deterministic` configuration parameter.
    """
    global _SEEN
    with tokenize_lock:
        seen_before, _SEEN = _SEEN, {}
        token = None
        try:
            _ENSURE_DETERMINISTIC.get()
        except LookupError:
            token = _ENSURE_DETERMINISTIC.set(ensure_deterministic)
        try:
            return _tokenize(*args, **kwargs)
        finally:
            if token:
                _ENSURE_DETERMINISTIC.reset(token)
            _SEEN = seen_before


def _maybe_raise_nondeterministic(msg: str) -> None:
    try:
        val = _ENSURE_DETERMINISTIC.get()
    except LookupError:
        val = None
    if val or val is None and config.get("tokenize.ensure-deterministic"):
        raise TokenizationError(msg)


_IDENTITY_DISPATCH = (
    int,
    float,
    str,
    bytes,
    type(None),
    slice,
    complex,
    type(Ellipsis),
    decimal.Decimal,
    datetime.date,
    datetime.time,
    datetime.datetime,
    datetime.timedelta,
    pathlib.PurePath,
)
normalize_token = Dispatch()
normalize_token.register(
    _IDENTITY_DISPATCH,
    identity,
)


@normalize_token.register((types.MappingProxyType, dict))
def normalize_dict(d):
    with tokenize_lock:
        if id(d) in _SEEN:
            return "__seen", _SEEN[id(d)][0]
        _SEEN[id(d)] = len(_SEEN), d
        try:
            return "dict", _normalize_seq_func(
                sorted(d.items(), key=lambda kv: str(kv[0]))
            )
        finally:
            _SEEN.pop(id(d), None)


@normalize_token.register(OrderedDict)
def normalize_ordered_dict(d):
    return _normalize_seq_func((type(d), list(d.items())))


@normalize_token.register(set)
def normalize_set(s):
    # Note: in some Python version / OS combinations, set order changes every
    # time you recreate the set (even within the same interpreter).
    # In most other cases, set ordering is consistent within the same interpreter.
    return "set", _normalize_seq_func(sorted(s, key=str))


def _normalize_seq_func(seq: Iterable[object]) -> tuple[object, ...]:
    def _inner_normalize_token(item):
        # Don't go through Dispatch. That's slow
        if isinstance(item, _IDENTITY_DISPATCH):
            return item
        return normalize_token(item)

    with tokenize_lock:
        if id(seq) in _SEEN:
            return "__seen", _SEEN[id(seq)][0]
        _SEEN[id(seq)] = len(_SEEN), seq
        try:
            return tuple(map(_inner_normalize_token, seq))
        finally:
            del _SEEN[id(seq)]


@normalize_token.register((tuple, list))
def normalize_seq(seq):
    return type(seq).__name__, _normalize_seq_func(seq)


@normalize_token.register(literal)
def normalize_literal(lit):
    return "literal", normalize_token(lit())


@normalize_token.register(Compose)
def normalize_compose(func):
    return _normalize_seq_func((func.first,) + func.funcs)


@normalize_token.register((partial, curry))
def normalize_partial(func):
    return _normalize_seq_func((func.func, func.args, func.keywords))


@normalize_token.register((types.MethodType, types.MethodWrapperType))
def normalize_bound_method(meth):
    return normalize_token(meth.__self__), meth.__name__


@normalize_token.register(types.BuiltinFunctionType)
def normalize_builtin_function_or_method(func):
    # Note: BuiltinMethodType is BuiltinFunctionType
    self = getattr(func, "__self__", None)
    if self is not None and not inspect.ismodule(self):
        return normalize_bound_method(func)
    else:
        return normalize_object(func)


@normalize_token.register(object)
def normalize_object(o):
    method = getattr(o, "__dask_tokenize__", None)
    if method is not None and not isinstance(o, type):
        return method()

    if type(o) is object:
        return _normalize_pure_object(o)

    if isinstance(o, type):
        copyreg._slotnames(o)

    if dataclasses.is_dataclass(o) and not isinstance(o, type):
        return _normalize_dataclass(o)

    try:
        return _normalize_pickle(o)
    except Exception:
        _maybe_raise_nondeterministic(
            f"Object {o!r} cannot be deterministically hashed. This likely "
            "indicates that the object cannot be serialized deterministically."
        )
        return uuid.uuid4().hex


_seen_objects = set()


def _normalize_pure_object(o: object) -> tuple[str, int]:
    _maybe_raise_nondeterministic(
        "object() cannot be deterministically hashed. See "
        "https://docs.dask.org/en/latest/custom-collections.html#implementing-deterministic-hashing "
        "for more information."
    )
    # Idempotent, but not deterministic. Make sure that the id is not reused.
    _seen_objects.add(o)
    return "object", id(o)


def _normalize_pickle(o: object) -> tuple:
    buffers: list[pickle.PickleBuffer] = []
    pik: int | None = None
    pik2: int | None = None
    for _ in range(3):
        buffers.clear()
        try:
            out = pickle.dumps(o, protocol=5, buffer_callback=buffers.append)
            if b"__main__" in out:
                # Use `cloudpickle` for objects defined in `__main__`
                buffers.clear()
                out = cloudpickle.dumps(o, protocol=5, buffer_callback=buffers.append)
            pickle.loads(out, buffers=buffers)
            pik2 = hash_buffer_hex(out)
        except Exception:
            buffers.clear()
            try:
                out = cloudpickle.dumps(o, protocol=5, buffer_callback=buffers.append)
                pickle.loads(out, buffers=buffers)
                pik2 = hash_buffer_hex(out)
            except Exception:
                break
        if pik and pik2 and pik == pik2:
            break
        pik = pik2
    else:
        _maybe_raise_nondeterministic("Failed to tokenize deterministically")
    if pik is None:
        _maybe_raise_nondeterministic("Failed to tokenize deterministically")
        pik = int(uuid.uuid4())
    return pik, [hash_buffer_hex(buf) for buf in buffers]


def _normalize_dataclass(obj):
    fields = [
        (field.name, normalize_token(getattr(obj, field.name, None)))
        for field in dataclasses.fields(obj)
    ]
    params = obj.__dataclass_params__
    params = [(attr, getattr(params, attr)) for attr in params.__slots__]

    return normalize_object(type(obj)), params, fields


@normalize_token.register_lazy("pandas")
def register_pandas():
    import pandas as pd

    @normalize_token.register(pd.RangeIndex)
    def normalize_range_index(x):
        return type(x), x.start, x.stop, x.step, x.dtype, x.name

    @normalize_token.register(pd.Index)
    def normalize_index(ind):
        values = ind.array
        return type(ind), ind.name, normalize_token(values)

    @normalize_token.register(pd.MultiIndex)
    def normalize_index(ind):
        codes = ind.codes
        return (
            [ind.name]
            + [normalize_token(x) for x in ind.levels]
            + [normalize_token(x) for x in codes]
        )

    @normalize_token.register(pd.Categorical)
    def normalize_categorical(cat):
        return [normalize_token(cat.codes), normalize_token(cat.dtype)]

    @normalize_token.register(pd.arrays.PeriodArray)
    @normalize_token.register(pd.arrays.DatetimeArray)
    @normalize_token.register(pd.arrays.TimedeltaArray)
    def normalize_period_array(arr):
        return [normalize_token(arr.asi8), normalize_token(arr.dtype)]

    @normalize_token.register(pd.arrays.IntervalArray)
    def normalize_interval_array(arr):
        return [
            normalize_token(arr.left),
            normalize_token(arr.right),
            normalize_token(arr.closed),
        ]

    @normalize_token.register(pd.Series)
    def normalize_series(s):
        return [
            s.name,
            s.dtype,
            normalize_token(s._values),
            normalize_token(s.index),
        ]

    @normalize_token.register(pd.DataFrame)
    def normalize_dataframe(df):
        mgr = df._mgr
        data = list(mgr.arrays) + [df.columns, df.index]
        return list(map(normalize_token, data))

    @normalize_token.register(pd.arrays.ArrowExtensionArray)
    def normalize_extension_array(arr):
        try:
            return (type(arr), normalize_token(arr._pa_array))
        except AttributeError:
            return (type(arr), normalize_token(arr._data))

    @normalize_token.register(pd.api.extensions.ExtensionArray)
    def normalize_extension_array(arr):
        import numpy as np

        return normalize_token(np.asarray(arr))

    # Dtypes
    @normalize_token.register(pd.api.types.CategoricalDtype)
    def normalize_categorical_dtype(dtype):
        return [normalize_token(dtype.categories), normalize_token(dtype.ordered)]

    @normalize_token.register(pd.api.extensions.ExtensionDtype)
    def normalize_period_dtype(dtype):
        return normalize_token(dtype.name)

    @normalize_token.register(type(pd.NA))
    def normalize_na(na):
        return pd.NA

    @normalize_token.register(pd.offsets.BaseOffset)
    def normalize_offset(offset):
        return offset.freqstr


@normalize_token.register_lazy("numba")
def register_numba():
    import numba

    @normalize_token.register(numba.core.serialize.ReduceMixin)
    def normalize_numba_ufunc(obj):
        return normalize_token((obj._reduce_class(), obj._reduce_states()))


@normalize_token.register_lazy("pyarrow")
def register_pyarrow():
    import pyarrow as pa

    @normalize_token.register(pa.DataType)
    def normalize_datatype(dt):
        return pickle.dumps(dt, protocol=4)

    @normalize_token.register(pa.Table)
    def normalize_table(dt):
        return (
            "pa.Table",
            normalize_token(dt.schema),
            normalize_token(dt.columns),
        )

    @normalize_token.register(pa.ChunkedArray)
    def normalize_chunked_array(arr):
        return (
            "pa.ChunkedArray",
            normalize_token(arr.type),
            normalize_token(arr.chunks),
        )

    @normalize_token.register(pa.Array)
    def normalize_chunked_array(arr):
        return (
            "pa.Array",
            normalize_token(arr.type),
            normalize_token(arr.buffers()),
        )

    @normalize_token.register(pa.Buffer)
    def normalize_chunked_array(buf):
        return ("pa.Buffer", hash_buffer_hex(buf))


@normalize_token.register_lazy("numpy")
def register_numpy():
    import numpy as np

    @normalize_token.register(np.ndarray)
    def normalize_array(x):
        if not x.shape:
            return (x.item(), x.dtype)
        if x.dtype.hasobject:
            try:
                try:
                    # string fast-path
                    data = hash_buffer_hex(
                        "-".join(x.flat).encode(
                            encoding="utf-8", errors="surrogatepass"
                        )
                    )
                except UnicodeDecodeError:
                    # bytes fast-path
                    data = hash_buffer_hex(b"-".join(x.flat))
            except (TypeError, UnicodeDecodeError):
                return normalize_object(x)
        else:
            try:
                data = hash_buffer_hex(x.ravel(order="K").view("i1"))
            except (BufferError, AttributeError, ValueError):
                data = hash_buffer_hex(x.copy().ravel(order="K").view("i1"))
        return (data, x.dtype, x.shape)

    @normalize_token.register(np.memmap)
    def normalize_mmap(mm):
        return hash_buffer_hex(np.ascontiguousarray(mm))

    @normalize_token.register(np.ufunc)
    def normalize_ufunc(func):
        try:
            return _normalize_pickle(func)
        except Exception:
            _maybe_raise_nondeterministic(
                f"Cannot tokenize numpy ufunc {func!r}. Please use functions "
                "of the dask.array.ufunc module instead. See also "
                "https://docs.dask.org/en/latest/array-numpy-compatibility.html"
            )
            return uuid.uuid4().hex

    @normalize_token.register(np.dtype)
    def normalize_dtype(dtype):
        return dtype.str


def _tokenize_deterministic(*args, **kwargs) -> str:
    # Utility to be strict about deterministic tokens
    return tokenize(*args, ensure_deterministic=True, **kwargs)
