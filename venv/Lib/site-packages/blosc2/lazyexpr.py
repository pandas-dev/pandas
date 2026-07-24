#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

# Avoid checking the name of type annotations at run time
from __future__ import annotations

import ast
import builtins
import concurrent.futures
import contextlib
import copy
import enum
import inspect
import linecache
import math
import os
import pathlib
import re
import sys
import textwrap
import threading
from abc import ABC, abstractmethod, abstractproperty
from collections.abc import MutableMapping
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from queue import Empty, Full, Queue
from typing import TYPE_CHECKING, Any

from numpy.exceptions import ComplexWarning

from . import exceptions

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

import ndindex
import numpy as np

import blosc2
from blosc2 import compute_chunks_blocks
from blosc2.info import InfoReporter

from .b2objects import (
    encode_b2object_payload,
    make_b2object_carrier,
    read_b2object_user_vlmeta,
    write_b2object_payload,
    write_b2object_user_vlmeta,
)
from .dsl_kernel import DSLKernel, DSLSyntaxError, DSLValidator, specialize_miniexpr_inputs
from .proxy import convert_dtype
from .utils import (
    check_smaller_shape,
    compute_smaller_slice,
    constructors,
    elementwise_funcs,
    format_expr_scalar,
    get_chunk_operands,
    get_chunks_idx,
    get_intersecting_chunks,
    infer_shape,
    linalg_attrs,
    linalg_funcs,
    npcumprod,
    npcumsum,
    populate_safe_numpy_globals,
    process_key,
    reducers,
    safe_numpy_globals,
    sliced_chunk_iter,
    try_miniexpr,
)

if not blosc2.IS_WASM:
    import numexpr

global safe_blosc2_globals
safe_blosc2_globals = {}


def ne_evaluate(expression, local_dict=None, **kwargs):
    """Safely evaluate expressions using numexpr when possible, falling back to numpy."""
    if local_dict is None:
        local_dict = {}
    # Get local vars dict from the stack frame
    _frame_depth = kwargs.pop("_frame_depth", 1)
    local_dict |= {
        k: v
        for k, v in dict(sys._getframe(_frame_depth).f_locals).items()
        if (
            (hasattr(v, "shape") or np.isscalar(v))
            and
            # Do not overwrite the local_dict with the expression variables
            not (k in local_dict or k in ("_where_x", "_where_y"))
        )
    }
    if blosc2.IS_WASM:
        global safe_numpy_globals
        populate_safe_numpy_globals(expression)
        if "out" in kwargs:
            out = kwargs.pop("out")
            out[:] = eval(expression, safe_numpy_globals, local_dict)
            return out
        res = eval(expression, safe_numpy_globals, local_dict)
        return np.asarray(res) if not hasattr(res, "shape") else res
    try:
        return numexpr.evaluate(expression, local_dict=local_dict, **kwargs)
    except ValueError as e:
        if e.args and e.args[0] == "NumExpr 2 does not support Unicode as a dtype.":
            pass
        else:
            raise  # unsafe expression
    except Exception:
        pass
    # Try with blosc2 funcs as presence of non-numexpr funcs probably caused failure
    # ne_evaluate will need safe_blosc2_globals for some functions (e.g. clip, logaddexp,
    # startswith, matmul) that are implemented incompletely in numexpr/miniexpr or not implemented at all
    global safe_blosc2_globals
    if len(safe_blosc2_globals) == 0:
        # First eval call, fill blosc2_safe_globals
        safe_blosc2_globals = {"blosc2": blosc2}
        # Add all first-level blosc2 functions
        safe_blosc2_globals.update(
            {
                name: getattr(blosc2, name)
                for name in dir(blosc2)
                if callable(getattr(blosc2, name)) and not name.startswith("_")
            }
        )
        # Expression strings can carry non-finite float literals ("nan",
        # "inf" — the repr of such scalars, or typed by the user); numexpr
        # has no such constants, so this python-eval fallback must.
        safe_blosc2_globals.update({"nan": math.nan, "inf": math.inf})
    res = eval(expression, safe_blosc2_globals, local_dict)
    if "out" in kwargs:
        out = kwargs.pop("out")
        out[:] = res  # will handle calc/decomp if res is lazyarray
        return out
    return res[()] if isinstance(res, blosc2.Operand) else res


def _get_result(expression, chunk_operands, ne_args, where=None, indices=None, _order=None):
    chunk_indices = None

    # Apply the where condition (in result) — fusion path, evaluate before shortcut
    if where is not None and len(where) == 2:
        # x = chunk_operands["_where_x"]
        # y = chunk_operands["_where_y"]
        # numexpr is a bit faster than np.where, and we can fuse operations in this case
        new_expr = f"where({expression}, _where_x, _where_y)"
        return ne_evaluate(new_expr, chunk_operands, **ne_args), None

    # If the expression is a simple operand reference (e.g. "key", "o0"),
    # grab it directly from chunk_operands instead of calling ne_evaluate.
    # This avoids ~150 µs of numexpr parsing/setup overhead per chunk.
    _expr = expression.strip("()")
    if _expr in chunk_operands:
        result = chunk_operands[_expr]
    else:
        result = ne_evaluate(expression, chunk_operands, **ne_args)
    if where is None:
        return result, None
    elif len(where) == 1:
        x = chunk_operands["_where_x"]
        if (indices is not None) or (_order is not None):
            # Return indices only makes sense when the where condition is a tuple with one element
            # and result is a boolean array
            if len(x.shape) > 1:
                raise ValueError("argsort() and sort() only support 1D arrays")
            if result.dtype != np.bool_:
                raise ValueError("argsort() and sort() only support bool conditions")
            if _order:
                # We need to cumulate all the fields in _order, as well as indices
                chunk_indices = indices[result]
                result = x[_order][result]
            else:
                chunk_indices = None
                result = indices[result]
            return result, chunk_indices
        else:
            return x[result], None
    raise ValueError("The where condition must be a tuple with one or two elements")


# Define empty ndindex tuple for function defaults
NDINDEX_EMPTY_TUPLE = ndindex.Tuple()

# All the dtypes that are supported by the expression evaluator
dtype_symbols = {
    "int8": np.int8,
    "int16": np.int16,
    "int32": np.int32,
    "int64": np.int64,
    "uint8": np.uint8,
    "uint16": np.uint16,
    "uint32": np.uint32,
    "uint64": np.uint64,
    "float32": np.float32,
    "float64": np.float64,
    "complex64": np.complex64,
    "complex128": np.complex128,
    "bool": np.bool_,
    "str": np.str_,
    "bytes": np.bytes_,
    "i1": np.int8,
    "i2": np.int16,
    "i4": np.int32,
    "i8": np.int64,
    "u1": np.uint8,
    "u2": np.uint16,
    "u4": np.uint32,
    "u8": np.uint64,
    "f4": np.float32,
    "f8": np.float64,
    "c8": np.complex64,
    "c16": np.complex128,
    "b1": np.bool_,
    "S": np.str_,
    "V": np.bytes_,
}
blosc2_funcs = constructors + linalg_funcs + elementwise_funcs + reducers
# functions that have to be evaluated before chunkwise lazyexpr machinery
eager_funcs = linalg_funcs + reducers + ["slice"] + ["." + attr for attr in linalg_attrs]
functions = blosc2_funcs
_TRANSIENT_MASK_CPARAMS = blosc2.CParams(codec=blosc2.Codec.LZ4, clevel=5, filters=[blosc2.Filter.SHUFFLE])
_constructor_call_patterns = {name: re.compile(rf"\b{re.escape(name)}\s*\(") for name in constructors}


def _has_constructor_call(expression: str, constructor: str) -> bool:
    return _constructor_call_patterns[constructor].search(expression) is not None


def _find_constructor_call(expression: str, constructor: str) -> re.Match | None:
    return _constructor_call_patterns[constructor].search(expression)


relational_ops = ["==", "!=", "<", "<=", ">", ">="]
logical_ops = ["&", "|", "^", "~"]
not_complex_ops = ["maximum", "minimum", "<", "<=", ">", ">="]
funcs_2args = (
    "arctan2",
    "contains",
    "pow",
    "power",
    "nextafter",
    "copysign",
    "hypot",
    "maximum",
    "minimum",
    "startswith",
    "endswith",
)


def get_expr_globals(expression):
    """Build a dictionary of functions needed for evaluating the expression."""
    _globals = {"np": np, "blosc2": blosc2, "nan": math.nan, "inf": math.inf}
    # Only check for functions that actually appear in the expression.
    for func in functions:
        if func in expression:
            if hasattr(blosc2, func):
                _globals[func] = getattr(blosc2, func)
            else:
                try:
                    _globals[func] = safe_numpy_globals[func]
                except KeyError as e:
                    raise AttributeError(f"Function {func} not found in blosc2 or numpy") from e

    # Lazily support bare numpy calls not covered by the Blosc2 function list.
    populate_safe_numpy_globals(expression)
    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError:
        return _globals
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        func = node.func.id
        if func in _globals:
            continue
        if hasattr(blosc2, func):
            _globals[func] = getattr(blosc2, func)
        elif func in safe_numpy_globals:
            _globals[func] = safe_numpy_globals[func]

    return _globals


if not hasattr(enum, "member"):
    # copy-pasted from Lib/enum.py
    class MyMember:
        """
        Forces item to become an Enum member during class creation.
        """

        def __init__(self, value):
            self.value = value
else:
    MyMember = enum.member  # only available after python 3.11


class ReduceOp(Enum):
    """
    Available reduce operations.
    """

    # wrap as enum.member so that Python doesn't treat some funcs
    # as class methods (rather than Enum members)
    SUM = MyMember(np.add)
    PROD = MyMember(np.multiply)
    MEAN = MyMember(np.mean)
    STD = MyMember(np.std)
    VAR = MyMember(np.var)
    # Computing a median from partial results is not straightforward because the median
    # is a positional statistic, which means it depends on the relative ordering of all
    # the data points. Unlike statistics such as the sum or mean, you can't compute a median
    # from partial results without knowing the entire dataset, and this is way too expensive
    # for arrays that cannot typically fit in-memory (e.g. disk-based NDArray).
    # MEDIAN = np.median
    MAX = MyMember(np.maximum)
    MIN = MyMember(np.minimum)
    ANY = MyMember(np.any)
    ALL = MyMember(np.all)
    ARGMAX = MyMember(np.argmax)
    ARGMIN = MyMember(np.argmin)
    CUMULATIVE_SUM = MyMember(npcumsum)
    CUMULATIVE_PROD = MyMember(npcumprod)


class LazyArrayEnum(Enum):
    """
    Available LazyArrays.
    """

    Expr = 0
    UDF = 1


class LazyArrayVLMeta(MutableMapping):
    """User metadata attached to a LazyArray."""

    def __init__(self, lazyarr: LazyArray):
        self.lazyarr = lazyarr

    def __getitem__(self, key):
        return self.lazyarr._get_user_vlmeta()[key]

    def __setitem__(self, key, value):
        data = self.lazyarr._get_user_vlmeta()
        data[key] = value
        self.lazyarr._sync_user_vlmeta()

    def __delitem__(self, key):
        data = self.lazyarr._get_user_vlmeta()
        del data[key]
        self.lazyarr._sync_user_vlmeta()

    def __iter__(self):
        return iter(self.lazyarr._get_user_vlmeta())

    def __len__(self):
        return len(self.lazyarr._get_user_vlmeta())

    def getall(self):
        return self.lazyarr._get_user_vlmeta().copy()

    def __repr__(self):
        return repr(self.getall())

    def __str__(self):
        return str(self.getall())


class LazyArray(ABC, blosc2.Operand):
    """Base class for lazy array expressions that compute data on demand."""

    def _get_user_vlmeta(self) -> dict[str, Any]:
        if not hasattr(self, "_vlmeta_user"):
            self._vlmeta_user = {}
        return self._vlmeta_user

    def _set_user_vlmeta(self, metadata: dict[str, Any], *, sync: bool = True) -> None:
        self._vlmeta_user = dict(metadata)
        if sync:
            self._sync_user_vlmeta()

    def _sync_user_vlmeta(self) -> None:
        array = getattr(self, "array", None)
        if array is not None:
            write_b2object_user_vlmeta(array, self._get_user_vlmeta())

    @property
    def vlmeta(self) -> LazyArrayVLMeta:
        """User variable-length metadata for this LazyArray."""
        if not hasattr(self, "_vlmeta_proxy"):
            self._vlmeta_proxy = LazyArrayVLMeta(self)
        return self._vlmeta_proxy

    def __enter__(self) -> LazyArray:
        """Enter a context manager and return this lazy array."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        """Exit a context manager.

        Lazy arrays do not currently keep explicit closeable resources, so this
        is a logical no-op kept for API consistency with :func:`blosc2.open`.
        """
        return False

    @abstractmethod
    def argsort(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        """
        Return an :ref:`LazyArray` containing the positions selected by self.

        The LazyArray must be of bool dtype (e.g. a condition).

        Parameters
        ----------
        order: str, list of str, optional
            Specifies which fields to compare first, second, etc. A single
            field can be specified as a string. Not all fields need to be
            specified, only the ones by which the array is to be sorted.

        Returns
        -------
        out: :ref:`LazyArray`
            The positions of the :ref:`LazyArray` self that are True.
        """
        pass

    @abstractmethod
    def sort(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        """
        Return a sorted :ref:`LazyArray`.

        This is only valid for LazyArrays with structured dtypes.

        Parameters
        ----------
        order: str, list of str, optional
            Specifies which fields to compare first, second, etc. A single
            field can be specified as a string. Not all fields need to be
            specified, only the ones by which the array is to be sorted.

        Returns
        -------
        out: :ref:`LazyArray`
            A sorted :ref:`LazyArray`.
        """
        pass

    @abstractmethod
    def compute(
        self,
        item: slice | list[slice] | None = None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs: Any,
    ) -> blosc2.NDArray:
        """
        Return a :ref:`NDArray` containing the evaluation of the :ref:`LazyArray`.

        Parameters
        ----------
        item: slice, list of slices, optional
            If provided, item is used to slice the operands *prior* to computation; not to retrieve specified slices of
            the evaluated result. This difference between slicing operands and slicing the final expression
            is important when reductions or a where clause are used in the expression.

        fp_accuracy: :class:`blosc2.FPAccuracy`, optional
            Specifies the floating-point accuracy to be used during computation.
            By default, :attr:`blosc2.FPAccuracy.DEFAULT` is used.

        kwargs: Any, optional
            Keyword arguments that are supported by the :func:`empty` constructor.
            These arguments will be set in the resulting :ref:`NDArray`.
            Additionally, the following special kwargs are supported:

            - ``strict_miniexpr`` (bool): controls whether miniexpr compilation/execution
              failures are raised instead of silently falling back to regular chunked eval
              for non-DSL expressions.  Setting it ``True`` also opts a DSL kernel out of the
              WebAssembly prefer-js default, keeping it on miniexpr.

            - ``jit`` (bool | None): enable (``True``) or disable (``False``) JIT compilation
              of the expression via miniexpr.  When ``None`` (default), JIT is only used
              for DSL kernels; plain expressions are evaluated by the bytecode interpreter.
              Setting ``jit=True`` forces auto-lift of plain expressions into JIT-compiled
              kernels.

            - ``jit_backend`` (str | None): select the JIT compiler backend.  Valid
              values are ``"tcc"`` (bundled Tiny C Compiler), ``"cc"`` (system C
              compiler, e.g. gcc or clang), and ``"js"`` (transpile the DSL kernel to
              JavaScript; browser/Pyodide only — see below).  ``None`` (default) defers
              to the miniexpr default (``"tcc"``), except under WebAssembly where — unless
              ``jit=False`` — it *prefers* ``"js"`` for transpilable float DSL kernels and
              falls back to miniexpr otherwise.  Since ``"js"`` is itself JIT-compiled by
              the JS engine, ``jit=True`` prefers it too; force miniexpr with
              ``jit_backend="tcc"``/``"cc"``.

            - ``"js"`` backend (WebAssembly/Pyodide only): transpiles a
              :func:`blosc2.dsl_kernel` to JavaScript so it runs at the browser engine's
              optimized native speed.  It tends to beat the WASM miniexpr JIT (~2x) for
              float kernels dominated by arithmetic and control flow, and is roughly a
              wash for transcendental-heavy or trivial kernels.  Outside WebAssembly,
              ``jit_backend="js"`` raises.  Forcing ``"tcc"``/``"cc"`` always uses miniexpr.

            - ``BLOSC_ME_JIT`` environment variable: when set to ``"1"``, ``"true"``,
              ``"on"``, ``"tcc"``, or ``"cc"``, it forces ``jit=True`` and overrides
              both the ``jit`` and ``jit_backend`` arguments — this lets you switch
              JIT on or change backends from the command line without touching code.
              Setting it to ``"tcc"`` or ``"cc"`` also selects that backend.

            - ``BLOSC_ME_JIT_TRACE`` environment variable: when set to ``"1"``,
              ``"true"``, or ``"on"``, prints a one-line diagnostic to stdout
              showing which compute engine was selected (``miniexpr`` or
              ``ne_evaluate``), the JIT mode and backend if applicable, and the
              expression being evaluated.

        Returns
        -------
        out: :ref:`NDArray`
            A :ref:`NDArray` containing the result of evaluating the
            :ref:`LazyUDF` or :ref:`LazyExpr`.

        Notes
        -----
        * If self is a LazyArray from an udf, the kwargs used to store the resulting
          array will be the ones passed to the constructor in :func:`lazyudf` (except the
          `urlpath`) updated with the kwargs passed when calling this method.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> dtype = np.float64
        >>> shape = [3, 3]
        >>> size = shape[0] * shape[1]
        >>> a = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
        >>> b = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
        >>> #  Convert numpy arrays to Blosc2 arrays
        >>> a1 = blosc2.asarray(a)
        >>> b1 = blosc2.asarray(b)
        >>> # Perform the mathematical operation
        >>> expr = a1 + b1
        >>> output = expr.compute()
        >>> f"Result of a + b (lazy evaluation): {output[:]}"
        Result of a + b (lazy evaluation):
                    [[ 0.    1.25  2.5 ]
                    [ 3.75  5.    6.25]
                    [ 7.5   8.75 10.  ]]
        """
        pass

    @abstractmethod
    def __getitem__(self, item: int | slice | Sequence[slice]) -> np.ndarray:
        """
        Return a numpy.ndarray containing the evaluation of the :ref:`LazyArray`.

        Parameters
        ----------
        item: int, slice or sequence of slices
            If provided, item is used to slice the operands *prior* to computation; not to retrieve specified slices of
            the evaluated result. This difference between slicing operands and slicing the final expression
            is important when reductions or a where clause are used in the expression.

        Returns
        -------
        out: np.ndarray
            An array with the data containing the evaluated slice.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> dtype = np.float64
        >>> shape = [30, 4]
        >>> size = shape[0] * shape[1]
        >>> a = np.linspace(0, 10, num=size, dtype=dtype).reshape(shape)
        >>> b = np.linspace(0, 10, num=size, dtype=dtype).reshape(shape)
        >>> #  Convert numpy arrays to Blosc2 arrays
        >>> a1 = blosc2.asarray(a)
        >>> b1 = blosc2.asarray(b)
        >>> # Perform the mathematical operation
        >>> expr = a1 + b1  # LazyExpr expression
        >>> expr[3]
        [2.01680672 2.18487395 2.35294118 2.5210084 ]
        >>> expr[2:4]
        [[1.34453782 1.51260504 1.68067227 1.8487395 ]
        [2.01680672 2.18487395 2.35294118 2.5210084 ]]
        """
        pass

    @abstractmethod
    def save(self, **kwargs: Any) -> None:
        """
        Save the :ref:`LazyArray` on disk.

        Parameters
        ----------
        kwargs: Any, optional
            Keyword arguments that are supported by the :func:`empty` constructor.
            The `urlpath` must always be provided.

        Returns
        -------
        out: None

        Notes
        -----
        * All the operands of the LazyArray must be Python scalars, or :class:`blosc2.Array` objects.
        * If an operand is a :ref:`Proxy`, keep in mind that Python-Blosc2 will only be able to reopen it as such
          if its source is a :ref:`SChunk`, :ref:`NDArray` or a :ref:`C2Array` (see :func:`blosc2.open` notes
          section for more info).
        * This is currently only supported for :ref:`LazyExpr` and :ref:`LazyUDF`
          (including kernels decorated with :func:`blosc2.dsl_kernel`).
        * User metadata can be attached via :attr:`vlmeta`. For in-memory LazyArrays
          this stays in memory; for persisted LazyArrays it is serialized and restored
          on reopen.

        Examples
        --------
        >>> import blosc2
        >>> import numpy as np
        >>> dtype = np.float64
        >>> shape = [3, 3]
        >>> size = shape[0] * shape[1]
        >>> a = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
        >>> b = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
        >>> # Define file paths for storing the arrays
        >>> a1 = blosc2.asarray(a, urlpath='a_array.b2nd', mode='w')
        >>> b1 = blosc2.asarray(b, urlpath='b_array.b2nd', mode='w')
        >>> # Perform the mathematical operation to create a LazyExpr expression
        >>> expr = a1 + b1
        >>> # Save the LazyExpr to disk
        >>> expr.save(urlpath='lazy_array.b2nd', mode='w')
        >>> # Open and load the LazyExpr from disk
        >>> disk_expr = blosc2.open('lazy_array.b2nd', mode='r')
        >>> disk_expr[:2]
        [[0.   1.25 2.5 ]
        [3.75 5.   6.25]]
        """
        pass

    # Provide a way to serialize the LazyArray
    def to_cframe(self) -> bytes:
        """
        Compute LazyArray and convert to cframe.

        Returns
        -------
        out: bytes
            The buffer containing the serialized :ref:`NDArray` instance.
        """
        return self.compute().to_cframe()

    @abstractproperty
    def chunks(self) -> tuple[int]:
        """
        Return :ref:`LazyArray` chunks.
        """
        pass

    @abstractproperty
    def blocks(self) -> tuple[int]:
        """
        Return :ref:`LazyArray` blocks.
        """
        pass

    def get_chunk(self, nchunk):
        """Get the `nchunk` of the expression, evaluating only that one."""
        # Create an empty array with the chunkshape and dtype; this is fast
        shape = self.shape
        chunks = self.chunks
        # Calculate the shape of the (chunk) slice_ (especially at the end of the array)
        chunks_idx, _ = get_chunks_idx(shape, chunks)
        coords = tuple(np.unravel_index(nchunk, chunks_idx))
        slice_ = tuple(
            slice(c * s, min((c + 1) * s, shape[i]))
            for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
        )
        loc_chunks = tuple(s.stop - s.start for s in slice_)
        out = blosc2.empty(shape=self.chunks, dtype=self.dtype, chunks=self.chunks, blocks=self.blocks)
        if loc_chunks == self.chunks:
            self.compute(item=slice_, out=out)
        else:
            _slice_ = tuple(slice(0, s) for s in loc_chunks)
            out[_slice_] = self.compute(item=slice_)
        return out.schunk.get_chunk(0)


def convert_inputs(inputs):
    if not inputs or len(inputs) == 0:
        return []
    inputs_ = []
    for obj in inputs:
        # CTable Column — unwrap to the backing NDArray so shape and identity match.
        with contextlib.suppress(AttributeError):
            obj = obj._raw_col
        if not isinstance(obj, np.ndarray | blosc2.Operand) and not np.isscalar(obj):
            try:
                obj = blosc2.SimpleProxy(obj)
            except Exception:
                print(
                    "Inputs not being np.ndarray, Array or Python scalar objects"
                    " should be convertible to SimpleProxy."
                )
                raise
        inputs_.append(obj)
    return inputs_


def compute_broadcast_shape(arrays):
    """
    Returns the shape of the outcome of an operation with the input arrays.
    """
    # When dealing with UDFs, one can arrive params that are not arrays
    shapes = [arr.shape for arr in arrays if hasattr(arr, "shape") and arr is not np]
    return np.broadcast_shapes(*shapes) if shapes else None


def _jit_from_env(jit, jit_backend):
    """Apply BLOSC_ME_JIT environment variable to jit/jit_backend defaults."""
    env_jit = os.environ.get("BLOSC_ME_JIT", "")
    if not env_jit:
        return jit, jit_backend
    env_jit_lower = env_jit.lower()
    # Env var always wins over both jit= and jit_backend= for easy CLI experimentation.
    if env_jit_lower in ("1", "true", "on", "tcc", "cc"):
        jit = True
    if env_jit_lower in ("tcc", "cc"):
        jit_backend = env_jit_lower
    return jit, jit_backend


# Define the patterns for validation
validation_patterns = [
    r"[\;]",  # Flow control characters
    r"(^|[^\w])__[\w]+__($|[^\w])",  # Dunder methods
    r"\.\b(?!real|imag|T|mT|(\d*[eE]?[+-]?\d+)|(\d*[eE]?[+-]?\d+j)|\d*j\b|(sum|prod|min|max|std|mean|var|any|all|where)"
    r"\s*\([^)]*\)|[a-zA-Z_]\w*\s*\([^)]*\))",  # Attribute patterns
]

# Compile the blacklist regex
_blacklist_re = re.compile("|".join(validation_patterns))

# Define valid method names
valid_methods = {
    "sum",
    "prod",
    "min",
    "max",
    "std",
    "mean",
    "var",
    "any",
    "all",
    "where",
    "reshape",
    "slice",
}
valid_methods |= {"int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"}
valid_methods |= {"float32", "float64", "complex64", "complex128"}
valid_methods |= {"bool", "str", "bytes"}
valid_methods |= {
    name for name in dir(blosc2.NDArray) if not name.startswith("_")
}  # allow attributes and methods


def validate_expr(expr: str) -> None:
    """
    Validate expression for forbidden syntax and valid method names.

    Parameters
    ----------
    expr : str
        The expression to validate.

    Returns
    -------
    None
    """
    # Remove whitespace and skip quoted strings
    no_whitespace = re.sub(r"\s+", "", expr)
    skip_quotes = re.sub(r"(\'[^\']*\')", "", no_whitespace)

    # Check for forbidden patterns
    forbiddens = _blacklist_re.search(skip_quotes)
    if forbiddens is not None:
        raise ValueError(f"'{expr}' is not a valid expression.")

    # Check for invalid characters not covered by the tokenizer
    invalid_chars = re.compile(r"[^\w\s+\-*/%()[].,=<>!&|~^]")
    if invalid_chars.search(skip_quotes) is not None:
        invalid_chars = invalid_chars.findall(skip_quotes)
        raise ValueError(f"Expression {expr} contains invalid characters: {invalid_chars}")

    # Check for invalid method names
    method_calls = re.findall(r"\.\b(\w+)\s*\(", skip_quotes)
    for method in method_calls:
        if method not in valid_methods:
            raise ValueError(f"Invalid method name: {method}")


def extract_and_replace_slices(expr, operands):
    """
    Return new expression and operands with op.slice(...) replaced by temporary operands.
    """
    # Copy shapes and operands
    shapes = {k: () if not hasattr(v, "shape") else v.shape for k, v in operands.items()}
    new_ops = operands.copy()  # copy dictionary

    # Parse the expression
    tree = ast.parse(expr, mode="eval")

    # Mapping of AST nodes to new variable names
    replacements = {}

    class SliceCollector(ast.NodeTransformer):
        def visit_Call(self, node):
            # Recursively visit children first
            self.generic_visit(node)

            # Detect method calls: obj.slice(...)
            if isinstance(node.func, ast.Attribute) and node.func.attr == "slice":
                obj = node.func.value

                # If the object is already replaced, keep the replacement
                base_name = None
                if isinstance(obj, ast.Name):
                    base_name = obj.id
                elif isinstance(obj, ast.Call) and obj in replacements:
                    base_name = replacements[obj]["base_var"]

                # Build the full slice chain expression as a string
                full_expr = ast.unparse(node)

                # Create a new temporary variable
                new_var = f"o{len(new_ops)}"

                # Infer shape
                try:
                    shape = infer_shape(full_expr, shapes)
                except Exception as e:
                    print(f"Shape inference failed for {full_expr}: {e}")
                    shape = ()

                # Determine dtype
                dtype = new_ops[base_name].dtype if base_name else None

                # Create placeholder array
                if isinstance(new_ops[base_name], blosc2.NDArray):
                    new_op = blosc2.ones((1,) * len(shape), dtype=dtype)
                else:
                    new_op = np.ones((1,) * len(shape), dtype=dtype)

                new_ops[new_var] = new_op
                shapes[new_var] = shape

                # Record replacement
                replacements[node] = {"new_var": new_var, "base_var": base_name}

                # Replace the AST node with the new variable
                return ast.Name(id=new_var, ctx=ast.Load())

            return node

    # Transform the AST
    transformer = SliceCollector()
    new_tree = transformer.visit(tree)
    ast.fix_missing_locations(new_tree)

    # Convert back to expression string
    new_expr = ast.unparse(new_tree)

    return new_expr, new_ops


def get_expr_operands(expression: str) -> set:
    """
    Given an expression in string form, return its operands.

    Parameters
    ----------
    expression : str
        The expression in string form.

    Returns
    -------
    set
        A set of operands found in the expression.
    """

    class OperandVisitor(ast.NodeVisitor):
        def __init__(self):
            self.operands = set()
            self.function_names = set()

        def visit_Name(self, node):
            if node.id == "np":
                # Skip NumPy namespace (e.g. np.int8, which will be treated separately)
                return
            if node.id not in self.function_names and node.id not in dtype_symbols:
                self.operands.add(node.id)
            self.generic_visit(node)

        def visit_Call(self, node):
            if isinstance(node.func, ast.Name):
                self.function_names.add(node.func.id)
            self.generic_visit(node)

    tree = ast.parse(expression)
    visitor = OperandVisitor()
    visitor.visit(tree)
    return set(visitor.operands)


def conserve_functions(  # noqa: C901
    expression: str,
    operands_old: dict[str, blosc2.Array],
    operands_new: dict[str, blosc2.Array],
) -> tuple[str, dict[str, blosc2.Array]]:
    """
    Given an expression in string form, return its operands.

    Parameters
    ----------
    expression : str
        The expression in string form.

    operands_old: dict[str : blosc2.ndarray | blosc2.LazyExpr]
        Dict of operands from expression prior to eval.

    operands_new: dict[str : blosc2.ndarray | blosc2.LazyExpr]
        Dict of operands from expression after eval.
    Returns
    -------
    newexpression
        A modified string expression with the functions/constructors conserved and
        true operands rebased and written in o- notation.
    newoperands
        Dict of the set of rebased operands.
    """

    operand_to_key = {id(v): k for k, v in operands_new.items()}
    for k, v in operands_old.items():  # extend operands_to_key with old operands
        if isinstance(
            v, blosc2.LazyExpr
        ):  # unroll operands in LazyExpr (only necessary when have reduced a lazyexpr)
            d = v.operands
        else:
            d = {k: v}
        for newk, newv in d.items():
            try:
                operand_to_key[id(newv)]
            except KeyError:
                newk = (
                    f"o{len(operands_new)}" if newk in operands_new else newk
                )  # possible that names coincide
                operand_to_key[id(newv)] = newk
                operands_new[newk] = newv

    class OperandVisitor(ast.NodeVisitor):
        def __init__(self):
            self.operandmap = {}
            self.operands = {}
            self.opcounter = 0
            self.function_names = set()

        def update_func(self, localop):
            k = operand_to_key[id(localop)]
            if k not in self.operandmap:
                newkey = f"o{self.opcounter}"
                self.operands[newkey] = operands_new[k]
                self.operandmap[k] = newkey
                self.opcounter += 1
                return newkey
            else:
                return self.operandmap[k]

        def visit_Name(self, node):
            if node.id == "np":  # Skip NumPy namespace (e.g. np.int8, which will be treated separately)
                return
            if node.id in self.function_names:  # Skip function names
                return
            elif node.id not in dtype_symbols:
                if node.id in ("nan", "inf") and node.id not in operands_old:
                    return  # non-finite float literal, not an operand
                localop = operands_old[node.id]
                if isinstance(localop, blosc2.LazyExpr):
                    newexpr = localop.expression
                    for (
                        opname,
                        v,
                    ) in localop.operands.items():  # expression operands already in terms of basic operands
                        # add illegal character ; to track changed operands and not overwrite later
                        newopname = ";" + self.update_func(v)
                        newexpr = re.sub(
                            rf"(?<=\s){opname}|(?<=\(){opname}", newopname, newexpr
                        )  # replace with newopname
                    # remove all instances of ; as all changes completed
                    node.id = newexpr.replace(";", "")
                else:
                    node.id = self.update_func(localop)
            self.generic_visit(node)

        def visit_Call(self, node):
            if isinstance(
                node.func, ast.Name
            ):  # visits Call first, then Name, so don't increment operandcounter yet
                self.function_names.add(node.func.id)
            self.generic_visit(node)

    tree = ast.parse(expression)
    visitor = OperandVisitor()
    visitor.visit(tree)
    newexpression, newoperands = ast.unparse(tree), visitor.operands
    return newexpression, newoperands


def convert_to_slice(expression):
    """
    Takes expression and converts all instances of [] to .slice(....)

    Parameters
    ----------
    expression: str

    Returns
    -------
    new_expr : str
    """

    new_expr = ""
    skip_to_char = 0
    for i, expr_i in enumerate(expression):
        if i < skip_to_char:
            continue
        if expr_i == "[":
            k = expression[i:].find("]")  # start checking from after [
            slice_convert = expression[i : i + k + 1]  # include [ and ]
            try:
                slicer = eval(f"np.s_{slice_convert}")
                slicer = (slicer,) if not isinstance(slicer, tuple) else slicer  # standardise to tuple
                if any(isinstance(el, str) for el in slicer):  # handle fields
                    raise ValueError("Cannot handle fields for slicing lazy expressions.")
                slicer = str(slicer)
                # use slice so that lazyexpr uses blosc arrays internally
                # (and doesn't decompress according to getitem syntax)
                new_expr += f".slice({slicer})"
                skip_to_char = i + k + 1
                continue
            except Exception:
                pass
        new_expr += expr_i  # if slice_convert is e.g. a list, not a slice, do nothing
    return new_expr


class TransformNumpyCalls(ast.NodeTransformer):
    def __init__(self):
        self.replacements = {}
        self.tmp_counter = 0

    def visit_Call(self, node):
        # Check if the call is a numpy type-casting call
        if (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id in ["np", "numpy"]
            and isinstance(node.args[0], ast.Constant)
        ):
            # Create a new temporary variable name
            tmp_var = f"tmp{self.tmp_counter}"
            self.tmp_counter += 1

            # Evaluate the type-casting call to create the new variable's value
            numpy_type = getattr(np, node.func.attr)
            self.replacements[tmp_var] = numpy_type(node.args[0].value)

            # Replace the call node with a variable node
            return ast.copy_location(ast.Name(id=tmp_var, ctx=ast.Load()), node)
        return self.generic_visit(node)


def extract_numpy_scalars(expr: str):
    # Parse the expression into an AST
    tree = ast.parse(expr, mode="eval")

    # Transform the AST
    transformer = TransformNumpyCalls()
    transformed_tree = transformer.visit(tree)

    # Generate the modified expression
    transformed_expr = ast.unparse(transformed_tree)

    return transformed_expr, transformer.replacements


def _isscalar(arr):
    return np.isscalar(arr) or (hasattr(arr, "shape") and arr.shape == ())


def validate_inputs(inputs: dict, out=None, reduce=False) -> tuple:  # noqa: C901
    """Validate the inputs for the expression."""
    if not inputs:
        if out is None:
            raise ValueError(
                "You really want to pass at least one input or one output for building a LazyArray."
                " Maybe you want blosc2.empty() instead?"
            )
        if isinstance(out, blosc2.NDArray):
            return out.shape, out.chunks, out.blocks, True
        else:
            return out.shape, None, None, True

    raw_inputs = [input_ for input_ in inputs.values() if (input_ is not np and not _isscalar(input_))]
    if not raw_inputs:
        # Scalar-only expressions have scalar output shape but can use miniexpr
        return (), None, None, True

    # This will raise an exception if the input shapes are not compatible
    shape = compute_broadcast_shape(raw_inputs)

    if not all(np.array_equal(shape, input.shape) for input in raw_inputs):
        # If inputs have different shapes (other than scalars), we cannot take the fast path
        return shape, None, None, False

    # More checks specific to NDArray inputs
    # NDInputs are either non-SimpleProxy with chunks or are SimpleProxy with src having chunks
    NDinputs = [
        input
        for input in raw_inputs
        if (hasattr(input, "chunks") and not isinstance(input, blosc2.SimpleProxy))
        or (isinstance(input, blosc2.SimpleProxy) and hasattr(input.src, "chunks"))
    ]
    if not NDinputs:
        # All inputs are NumPy arrays, so we cannot take the fast path
        if raw_inputs and hasattr(raw_inputs[0], "shape"):
            shape = raw_inputs[0].shape
        else:
            shape = None
        return shape, None, None, False

    # Check if we can take the fast path
    # For this we need that the chunks and blocks for all inputs (and a possible output)
    # are the same
    fast_path = True
    first_input = NDinputs[0]
    # Check the out NDArray (if present) first
    if isinstance(out, blosc2.NDArray) and not reduce:
        if first_input.shape != out.shape:
            return None, None, None, False
        if first_input.chunks != out.chunks:
            fast_path = False
        if first_input.blocks != out.blocks:
            fast_path = False
        if 0 in out.chunks:  # fast_eval has zero division error for 0 shapes
            fast_path = False
    # Then, the rest of the operands
    for input_ in NDinputs:
        if first_input.chunks != input_.chunks:
            fast_path = False
        if first_input.blocks != input_.blocks:
            fast_path = False
        if 0 in input_.chunks:  # fast_eval has zero division error for 0 shapes
            fast_path = False

    return first_input.shape, first_input.chunks, first_input.blocks, fast_path


def is_full_slice(item):
    """Check whether the slice represented by item is a full slice."""
    if item == ():
        # This is the case when the user does not pass any slice in compute() method
        return True
    if isinstance(item, tuple):
        return all((isinstance(i, slice) and i == slice(None, None, None)) or i == Ellipsis for i in item)
    elif isinstance(item, int | bool):
        return False
    else:
        return item in (slice(None, None, None), Ellipsis)


def do_slices_intersect(slice1: list | tuple, slice2: list | tuple) -> bool:
    """
    Check whether two slices intersect.

    Parameters
    ----------
    slice1: list of slices
        The first slice
    slice2: list of slices
        The second slice

    Returns
    -------
    bool
        Whether the slices intersect
    """

    # Pad the shorter slice list with full slices (:)
    while len(slice1) < len(slice2):
        slice1.append(slice(None))
    while len(slice2) < len(slice1):
        slice2.append(slice(None))

    # Check each dimension for intersection
    for s1, s2 in zip(slice1, slice2, strict=True):
        if s1 is Ellipsis or s2 is Ellipsis:
            return True
        if s1.start >= s2.stop:
            return False
        if s1.stop <= s2.start:
            return False

    return True


def get_chunk(arr, info, nchunk):
    reduc, aligned, low_mem, chunks_idx = info

    if low_mem:
        # We don't want to uncompress the chunk, so keep it compressed and
        # decompress it just before execution.  This is normally slower, but
        # can be useful in scarce memory situations.
        return arr.schunk.get_chunk(nchunk)

    # First check if the chunk is a special zero chunk.
    # Using lazychunks is very effective here because we only need to read the header.
    if reduc:
        # Reductions can treat zero scalars as zero chunks
        chunk = arr.schunk.get_lazychunk(nchunk)
        special = blosc2.SpecialValue((chunk[31] & 0x70) >> 4)
        if special == blosc2.SpecialValue.ZERO:
            return np.zeros((), dtype=arr.dtype)

    shape, chunks = arr.shape, arr.chunks
    coords = tuple(np.unravel_index(nchunk, chunks_idx))
    slice_ = tuple(
        # slice(c * s, min((c + 1) * s, shape))  # uncomment to make code hang here
        slice(c * s, min((c + 1) * s, shape[i]))
        for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
    )
    chunks_ = tuple(s.stop - s.start for s in slice_)

    if aligned:
        # Decompress the whole chunk and return it
        buff = arr.schunk.decompress_chunk(nchunk)
        bsize = arr.dtype.itemsize * math.prod(chunks_)
        return np.frombuffer(buff[:bsize], dtype=arr.dtype).reshape(chunks_)

    return arr[slice_]


def _stoppable_put(queue, item, stop):
    """Put *item* into the bounded *queue*, giving up when *stop* gets set.

    Returns False when aborted, so the producer can exit instead of blocking
    forever on a queue whose consumer is gone.
    """
    while not stop.is_set():
        try:
            queue.put(item, timeout=0.1)
            return True
        except Full:
            continue
    return False


def read_chunks_worker(arrs, info, queue, stop):
    """Read the chunks of all operands concurrently and feed them into *queue*.

    For each chunk index, the reads are submitted to a thread pool (one task per
    operand) so that file reads and decompression overlap; the bounded queue in
    :func:`sync_read_chunks` provides the prefetch ahead of the consumer.
    """
    shape, chunks_ = arrs[0].shape, arrs[0].chunks
    max_workers = max(1, min(len(arrs), int(getattr(blosc2, "nthreads", 1) or 1)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        my_chunk_iter = range(arrs[0].schunk.nchunks)
        if len(info) == 5:
            if info[-1] is not None:
                my_chunk_iter = sliced_chunk_iter(chunks_, (), shape, axis=info[-1], nchunk=True)
            info = info[:4]
        for i, nchunk in enumerate(my_chunk_iter):
            futures = [executor.submit(get_chunk, arr, info, nchunk) for arr in arrs]
            # result() keeps operand order and propagates the first exception raised
            if not _stoppable_put(queue, (i, [future.result() for future in futures]), stop):
                return

    _stoppable_put(queue, None, stop)  # signal the end of the chunks


def sync_read_chunks(arrs, info):
    queue_size = 2  # maximum number of chunks in the queue
    queue = Queue(maxsize=queue_size)
    # Signals the producer to bail out when the consumer goes away (e.g. an
    # exception during evaluation closes this generator early); without it,
    # the producer can block forever on a full queue and deadlock the
    # thread.join() below during generator finalization.
    stop = threading.Event()
    worker_exc = None

    def _run_reader():
        nonlocal worker_exc
        try:
            read_chunks_worker(arrs, info, queue, stop)
        except BaseException as exc:
            worker_exc = exc
            _stoppable_put(queue, None, stop)

    # Start the file reading in a separate thread
    thread = threading.Thread(target=_run_reader)
    thread.start()

    try:
        # Read the chunks synchronously from the queue
        while True:
            try:
                chunks = queue.get(timeout=1)  # Wait for the next chunk
                if chunks is None:  # End of chunks
                    if worker_exc is not None:
                        raise worker_exc
                    break
                yield chunks
            except Empty:
                if not thread.is_alive():
                    if worker_exc is not None:
                        raise worker_exc from None
                    break
                continue
    finally:
        stop.set()
        thread.join()


def read_nchunk(arrs, info):
    for _, chunks in sync_read_chunks(arrs, info):
        yield chunks


iter_chunks = None


def fill_chunk_operands(
    operands, slice_, chunks_, full_chunk, aligned, nchunk, iter_disk, chunk_operands, reduc=False, axis=None
):
    """Retrieve the chunk operands for evaluating an expression.

    This function provides an optimized path for full chunks and a slower path for partial chunks.
    """
    global iter_chunks

    if iter_disk:
        # Use an environment variable to control the memory usage
        low_mem = os.environ.get("BLOSC_LOW_MEM", False)
        # This method is only useful when all operands are NDArray and shows better
        # performance only when at least one of them is persisted on disk
        if iter_chunks is None:
            # Initialize the iterator for reading the chunks
            # Take any operand (all should have the same shape and chunks)
            key, arr = next(iter(operands.items()))
            chunks_idx, _ = get_chunks_idx(arr.shape, arr.chunks)
            info = (reduc, aligned[key], low_mem, chunks_idx, axis)
            iter_chunks = read_nchunk(list(operands.values()), info)
        # Run the asynchronous file reading function from a synchronous context
        chunks = next(iter_chunks)

        for i, (key, value) in enumerate(operands.items()):
            # Chunks are already decompressed, so we can use them directly
            if not low_mem:
                if full_chunk:
                    chunk_operands[key] = chunks[i]
                else:
                    chunk_operands[key] = value[slice_]
                continue
            # Otherwise, we need to decompress them
            if aligned[key]:
                buff = blosc2.decompress2(chunks[i])
                bsize = value.dtype.itemsize * math.prod(chunks_)
                chunk_operands[key] = np.frombuffer(buff[:bsize], dtype=value.dtype).reshape(chunks_)
            else:
                chunk_operands[key] = value[slice_]
        return

    # Get the starts and stops for the slice
    starts = [s.start if s.start is not None else 0 for s in slice_]
    stops = [s.stop if s.stop is not None else sh for s, sh in zip(slice_, chunks_, strict=True)]

    for key, value in operands.items():
        if np.isscalar(value):
            chunk_operands[key] = value
            continue
        if value.shape == ():
            chunk_operands[key] = value[()]
            continue
        if not full_chunk or not isinstance(value, blosc2.NDArray):
            # The chunk is not a full one, or has padding, or is not a blosc2.NDArray,
            # so we need to go the slow path
            chunk_operands[key] = value[slice_]
            continue

        # If key is in operands, we can reuse the buffer
        if (
            key in chunk_operands
            and chunks_ == chunk_operands[key].shape
            and isinstance(value, blosc2.NDArray)
        ):
            value.get_slice_numpy(chunk_operands[key], (starts, stops))
            continue

        if aligned[key]:
            # Decompress the whole chunk and store it
            buff = value.schunk.decompress_chunk(nchunk)
            bsize = value.dtype.itemsize * math.prod(chunks_)
            chunk_operands[key] = np.frombuffer(buff[:bsize], dtype=value.dtype).reshape(chunks_)
        else:
            chunk_operands[key] = value[slice_]


def _apply_jit_backend_pragma(expression: str, inputs: dict, jit_backend: str | None) -> str:
    if jit_backend is None:
        return expression
    if jit_backend == "js":
        # "js" is handled earlier (DSL kernels -> JS bridge); it never carries a C pragma.
        return expression
    if jit_backend not in ("tcc", "cc"):
        raise ValueError("jit_backend must be one of: None, 'tcc', 'cc', 'js'")

    pragma = f"# me:compiler={jit_backend}\n"
    stripped = expression.lstrip()
    if stripped.startswith("def "):
        if "# me:compiler=" in expression:
            return expression
        return pragma + expression
    params = ", ".join(k for k, v in inputs.items() if hasattr(v, "dtype"))
    return f"{pragma}def __me_auto({params}):\n    return {expression}"


def _inject_dummy_param_for_zero_input_dsl(expression: str, param_name: str) -> str:
    pattern = re.compile(r"^(\s*def\s+[A-Za-z_]\w*)\(\s*\)(\s*:)", re.MULTILINE)
    rewritten, nsubs = pattern.subn(rf"\1({param_name})\2", expression, count=1)
    if nsubs == 0:
        raise ValueError("Could not inject dummy DSL parameter for zero-input kernel")
    return rewritten


def _is_dsl_kernel_expression(expression) -> bool:
    return isinstance(expression, DSLKernel) and expression.dsl_source is not None


def _as_js_udf(expression, shape=None):
    """For jit_backend="js": transpile a DSL kernel to JS and return a plain per-block
    callable (so the normal UDF path runs it). Browser/Pyodide only.

    *shape* (the whole-array output shape) is forwarded to the transpiler so kernels using
    index/shape symbols can reconstruct global coordinates per block."""
    if not _is_dsl_kernel_expression(expression):
        raise ValueError('jit_backend="js" requires a blosc2.dsl_kernel-decorated kernel')
    if not blosc2.IS_WASM:
        raise RuntimeError('jit_backend="js" is only available under WebAssembly/Pyodide')
    from .dsl_js import js_kernel

    return js_kernel(expression, shape=shape)


def _js_dtypes_ok(operands, kwargs) -> bool:
    """True only if the JS bridge (which computes in float64) is safe for these operands.

    The output dtype must be floating: integer/complex *output* goes to miniexpr (the bridge
    can't reproduce integer division/overflow/truncation semantics, and float64 can't hold
    int64 exactly).  Given a floating output, integer *inputs* are fine -- the bridge converts
    every operand to float64, which is exactly what miniexpr does when promoting integer inputs
    for a float result (so any values above 2**53 lose precision identically).  Complex inputs
    are rejected (the bridge is real-only)."""
    dt = kwargs.get("dtype")
    if dt is None:
        # Inferred output: only safe when all operands are float (so the output is float too).
        return all(
            np.issubdtype(op.dtype, np.floating)
            for op in operands.values()
            if isinstance(op, blosc2.NDArray)
        )
    if not np.issubdtype(np.dtype(dt), np.floating):
        return False
    return all(
        np.issubdtype(op.dtype, np.floating) or np.issubdtype(op.dtype, np.integer)
        for op in operands.values()
        if isinstance(op, blosc2.NDArray)
    )


def _maybe_js_backend(expression, jit, jit_backend, reduce_args, operands, kwargs, shape=None):
    """Resolve the JS backend for a DSL kernel.

    - ``jit_backend="js"`` (explicit): transpile to the JS bridge, or raise if it can't.
    - ``jit_backend=None`` under WebAssembly, unless ``jit=False``: *prefer* JS (it is a
      JIT too, and the fastest one here) for transpilable float DSL kernels, silently
      falling back to miniexpr for anything it can't do (non-float dtypes, reductions, or
      unsupported DSL constructs).  ``jit=True`` and ``jit=None`` both prefer JS; only
      ``jit=False`` (interpreter), ``strict_miniexpr=True``, or an explicit ``jit_backend``
      opts out.

    *shape* is the whole-array output shape, forwarded to the transpiler for kernels that
    use index/shape symbols (``_i0``/``_n0``/``_flat_idx``); without it such kernels fall
    back to miniexpr.

    Returns ``(expression, jit, jit_backend)`` — expression becomes a plain per-block
    callable when JS is chosen, else everything passes through unchanged.
    """
    if jit_backend == "js":
        if reduce_args:
            raise ValueError('jit_backend="js" does not support reductions')
        out_dtype = kwargs.get("dtype")
        if out_dtype is not None and not np.issubdtype(np.dtype(out_dtype), np.floating):
            # The JS bridge computes in float64 and cannot reproduce integer/complex output
            # semantics (division/overflow/truncation); keep those on miniexpr.
            raise ValueError(
                'jit_backend="js" requires a floating-point output dtype '
                f"(got {np.dtype(out_dtype)}); drop jit_backend to use miniexpr"
            )
        return _as_js_udf(expression, shape), None, None
    prefer_js = (
        jit is not False  # jit=True/None prefer the best JIT (js); only jit=False forces interpreter
        and jit_backend is None
        and not kwargs.get("strict_miniexpr")  # explicit strict_miniexpr=True keeps miniexpr
        and blosc2.IS_WASM
        and _is_dsl_kernel_expression(expression)
        and operands  # at least one operand: the zero-input DSL path stays on miniexpr
        and not reduce_args
        and _js_dtypes_ok(operands, kwargs)
    )
    if not prefer_js:
        return expression, jit, jit_backend
    try:
        bridge = _as_js_udf(expression, shape)  # transpiles; raises on any unsupported construct
    except Exception:
        return expression, jit, jit_backend  # fall back to miniexpr, no regression
    return bridge, None, None


def _format_dsl_parse_error_hint(expr_text: str, backend_msg: str):
    marker = "parse_error_pos="
    pos0 = backend_msg.find(marker)
    if pos0 < 0:
        return None
    pos0 += len(marker)
    pos1 = pos0
    while pos1 < len(backend_msg) and backend_msg[pos1].isdigit():
        pos1 += 1
    if pos1 == pos0:
        return None
    err_pos = int(backend_msg[pos0:pos1])
    if err_pos < 0:
        return None
    if err_pos > len(expr_text):
        err_pos = len(expr_text)
    line_no = expr_text.count("\n", 0, err_pos) + 1
    line_start = expr_text.rfind("\n", 0, err_pos) + 1
    col_no = err_pos - line_start + 1
    dump = DSLValidator(expr_text)._format_source_with_pointer(line_no, col_no)
    return f"Parse error location (line {line_no}, col {col_no}, offset {err_pos}):\n{dump}"


def _dsl_miniexpr_required_message(reason: str | None = None) -> str:
    message = ""
    if reason:
        message = f"{message}{reason}"
    return message


def _raise_dsl_miniexpr_required(reason: str | None = None) -> None:
    raise RuntimeError(_dsl_miniexpr_required_message(reason))


def fast_eval(  # noqa: C901
    expression: str | Callable[[tuple, np.ndarray, tuple[int]], None],
    operands: dict,
    getitem: bool,
    **kwargs,
) -> blosc2.NDArray | np.ndarray:
    """Evaluate the expression in chunks of operands using a fast path.

    Parameters
    ----------
    expression: str or callable
        The expression or udf to evaluate.
    operands: dict
        A dictionary containing the operands for the expression.
    getitem: bool, optional
        Indicates whether the expression is being evaluated for a getitem operation or compute().
        Default is False.
    kwargs: Any, optional
        Additional keyword arguments supported by the :func:`empty` constructor.

    Returns
    -------
    :ref:`NDArray` or np.ndarray
        The output array.
    """
    global try_miniexpr

    # Use a local copy so we don't modify the global
    use_miniexpr = try_miniexpr

    is_dsl = _is_dsl_kernel_expression(expression)
    expr_string = expression.dsl_source if is_dsl else expression
    dsl_disable_reason = None

    # Disable miniexpr for UDFs (callable expressions), except DSL kernels
    if callable(expression) and not is_dsl:
        use_miniexpr = False

    out = kwargs.pop("_output", None)
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    fp_accuracy = kwargs.pop("fp_accuracy", blosc2.FPAccuracy.DEFAULT)
    jit = kwargs.pop("jit", None)
    jit_backend = kwargs.pop("jit_backend", None)
    strict_miniexpr = kwargs.pop("strict_miniexpr", None)
    _in_place = kwargs.pop("in_place", None)
    dtype = kwargs.pop("dtype", None)
    requested_shape = kwargs.pop("shape", None)
    where: dict | None = kwargs.pop("_where_args", None)
    # Optional candidate-block bitmap (1-D miniexpr only): uint8, one byte per
    # global block (0 = skip → false output, non-zero = evaluate).  Used by the
    # CTable SUMMARY-index path to skip non-candidate blocks inside miniexpr.
    candidate_blocks = kwargs.pop("_candidate_blocks", None)
    # Opt-in: skip update_data for chunks that contain no candidate block.  The
    # result array is left uninitialized in those chunks, so the caller MUST read
    # only candidate chunks (scoped extraction).  Used by the CTable SUMMARY path.
    prune_chunks = kwargs.pop("_prune_chunks", False)
    if strict_miniexpr is None:
        # Be strict by default for DSL kernels to avoid silently losing DSL fast-path regressions.
        strict_miniexpr = bool(is_dsl)
    if where is not None and len(where) != 2:
        # miniexpr does not support cardinality-changing where (len==1);
        # where(cond, x, y) with two args is element-wise and IS supported.
        use_miniexpr = False
        if is_dsl:
            dsl_disable_reason = "DSL kernels cannot be run without miniexpr."
    if isinstance(out, blosc2.NDArray):
        # If 'out' has been passed, and is a NDArray, use it as the base array
        basearr = out
    elif isinstance(out, np.ndarray):
        # If 'out' is a NumPy array, create a NDArray with the same shape and dtype
        basearr = blosc2.empty(out.shape, dtype=out.dtype, **kwargs)
    else:
        # Otherwise, find the operand with the 'chunks' attribute and the longest shape
        operands_with_chunks = [o for o in operands.values() if hasattr(o, "chunks")]
        if operands_with_chunks and any(arr.chunks != () for arr in operands_with_chunks):
            basearr = max(operands_with_chunks, key=lambda x: len(x.shape))
        else:
            if requested_shape is None:
                raise ValueError("Cannot infer output shape without operands; pass `shape` explicitly")
            if dtype is None:
                raise ValueError("Cannot infer output dtype without operands; pass `dtype` explicitly")
            basearr = blosc2.empty(
                requested_shape,
                dtype=dtype,
                chunks=kwargs.get("chunks"),
                blocks=kwargs.get("blocks"),
            )

    # Get the shape of the base array
    shape = basearr.shape
    chunks = kwargs.pop("chunks", None)
    if chunks is None:
        chunks = basearr.chunks
    blocks = kwargs.pop("blocks", None)
    if blocks is None:
        blocks = basearr.blocks
    # Check whether the partitions are aligned and behaved
    aligned = {
        k: False if not hasattr(k, "chunks") else blosc2.are_partitions_aligned(k.shape, k.chunks, k.blocks)
        for k in operands
    }
    behaved = blosc2.are_partitions_behaved(shape, chunks, blocks)

    # Check that all operands are NDArray for fast path
    all_ndarray = all(isinstance(value, blosc2.NDArray) and value.shape != () for value in operands.values())
    # Check that there is some NDArray that is persisted in the disk
    any_persisted = any(
        (isinstance(value, blosc2.NDArray) and value.shape != () and value.schunk.urlpath is not None)
        for value in operands.values()
    )
    if not blosc2.IS_WASM:
        iter_disk = all_ndarray and any_persisted
    else:
        # WebAssembly does not support threading, so we cannot use the iter_disk option
        iter_disk = False

    expr_string_miniexpr = expr_string
    operands_miniexpr = operands
    if use_miniexpr and isinstance(expr_string, str):
        try:
            expr_string_miniexpr, operands_miniexpr = specialize_miniexpr_inputs(expr_string, operands)
        except Exception:
            # If specialization fails, keep original expression/operands and let normal checks decide.
            expr_string_miniexpr = expr_string
            operands_miniexpr = operands

    # Check whether we can use miniexpr
    if use_miniexpr:
        if is_dsl and isinstance(expression, DSLKernel) and not operands_miniexpr:
            # Scalar specialization may remove all kernel inputs at runtime (e.g. `f(start)` with start=3),
            # so inject a dummy array operand for miniexpr even if the original DSL signature had parameters.
            dummy_name = "__me_dummy0"
            dummy_dtype = dtype if dtype is not None else np.uint8
            expr_string_miniexpr = _inject_dummy_param_for_zero_input_dsl(expr_string_miniexpr, dummy_name)
            operands_miniexpr = {
                dummy_name: blosc2.zeros(shape, dtype=dummy_dtype, chunks=chunks, blocks=blocks)
            }
        if math.prod(shape) <= 1:
            # Avoid miniexpr for scalar-like outputs; current prefilter path is unstable here.
            use_miniexpr = False
            if is_dsl and dsl_disable_reason is None:
                dsl_disable_reason = "scalar-like outputs are not supported by the DSL miniexpr path."
        if (
            isinstance(expr_string_miniexpr, str)
            and
            # Prefix scans are stateful across chunks and not safe for miniexpr prefilter execution.
            any(tok in expr_string_miniexpr for tok in ("cumsum(", "cumprod(", "cumulative_sum("))
        ):
            use_miniexpr = False
            if is_dsl and dsl_disable_reason is None:
                dsl_disable_reason = "cumulative scans are not supported by the DSL miniexpr path."
        if isinstance(expr_string_miniexpr, str):
            expr_string_miniexpr = _apply_jit_backend_pragma(
                expr_string_miniexpr, operands_miniexpr, jit_backend
            )
        all_ndarray_miniexpr = all(
            isinstance(value, blosc2.NDArray) and value.shape != () for value in operands_miniexpr.values()
        )
        # Require aligned NDArray operands with identical chunk/block grid.
        same_shape = all(hasattr(op, "shape") and op.shape == shape for op in operands_miniexpr.values())
        same_chunks = all(hasattr(op, "chunks") and op.chunks == chunks for op in operands_miniexpr.values())
        same_blocks = all(hasattr(op, "blocks") and op.blocks == blocks for op in operands_miniexpr.values())
        if not (same_shape and same_chunks and same_blocks):
            use_miniexpr = False
            if is_dsl and dsl_disable_reason is None:
                dsl_disable_reason = "all DSL operands must share shape/chunks/blocks."
        if not (all_ndarray_miniexpr and out is None):
            use_miniexpr = False
            if is_dsl and dsl_disable_reason is None:
                dsl_disable_reason = (
                    "DSL kernels require NDArray inputs and do not support the `out` argument."
                )
        has_complex = any(
            isinstance(op, blosc2.NDArray) and blosc2.isdtype(op.dtype, "complex floating")
            for op in operands_miniexpr.values()
        )
        if isinstance(expr_string_miniexpr, str) and has_complex:
            if sys.platform == "win32" or blosc2.IS_WASM:
                # On Windows and WebAssembly, miniexpr has issues with complex numbers
                use_miniexpr = False
                if is_dsl and dsl_disable_reason is None:
                    dsl_disable_reason = "complex DSL kernels are disabled on Windows and WebAssembly."
            if any(tok in expr_string_miniexpr for tok in ("!=", "==", "<=", ">=", "<", ">")):
                use_miniexpr = False
                if is_dsl and dsl_disable_reason is None:
                    dsl_disable_reason = "complex comparisons are not supported by miniexpr."

    if is_dsl and not use_miniexpr:
        _raise_dsl_miniexpr_required(dsl_disable_reason)

    if os.environ.get("BLOSC_ME_JIT_TRACE", "").lower() in ("1", "true", "on"):
        engine = (
            "miniexpr" if use_miniexpr else ("ne_evaluate" if isinstance(expr_string, str) else "python-udf")
        )
        jit_info = f"jit={jit}, backend={jit_backend}" if use_miniexpr else ""
        expr_short = str(expr_string)[:120].replace("\n", " ")
        print(f"[blosc2] engine={engine} {jit_info} expr={expr_short}", flush=True)

    if use_miniexpr:
        cparams = kwargs.pop("cparams", blosc2.CParams())
        # All values will be overwritten, so we can use an uninitialized array
        res_eval = blosc2.uninit(shape, dtype, chunks=chunks, blocks=blocks, cparams=cparams, **kwargs)
        prefilter_set = False
        try:
            # Fuse where(cond, x, y) into the expression for miniexpr
            _pref_expr = expr_string_miniexpr
            _pref_ops = operands_miniexpr
            if where is not None and len(where) == 2:
                _pref_expr = f"where({_pref_expr}, _where_x, _where_y)"
            # _cb_anchor keeps the contiguous bitmap alive for the whole
            # update_data loop (the prefilter dereferences it per block).  Pass
            # candidate_blocks only when present so external wrappers of
            # _set_pref_expr with the legacy signature keep working.
            _pref_kwargs = {"fp_accuracy": fp_accuracy, "jit": jit}
            if candidate_blocks is not None:
                _pref_kwargs["candidate_blocks"] = candidate_blocks
            _cb_anchor = res_eval._set_pref_expr(_pref_expr, _pref_ops, **_pref_kwargs)
            prefilter_set = True
            # print("expr->miniexpr:", expression, fp_accuracy)
            # Data to compress is fetched from operands, so it can be uninitialized here
            data = np.empty(res_eval.schunk.chunksize, dtype=np.uint8)
            # Exercise prefilter for each chunk.  With a candidate bitmap the
            # prefilter zeroes non-candidate blocks without decompressing inputs
            # or running miniexpr, which is where the savings come from.  With
            # _prune_chunks, whole chunks that contain no candidate block are
            # skipped entirely (left uninitialized — caller reads only candidates).
            chunk_has_candidate = None
            if prune_chunks and _cb_anchor is not None and len(chunks) == 1 and len(blocks) == 1:
                bpc = -(-int(chunks[0]) // int(blocks[0]))  # blocks per chunk (ceil)
                nct = res_eval.schunk.nchunks
                if bpc > 0 and _cb_anchor.shape[0] >= nct * bpc:
                    chunk_has_candidate = _cb_anchor[: nct * bpc].reshape(nct, bpc).any(axis=1)
            for nchunk in range(res_eval.schunk.nchunks):
                if chunk_has_candidate is not None and not chunk_has_candidate[nchunk]:
                    continue
                res_eval.schunk.update_data(nchunk, data, copy=False)
            # _cb_anchor (if any) stays referenced until fast_eval returns, which
            # is after the finally below removes the prefilter — so the buffer
            # safely outlives every prefilter dereference.
        except Exception as e:
            use_miniexpr = False
            if is_dsl:
                reason = (
                    f"miniexpr compilation or execution failed for this DSL kernel:\n{expression.dsl_source}"
                )
                backend_error = str(e)
                parse_hint = None
                if isinstance(expr_string_miniexpr, str):
                    parse_hint = _format_dsl_parse_error_hint(expr_string_miniexpr, backend_error)
                reason = f"{reason}\nBackend error: {backend_error}"
                if parse_hint is not None:
                    reason = f"{reason}\n{parse_hint}"
                raise RuntimeError(_dsl_miniexpr_required_message(reason)) from e
            if strict_miniexpr:
                raise RuntimeError("miniexpr evaluation failed while strict_miniexpr=True") from e
        finally:
            if prefilter_set:
                res_eval.schunk.remove_prefilter("miniexpr")
            global iter_chunks
            # Ensure any background reading thread is closed
            iter_chunks = None

        if not use_miniexpr:
            # If miniexpr failed, fallback to regular evaluation
            # (continue to the manual chunked evaluation below)
            pass
        else:
            if getitem:
                return res_eval[()]
            return res_eval

    chunk_operands = {}
    # Check which chunks intersect with _slice
    all_chunks = get_intersecting_chunks((), shape, chunks)  # if _slice is (), returns all chunks
    for nchunk, chunk_slice in enumerate(all_chunks):
        cslice = chunk_slice.raw
        offset = tuple(s.start for s in cslice)  # offset for the udf
        chunks_ = tuple(s.stop - s.start for s in cslice)

        full_chunk = chunks_ == chunks  # slice is same as chunk
        fill_chunk_operands(
            operands, cslice, chunks_, full_chunk, aligned, nchunk, iter_disk, chunk_operands
        )

        # Since ne_evaluate() can return a dtype larger than the one in computed in the expression,
        # we cannot take this fast path
        # if isinstance(out, np.ndarray) and not where:
        #     # Fast path: put the result straight in the output array (avoiding a memory copy)
        #     if callable(expression):
        #         expression(tuple(chunk_operands.values()), out[slice_], offset=offset)
        #     else:
        #         ne_evaluate(expression, chunk_operands, out=out[slice_])
        #     continue
        if out is None:
            # We can enter here when using any of the compute() or __getitem__() methods
            if getitem:
                out = np.empty(shape, dtype=dtype)
            else:
                out = blosc2.empty(shape, chunks=chunks, blocks=blocks, dtype=dtype, **kwargs)

        if callable(expression):
            if _is_dsl_kernel_expression(expression):
                _raise_dsl_miniexpr_required(
                    "internal fallback attempted to execute the DSL kernel directly in Python."
                )
            if _in_place:
                expression(tuple(chunk_operands.values()), out, offset=offset)
                continue
            result = np.empty(chunks_, dtype=out.dtype)
            expression(tuple(chunk_operands.values()), result, offset=offset)
        else:
            if where is None:
                result = ne_evaluate(expression, chunk_operands, **ne_args)
            else:
                # Apply the where condition (in result)
                if len(where) == 2:
                    new_expr = f"where({expression}, _where_x, _where_y)"
                    result = ne_evaluate(new_expr, chunk_operands, **ne_args)
                else:
                    # We do not support one or zero operands in the fast path yet
                    raise ValueError("Fast path: the where condition must be a tuple with two elements")

        # Store the result in the output array
        if getitem:
            try:
                out[cslice] = result
            except ComplexWarning:
                # The result is a complex number, so we need to convert it to real.
                # This is a workaround for rigidness of NumExpr with type casting.
                result = result.real.astype(out.dtype)
                out[cslice] = result
        else:
            if behaved and result.shape == chunks_ and result.dtype == out.dtype:
                # Fast path only works for results that are full chunks
                out.schunk.update_data(nchunk, result, copy=False)
            else:
                out[cslice] = result

    return out


def compute_start_index(shape, slice_obj):
    """
    Compute the index of the starting element of a slice in an n-dimensional array.

    Parameters
    ----------
    shape : tuple
        The shape of the n-dimensional array.
    slice_obj : tuple of slices
        The slice object representing the slice of the array.

    Returns
    -------
    start_index : int
        The index of the starting element of the slice.
    """
    if not isinstance(slice_obj, tuple):
        slice_obj = (slice_obj,)

    start_index = 0
    stride = 1

    for dim, sl in reversed(list(enumerate(slice_obj))):
        if isinstance(sl, slice):
            start = sl.start if sl.start is not None else 0
        elif sl is Ellipsis:
            start = 0
        else:
            start = sl

        start_index += start * stride
        stride *= shape[dim]

    return start_index


def slices_eval(  # noqa: C901
    expression: str | Callable[[tuple, np.ndarray, tuple[int]], None],
    operands: dict,
    getitem: bool,
    _slice=NDINDEX_EMPTY_TUPLE,
    shape=None,
    **kwargs,
) -> blosc2.NDArray | np.ndarray:
    """Evaluate the expression in chunks of operands.

    This function can handle operands with different chunk shapes and
    can evaluate only a slice of the output array if needed.

    This is also flexible enough to work with operands of different shapes.

    Parameters
    ----------
    expression: str or callable
        The expression or user-defined (udf) to evaluate.
    operands: dict
        A dictionary containing the operands for the expression.
    getitem: bool, optional
        Indicates whether the expression is being evaluated for a getitem operation or compute().
        Default is False.
    _slice: ndindex.Tuple sequence of slices and ints. Default = ndindex.Tuple(), optional
        If provided, only the chunks that intersect with this slice
        will be evaluated.
    shape: tuple | None
        The shape of the full (unsliced result). Typically passed on from parent LazyArray.
        If None, a guess is made from broadcasting the operands.
    kwargs: Any, optional
        Additional keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    :ref:`NDArray` or np.ndarray
        The output array.
    """
    out: blosc2.NDArray | None = kwargs.pop("_output", None)
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    chunks = kwargs.get("chunks")
    where: dict | None = kwargs.pop("_where_args", None)
    use_index = kwargs.pop("_use_index", True)
    _indices = kwargs.pop("_indices", False)
    if _indices and (not where or len(where) != 1):
        raise NotImplementedError("Indices can only be used with one where condition")
    _order = kwargs.pop("_order", None)
    if _order is not None and not isinstance(_order, list):
        # Always use a list for _order
        _order = [_order]

    dtype = kwargs.pop("dtype", None)
    _in_place = kwargs.pop("in_place", False)
    shape_slice = None
    need_final_slice = False

    # keep orig_slice
    _slice = _slice.raw
    orig_slice = _slice

    # Compute the shape and chunks of the output array, including broadcasting
    if shape is None:  # lazyudf provides shape kwarg
        shape = compute_broadcast_shape(operands.values())

    if _slice != ():
        # Check whether _slice contains an integer, or any step that are not None or 1
        if any((isinstance(s, int)) for s in _slice):
            need_final_slice = True
        _slice = tuple(slice(i, i + 1, 1) if isinstance(i, int) else i for i in _slice)
        # shape_slice in general not equal to final shape:
        # dummy dims (due to ints) will be dealt with by taking final_slice
        shape_slice = ndindex.ndindex(_slice).newshape(shape)
        mask_slice = np.array([isinstance(i, int) for i in orig_slice], dtype=np.bool_)
    if out is not None:
        shape_ = shape_slice if shape_slice is not None else shape
        if shape_ != out.shape and not _in_place:
            raise ValueError("Provided output shape does not match the slice shape.")

    if chunks is None:  # Guess chunk shape
        # Either out, or operand with `chunks`, can be used to get the chunks
        operands_ = [o for o in operands.values() if hasattr(o, "chunks") and o.shape == shape]
        if out is not None and hasattr(out, "chunks"):
            chunks = out.chunks
        elif len(operands_) > 0:
            # Use the first operand with chunks to get the necessary chunking information
            chunks = operands_[0].chunks
        else:
            # Typically, we enter here when using UDFs, and out is a NumPy array.
            # Use operands to get the shape and chunks
            # operand will be a 'fake' NDArray just to get the necessary chunking information
            fp_accuracy = kwargs.pop("fp_accuracy", None)
            temp = blosc2.empty(shape, dtype=dtype)
            if fp_accuracy is not None:
                kwargs["fp_accuracy"] = fp_accuracy
            chunks = temp.chunks
            del temp

    # The starting point for the indices of the inputs
    leninputs = compute_start_index(shape, orig_slice) if orig_slice != () else 0
    lenout = 0
    behaved = False
    indices_ = None
    chunk_indices = None
    dtype_ = np.int64 if _indices else dtype
    if _order is not None:
        # Get the dtype of the array to sort
        dtype_ = operands["_where_x"].dtype
        # Now, use only the fields that are necessary for the sorting
        if dtype_.fields is not None and all(f in dtype_.fields for f in _order):
            dtype_ = np.dtype([(f, dtype_[f]) for f in _order])
        else:
            dtype_ = np.dtype(np.int64)

    # Iterate over the operands and get the chunks
    chunk_operands = {}
    # Check which chunks intersect with _slice (handles zero chunks internally)
    intersecting_chunks = get_intersecting_chunks(
        _slice, shape, chunks
    )  # if _slice is (), returns all chunks
    ratio = (
        np.ceil(np.asarray(shape) / np.asarray(chunks)).astype(np.int64)
        if 0 not in chunks
        else np.asarray(shape)
    )
    index_plan = None
    if where is not None and len(where) == 1 and use_index and _slice == ():
        from . import indexing

        _cache_array = where["_where_x"]
        _cache_tokens = [indexing.SELF_TARGET_NAME]

        # --- Ordered path ---
        if _order is not None:
            ordered_plan = indexing.plan_ordered_query(expression, operands, where, _order)
            if ordered_plan.usable:
                cached_coords = indexing.get_cached_coords(_cache_array, expression, _cache_tokens, _order)
                if cached_coords is not None:
                    return cached_coords
                ordered_positions = indexing.ordered_query_indices(expression, operands, where, _order)
                if ordered_positions is not None:
                    indexing.store_cached_coords(
                        _cache_array, expression, _cache_tokens, _order, ordered_positions
                    )
                    return ordered_positions
            elif indexing.is_expression_order(where["_where_x"], _order):
                raise ValueError("expression order requires a matching full expression index")

        # --- Argsort-only path (.argsort().compute()) ---
        if _indices and _order is None:
            cached_coords = indexing.get_cached_coords(_cache_array, expression, _cache_tokens, None)
            if cached_coords is not None:
                return cached_coords

        # --- Value-returning path (arr[cond][:]) — cache check before plan_query ---
        _cache_urlpath = getattr(_cache_array, "urlpath", None) or getattr(
            getattr(_cache_array, "ndarr", None), "urlpath", None
        )
        if not _indices and _order is None:
            cached_coords = indexing.get_cached_coords(_cache_array, expression, _cache_tokens, None)
            if cached_coords is not None:
                cached_plan = indexing.IndexPlan(
                    usable=True, reason="cache-hit", base=_cache_array, exact_positions=cached_coords
                )
                return indexing.evaluate_full_query(where, cached_plan)

        index_plan = indexing.plan_query(expression, operands, where, use_index=use_index)

        if _indices and _order is None and index_plan.usable:
            if index_plan.exact_positions is not None:
                coords = np.asarray(index_plan.exact_positions, dtype=np.int64)
                indexing.store_cached_coords(_cache_array, expression, _cache_tokens, None, coords)
                return coords
            if index_plan.bucket_masks is not None:
                _, coords = indexing.evaluate_bucket_query(
                    expression, operands, ne_args, where, index_plan, return_positions=True
                )
                indexing.store_cached_coords(_cache_array, expression, _cache_tokens, None, coords)
                return coords
            if index_plan.candidate_units is not None and index_plan.segment_len is not None:
                _, coords = indexing.evaluate_segment_query(
                    expression, operands, ne_args, where, index_plan, return_positions=True
                )
                indexing.store_cached_coords(_cache_array, expression, _cache_tokens, None, coords)
                return coords
        if index_plan.usable and not (_indices or _order):
            if index_plan.exact_positions is not None:
                coords = np.asarray(index_plan.exact_positions, dtype=np.int64)
                indexing.store_cached_coords(_cache_array, expression, _cache_tokens, None, coords)
                return indexing.evaluate_full_query(where, index_plan)
            if index_plan.bucket_masks is not None:
                result, coords = indexing.evaluate_bucket_query(
                    expression, operands, ne_args, where, index_plan, return_positions=True
                )
                indexing.store_cached_coords(_cache_array, expression, _cache_tokens, None, coords)
                return result
            if index_plan.candidate_units is not None and index_plan.segment_len is not None:
                result, coords = indexing.evaluate_segment_query(
                    expression, operands, ne_args, where, index_plan, return_positions=True
                )
                indexing.store_cached_coords(_cache_array, expression, _cache_tokens, None, coords)
                return result

    for chunk_slice in intersecting_chunks:
        # Check whether current cslice intersects with _slice
        cslice = chunk_slice.raw
        nchunk = (
            builtins.sum(c.start // chunks[i] * np.prod(ratio[i + 1 :]) for i, c in enumerate(cslice))
            if 0 not in chunks
            else 0
        )
        if cslice != () and _slice != ():
            # get intersection of chunk and target
            cslice = step_handler(cslice, _slice)
        offset = tuple(s.start for s in cslice)  # offset for the udf
        cslice_shape = tuple(s.stop - s.start for s in cslice)
        len_chunk = math.prod(cslice_shape)
        if (
            index_plan is not None
            and index_plan.usable
            and index_plan.level == "chunk"
            and not index_plan.candidate_units[nchunk]
        ):
            if _indices or _order:
                leninputs += len_chunk
            continue
        # get local index of part of out that is to be updated
        cslice_subidx = (
            ndindex.ndindex(cslice).as_subindex(_slice).raw
        )  # in the case _slice=(), just gives cslice

        get_chunk_operands(operands, cslice, chunk_operands, shape)

        if out is None:
            shape_ = shape_slice if shape_slice is not None else shape
            if where is not None and len(where) < 2:
                # The result is a linear array
                shape_ = math.prod(shape_)
            if getitem or _order:
                out = np.empty(shape_, dtype=dtype_)
                if _order:
                    indices_ = np.empty(shape_, dtype=np.int64)
            else:
                # if "chunks" not in kwargs and (where is None or len(where) == 2):
                # Let's use the same chunks as the first operand (it could have been automatic too)
                # out = blosc2.empty(shape_, chunks=chunks, dtype=dtype_, **kwargs)
                # out = blosc2.empty(shape_, dtype=dtype_, **kwargs)
                if "chunks" in kwargs and (where is not None and len(where) < 2 and len(shape_) > 1):
                    # Remove the chunks argument if the where condition is not a tuple with two elements
                    kwargs.pop("chunks")
                fp_accuracy = kwargs.pop("fp_accuracy", None)
                out = blosc2.empty(shape_, dtype=dtype_, **kwargs)
                if fp_accuracy is not None:
                    kwargs["fp_accuracy"] = fp_accuracy
                # Check if the in out partitions are well-behaved (i.e. no padding)
                behaved = blosc2.are_partitions_behaved(out.shape, out.chunks, out.blocks)
        # Evaluate the expression using chunks of operands

        if callable(expression):
            if _is_dsl_kernel_expression(expression):
                _raise_dsl_miniexpr_required(
                    "internal sliced fallback attempted to execute the DSL kernel directly in Python."
                )
            if _in_place:  # presumably the user knows what they're doing
                # edit out in-place
                expression(tuple(chunk_operands.values()), out, offset=offset)
            else:
                result = np.empty(cslice_shape, dtype=out.dtype)  # raises error if out is None
                # cslice should be equal to cslice_subidx
                # Call the udf directly and use result as the output array
                expression(tuple(chunk_operands.values()), result, offset=offset)
                out[cslice_subidx] = result
            continue

        if _indices or _order:
            indices = np.arange(leninputs, leninputs + len_chunk, dtype=np.int64).reshape(cslice_shape)
            leninputs += len_chunk
            result, chunk_indices = _get_result(expression, chunk_operands, ne_args, where, indices, _order)
        else:
            result, _ = _get_result(expression, chunk_operands, ne_args, where)
        # Enforce contiguity of result (necessary to fill the out array)
        # but avoid copy if already contiguous
        result = np.require(result, requirements="C")

        if where is None or len(where) == 2:
            if behaved and result.shape == out.chunks and result.dtype == out.dtype:
                # Fast path: only use it when the output chunk index is valid
                # (operand and output may have different chunk layouts when slicing)
                if nchunk < out.schunk.nchunks:
                    out.schunk.update_data(nchunk, result, copy=False)
                else:
                    out[cslice_subidx] = result
            else:
                try:
                    out[cslice_subidx] = result
                except ComplexWarning:
                    # The result is a complex number, so we need to convert it to real.
                    # This is a workaround for rigidness of numpy with type casting.
                    result = result.real.astype(out.dtype)
                    out[cslice_subidx] = result
        elif len(where) == 1:
            lenres = len(result)
            out[lenout : lenout + lenres] = result
            if _order is not None:
                indices_[lenout : lenout + lenres] = chunk_indices
            lenout += lenres
        else:
            raise ValueError("The where condition must be a tuple with one or two elements")

    if where is not None and len(where) < 2:  # Don't need to take final_slice since filled up from 0 index
        if _order is not None:
            # argsort the result following _order
            new_order = np.argsort(out[:lenout])
            # And get the corresponding indices in array
            out = indices_[new_order]
        # Cap the output array to the actual length
        if isinstance(out, np.ndarray):
            out = out[:lenout]
        else:
            out.resize((lenout,))

    else:  # Need to take final_slice since filled up array according to slice_ for each chunk
        if need_final_slice:  # only called if out was None
            if isinstance(out, np.ndarray):
                squeeze_axis = np.where(mask_slice)[0]
                squeeze_axis = np.squeeze(squeeze_axis)  # handle 1d mask_slice
                out = np.squeeze(out, squeeze_axis)
            elif isinstance(out, blosc2.NDArray):
                # It *seems* better to choose an automatic chunks and blocks for the output array
                # out = out.slice(_slice, chunks=out.chunks, blocks=out.blocks)
                out = out.squeeze(np.where(mask_slice)[0])
            else:
                raise ValueError("The output array is not a NumPy array or a NDArray")

    return out


def slices_eval_getitem(
    expression: str,
    operands: dict,
    _slice=NDINDEX_EMPTY_TUPLE,
    **kwargs,
) -> np.ndarray:
    """Evaluate the expression in slices of operands.

    This function can handle operands with different chunk shapes and
    can evaluate only a slice of the output array if needed.

    This is a special (and much simplified) version of slices_eval() that
    only works for the case we are returning a NumPy array, where is
    either None or has two args, and expression is not callable.

    One inconvenient of this function is that it tries to evaluate
    the whole slice in one go.  For small slices, this is good, as it
    is normally way more efficient.  However, for larger slices this
    can require large amounts of memory per operand.

    Parameters
    ----------
    expression: str or callable
        The expression or user-defined (udf) to evaluate.
    operands: dict
        A dictionary containing the operands for the expression.
    _slice: ndindex.Tuple sequence of slices and ints. Default = ndindex.Tuple(), optional
        If provided, this slice will be evaluated.
    kwargs: Any, optional
        Additional keyword arguments that are supported by the :func:`empty` constructor.

    Returns
    -------
    :ref:`NDArray` or np.ndarray
        The output array.
    """
    out: np.ndarray | None = kwargs.pop("_output", None)
    ne_args: dict = kwargs.pop("_ne_args", {})
    _in_place = kwargs.pop("in_place", False)
    if ne_args is None:
        ne_args = {}
    where: dict | None = kwargs.pop("_where_args", None)

    dtype = kwargs.pop("dtype", None)
    shape = kwargs.pop("shape", None)
    if shape is None:
        if out is None:
            # Compute the shape and chunks of the output array, including broadcasting
            shape = compute_broadcast_shape(operands.values())
        else:
            shape = out.shape

    # compute the shape of the output array
    _slice = _slice.raw
    _slice_bcast = tuple(slice(i, i + 1) if isinstance(i, int) else i for i in _slice)
    slice_shape = ndindex.ndindex(_slice_bcast).newshape(shape)  # includes dummy dimensions

    # Get the slice of each operand
    slice_operands = {}
    for key, value in operands.items():
        if np.isscalar(value):
            slice_operands[key] = value
            continue
        if value.shape == ():
            slice_operands[key] = value[()]
            continue
        if check_smaller_shape(value.shape, shape, slice_shape, _slice_bcast):
            # We need to fetch the part of the value that broadcasts with the operand
            smaller_slice = compute_smaller_slice(shape, value.shape, _slice)
            slice_operands[key] = value[smaller_slice]
            continue

        slice_operands[key] = value[_slice]

    # Evaluate the expression using slices of operands
    if callable(expression):
        if _is_dsl_kernel_expression(expression):
            _raise_dsl_miniexpr_required(
                "internal getitem fallback attempted to execute the DSL kernel directly in Python."
            )
        offset = tuple(0 if s is None else s.start for s in _slice_bcast)  # offset for the udf
        if _in_place:
            expression(tuple(slice_operands.values()), out, offset=offset)
            return out
        else:
            result = np.empty(slice_shape, dtype=dtype)
            expression(tuple(slice_operands.values()), result, offset=offset)
    else:
        result, _ = _get_result(expression, slice_operands, ne_args, where)

    if out is None:  # avoid copying unnecessarily
        try:
            return result.astype(dtype, copy=False)
        except ComplexWarning:
            # The result is a complex number, so we need to convert it to real.
            # This is a workaround for rigidness of numpy with type casting.
            return result.real.astype(dtype, copy=False)
    else:
        # out should always have maximal shape
        out[_slice] = result
        return out


def infer_reduction_dtype(dtype, operation):
    # It may change in the future, but mostly array-api compliant
    my_float = np.result_type(
        dtype, np.float32 if dtype in (np.float32, np.complex64) else blosc2.DEFAULT_FLOAT
    )
    if operation in {ReduceOp.SUM, ReduceOp.PROD, ReduceOp.CUMULATIVE_SUM, ReduceOp.CUMULATIVE_PROD}:
        if np.issubdtype(dtype, np.bool_):
            return np.int64
        if np.issubdtype(dtype, np.unsignedinteger):
            return np.result_type(dtype, np.uint64)
        return np.result_type(dtype, np.int64 if np.issubdtype(dtype, np.integer) else my_float)
    elif operation in {ReduceOp.MEAN, ReduceOp.STD, ReduceOp.VAR}:
        return my_float
    elif operation in {ReduceOp.MIN, ReduceOp.MAX}:
        return dtype
    elif operation in {ReduceOp.ANY, ReduceOp.ALL}:
        return np.bool_
    elif operation in {ReduceOp.ARGMAX, ReduceOp.ARGMIN}:
        return np.int64
    else:
        raise ValueError(f"Unsupported operation: {operation}")


def step_handler(cslice, _slice):
    out = ()
    for s1, s2 in zip(cslice, _slice, strict=True):
        s1start, s1stop = s1.start, s1.stop
        s2start, s2stop, s2step = s2.start, s2.stop, s2.step
        # assume s1step = 1
        newstart = builtins.max(s1start, s2start)
        newstop = builtins.min(s1stop, s2stop)
        rem = (newstart - s2start) % s2step
        if rem != 0:  # only pass through here if s2step is not 1
            newstart += s2step - rem
            # true_stop = start + n*step + 1 -> stop = start + n * step + 1 + residual
            # so n = (stop - start - 1) // step
            newstop = newstart + (newstop - newstart - 1) // s2step * s2step + 1
        out += (slice(newstart, newstop, s2step),)
    return out


def reduce_slices(  # noqa: C901
    expression: str | Callable[[tuple, np.ndarray, tuple[int]], None],
    operands: dict,
    reduce_args,
    _slice=NDINDEX_EMPTY_TUPLE,
    **kwargs,
) -> blosc2.NDArray | np.ndarray:
    """Evaluate the expression in chunks of operands.

    This function can handle operands with different chunk shapes.
    Also, it can be used when only a slice of the output array is needed.

    Parameters
    ----------
    expression: str or callable
        The expression or user-defined function (udf) to evaluate.
    operands: dict
        A dictionary containing the operands for the operands.
    reduce_args: dict
        A dictionary with arguments to be passed to the reduction function.
    _slice: ndindex.Tuple sequence of slices and ints. Default = ndindex.Tuple(), optional
        If provided, only the chunks that intersect with this slice
        will be evaluated.
    kwargs: Any, optional
        Additional keyword arguments supported by the :func:`empty` constructor.

    Returns
    -------
    :ref:`NDArray` or np.ndarray
        The resulting output array.
    """
    global try_miniexpr

    # Use a local copy so we don't modify the global
    use_miniexpr = try_miniexpr  # & False

    if blosc2.IS_WASM:
        # Reduction miniexpr on wasm is currently unstable for scalar reductions (axis=None).
        # Keep wasm reduction evaluation on the regular chunked path until stabilized.
        use_miniexpr = False

    out = kwargs.pop("_output", None)
    res_out_ = None  # temporary required to store max/min for argmax/argmin
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    fp_accuracy = kwargs.pop("fp_accuracy", blosc2.FPAccuracy.DEFAULT)
    jit = kwargs.pop("jit", None)
    jit_backend = kwargs.pop("jit_backend", None)
    where: dict | None = kwargs.pop("_where_args", None)
    reduce_op = reduce_args.pop("op")
    reduce_op_str = reduce_args.pop("op_str", None)
    axis = reduce_args["axis"]
    keepdims = reduce_args.get("keepdims", False)
    include_initial = reduce_args.pop("include_initial", False)
    dtype = reduce_args.get("dtype", None)
    if dtype is None:
        dtype = kwargs.pop("dtype", None)
        dtype = infer_reduction_dtype(dtype, reduce_op)
    else:
        del kwargs["dtype"]

    # Compute the shape and chunks of the output array, including broadcasting
    shape = compute_broadcast_shape(operands.values())

    # Validate axis against operand dimensions before any computation.
    if axis is not None and not np.isscalar(axis):
        ndim = len(shape)
        for ax in axis:
            if ax < -ndim or ax >= ndim:
                raise np.exceptions.AxisError(ax, ndim)

    _slice = _slice.raw
    shape_slice = shape
    mask_slice = np.array([isinstance(i, int) for i in _slice], dtype=np.bool_)
    if out is None and _slice != ():
        _slice = tuple(slice(i, i + 1, 1) if isinstance(i, int) else i for i in _slice)
        shape_slice = ndindex.ndindex(_slice).newshape(shape)
        # shape_slice in general not equal to final shape:
        # dummy dims (due to ints) will be dealt with by taking final_slice

    # after slicing, we reduce to calculate shape of output
    if axis is None:
        axis = tuple(range(len(shape_slice)))
    elif np.isscalar(axis):
        axis = (axis,)
    axis = tuple(a if a >= 0 else a + len(shape_slice) for a in axis)
    if np.any(mask_slice):
        add_idx = np.cumsum(mask_slice)
        axis = tuple(a + add_idx[a] for a in axis)  # axis now refers to new shape with dummy dims
    if reduce_args["axis"] is not None:
        # conserve as integer if was not tuple originally
        reduce_args["axis"] = axis[0] if np.isscalar(reduce_args["axis"]) else axis
    if reduce_op in {ReduceOp.CUMULATIVE_SUM, ReduceOp.CUMULATIVE_PROD}:
        reduced_shape = (np.prod(shape_slice),) if reduce_args["axis"] is None else shape_slice
        # if reduce_args["axis"] is None, have to have 1D input array; otherwise, ensure positive scalar
        reduce_args["axis"] = 0 if reduce_args["axis"] is None else axis[0]
        if include_initial:
            reduced_shape = tuple(
                s + 1 if i == reduce_args["axis"] else s for i, s in enumerate(shape_slice)
            )
    else:
        if keepdims:
            reduced_shape = tuple(1 if i in axis else s for i, s in enumerate(shape_slice))
        else:
            reduced_shape = tuple(s for i, s in enumerate(shape_slice) if i not in axis)
            mask_slice = mask_slice[[i for i in range(len(mask_slice)) if i not in axis]]

    if out is not None and reduced_shape != out.shape:
        raise ValueError("Provided output shape does not match the reduced shape.")

    # Choose the array with the largest shape as the reference for chunks
    # Note: we could have expr = blosc2.lazyexpr('numpy_array + 1') (i.e. no choice for chunks)
    blosc2_arrs = tuple(o for o in operands.values() if hasattr(o, "chunks"))
    fast_path = False
    all_ndarray = False
    any_persisted = False
    chunks = None
    blocks = None
    if blosc2_arrs:  # fast path only relevant if there are blosc2 arrays
        operand = max(blosc2_arrs, key=lambda x: len(x.shape))

        # Check if the partitions are aligned (i.e. all operands have the same shape,
        # chunks and blocks, and have no padding). This will allow us to take the fast path.
        same_shape = all(operand.shape == o.shape for o in operands.values() if hasattr(o, "shape"))
        same_chunks = all(operand.chunks == o.chunks for o in operands.values() if hasattr(o, "chunks"))
        same_blocks = all(operand.blocks == o.blocks for o in operands.values() if hasattr(o, "blocks"))
        fast_path = same_shape and same_chunks and same_blocks and (0 not in operand.chunks)
        aligned = dict.fromkeys(operands.keys(), False)
        iter_disk = False
        if fast_path:
            chunks = operand.chunks
            blocks = operand.blocks
            # Check that all operands are NDArray for fast path
            all_ndarray = all(
                isinstance(value, blosc2.NDArray) and value.shape != () for value in operands.values()
            )
            # Check that there is some NDArray that is persisted in the disk
            any_persisted = any(
                (
                    isinstance(value, blosc2.NDArray)
                    and value.shape != ()
                    and value.schunk.urlpath is not None
                )
                for value in operands.values()
            )
            if not blosc2.IS_WASM:
                iter_disk = all_ndarray and any_persisted
                # Experiments say that iter_disk is faster than the regular path for reductions
                # even when all operands are in memory, so no need to check any_persisted
                # New benchmarks are saying the contrary (> 10% slower), so this needs more
                # investigation
                # iter_disk = all_ndarray
            else:
                # WebAssembly does not support threading, so we cannot use the iter_disk option
                iter_disk = False
        else:
            for arr in blosc2_arrs:
                if arr.shape == shape:
                    chunks = arr.chunks
                    break
    if chunks is None:  # have to calculate chunks (this is cheap as empty just creates a thin metalayer)
        temp = blosc2.empty(shape, dtype=dtype)
        chunks = temp.chunks
        del temp

    # miniexpr reduction path only supported for some cases so far
    if not (fast_path and all_ndarray and reduced_shape == () and _slice == ()):
        use_miniexpr = False

    # Some reductions are not supported yet in miniexpr
    if reduce_op in (ReduceOp.ARGMAX, ReduceOp.ARGMIN, ReduceOp.CUMULATIVE_PROD, ReduceOp.CUMULATIVE_SUM):
        use_miniexpr = False

    # Check whether we can use miniexpr
    if use_miniexpr and isinstance(expression, str):
        has_complex = any(
            isinstance(op, blosc2.NDArray) and blosc2.isdtype(op.dtype, "complex floating")
            for op in operands.values()
        )
        if has_complex and (sys.platform == "win32" or blosc2.IS_WASM):
            # On Windows and WebAssembly, miniexpr has issues with complex numbers
            use_miniexpr = False
        if has_complex and any(tok in expression for tok in ("!=", "==", "<=", ">=", "<", ">")):
            use_miniexpr = False
        if where is not None and len(where) != 2:
            use_miniexpr = False

    if use_miniexpr:
        # Experiments say that not splitting is best (at least on Apple Silicon M4 Pro)
        cparams = kwargs.pop("cparams", blosc2.CParams(splitmode=blosc2.SplitMode.NEVER_SPLIT))
        # Create a fake NDArray just to drive the miniexpr evaluation (values won't be used)
        res_eval = blosc2.uninit(shape, dtype, chunks=chunks, blocks=blocks, cparams=cparams, **kwargs)
        # Compute the number of blocks in the result
        nblocks = res_eval.nbytes // res_eval.blocksize
        # Initialize aux_reduc based on the reduction operation
        # Padding blocks won't be written, so initial values matter for the final reduction
        if reduce_op in {ReduceOp.SUM, ReduceOp.ANY, ReduceOp.CUMULATIVE_SUM}:
            aux_reduc = np.zeros(nblocks, dtype=dtype)
        elif reduce_op in {ReduceOp.PROD, ReduceOp.ALL, ReduceOp.CUMULATIVE_PROD}:
            aux_reduc = np.ones(nblocks, dtype=dtype)
        elif reduce_op == ReduceOp.MIN:
            if np.issubdtype(dtype, np.integer):
                aux_reduc = np.full(nblocks, np.iinfo(dtype).max, dtype=dtype)
            else:
                aux_reduc = np.full(nblocks, np.inf, dtype=dtype)
        elif reduce_op == ReduceOp.MAX:
            if np.issubdtype(dtype, np.integer):
                aux_reduc = np.full(nblocks, np.iinfo(dtype).min, dtype=dtype)
            else:
                aux_reduc = np.full(nblocks, -np.inf, dtype=dtype)
        else:
            # For other operations, zeros should be safe
            aux_reduc = np.zeros(nblocks, dtype=dtype)
        prefilter_set = False
        expression_miniexpr = None
        try:
            if where is not None:
                expression_miniexpr = f"{reduce_op_str}(where({expression}, _where_x, _where_y))"
            else:
                expression_miniexpr = f"{reduce_op_str}({expression})"
            expression_miniexpr = _apply_jit_backend_pragma(expression_miniexpr, operands, jit_backend)
            res_eval._set_pref_expr(expression_miniexpr, operands, fp_accuracy, aux_reduc, jit=jit)
            prefilter_set = True
            # print("expr->miniexpr:", expression, reduce_op, fp_accuracy)
            # Data won't even try to be compressed, so buffers can be unitialized and reused
            data = np.empty(res_eval.schunk.chunksize, dtype=np.uint8)
            chunk_data = np.empty(res_eval.schunk.chunksize + blosc2.MAX_OVERHEAD, dtype=np.uint8)
            # Exercise prefilter for each chunk
            for nchunk in range(res_eval.schunk.nchunks):
                res_eval.schunk._prefilter_data(nchunk, data, chunk_data)
        except Exception as e:
            use_miniexpr = False
            if callable(expression) and _is_dsl_kernel_expression(expression):
                reason = "miniexpr compilation or execution failed for this DSL kernel."
                backend_error = str(e)
                parse_hint = None
                if isinstance(expression_miniexpr, str):
                    parse_hint = _format_dsl_parse_error_hint(expression_miniexpr, backend_error)
                reason = f"{reason}\nBackend error: {backend_error}"
                if parse_hint is not None:
                    reason = f"{reason}\n{parse_hint}"
                raise RuntimeError(_dsl_miniexpr_required_message(reason)) from e
        finally:
            if prefilter_set:
                res_eval.schunk.remove_prefilter("miniexpr")
            global iter_chunks
            # Ensure any background reading thread is closed
            iter_chunks = None

        if not use_miniexpr:
            # If miniexpr failed, fallback to regular evaluation
            # (continue to the manual chunked evaluation below)
            pass
        else:
            if reduce_op in {ReduceOp.ANY, ReduceOp.ALL}:
                result = reduce_op.value(aux_reduc, **reduce_args)
            else:
                # The accumulator is always 1-D (one slot per output block).
                # The original axis may refer to dimensions that no longer
                # exist after per-block reduction.  Use axis=0 to combine
                # all block results.
                result = reduce_op.value.reduce(aux_reduc, axis=0)
            return result

    # Iterate over the operands and get the chunks
    chunk_operands = {}
    # Check which chunks intersect with _slice
    if np.isscalar(reduce_args["axis"]):  # iterate over chunks incrementing along reduction axis
        intersecting_chunks = get_intersecting_chunks(_slice, shape, chunks, axis=reduce_args["axis"])
    else:  # iterate over chunks incrementing along last axis
        intersecting_chunks = get_intersecting_chunks(_slice, shape, chunks)
    out_init = False
    res_out_init = False
    ratio = (
        np.ceil(np.asarray(shape) / np.asarray(chunks)).astype(np.int64)
        if 0 not in chunks
        else np.asarray(shape)
    )

    for chunk_slice in intersecting_chunks:
        cslice = chunk_slice.raw
        nchunk = (
            builtins.sum(c.start // chunks[i] * np.prod(ratio[i + 1 :]) for i, c in enumerate(cslice))
            if 0 not in chunks
            else 0
        )
        # Check whether current cslice intersects with _slice
        if cslice != () and _slice != ():
            # get intersection of chunk and target
            cslice = step_handler(cslice, _slice)
        offset = tuple(s.start for s in cslice)  # offset for the udf
        starts = [s.start if s.start is not None else 0 for s in cslice]
        unit_steps = np.all([s.step == 1 for s in cslice])
        cslice_shape = tuple(s.stop - s.start for s in cslice)
        # get local index of part of out that is to be updated
        cslice_subidx = ndindex.ndindex(cslice).as_subindex(_slice).raw  # if _slice is (), just gives cslice
        if _slice == () and fast_path and unit_steps:
            # Fast path
            full_chunk = cslice_shape == chunks
            fill_chunk_operands(
                operands,
                cslice,
                cslice_shape,
                full_chunk,
                aligned,
                nchunk,
                iter_disk,
                chunk_operands,
                reduc=True,
                axis=reduce_args["axis"] if np.isscalar(reduce_args["axis"]) else None,
            )
        else:
            get_chunk_operands(operands, cslice, chunk_operands, shape)

        if reduce_op in {ReduceOp.CUMULATIVE_PROD, ReduceOp.CUMULATIVE_SUM}:
            reduced_slice = (
                tuple(
                    slice(sl.start + 1, sl.stop + 1, sl.step) if i == reduce_args["axis"] else sl
                    for i, sl in enumerate(cslice_subidx)
                )
                if include_initial
                else cslice_subidx
            )
        else:
            reduced_slice = (
                tuple(slice(None) if i in axis else sl for i, sl in enumerate(cslice_subidx))
                if keepdims
                else tuple(sl for i, sl in enumerate(cslice_subidx) if i not in axis)
            )

        # Evaluate and reduce the expression using chunks of operands

        if callable(expression):
            if _is_dsl_kernel_expression(expression):
                _raise_dsl_miniexpr_required(
                    "internal reduction fallback attempted to execute the DSL kernel directly in Python."
                )
            # TODO: Implement the reductions for UDFs (and test them)
            result = np.empty(cslice_shape, dtype=out.dtype)
            expression(tuple(chunk_operands.values()), result, offset=offset)
            # Reduce the result
            result = reduce_op.value.reduce(result, **reduce_args)
            # Update the output array with the result
            out[reduced_slice] = reduce_op.value(out[reduced_slice], result)
            continue

        result, _ = _get_result(expression, chunk_operands, ne_args, where)
        # Enforce contiguity of result (necessary to fill the out array)
        # but avoid copy if already contiguous
        result = np.require(result, requirements="C")

        # Reduce the result
        if result.shape == ():
            if reduce_op == ReduceOp.SUM and result[()] == 0:
                # Avoid a reduction when result is a zero scalar. Faster for sparse data.
                continue
            # Note that cslice_shape refers to slice of operand chunks, not reduced_slice
            result = np.full(cslice_shape, result[()])
        if reduce_op in {ReduceOp.ANY, ReduceOp.ALL, ReduceOp.CUMULATIVE_SUM, ReduceOp.CUMULATIVE_PROD}:
            result = reduce_op.value(result, **reduce_args)
        elif reduce_op in {ReduceOp.ARGMAX, ReduceOp.ARGMIN}:
            # offset for start of slice
            slice_ref = (
                starts
                if _slice == ()
                else [
                    (s - sl.start - np.sign(sl.step)) // sl.step + 1
                    for s, sl in zip(starts, _slice, strict=True)
                ]
            )
            result_idx = reduce_op.value(result, **reduce_args)
            if reduce_args["axis"] is None:  # indexing into flattened array
                result = result[np.unravel_index(result_idx, shape=result.shape)]
                idx_within_cslice = np.unravel_index(result_idx, shape=cslice_shape)
                result_idx = np.ravel_multi_index(
                    tuple(o + i for o, i in zip(slice_ref, idx_within_cslice, strict=True)), shape_slice
                )
            else:  # axis is an integer
                result = np.take_along_axis(
                    result,
                    np.expand_dims(result_idx, axis=reduce_args["axis"]) if not keepdims else result_idx,
                    axis=reduce_args["axis"],
                )
                result = result if keepdims else result.squeeze(axis=reduce_args["axis"])
                result_idx += slice_ref[reduce_args["axis"]]
        else:
            result = reduce_op.value.reduce(result, **reduce_args)

        if not out_init:
            # if cumsum/cumprod and arrays large, return blosc2 array with same chunks
            chunks_out = (
                chunks
                if np.prod(reduced_shape) * np.dtype(dtype).itemsize > 4 * blosc2.MAX_FAST_PATH_SIZE
                else None
            )
            chunks_out = chunks_out if _slice == () else None
            out_ = convert_none_out(result.dtype, reduce_op, reduced_shape, chunks=chunks_out)
            if out is not None:
                out[:] = out_
                del out_
            else:
                out = out_
            behaved = (
                False
                if not hasattr(out, "chunks")
                else blosc2.are_partitions_behaved(out.shape, out.chunks, out.blocks)
            )
            out_init = True

        # res_out only used be argmin/max and cumulative_sum/prod which only accept axis=int argument
        if (not res_out_init) or (
            np.isscalar(reduce_args["axis"]) and cslice_subidx[reduce_args["axis"]].start == 0
        ):  # starting reduction again along axis
            res_out_ = _get_res_out(result.shape, reduce_args["axis"], dtype, reduce_op)
            res_out_init = True

        # Update the output array with the result
        if reduce_op == ReduceOp.ANY:
            out[reduced_slice] += result
        elif reduce_op == ReduceOp.ALL:
            out[reduced_slice] *= result
        elif res_out_ is not None:
            # need lowest index for which optimum attained
            if reduce_op in {ReduceOp.ARGMAX, ReduceOp.ARGMIN}:
                cond = (res_out_ == result) & (result_idx < out[reduced_slice])
                cond |= res_out_ < result if reduce_op == ReduceOp.ARGMAX else res_out_ > result
                out[reduced_slice] = np.where(cond, result_idx, out[reduced_slice])
                res_out_ = np.where(cond, result, res_out_)
            else:  # CUMULATIVE_SUM or CUMULATIVE_PROD
                idx_lastval = tuple(
                    slice(-1, None) if i == reduce_args["axis"] else slice(None, None)
                    for i, c in enumerate(reduced_slice)
                )
                if reduce_op == ReduceOp.CUMULATIVE_SUM:
                    result += res_out_
                else:  # CUMULATIVE_PROD
                    result *= res_out_
                res_out_ = result[idx_lastval]
                if behaved and result.shape == out.chunks and result.dtype == out.dtype and _slice == ():
                    # Fast path
                    # TODO: Check this only works when slice is () as nchunk is incorrect  for out otherwise
                    out.schunk.update_data(nchunk, result, copy=False)
                else:
                    out[reduced_slice] = result
        else:
            out[reduced_slice] = reduce_op.value(out[reduced_slice], result)

    # No longer need res_out_
    del res_out_

    if out is None:
        if reduce_op in (ReduceOp.MIN, ReduceOp.MAX, ReduceOp.ARGMIN, ReduceOp.ARGMAX):
            raise ValueError("zero-size array in (arg-)min/max reduction operation is not supported")
        if dtype is None:
            # We have no hint here, so choose a default dtype
            dtype = np.float64
        out = convert_none_out(dtype, reduce_op, reduced_shape)

    if reduced_shape == ():
        # convert_none_out() may allocate shape (1,) as an internal buffer for scalar reductions.
        # Collapse it to a numpy scalar while handling both 0-d and 1-d singleton arrays.
        if isinstance(out, np.ndarray):
            out = out[()] if out.ndim == 0 else out[0]
        else:
            out = out[()]
    final_mask = tuple(np.where(mask_slice)[0])
    if np.any(mask_slice):  # remove dummy dims
        out = np.squeeze(out, axis=final_mask)
    # Check if the output array needs to be converted into a blosc2.NDArray
    if kwargs != {} and not np.isscalar(out):
        out = blosc2.asarray(out, **kwargs)
    return out


def _get_res_out(reduced_shape, axis, dtype, reduce_op):
    reduced_shape = (1,) if reduced_shape == () else reduced_shape
    # Get res_out to hold running sums along axes for chunks when doing cumulative sums/prods with axis not None
    if reduce_op in {ReduceOp.CUMULATIVE_SUM, ReduceOp.CUMULATIVE_PROD}:
        temp_shape = tuple(1 if i == axis else s for i, s in enumerate(reduced_shape))
        res_out_ = (
            np.zeros(temp_shape, dtype=dtype)
            if reduce_op == ReduceOp.CUMULATIVE_SUM
            else np.ones(temp_shape, dtype=dtype)
        )
    elif reduce_op in {ReduceOp.ARGMIN, ReduceOp.ARGMAX}:
        temp_shape = reduced_shape
        res_out_ = np.ones(temp_shape, dtype=dtype)
        if np.issubdtype(dtype, np.integer):
            res_out_ *= np.iinfo(dtype).max if reduce_op == ReduceOp.ARGMIN else np.iinfo(dtype).min
        elif np.issubdtype(dtype, np.bool):
            res_out_ = res_out_ if reduce_op == ReduceOp.ARGMIN else np.zeros(temp_shape, dtype=dtype)
        else:
            res_out_ *= np.inf if reduce_op == ReduceOp.ARGMIN else -np.inf
    else:
        res_out_ = None
    return res_out_


def convert_none_out(dtype, reduce_op, reduced_shape, chunks=None):
    reduced_shape = (1,) if reduced_shape == () else reduced_shape
    # out will be a proper numpy.ndarray
    if reduce_op in {ReduceOp.SUM, ReduceOp.CUMULATIVE_SUM, ReduceOp.PROD, ReduceOp.CUMULATIVE_PROD}:
        if reduce_op in (ReduceOp.CUMULATIVE_SUM, ReduceOp.CUMULATIVE_PROD) and chunks is not None:
            out = (
                blosc2.zeros(reduced_shape, dtype=dtype, chunks=chunks)
                if reduce_op == ReduceOp.CUMULATIVE_SUM
                else blosc2.ones(reduced_shape, dtype=dtype, chunks=chunks)
            )
        else:
            out = (
                np.zeros(reduced_shape, dtype=dtype)
                if reduce_op in {ReduceOp.SUM, ReduceOp.CUMULATIVE_SUM}
                else np.ones(reduced_shape, dtype=dtype)
            )
    elif reduce_op == ReduceOp.MIN:
        if np.issubdtype(dtype, np.integer):
            out = np.iinfo(dtype).max * np.ones(reduced_shape, dtype=dtype)
        else:
            out = np.inf * np.ones(reduced_shape, dtype=dtype)
    elif reduce_op == ReduceOp.MAX:
        if np.issubdtype(dtype, np.integer):
            out = np.iinfo(dtype).min * np.ones(reduced_shape, dtype=dtype)
        else:
            out = -np.inf * np.ones(reduced_shape, dtype=dtype)
    elif reduce_op == ReduceOp.ANY:
        out = np.zeros(reduced_shape, dtype=np.bool_)
    elif reduce_op == ReduceOp.ALL:
        out = np.ones(reduced_shape, dtype=np.bool_)
    elif reduce_op in {ReduceOp.ARGMIN, ReduceOp.ARGMAX}:
        out = np.zeros(reduced_shape, dtype=blosc2.DEFAULT_INDEX)
    return out


def _validate_chunked_eval_inputs(operands: dict, out, shape, reduce_args: dict) -> bool:
    if operands:
        _, _, _, fast_path = validate_inputs(operands, out, reduce=reduce_args != {})
        return fast_path
    if shape is None and out is None:
        raise ValueError(
            "For UDFs with no inputs, provide `shape` (or an output array) to indicate result shape"
        )
    return False


def _eval_zero_input_dsl_if_needed(
    expression,
    operands: dict,
    where,
    getitem: bool,
    item,
    shape,
    jit,
    jit_backend,
    kwargs: dict,
):
    use_zero_input_dsl_fast_eval = (
        not operands
        and isinstance(expression, DSLKernel)
        and expression.dsl_source is not None
        and where is None
    )
    if not use_zero_input_dsl_fast_eval:
        return False, None

    full_res = fast_eval(
        expression,
        operands,
        getitem=False,
        shape=shape,
        jit=jit,
        jit_backend=jit_backend,
        **kwargs,
    )
    if getitem:
        return True, full_res[item.raw]
    return True, full_res


def chunked_eval(
    expression: str | Callable[[tuple, np.ndarray, tuple[int]], None], operands: dict, item=(), **kwargs
):
    """
    Evaluate the expression in chunks of operands.

    This chooses the best algorithm exploring different paths depending on the input operands.

    Parameters
    ----------
    expression: str or callable
        The expression or user-defined function (udf) to evaluate.
    operands: dict
        A dictionary containing the operands for the expression.
    item: int, sequence of ints, slice, sequence of slices or None, optional
        The slice(s) of the operands to be used in computation. Note that step parameter is not honored yet.
        Item is used to slice the operands PRIOR to computation.
    kwargs: Any, optional
        Additional keyword arguments supported by the :func:`empty` constructor.  In addition,
        the following keyword arguments are supported:
        _getitem: bool, optional
            Indicates whether the expression is being evaluated for a getitem operation.
            Default is False.
        _output: blosc2.Array, optional
            The output array to store the result.
        _ne_args: dict, optional
            Additional arguments to be passed to `numexpr.evaluate()` function.
        _where_args: dict, optional
            Additional arguments for conditional evaluation.
    """
    try:
        # standardise slice to be ndindex.Tuple
        item = () if item == slice(None, None, None) else item
        item = item if isinstance(item, tuple) else (item,)
        item = tuple(
            slice(s.start, s.stop, 1 if s.step is None else s.step) if isinstance(s, slice) else s
            for s in item
        )
        item = ndindex.ndindex(item)
        shape = kwargs.pop("shape", None)
        if item.raw != () and shape is not None:
            item = item.expand(shape)  # converts to standard tuple form

        getitem = kwargs.pop("_getitem", False)
        out = kwargs.get("_output")
        # Execution policy for miniexpr JIT paths only; never forward to array constructors.
        jit = kwargs.pop("jit", None)
        jit_backend = kwargs.pop("jit_backend", None)

        where: dict | None = kwargs.get("_where_args")
        if where:
            # Make the where arguments part of the operands
            operands = {**operands, **where}

        reduce_args = kwargs.pop("_reduce_args", {})
        # Resolve the JS backend: explicit jit_backend="js", or prefer-js-with-fallback under
        # WebAssembly when the user left jit_backend unset (see _maybe_js_backend).
        expression, jit, jit_backend = _maybe_js_backend(
            expression, jit, jit_backend, reduce_args, operands, kwargs, shape=shape
        )

        fast_path = _validate_chunked_eval_inputs(operands, out, shape, reduce_args)

        # Activate last read cache for NDField instances
        for op in operands:
            if isinstance(operands[op], blosc2.NDField):
                operands[op].ndarr.keep_last_read = True

        if reduce_args:
            # Eval and reduce the expression in a single step
            return reduce_slices(
                expression,
                operands,
                reduce_args=reduce_args,
                _slice=item,
                jit=jit,
                jit_backend=jit_backend,
                **kwargs,
            )

        handled, result = _eval_zero_input_dsl_if_needed(
            expression, operands, where, getitem, item, shape, jit, jit_backend, kwargs
        )
        if handled:
            return result

        full_slice = is_full_slice(item.raw)
        if not full_slice or (where is not None and len(where) < 2):
            # The fast path is possible under a few conditions
            if getitem and (where is None or len(where) == 2):
                # Compute the size of operands for the fast path
                unit_steps = np.all([s.step == 1 for s in item.raw if isinstance(s, slice)])
                # shape of slice, if non-unit steps have to decompress full array into memory
                shape_operands = item.newshape(shape) if unit_steps else shape
                _dtype = np.dtype(kwargs.get("dtype", np.float64))
                size_operands = math.prod(shape_operands) * len(operands) * _dtype.itemsize
                # Only take the fast path if the size of operands is relatively small
                if size_operands < blosc2.MAX_FAST_PATH_SIZE:
                    return slices_eval_getitem(expression, operands, _slice=item, shape=shape, **kwargs)
            return slices_eval(expression, operands, getitem=getitem, _slice=item, shape=shape, **kwargs)

        fast_path = full_slice and fast_path
        if fast_path:  # necessarily item is ()
            if getitem:
                # When using getitem, taking the fast path is always possible
                return fast_eval(
                    expression,
                    operands,
                    getitem=True,
                    jit=jit,
                    jit_backend=jit_backend,
                    shape=shape,
                    **kwargs,
                )
            elif (kwargs.get("chunks") is None and kwargs.get("blocks") is None) and (
                out is None or isinstance(out, blosc2.NDArray)
            ):
                # If not, the conditions to use the fast path are a bit more restrictive
                # e.g. the user cannot specify chunks or blocks, or an output that is not
                # a blosc2.NDArray
                return fast_eval(
                    expression,
                    operands,
                    getitem=False,
                    jit=jit,
                    jit_backend=jit_backend,
                    shape=shape,
                    **kwargs,
                )
            elif _is_dsl_kernel_expression(expression) and (out is None or isinstance(out, blosc2.NDArray)):
                # DSL kernels require miniexpr and must not fall back to Python execution.
                return fast_eval(
                    expression,
                    operands,
                    getitem=False,
                    jit=jit,
                    jit_backend=jit_backend,
                    shape=shape,
                    **kwargs,
                )

        # End up here by default
        return slices_eval(expression, operands, getitem=getitem, _slice=item, shape=shape, **kwargs)

    finally:
        global iter_chunks
        # Ensure any background reading thread is closed
        iter_chunks = None


def fuse_operands(operands1, operands2):
    new_operands = {}
    dup_operands = {}
    new_pos = len(operands1)
    operand_to_key = {id(v): k for k, v in operands1.items()}
    for k2, v2 in operands2.items():
        try:
            k1 = operand_to_key[id(v2)]
            # The operand is duplicated; keep track of it
            dup_operands[k2] = k1
        except KeyError:
            # The value is not among operands1, so rebase it
            new_op = f"o{new_pos}"
            new_pos += 1
            new_operands[new_op] = v2
    return new_operands, dup_operands


def fuse_expressions(expr, new_base, dup_op):
    new_expr = ""
    skip_to_char = 0
    old_base = 0
    prev_pos = {}
    for i, expr_i in enumerate(expr):
        if i < skip_to_char:
            continue
        if expr_i == "o":
            if i > 0 and expr[i - 1] not in {" ", "("}:
                # Not a variable
                new_expr += expr_i
                continue
            # This is a variable.  Find the end of it.
            j = i + 1
            for k in range(len(expr[j:])):
                if expr[j + k] in " )[,":  # Added comma to the list of delimiters
                    j = k
                    break
            if expr[i + j] == ")":
                j -= 1
            # Extract only the numeric part, handling cases where there might be a comma
            operand_str = expr[i + 1 : i + j + 1]
            # Split by comma and take the first part (the operand index)
            operand_num_str = operand_str.split(",")[0]
            old_pos = int(operand_num_str)
            old_op = f"o{old_pos}"
            if old_op not in dup_op:
                if old_pos in prev_pos:
                    # Keep track of duplicated old positions inside expr
                    new_pos = prev_pos[old_pos]
                else:
                    new_pos = old_base + new_base
                    old_base += 1
                new_expr += f"o{new_pos}"
                prev_pos[old_pos] = new_pos
            else:
                new_expr += dup_op[old_op]
            skip_to_char = i + j + 1
        else:
            new_expr += expr_i
    return new_expr


def check_dtype(op, value1, value2):
    if op in ("contains", "startswith", "endswith"):
        return np.dtype(np.bool_)

    v1_dtype = blosc2.result_type(value1)
    v2_dtype = v1_dtype if value2 is None else blosc2.result_type(value2)
    if op in not_complex_ops and (v1_dtype == np.complex128 or v2_dtype == np.complex128):
        # Ensure that throw exception for functions which don't support complex args
        raise ValueError(f"Invalid operand type for {op}: {v1_dtype, v2_dtype}")
    if op in relational_ops:
        return np.dtype(np.bool_)
    if op in logical_ops:
        # Ensure that both operands are booleans or ints
        if v1_dtype not in (np.bool_, np.int32, np.int64):
            raise ValueError(f"Invalid operand type for {op}: {v1_dtype}")
        if v2_dtype not in (np.bool_, np.int32, np.int64):
            raise ValueError(f"Invalid operand type for {op}: {v2_dtype}")

    if op == "/":
        if v1_dtype == np.int32 and v2_dtype == np.int32:
            return blosc2.float32
        if np.issubdtype(v1_dtype, np.integer) and np.issubdtype(v2_dtype, np.integer):
            return blosc2.float64

    # Follow NumPy rules for scalar-array operations
    return blosc2.result_type(value1, value2)


def result_type(
    *arrays_and_dtypes: blosc2.NDArray | int | float | complex | bool | str | blosc2.dtype,
) -> blosc2.dtype:
    """
    Returns the dtype that results from applying type promotion rules (see Type Promotion Rules) to the arguments.

    Parameters
    ----------
    arrays_and_dtypes: Sequence[NDarray | int | float | complex | bool | blosc2.dtype])
        An arbitrary number of input arrays, scalars, and/or dtypes.

    Returns
    -------
    out: blosc2.dtype
        The dtype resulting from an operation involving the input arrays, scalars, and/or dtypes.
    """
    # Follow NumPy rules for scalar-array operations
    # Create small arrays with the same dtypes and let NumPy's type promotion determine the result type
    arrs = [
        (np.array(value).dtype if isinstance(value, str | bytes) else value)
        if (np.isscalar(value) or not hasattr(value, "dtype"))
        else np.array([0], dtype=convert_dtype(value.dtype))
        for value in arrays_and_dtypes
    ]
    return np.result_type(*arrs)


def can_cast(from_: blosc2.dtype | blosc2.NDArray, to: blosc2.dtype) -> bool:
    """
    Determines if one data type can be cast to another data type according to (NumPy) type promotion rules.

    Parameters
    ----------
    from_: dtype | NDArray
        Input data type or array from which to cast.

    to: dtype
    Desired data type.

    Returns
    -------
    out:bool
        True if the cast can occur according to type promotion rules; otherwise, False.
    """
    arrs = np.array([0], dtype=from_.dtype) if hasattr(from_, "shape") else from_
    return np.result_type(arrs)


class LazyExpr(LazyArray):
    """Class for hosting lazy expressions.

    This is not meant to be called directly from user space.

    Once the lazy expression is created, it can be evaluated via :func:`LazyExpr.compute`.
    """

    def __init__(self, new_op):  # noqa: C901
        if new_op is None:
            self.expression = ""
            self.operands = {}
            return
        value1, op, value2 = new_op
        value1 = value1._raw_col if hasattr(value1, "_raw_col") else value1
        value2 = value2._raw_col if hasattr(value2, "_raw_col") else value2
        dtype_ = check_dtype(op, value1, value2)  # perform some checks
        # Check that operands are proper Operands, LazyArray or scalars; if not, convert to NDArray objects
        value1 = (
            blosc2.SimpleProxy(value1)
            if not (isinstance(value1, blosc2.Operand | np.ndarray) or np.isscalar(value1))
            else value1
        )
        # Reset values represented as np.int64 etc. to be set as Python natives
        value1 = value1.item() if np.isscalar(value1) and hasattr(value1, "item") else value1
        if value2 is None:
            if isinstance(value1, LazyExpr):
                self.expression = value1.expression if op is None else f"{op}({value1.expression})"
                # handle constructors which can give empty operands
                self._dtype = (
                    value1.dtype
                    if op is None
                    else _numpy_eval_expr(f"{op}(o0)", {"o0": value1}, prefer_blosc=False).dtype
                )
                self.operands = value1.operands
            else:
                if np.isscalar(value1):
                    value1 = ne_evaluate(f"{op}({value1!r})")
                    op = None
                self.operands = {"o0": value1}
                self.expression = "o0" if op is None else f"{op}(o0)"
            return
        value2 = (
            blosc2.SimpleProxy(value2)
            if not (isinstance(value2, blosc2.Operand | np.ndarray) or np.isscalar(value2))
            else value2
        )

        # Reset values represented as np.int64 etc. to be set as Python natives,
        # BUT preserve numpy integer scalars that require explicit typing (unsigned or
        # 64-bit) so that dtype-sensitive backends (numexpr) don't downcast them to int32.
        def _to_native_if_safe(v):
            if not (np.isscalar(v) and hasattr(v, "item")):
                return v
            dt = np.dtype(type(v))
            # Keep typed when unsigned or itemsize >= 8 to avoid silent int32 truncation.
            if np.issubdtype(dt, np.unsignedinteger) or dt.itemsize >= 8:
                return v
            return v.item()

        value2 = _to_native_if_safe(value2)

        # Non-finite Python floats repr as bare names ("nan", "inf") that no
        # expression evaluator defines; retype them so they take the
        # typed-scalar branches below and ride as named operands instead of
        # expression text.
        def _retype_nonfinite(v):
            if isinstance(v, float) and not math.isfinite(v):
                return np.float64(v)
            return v

        value1 = _retype_nonfinite(value1)
        value2 = _retype_nonfinite(value2)

        if isinstance(value1, LazyExpr) or isinstance(value2, LazyExpr):
            if isinstance(value1, LazyExpr):
                newexpr = value1.update_expr(new_op)
            else:
                newexpr = value2.update_expr(new_op)
            self.expression = newexpr.expression
            self.operands = newexpr.operands
            self._dtype = newexpr.dtype
            return
        elif op in funcs_2args:
            if np.isscalar(value1) and np.isscalar(value2):
                self.expression = "o0"
                svalue1 = format_expr_scalar(value1)
                svalue2 = format_expr_scalar(value2)
                self.operands = {"o0": ne_evaluate(f"{op}({svalue1}, {svalue2})")}  # eager evaluation
            elif np.isscalar(value2):
                if isinstance(value2, np.floating) and not math.isfinite(value2):
                    # nan/inf have no expression literal — keep as named operand
                    self.operands = {"o0": value1, "o1": value2}
                    self.expression = f"{op}(o0, o1)"
                else:
                    self.operands = {"o0": value1}
                    self.expression = f"{op}(o0, {format_expr_scalar(value2)})"
            elif np.isscalar(value1):
                if isinstance(value1, np.floating) and not math.isfinite(value1):
                    self.operands = {"o0": value2, "o1": value1}
                    self.expression = f"{op}(o1, o0)"
                else:
                    self.operands = {"o0": value2}
                    self.expression = f"{op}({format_expr_scalar(value1)}, o0)"
            else:
                self.operands = {"o0": value1, "o1": value2}
                self.expression = f"{op}(o0, o1)"
            return

        self._dtype = dtype_
        if np.isscalar(value1) and np.isscalar(value2):
            self.expression = "o0"
            self.operands = {"o0": ne_evaluate(f"({value1!r} {op} {value2!r})")}  # eager evaluation
        elif np.isscalar(value2):
            if hasattr(value2, "dtype"):  # typed numpy scalar — keep as named operand
                self.operands = {"o0": value1, "o1": value2}
                self.expression = f"(o0 {op} o1)"
            else:
                self.operands = {"o0": value1}
                self.expression = f"(o0 {op} {value2!r})"
        elif hasattr(value2, "shape") and value2.shape == ():
            self.operands = {"o0": value1}
            self.expression = f"(o0 {op} {value2[()]})"
        elif np.isscalar(value1):
            if hasattr(value1, "dtype"):  # typed numpy scalar — keep as named operand
                self.operands = {"o0": value2, "o1": value1}
                self.expression = f"(o1 {op} o0)"
            else:
                self.operands = {"o0": value2}
                self.expression = f"({value1!r} {op} o0)"
        elif hasattr(value1, "shape") and value1.shape == ():
            self.operands = {"o0": value2}
            self.expression = f"({value1[()]} {op} o0)"
        else:
            if value1 is value2:
                self.operands = {"o0": value1}
                self.expression = f"(o0 {op} o0)"
            else:
                # This is the very first time that a LazyExpr is formed from two operands
                # that are not LazyExpr themselves
                self.operands = {"o0": value1, "o1": value2}
                self.expression = f"(o0 {op} o1)"

    def update_expr(self, new_op):
        prev_flag = blosc2._disable_overloaded_equal
        # We use a lot of the original NDArray.__eq__ as 'is', so deactivate the overloaded one
        blosc2._disable_overloaded_equal = True
        # One of the two operands are LazyExpr instances
        try:
            value1, op, value2 = new_op
            value1 = value1._raw_col if hasattr(value1, "_raw_col") else value1
            value2 = value2._raw_col if hasattr(value2, "_raw_col") else value2
            dtype_ = check_dtype(op, value1, value2)  # conserve dtype
            # The new expression and operands
            expression = None
            new_operands = {}
            # where() handling requires evaluating the expression prior to merge.
            # This is different from reductions, where the expression is evaluated
            # and returned a NumPy array (for usability convenience).
            # We do things like this to enable the fusion of operations like
            # `a.where(0, 1).sum()`.
            # Another possibility would have been to always evaluate where() and produce
            # an NDArray, but that would have been less efficient for the case above.
            if hasattr(value1, "_where_args"):
                value1 = value1.compute()
            if hasattr(value2, "_where_args"):
                value2 = value2.compute()

            if not isinstance(value1, LazyExpr) and not isinstance(value2, LazyExpr):
                # We converted some of the operands to NDArray (where() handling above)
                new_operands = {"o0": value1, "o1": value2}
                expression = "op(o0, o1)" if op in funcs_2args else f"(o0 {op} o1)"
                return self._new_expr(expression, new_operands, guess=False, out=None, where=None)
            elif isinstance(value1, LazyExpr) and isinstance(value2, LazyExpr):
                # Expression fusion
                # Fuse operands in expressions and detect duplicates
                new_operands, dup_op = fuse_operands(value1.operands, value2.operands)
                # Take expression 2 and rebase the operands while removing duplicates
                new_expr = fuse_expressions(value2.expression, len(value1.operands), dup_op)
                expression = (
                    f"{op}({value1.expression}, {new_expr})"
                    if op in funcs_2args
                    else f"({value1.expression} {op} {new_expr})"
                )
                def_operands = value1.operands
            elif isinstance(value1, LazyExpr):
                if np.isscalar(value2):
                    v2 = format_expr_scalar(value2)
                elif hasattr(value2, "shape") and value2.shape == ():
                    v2 = format_expr_scalar(value2[()])
                else:
                    operand_to_key = {id(v): k for k, v in value1.operands.items()}
                    try:
                        v2 = operand_to_key[id(value2)]
                    except KeyError:
                        v2 = f"o{len(value1.operands)}"
                        new_operands = {v2: value2}
                if op == "~":
                    expression = f"({op}{value1.expression})"
                else:
                    expression = (
                        f"{op}({value1.expression}, {v2})"
                        if op in funcs_2args
                        else f"({value1.expression} {op} {v2})"
                    )
                def_operands = value1.operands
            else:
                if np.isscalar(value1):
                    v1 = format_expr_scalar(value1)
                elif hasattr(value1, "shape") and value1.shape == ():
                    v1 = format_expr_scalar(value1[()])
                else:
                    operand_to_key = {id(v): k for k, v in value2.operands.items()}
                    try:
                        v1 = operand_to_key[id(value1)]
                    except KeyError:
                        v1 = f"o{len(value2.operands)}"
                        new_operands = {v1: value1}
                if op == "[]":  # syntactic sugar for slicing
                    expression = f"({v1}[{value2.expression}])"
                else:
                    expression = (
                        f"{op}({v1}, {value2.expression})"
                        if op in funcs_2args
                        else f"({v1} {op} {value2.expression})"
                    )
                    def_operands = value2.operands
            # Return a new expression
            operands = def_operands | new_operands
            expr = self._new_expr(expression, operands, guess=False, out=None, where=None)
            expr._dtype = dtype_  # override dtype with preserved dtype
            return expr
        finally:
            blosc2._disable_overloaded_equal = prev_flag

    @property
    def dtype(self):
        # Honor self._dtype; it can be set during the building of the expression
        if hasattr(self, "_dtype"):
            # In some situations, we already know the dtype
            return self._dtype
        if (
            hasattr(self, "_dtype_")
            and hasattr(self, "_expression_")
            and self._expression_ == self.expression
        ):
            # Use the cached dtype
            return self._dtype_

        # Return None if there is a missing operand (e.g. a removed file on disk)
        if any(v is None for v in self.operands.values()):
            return None

        _out = _numpy_eval_expr(self.expression, self.operands, prefer_blosc=False)
        self._dtype_ = _out.dtype
        self._expression_ = self.expression
        return self._dtype_

    @property
    def ndim(self) -> int:
        return len(self.shape)

    @property
    def shape(self):
        # Honor self._shape; it can be set during the building of the expression
        if hasattr(self, "_shape"):
            return self._shape
        if (
            hasattr(self, "_shape_")
            and hasattr(self, "_expression_")
            and self._expression_ == self.expression
        ):
            # Use the cached shape
            return self._shape_

        # Return None if there is a missing operand (e.g. a removed file on disk)
        if any(v is None for v in self.operands.values()):
            return None

        # Operands shape can change, so we always need to recompute this
        if any(_has_constructor_call(self.expression, constructor) for constructor in constructors):
            # might have an expression with pure constructors
            opshapes = {k: v if not hasattr(v, "shape") else v.shape for k, v in self.operands.items()}
            _shape = infer_shape(self.expression, opshapes)  # infer shape, includes constructors
        else:
            _shape, chunks, blocks, fast_path = validate_inputs(self.operands, getattr(self, "_out", None))
            if fast_path:
                # fast_path ensure that all the operands have the same partitions
                self._chunks = chunks
                self._blocks = blocks

        self._shape_ = _shape
        self._expression_ = self.expression
        return _shape

    @property
    def chunks(self):
        if hasattr(self, "_chunks"):
            return self._chunks
        shape, self._chunks, self._blocks, fast_path = validate_inputs(
            self.operands, getattr(self, "_out", None)
        )
        if not hasattr(self, "_shape"):
            self._shape = shape
        if self._shape != shape:  # validate inputs only works for elementwise funcs so returned shape might
            fast_path = False  # be incompatible with true output shape
        if not fast_path:
            # Not using the fast path, so we need to compute the chunks/blocks automatically
            self._chunks, self._blocks = compute_chunks_blocks(self.shape, None, None, dtype=self.dtype)
        return self._chunks

    @property
    def blocks(self):
        if hasattr(self, "_blocks"):
            return self._blocks
        shape, self._chunks, self._blocks, fast_path = validate_inputs(
            self.operands, getattr(self, "_out", None)
        )
        if not hasattr(self, "_shape"):
            self._shape = shape
        if self._shape != shape:  # validate inputs only works for elementwise funcs so returned shape might
            fast_path = False  # be incompatible with true output shape
        if not fast_path:
            # Not using the fast path, so we need to compute the chunks/blocks automatically
            self._chunks, self._blocks = compute_chunks_blocks(self.shape, None, None, dtype=self.dtype)
        return self._blocks

    def where(self, value1=None, value2=None):
        """
        Select value1 or value2 values based on the condition of the current expression.

        Parameters
        ----------
        value1: array_like, optional
            The value to select when the condition is True.
        value2: array_like, optional
            The value to select when the condition is False.

        Returns
        -------
        out: LazyExpr
            A new expression with the where condition applied.
        """
        if not np.issubdtype(self.dtype, np.bool_):
            raise ValueError("where() can only be used with boolean expressions")
        # This just acts as a 'decorator' for the existing expression
        if value1 is not None and value2 is not None:
            # Guess the outcome dtype for value1 and value2
            dtype = blosc2.result_type(value1, value2)
            args = {"_where_x": value1, "_where_y": value2}
        elif value1 is not None:
            if hasattr(value1, "dtype"):
                dtype = value1.dtype
            else:
                dtype = np.asarray(value1).dtype
            args = {"_where_x": value1}
        elif value2 is not None:
            raise ValueError("where() requires value1 when using value2")
        else:
            args = {}
            dtype = None

        # Create a new expression
        new_expr = blosc2.LazyExpr(new_op=(self, None, None))
        new_expr.expression = self.expression
        new_expr.operands = self.operands
        new_expr._where_args = args
        new_expr._dtype = dtype
        return new_expr

    @staticmethod
    def _normalize_where(where):
        if where is None:
            return None
        raw_col = getattr(where, "_raw_col", None)
        if raw_col is not None:
            where = raw_col
        if isinstance(where, np.ndarray) and where.dtype == np.bool_:
            where = blosc2.asarray(where)
        if not (
            isinstance(where, (blosc2.NDArray, blosc2.LazyExpr))
            and getattr(where, "dtype", None) == np.bool_
        ):
            raise TypeError(f"Expected boolean blosc2.NDArray or LazyExpr, got {type(where).__name__}")
        return where

    def _where_selected(self, where):
        where = self._normalize_where(where)
        return self if where is None else where.where(self)

    def _where_identity_expr(self, where, identity):
        where = self._normalize_where(where)
        return self if where is None else where.where(self, identity)

    def _reduction_identity(self, op):
        dtype = np.dtype(self.dtype)
        if dtype.kind == "b":
            return op == "min"
        if dtype.kind in "iu":
            info = np.iinfo(dtype)
            return info.max if op == "min" else info.min
        if dtype.kind == "f":
            return np.inf if op == "min" else -np.inf
        raise TypeError(f"where= for {op} is not supported for dtype {dtype!r}")

    def sum(
        self,
        axis=None,
        dtype=None,
        keepdims=False,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        if where is not None:
            return self._where_identity_expr(where, 0).sum(
                axis=axis, dtype=dtype, keepdims=keepdims, fp_accuracy=fp_accuracy, **kwargs
            )
        reduce_args = {
            "op": ReduceOp.SUM,
            "op_str": "sum",
            "axis": axis,
            "dtype": dtype,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def prod(
        self,
        axis=None,
        dtype=None,
        keepdims=False,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        if where is not None:
            return self._where_identity_expr(where, 1).prod(
                axis=axis, dtype=dtype, keepdims=keepdims, fp_accuracy=fp_accuracy, **kwargs
            )
        reduce_args = {
            "op": ReduceOp.PROD,
            "op_str": "prod",
            "axis": axis,
            "dtype": dtype,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def get_num_elements(self, axis, item):
        if hasattr(self, "_where_args") and len(self._where_args) == 1:
            # We have a where condition, so we need to count the number of elements
            # fulfilling the condition
            orig_where_args = self._where_args
            self._where_args = {"_where_x": blosc2.ones(self.shape, dtype=np.int8)}
            num_elements = self.sum(axis=axis, dtype=np.int64, item=item)
            self._where_args = orig_where_args
            return num_elements
        # Compute the number of elements in the array
        shape = self.shape
        if np.isscalar(axis):
            axis = (axis,)
        if item != ():
            # Compute the shape of the slice
            shape = ndindex.ndindex(item).newshape(shape)
        axis = tuple(range(len(shape))) if axis is None else axis
        axis = tuple(a if a >= 0 else a + len(shape) for a in axis)  # handle negative indexing
        return math.prod([shape[i] for i in axis])

    def mean(
        self,
        axis=None,
        dtype=None,
        keepdims=False,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        where = self._normalize_where(where)
        expr = self if where is None else where.where(self, 0)
        item = kwargs.pop("item", ())
        total_sum = expr.sum(
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            item=item,
            fp_accuracy=fp_accuracy,
        )
        num_elements = (
            self.get_num_elements(axis, item)
            if where is None
            else where.where(blosc2.ones(self.shape, dtype=np.int64), 0).sum(axis=axis, dtype=np.int64)
        )
        if np.isscalar(num_elements) and num_elements == 0:
            raise ValueError("mean of an empty array is not defined")
        out = total_sum / num_elements
        out2 = kwargs.pop("out", None)
        if out2 is not None:
            out2[:] = out
            return out2
        if kwargs != {} and not np.isscalar(out):
            out = blosc2.asarray(out, **kwargs)
        return out

    def std(
        self,
        axis=None,
        dtype=None,
        keepdims=False,
        ddof=0,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        where = self._normalize_where(where)
        item = kwargs.pop("item", ())
        if item == ():  # fast path
            mean_value = self.mean(
                axis=axis, dtype=dtype, keepdims=True, where=where, fp_accuracy=fp_accuracy
            )
            expr = (self - mean_value) ** 2
        else:
            mean_value = self.mean(
                axis=axis, dtype=dtype, keepdims=True, where=where, item=item, fp_accuracy=fp_accuracy
            )
            # TODO: Not optimal because we load the whole slice in memory. Would have to write
            #  a bespoke std function that executed within slice_eval to avoid this probably.
            expr = (self.slice(item) - mean_value) ** 2
        out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, where=where, fp_accuracy=fp_accuracy)
        if ddof != 0:
            num_elements = (
                self.get_num_elements(axis, item)
                if where is None
                else where.where(blosc2.ones(self.shape, dtype=np.int64), 0).sum(axis=axis, dtype=np.int64)
            )
            out = np.sqrt(out * num_elements / (num_elements - ddof))
        else:
            out = np.sqrt(out)
        out2 = kwargs.pop("out", None)
        if out2 is not None:
            out2[:] = out
            return out2
        if kwargs != {} and not np.isscalar(out):
            out = blosc2.asarray(out, **kwargs)
        return out

    def var(
        self,
        axis=None,
        dtype=None,
        keepdims=False,
        ddof=0,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        where = self._normalize_where(where)
        item = kwargs.pop("item", ())
        if item == ():  # fast path
            mean_value = self.mean(
                axis=axis, dtype=dtype, keepdims=True, where=where, fp_accuracy=fp_accuracy
            )
            expr = (self - mean_value) ** 2
        else:
            mean_value = self.mean(
                axis=axis, dtype=dtype, keepdims=True, where=where, item=item, fp_accuracy=fp_accuracy
            )
            # TODO: Not optimal because we load the whole slice in memory. Would have to write
            #  a bespoke var function that executed within slice_eval to avoid this probably.
            expr = (self.slice(item) - mean_value) ** 2
        out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, where=where, fp_accuracy=fp_accuracy)
        if ddof != 0:
            num_elements = (
                self.get_num_elements(axis, item)
                if where is None
                else where.where(blosc2.ones(self.shape, dtype=np.int64), 0).sum(axis=axis, dtype=np.int64)
            )
            out = out * num_elements / (num_elements - ddof)
        out2 = kwargs.pop("out", None)
        if out2 is not None:
            out2[:] = out
            return out2
        if kwargs != {} and not np.isscalar(out):
            out = blosc2.asarray(out, **kwargs)
        return out

    def min(
        self,
        axis=None,
        keepdims=False,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        if where is not None:
            identity = self._reduction_identity("min")
            return self._where_identity_expr(where, identity).min(
                axis=axis, keepdims=keepdims, fp_accuracy=fp_accuracy, **kwargs
            )
        reduce_args = {
            "op": ReduceOp.MIN,
            "op_str": "min",
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def max(
        self,
        axis=None,
        keepdims=False,
        where=None,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        if where is not None:
            identity = self._reduction_identity("max")
            return self._where_identity_expr(where, identity).max(
                axis=axis, keepdims=keepdims, fp_accuracy=fp_accuracy, **kwargs
            )
        reduce_args = {
            "op": ReduceOp.MAX,
            "op_str": "max",
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def any(
        self,
        axis=None,
        keepdims=False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        reduce_args = {
            "op": ReduceOp.ANY,
            "op_str": "any",
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def all(
        self,
        axis=None,
        keepdims=False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        reduce_args = {
            "op": ReduceOp.ALL,
            "op_str": "all",
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def argmax(
        self,
        axis=None,
        keepdims=False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        reduce_args = {
            "op": ReduceOp.ARGMAX,
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def argmin(
        self,
        axis=None,
        keepdims=False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        reduce_args = {
            "op": ReduceOp.ARGMIN,
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def cumulative_sum(
        self,
        axis=None,
        include_initial: bool = False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        reduce_args = {
            "op": ReduceOp.CUMULATIVE_SUM,
            "axis": axis,
            "include_initial": include_initial,
        }
        if self.ndim != 1 and axis is None:
            raise ValueError("axis must be specified for cumulative_sum of non-1D array.")
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def cumulative_prod(
        self,
        axis=None,
        include_initial: bool = False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        reduce_args = {
            "op": ReduceOp.CUMULATIVE_PROD,
            "axis": axis,
            "include_initial": include_initial,
        }
        if self.ndim != 1 and axis is None:
            raise ValueError("axis must be specified for cumulative_prod of non-1D array.")
        return self.compute(_reduce_args=reduce_args, fp_accuracy=fp_accuracy, **kwargs)

    def _eval_constructor(self, expression, constructor, operands):
        """Evaluate a constructor function inside a string expression."""

        def find_args(expr):
            idx = expr.find("(") + 1
            count = 1
            for i, c in enumerate(expr[idx:], start=idx):
                if c == "(":
                    count += 1
                elif c == ")":
                    count -= 1
                if count == 0:
                    return expr[idx:i], i + 1
            raise ValueError("Unbalanced parenthesis in expression")

        # Find the index of the first constructor call.
        match = _find_constructor_call(expression, constructor)
        if match is None:
            raise ValueError(f"Constructor '{constructor}' not found in expression: {expression}")
        idx = match.start()
        # Find the arguments of the constructor function
        try:
            args, idx2 = find_args(expression[idx + len(constructor) :])
        except ValueError as err:
            raise ValueError(f"Unbalanced parenthesis in expression: {expression}") from err
        idx2 = idx + len(constructor) + idx2

        # Give a chance to a possible .reshape() method
        if expression[idx2 : idx2 + len(".reshape(")] == ".reshape(":
            args2, idx3 = find_args(expression[idx2 + len("reshape(") :])
            # Remove a possible shape= from the reshape call (due to rewriting the expression
            # via extract_numpy_scalars(), other variants like .reshape(shape = shape_) work too)
            args2 = args2.replace("shape=", "")
            args = f"{args}, shape={args2}"
            idx2 += len(".reshape") + idx3

        # Evaluate the constructor function
        constructor_func = getattr(blosc2, constructor)
        _globals = {constructor: constructor_func}
        # Add the blosc2 constructors and dtype symbols to the globals
        _globals |= {k: getattr(blosc2, k) for k in constructors}
        _globals |= dtype_symbols
        evalcons = f"{constructor}({args})"

        # Internal constructors will be cached for avoiding multiple computations
        if not hasattr(self, "cons_cache"):
            self.cons_cache = {}
        if evalcons in self.cons_cache:
            return self.cons_cache[evalcons], expression[idx:idx2]
        value = eval(evalcons, _globals, operands)
        self.cons_cache[evalcons] = value

        return value, expression[idx:idx2]

    @staticmethod
    def _is_full_slice(lazy_item):
        """Return True if *lazy_item* is a no-op full slice (() or slice(None))."""
        if isinstance(lazy_item, slice):
            return lazy_item == slice(None)
        if isinstance(lazy_item, tuple):
            return lazy_item == () or all(isinstance(s, slice) and s == slice(None) for s in lazy_item)
        return False

    @staticmethod
    def _collect_flat_indices_from_bool_ndarray(bool_ndarray):
        """Collect flat indices of True positions from a compressed boolean NDArray.

        Uses :meth:`~blosc2.NDArray.iterchunks_info` to skip chunks that are
        special values (e.g. all-False ``ZERO``), avoiding decompression and
        scanning for those chunks.

        Parameters
        ----------
        bool_ndarray: blosc2.NDArray
            A 1D NDArray with boolean dtype.

        Returns
        -------
        np.ndarray
            Flat indices of True positions (int64).
        """
        chunk_len = bool_ndarray.chunks[0]
        all_indices = []

        for info in bool_ndarray.iterchunks_info():
            # Skip special-value chunks that are entirely False
            if info.special == blosc2.SpecialValue.ZERO:
                continue
            if info.special == blosc2.SpecialValue.VALUE:
                if not info.repeated_value:  # repeated_value is False/0
                    continue
                # repeated_value is True: all elements in this chunk are True
                offset = info.nchunk * chunk_len
                all_indices.append(np.arange(offset, offset + chunk_len, dtype=np.int64))
                continue

            # Normal chunk: decompress and scan for True positions
            raw = bool_ndarray.schunk.decompress_chunk(info.nchunk)
            arr = np.frombuffer(raw, dtype=np.bool_)
            # Truncate to the logical chunk size (buffer may include padding)
            if len(arr) > chunk_len:
                arr = arr[:chunk_len]
            idx = np.flatnonzero(arr)
            if len(idx) > 0:
                offset = info.nchunk * chunk_len
                all_indices.append(idx + offset)

        if not all_indices:
            return np.array([], dtype=np.int64)
        return np.concatenate(all_indices)

    def _where_getitem_fastpath(self, item, kwargs):
        """Fast path for where(cond, x) full-slice getitem calls.

        Returns ``None`` when the fast path does not apply.
        """
        from . import indexing

        simple_operand_expr = self.expression.strip("() ") in self.operands
        if not (
            hasattr(self, "_where_args")
            and len(self._where_args) == 1
            and not hasattr(self, "_indices")
            and not hasattr(self, "_order")
            and "_reduce_args" not in kwargs
            and isinstance(self._where_args["_where_x"], blosc2.NDArray)
            and self._is_full_slice(item)
            and not simple_operand_expr
        ):
            return None

        # Preserve index/caching behavior for indexed queries.
        if kwargs.get("_use_index", True) and indexing.will_use_index(self):
            return None

        cond_expr = blosc2.LazyExpr._new_expr(self.expression, self.operands, guess=False)
        if not blosc2.isdtype(cond_expr.dtype, "bool"):
            return None

        target = self._where_args["_where_x"]
        if cond_expr.ndim != 1 or target.ndim != 1:
            return None

        cache_tokens = [indexing.SELF_TARGET_NAME]
        cached_coords = indexing.get_cached_coords(target, self.expression, cache_tokens, None)
        if cached_coords is not None:
            cached_plan = indexing.IndexPlan(
                usable=True, reason="cache-hit", base=target, exact_positions=cached_coords
            )
            return indexing.evaluate_full_query(self._where_args, cached_plan)

        # Evaluate the condition using the miniexpr prefilter (fastest first pass)
        mask = cond_expr.compute((), cparams=_TRANSIENT_MASK_CPARAMS)

        # Collect flat indices by iterating the compressed bool chunks,
        # avoiding a full-mask decompression + count_nonzero + flatnonzero.
        flat_indices = self._collect_flat_indices_from_bool_ndarray(mask)
        indexing.store_cached_coords(target, self.expression, cache_tokens, None, flat_indices)
        plan = indexing.IndexPlan(usable=True, reason="mask-scan", base=target, exact_positions=flat_indices)
        return indexing.evaluate_full_query(self._where_args, plan)

    def _compute_expr(self, item, kwargs):
        if any(method in self.expression for method in eager_funcs):
            # We have reductions in the expression (probably coming from a string lazyexpr)
            # Also includes slice
            _globals = get_expr_globals(self.expression)
            lazy_expr = eval(self.expression, _globals, self.operands)
            if not isinstance(lazy_expr, blosc2.LazyExpr):
                key, mask = process_key(item, lazy_expr.shape)
                # An immediate evaluation happened (e.g. all operands are numpy arrays)
                if hasattr(self, "_where_args"):
                    # We need to apply the where() operation
                    if len(self._where_args) == 1:
                        # We have a single argument
                        where_x = self._where_args["_where_x"]
                        return (where_x[:][lazy_expr])[key]
                    if len(self._where_args) == 2:
                        # We have two arguments
                        where_x = self._where_args["_where_x"]
                        where_y = self._where_args["_where_y"]
                        return np.where(lazy_expr, where_x, where_y)[key]
                out = kwargs.get("_output", None)
                if out is not None:
                    # This is not exactly optimized, but it works for now
                    out[:] = lazy_expr[key]
                    return out
                arr = lazy_expr[key]
                if builtins.sum(mask) > 0:
                    # Correct shape to adjust to NumPy convention.
                    new_shape = tuple(arr.shape[i] for i in range(len(mask)) if not mask[i])
                    arr = np.reshape(arr, new_shape)
                return arr

            return chunked_eval(lazy_expr.expression, lazy_expr.operands, item, **kwargs)

        if any(_has_constructor_call(self.expression, constructor) for constructor in constructors):
            expression = self.expression
            newexpr = expression
            newops = self.operands.copy()
            # We have constructors in the expression (probably coming from a string lazyexpr)
            # Let's replace the constructors with the actual NDArray objects
            for constructor in constructors:
                if not _has_constructor_call(newexpr, constructor):
                    continue
                while _has_constructor_call(newexpr, constructor):
                    # Get the constructor function and replace it by an NDArray object in the operands
                    # Find the constructor call and its arguments
                    value, constexpr = self._eval_constructor(newexpr, constructor, newops)
                    # Add the new operand to the operands; its name will be temporary
                    newop = f"_c{len(newops)}"
                    newops[newop] = value
                    # Replace the constructor call by the new operand
                    newexpr = newexpr.replace(constexpr, newop)

            _globals = get_expr_globals(newexpr)
            lazy_expr = eval(newexpr, _globals, newops)
            if isinstance(lazy_expr, blosc2.NDArray):
                # Almost done (probably the expression is made of only constructors)
                # We only have to define the trivial expression ("o0")
                lazy_expr = blosc2.LazyExpr(new_op=(lazy_expr, None, None))

            return chunked_eval(lazy_expr.expression, lazy_expr.operands, item, **kwargs)

        # Optimization: for where(cond, x) (1-arg) with a boolean condition,
        # stream matching values chunk-by-chunk without materializing the full mask.
        fastpath_result = self._where_getitem_fastpath(item, kwargs)
        if fastpath_result is not None:
            return fastpath_result

        return chunked_eval(self.expression, self.operands, item, **kwargs)

    # TODO: argsort and sort are repeated in LazyUDF; refactor
    def argsort(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        if self.dtype.fields is None:
            raise NotImplementedError("argsort() can only be used with structured arrays")
        if not hasattr(self, "_where_args") or len(self._where_args) != 1:
            raise ValueError("argsort() can only be used with conditions")
        # Build a new lazy array
        lazy_expr = copy.copy(self)
        # ... and assign the new attributes
        lazy_expr._indices = True
        if order:
            lazy_expr._order = order
        # dtype changes to int64
        lazy_expr._dtype = np.dtype(np.int64)
        return lazy_expr

    def sort(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        if self.dtype.fields is None:
            raise NotImplementedError("sort() can only be used with structured arrays")
        if not hasattr(self, "_where_args") or len(self._where_args) != 1:
            raise ValueError("sort() can only be used with conditions")
        # Build a new lazy expression
        lazy_expr = copy.copy(self)
        # ... and assign the new attributes
        if order:
            lazy_expr._order = order
        return lazy_expr

    def will_use_index(self) -> bool:
        """Return whether the current lazy query can use an index."""
        from . import indexing

        return indexing.will_use_index(self)

    def explain(self) -> dict:
        """Explain how this lazy query will be executed.

        Returns a dictionary describing the planner decision for the current
        query. Typical fields include whether an index will be used, the chosen
        index kind and level, candidate counts, and the lookup path selected
        for ``full`` indexes.

        Returns:
            dict: Query planning metadata for the current expression.

        Examples:
            >>> import numpy as np
            >>> import blosc2
            >>> arr = blosc2.asarray(np.arange(10))
            >>> _ = arr.create_index(kind=blosc2.IndexKind.FULL)
            >>> expr = blosc2.lazyexpr("(a >= 3) & (a < 6)", {"a": arr}).where(arr)
            >>> info = expr.explain()
            >>> info["will_use_index"]
            True
            >>> info["kind"]
            'full'
        """
        from . import indexing

        return indexing.explain_query(self)

    def compute(
        self,
        item=(),
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        jit=None,
        jit_backend: str | None = None,
        **kwargs,
    ) -> blosc2.NDArray:
        # When NumPy ufuncs are called, the user may add an `out` parameter to kwargs
        if "out" in kwargs:  # use provided out preferentially
            kwargs["_output"] = kwargs.pop("out")
        elif hasattr(self, "_output"):
            kwargs["_output"] = self._output

        if "ne_args" in kwargs:
            kwargs["_ne_args"] = kwargs.pop("ne_args")
        if hasattr(self, "_ne_args"):
            kwargs["_ne_args"] = self._ne_args
        if hasattr(self, "_where_args"):
            kwargs["_where_args"] = self._where_args
        kwargs.setdefault("fp_accuracy", fp_accuracy)
        jit, jit_backend = _jit_from_env(jit, jit_backend)
        if jit is not None:
            kwargs["jit"] = jit
        if jit_backend is not None:
            kwargs["jit_backend"] = jit_backend
        kwargs["dtype"] = self.dtype
        kwargs["shape"] = self.shape
        if hasattr(self, "_indices"):
            kwargs["_indices"] = self._indices
        if hasattr(self, "_order"):
            kwargs["_order"] = self._order
        result = self._compute_expr(item, kwargs)
        if "_order" in kwargs and "_indices" not in kwargs:
            # We still need to apply the index in result
            x = self._where_args["_where_x"]
            result = x[result]  # always a numpy array; TODO: optimize this for _getitem not in kwargs
        if (
            "_getitem" not in kwargs
            and "_output" not in kwargs
            and "_reduce_args" not in kwargs
            and not isinstance(result, blosc2.NDArray)
        ):
            # Get rid of all the extra kwargs that are not accepted by blosc2.asarray
            kwargs_not_accepted = {
                "_where_args",
                "_indices",
                "_order",
                "_ne_args",
                "_use_index",
                "dtype",
                "shape",
                "fp_accuracy",
            }
            kwargs = {key: value for key, value in kwargs.items() if key not in kwargs_not_accepted}
            result = blosc2.asarray(result, **kwargs)
        return result

    def __getitem__(self, item):
        kwargs = {"_getitem": True}
        result = self.compute(item, **kwargs)
        # Squeeze single-element dimensions when indexing with integers
        # See e.g. examples/ndarray/animated_plot.py
        if isinstance(item, int) or (hasattr(item, "__iter__") and any(isinstance(i, int) for i in item)):
            result = result.squeeze(axis=tuple(i for i in range(result.ndim) if result.shape[i] == 1))
        return result

    def slice(self, item):
        return self.compute(item)  # should do a slice since _getitem = False

    def __str__(self):
        return f"{self.expression}"

    @property
    def info(self):
        return InfoReporter(self)

    @property
    def info_items(self):
        items = []
        items += [("type", f"{self.__class__.__name__}")]
        items += [("expression", self.expression)]
        opsinfo = {
            key: str(value) if value.schunk.urlpath is None else value.schunk.urlpath
            for key, value in self.operands.items()
        }
        items += [("operands", opsinfo)]
        items += [("shape", self.shape)]
        items += [("dtype", self.dtype)]
        return items

    def save(self, urlpath=None, **kwargs):
        if urlpath is None:
            raise ValueError("To save a LazyArray you must provide an urlpath")

        kwargs["urlpath"] = urlpath
        kwargs["mode"] = "w"  # always overwrite the file in urlpath
        self._to_b2object_carrier(**kwargs)

    def to_cframe(self) -> bytes:
        return self._to_b2object_carrier().to_cframe()

    def _to_b2object_carrier(self, **kwargs):
        expression = self.expression_tosave if hasattr(self, "expression_tosave") else self.expression
        operands_ = self.operands_tosave if hasattr(self, "operands_tosave") else self.operands
        validate_expr(expression)

        payload = {"kind": "lazyexpr", "version": 1, "expression": expression, "operands": {}}
        carrier_urlpath = kwargs.get("urlpath")
        carrier_parent = Path(carrier_urlpath).parent if carrier_urlpath is not None else None
        for key, value in operands_.items():
            if isinstance(value, blosc2.C2Array):
                payload["operands"][key] = encode_b2object_payload(value)
                continue
            if isinstance(value, blosc2.Proxy):
                value = value._cache
            ref = getattr(value, "_blosc2_ref", None)
            if isinstance(ref, blosc2.Ref):
                payload["operands"][key] = ref.to_dict()
                continue
            if not hasattr(value, "schunk"):
                raise ValueError(
                    "To save a LazyArray, all operands must be blosc2.NDArray or blosc2.C2Array objects"
                )
            if value.schunk.urlpath is None:
                raise ValueError("To save a LazyArray, all operands must be stored on disk/network")
            operand_urlpath = Path(value.schunk.urlpath)
            if carrier_parent is not None and not operand_urlpath.is_absolute():
                ref_urlpath = operand_urlpath.as_posix()
            elif carrier_parent is not None:
                try:
                    ref_urlpath = operand_urlpath.relative_to(carrier_parent).as_posix()
                except ValueError:
                    ref_urlpath = operand_urlpath.as_posix()
            else:
                ref_urlpath = operand_urlpath.as_posix()
            payload["operands"][key] = {"kind": "urlpath", "version": 1, "urlpath": ref_urlpath}
        array = make_b2object_carrier(
            "lazyexpr",
            self.shape,
            self.dtype,
            chunks=self.chunks,
            blocks=self.blocks,
            **kwargs,
        )
        write_b2object_payload(array, payload)
        write_b2object_user_vlmeta(array, self._get_user_vlmeta())
        return array

    @classmethod
    def _new_expr(cls, expression, operands, guess, out=None, where=None, ne_args=None):
        # Validate the expression
        validate_expr(expression)
        expression = convert_to_slice(expression)
        chunks, blocks = None, None
        if guess:
            # The expression has been validated, so we can evaluate it
            # in guessing mode to avoid computing reductions
            # Extract possible numpy scalars
            _expression, local_vars = extract_numpy_scalars(expression)
            _operands = operands | local_vars
            # Check that operands are proper Operands, LazyArray or scalars; if not, convert to NDArray objects
            for op, val in _operands.items():
                if not (isinstance(val, blosc2.Operand | np.ndarray) or np.isscalar(val)):
                    _operands[op] = blosc2.SimpleProxy(val)
            # for scalars just return value (internally converts to () if necessary)
            opshapes = {k: v if not hasattr(v, "shape") else v.shape for k, v in _operands.items()}
            _shape = infer_shape(_expression, opshapes)  # infer shape, includes constructors
            # have to handle slices since a[10] on a dummy variable of shape (1,1) doesn't work
            desliced_expr, desliced_ops = extract_and_replace_slices(_expression, _operands)
            # substitutes with dummy operands (cheap for reductions) and
            # defaults to blosc2 functions (cheap for constructors)
            new_expr = _numpy_eval_expr(desliced_expr, desliced_ops, prefer_blosc=True)
            _dtype = new_expr.dtype if hasattr(new_expr, "dtype") else np.dtype(type(new_expr))
            if isinstance(new_expr, blosc2.LazyExpr):
                # DO NOT restore the original expression and operands
                # Instead rebase operands and restore only constructors
                expression_, operands_ = conserve_functions(
                    _expression, _operands, new_expr.operands | local_vars
                )
            elif _shape == () and not _operands:  # passed scalars
                expression_ = "o0"
                operands_ = {"o0": ne_evaluate(_expression)}
            else:
                # An immediate evaluation happened
                # (e.g. all operands are numpy arrays or constructors)
                # or passed "a", "a[:10]", 'sum(a)'
                expression_, operands_ = conserve_functions(_expression, _operands, local_vars)
                if hasattr(new_expr, "chunks") and new_expr.chunks != (1,) * len(_shape):
                    # for constructors with chunks in kwargs, chunks will be specified
                    # for general expression new_expr is just with dummy scalar variables (so ignore)
                    chunks = new_expr.chunks
                    blocks = new_expr.blocks
            new_expr = cls(None)
            new_expr.expression = f"({expression_})"  # force parenthesis
            new_expr.operands = operands_
            new_expr.expression_tosave = expression
            new_expr.operands_tosave = operands
            # Cache the dtype and shape (should be immutable)
            new_expr._dtype = _dtype
            new_expr._shape = _shape
            if chunks is not None and blocks is not None:
                new_expr._chunks, new_expr._blocks = chunks, blocks
        else:
            # Create a new LazyExpr object
            new_expr = cls(None)
            new_expr.expression = expression
            new_expr.operands = operands
        if out is not None:
            new_expr._output = out
        if where is not None:
            new_expr._where_args = where
        new_expr._ne_args = ne_args
        return new_expr


class LazyUDF(LazyArray):
    def __init__(
        self, func, inputs, dtype, shape=None, chunked_eval=True, jit=None, jit_backend=None, **kwargs
    ):
        # After this, all the inputs should be np.ndarray or NDArray objects
        self.inputs = convert_inputs(inputs)
        # Get res shape
        if shape is None:
            self._shape = compute_broadcast_shape(self.inputs)
            if self._shape is None:
                raise ValueError(
                    "If all inputs are scalars, pass a `shape` argument to indicate the output shape"
                )
        else:
            self._shape = shape

        self.kwargs = kwargs
        self.kwargs["dtype"] = dtype
        self.kwargs["shape"] = self._shape
        self.kwargs["jit"] = jit
        self.kwargs["jit_backend"] = jit_backend
        in_place = kwargs.get("in_place", False)
        self.kwargs["in_place"] = in_place
        self._dtype = dtype
        self.func = func
        if isinstance(self.func, DSLKernel) and self.func.dsl_error is not None:
            udf_name = getattr(self.func.func, "__name__", self.func.__name__)
            raise DSLSyntaxError(f"Invalid DSL kernel '{udf_name}'.\n{self.func.dsl_error}") from None

        # Prepare internal array for __getitem__
        # Deep copy the kwargs to avoid modifying them
        kwargs_getitem = copy.deepcopy(self.kwargs)
        # Cannot use multithreading when applying a postfilter, dparams['nthreads'] ignored
        dparams = kwargs_getitem.get("dparams", {})
        if isinstance(dparams, dict):
            dparams["nthreads"] = 1
        else:
            raise TypeError("dparams should be a dictionary")
        kwargs_getitem["dparams"] = dparams

        if isinstance(self.func, DSLKernel) and self.func.input_names:
            # DSL kernels are using input names that are extracted from params as a list,
            # and we need to use them for matching variables in miniexpr
            # (instead of the 'o{%d}' notation).
            self.inputs_dict = dict(zip(self.func.input_names, self.inputs, strict=True))
        else:
            self.inputs_dict = {f"o{i}": obj for i, obj in enumerate(self.inputs)}

    @property
    def dtype(self):
        return self._dtype

    @property
    def ndim(self) -> int:
        return len(self.shape)

    @property
    def shape(self):
        return self._shape

    @property
    def info(self):
        return InfoReporter(self)

    @property
    def info_items(self):
        inputs = {}
        for key, value in self.inputs_dict.items():
            if isinstance(value, blosc2.Array):
                inputs[key] = f"<{value.__class__.__name__}> {value.shape} {value.dtype}"
            else:
                inputs[key] = str(value)
        return [
            ("type", f"{self.__class__.__name__}"),
            ("inputs", inputs),
            ("shape", self.shape),
            ("dtype", self.dtype),
        ]

    @property
    def chunks(self):
        if hasattr(self, "_chunks"):
            return self._chunks
        if not self.inputs_dict:
            req_chunks = self.kwargs.get("chunks")
            req_blocks = self.kwargs.get("blocks")
            self._chunks, self._blocks = compute_chunks_blocks(
                self.shape, req_chunks, req_blocks, dtype=self.dtype
            )
            return self._chunks
        shape, self._chunks, self._blocks, fast_path = validate_inputs(
            self.inputs_dict, getattr(self, "_out", None)
        )
        if not hasattr(self, "_shape"):
            self._shape = shape
        if self._shape != shape:  # validate inputs only works for elementwise funcs so returned shape might
            fast_path = False  # be incompatible with true output shape
        if not fast_path:
            # Not using the fast path, so we need to compute the chunks/blocks automatically
            self._chunks, self._blocks = compute_chunks_blocks(self.shape, None, None, dtype=self.dtype)
        return self._chunks

    @property
    def blocks(self):
        if hasattr(self, "_blocks"):
            return self._blocks
        if not self.inputs_dict:
            req_chunks = self.kwargs.get("chunks")
            req_blocks = self.kwargs.get("blocks")
            self._chunks, self._blocks = compute_chunks_blocks(
                self.shape, req_chunks, req_blocks, dtype=self.dtype
            )
            return self._blocks
        shape, self._chunks, self._blocks, fast_path = validate_inputs(
            self.inputs_dict, getattr(self, "_out", None)
        )
        if not hasattr(self, "_shape"):
            self._shape = shape
        if self._shape != shape:  # validate inputs only works for elementwise funcs so returned shape might
            fast_path = False  # be incompatible with true output shape
        if not fast_path:
            # Not using the fast path, so we need to compute the chunks/blocks automatically
            self._chunks, self._blocks = compute_chunks_blocks(self.shape, None, None, dtype=self.dtype)
        return self._blocks

    # TODO: argsort and sort are repeated in LazyExpr; refactor
    def argsort(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        if self.dtype.fields is None:
            raise NotImplementedError("argsort() can only be used with structured arrays")
        if not hasattr(self, "_where_args") or len(self._where_args) != 1:
            raise ValueError("argsort() can only be used with conditions")
        # Build a new lazy array
        lazy_expr = copy.copy(self)
        # ... and assign the new attributes
        lazy_expr._indices = True
        if order:
            lazy_expr._order = order
        # dtype changes to int64
        lazy_expr._dtype = np.dtype(np.int64)
        return lazy_expr

    def sort(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        if self.dtype.fields is None:
            raise NotImplementedError("sort() can only be used with structured arrays")
        if not hasattr(self, "_where_args") or len(self._where_args) != 1:
            raise ValueError("sort() can only be used with conditions")
        # Build a new lazy expression
        lazy_expr = copy.copy(self)
        # ... and assign the new attributes
        if order:
            lazy_expr._order = order
        return lazy_expr

    def compute(
        self,
        item=(),
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        jit=None,
        jit_backend=None,
        **kwargs,
    ):
        # Get kwargs
        if kwargs is None:
            kwargs = {}
        # Do copy to avoid modifying the original parameters
        aux_kwargs = copy.deepcopy(self.kwargs)

        # Update is not recursive
        aux_cparams = aux_kwargs.get("cparams", {})
        if isinstance(aux_cparams, blosc2.CParams):
            # Convert to dictionary
            aux_cparams = asdict(aux_cparams)
        cparams = kwargs.get("cparams", {})
        if isinstance(cparams, blosc2.CParams):
            # Convert to dictionary
            cparams = asdict(cparams)
        aux_cparams.update(cparams)
        aux_kwargs["cparams"] = aux_cparams

        aux_dparams = aux_kwargs.get("dparams", {})
        if isinstance(aux_dparams, blosc2.DParams):
            # Convert to dictionary
            aux_dparams = asdict(aux_dparams)
        dparams = kwargs.get("dparams", {})
        if isinstance(dparams, blosc2.DParams):
            # Convert to dictionary
            dparams = asdict(dparams)
        aux_dparams.update(dparams)
        aux_kwargs["dparams"] = aux_dparams

        _ = kwargs.pop("cparams", None)
        _ = kwargs.pop("dparams", None)
        if jit is not None:
            aux_kwargs["jit"] = jit
        if jit_backend is not None:
            aux_kwargs["jit_backend"] = jit_backend
        urlpath = kwargs.get("urlpath")
        if urlpath is not None and urlpath == aux_kwargs.get(
            "urlpath",
        ):
            raise ValueError("Cannot use the same urlpath for LazyArray and eval NDArray")
        _ = aux_kwargs.pop("urlpath", None)

        if "out" in kwargs:  # use provided out preferentially
            aux_kwargs["_output"] = kwargs.pop("out")
        elif hasattr(self, "_output"):
            aux_kwargs["_output"] = self._output
        aux_kwargs.update(kwargs)

        # aux_kwargs includes self.shape and self.dtype
        return chunked_eval(self.func, self.inputs_dict, item, _getitem=False, **aux_kwargs)

    def __getitem__(self, item):
        return chunked_eval(self.func, self.inputs_dict, item, _getitem=True, **self.kwargs)

    def save(self, urlpath=None, **kwargs):
        """
        Save the :ref:`LazyUDF` on disk.

        Parameters
        ----------
        urlpath: str
            The path to the file where the LazyUDF will be stored.
        kwargs: Any, optional
            Keyword arguments that are supported by the :func:`empty` constructor.

        Returns
        -------
        out: None

        Notes
        -----
        * All operands must be :ref:`NDArray` or :ref:`C2Array` objects stored on
          disk or a remote server (i.e. they must have a ``urlpath``).
        * When the :ref:`LazyUDF` wraps a :func:`blosc2.dsl_kernel`-decorated
          function, the DSL source is preserved verbatim in the saved metadata.
          On reload via :func:`blosc2.open`, the function is restored as a full
          :class:`~blosc2.dsl_kernel.DSLKernel` so the miniexpr JIT fast path
          remains available without any extra work from the caller.
        """
        if urlpath is None:
            raise ValueError("To save a LazyArray you must provide an urlpath")

        kwargs["urlpath"] = urlpath
        kwargs["mode"] = "w"  # always overwrite the file in urlpath
        try:
            self._to_b2object_carrier(**kwargs)
        except (TypeError, ValueError):
            meta = kwargs.get("meta", {})
            meta["LazyArray"] = LazyArrayEnum.UDF.value
            kwargs["meta"] = meta

            # Create an empty array; useful for providing the shape and dtype of the outcome
            array = blosc2.empty(shape=self.shape, dtype=self.dtype, **kwargs)

            # Save the expression and operands in the metadata
            operands = {}
            operands_ = self.inputs_dict
            for i, (_key, value) in enumerate(operands_.items()):
                pos_key = f"o{i}"  # always use positional keys for consistent loading
                if isinstance(value, blosc2.C2Array):
                    operands[pos_key] = {
                        "path": str(value.path),
                        "urlbase": value.urlbase,
                    }
                    continue
                if isinstance(value, blosc2.Proxy):
                    # Take the required info from the Proxy._cache container
                    value = value._cache
                if not hasattr(value, "schunk"):
                    raise ValueError(
                        "To save a LazyArray, all operands must be blosc2.NDArray or blosc2.C2Array objects"
                    ) from None
                if value.schunk.urlpath is None:
                    raise ValueError(
                        "To save a LazyArray, all operands must be stored on disk/network"
                    ) from None
                operands[pos_key] = value.schunk.urlpath
            udf_func = self.func.func if isinstance(self.func, DSLKernel) else self.func
            udf_name = getattr(udf_func, "__name__", self.func.__name__)
            try:
                udf_source = textwrap.dedent(inspect.getsource(udf_func)).lstrip()
            except Exception:
                udf_source = None
            meta = {
                "UDF": udf_source,
                "operands": operands,
                "name": udf_name,
            }
            if isinstance(self.func, DSLKernel) and self.func.dsl_source is not None:
                meta["dsl_source"] = self.func.dsl_source
            array.schunk.vlmeta["_LazyArray"] = meta
            write_b2object_user_vlmeta(array, self._get_user_vlmeta())

    def to_cframe(self) -> bytes:
        return self._to_b2object_carrier().to_cframe()

    def _to_b2object_carrier(self, **kwargs):
        payload = encode_b2object_payload(self)
        if payload is None:
            raise TypeError("Persistent Blosc2 object payload is not supported for this LazyUDF")

        carrier_urlpath = kwargs.get("urlpath")
        carrier_parent = Path(carrier_urlpath).parent if carrier_urlpath is not None else None
        for operand_payload in payload["operands"].values():
            if operand_payload["kind"] not in {"urlpath", "dictstore_key"}:
                continue
            operand_urlpath = Path(operand_payload["urlpath"])
            if carrier_parent is not None and not operand_urlpath.is_absolute():
                ref_urlpath = operand_urlpath.as_posix()
            elif carrier_parent is not None:
                try:
                    ref_urlpath = operand_urlpath.relative_to(carrier_parent).as_posix()
                except ValueError:
                    ref_urlpath = operand_urlpath.as_posix()
            else:
                ref_urlpath = operand_urlpath.as_posix()
            operand_payload["urlpath"] = ref_urlpath

        array = make_b2object_carrier(
            "lazyudf",
            self.shape,
            self.dtype,
            chunks=self.chunks,
            blocks=self.blocks,
            **kwargs,
        )
        write_b2object_payload(array, payload)
        write_b2object_user_vlmeta(array, self._get_user_vlmeta())
        return array


def _numpy_eval_expr(expression, operands, prefer_blosc=False):
    npops = {
        key: np.ones(np.ones(len(value.shape), dtype=int), dtype=value.dtype)
        if hasattr(value, "shape")
        else value
        for key, value in operands.items()
    }
    if prefer_blosc:
        # convert blosc arrays to small dummies
        ops = {
            key: blosc2.ones((1,) * len(value.shape), dtype=value.dtype)
            if hasattr(value, "chunks")
            else value  # some of these could be numpy arrays
            for key, value in operands.items()
        }
        # change numpy arrays
        ops = {
            key: np.ones((1,) * len(value.shape), dtype=value.dtype)
            if isinstance(value, np.ndarray)
            else value
            for key, value in ops.items()
        }
    else:  # wasm pathway assumes numpy arrs
        ops = npops

    # Create a globals dict with blosc2 version of functions preferentially
    # (default to numpy func if not implemented in blosc2)
    if prefer_blosc:
        _globals = get_expr_globals(expression)
        _globals |= dtype_symbols
    else:
        _globals = safe_numpy_globals
    try:
        _out = eval(expression, _globals, ops)
    except RuntimeWarning:
        # Sometimes, numpy gets a RuntimeWarning when evaluating expressions
        # with synthetic operands (1's). Let's try with numexpr, which is not so picky
        # about this.
        ops = npops if blosc2.IS_WASM else ops
        _out = ne_evaluate(expression, local_dict=ops)
    return _out


def lazyudf(
    func: Callable[[tuple, np.ndarray, tuple[int]], None],
    inputs: Sequence[Any] | None,
    dtype: np.dtype | None = None,
    shape: tuple | list | None = None,
    chunked_eval: bool = True,
    jit: bool | None = None,
    jit_backend: str | None = None,
    **kwargs: Any,
) -> LazyUDF:
    """
    Get a LazyUDF from a python user-defined function.

    Parameters
    ----------
    func: Python function
        The user-defined function to apply to each block. This function will
        always receive the following parameters:

        - `inputs_tuple`: A tuple containing the corresponding slice for the block of each input
          in :paramref:`inputs`.
        - `output`: The buffer to be filled as a multidimensional numpy.ndarray.
        - `offset`: The multidimensional offset corresponding to the start of the block
          being computed. Example signature::

            def myudf(inputs_tuple, output, offset):
                x, y = inputs_tuple
                ...
                output[:] = result
    inputs: Sequence[Any] or None
        The sequence of inputs. Besides objects compliant with the blosc2.Array protocol,
        any other object is supported too, and it will be passed as-is to the
        user-defined function. If not needed, this can be empty, but `shape` must
        be provided.
    dtype: np.dtype, optional
        The resulting ndarray dtype in NumPy format.  When omitted and *func*
        is a :class:`DSLKernel`, the dtype is inferred by NumPy type promotion
        of the input dtypes.  For type-changing kernels (comparisons, casts)
        pass *dtype* explicitly.  Required for plain Python UDFs.
    shape: tuple, optional
        The shape of the resulting array. If None, the shape will be guessed from inputs.
    chunked_eval: bool, optional
        Whether to evaluate the function in chunks or not (blocks).
    jit: bool or None, optional
        JIT policy for miniexpr-backed execution:
        ``None`` uses default behavior (currently, JIT is tried out), ``True`` prefers JIT, ``False`` disables JIT.
    jit_backend: {"tcc", "cc", "js"} or None, optional
        JIT backend selection. ``None`` uses backend defaults (miniexpr "tcc"), except under
        WebAssembly where — unless ``jit=False`` — it *prefers* ``"js"`` for transpilable
        float DSL kernels and falls back to miniexpr otherwise (``jit=True`` prefers ``"js"``
        too, since it is JIT-compiled by the JS engine). ``"tcc"`` forces libtcc, ``"cc"``
        forces the C compiler backend, and ``"js"`` transpiles a :func:`blosc2.dsl_kernel`
        to JavaScript (browser/Pyodide only; raises elsewhere).
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.
        These arguments will be used by the :meth:`LazyArray.__getitem__` and
        :meth:`LazyArray.compute` methods. The
        last one will ignore the `urlpath` parameter passed in this function.
        In addition, one may provide ``in_place``, a bool (default False), which indicates whether
        the function should modify the output directly (rather than chunks of the output, which
        are later written to output). Example::

            def inplace_udf(inputs_tuple, output, offset):
                x, y = inputs_tuple
                ...
                out[3] += 1

    Returns
    -------
    out: :ref:`LazyUDF`
        A :ref:`LazyUDF` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> dtype = np.float64
    >>> shape = [3, 3]
    >>> size = shape[0] * shape[1]
    >>> a = np.linspace(0, 10, num=size, dtype=dtype).reshape(shape)
    >>> b = np.linspace(10, 20, num=size, dtype=dtype).reshape(shape)
    >>> a1 = blosc2.asarray(a)
    >>> b1 = blosc2.asarray(b)
    >>> # Define a user-defined function that will be applied to each block of data
    >>> def my_function(inputs_tuple, output, offset):
    >>>     a, b = inputs_tuple
    >>>     output[:] = a + b
    >>> # Create a LazyUDF object using the user-defined function
    >>> lazy_udf = blosc2.lazyudf(my_function, [a1, b1], dtype)
    >>> type(lazy_udf)
    <class 'blosc2.lazyexpr.LazyUDF'>
    >>> f"Result of LazyUDF evaluation: {lazy_udf[:]}"
    Result of LazyUDF evaluation:
            [[10.  12.5 15. ]
            [17.5 20.  22.5]
            [25.  27.5 30. ]]
    """
    if isinstance(func, DSLKernel) and func.dsl_error is not None:
        udf_name = getattr(func.func, "__name__", func.__name__)
        raise DSLSyntaxError(f"Invalid DSL kernel '{udf_name}'.\n{func.dsl_error}") from None
    if dtype is None:
        if isinstance(func, DSLKernel):
            dep_dtypes = [arr.dtype for arr in (inputs or []) if hasattr(arr, "dtype")]
            if not dep_dtypes:
                raise TypeError(
                    "Cannot infer dtype for DSL kernel with no array inputs; pass dtype= explicitly."
                )
            dtype = np.result_type(*dep_dtypes)
        else:
            raise TypeError("dtype is required for non-DSL UDFs.")
    return LazyUDF(func, inputs, dtype, shape, chunked_eval, jit, jit_backend, **kwargs)


def seek_operands(names, local_dict=None, global_dict=None, _frame_depth: int = 2):
    """
    Get the arguments based on the names.
    """
    call_frame = sys._getframe(_frame_depth)

    clear_local_dict = False
    if local_dict is None:
        local_dict = call_frame.f_locals
        clear_local_dict = True
    try:
        frame_globals = call_frame.f_globals
        if global_dict is None:
            global_dict = frame_globals

        # If `call_frame` is the top frame of the interpreter we can't clear its
        # `local_dict`, because it is actually the `global_dict`.
        clear_local_dict = clear_local_dict and frame_globals is not local_dict

        op_dict = {}
        for name in names:
            try:
                a = local_dict[name]
            except KeyError:
                a = global_dict[name]
            op_dict[name] = a
    finally:
        # If we generated local_dict via an explicit reference to f_locals,
        # clear the dict to prevent creating extra ref counts in the caller's scope
        if clear_local_dict and hasattr(local_dict, "clear"):
            local_dict.clear()

    return op_dict


def lazyexpr(
    expression: str | bytes | LazyArray | blosc2.NDArray,
    operands: dict | None = None,
    out: blosc2.Array = None,
    where: tuple | list | None = None,
    local_dict: dict | None = None,
    global_dict: dict | None = None,
    ne_args: dict | None = None,
    _frame_depth: int = 2,
) -> LazyExpr:
    """
    Get a LazyExpr from an expression.

    Parameters
    ----------
    expression: str or bytes or LazyExpr or NDArray
        The expression to evaluate. This can be any valid expression that numexpr
        can ingest. If a LazyExpr is passed, the expression will be
        updated with the new operands.
    operands: dict[blosc2.Array], optional
        The dictionary with operands. Supported values are Python scalars,
        or any instance that is blosc2.Array compliant.
        If None, the operands will be seeked in the local and global dictionaries.
    out: blosc2.Array, optional
        The output array where the result will be stored. If not provided,
        a new NumPy array will be created and returned.
    where: tuple, list, optional
        A sequence of arguments for the where clause in the expression.
    local_dict: dict, optional
        The local dictionary to use when looking for operands in the expression.
        If not provided, the local dictionary of the caller will be used.
    global_dict: dict, optional
        The global dictionary to use when looking for operands in the expression.
        If not provided, the global dictionary of the caller will be used.
    ne_args: dict, optional
        Additional arguments to be passed to `numexpr.evaluate()` function.
    _frame_depth: int, optional
        The depth of the frame to use when looking for operands in the expression.
        The default value is 2.


    Returns
    -------
    out: :ref:`LazyExpr`
        A :ref:`LazyExpr` is returned.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> dtype = np.float64
    >>> shape = [3, 3]
    >>> size = shape[0] * shape[1]
    >>> a = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
    >>> b = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
    >>> a1 = blosc2.asarray(a)
    >>> a1[:]
    [[0.    0.625 1.25 ]
    [1.875 2.5   3.125]
    [3.75  4.375 5.   ]]
    >>> b1 = blosc2.asarray(b)
    >>> expr = 'a * b + 2'
    >>> operands = { 'a': a1, 'b': b1 }
    >>> lazy_expr = blosc2.lazyexpr(expr, operands=operands)
    >>> f"Lazy expression created: {lazy_expr}"
    Lazy expression created: a * b + 2
    >>> lazy_expr[:]
    [[ 2.        2.390625  3.5625  ]
    [ 5.515625  8.25     11.765625]
    [16.0625   21.140625 27.      ]]
    """
    if isinstance(expression, LazyExpr):
        if operands is not None:
            expression.operands.update(operands)
        if out is not None:
            expression._output = out
        expression._ne_args = ne_args
        if where is not None:
            where_args = {"_where_x": where[0], "_where_y": where[1]}
            expression._where_args = where_args
        return expression
    elif isinstance(expression, blosc2.NDArray):
        operands = {"o0": expression}
        return LazyExpr._new_expr("o0", operands, guess=False, out=out, where=where, ne_args=ne_args)

    if operands is None:
        # Try to get operands from variables in the stack
        operand_set = get_expr_operands(expression)
        # If no operands are found, raise an error
        if operand_set:
            # Look for operands in the stack
            operands = seek_operands(operand_set, local_dict, global_dict, _frame_depth=_frame_depth)
        else:
            # No operands found in the expression. Maybe a constructor?
            constructor = any(_has_constructor_call(expression, constructor) for constructor in constructors)
            if not constructor:
                raise ValueError("No operands nor constructors found in the expression")
            # _new_expr will take care of the constructor, but needs an empty dict in operands
            operands = {}

    return LazyExpr._new_expr(expression, operands, guess=True, out=out, where=where, ne_args=ne_args)


def _reconstruct_lazyudf(expr, lazyarray, operands_dict, array):
    """Reconstruct a LazyUDF (including DSL kernels) from saved metadata."""
    local_ns = {}
    name = lazyarray["name"]
    filename = f"<{name}>"  # any unique name
    SAFE_GLOBALS = {
        "__builtins__": {k: v for k, v in builtins.__dict__.items() if k != "__import__"},
        "np": np,
        "blosc2": blosc2,
    }
    if blosc2._HAS_NUMBA:
        import numba

        SAFE_GLOBALS["numba"] = numba

    # Register the source so inspect can find it
    linecache.cache[filename] = (len(expr), None, expr.splitlines(True), filename)

    exec(compile(expr, filename, "exec"), SAFE_GLOBALS, local_ns)
    func = local_ns[name]
    # If the saved LazyUDF was a DSL kernel, re-wrap and restore the dsl_source
    if "dsl_source" in lazyarray:
        if not isinstance(func, DSLKernel):
            func = DSLKernel(func)
        if func.dsl_source is None:
            # Re-extraction from linecache failed; use the saved verbatim dsl_source
            func.dsl_source = lazyarray["dsl_source"]
    # TODO: make more robust for general kwargs (not just cparams)
    return blosc2.lazyudf(
        func,
        tuple(operands_dict[f"o{n}"] for n in range(len(operands_dict))),
        shape=array.shape,
        dtype=array.dtype,
        cparams=array.cparams,
    )


def open_lazyarray(array):
    value = array.schunk.meta["LazyArray"]
    lazyarray = array.schunk.vlmeta["_LazyArray"]
    if value == LazyArrayEnum.Expr.value:
        expr = lazyarray["expression"]
    elif value == LazyArrayEnum.UDF.value:
        expr = lazyarray["UDF"]
    else:
        raise ValueError("Argument `array` is not LazyExpr or LazyUDF instance.")

    operands = lazyarray["operands"]
    parent_path = Path(array.schunk.urlpath).parent
    operands_dict = {}
    missing_ops = {}
    for key, v in operands.items():
        if isinstance(v, str):
            v = parent_path / v
            try:
                op = blosc2.open(v, mode="r")
            except FileNotFoundError:
                missing_ops[key] = v
            else:
                operands_dict[key] = op
        elif isinstance(v, dict):
            # C2Array
            operands_dict[key] = blosc2.C2Array(
                pathlib.Path(v["path"]).as_posix(),
                urlbase=v["urlbase"],
            )
        else:
            raise TypeError("Error when retrieving the operands")

    if missing_ops:
        exc = exceptions.MissingOperands(expr, missing_ops)
        exc.expr = expr
        exc.missing_ops = missing_ops
        raise exc

    # LazyExpr
    if value == LazyArrayEnum.Expr.value:
        new_expr = LazyExpr._new_expr(expr, operands_dict, guess=True, out=None, where=None)
    elif value == LazyArrayEnum.UDF.value:
        new_expr = _reconstruct_lazyudf(expr, lazyarray, operands_dict, array)

    # Make the array info available for the user (only available when opened from disk)
    new_expr.array = array
    # We want to expose schunk too, so that .info() can be used on the LazyArray
    new_expr.schunk = array.schunk
    new_expr._set_user_vlmeta(read_b2object_user_vlmeta(array), sync=False)
    return new_expr


# Mimim numexpr's evaluate function
def evaluate(
    ex: str,
    local_dict: dict | None = None,
    global_dict: dict | None = None,
    out: blosc2.Array = None,
    **kwargs: Any,
) -> blosc2.Array:
    """
    Evaluate a string expression using the Blosc2 compute engine.

    This is a drop-in replacement for `numexpr.evaluate()`, but using the
    Blosc2 compute engine.  This allows for:

    1) Use more functionality (e.g. reductions) than numexpr.
    2) Follow casting rules of NumPy more closely.
    3) Use both NumPy arrays and Blosc2 NDArrays in the same expression.

    As NDArrays can be on-disk, the expression can be evaluated without loading
    the whole array into memory (i.e. using an out-of-core approach).

    Parameters
    ----------
    ex: str
        The expression to evaluate.
    local_dict: dict, optional
        The local dictionary to use when looking for operands in the expression.
        If not provided, the local dictionary of the caller will be used.
    global_dict: dict, optional
        The global dictionary to use when looking for operands in the expression.
        If not provided, the global dictionary of the caller will be used.
    out: blosc2.Array, optional
        The output array where the result will be stored. If not provided,
        a new NumPy array will be created and returned.
    kwargs: Any, optional
        Additional arguments to be passed to `numexpr.evaluate()` function.

    Returns
    -------
    out: blosc2.Array
        The result of the expression evaluation.  If out is provided, the result
        will be stored in out and returned at the same time.

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> dtype = np.float64
    >>> shape = [3, 3]
    >>> size = shape[0] * shape[1]
    >>> a = np.linspace(0, 5, num=size, dtype=dtype).reshape(shape)
    >>> b = blosc2.linspace(0, 5, num=size, dtype=dtype, shape=shape)
    >>> expr = 'a * b + 2'
    >>> out = blosc2.evaluate(expr)
    >>> out
    [[ 2.        2.390625  3.5625  ]
    [ 5.515625  8.25     11.765625]
    [16.0625   21.140625 27.      ]]
    """
    lexpr = lazyexpr(
        ex, local_dict=local_dict, global_dict=global_dict, out=out, ne_args=kwargs, _frame_depth=3
    )
    if out is not None:
        # The user specified an output array
        return lexpr.compute()
    # The user did not specify an output array, so return a NumPy array
    return lexpr[()]


if __name__ == "__main__":
    from time import time

    # Create initial containers
    na1 = np.linspace(0, 10, 10_000_000, dtype=np.float64)
    a1 = blosc2.asarray(na1)
    na2 = np.copy(na1)
    a2 = blosc2.asarray(na2)
    na3 = np.copy(na1)
    a3 = blosc2.asarray(na3)
    na4 = np.copy(na1)
    a4 = blosc2.asarray(na4)
    # Interesting slice
    # sl = None
    sl = slice(0, 10_000)
    # Create a simple lazy expression
    expr = a1 + a2
    print(expr)
    t0 = time()
    nres = na1 + na2
    print(f"Elapsed time (numpy, [:]): {time() - t0:.3f} s")
    t0 = time()
    nres = ne_evaluate("na1 + na2")
    print(f"Elapsed time (numexpr, [:]): {time() - t0:.3f} s")
    nres = nres[sl] if sl is not None else nres
    t0 = time()
    res = expr.compute(item=sl)
    print(f"Elapsed time (evaluate): {time() - t0:.3f} s")
    res = res[sl] if sl is not None else res[:]
    t0 = time()
    res2 = expr[sl]
    print(f"Elapsed time (getitem): {time() - t0:.3f} s")
    np.testing.assert_allclose(res, nres)
    np.testing.assert_allclose(res2, nres)

    # Complex lazy expression
    expr = blosc2.tan(a1) * (blosc2.sin(a2) * blosc2.sin(a2) + blosc2.cos(a3)) + (blosc2.sqrt(a4) * 2)
    # expr = blosc2.sin(a1) + 2 * a1 + 1
    expr += 2
    print(expr)
    t0 = time()
    nres = np.tan(na1) * (np.sin(na2) * np.sin(na2) + np.cos(na3)) + (np.sqrt(na4) * 2) + 2
    # nres = np.sin(na1[:]) + 2 * na1[:] + 1 + 2
    print(f"Elapsed time (numpy, [:]): {time() - t0:.3f} s")
    t0 = time()
    nres = ne_evaluate("tan(na1) * (sin(na2) * sin(na2) + cos(na3)) + (sqrt(na4) * 2) + 2")
    print(f"Elapsed time (numexpr, [:]): {time() - t0:.3f} s")
    nres = nres[sl] if sl is not None else nres
    t0 = time()
    res = expr.compute(sl)
    print(f"Elapsed time (evaluate): {time() - t0:.3f} s")
    res = res[sl] if sl is not None else res[:]
    t0 = time()
    res2 = expr[sl]
    print(f"Elapsed time (getitem): {time() - t0:.3f} s")
    np.testing.assert_allclose(res, nres)
    np.testing.assert_allclose(res2, nres)
    print("Everything is working fine")
