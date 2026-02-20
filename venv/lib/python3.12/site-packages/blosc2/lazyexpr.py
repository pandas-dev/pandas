#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

# Avoid checking the name of type annotations at run time
from __future__ import annotations

import ast
import asyncio
import builtins
import concurrent.futures
import copy
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
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from queue import Empty, Queue
from typing import TYPE_CHECKING, Any

from numpy.exceptions import ComplexWarning

from . import exceptions

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

import ndindex
import numpy as np

import blosc2

if blosc2._HAS_NUMBA:
    import numba
from blosc2 import compute_chunks_blocks
from blosc2.info import InfoReporter

from .proxy import _convert_dtype
from .utils import (
    NUMPY_GE_2_0,
    constructors,
    elementwise_funcs,
    get_chunks_idx,
    get_intersecting_chunks,
    infer_shape,
    linalg_attrs,
    linalg_funcs,
    npvecdot,
    process_key,
    reducers,
)

if not blosc2.IS_WASM:
    import numexpr

global safe_blosc2_globals
safe_blosc2_globals = {}
global safe_numpy_globals
# Use numpy eval when running in WebAssembly
safe_numpy_globals = {"np": np}
# Add all first-level numpy functions
safe_numpy_globals.update(
    {name: getattr(np, name) for name in dir(np) if callable(getattr(np, name)) and not name.startswith("_")}
)

if not NUMPY_GE_2_0:  # handle non-array-api compliance
    safe_numpy_globals["acos"] = np.arccos
    safe_numpy_globals["acosh"] = np.arccosh
    safe_numpy_globals["asin"] = np.arcsin
    safe_numpy_globals["asinh"] = np.arcsinh
    safe_numpy_globals["atan"] = np.arctan
    safe_numpy_globals["atanh"] = np.arctanh
    safe_numpy_globals["atan2"] = np.arctan2
    safe_numpy_globals["permute_dims"] = np.transpose
    safe_numpy_globals["pow"] = np.power
    safe_numpy_globals["bitwise_left_shift"] = np.left_shift
    safe_numpy_globals["bitwise_right_shift"] = np.right_shift
    safe_numpy_globals["bitwise_invert"] = np.bitwise_not
    safe_numpy_globals["concat"] = np.concatenate
    safe_numpy_globals["matrix_transpose"] = np.transpose
    safe_numpy_globals["vecdot"] = npvecdot

# Set this to False if miniexpr should not be tried out
try_miniexpr = True
if blosc2.IS_WASM:
    try_miniexpr = False

_MINIEXPR_WINDOWS_OVERRIDE = os.environ.get("BLOSC2_ENABLE_MINIEXPR_WINDOWS", "").strip().lower()
_MINIEXPR_WINDOWS_OVERRIDE = _MINIEXPR_WINDOWS_OVERRIDE not in ("", "0", "false", "no", "off")


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
        if "out" in kwargs:
            out = kwargs.pop("out")
            out[:] = eval(expression, safe_numpy_globals, local_dict)
            return out
        return eval(expression, safe_numpy_globals, local_dict)
    try:
        return numexpr.evaluate(expression, local_dict=local_dict, **kwargs)
    except ValueError as e:
        raise e  # unsafe expression
    except Exception:  # non_numexpr functions present
        global safe_blosc2_globals
        res = eval(expression, safe_blosc2_globals, local_dict)
        if "out" in kwargs:
            out = kwargs.pop("out")
            out[:] = res  # will handle calc/decomp if res is lazyarray
            return out
        return res[()] if isinstance(res, blosc2.Operand) else res


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
# Gather all callable functions in numpy
numpy_funcs = {
    name
    for name, member in inspect.getmembers(np, callable)
    if not name.startswith("_") and not isinstance(member, np.ufunc)
}
numpy_ufuncs = {name for name, member in inspect.getmembers(np, lambda x: isinstance(x, np.ufunc))}
# Add these functions to the list of available functions
# (will be evaluated via the array interface)
additional_funcs = sorted((numpy_funcs | numpy_ufuncs) - set(blosc2_funcs))
functions = blosc2_funcs + additional_funcs

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
)


def get_expr_globals(expression):
    """Build a dictionary of functions needed for evaluating the expression."""
    _globals = {"np": np, "blosc2": blosc2}
    # Only check for functions that actually appear in the expression
    # This avoids many unnecessary string searches
    for func in functions:
        if func in expression:
            # Try blosc2 first
            if hasattr(blosc2, func):
                _globals[func] = getattr(blosc2, func)
            # Fall back to numpy
            elif hasattr(np, func):
                _globals[func] = getattr(np, func)
            # Function not found in either module
            else:
                raise AttributeError(f"Function {func} not found in blosc2 or numpy")

    return _globals


class ReduceOp(Enum):
    """
    Available reduce operations.
    """

    SUM = np.add
    PROD = np.multiply
    MEAN = np.mean
    STD = np.std
    VAR = np.var
    # Computing a median from partial results is not straightforward because the median
    # is a positional statistic, which means it depends on the relative ordering of all
    # the data points. Unlike statistics such as the sum or mean, you can't compute a median
    # from partial results without knowing the entire dataset, and this is way too expensive
    # for arrays that cannot typically fit in-memory (e.g. disk-based NDArray).
    # MEDIAN = np.median
    MAX = np.maximum
    MIN = np.minimum
    ANY = np.any
    ALL = np.all
    ARGMAX = np.argmax
    ARGMIN = np.argmin


class LazyArrayEnum(Enum):
    """
    Available LazyArrays.
    """

    Expr = 0
    UDF = 1


class LazyArray(ABC, blosc2.Operand):
    @abstractmethod
    def indices(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        """
        Return an :ref:`LazyArray` containing the indices where self is True.

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
            The indices of the :ref:`LazyArray` self that are True.
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

        fp_accuracy: :ref:`blosc2.FPAccuracy`, optional
            Specifies the floating-point accuracy to be used during computation.
            By default, :ref:`blosc2.FPAccuracy.DEFAULT` is used.

        kwargs: Any, optional
            Keyword arguments that are supported by the :func:`empty` constructor.
            These arguments will be set in the resulting :ref:`NDArray`.
            Additionally, the following special kwargs are supported:

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
        * All the operands of the LazyArray must be Python scalars, or :ref:`blosc2.Array` objects.
        * If an operand is a :ref:`Proxy`, keep in mind that Python-Blosc2 will only be able to reopen it as such
          if its source is a :ref:`SChunk`, :ref:`NDArray` or a :ref:`C2Array` (see :func:`blosc2.open` notes
          section for more info).
        * This is currently only supported for :ref:`LazyExpr`.

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
        >>> disk_expr = blosc2.open('lazy_array.b2nd')
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
        if not isinstance(obj, (np.ndarray, blosc2.Operand)) and not np.isscalar(obj):
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


def check_smaller_shape(value_shape, shape, slice_shape, slice_):
    """Check whether the shape of the value is smaller than the shape of the array.

    This follows the NumPy broadcasting rules.
    """
    # slice_shape must be as long as shape
    if len(slice_shape) != len(slice_):
        raise ValueError("slice_shape must be as long as slice_")
    no_nones_shape = tuple(sh for sh, s in zip(slice_shape, slice_, strict=True) if s is not None)
    no_nones_slice = tuple(s for sh, s in zip(slice_shape, slice_, strict=True) if s is not None)
    is_smaller_shape = any(
        s > (1 if i >= len(value_shape) else value_shape[i]) for i, s in enumerate(no_nones_shape)
    )
    slice_past_bounds = any(
        s.stop > (1 if i >= len(value_shape) else value_shape[i]) for i, s in enumerate(no_nones_slice)
    )
    return len(value_shape) < len(shape) or is_smaller_shape or slice_past_bounds


def _compute_smaller_slice(larger_shape, smaller_shape, larger_slice):
    smaller_slice = []
    diff_dims = len(larger_shape) - len(smaller_shape)

    for i in range(len(larger_shape)):
        if i < diff_dims:
            # For leading dimensions of the larger array that the smaller array doesn't have,
            # we don't add anything to the smaller slice
            pass
        else:
            # For dimensions that both arrays have, the slice for the smaller array should be
            # the same as the larger array unless the smaller array's size along that dimension
            # is 1, in which case we use None to indicate the full slice
            if smaller_shape[i - diff_dims] != 1:
                smaller_slice.append(larger_slice[i])
            else:
                smaller_slice.append(slice(0, larger_shape[i]))

    return tuple(smaller_slice)


# A more compact version of the function above, albeit less readable
def compute_smaller_slice(larger_shape, smaller_shape, larger_slice):
    """
    Returns the slice of the smaller array that corresponds to the slice of the larger array.
    """
    j_small = len(smaller_shape) - 1
    j_large = len(larger_shape) - 1
    smaller_shape_nones = []
    larger_shape_nones = []
    for s in reversed(larger_slice):
        if s is None:
            smaller_shape_nones.append(1)
            larger_shape_nones.append(1)
        else:
            if j_small >= 0:
                smaller_shape_nones.append(smaller_shape[j_small])
                j_small -= 1
            if j_large >= 0:
                larger_shape_nones.append(larger_shape[j_large])
                j_large -= 1
    smaller_shape_nones.reverse()
    larger_shape_nones.reverse()
    diff_dims = len(larger_shape_nones) - len(smaller_shape_nones)
    return tuple(
        None
        if larger_slice[i] is None
        else (
            larger_slice[i] if smaller_shape_nones[i - diff_dims] != 1 else slice(0, larger_shape_nones[i])
        )
        for i in range(diff_dims, len(larger_shape_nones))
    )


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

    inputs = [input for input in inputs.values() if hasattr(input, "shape") and input is not np]
    # This will raise an exception if the input shapes are not compatible
    shape = compute_broadcast_shape(inputs)

    if not all(np.array_equal(shape, input.shape) for input in inputs):
        # If inputs have different shapes, we cannot take the fast path
        return shape, None, None, False

    # More checks specific of NDArray inputs
    # NDInputs are either non-SimpleProxy with chunks or are SimpleProxy with src having chunks
    NDinputs = [
        input
        for input in inputs
        if (hasattr(input, "chunks") and not isinstance(input, blosc2.SimpleProxy))
        or (isinstance(input, blosc2.SimpleProxy) and hasattr(input.src, "chunks"))
    ]
    if not NDinputs:
        # All inputs are NumPy arrays, so we cannot take the fast path
        if inputs and hasattr(inputs[0], "shape"):
            shape = inputs[0].shape
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


async def async_read_chunks(arrs, info, queue):
    loop = asyncio.get_event_loop()
    nchunks = arrs[0].schunk.nchunks

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for nchunk in range(nchunks):
            futures = [
                (index, loop.run_in_executor(executor, get_chunk, arr, info, nchunk))
                for index, arr in enumerate(arrs)
            ]
            chunks = await asyncio.gather(*(future for index, future in futures), return_exceptions=True)
            chunks_sorted = []
            for chunk in chunks:
                if isinstance(chunk, Exception):
                    # Handle the exception (e.g., log it, raise a custom exception, etc.)
                    print(f"Exception occurred: {chunk}")
                    raise chunk
                chunks_sorted.append(chunk)
            queue.put((nchunk, chunks_sorted))  # use non-async queue.put()

    queue.put(None)  # signal the end of the chunks


def async_read_chunks_thread(arrs, info, queue):
    asyncio.run(async_read_chunks(arrs, info, queue))


def sync_read_chunks(arrs, info):
    queue_size = 2  # maximum number of chunks in the queue
    queue = Queue(maxsize=queue_size)

    # Start the async file reading in a separate thread
    thread = threading.Thread(target=async_read_chunks_thread, args=(arrs, info, queue))
    thread.start()

    # Read the chunks synchronously from the queue
    while True:
        try:
            chunks = queue.get(timeout=1)  # Wait for the next chunk
            if chunks is None:  # End of chunks
                break
            yield chunks
        except Empty:
            continue


def read_nchunk(arrs, info):
    for _, chunks in sync_read_chunks(arrs, info):
        yield chunks


iter_chunks = None


def fill_chunk_operands(
    operands, slice_, chunks_, full_chunk, aligned, nchunk, iter_disk, chunk_operands, reduc=False
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
        if nchunk == 0:
            # Initialize the iterator for reading the chunks
            # Take any operand (all should have the same shape and chunks)
            key, arr = next(iter(operands.items()))
            chunks_idx, _ = get_chunks_idx(arr.shape, arr.chunks)
            info = (reduc, aligned[key], low_mem, chunks_idx)
            iter_chunks = read_nchunk(list(operands.values()), info)
        # Run the asynchronous file reading function from a synchronous context
        chunks = next(iter_chunks)

        for i, (key, value) in enumerate(operands.items()):
            # Chunks are already decompressed, so we can use them directly
            if not low_mem:
                chunk_operands[key] = chunks[i]
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

    # Disable miniexpr for UDFs (callable expressions)
    if callable(expression):
        use_miniexpr = False

    out = kwargs.pop("_output", None)
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    fp_accuracy = kwargs.pop("fp_accuracy", blosc2.FPAccuracy.DEFAULT)
    dtype = kwargs.pop("dtype", None)
    where: dict | None = kwargs.pop("_where_args", None)
    if where is not None:
        # miniexpr does not support where(); use the regular path.
        use_miniexpr = False
    if isinstance(out, blosc2.NDArray):
        # If 'out' has been passed, and is a NDArray, use it as the base array
        basearr = out
    elif isinstance(out, np.ndarray):
        # If 'out' is a NumPy array, create a NDArray with the same shape and dtype
        basearr = blosc2.empty(out.shape, dtype=out.dtype, **kwargs)
    else:
        # Otherwise, find the operand with the 'chunks' attribute and the longest shape
        operands_with_chunks = [o for o in operands.values() if hasattr(o, "chunks")]
        basearr = max(operands_with_chunks, key=lambda x: len(x.shape))

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

    # Check whether we can use miniexpr
    if use_miniexpr:
        # Require aligned NDArray operands with identical chunk/block grid.
        same_shape = all(hasattr(op, "shape") and op.shape == shape for op in operands.values())
        same_chunks = all(hasattr(op, "chunks") and op.chunks == chunks for op in operands.values())
        same_blocks = all(hasattr(op, "blocks") and op.blocks == blocks for op in operands.values())
        if not (same_shape and same_chunks and same_blocks):
            use_miniexpr = False
        if not (all_ndarray and out is None):
            use_miniexpr = False
        has_complex = any(
            isinstance(op, blosc2.NDArray) and blosc2.isdtype(op.dtype, "complex floating")
            for op in operands.values()
        )
        if isinstance(expression, str) and has_complex:
            if sys.platform == "win32":
                # On Windows, miniexpr has issues with complex numbers
                use_miniexpr = False
            if any(tok in expression for tok in ("!=", "==", "<=", ">=", "<", ">")):
                use_miniexpr = False
        if sys.platform == "win32" and use_miniexpr and not _MINIEXPR_WINDOWS_OVERRIDE:
            # Work around Windows miniexpr issues for integer outputs and dtype conversions.
            if blosc2.isdtype(dtype, "integral"):
                use_miniexpr = False
            else:
                dtype_mismatch = any(
                    isinstance(op, blosc2.NDArray) and op.dtype != dtype for op in operands.values()
                )
                if dtype_mismatch:
                    use_miniexpr = False

    if use_miniexpr:
        cparams = kwargs.pop("cparams", blosc2.CParams())
        # All values will be overwritten, so we can use an uninitialized array
        res_eval = blosc2.uninit(shape, dtype, chunks=chunks, blocks=blocks, cparams=cparams, **kwargs)
        try:
            res_eval._set_pref_expr(expression, operands, fp_accuracy=fp_accuracy)
            # print("expr->miniexpr:", expression, fp_accuracy)
            # Data to compress is fetched from operands, so it can be uninitialized here
            data = np.empty(res_eval.schunk.chunksize, dtype=np.uint8)
            # Exercise prefilter for each chunk
            for nchunk in range(res_eval.schunk.nchunks):
                res_eval.schunk.update_data(nchunk, data, copy=False)
        except Exception:
            use_miniexpr = False
        finally:
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
    _indices = kwargs.pop("_indices", False)
    if _indices and (not where or len(where) != 1):
        raise NotImplementedError("Indices can only be used with one where condition")
    _order = kwargs.pop("_order", None)
    if _order is not None and not isinstance(_order, list):
        # Always use a list for _order
        _order = [_order]

    dtype = kwargs.pop("dtype", None)
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
        if shape_ != out.shape:
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
        dtype_ = np.dtype([(f, dtype_[f]) for f in _order])

    # Iterate over the operands and get the chunks
    chunk_operands = {}
    # Check which chunks intersect with _slice (handles zero chunks internally)
    intersecting_chunks = get_intersecting_chunks(
        _slice, shape, chunks
    )  # if _slice is (), returns all chunks

    for nchunk, chunk_slice in enumerate(intersecting_chunks):
        # get intersection of chunk and target
        if _slice != ():
            cslice = step_handler(chunk_slice.raw, _slice)
        else:
            cslice = chunk_slice.raw

        cslice_shape = tuple(s.stop - s.start for s in cslice)
        len_chunk = math.prod(cslice_shape)
        # get local index of part of out that is to be updated
        cslice_subidx = (
            ndindex.ndindex(cslice).as_subindex(_slice).raw
        )  # in the case _slice=(), just gives cslice

        # Get the starts and stops for the slice
        starts = [s.start if s.start is not None else 0 for s in cslice]
        stops = [s.stop if s.stop is not None else sh for s, sh in zip(cslice, cslice_shape, strict=True)]
        unit_steps = np.all([s.step == 1 for s in cslice])

        # Get the slice of each operand
        for key, value in operands.items():
            if np.isscalar(value):
                chunk_operands[key] = value
                continue
            if value.shape == ():
                chunk_operands[key] = value[()]
                continue
            if check_smaller_shape(value.shape, shape, cslice_shape, cslice):
                # We need to fetch the part of the value that broadcasts with the operand
                smaller_slice = compute_smaller_slice(shape, value.shape, cslice)
                chunk_operands[key] = value[smaller_slice]
                continue
            # If key is in operands, we can reuse the buffer
            if (
                key in chunk_operands
                and cslice_shape == chunk_operands[key].shape
                and isinstance(value, blosc2.NDArray)
                and unit_steps
            ):
                value.get_slice_numpy(chunk_operands[key], (starts, stops))
                continue

            chunk_operands[key] = value[cslice]

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
            result = np.empty(cslice_shape, dtype=out.dtype)  # raises error if out is None
            # cslice should be equal to cslice_subidx
            # Call the udf directly and use result as the output array
            offset = tuple(s.start for s in cslice)
            expression(tuple(chunk_operands.values()), result, offset=offset)
            out[cslice_subidx] = result
            continue

        if where is None:
            result = ne_evaluate(expression, chunk_operands, **ne_args)
        else:
            # Apply the where condition (in result)
            if len(where) == 2:
                # x = chunk_operands["_where_x"]
                # y = chunk_operands["_where_y"]
                # result = np.where(result, x, y)
                # numexpr is a bit faster than np.where, and we can fuse operations in this case
                new_expr = f"where({expression}, _where_x, _where_y)"
                result = ne_evaluate(new_expr, chunk_operands, **ne_args)
            elif len(where) == 1:
                result = ne_evaluate(expression, chunk_operands, **ne_args)
                if _indices or _order:
                    # Return indices only makes sense when the where condition is a tuple with one element
                    # and result is a boolean array
                    x = chunk_operands["_where_x"]
                    if len(x.shape) > 1:
                        raise ValueError("indices() and sort() only support 1D arrays")
                    if result.dtype != np.bool_:
                        raise ValueError("indices() and sort() only support bool conditions")
                    indices = np.arange(leninputs, leninputs + len_chunk, dtype=np.int64).reshape(
                        cslice_shape
                    )
                    if _order:
                        # We need to cumulate all the fields in _order, as well as indices
                        chunk_indices = indices[result]
                        result = x[_order][result]
                    else:
                        result = indices[result]
                    leninputs += len_chunk
                else:
                    x = chunk_operands["_where_x"]
                    result = x[result]
            else:
                raise ValueError("The where condition must be a tuple with one or two elements")
        # Enforce contiguity of result (necessary to fill the out array)
        # but avoid copy if already contiguous
        result = np.require(result, requirements="C")

        if where is None or len(where) == 2:
            if behaved and result.shape == out.chunks and result.dtype == out.dtype:
                # Fast path
                # TODO: Check this only works when slice is ()
                out.schunk.update_data(nchunk, result, copy=False)
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
    _slice_bcast = tuple(slice(i, i + 1) if isinstance(i, int) else i for i in _slice.raw)
    slice_shape = ndindex.ndindex(_slice_bcast).newshape(shape)  # includes dummy dimensions
    _slice = _slice.raw

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
        offset = tuple(0 if s is None else s.start for s in _slice_bcast)  # offset for the udf
        result = np.empty(slice_shape, dtype=dtype)
        expression(tuple(slice_operands.values()), result, offset=offset)
    else:
        if where is None:
            result = ne_evaluate(expression, slice_operands, **ne_args)
        else:
            # Apply the where condition (in result)
            new_expr = f"where({expression}, _where_x, _where_y)"
            result = ne_evaluate(new_expr, slice_operands, **ne_args)

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
    if operation in {ReduceOp.SUM, ReduceOp.PROD}:
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

    out = kwargs.pop("_output", None)
    res_out_ = None  # temporary required to store max/min for argmax/argmin
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    fp_accuracy = kwargs.pop("fp_accuracy", blosc2.FPAccuracy.DEFAULT)
    where: dict | None = kwargs.pop("_where_args", None)
    reduce_op = reduce_args.pop("op")
    reduce_op_str = reduce_args.pop("op_str", None)
    axis = reduce_args["axis"]
    keepdims = reduce_args["keepdims"]
    dtype = reduce_args.get("dtype", None)
    if dtype is None:
        dtype = kwargs.pop("dtype", None)
        dtype = infer_reduction_dtype(dtype, reduce_op)
    else:
        del kwargs["dtype"]

    # Compute the shape and chunks of the output array, including broadcasting
    shape = compute_broadcast_shape(operands.values())

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
    if not (fast_path and all_ndarray and reduced_shape == ()):
        use_miniexpr = False

    # Some reductions are not supported yet in miniexpr
    if reduce_op in (ReduceOp.ARGMAX, ReduceOp.ARGMIN):
        use_miniexpr = False

    # Check whether we can use miniexpr
    if use_miniexpr and isinstance(expression, str):
        has_complex = any(
            isinstance(op, blosc2.NDArray) and blosc2.isdtype(op.dtype, "complex floating")
            for op in operands.values()
        )
        if has_complex and sys.platform == "win32":
            # On Windows, miniexpr has issues with complex numbers
            use_miniexpr = False
        if sys.platform == "win32" and use_miniexpr and not _MINIEXPR_WINDOWS_OVERRIDE:
            if blosc2.isdtype(dtype, "integral"):
                use_miniexpr = False
            else:
                dtype_mismatch = any(
                    isinstance(op, blosc2.NDArray) and op.dtype != dtype for op in operands.values()
                )
                if dtype_mismatch:
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
        if reduce_op == ReduceOp.SUM or reduce_op == ReduceOp.ANY:
            aux_reduc = np.zeros(nblocks, dtype=dtype)
        elif reduce_op == ReduceOp.PROD or reduce_op == ReduceOp.ALL:
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
        try:
            if where is not None:
                expression_miniexpr = f"{reduce_op_str}(where({expression}, _where_x, _where_y))"
            else:
                expression_miniexpr = f"{reduce_op_str}({expression})"
            res_eval._set_pref_expr(expression_miniexpr, operands, fp_accuracy, aux_reduc)
            # print("expr->miniexpr:", expression, reduce_op, fp_accuracy)
            # Data won't even try to be compressed, so buffers can be unitialized and reused
            data = np.empty(res_eval.schunk.chunksize, dtype=np.uint8)
            chunk_data = np.empty(res_eval.schunk.chunksize + blosc2.MAX_OVERHEAD, dtype=np.uint8)
            # Exercise prefilter for each chunk
            for nchunk in range(res_eval.schunk.nchunks):
                res_eval.schunk._prefilter_data(nchunk, data, chunk_data)
        except Exception:
            use_miniexpr = False
        finally:
            res_eval.schunk.remove_prefilter("miniexpr")
            global iter_chunks
            # Ensure any background reading thread is closed
            iter_chunks = None

        if not use_miniexpr:
            # If miniexpr failed, fallback to regular evaluation
            # (continue to the manual chunked evaluation below)
            pass
        else:
            if reduce_op == ReduceOp.ANY:
                result = np.any(aux_reduc, **reduce_args)
            elif reduce_op == ReduceOp.ALL:
                result = np.all(aux_reduc, **reduce_args)
            else:
                result = reduce_op.value.reduce(aux_reduc, **reduce_args)
            return result

    # Iterate over the operands and get the chunks
    chunk_operands = {}
    # Check which chunks intersect with _slice
    # if chunks has 0 we loop once but fast path is false as gives error (schunk has no chunks)
    intersecting_chunks = get_intersecting_chunks(_slice, shape, chunks)
    out_init = False

    for nchunk, chunk_slice in enumerate(intersecting_chunks):
        cslice = chunk_slice.raw
        offset = tuple(s.start for s in cslice)  # offset for the udf
        # Check whether current cslice intersects with _slice
        if cslice != () and _slice != ():
            # get intersection of chunk and target
            cslice = step_handler(cslice, _slice)
        chunks_ = tuple(s.stop - s.start for s in cslice)
        unit_steps = np.all([s.step == 1 for s in cslice])
        # Starts for slice
        starts = [s.start if s.start is not None else 0 for s in cslice]
        if _slice == () and fast_path and unit_steps:
            # Fast path
            full_chunk = chunks_ == chunks
            fill_chunk_operands(
                operands, cslice, chunks_, full_chunk, aligned, nchunk, iter_disk, chunk_operands, reduc=True
            )
        else:
            # Get the stops for the slice
            stops = [s.stop if s.stop is not None else sh for s, sh in zip(cslice, chunks_, strict=True)]
            # Get the slice of each operand
            for key, value in operands.items():
                if np.isscalar(value):
                    chunk_operands[key] = value
                    continue
                if value.shape == ():
                    chunk_operands[key] = value[()]
                    continue
                if check_smaller_shape(value.shape, shape, chunks_, cslice):
                    # We need to fetch the part of the value that broadcasts with the operand
                    smaller_slice = compute_smaller_slice(operand.shape, value.shape, cslice)
                    chunk_operands[key] = value[smaller_slice]
                    continue
                # If key is in operands, we can reuse the buffer
                if (
                    key in chunk_operands
                    and chunks_ == chunk_operands[key].shape
                    and isinstance(value, blosc2.NDArray)
                    and unit_steps
                ):
                    value.get_slice_numpy(chunk_operands[key], (starts, stops))
                    continue
                chunk_operands[key] = value[cslice]

        # get local index of part of out that is to be updated
        cslice_subidx = ndindex.ndindex(cslice).as_subindex(_slice).raw  # if _slice is (), just gives cslice
        if keepdims:
            reduced_slice = tuple(slice(None) if i in axis else sl for i, sl in enumerate(cslice_subidx))
        else:
            reduced_slice = tuple(sl for i, sl in enumerate(cslice_subidx) if i not in axis)

        # Evaluate and reduce the expression using chunks of operands

        if callable(expression):
            # TODO: Implement the reductions for UDFs (and test them)
            result = np.empty(chunks_, dtype=out.dtype)
            expression(tuple(chunk_operands.values()), result, offset=offset)
            # Reduce the result
            result = reduce_op.value.reduce(result, **reduce_args)
            # Update the output array with the result
            out[reduced_slice] = reduce_op.value(out[reduced_slice], result)
            continue

        if where is None:
            if expression == "o0" or expression == "(o0)":
                # We don't have an actual expression, so avoid a copy except to make contiguous
                result = np.require(chunk_operands["o0"], requirements="C")
            else:
                result = ne_evaluate(expression, chunk_operands, **ne_args)
        else:
            # Apply the where condition (in result)
            if len(where) == 2:
                new_expr = f"where({expression}, _where_x, _where_y)"
                result = ne_evaluate(new_expr, chunk_operands, **ne_args)
            elif len(where) == 1:
                result = ne_evaluate(expression, chunk_operands, **ne_args)
                x = chunk_operands["_where_x"]
                result = x[result]
            else:
                raise NotImplementedError(
                    "A where condition with no params in combination with reductions is not supported yet"
                )

        # Reduce the result
        if result.shape == ():
            if reduce_op == ReduceOp.SUM and result[()] == 0:
                # Avoid a reduction when result is a zero scalar. Faster for sparse data.
                continue
            # Note that chunks_ refers to slice of operand chunks, not reduced_slice
            result = np.full(chunks_, result[()])
        if reduce_op == ReduceOp.ANY:
            result = np.any(result, **reduce_args)
        elif reduce_op == ReduceOp.ALL:
            result = np.all(result, **reduce_args)
        elif reduce_op == ReduceOp.ARGMAX or reduce_op == ReduceOp.ARGMIN:
            # offset for start of slice
            slice_ref = (
                starts
                if _slice == ()
                else [
                    (s - sl.start - np.sign(sl.step)) // sl.step + 1
                    for s, sl in zip(starts, _slice, strict=True)
                ]
            )
            result_idx = (
                np.argmin(result, **reduce_args)
                if reduce_op == ReduceOp.ARGMIN
                else np.argmax(result, **reduce_args)
            )
            if reduce_args["axis"] is None:  # indexing into flattened array
                result = result[np.unravel_index(result_idx, shape=result.shape)]
                idx_within_cslice = np.unravel_index(result_idx, shape=chunks_)
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
            out_, res_out_ = convert_none_out(result.dtype, reduce_op, reduced_shape)
            if out is not None:
                out[:] = out_
                del out_
            else:
                out = out_
            out_init = True

        # Update the output array with the result
        if reduce_op == ReduceOp.ANY:
            out[reduced_slice] += result
        elif reduce_op == ReduceOp.ALL:
            out[reduced_slice] *= result
        elif res_out_ is not None:  # i.e. ReduceOp.ARGMAX or ReduceOp.ARGMIN
            # need lowest index for which optimum attained
            cond = (res_out_[reduced_slice] == result) & (result_idx < out[reduced_slice])
            if reduce_op == ReduceOp.ARGMAX:
                cond |= res_out_[reduced_slice] < result
            else:  # ARGMIN
                cond |= res_out_[reduced_slice] > result
            if reduced_slice == ():
                out = np.where(cond, result_idx, out[reduced_slice])
                res_out_ = np.where(cond, result, res_out_[reduced_slice])
            else:
                out[reduced_slice] = np.where(cond, result_idx, out[reduced_slice])
                res_out_[reduced_slice] = np.where(cond, result, res_out_[reduced_slice])
        else:
            if reduced_slice == ():
                out = reduce_op.value(out, result)
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
        out, _ = convert_none_out(dtype, reduce_op, reduced_shape)

    final_mask = tuple(np.where(mask_slice)[0])
    if np.any(mask_slice):  # remove dummy dims
        out = np.squeeze(out, axis=final_mask)
    # Check if the output array needs to be converted into a blosc2.NDArray
    if kwargs != {} and not np.isscalar(out):
        out = blosc2.asarray(out, **kwargs)
    return out


def convert_none_out(dtype, reduce_op, reduced_shape):
    out = None
    # out will be a proper numpy.ndarray
    if reduce_op == ReduceOp.SUM:
        out = np.zeros(reduced_shape, dtype=dtype)
    elif reduce_op == ReduceOp.PROD:
        out = np.ones(reduced_shape, dtype=dtype)
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
    elif reduce_op == ReduceOp.ARGMIN:
        if np.issubdtype(dtype, np.integer):
            res_out_ = np.iinfo(dtype).max * np.ones(reduced_shape, dtype=dtype)
        else:
            res_out_ = np.inf * np.ones(reduced_shape, dtype=dtype)
        out = (np.zeros(reduced_shape, dtype=blosc2.DEFAULT_INDEX), res_out_)
    elif reduce_op == ReduceOp.ARGMAX:
        if np.issubdtype(dtype, np.integer):
            res_out_ = np.iinfo(dtype).min * np.ones(reduced_shape, dtype=dtype)
        else:
            res_out_ = -np.inf * np.ones(reduced_shape, dtype=dtype)
        out = (np.zeros(reduced_shape, dtype=blosc2.DEFAULT_INDEX), res_out_)
    return out if isinstance(out, tuple) else (out, None)


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

        where: dict | None = kwargs.get("_where_args")
        if where:
            # Make the where arguments part of the operands
            operands = {**operands, **where}

        reduce_args = kwargs.pop("_reduce_args", {})
        _, _, _, fast_path = validate_inputs(operands, out, reduce=reduce_args != {})

        # Activate last read cache for NDField instances
        for op in operands:
            if isinstance(operands[op], blosc2.NDField):
                operands[op].ndarr.keep_last_read = True

        if reduce_args:
            # Eval and reduce the expression in a single step
            return reduce_slices(expression, operands, reduce_args=reduce_args, _slice=item, **kwargs)

        if not is_full_slice(item.raw) or (where is not None and len(where) < 2):
            # The fast path is possible under a few conditions
            if getitem and (where is None or len(where) == 2):
                # Compute the size of operands for the fast path
                unit_steps = np.all([s.step == 1 for s in item.raw if isinstance(s, slice)])
                # shape of slice, if non-unit steps have to decompress full array into memory
                shape_operands = item.newshape(shape) if unit_steps else shape
                _dtype = kwargs.get("dtype", np.float64)
                size_operands = math.prod(shape_operands) * len(operands) * _dtype.itemsize
                # Only take the fast path if the size of operands is relatively small
                if size_operands < blosc2.MAX_FAST_PATH_SIZE:
                    return slices_eval_getitem(expression, operands, _slice=item, shape=shape, **kwargs)
            return slices_eval(expression, operands, getitem=getitem, _slice=item, shape=shape, **kwargs)

        fast_path = is_full_slice(item.raw) and fast_path
        if fast_path:  # necessarily item is ()
            if getitem:
                # When using getitem, taking the fast path is always possible
                return fast_eval(expression, operands, getitem=True, **kwargs)
            elif (kwargs.get("chunks") is None and kwargs.get("blocks") is None) and (
                out is None or isinstance(out, blosc2.NDArray)
            ):
                # If not, the conditions to use the fast path are a bit more restrictive
                # e.g. the user cannot specify chunks or blocks, or an output that is not
                # a blosc2.NDArray
                return fast_eval(expression, operands, getitem=False, **kwargs)

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
    if op == "contains":
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
    *arrays_and_dtypes: blosc2.NDArray | int | float | complex | bool | blosc2.dtype,
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
        value
        if (np.isscalar(value) or not hasattr(value, "dtype"))
        else np.array([0], dtype=_convert_dtype(value.dtype))
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
        dtype_ = check_dtype(op, value1, value2)  # perform some checks
        # Check that operands are proper Operands, LazyArray or scalars; if not, convert to NDArray objects
        value1 = (
            blosc2.SimpleProxy(value1)
            if not (isinstance(value1, (blosc2.Operand, np.ndarray)) or np.isscalar(value1))
            else value1
        )
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
                    value1 = ne_evaluate(f"{op}({value1})")
                    op = None
                self.operands = {"o0": value1}
                self.expression = "o0" if op is None else f"{op}(o0)"
            return
        value2 = (
            blosc2.SimpleProxy(value2)
            if not (isinstance(value2, (blosc2.Operand, np.ndarray)) or np.isscalar(value2))
            else value2
        )
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
                self.operands = {"o0": ne_evaluate(f"{op}({value1}, {value2})")}  # eager evaluation
            elif np.isscalar(value2):
                self.operands = {"o0": value1}
                self.expression = f"{op}(o0, {value2})"
            elif np.isscalar(value1):
                self.operands = {"o0": value2}
                self.expression = f"{op}({value1}, o0)"
            else:
                self.operands = {"o0": value1, "o1": value2}
                self.expression = f"{op}(o0, o1)"
            return

        self._dtype = dtype_
        if np.isscalar(value1) and np.isscalar(value2):
            self.expression = "o0"
            self.operands = {"o0": ne_evaluate(f"({value1} {op} {value2})")}  # eager evaluation
        elif np.isscalar(value2):
            self.operands = {"o0": value1}
            self.expression = f"(o0 {op} {value2})"
        elif hasattr(value2, "shape") and value2.shape == ():
            self.operands = {"o0": value1}
            self.expression = f"(o0 {op} {value2[()]})"
        elif np.isscalar(value1):
            self.operands = {"o0": value2}
            self.expression = f"({value1} {op} o0)"
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

    def update_expr(self, new_op):  # noqa: C901
        prev_flag = blosc2._disable_overloaded_equal
        # We use a lot of the original NDArray.__eq__ as 'is', so deactivate the overloaded one
        blosc2._disable_overloaded_equal = True
        # One of the two operands are LazyExpr instances
        try:
            value1, op, value2 = new_op
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
                    v2 = value2
                elif hasattr(value2, "shape") and value2.shape == ():
                    v2 = value2[()]
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
                    v1 = value1
                elif hasattr(value1, "shape") and value1.shape == ():
                    v1 = value1[()]
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
        if any(constructor in self.expression for constructor in constructors):
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

    def sum(
        self,
        axis=None,
        dtype=None,
        keepdims=False,
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
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
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
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
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        item = kwargs.pop("item", ())
        total_sum = self.sum(
            axis=axis,
            dtype=dtype,
            keepdims=keepdims,
            item=item,
            fp_accuracy=fp_accuracy,
        )
        num_elements = self.get_num_elements(axis, item)
        if num_elements == 0:
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
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        item = kwargs.pop("item", ())
        if item == ():  # fast path
            mean_value = self.mean(axis=axis, dtype=dtype, keepdims=True, fp_accuracy=fp_accuracy)
            expr = (self - mean_value) ** 2
        else:
            mean_value = self.mean(axis=axis, dtype=dtype, keepdims=True, item=item, fp_accuracy=fp_accuracy)
            # TODO: Not optimal because we load the whole slice in memory. Would have to write
            #  a bespoke std function that executed within slice_eval to avoid this probably.
            expr = (self.slice(item) - mean_value) ** 2
        out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, fp_accuracy=fp_accuracy)
        if ddof != 0:
            num_elements = self.get_num_elements(axis, item)
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
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
        item = kwargs.pop("item", ())
        if item == ():  # fast path
            mean_value = self.mean(axis=axis, dtype=dtype, keepdims=True, fp_accuracy=fp_accuracy)
            expr = (self - mean_value) ** 2
        else:
            mean_value = self.mean(axis=axis, dtype=dtype, keepdims=True, item=item, fp_accuracy=fp_accuracy)
            # TODO: Not optimal because we load the whole slice in memory. Would have to write
            #  a bespoke var function that executed within slice_eval to avoid this probably.
            expr = (self.slice(item) - mean_value) ** 2
        out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, fp_accuracy=fp_accuracy)
        if ddof != 0:
            num_elements = self.get_num_elements(axis, item)
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
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
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
        fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT,
        **kwargs,
    ):
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

        # Find the index of the first parenthesis after the constructor
        idx = expression.find(f"{constructor}")
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

    def _compute_expr(self, item, kwargs):  # noqa : C901
        # ne_evaluate will need safe_blosc2_globals for some functions (e.g. clip, logaddexp)
        # that are implemented in python-blosc2 not in numexpr
        global safe_blosc2_globals
        if len(safe_blosc2_globals) == 0:
            # First eval call, fill blosc2_safe_globals for ne_evaluate
            safe_blosc2_globals = {"blosc2": blosc2}
            # Add all first-level blosc2 functions
            safe_blosc2_globals.update(
                {
                    name: getattr(blosc2, name)
                    for name in dir(blosc2)
                    if callable(getattr(blosc2, name)) and not name.startswith("_")
                }
            )

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
                    # Correct shape to adjust to NumPy convention
                    arr.shape = tuple(arr.shape[i] for i in range(len(mask)) if not mask[i])
                return arr

            return chunked_eval(lazy_expr.expression, lazy_expr.operands, item, **kwargs)

        if any(constructor in self.expression for constructor in constructors):
            expression = self.expression
            newexpr = expression
            newops = self.operands.copy()
            # We have constructors in the expression (probably coming from a string lazyexpr)
            # Let's replace the constructors with the actual NDArray objects
            for constructor in constructors:
                if constructor not in newexpr:
                    continue
                while constructor in newexpr:
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

        return chunked_eval(self.expression, self.operands, item, **kwargs)

    # TODO: indices and sort are repeated in LazyUDF; refactor
    def indices(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        if self.dtype.fields is None:
            raise NotImplementedError("indices() can only be used with structured arrays")
        if not hasattr(self, "_where_args") or len(self._where_args) != 1:
            raise ValueError("indices() can only be used with conditions")
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
        self, item=(), fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT, **kwargs
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

        expression = self.expression_tosave if hasattr(self, "expression_tosave") else self.expression
        operands_ = self.operands_tosave if hasattr(self, "operands_tosave") else self.operands
        # Validate expression
        validate_expr(expression)

        meta = kwargs.get("meta", {})
        meta["LazyArray"] = LazyArrayEnum.Expr.value
        kwargs["urlpath"] = urlpath
        kwargs["meta"] = meta
        kwargs["mode"] = "w"  # always overwrite the file in urlpath

        # Create an empty array; useful for providing the shape and dtype of the outcome
        array = blosc2.empty(shape=self.shape, dtype=self.dtype, **kwargs)

        # Save the expression and operands in the metadata
        operands = {}
        for key, value in operands_.items():
            if isinstance(value, blosc2.C2Array):
                operands[key] = {
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
                )
            if value.schunk.urlpath is None:
                raise ValueError("To save a LazyArray, all operands must be stored on disk/network")
            operands[key] = value.schunk.urlpath
        array.schunk.vlmeta["_LazyArray"] = {
            "expression": expression,
            "UDF": None,
            "operands": operands,
        }

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
                if not (isinstance(val, (blosc2.Operand, np.ndarray)) or np.isscalar(val)):
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
    def __init__(self, func, inputs, dtype, shape=None, chunked_eval=True, **kwargs):
        # After this, all the inputs should be np.ndarray or NDArray objects
        self.inputs = convert_inputs(inputs)
        self.chunked_eval = True  # chunked_eval
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
        self._dtype = dtype
        self.func = func

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

        # TODO: enable parallelism using python 3.14t
        # self.res_getitem = blosc2.empty(self._shape, self._dtype, **kwargs_getitem)
        # # Register a postfilter for getitem
        # if 0 not in self._shape:
        #     self.res_getitem._set_postf_udf(self.func, id(self.inputs))

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

    # TODO: indices and sort are repeated in LazyExpr; refactor
    def indices(self, order: str | list[str] | None = None) -> blosc2.LazyArray:
        if self.dtype.fields is None:
            raise NotImplementedError("indices() can only be used with structured arrays")
        if not hasattr(self, "_where_args") or len(self._where_args) != 1:
            raise ValueError("indices() can only be used with conditions")
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

    def compute(self, item=(), fp_accuracy: blosc2.FPAccuracy = blosc2.FPAccuracy.DEFAULT, **kwargs):
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

        if self.chunked_eval:
            # aux_kwargs includes self.shape and self.dtype
            return chunked_eval(self.func, self.inputs_dict, item, _getitem=False, **aux_kwargs)

        # TODO: Implement multithreading
        # # Cannot use multithreading when applying a prefilter, save nthreads to set them
        # # after the evaluation
        # cparams = aux_kwargs.get("cparams", {})
        # if isinstance(cparams, dict):
        #     self._cnthreads = cparams.get("nthreads", blosc2.cparams_dflts["nthreads"])
        #     cparams["nthreads"] = 1
        # else:
        #     raise ValueError("cparams should be a dictionary")
        # aux_kwargs["cparams"] = cparams

        # res_eval = blosc2.empty(self.shape, self.dtype, **aux_kwargs)
        # # Register a prefilter for eval
        # res_eval._set_pref_udf(self.func, id(self.inputs))

        # This line would NOT allocate physical RAM on any modern OS:
        # aux = np.empty(res_eval.shape, res_eval.dtype)
        # Physical allocation happens here (when writing):
        # res_eval[...] = aux
        # res_eval.schunk.remove_prefilter(self.func.__name__)
        # res_eval.schunk.cparams.nthreads = self._cnthreads

        # return res_eval
        return None

    def __getitem__(self, item):
        if self.chunked_eval:
            # It is important to pass kwargs here, because chunks can be used internally
            # self.kwargs includes self.shape and self.dtype
            return chunked_eval(self.func, self.inputs_dict, item, _getitem=True, **self.kwargs)
        # return self.res_getitem[item] # TODO: implement multithreading
        return None

    def save(self, urlpath=None, **kwargs):
        if urlpath is None:
            raise ValueError("To save a LazyArray you must provide an urlpath")

        meta = kwargs.get("meta", {})
        meta["LazyArray"] = LazyArrayEnum.UDF.value
        kwargs["urlpath"] = urlpath
        kwargs["meta"] = meta
        kwargs["mode"] = "w"  # always overwrite the file in urlpath

        # Create an empty array; useful for providing the shape and dtype of the outcome
        array = blosc2.empty(shape=self.shape, dtype=self.dtype, **kwargs)

        # Save the expression and operands in the metadata
        operands = {}
        operands_ = self.inputs_dict
        for key, value in operands_.items():
            if isinstance(value, blosc2.C2Array):
                operands[key] = {
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
                )
            if value.schunk.urlpath is None:
                raise ValueError("To save a LazyArray, all operands must be stored on disk/network")
            operands[key] = value.schunk.urlpath
        array.schunk.vlmeta["_LazyArray"] = {
            "UDF": textwrap.dedent(inspect.getsource(self.func)).lstrip(),
            "operands": operands,
            "name": self.func.__name__,
        }


def _numpy_eval_expr(expression, operands, prefer_blosc=False):
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
    else:
        ops = {
            key: np.ones(np.ones(len(value.shape), dtype=int), dtype=value.dtype)
            if hasattr(value, "shape")
            else value
            for key, value in operands.items()
        }

    if "contains" in expression:
        _out = ne_evaluate(expression, local_dict=ops)
    else:
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
            _out = ne_evaluate(expression, local_dict=ops)
    return _out


def lazyudf(
    func: Callable[[tuple, np.ndarray, tuple[int]], None],
    inputs: Sequence[Any] | None,
    dtype: np.dtype,
    shape: tuple | list | None = None,
    chunked_eval: bool = True,
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
        - `offset`: The multidimensional offset corresponding to the start of the block being computed.
    inputs: Sequence[Any] or None
        The sequence of inputs. Besides objects compliant with the blosc2.Array protocol,
        any other object is supported too, and it will be passed as-is to the
        user-defined function. If not needed, this can be empty, but `shape` must
        be provided.
    dtype: np.dtype
        The resulting ndarray dtype in NumPy format.
    shape: tuple, optional
        The shape of the resulting array. If None, the shape will be guessed from inputs.
    chunked_eval: bool, optional
        Whether to evaluate the function in chunks or not (blocks).
    kwargs: Any, optional
        Keyword arguments that are supported by the :func:`empty` constructor.
        These arguments will be used by the :meth:`LazyArray.__getitem__` and
        :meth:`LazyArray.compute` methods. The
        last one will ignore the `urlpath` parameter passed in this function.

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
    return LazyUDF(func, inputs, dtype, shape, chunked_eval, **kwargs)


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
            constructor = any(constructor in expression for constructor in constructors)
            if not constructor:
                raise ValueError("No operands nor constructors found in the expression")
            # _new_expr will take care of the constructor, but needs an empty dict in operands
            operands = {}

    return LazyExpr._new_expr(expression, operands, guess=True, out=out, where=where, ne_args=ne_args)


def _open_lazyarray(array):
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
                op = blosc2.open(v)
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
        local_ns = {}
        name = lazyarray["name"]
        filename = f"<{name}>"  # any unique name
        SAFE_GLOBALS = {
            "__builtins__": {
                name: value for name, value in builtins.__dict__.items() if name != "__import__"
            },
            "np": np,
            "blosc2": blosc2,
        }
        if blosc2._HAS_NUMBA:
            SAFE_GLOBALS["numba"] = numba

        # Register the source so inspect can find it
        linecache.cache[filename] = (len(expr), None, expr.splitlines(True), filename)

        exec(compile(expr, filename, "exec"), SAFE_GLOBALS, local_ns)
        func = local_ns[name]
        # TODO: make more robust for general kwargs (not just cparams)
        new_expr = blosc2.lazyudf(
            func,
            tuple(operands_dict[f"o{n}"] for n in range(len(operands_dict))),
            shape=array.shape,
            dtype=array.dtype,
            cparams=array.cparams,
        )

    # Make the array info available for the user (only available when opened from disk)
    new_expr.array = array
    # We want to expose schunk too, so that .info() can be used on the LazyArray
    new_expr.schunk = array.schunk
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
