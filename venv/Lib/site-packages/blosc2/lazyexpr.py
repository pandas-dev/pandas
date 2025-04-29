#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################
# Avoid checking the name of type annotations at run time
from __future__ import annotations

import ast
import asyncio
import concurrent.futures
import copy
import inspect
import math
import os
import pathlib
import re
import sys
import threading
from abc import ABC, abstractmethod
from dataclasses import asdict
from enum import Enum
from pathlib import Path
from queue import Empty, Queue
from typing import TYPE_CHECKING, Any

from numpy.exceptions import ComplexWarning

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

import ndindex
import numpy as np

import blosc2
from blosc2 import compute_chunks_blocks
from blosc2.info import InfoReporter
from blosc2.ndarray import _check_allowed_dtypes, get_chunks_idx, is_inside_new_expr

if not blosc2.IS_WASM:
    import numexpr


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
        # Use numpy eval when running in WebAssembly
        safe_globals = {"np": np}
        # Add all first-level numpy functions
        safe_globals.update(
            {
                name: getattr(np, name)
                for name in dir(np)
                if callable(getattr(np, name)) and not name.startswith("_")
            }
        )
        if "out" in kwargs:
            out = kwargs.pop("out")
            out[:] = eval(expression, safe_globals, local_dict)
            return out
        return eval(expression, safe_globals, local_dict)
    return numexpr.evaluate(expression, local_dict=local_dict, **kwargs)


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

# All the available constructors and reducers necessary for the (string) expression evaluator
constructors = ("arange", "linspace", "fromiter", "zeros", "ones", "empty", "full", "frombuffer")
# Note that, as reshape is accepted as a method too, it should always come last in the list
constructors += ("reshape",)
reducers = ("sum", "prod", "min", "max", "std", "mean", "var", "any", "all")

functions = [
    "sin",
    "cos",
    "tan",
    "sqrt",
    "sinh",
    "cosh",
    "tanh",
    "arcsin",
    "arccos",
    "arctan",
    "arctan2",
    "arcsinh",
    "arccosh",
    "arctanh",
    "exp",
    "expm1",
    "log",
    "log10",
    "log1p",
    "conj",
    "real",
    "imag",
    "contains",
    "abs",
    "sum",
    "prod",
    "mean",
    "std",
    "var",
    "min",
    "max",
    "any",
    "all",
    "pow" if np.__version__.startswith("2.") else "power",
    "where",
]

# Gather all callable functions in numpy
numpy_funcs = {
    name
    for name, member in inspect.getmembers(np, callable)
    if not name.startswith("_") and not isinstance(member, np.ufunc)
}
numpy_ufuncs = {name for name, member in inspect.getmembers(np, lambda x: isinstance(x, np.ufunc))}
# Add these functions to the list of available functions
# (will be evaluated via the array interface)
additional_funcs = sorted((numpy_funcs | numpy_ufuncs) - set(functions))
functions += additional_funcs

functions += constructors

relational_ops = ["==", "!=", "<", "<=", ">", ">="]
logical_ops = ["&", "|", "^", "~"]


def get_expr_globals(expression):
    """Build a dictionary of functions needed for evaluating the expression."""
    _globals = {}

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


class LazyArrayEnum(Enum):
    """
    Available LazyArrays.
    """

    Expr = 0
    UDF = 1


class LazyArray(ABC):
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
    def compute(self, item: slice | list[slice] | None = None, **kwargs: Any) -> blosc2.NDArray:
        """
        Return an :ref:`NDArray` containing the evaluation of the :ref:`LazyArray`.

        Parameters
        ----------
        item: slice, list of slices, optional
            If not None, only the chunks that intersect with the slices
            in items will be evaluated.

        kwargs: Any, optional
            Keyword arguments that are supported by the :func:`empty` constructor.
            These arguments will be set in the resulting :ref:`NDArray`.

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
    def __getitem__(self, item: int | slice | Sequence[slice]) -> blosc2.NDArray:
        """
        Return a NumPy.ndarray containing the evaluation of the :ref:`LazyArray`.

        Parameters
        ----------
        item: int, slice or sequence of slices
            The slice(s) to be retrieved. Note that step parameter is not yet honored.

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
        * All the operands of the LazyArray must be Python scalars, :ref:`NDArray`, :ref:`C2Array` or :ref:`Proxy`.
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

    @property
    @abstractmethod
    def dtype(self) -> np.dtype:
        """
        Get the data type of the :ref:`LazyArray`.

        Returns
        -------
        out: np.dtype
            The data type of the :ref:`LazyArray`.
        """
        pass

    @property
    @abstractmethod
    def shape(self) -> tuple[int]:
        """
        Get the shape of the :ref:`LazyArray`.

        Returns
        -------
        out: tuple
                The shape of the :ref:`LazyArray`.
        """
        pass

    @property
    @abstractmethod
    def ndim(self) -> int:
        """
        Get the number of dimensions of the :ref:`LazyArray`.

        Returns
        -------
        out: int
            The number of dimensions of the :ref:`LazyArray`.
        """
        pass

    @property
    @abstractmethod
    def info(self) -> InfoReporter:
        """
        Get information about the :ref:`LazyArray`.

        Returns
        -------
        out: InfoReporter
            A printable class with information about the :ref:`LazyArray`.
        """
        pass

    # Provide minimal __array_interface__ to allow NumPy to work with this object
    @property
    def __array_interface__(self):
        return {
            "shape": self.shape,
            "typestr": self.dtype.str,
            "data": self[()],
            "version": 3,
        }


def convert_inputs(inputs):
    if not inputs or len(inputs) == 0:
        return []
    inputs_ = []
    for obj in inputs:
        if not isinstance(
            obj, np.ndarray | blosc2.NDArray | blosc2.NDField | blosc2.C2Array
        ) and not np.isscalar(obj):
            try:
                obj = np.asarray(obj)
            except Exception:
                print(
                    "Inputs not being np.ndarray, NDArray, NDField, C2Array or Python scalar objects"
                    " should be convertible to np.ndarray."
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


def check_smaller_shape(value, shape, slice_shape):
    """Check whether the shape of the value is smaller than the shape of the array.

    This follows the NumPy broadcasting rules.
    """
    is_smaller_shape = any(
        s > (1 if i >= len(value.shape) else value.shape[i]) for i, s in enumerate(slice_shape)
    )
    return len(value.shape) < len(shape) or is_smaller_shape


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
    diff_dims = len(larger_shape) - len(smaller_shape)
    return tuple(
        larger_slice[i] if smaller_shape[i - diff_dims] != 1 else slice(0, larger_shape[i])
        for i in range(diff_dims, len(larger_shape))
    )


# Define the patterns for validation
validation_patterns = [
    r"[\;\[\:]",  # Flow control characters
    r"(^|[^\w])__[\w]+__($|[^\w])",  # Dunder methods
    r"\.\b(?!real|imag|(\d*[eE]?[+-]?\d+)|(\d*[eE]?[+-]?\d+j)|\d*j\b|(sum|prod|min|max|std|mean|var|any|all|where)"
    r"\s*\([^)]*\)|[a-zA-Z_]\w*\s*\([^)]*\))",  # Attribute patterns
]

# Compile the blacklist regex
_blacklist_re = re.compile("|".join(validation_patterns))

# Define valid method names
valid_methods = {"sum", "prod", "min", "max", "std", "mean", "var", "any", "all", "where", "reshape"}
valid_methods |= {"int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"}
valid_methods |= {"float32", "float64", "complex64", "complex128"}
valid_methods |= {"bool", "str", "bytes"}


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
    if _blacklist_re.search(skip_quotes) is not None:
        raise ValueError(f"'{expr}' is not a valid expression.")

    # Check for invalid characters not covered by the tokenizer
    invalid_chars = re.compile(r"[^\w\s+\-*/%().,=<>!&|~^]")
    if invalid_chars.search(skip_quotes) is not None:
        invalid_chars = invalid_chars.findall(skip_quotes)
        raise ValueError(f"Expression {expr} contains invalid characters: {invalid_chars}")

    # Check for invalid method names
    method_calls = re.findall(r"\.\b(\w+)\s*\(", skip_quotes)
    for method in method_calls:
        if method not in valid_methods:
            raise ValueError(f"Invalid method name: {method}")


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
    if len(inputs) == 0:
        if out is None:
            raise ValueError(
                "You really want to pass at least one input or one output for building a LazyArray."
                "  Maybe you want blosc2.empty() instead?"
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
    NDinputs = [input for input in inputs if hasattr(input, "chunks")]
    if len(NDinputs) == 0:
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
            raise ValueError("Output shape does not match the first input shape")
        if first_input.chunks != out.chunks:
            fast_path = False
        if first_input.blocks != out.blocks:
            fast_path = False
    # Then, the rest of the operands
    for input_ in NDinputs:
        if first_input.chunks != input_.chunks:
            fast_path = False
        if first_input.blocks != input_.blocks:
            fast_path = False

    return first_input.shape, first_input.chunks, first_input.blocks, fast_path


def is_full_slice(item):
    """Check whether the slice represented by item is a full slice."""
    if item is None:
        # This is the case when the user does not pass any slice in compute() method
        return True
    if isinstance(item, tuple):
        return all((isinstance(i, slice) and i == slice(None, None, None)) or i == Ellipsis for i in item)
    elif isinstance(item, int | bool):
        return False
    else:
        return item == slice(None, None, None) or item == Ellipsis


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


def fill_chunk_operands(  # noqa: C901
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
            arr = next(iter(operands.values()))
            chunks_idx, _ = get_chunks_idx(arr.shape, arr.chunks)
            info = (reduc, aligned, low_mem, chunks_idx)
            iter_chunks = read_nchunk(list(operands.values()), info)
        # Run the asynchronous file reading function from a synchronous context
        chunks = next(iter_chunks)

        for i, (key, value) in enumerate(operands.items()):
            # Chunks are already decompressed, so we can use them directly
            if not low_mem:
                chunk_operands[key] = chunks[i]
                continue
            # Otherwise, we need to decompress them
            special = blosc2.SpecialValue((chunks[i][31] & 0x70) >> 4)
            if special == blosc2.SpecialValue.ZERO:
                # The chunk is a special zero chunk, so we can treat it as a scalar
                chunk_operands[key] = np.zeros((), dtype=value.dtype)
                continue
            if aligned:
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

        if isinstance(value, np.ndarray | blosc2.C2Array):
            chunk_operands[key] = value[slice_]
            continue

        if not full_chunk or not isinstance(value, blosc2.NDArray):
            # The chunk is not a full one, or has padding, or is not a blosc2.NDArray,
            # so we need to go the slow path
            chunk_operands[key] = value[slice_]
            continue

        # First check if the chunk is a special zero chunk.
        # Using lazychunks is very effective here because we only need to read the header.
        chunk = value.schunk.get_lazychunk(nchunk)
        special = blosc2.SpecialValue((chunk[31] & 0x70) >> 4)
        if special == blosc2.SpecialValue.ZERO:
            # The chunk is a special zero chunk, so we can treat it as a scalar
            chunk_operands[key] = np.zeros((), dtype=value.dtype)
            continue

        # If key is in operands, we can reuse the buffer
        if (
            key in chunk_operands
            and chunks_ == chunk_operands[key].shape
            and isinstance(value, blosc2.NDArray)
        ):
            value.get_slice_numpy(chunk_operands[key], (starts, stops))
            continue

        if aligned:
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
    out = kwargs.pop("_output", None)
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    dtype = kwargs.pop("dtype", None)
    where: dict | None = kwargs.pop("_where_args", None)
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
    aligned = blosc2.are_partitions_aligned(shape, chunks, blocks)
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

    # Iterate over the chunks and evaluate the expression
    chunks_idx, nchunks = get_chunks_idx(shape, chunks)
    chunk_operands = {}
    for nchunk in range(nchunks):
        coords = tuple(np.unravel_index(nchunk, chunks_idx))
        slice_ = tuple(
            slice(c * s, min((c + 1) * s, shape[i]))
            for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
        )
        offset = tuple(s.start for s in slice_)  # offset for the udf
        chunks_ = tuple(s.stop - s.start for s in slice_)

        full_chunk = chunks_ == chunks
        fill_chunk_operands(
            operands, slice_, chunks_, full_chunk, aligned, nchunk, iter_disk, chunk_operands
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

            if out is None:
                # We can enter here when using any of the compute() or __getitem__() methods
                if getitem:
                    out = np.empty(shape, dtype=dtype)
                else:
                    out = blosc2.empty(shape, chunks=chunks, blocks=blocks, dtype=dtype, **kwargs)

        # Store the result in the output array
        if getitem:
            try:
                out[slice_] = result
            except ComplexWarning:
                # The result is a complex number, so we need to convert it to real.
                # This is a workaround for rigidness of NumExpr with type casting.
                result = result.real.astype(out.dtype)
                out[slice_] = result
        else:
            if behaved and result.shape == chunks_ and result.dtype == out.dtype:
                # Fast path only works for results that are full chunks
                out.schunk.update_data(nchunk, result, copy=False)
            else:
                out[slice_] = result

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
    _slice=None,
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
        Indicates whether the expression is being evaluated for a getitem operation.
    _slice: slice, list of slices, optional
        If provided, only the chunks that intersect with this slice
        will be evaluated.
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
    if out is None:
        # Compute the shape and chunks of the output array, including broadcasting
        shape = compute_broadcast_shape(operands.values())
    else:
        shape = out.shape

    # We need to keep the original _slice arg, for allowing a final getitem (if necessary)
    orig_slice = _slice

    if chunks is None:
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
            temp = blosc2.empty(shape, dtype=dtype)
            chunks = temp.chunks
            del temp

    # Get the indexes for chunks
    chunks_idx, nchunks = get_chunks_idx(shape, chunks)
    # The starting point for the indices of the inputs
    leninputs = compute_start_index(shape, orig_slice) if orig_slice is not None else 0
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
    for nchunk in range(nchunks):
        coords = tuple(np.unravel_index(nchunk, chunks_idx))
        # Calculate the shape of the (chunk) slice_ (specially at the end of the array)
        slice_ = tuple(
            slice(c * s, min((c + 1) * s, shape[i]))
            for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
        )
        # Check whether current slice_ intersects with _slice
        if _slice is not None and _slice != ():
            # Ensure that _slice is of type slice
            key = ndindex.ndindex(_slice).expand(shape).raw
            _slice = tuple(k if isinstance(k, slice) else slice(k, k + 1, None) for k in key)
            # Ensure that slices do not have any None as start or stop
            _slice = tuple(slice(s.start or 0, s.stop or shape[i], s.step) for i, s in enumerate(_slice))
            slice_ = tuple(slice(s.start or 0, s.stop or shape[i], s.step) for i, s in enumerate(slice_))
            intersects = do_slices_intersect(_slice, slice_)
            if not intersects:
                continue
            # Compute the part of the slice_ that intersects with _slice
            slice_ = tuple(
                slice(max(s1.start, s2.start), min(s1.stop, s2.stop))
                for s1, s2 in zip(slice_, _slice, strict=True)
            )
        slice_shape = tuple(s.stop - s.start for s in slice_)
        len_chunk = math.prod(slice_shape)

        # Get the starts and stops for the slice
        starts = [s.start if s.start is not None else 0 for s in slice_]
        stops = [s.stop if s.stop is not None else sh for s, sh in zip(slice_, slice_shape, strict=True)]

        # Get the slice of each operand
        for key, value in operands.items():
            if np.isscalar(value):
                chunk_operands[key] = value
                continue
            if value.shape == ():
                chunk_operands[key] = value[()]
                continue
            if check_smaller_shape(value, shape, slice_shape):
                # We need to fetch the part of the value that broadcasts with the operand
                smaller_slice = compute_smaller_slice(shape, value.shape, slice_)
                chunk_operands[key] = value[smaller_slice]
                continue
            # If key is in operands, we can reuse the buffer
            if (
                key in chunk_operands
                and slice_shape == chunk_operands[key].shape
                and isinstance(value, blosc2.NDArray)
            ):
                value.get_slice_numpy(chunk_operands[key], (starts, stops))
                continue

            chunk_operands[key] = value[slice_]

        # Evaluate the expression using chunks of operands

        if callable(expression):
            result = np.empty(slice_shape, dtype=out.dtype)
            # Call the udf directly and use result as the output array
            offset = tuple(s.start for s in slice_)
            expression(tuple(chunk_operands.values()), result, offset=offset)
            out[slice_] = result
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
                        slice_shape
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

        if out is None:
            shape_ = shape
            if where is not None and len(where) < 2:
                # The result is a linear array
                shape_ = math.prod(shape)
            if getitem or _order:
                out = np.empty(shape_, dtype=dtype_)
                if _order:
                    indices_ = np.empty(shape_, dtype=np.int64)
            else:
                if "chunks" not in kwargs and (where is None or len(where) == 2):
                    # Let's use the same chunks as the first operand (it could have been automatic too)
                    out = blosc2.empty(shape_, chunks=chunks, dtype=dtype_, **kwargs)
                elif "chunks" in kwargs and (where is not None and len(where) < 2 and len(shape_) > 1):
                    # Remove the chunks argument if the where condition is not a tuple with two elements
                    kwargs.pop("chunks")
                    out = blosc2.empty(shape_, dtype=dtype_, **kwargs)
                else:
                    out = blosc2.empty(shape_, dtype=dtype_, **kwargs)
                # Check if the in out partitions are well-behaved (i.e. no padding)
                behaved = blosc2.are_partitions_behaved(out.shape, out.chunks, out.blocks)

        if where is None or len(where) == 2:
            if behaved and result.shape == out.chunks and result.dtype == out.dtype:
                # Fast path
                out.schunk.update_data(nchunk, result, copy=False)
            else:
                out[slice_] = result
        elif len(where) == 1:
            lenres = len(result)
            out[lenout : lenout + lenres] = result
            if _order is not None:
                indices_[lenout : lenout + lenres] = chunk_indices
            lenout += lenres
        else:
            raise ValueError("The where condition must be a tuple with one or two elements")

    if orig_slice is not None:
        if isinstance(out, np.ndarray):
            out = out[orig_slice]
            if _order is not None:
                indices_ = indices_[orig_slice]
        elif isinstance(out, blosc2.NDArray):
            # It *seems* better to choose an automatic chunks and blocks for the output array
            # out = out.slice(orig_slice, chunks=out.chunks, blocks=out.blocks)
            out = out.slice(orig_slice)
        else:
            raise ValueError("The output array is not a NumPy array or a NDArray")

    if where is not None and len(where) < 2:
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

    return out


def infer_reduction_dtype(dtype, operation):
    # It may change in the future, but for now, this mimics NumPy's (2.1) behavior pretty well
    if operation in {ReduceOp.SUM, ReduceOp.PROD}:
        return np.result_type(dtype, np.int64 if np.issubdtype(dtype, np.integer) else np.float64)
    elif operation in {ReduceOp.MEAN, ReduceOp.STD, ReduceOp.VAR}:
        return np.float64
    elif operation in {ReduceOp.MIN, ReduceOp.MAX}:
        return dtype
    elif operation in {ReduceOp.ANY, ReduceOp.ALL}:
        return np.bool_
    else:
        raise ValueError(f"Unsupported operation: {operation}")


def reduce_slices(  # noqa: C901
    expression: str | Callable[[tuple, np.ndarray, tuple[int]], None],
    operands: dict,
    reduce_args,
    _slice=None,
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
    _slice: slice, list of slices, optional
        If provided, only the chunks that intersect with this slice
        will be evaluated.
    kwargs: Any, optional
        Additional keyword arguments supported by the :func:`empty` constructor.

    Returns
    -------
    :ref:`NDArray` or np.ndarray
        The resulting output array.
    """
    out = kwargs.pop("_output", None)
    ne_args: dict = kwargs.pop("_ne_args", {})
    if ne_args is None:
        ne_args = {}
    where: dict | None = kwargs.pop("_where_args", None)
    reduce_op = reduce_args.pop("op")
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

    if axis is None:
        axis = tuple(range(len(shape)))
    elif not isinstance(axis, tuple):
        axis = (axis,)
    if keepdims:
        reduced_shape = tuple(1 if i in axis else s for i, s in enumerate(shape))
    else:
        reduced_shape = tuple(s for i, s in enumerate(shape) if i not in axis)

    if is_inside_new_expr():
        # We already have the dtype and reduced_shape, so return immediately
        # Use a blosc2 container, as it consumes less memory in general
        return blosc2.zeros(reduced_shape, dtype=dtype)

    # Choose the array with the largest shape as the reference for chunks
    operand = max((o for o in operands.values() if hasattr(o, "chunks")), key=lambda x: len(x.shape))
    chunks = operand.chunks

    # Check if the partitions are aligned (i.e. all operands have the same shape,
    # chunks and blocks, and have no padding). This will allow us to take the fast path.
    same_shape = all(operand.shape == o.shape for o in operands.values() if hasattr(o, "shape"))
    same_chunks = all(operand.chunks == o.chunks for o in operands.values() if hasattr(o, "chunks"))
    same_blocks = all(operand.blocks == o.blocks for o in operands.values() if hasattr(o, "blocks"))
    fast_path = same_shape and same_chunks and same_blocks
    aligned, iter_disk = False, False
    if fast_path:
        # Check that all operands are NDArray for fast path
        all_ndarray = all(
            isinstance(value, blosc2.NDArray) and value.shape != () for value in operands.values()
        )
        # Check that there is some NDArray that is persisted in the disk
        any_persisted = any(
            (isinstance(value, blosc2.NDArray) and value.shape != () and value.schunk.urlpath is not None)
            for value in operands.values()
        )
        if not blosc2.IS_WASM:
            iter_disk = all_ndarray and any_persisted
            # Experiments say that iter_disk is faster than the regular path for reductions
            # even when all operands are in memory, so no need to check any_persisted
            # New benchs are saying the contrary (> 10% slower), so this needs more investigation
            # iter_disk = all_ndarray
        else:
            # WebAssembly does not support threading, so we cannot use the iter_disk option
            iter_disk = False

    # Iterate over the operands and get the chunks
    chunks_idx, nchunks = get_chunks_idx(shape, chunks)
    chunk_operands = {}
    out_init = False
    for nchunk in range(nchunks):
        coords = tuple(np.unravel_index(nchunk, chunks_idx))
        # Calculate the shape of the (chunk) slice_ (specially at the end of the array)
        slice_ = tuple(
            slice(c * s, min((c + 1) * s, shape[i]))
            for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
        )
        if keepdims:
            reduced_slice = tuple(slice(None) if i in axis else sl for i, sl in enumerate(slice_))
        else:
            reduced_slice = tuple(sl for i, sl in enumerate(slice_) if i not in axis)
        offset = tuple(s.start for s in slice_)  # offset for the udf
        # Check whether current slice_ intersects with _slice
        if _slice is not None and _slice != ():
            # Ensure that slices do not have any None as start or stop
            _slice = tuple(slice(s.start or 0, s.stop or shape[i], s.step) for i, s in enumerate(_slice))
            slice_ = tuple(slice(s.start or 0, s.stop or shape[i], s.step) for i, s in enumerate(slice_))
            intersects = do_slices_intersect(_slice, slice_)
            if not intersects:
                continue
            # Compute the part of the slice_ that intersects with _slice
            slice_ = tuple(
                slice(max(s1.start, s2.start), min(s1.stop, s2.stop))
                for s1, s2 in zip(slice_, _slice, strict=True)
            )
        chunks_ = tuple(s.stop - s.start for s in slice_)

        if _slice in (None, ()) and fast_path:
            # Fast path
            full_chunk = chunks_ == chunks
            fill_chunk_operands(
                operands, slice_, chunks_, full_chunk, aligned, nchunk, iter_disk, chunk_operands, reduc=True
            )
        else:
            # Get the starts and stops for the slice
            starts = [s.start if s.start is not None else 0 for s in slice_]
            stops = [s.stop if s.stop is not None else sh for s, sh in zip(slice_, chunks_, strict=True)]
            # Get the slice of each operand
            for key, value in operands.items():
                if np.isscalar(value):
                    chunk_operands[key] = value
                    continue
                if value.shape == ():
                    chunk_operands[key] = value[()]
                    continue
                if check_smaller_shape(value, shape, chunks_):
                    # We need to fetch the part of the value that broadcasts with the operand
                    smaller_slice = compute_smaller_slice(operand.shape, value.shape, slice_)
                    chunk_operands[key] = value[smaller_slice]
                    continue
                # If key is in operands, we can reuse the buffer
                if (
                    key in chunk_operands
                    and chunks_ == chunk_operands[key].shape
                    and isinstance(value, blosc2.NDArray)
                ):
                    value.get_slice_numpy(chunk_operands[key], (starts, stops))
                    continue
                chunk_operands[key] = value[slice_]

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
            if expression == "o0":
                # We don't have an actual expression, so avoid a copy
                result = chunk_operands["o0"]
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
            chunks_ = tuple(s.stop - s.start for s in slice_)
            result = np.full(chunks_, result[()])
        if reduce_op == ReduceOp.ANY:
            result = np.any(result, **reduce_args)
        elif reduce_op == ReduceOp.ALL:
            result = np.all(result, **reduce_args)
        else:
            result = reduce_op.value.reduce(result, **reduce_args)

        if not out_init:
            if out is None:
                out = convert_none_out(result.dtype, reduce_op, reduced_shape)
            else:
                out2 = convert_none_out(result.dtype, reduce_op, reduced_shape)
                out[:] = out2
                del out2
            out_init = True

        # Update the output array with the result
        if reduce_op == ReduceOp.ANY:
            out[reduced_slice] += result
        elif reduce_op == ReduceOp.ALL:
            out[reduced_slice] *= result
        else:
            if reduced_slice == ():
                out = reduce_op.value(out, result)
            else:
                out[reduced_slice] = reduce_op.value(out[reduced_slice], result)

    if out is None:
        if reduce_op in (ReduceOp.MIN, ReduceOp.MAX):
            raise ValueError("zero-size array in min/max reduction operation is not supported")
        if dtype is None:
            # We have no hint here, so choose a default dtype
            dtype = np.float64
        out = convert_none_out(dtype, reduce_op, reduced_shape)

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
    return out


def chunked_eval(  # noqa: C901
    expression: str | Callable[[tuple, np.ndarray, tuple[int]], None], operands: dict, item=None, **kwargs
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
    item: int, slice or sequence of slices, optional
        The slice(s) to be retrieved. Note that step parameter is not honored yet.
    kwargs: Any, optional
        Additional keyword arguments supported by the :func:`empty` constructor.  In addition,
        the following keyword arguments are supported:
        _getitem: bool, optional
            Indicates whether the expression is being evaluated for a getitem operation.
            Default is False.
        _output: NDArray or np.ndarray, optional
            The output array to store the result.
        _ne_args: dict, optional
            Additional arguments to be passed to `numexpr.evaluate()` function.
        _where_args: dict, optional
            Additional arguments for conditional evaluation.
    """
    try:
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

        if not is_full_slice(item) or (where is not None and len(where) < 2):
            # The fast path is not possible when using partial slices or where returning
            # a variable number of elements
            return slices_eval(expression, operands, getitem=getitem, _slice=item, **kwargs)

        if fast_path:
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

        res = slices_eval(expression, operands, getitem=getitem, _slice=item, **kwargs)

    finally:
        # Deactivate cache for NDField instances
        for op in operands:
            if isinstance(operands[op], blosc2.NDField):
                operands[op].ndarr.keep_last_read = False

    return res


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
            if i > 0 and (expr[i - 1] != " " and expr[i - 1] != "("):
                # Not a variable
                new_expr += expr_i
                continue
            # This is a variable.  Find the end of it.
            j = i + 1
            for k in range(len(expr[j:])):
                if expr[j + k] in " )[":
                    j = k
                    break
            if expr[i + j] == ")":
                j -= 1
            old_pos = int(expr[i + 1 : i + j + 1])
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


def infer_dtype(op, value1, value2):
    if op in relational_ops:
        return np.dtype(np.bool_)
    if op in logical_ops:
        # Ensure that both operands are boolean
        if value1.dtype != np.bool_:
            raise ValueError(f"Invalid operand type for {op}: {value1.dtype}")
        if op != "~" and value2.dtype != np.bool_:
            raise ValueError(f"Invalid operand type for {op}: {value2.dtype}")
        return np.dtype(np.bool_)
    dtype1 = value1.dtype if hasattr(value1, "dtype") else np.array(value1).dtype
    dtype2 = value2.dtype if hasattr(value2, "dtype") else np.array(value2).dtype
    return np.result_type(dtype1, dtype2)


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
        if value2 is None:
            if isinstance(value1, LazyExpr):
                self.expression = f"{op}({value1.expression})"
                self.operands = value1.operands
            else:
                self.operands = {"o0": value1}
                self.expression = "o0" if op is None else f"{op}(o0)"
            return
        elif op in ("arctan2", "contains", "pow", "power"):
            if np.isscalar(value1) and np.isscalar(value2):
                self.expression = f"{op}(o0, o1)"
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

        self._dtype = infer_dtype(op, value1, value2)
        if np.isscalar(value1) and np.isscalar(value2):
            self.expression = f"({value1} {op} {value2})"
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
            elif isinstance(value1, LazyExpr) or isinstance(value2, LazyExpr):
                if isinstance(value1, LazyExpr):
                    self.expression = value1.expression
                    self.operands = {"o0": value2}
                else:
                    self.expression = value2.expression
                    self.operands = {"o0": value1}
                newexpr = self.update_expr(new_op)
                self.expression = newexpr.expression
                self.operands = newexpr.operands
            else:
                # This is the very first time that a LazyExpr is formed from two operands
                # that are not LazyExpr themselves
                self.operands = {"o0": value1, "o1": value2}
                self.expression = f"(o0 {op} o1)"

    def get_chunk(self, nchunk):
        """Get the `nchunk` of the expression, evaluating only that one."""
        # Create an empty array with the same shape and dtype; this is fast
        out = blosc2.empty(shape=self.shape, dtype=self.dtype, chunks=self.chunks, blocks=self.blocks)
        shape = out.shape
        chunks = out.chunks
        # Calculate the shape of the (chunk) slice_ (specially at the end of the array)
        chunks_idx, _ = get_chunks_idx(shape, chunks)
        coords = tuple(np.unravel_index(nchunk, chunks_idx))
        slice_ = tuple(
            slice(c * s, min((c + 1) * s, shape[i]))
            for i, (c, s) in enumerate(zip(coords, chunks, strict=True))
        )
        # TODO: we need more metadata for treating reductions
        # We want to fill a single chunk, so we need to evaluate the expression on out
        expr = lazyexpr(self, out=out)
        # The evals below produce arrays with different chunks and blocks;
        # we choose the ones for LazyExpr main class
        expr.compute(item=slice_)
        # out = expr.compute(item=slice_)
        return out.schunk.get_chunk(nchunk)

    def update_expr(self, new_op):  # noqa: C901
        # We use a lot of the original NDArray.__eq__ as 'is', so deactivate the overloaded one
        blosc2._disable_overloaded_equal = True
        # One of the two operands are LazyExpr instances
        value1, op, value2 = new_op
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

        self._dtype = infer_dtype(op, value1, value2)
        if not isinstance(value1, LazyExpr) and not isinstance(value2, LazyExpr):
            # We converted some of the operands to NDArray (where() handling above)
            new_operands = {"o0": value1, "o1": value2}
            expression = f"(o0 {op} o1)"
            return self._new_expr(expression, new_operands, guess=False, out=None, where=None)
        elif isinstance(value1, LazyExpr) and isinstance(value2, LazyExpr):
            # Expression fusion
            # Fuse operands in expressions and detect duplicates
            new_operands, dup_op = fuse_operands(value1.operands, value2.operands)
            # Take expression 2 and rebase the operands while removing duplicates
            new_expr = fuse_expressions(value2.expression, len(value1.operands), dup_op)
            expression = f"({self.expression} {op} {new_expr})"
        elif isinstance(value1, LazyExpr):
            if op == "~":
                expression = f"({op}{self.expression})"
            elif np.isscalar(value2):
                expression = f"({self.expression} {op} {value2})"
            elif hasattr(value2, "shape") and value2.shape == ():
                expression = f"({self.expression} {op} {value2[()]})"
            else:
                operand_to_key = {id(v): k for k, v in value1.operands.items()}
                try:
                    op_name = operand_to_key[id(value2)]
                except KeyError:
                    op_name = f"o{len(self.operands)}"
                    new_operands = {op_name: value2}
                expression = f"({self.expression} {op} {op_name})"
            self.operands = value1.operands
        else:
            if np.isscalar(value1):
                expression = f"({value1} {op} {self.expression})"
            elif hasattr(value1, "shape") and value1.shape == ():
                expression = f"({value1[()]} {op} {self.expression})"
            else:
                operand_to_key = {id(v): k for k, v in value2.operands.items()}
                try:
                    op_name = operand_to_key[id(value1)]
                except KeyError:
                    op_name = f"o{len(value2.operands)}"
                    new_operands = {op_name: value1}
                if op == "[]":  # syntactic sugar for slicing
                    expression = f"({op_name}[{self.expression}])"
                else:
                    expression = f"({op_name} {op} {self.expression})"
                self.operands = value2.operands
        blosc2._disable_overloaded_equal = False
        # Return a new expression
        operands = self.operands | new_operands
        return self._new_expr(expression, operands, guess=False, out=None, where=None)

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
        operands = {
            key: np.ones(np.ones(len(value.shape), dtype=int), dtype=value.dtype)
            for key, value in self.operands.items()
        }
        if "contains" in self.expression:
            _out = ne_evaluate(self.expression, local_dict=operands)
        else:
            # Create a globals dict with the functions of numpy
            globals_dict = {f: getattr(np, f) for f in functions if f not in ("contains", "pow")}
            try:
                _out = eval(self.expression, globals_dict, operands)
            except RuntimeWarning:
                # Sometimes, numpy gets a RuntimeWarning when evaluating expressions
                # with synthetic operands (1's). Let's try with numexpr, which is not so picky
                # about this.
                _out = ne_evaluate(self.expression, local_dict=operands)
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
        # Operands shape can change, so we always need to recompute this
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
        self._shape, self._chunks, self._blocks, fast_path = validate_inputs(
            self.operands, getattr(self, "_out", None)
        )
        if not fast_path:
            # Not using the fast path, so we need to compute the chunks/blocks automatically
            self._chunks, self._blocks = compute_chunks_blocks(self.shape, None, None, dtype=self.dtype)
        return self._chunks

    @property
    def blocks(self):
        if hasattr(self, "_blocks"):
            return self._blocks
        self._shape, self._chunks, self._blocks, fast_path = validate_inputs(
            self.operands, getattr(self, "_out", None)
        )
        if not fast_path:
            # Not using the fast path, so we need to compute the chunks/blocks automatically
            self._chunks, self._blocks = compute_chunks_blocks(self.shape, None, None, dtype=self.dtype)
        return self._blocks

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Handle operations at the array level
        if method != "__call__":
            return NotImplemented

        ufunc_map = {
            np.add: "+",
            np.subtract: "-",
            np.multiply: "*",
            np.divide: "/",
            np.true_divide: "/",
            np.power: "**",
            np.less: "<",
            np.less_equal: "<=",
            np.greater: ">",
            np.greater_equal: ">=",
            np.equal: "==",
            np.not_equal: "!=",
            np.bitwise_and: "&",
            np.bitwise_or: "|",
            np.bitwise_xor: "^",
        }

        ufunc_map_1param = {
            np.sqrt: "sqrt",
            np.sin: "sin",
            np.cos: "cos",
            np.tan: "tan",
            np.arcsin: "arcsin",
            np.arccos: "arccos",
            np.arctan: "arctan",
            np.sinh: "sinh",
            np.cosh: "cosh",
            np.tanh: "tanh",
            np.arcsinh: "arcsinh",
            np.arccosh: "arccosh",
            np.arctanh: "arctanh",
            np.exp: "exp",
            np.expm1: "expm1",
            np.log: "log",
            np.log10: "log10",
            np.log1p: "log1p",
            np.abs: "abs",
            np.conj: "conj",
            np.real: "real",
            np.imag: "imag",
            np.bitwise_not: "~",
        }

        if ufunc in ufunc_map:
            value = inputs[0] if inputs[1] is self else inputs[1]
            _check_allowed_dtypes(value)
            return blosc2.LazyExpr(new_op=(value, ufunc_map[ufunc], self))

        if ufunc in ufunc_map_1param:
            value = inputs[0]
            _check_allowed_dtypes(value)
            return blosc2.LazyExpr(new_op=(value, ufunc_map_1param[ufunc], None))

        return NotImplemented

    def __neg__(self):
        return self.update_expr(new_op=(0, "-", self))

    def __add__(self, value):
        return self.update_expr(new_op=(self, "+", value))

    def __iadd__(self, other):
        return self.update_expr(new_op=(self, "+", other))

    def __radd__(self, value):
        return self.update_expr(new_op=(value, "+", self))

    def __sub__(self, value):
        return self.update_expr(new_op=(self, "-", value))

    def __isub__(self, value):
        return self.update_expr(new_op=(self, "-", value))

    def __rsub__(self, value):
        return self.update_expr(new_op=(value, "-", self))

    def __mul__(self, value):
        return self.update_expr(new_op=(self, "*", value))

    def __imul__(self, value):
        return self.update_expr(new_op=(self, "*", value))

    def __rmul__(self, value):
        return self.update_expr(new_op=(value, "*", self))

    def __truediv__(self, value):
        return self.update_expr(new_op=(self, "/", value))

    def __itruediv__(self, value):
        return self.update_expr(new_op=(self, "/", value))

    def __rtruediv__(self, value):
        return self.update_expr(new_op=(value, "/", self))

    def __and__(self, value):
        return self.update_expr(new_op=(self, "&", value))

    def __rand__(self, value):
        return self.update_expr(new_op=(value, "&", self))

    def __or__(self, value):
        return self.update_expr(new_op=(self, "|", value))

    def __ror__(self, value):
        return self.update_expr(new_op=(value, "|", self))

    def __invert__(self):
        return self.update_expr(new_op=(self, "~", None))

    def __pow__(self, value):
        return self.update_expr(new_op=(self, "**", value))

    def __rpow__(self, value):
        return self.update_expr(new_op=(value, "**", self))

    def __ipow__(self, value):
        return self.update_expr(new_op=(self, "**", value))

    def __lt__(self, value):
        return self.update_expr(new_op=(self, "<", value))

    def __le__(self, value):
        return self.update_expr(new_op=(self, "<=", value))

    def __eq__(self, value):
        return self.update_expr(new_op=(self, "==", value))

    def __ne__(self, value):
        return self.update_expr(new_op=(self, "!=", value))

    def __gt__(self, value):
        return self.update_expr(new_op=(self, ">", value))

    def __ge__(self, value):
        return self.update_expr(new_op=(self, ">=", value))

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
            dtype = np.result_type(value1, value2)
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

    def sum(self, axis=None, dtype=None, keepdims=False, **kwargs):
        reduce_args = {
            "op": ReduceOp.SUM,
            "axis": axis,
            "dtype": dtype,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, **kwargs)

    def prod(self, axis=None, dtype=None, keepdims=False, **kwargs):
        reduce_args = {
            "op": ReduceOp.PROD,
            "axis": axis,
            "dtype": dtype,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, **kwargs)

    def get_num_elements(self, axis, item):
        if hasattr(self, "_where_args") and len(self._where_args) == 1:
            # We have a where condition, so we need to count the number of elements
            # fulfilling the condition
            orig_where_args = self._where_args
            self._where_args = {"_where_x": blosc2.ones(self.shape, dtype=np.int8)}
            num_elements = self.sum(axis=axis, dtype=np.int64, item=item)
            self._where_args = orig_where_args
            return num_elements
        if np.isscalar(axis):
            axis = (axis,)
        # Compute the number of elements in the array
        shape = self.shape
        if item is not None:
            # Compute the shape of the slice
            if not isinstance(item, tuple):
                item = (item,)
            # Ensure that the limits in item slices are not None
            item = tuple(slice(s.start or 0, s.stop or self.shape[i], s.step) for i, s in enumerate(item))
            # Compute the intersection of the slice with the shape
            item = tuple(slice(s1.start, min(s1.stop, s2)) for s1, s2 in zip(item, shape, strict=True))
            if axis is None:
                shape = [s.stop - s.start for s in item]
            else:
                shape = [s.stop - s.start for i, s in enumerate(item) if i in axis]
        return math.prod(shape) if axis is None else math.prod([shape[i] for i in axis])

    def mean(self, axis=None, dtype=None, keepdims=False, **kwargs):
        item = kwargs.pop("item", None)
        total_sum = self.sum(axis=axis, dtype=dtype, keepdims=keepdims, item=item)
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

    def std(self, axis=None, dtype=None, keepdims=False, ddof=0, **kwargs):
        item = kwargs.pop("item", None)
        mean_value = self.mean(axis=axis, dtype=dtype, keepdims=True, item=item)
        expr = (self - mean_value) ** 2
        out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, item=item)
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

    def var(self, axis=None, dtype=None, keepdims=False, ddof=0, **kwargs):
        item = kwargs.pop("item", None)
        mean_value = self.mean(axis=axis, dtype=dtype, keepdims=True, item=item)
        expr = (self - mean_value) ** 2
        if ddof != 0:
            out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, item=item)
            num_elements = self.get_num_elements(axis, item)
            out = out * num_elements / (num_elements - ddof)
        else:
            out = expr.mean(axis=axis, dtype=dtype, keepdims=keepdims, item=item)
        out2 = kwargs.pop("out", None)
        if out2 is not None:
            out2[:] = out
            return out2
        if kwargs != {} and not np.isscalar(out):
            out = blosc2.asarray(out, **kwargs)
        return out

    def min(self, axis=None, keepdims=False, **kwargs):
        reduce_args = {
            "op": ReduceOp.MIN,
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, **kwargs)

    def max(self, axis=None, keepdims=False, **kwargs):
        reduce_args = {
            "op": ReduceOp.MAX,
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, **kwargs)

    def any(self, axis=None, keepdims=False, **kwargs):
        reduce_args = {
            "op": ReduceOp.ANY,
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, **kwargs)

    def all(self, axis=None, keepdims=False, **kwargs):
        reduce_args = {
            "op": ReduceOp.ALL,
            "axis": axis,
            "keepdims": keepdims,
        }
        return self.compute(_reduce_args=reduce_args, **kwargs)

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

    def _compute_expr(self, item, kwargs):  # noqa: C901
        if any(method in self.expression for method in reducers):
            # We have reductions in the expression (probably coming from a string lazyexpr)
            _globals = get_expr_globals(self.expression)
            lazy_expr = eval(self.expression, _globals, self.operands)
            if not isinstance(lazy_expr, blosc2.LazyExpr):
                # An immediate evaluation happened (e.g. all operands are numpy arrays)
                if hasattr(self, "_where_args"):
                    # We need to apply the where() operation
                    if len(self._where_args) == 1:
                        # We have a single argument
                        where_x = self._where_args["_where_x"]
                        return where_x[:][lazy_expr]
                    if len(self._where_args) == 2:
                        # We have two arguments
                        where_x = self._where_args["_where_x"]
                        where_y = self._where_args["_where_y"]
                        return np.where(lazy_expr, where_x, where_y)
                if hasattr(self, "_output"):
                    # This is not exactly optimized, but it works for now
                    self._output[:] = lazy_expr
                return lazy_expr

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

            _globals = {func: getattr(blosc2, func) for func in functions if func in newexpr}
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

    def compute(self, item=None, **kwargs) -> blosc2.NDArray:
        # When NumPy ufuncs are called, the user may add an `out` parameter to kwargs
        if "out" in kwargs:
            kwargs["_output"] = kwargs.pop("out")
        if hasattr(self, "_output"):
            kwargs["_output"] = self._output
        if "ne_args" in kwargs:
            kwargs["_ne_args"] = kwargs.pop("ne_args")
        if hasattr(self, "_ne_args"):
            kwargs["_ne_args"] = self._ne_args
        if hasattr(self, "_where_args"):
            kwargs["_where_args"] = self._where_args
        kwargs["dtype"] = self.dtype
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
            kwargs_not_accepted = {"_where_args", "_indices", "_order", "_ne_args", "dtype"}
            kwargs = {key: value for key, value in kwargs.items() if key not in kwargs_not_accepted}
            result = blosc2.asarray(result, **kwargs)
        return result

    def __getitem__(self, item):
        kwargs = {"_getitem": True}
        return self.compute(item, **kwargs)

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

    def _save(self, urlpath=None, **kwargs):
        if urlpath is None:
            raise ValueError("To save a LazyArray you must provide an urlpath")

        # Validate expression
        validate_expr(self.expression)

        meta = kwargs.get("meta", {})
        meta["LazyArray"] = LazyArrayEnum.Expr.value
        kwargs["urlpath"] = urlpath
        kwargs["meta"] = meta
        kwargs["mode"] = "w"  # always overwrite the file in urlpath

        # Create an empty array; useful for providing the shape and dtype of the outcome
        array = blosc2.empty(shape=self.shape, dtype=self.dtype, **kwargs)

        # Save the expression and operands in the metadata
        operands = {}
        for key, value in self.operands.items():
            if isinstance(value, blosc2.C2Array):
                operands[key] = {
                    "path": str(value.path),
                    "urlbase": value.urlbase,
                }
                continue
            if key in {"numpy", "np"}:
                # Provide access to cast funcs like int8 et al.
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
            "expression": self.expression,
            "UDF": None,
            "operands": operands,
        }

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
        if guess:
            # The expression has been validated, so we can evaluate it
            # in guessing mode to avoid computing reductions
            # Extract possible numpy scalars
            _expression, local_vars = extract_numpy_scalars(expression)
            # Let's include numpy and blosc2 as operands so that some functions can be used
            # Most in particular, castings like np.int8 et al. can be very useful to allow
            # for desired data types in the output.
            _operands = operands | local_vars
            _globals = get_expr_globals(expression)
            _globals |= dtype_symbols
            new_expr = eval(_expression, _globals, _operands)
            _dtype = new_expr.dtype
            _shape = new_expr.shape
            if isinstance(new_expr, blosc2.LazyExpr):
                # Restore the original expression and operands
                new_expr.expression = _expression
                new_expr.expression_tosave = expression
                new_expr.operands = _operands
                new_expr.operands_tosave = operands
            else:
                # An immediate evaluation happened (e.g. all operands are numpy arrays)
                new_expr = cls(None)
                new_expr.expression = expression
                new_expr.operands = operands
            # Cache the dtype and shape (should be immutable)
            new_expr._dtype = _dtype
            new_expr._shape = _shape
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
        self.chunked_eval = chunked_eval
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

        self.res_getitem = blosc2.empty(self._shape, self._dtype, **kwargs_getitem)
        # Register a postfilter for getitem
        self.res_getitem._set_postf_udf(self.func, id(self.inputs))

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
        items = []
        items += [("type", f"{self.__class__.__name__}")]
        inputs = {}
        for key, value in self.inputs_dict.items():
            if isinstance(value, np.ndarray | blosc2.NDArray | blosc2.C2Array):
                inputs[key] = f"<{value.__class__.__name__}> {value.shape} {value.dtype}"
            else:
                inputs[key] = str(value)
        items += [("inputs", inputs)]
        items += [("shape", self.shape)]
        items += [("dtype", self.dtype)]
        return items

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

    def compute(self, item=None, **kwargs):
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
        dparams = aux_kwargs.get("dparams", {})
        if isinstance(dparams, blosc2.DParams):
            # Convert to dictionary
            dparams = asdict(dparams)
        dparams.update(kwargs.get("dparams", {}))
        aux_kwargs["dparams"] = dparams
        _ = kwargs.pop("cparams", None)
        _ = kwargs.pop("dparams", None)
        urlpath = kwargs.get("urlpath")
        if urlpath is not None and urlpath == aux_kwargs.get(
            "urlpath",
        ):
            raise ValueError("Cannot use the same urlpath for LazyArray and eval NDArray")
        _ = aux_kwargs.pop("urlpath", None)
        aux_kwargs.update(kwargs)

        if item is None:
            if self.chunked_eval:
                res_eval = blosc2.empty(self.shape, self.dtype, **aux_kwargs)
                chunked_eval(self.func, self.inputs_dict, None, _getitem=False, _output=res_eval)
                return res_eval

            # Cannot use multithreading when applying a prefilter, save nthreads to set them
            # after the evaluation
            cparams = aux_kwargs.get("cparams", {})
            if isinstance(cparams, dict):
                self._cnthreads = cparams.get("nthreads", blosc2.cparams_dflts["nthreads"])
                cparams["nthreads"] = 1
            else:
                raise ValueError("cparams should be a dictionary")
            aux_kwargs["cparams"] = cparams

            res_eval = blosc2.empty(self.shape, self.dtype, **aux_kwargs)
            # Register a prefilter for eval
            res_eval._set_pref_udf(self.func, id(self.inputs))

            aux = np.empty(res_eval.shape, res_eval.dtype)
            res_eval[...] = aux
            res_eval.schunk.remove_prefilter(self.func.__name__)
            res_eval.schunk.cparams.nthreads = self._cnthreads

            return res_eval
        else:
            # Get only a slice
            np_array = self.__getitem__(item)
            return blosc2.asarray(np_array, **aux_kwargs)

    def __getitem__(self, item):
        if self.chunked_eval:
            output = np.empty(self.shape, self.dtype)
            # It is important to pass kwargs here, because chunks can be used internally
            chunked_eval(self.func, self.inputs_dict, item, _getitem=True, _output=output, **self.kwargs)
            return output[item]
        return self.res_getitem[item]

    def save(self, **kwargs):
        raise NotImplementedError("For safety reasons, this is not implemented for UDFs")


def lazyudf(
    func: Callable[[tuple, np.ndarray, tuple[int]], None],
    inputs: tuple | list | None,
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
    inputs: tuple or list or None
        The sequence of inputs. Supported inputs are:
        NumPy.ndarray, :ref:`NDArray`, :ref:`NDField`, :ref:`C2Array`.
        Any other object is supported too, and will be passed as is to the user-defined function.
        If not needed, this can be empty, but `shape` must be provided.
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
    expression: str | bytes | LazyExpr | blosc2.NDArray,
    operands: dict | None = None,
    out: blosc2.NDArray | np.ndarray = None,
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
    expression: str or bytes or LazyExpr
        The expression to evaluate. This can be any valid expression that can be
        ingested by numexpr. If a LazyExpr is passed, the expression will be
        updated with the new operands.
    operands: dict
        The dictionary with operands. Supported values are NumPy.ndarray,
        Python scalars, :ref:`NDArray`, :ref:`NDField` or :ref:`C2Array` instances.
        If None, the operands will be seeked in the local and global dictionaries.
    out: NDArray or np.ndarray, optional
        The output array where the result will be stored. If not provided,
        a new array will be created.
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
    >>> expr = 'a1 * b1 + 2'
    >>> operands = { 'a': a1, 'b': b1 }
    >>> lazy_expr = blosc2.lazyexpr(expr, operands=operands)
    >>> f"Lazy expression created: {lazy_expr}"
    Lazy expression created: a1 * b1 + 2
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
    if value == LazyArrayEnum.UDF.value:
        raise NotImplementedError("For safety reasons, persistent UDFs are not supported")

    # LazyExpr
    lazyarray = array.schunk.vlmeta["_LazyArray"]
    operands = lazyarray["operands"]
    parent_path = Path(array.schunk.urlpath).parent
    operands_dict = {}
    for key, value in operands.items():
        if isinstance(value, str):
            value = parent_path / value
            op = blosc2.open(value)
            operands_dict[key] = op
        elif isinstance(value, dict):
            # C2Array
            operands_dict[key] = blosc2.C2Array(
                pathlib.Path(value["path"]).as_posix(),
                urlbase=value["urlbase"],
            )
        else:
            raise TypeError("Error when retrieving the operands")

    expr = lazyarray["expression"]
    new_expr = LazyExpr._new_expr(expr, operands_dict, guess=True, out=None, where=None)
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
    out: np.ndarray | blosc2.NDArray = None,
    **kwargs: Any,
) -> np.ndarray | blosc2.NDArray:
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
    out: NDArray or NumPy array, optional
        The output array where the result will be stored. If not provided,
        a new NumPy array will be created and returned.
    kwargs: Any, optional
        Additional arguments to be passed to `numexpr.evaluate()` function.

    Returns
    -------
    out: NumPy or NDArray
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
