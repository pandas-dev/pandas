#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

import ast
import builtins
import contextlib
import inspect
import math
import sys
import warnings
from itertools import product

import ndindex
import numpy as np
from ndindex.subindex_helpers import ceiling
from numpy import broadcast_shapes

import blosc2

# Set this to False if miniexpr should not be tried out
try_miniexpr = not blosc2.IS_WASM or getattr(blosc2, "_WASM_MINIEXPR_ENABLED", False)


def _toggle_miniexpr(FLAG):
    global try_miniexpr
    try_miniexpr = FLAG
    for module_name in ("blosc2.lazyexpr", "blosc2.linalg"):
        module = sys.modules.get(module_name)
        if module is not None:
            module.try_miniexpr = FLAG


# NumPy version and a convenient boolean flag
NUMPY_GE_2_0 = np.__version__ >= "2.0"
# handle different numpy versions
if NUMPY_GE_2_0:  # array-api compliant
    nplshift = np.bitwise_left_shift
    nprshift = np.bitwise_right_shift
    npbinvert = np.bitwise_invert
    npvecdot = np.vecdot
    nptranspose = np.permute_dims
    if hasattr(np, "cumulative_sum"):
        npcumsum = np.cumulative_sum
        npcumprod = np.cumulative_prod
    else:
        npcumsum = np.cumsum
        npcumprod = np.cumprod
else:  # not array-api compliant
    nplshift = np.left_shift
    nprshift = np.right_shift
    npbinvert = np.bitwise_not
    nptranspose = np.transpose
    npcumsum = np.cumsum
    npcumprod = np.cumprod

    def npvecdot(a, b, axis=-1):
        return np.einsum("...i,...i->...", np.moveaxis(np.conj(a), axis, -1), np.moveaxis(b, axis, -1))


def _string_contains(a, b):
    return np.char.find(a, b) >= 0


def _string_startswith(a, b):
    return np.char.startswith(a, b)


def _string_lower(a):
    return np.char.lower(a)


def _string_upper(a):
    return np.char.upper(a)


def _string_endswith(a, b):
    return np.char.endswith(a, b)


def format_expr_scalar(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, str | bytes):
        return repr(value)
    return value


elementwise_funcs = [
    "abs",
    "acos",
    "acosh",
    "add",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctan2",
    "arctanh",
    "asin",
    "asinh",
    "atan",
    "atan2",
    "atanh",
    "bitwise_and",
    "bitwise_invert",
    "bitwise_left_shift",
    "bitwise_or",
    "bitwise_right_shift",
    "bitwise_xor",
    "broadcast_to",
    "ceil",
    "clip",
    "conj",
    "contains",
    "copysign",
    "cos",
    "cosh",
    "divide",
    "endswith",
    "equal",
    "exp",
    "expm1",
    "floor",
    "floor_divide",
    "greater",
    "greater_equal",
    "hypot",
    "imag",
    "isfinite",
    "isinf",
    "isnan",
    "less_equal",
    "less",
    "log",
    "log1p",
    "log2",
    "log10",
    "logaddexp",
    "logical_and",
    "logical_not",
    "logical_or",
    "logical_xor",
    "lower",
    "maximum",
    "minimum",
    "multiply",
    "negative",
    "nextafter",
    "not_equal",
    "positive",
    "pow",
    "real",
    "reciprocal",
    "remainder",
    "round",
    "sign",
    "signbit",
    "sin",
    "sinh",
    "sqrt",
    "square",
    "startswith",
    "subtract",
    "tan",
    "tanh",
    "trunc",
    "upper",
    "where",
]

linalg_funcs = [
    "concat",
    "diagonal",
    "expand_dims",
    "matmul",
    "matrix_transpose",
    "outer",
    "permute_dims",
    "squeeze",
    "stack",
    "tensordot",
    "transpose",
    "vecdot",
]

linalg_attrs = ["T", "mT"]
reducers = [
    "sum",
    "prod",
    "min",
    "max",
    "std",
    "mean",
    "var",
    "any",
    "all",
    "count_nonzero",
    "argmax",
    "argmin",
    "cumulative_sum",
    "cumulative_prod",
]

# All the available constructors and reducers necessary for the (string) expression evaluator
constructors = [
    "asarray",
    "arange",
    "copy",
    "linspace",
    "fromiter",
    "zeros",
    "ones",
    "empty",
    "full",
    "frombuffer",
    "full_like",
    "zeros_like",
    "ones_like",
    "empty_like",
    "eye",
    "nans",
    "ndarray_from_cframe",
    "uninit",
    "meshgrid",
]

# Note that, as reshape is accepted as a method too, it should always come last in the list
constructors += ["reshape"]


_NUMPY_ALIASES = {
    "acos": np.arccos,
    "acosh": np.arccosh,
    "asin": np.arcsin,
    "asinh": np.arcsinh,
    "atan": np.arctan,
    "atanh": np.arctanh,
    "atan2": np.arctan2,
    "concat": getattr(np, "concat", np.concatenate),
    "contains": _string_contains,
    "cumulative_prod": npcumprod,
    "cumulative_sum": npcumsum,
    "endswith": _string_endswith,
    "lower": _string_lower,
    "matrix_transpose": getattr(np, "matrix_transpose", np.transpose),
    "permute_dims": nptranspose,
    "pow": np.power,
    "startswith": _string_startswith,
    "upper": _string_upper,
    "vecdot": npvecdot,
}
if not NUMPY_GE_2_0:  # handle non-array-api compliance
    _NUMPY_ALIASES.update(
        {
            "bitwise_invert": np.bitwise_not,
            "bitwise_left_shift": np.left_shift,
            "bitwise_right_shift": np.right_shift,
        }
    )

# Use numpy eval when running in WebAssembly.  Keep this intentionally small:
# scanning every callable in numpy triggers lazy imports such as numpy.f2py and
# numpy.testing during ``import blosc2``.
safe_numpy_globals = {"np": np, "nan": np.nan, "inf": np.inf, **_NUMPY_ALIASES}
for _name in set(elementwise_funcs + linalg_funcs + reducers + constructors):
    if _name not in safe_numpy_globals and not _name.startswith("_"):
        with contextlib.suppress(AttributeError):
            _value = getattr(np, _name)
            if callable(_value):
                safe_numpy_globals[_name] = _value


def populate_safe_numpy_globals(expression: str) -> None:
    """Add bare numpy call names used by *expression* to safe_numpy_globals."""
    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError:
        return
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        name = node.func.id
        if name in safe_numpy_globals or name.startswith("_"):
            continue
        with contextlib.suppress(AttributeError):
            value = getattr(np, name)
            if callable(value):
                safe_numpy_globals[name] = value


# --- Shape utilities ---
def linalg_shape(func_name, args, kwargs):  # noqa: C901
    # --- Linear algebra and tensor manipulation ---
    a = args[0] if args else None
    if a is None or any(s is None for s in a):
        return None
    b = args[1] if len(args) > 1 else None
    axis = kwargs.get("axis", None)
    axes = kwargs.get("axes", None)
    offset = kwargs.get("offset", 0)

    # --- concat ---
    if func_name == "concat":
        shapes = args[0]
        if axis is None and len(args) > 1:
            axis = args[1]

        # Coerce axis to int if tuple single-element
        axis = 0 if axis is None else axis
        # normalize negative axis
        axis = axis + len(shapes[0]) if axis < 0 else axis
        concat_dim = builtins.sum(s[axis] for s in shapes)
        return tuple(s if i != axis else concat_dim for i, s in enumerate(shapes[0]))

    # --- diagonal ---
    elif func_name == "diagonal":
        axis1 = len(a) - 2
        axis2 = len(a) - 1
        new_shape = [d for i, d in enumerate(a) if i not in (axis1, axis2)]
        d1, d2 = a[axis1], a[axis2]
        diag_len = builtins.max(0, min(d1, d2) - abs(offset))
        new_shape.append(diag_len)
        return tuple(new_shape)

    # --- expand_dims ---
    elif func_name == "expand_dims":
        # positional axis may be second positional argument
        if axis is None and len(args) > 1:
            axis = args[1]
        if axis is None:
            axis = 0
        axis = [axis] if isinstance(axis, int) else axis
        new_shape = list(a)
        for ax in sorted(axis):
            ax = ax if ax >= 0 else len(new_shape) + ax + 1
            new_shape.insert(ax, 1)
        return tuple(new_shape)

    # --- matmul ---
    elif func_name == "matmul":
        if b is None:
            return None
        x1_is_vector = False
        x2_is_vector = False
        if len(a) == 1:
            a = (1,) + a  # (N,) -> (1, N)
            x1_is_vector = True
        if len(b) == 1:
            b += (1,)  # (M,) -> (M, 1)
            x2_is_vector = True
        batch = broadcast_shapes(a[:-2], b[:-2])
        shape = batch
        if not x1_is_vector:
            shape += (a[-2],)
        if not x2_is_vector:
            shape += (b[-1],)
        return shape

    # --- matrix_transpose ---
    elif func_name == "matrix_transpose":
        if len(a) < 2:
            return a
        return a[:-2] + (a[-1], a[-2])

    # --- outer ---
    elif func_name == "outer":
        if b is None:
            return None
        return a + b

    # --- permute_dims ---
    elif func_name == "permute_dims":
        if axes is None and len(args) > 1:
            axes = args[1]
        if axes is None:
            axes = tuple(reversed(range(len(a))))
        return tuple(a[i] for i in axes)

    # --- squeeze ---
    elif func_name == "squeeze":
        if axis is None and len(args) > 1:
            axis = args[1]
        if axis is None:
            return tuple(d for d in a if d != 1)
        if isinstance(axis, int):
            axis = (axis,)
        axis = tuple(ax if ax >= 0 else len(a) + ax for ax in axis)
        return tuple(d for i, d in enumerate(a) if i not in axis or d != 1)

    # --- stack ---
    elif func_name == "stack":
        # detect axis as last positional if candidate
        elems = args[0]
        if axis is None and len(args) > 1:
            axis = args[1]
        if axis is None:
            axis = 0
        if axis < 0:
            axis += len(elems[0]) + 1
        return elems[0][:axis] + (len(elems),) + elems[0][axis:]

    # --- tensordot ---
    elif func_name == "tensordot":
        if axes is None and len(args) > 2:
            axes = args[2]
        if axes is None:
            axes = 2
        if b is None:
            return None
        if isinstance(axes, int):
            a_rest = a[:-axes]
            b_rest = b[axes:]
        else:
            a_axes, b_axes = axes
            a_rest = tuple(d for i, d in enumerate(a) if i not in a_axes)
            b_rest = tuple(d for i, d in enumerate(b) if i not in b_axes)
        return a_rest + b_rest

    # --- transpose ---
    elif func_name in ("transpose", "T", "mT"):
        return a[:-2] + (a[-1], a[-2])

    # --- vecdot ---
    elif func_name == "vecdot":
        if axis is None and len(args) > 2:
            axis = args[2]
        if axis is None:
            axis = -1
        if b is None:
            return None
        a_axis = axis + len(a) if axis < 0 else axis
        b_axis = axis + len(b) if axis < 0 else axis
        a_rem = tuple(d for i, d in enumerate(a) if i != a_axis)
        b_rem = tuple(d for i, d in enumerate(b) if i != b_axis)
        return broadcast_shapes(a_rem, b_rem)
    else:
        return None


def reduce_shape(shape, axis, keepdims):
    """Reduce shape along given axis or axes (collapse dimensions)."""
    if shape is None:
        return None  # unknown shape

    # full reduction
    if axis is None:
        return (1,) * len(shape) if keepdims else ()

    # normalize to tuple
    if isinstance(axis, int):
        axes = (axis,)
    else:
        axes = tuple(axis)

    # normalize negative axes
    axes = tuple(a + len(shape) if a < 0 else a for a in axes)

    if keepdims:
        return tuple(d if i not in axes else 1 for i, d in enumerate(shape))
    else:
        return tuple(d for i, d in enumerate(shape) if i not in axes)


def slice_shape(shape, slices):
    """Infer shape after slicing."""
    if shape is None:
        return None
    result = []
    for dim, sl in zip(shape, slices, strict=False):
        if isinstance(sl, int):  # indexing removes the axis
            continue
        if isinstance(sl, slice):
            start = sl.start or 0
            stop = sl.stop if sl.stop is not None else dim
            step = sl.step or 1
            length = max(0, (stop - start + (step - 1)) // step)
            result.append(length)
        else:
            raise ValueError(f"Unsupported slice type: {sl}")
    result.extend(shape[len(slices) :])  # untouched trailing dims
    return tuple(result)


def elementwise(*args):
    """All args must broadcast elementwise."""
    if None in args:
        return None
    return broadcast_shapes(*args)


def cumulative_shape(x, axis=None, include_initial=False, out=None):
    if axis is None:
        if len(x) == 1:
            axis = 0
        else:
            raise ValueError("axis can only be None for 1D arrays")
    return tuple(d + 1 if (i == axis and include_initial) else d for i, d in enumerate(x))


# --- Function registry ---
REDUCTIONS = {  # ignore out arg
    func: cumulative_shape
    if func in {"cumulative_sum", "cumulative_prod"}
    else lambda x, axis=None, keepdims=False, out=None: reduce_shape(x, axis, keepdims)
    for func in reducers
    # any unknown function will default to elementwise
}


# --- AST Shape Inferencer ---
class ShapeInferencer(ast.NodeVisitor):
    def __init__(self, shapes):
        self.shapes = shapes

    def visit_Name(self, node):
        if node.id not in self.shapes:
            if node.id in ("nan", "inf"):  # non-finite float literals: scalars
                return ()
            raise ValueError(f"Unknown symbol: {node.id}")
        s = self.shapes[node.id]
        if isinstance(s, tuple):
            return s
        else:  # passed a scalar value
            return ()

    def visit_Attribute(self, node):
        obj_shape = self.visit(node.value)
        attr = node.attr
        if attr == "reshape":
            if node.args:
                shape_arg = node.args[-1]
                if isinstance(shape_arg, ast.Tuple):
                    return tuple(self._lookup_value(e) for e in shape_arg.elts)
            return ()
        elif attr in ("T", "mT"):
            return linalg_shape(attr, (obj_shape,), {})
        return None

    def visit_Call(self, node):  # noqa : C901
        # Extract full function name (support np.func, blosc2.func)
        func_name = getattr(node.func, "id", None)
        attr_name = getattr(node.func, "attr", None)
        module_name = getattr(getattr(node.func, "value", None), "id", None)

        # Handle namespaced calls like np.func or blosc2.func
        if module_name in ("np", "blosc2"):
            qualified_name = f"{module_name}.{attr_name}"
        else:
            qualified_name = attr_name or func_name

        base_name = qualified_name.split(".")[-1]

        # --- Recursive method-chain support ---
        obj_shape = None
        if isinstance(node.func, ast.Attribute) and module_name not in (
            "np",
            "blosc2",
        ):  # check if genuine method and not module func
            obj_shape = self.visit(node.func.value)

        args = [self.visit(arg) for arg in node.args]
        # If it's a method call, prepend the object shape
        if obj_shape is not None and attr_name == base_name:
            args.insert(0, obj_shape)

        # --- Parse keyword args ---
        kwargs = {}
        for kw in node.keywords:
            kwargs[kw.arg] = self._lookup_value(kw.value)

        # ------- handle linear algebra ---------------
        if base_name in linalg_funcs:
            return linalg_shape(base_name, args, kwargs)

        # ------- handle constructors ---------------
        if base_name in constructors:
            # shape kwarg directly provided
            if "shape" in kwargs:
                val = kwargs["shape"]
                return val if isinstance(val, tuple) else (val,)

            # ---- array constructors like zeros, ones, full, etc. ----
            elif base_name in (
                "zeros",
                "ones",
                "empty",
                "full",
                "full_like",
                "zeros_like",
                "empty_like",
                "ones_like",
                "nans",
            ):
                if node.args:
                    shape_arg = node.args[0]
                    if isinstance(shape_arg, ast.Tuple):
                        shape = tuple(self._lookup_value(e) for e in shape_arg.elts)
                    elif isinstance(shape_arg, ast.Constant):
                        shape = (shape_arg.value,)
                    else:
                        shape = self._lookup_value(shape_arg)
                        shape = shape if isinstance(shape, tuple) else (shape,)
                    return shape

            # ---- arange ----
            elif base_name == "arange":
                start = self._lookup_value(node.args[0]) if node.args else 0
                stop = self._lookup_value(node.args[1]) if len(node.args) > 1 else None
                step = self._lookup_value(node.args[2]) if len(node.args) > 2 else 1
                shape = self._lookup_value(node.args[4]) if len(node.args) > 4 else kwargs.get("shape")

                if shape is not None:
                    return shape if isinstance(shape, tuple) else (shape,)

                # Fallback to numeric difference if possible
                if stop is None:
                    stop, start = start, 0
                try:
                    NUM = max(math.ceil((stop - start) / step), 0)
                except Exception:
                    # symbolic or non-numeric: unknown 1D
                    return ((),)
                return (NUM,)

            # ---- linspace ----
            elif base_name == "linspace":
                num = self._lookup_value(node.args[2]) if len(node.args) > 2 else kwargs.get("num")
                shape = self._lookup_value(node.args[5]) if len(node.args) > 5 else kwargs.get("shape")
                if shape is not None:
                    return shape if isinstance(shape, tuple) else (shape,)
                if num is not None:
                    return (num,)
                raise ValueError("linspace requires either shape or num argument")

            elif base_name in {"frombuffer", "fromiter"}:
                count = kwargs.get("count")
                return (count,) if count else ()

            elif base_name == "eye":
                N = self._lookup_value(node.args[0])
                M = self._lookup_value(node.args[1]) if len(node.args) > 1 else kwargs.get("M")
                return (N, N) if M is None else (N, M)

            elif base_name == "reshape":
                if node.args:
                    shape_arg = node.args[-1]
                    if isinstance(shape_arg, ast.Tuple):
                        return tuple(self._lookup_value(e) for e in shape_arg.elts)
                return ()

            else:
                raise ValueError(f"Unrecognized constructor or missing shape argument for {func_name}")

        # --- Special-case .slice((slice(...), ...)) ---
        if attr_name == "slice":
            if not node.args:
                raise ValueError(".slice() requires an argument")
            slice_arg = node.args[0]
            if isinstance(slice_arg, ast.Tuple):
                slices = [self._eval_slice(s) for s in slice_arg.elts]
            else:
                slices = [self._eval_slice(slice_arg)]
            return slice_shape(obj_shape, slices)

        if base_name in REDUCTIONS:
            return REDUCTIONS[base_name](*args, **kwargs)

        shapes = [s for s in args if s is not None]
        if base_name not in elementwise_funcs:
            warnings.warn(
                f"Function shape parser not implemented for {base_name}.", UserWarning, stacklevel=2
            )
        # default to elementwise but print warning that function not defined explicitly
        return elementwise(*shapes) if shapes else ()

    def visit_Compare(self, node):
        shapes = [self.visit(node.left)] + [self.visit(c) for c in node.comparators]
        return elementwise(*shapes)

    def visit_Constant(self, node):
        return () if not hasattr(node.value, "shape") else node.value.shape

    def visit_Tuple(self, node):
        return tuple(self.visit(arg) for arg in node.elts)

    def visit_List(self, node):
        return self.visit_Tuple(node)

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        return elementwise(left, right)

    def visit_UnaryOp(self, node):
        return self.visit(node.operand)

    def _eval_slice(self, node):
        if isinstance(node, ast.Slice):
            return slice(
                node.lower.value if node.lower else None,
                node.upper.value if node.upper else None,
                node.step.value if node.step else None,
            )
        elif isinstance(node, ast.Call) and getattr(node.func, "id", None) == "slice":
            # handle explicit slice() constructor
            args = [a.value if isinstance(a, ast.Constant) else None for a in node.args]
            return slice(*args)
        elif isinstance(node, ast.Constant):
            return node.value
        else:
            raise ValueError(f"Unsupported slice expression: {ast.dump(node)}")

    def _lookup_value(self, node):  # noqa : C901
        """Look up a value in self.shapes if node is a variable name, else constant value."""
        # Name -> lookup in shapes mapping
        if isinstance(node, ast.Name):
            return self.shapes.get(node.id, None)

        # Constant -> return its value
        if isinstance(node, ast.Constant):
            return node.value

        # Tuple of constants / expressions
        if isinstance(node, ast.Tuple):
            vals = []
            for e in node.elts:
                v = self._lookup_value(e)
                vals.append(v)
            return tuple(vals)

        # Unary operations (e.g. -1)
        if isinstance(node, ast.UnaryOp):
            # handle negative constants like -1
            if isinstance(node.op, ast.USub):
                val = self._lookup_value(node.operand)
                if isinstance(val, int | float):
                    return -val
            # handle + (USub) if needed
            if isinstance(node.op, ast.UAdd):
                return self._lookup_value(node.operand)
            return None

        # Simple binary ops with constant operands (e.g. 1+2)
        if isinstance(node, ast.BinOp):
            left = self._lookup_value(node.left)
            right = self._lookup_value(node.right)
            if left is None or right is None:
                return None
            try:
                if isinstance(node.op, ast.Add):
                    return left + right
                if isinstance(node.op, ast.Sub):
                    return left - right
                if isinstance(node.op, ast.Mult):
                    return left * right
                if isinstance(node.op, ast.FloorDiv):
                    return left // right
                if isinstance(node.op, ast.Div):
                    return left / right
                if isinstance(node.op, ast.Mod):
                    return left % right
            except Exception:
                return None
            return None

        # fallback
        return None


# --- Public API ---
def infer_shape(expr, shapes):
    tree = ast.parse(expr, mode="eval")
    inferencer = ShapeInferencer(shapes)
    return inferencer.visit(tree.body)


class MyChunkRange:
    def __init__(self, start, stop, step=1, n=1):
        self.start = start
        self.stop = stop
        self.step = step
        self.n = n

    def __iter__(self):
        for k in range(math.ceil((self.stop - self.start) / self.step)):
            yield (self.start + k * self.step) // self.n


def slice_to_chunktuple(s, n):
    # Adapted from _slice_iter in ndindex.ChunkSize.as_subchunks.
    start, stop, step = s.start, s.stop, s.step
    if step < 0:
        temp = stop
        stop = start + 1
        start = temp + 1
        step = -step  # get positive steps
    if step > n:
        return MyChunkRange(start, stop, step, n)
    else:
        return range(start // n, ceiling(stop, n))


def get_selection(ctuple, ptuple, chunks):
    # we assume that at least one element of chunk intersects with the slice
    # (as a consequence of only looping over intersecting chunks)
    # ptuple is global slice, ctuple is chunk coords (in units of chunks)
    pselection = ()
    for i, s, csize in zip(ctuple, ptuple, chunks, strict=True):
        # we need to advance to first element within chunk that intersects with slice, not
        # necessarily the first element of chunk
        # i * csize = s.start + n*step + k, already added n+1 elements, k in [1, step]
        if s.step > 0:
            np1 = (i * csize - s.start + s.step - 1) // s.step  # gives (n + 1)
            # can have n = -1 if s.start > i * csize, but never < -1 since have to intersect with chunk
            pselection += (
                slice(
                    builtins.max(
                        s.start, s.start + np1 * s.step
                    ),  # start+(n+1)*step gives i*csize if k=step
                    builtins.min(csize * (i + 1), s.stop),
                    s.step,
                ),
            )
        else:
            # (i + 1) * csize = s.start + n*step + k, already added n+1 elements, k in [step+1, 0]
            np1 = ((i + 1) * csize - s.start + s.step) // s.step  # gives (n + 1)
            # can have n = -1 if s.start < (i + 1) * csize, but never < -1 since have to intersect with chunk
            pselection += (
                slice(
                    builtins.min(s.start, s.start + np1 * s.step),  # start+n*step gives (i+1)*csize if k=0
                    builtins.max(csize * i - 1, s.stop),  # want to include csize * i
                    s.step,
                ),
            )

    # selection relative to coordinates of out (necessarily out_step = 1 as we work through out chunk-by-chunk of self)
    # when added n + 1 elements
    # ps.start = pt.start + step * (n+1) => n = (ps.start - pt.start - sign) // step
    # hence, out_start = n + 1
    # ps.stop = pt.start + step * (out_stop - 1) + k,  k in [step, -1] or [1, step]
    # => out_stop = (ps.stop - pt.start - sign) // step + 1
    out_pselection = ()
    for ps, pt in zip(pselection, ptuple, strict=True):
        sign_ = np.sign(pt.step)
        n = (ps.start - pt.start - sign_) // pt.step
        out_start = n + 1
        # ps.stop always positive except for case where get full array (it is then -1 since desire 0th element)
        out_stop = None if ps.stop == -1 else (ps.stop - pt.start - sign_) // pt.step + 1
        out_pselection += (
            slice(
                out_start,
                out_stop,
                1,
            ),
        )

    loc_selection = tuple(  # is s.stop is None, get whole chunk so s.start - 0
        slice(0, s.stop - s.start, s.step)
        if s.step > 0
        else slice(s.start if s.stop == -1 else s.start - s.stop, None, s.step)
        for s in pselection
    )  # local coords of loaded part of chunk

    return out_pselection, pselection, loc_selection


def get_local_slice(prior_selection, post_selection, chunk_bounds):
    chunk_begin, chunk_end = chunk_bounds
    # +1 for negative steps as have to include start (exclude stop)
    locbegin = np.hstack(
        (
            [s.start if s.step > 0 else s.stop + 1 for s in prior_selection],
            chunk_begin,
            [s.start if s.step > 0 else s.stop + 1 for s in post_selection],
        ),
        casting="unsafe",
        dtype="int64",
    )
    locend = np.hstack(
        (
            [s.stop if s.step > 0 else s.start + 1 for s in prior_selection],
            chunk_end,
            [s.stop if s.step > 0 else s.start + 1 for s in post_selection],
        ),
        casting="unsafe",
        dtype="int64",
    )
    return locbegin, locend


def sliced_chunk_iter(chunks, idx, shape, axis=None, nchunk=False):
    """
    If nchunk is True, retrun at iterator over the number of the chunk.
    """
    ratio = np.ceil(np.asarray(shape) / np.asarray(chunks)).astype(np.int64)
    idx = ndindex.ndindex(idx).expand(shape)
    if axis is not None:
        idx = tuple(a for i, a in enumerate(idx.args) if i != axis) + (idx.args[axis],)
        chunks_ = tuple(a for i, a in enumerate(chunks) if i != axis) + (chunks[axis],)
    else:
        chunks_ = chunks
    idx_iter = iter(idx)  # iterate over tuple of slices in order
    chunk_iter = iter(chunks_)  # iterate over chunk_shape in order

    iters = []
    while True:
        try:
            i = next(idx_iter)  # slice along axis
            n = next(chunk_iter)  # chunklen along dimension
        except StopIteration:
            break
        if not isinstance(i, ndindex.Slice):
            raise ValueError("Only slices may be used with axis arg")

        def _slice_iter(s, n):
            a, N, m = s.args
            if m > n:
                yield from ((a + k * m) // n for k in range(ceiling(N - a, m)))
            else:
                yield from range(a // n, ceiling(N, n))

        iters.append(_slice_iter(i, n))

    def _indices(iters):
        my_list = [ndindex.Slice(None, None)] * len(chunks)
        for p in product(*iters):
            # p increments over arg axis first before other axes
            # p = (...., -1, axis)
            if axis is None:
                my_list = [
                    ndindex.Slice(cs * ci, min(cs * (ci + 1), n), 1)
                    for n, cs, ci in zip(shape, chunks, p, strict=True)
                ]
            else:
                my_list[:axis] = [
                    ndindex.Slice(cs * ci, min(cs * (ci + 1), n), 1)
                    for n, cs, ci in zip(shape[:axis], chunks[:axis], p[:axis], strict=True)
                ]
                n, cs, ci = shape[axis], chunks[axis], p[-1]
                my_list[axis] = ndindex.Slice(cs * ci, min(cs * (ci + 1), n), 1)
                my_list[axis + 1 :] = [
                    ndindex.Slice(cs * ci, min(cs * (ci + 1), n), 1)
                    for n, cs, ci in zip(shape[axis + 1 :], chunks[axis + 1 :], p[axis:-1], strict=True)
                ]
            if nchunk:
                yield builtins.sum(
                    c.start // chunks[i] * np.prod(ratio[i + 1 :]) for i, c in enumerate(my_list)
                )
            else:
                yield ndindex.Tuple(*my_list)

    yield from _indices(iters)


def get_intersecting_chunks(idx, shape, chunks, axis=None):
    if len(chunks) != len(shape):
        raise ValueError("chunks must be same length as shape!")
    if 0 in chunks:  # chunk is whole array so just return full tuple to do loop once
        return (ndindex.ndindex(...).expand(shape),)
    chunk_size = ndindex.ChunkSize(chunks)
    if axis is None:
        return chunk_size.as_subchunks(idx, shape)  # if _slice is (), returns all chunks

    # special algorithm to iterate over axis first (adapted from ndindex source)
    return sliced_chunk_iter(chunks, idx, shape, axis)


def get_chunks_idx(shape, chunks):
    chunks_idx = tuple(math.ceil(s / c) for s, c in zip(shape, chunks, strict=True))
    nchunks = math.prod(chunks_idx)
    return chunks_idx, nchunks


def process_key(key, shape):
    key = ndindex.ndindex(key).expand(shape).raw
    mask = tuple(
        isinstance(k, int) for k in key
    )  # mask to track dummy dims introduced by int -> slice(k, k+1)
    key = tuple(slice(k, k + 1, None) if isinstance(k, int) else k for k in key)  # key is slice, None, int
    return key, mask


def is_inside_ne_evaluate() -> bool:
    """
    Whether the current code is being executed from an ne_evaluate call
    """
    # Get the current call stack
    stack = inspect.stack()
    return builtins.any(frame_info.function in {"ne_evaluate"} for frame_info in stack)


def incomplete_lazyfunc(func) -> None:
    """Decorator for lazy functions with incomplete numexpr/miniexpr coverage.

    This function will force eager execution when called from ne_evaluate.

    Returns
    -------
    out: None

    Examples
    --------
    .. code-block:: python

        @incomplete_lazyfunc()
        def filler(inputs_tuple, output, offset):
            output[:] = inputs_tuple[0] - inputs_tuple[1]

    """

    def wrapper(*args, **kwargs):
        if is_inside_ne_evaluate():  # haven't been able to use miniexpr so use numpy
            return safe_numpy_globals[func.__name__](*args, **kwargs)
        return func(*args, **kwargs)

    return wrapper


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


def get_chunk_operands(operands, cslice, chunk_operands, shape):
    # Get the starts and stops for the slice
    cslice_shape = tuple(s.stop - s.start for s in cslice)
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
