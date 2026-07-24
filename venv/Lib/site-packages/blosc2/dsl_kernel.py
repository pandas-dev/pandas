#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################

from __future__ import annotations

import ast
import contextlib
import inspect
import os
import textwrap
import tokenize
from io import StringIO
from typing import ClassVar

_PRINT_DSL_KERNEL = os.environ.get("PRINT_DSL_KERNEL", "").strip().lower()
_PRINT_DSL_KERNEL = _PRINT_DSL_KERNEL not in ("", "0", "false", "no", "off")
_DSL_USAGE_DOC_URL = "https://github.com/Blosc/python-blosc2/blob/main/doc/reference/dsl_syntax.md"


class DSLSyntaxError(ValueError):
    """Raised when a @dsl_kernel function uses unsupported DSL syntax."""


def _normalize_miniexpr_scalar(value):
    # NumPy scalar-like values expose .item(); plain Python scalars do not.
    # Do not call .item() on non-scalar arrays; for Blosc2 arrays this can be expensive.
    if hasattr(value, "shape") and value.shape != ():
        raise TypeError("Unsupported scalar type for miniexpr specialization")
    if hasattr(value, "item") and callable(value.item):
        with contextlib.suppress(Exception):
            value = value.item()
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int | float | str):
        return value
    raise TypeError("Unsupported scalar type for miniexpr specialization")


def _line_starts(text: str) -> list[int]:
    starts = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            starts.append(i + 1)
    return starts


def _to_abs(line_starts: list[int], line: int, col: int) -> int:
    return line_starts[line - 1] + col


def _find_def_signature_span(text: str):
    tokens = list(tokenize.generate_tokens(StringIO(text).readline))
    for i, tok in enumerate(tokens):
        if tok.type != tokenize.NAME or tok.string != "def":
            continue
        lparen = None
        rparen = None
        colon = None
        depth = 0
        for j in range(i + 1, len(tokens)):
            t = tokens[j]
            if lparen is None:
                if t.type == tokenize.OP and t.string == "(":
                    lparen = t
                    depth = 1
                continue
            if t.type == tokenize.OP and t.string == "(":
                depth += 1
                continue
            if t.type == tokenize.OP and t.string == ")":
                depth -= 1
                if depth == 0:
                    rparen = t
                    continue
            if rparen is not None and t.type == tokenize.OP and t.string == ":":
                colon = t
                break
        if lparen is not None and rparen is not None:
            return lparen, rparen, colon
    return None, None, None


def _remove_scalar_params_preserving_source(text: str, scalar_replacements: dict[str, int | float]):
    if not scalar_replacements:
        return text, 0

    lparen, rparen, colon = _find_def_signature_span(text)
    if lparen is None or rparen is None:
        return text, 0

    try:
        tree = ast.parse(text)
    except Exception:
        return text, 0

    func = next((n for n in tree.body if isinstance(n, ast.FunctionDef)), None)
    if func is None:
        return text, 0

    kept = [a.arg for a in (func.args.posonlyargs + func.args.args) if a.arg not in scalar_replacements]
    line_starts = _line_starts(text)
    pstart = _to_abs(line_starts, lparen.end[0], lparen.end[1])
    pend = _to_abs(line_starts, rparen.start[0], rparen.start[1])
    updated = f"{text[:pstart]}{', '.join(kept)}{text[pend:]}"
    body_start = 0
    if colon is not None:
        # Signature shrink can move ':' to an earlier column, so recompute
        # on the rewritten text to avoid skipping first-line body tokens.
        _, _, updated_colon = _find_def_signature_span(updated)
        if updated_colon is not None:
            body_start = _to_abs(_line_starts(updated), updated_colon.end[0], updated_colon.end[1])
    return updated, body_start


def _replace_scalar_names_preserving_source(
    text: str, scalar_replacements: dict[str, int | float], body_start: int
):
    if not scalar_replacements:
        return text

    line_starts = _line_starts(text)
    tokens = list(tokenize.generate_tokens(StringIO(text).readline))
    significant = {
        tokenize.NAME,
        tokenize.NUMBER,
        tokenize.STRING,
        tokenize.OP,
        tokenize.INDENT,
        tokenize.DEDENT,
    }
    assign_ops = {"=", "+=", "-=", "*=", "/=", "//=", "%=", "&=", "|=", "^=", "<<=", ">>=", ":="}
    edits = []
    for i, tok in enumerate(tokens):
        if tok.type != tokenize.NAME or tok.string not in scalar_replacements:
            continue
        start_abs = _to_abs(line_starts, tok.start[0], tok.start[1])
        if start_abs < body_start:
            continue

        prev_sig = None
        for j in range(i - 1, -1, -1):
            if tokens[j].type in significant:
                prev_sig = tokens[j]
                break
        if prev_sig is not None and prev_sig.type == tokenize.OP and prev_sig.string == ".":
            continue

        next_sig = None
        for j in range(i + 1, len(tokens)):
            if tokens[j].type in significant:
                next_sig = tokens[j]
                break
        if next_sig is not None and next_sig.type == tokenize.OP and next_sig.string in assign_ops:
            continue

        end_abs = _to_abs(line_starts, tok.end[0], tok.end[1])
        edits.append((start_abs, end_abs, repr(scalar_replacements[tok.string])))

    if not edits:
        return text

    out = text
    for start, end, repl in sorted(edits, key=lambda e: e[0], reverse=True):
        out = f"{out[:start]}{repl}{out[end:]}"
    return out


def _fold_numeric_cast_calls_preserving_source(text: str, body_start: int):
    """Fold float(<number>) and int(<number>) calls into literals.

    miniexpr parses DSL function calls in a restricted way, and scalar specialization can
    produce calls like float(200) that fail to parse.  Fold those into literals while
    preserving source formatting/comments elsewhere.
    """
    try:
        tree = ast.parse(text)
    except Exception:
        return text

    line_starts = _line_starts(text)
    edits = []

    def _numeric_literal_value(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int | float | bool):
            return node.value
        if (
            isinstance(node, ast.UnaryOp)
            and isinstance(node.op, ast.UAdd | ast.USub)
            and isinstance(node.operand, ast.Constant)
            and isinstance(node.operand.value, int | float | bool)
        ):
            value = node.operand.value
            return value if isinstance(node.op, ast.UAdd) else -value
        return None

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if node.keywords or len(node.args) != 1:
            continue
        if not isinstance(node.func, ast.Name) or node.func.id not in {"float", "int"}:
            continue

        arg = node.args[0]
        value = _numeric_literal_value(arg)
        if value is None:
            continue

        start_abs = _to_abs(line_starts, node.lineno, node.col_offset)
        if start_abs < body_start:
            continue
        end_abs = _to_abs(line_starts, node.end_lineno, node.end_col_offset)

        if node.func.id == "float":
            repl = repr(float(value))
        else:
            repl = repr(int(value))
        edits.append((start_abs, end_abs, repl))

    if not edits:
        return text

    out = text
    for start, end, repl in sorted(edits, key=lambda e: e[0], reverse=True):
        out = f"{out[:start]}{repl}{out[end:]}"
    return out


def specialize_miniexpr_inputs(expr_string: str, operands: dict):
    """Inline scalar operands as constants for miniexpr compilation."""
    scalar_replacements = {}
    array_operands = {}
    for name, value in operands.items():
        if hasattr(value, "shape") and value.shape == ():
            scalar_replacements[name] = _normalize_miniexpr_scalar(value[()])
            continue
        if isinstance(value, int | float | bool | str) or (
            not hasattr(value, "shape") and hasattr(value, "item") and callable(value.item)
        ):
            try:
                scalar_replacements[name] = _normalize_miniexpr_scalar(value)
                continue
            except TypeError:
                pass
        array_operands[name] = value

    if not scalar_replacements:
        return expr_string, operands

    rewritten, body_start = _remove_scalar_params_preserving_source(expr_string, scalar_replacements)
    rewritten = _replace_scalar_names_preserving_source(rewritten, scalar_replacements, body_start)
    rewritten = _fold_numeric_cast_calls_preserving_source(rewritten, body_start)
    return rewritten, array_operands


def specialize_dsl_miniexpr_inputs(expr_string: str, operands: dict):
    """Backward-compatible alias for DSL-specific callers."""
    return specialize_miniexpr_inputs(expr_string, operands)


class DSLValidator:
    _binop_map: ClassVar[dict[type[ast.operator], str]] = {
        ast.Add: "+",
        ast.Sub: "-",
        ast.Mult: "*",
        ast.Div: "/",
        ast.FloorDiv: "//",
        ast.Mod: "%",
        ast.Pow: "**",
        ast.BitAnd: "&",
        ast.BitOr: "|",
        ast.BitXor: "^",
        ast.LShift: "<<",
        ast.RShift: ">>",
    }
    _cmp_map: ClassVar[dict[type[ast.cmpop], str]] = {
        ast.Eq: "==",
        ast.NotEq: "!=",
        ast.Lt: "<",
        ast.LtE: "<=",
        ast.Gt: ">",
        ast.GtE: ">=",
    }

    def __init__(self, source: str, line_base: int = 0, input_names: list[str] | None = None):
        self._source = source
        self._line_base = line_base
        self._inputs = set(input_names or ())

    def validate(self, func_node: ast.FunctionDef):
        self._args(func_node)
        if not func_node.body:
            self._err(func_node, "DSL kernel must have a body")
        self._one_per_line(func_node.body)
        for stmt in func_node.body:
            self._stmt(stmt)

    def _one_per_line(self, body: list[ast.stmt]):
        # G1: miniexpr parses one statement per line; `;`-joined siblings share a lineno.
        prev = None
        for stmt in body:
            if prev is not None and stmt.lineno == prev:
                self._err(
                    stmt,
                    "Only one statement per line is supported in DSL kernels; "
                    "split ';'-joined statements onto separate lines",
                )
            prev = stmt.lineno

    def _err(self, node: ast.AST, msg: str, *, line: int | None = None, col: int | None = None):
        if line is None:
            line = getattr(node, "lineno", 0)
        if col is None:
            col = getattr(node, "col_offset", 0) + 1
        line -= self._line_base
        location = f"{msg} at line {line}, column {col}"
        dump = self._format_source_with_pointer(line, col)
        raise DSLSyntaxError(f"{location}\n\nDSL kernel source:\n{dump}\n\nSee: {_DSL_USAGE_DOC_URL}")

    def _format_source_with_pointer(self, line: int, col: int) -> str:
        lines = self._source.splitlines()
        if not lines:
            return "<empty>"
        width = len(str(len(lines)))
        out = []
        for lineno, text in enumerate(lines, start=1):
            out.append(f"{lineno:>{width}} | {text}")
            if lineno == line:
                pointer = " " * max(col - 1, 0)
                out.append(f"{' ' * width} | {pointer}^")
        return "\n".join(out)

    def _args(self, func_node: ast.FunctionDef):
        args = func_node.args
        if args.vararg or args.kwarg or args.kwonlyargs:
            self._err(args, "DSL kernel does not support *args/**kwargs/kwonly args")
        if args.defaults or args.kw_defaults:
            self._err(args, "DSL kernel does not support default arguments")

    def _check_input_assign(self, target: ast.Name):
        # G2: miniexpr forbids reassigning an input parameter (inputs alias operand buffers).
        if target.id in self._inputs:
            self._err(
                target,
                f"Cannot assign to input parameter '{target.id}'; "
                f"copy it into a local variable first (e.g. 'tmp = {target.id}')",
            )

    def _stmt(self, node: ast.stmt):  # noqa: C901
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                self._err(node, "Only simple assignments are supported in DSL kernels")
            self._check_input_assign(node.targets[0])
            self._expr(node.value)
            return
        if isinstance(node, ast.AugAssign):
            if not isinstance(node.target, ast.Name):
                self._err(node, "Only simple augmented assignments are supported")
            self._check_input_assign(node.target)
            self._binop(node.op)
            self._expr(node.value)
            return
        if isinstance(node, ast.Return):
            if node.value is None:
                self._err(node, "DSL kernel return must have a value")
            self._expr(node.value)
            return
        if isinstance(node, ast.Expr):
            self._expr(node.value)
            return
        if isinstance(node, ast.If):
            self._expr(node.test)
            if not node.body:
                self._err(node, "Empty if blocks are not supported in DSL kernels")
            self._one_per_line(node.body)
            self._one_per_line(node.orelse)
            for stmt in node.body:
                self._stmt(stmt)
            for stmt in node.orelse:
                self._stmt(stmt)
            return
        if isinstance(node, ast.For):
            if node.orelse:
                self._err(node, "for/else is not supported in DSL kernels")
            if not isinstance(node.target, ast.Name):
                self._err(node, "DSL for-loop target must be a simple name")
            if not isinstance(node.iter, ast.Call):
                self._err(node, "DSL for-loop must iterate over range()")
            func_name = self._call_name(node.iter.func)
            if func_name != "range":
                self._err(node, "DSL for-loop must iterate over range()")
            if node.iter.keywords or not (1 <= len(node.iter.args) <= 3):
                self._err(node, "DSL range() must take 1 to 3 positional arguments")
            for arg in node.iter.args:
                self._expr(arg)
            if not node.body:
                self._err(node, "Empty for-loop bodies are not supported in DSL kernels")
            self._one_per_line(node.body)
            for stmt in node.body:
                self._stmt(stmt)
            return
        if isinstance(node, ast.While):
            if node.orelse:
                self._err(node, "while/else is not supported in DSL kernels")
            self._expr(node.test)
            if not node.body:
                self._err(node, "Empty while-loop bodies are not supported in DSL kernels")
            self._one_per_line(node.body)
            for stmt in node.body:
                self._stmt(stmt)
            return
        if isinstance(node, ast.Break | ast.Continue):
            return
        self._err(node, f"Unsupported DSL statement: {type(node).__name__}")

    def _expr(self, node: ast.AST):  # noqa: C901
        if isinstance(node, ast.Name):
            return
        if isinstance(node, ast.Constant):
            val = node.value
            if isinstance(val, bool | int | float | str):
                return
            self._err(node, "Unsupported constant in DSL expression")
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.UAdd | ast.USub | ast.Not):
                self._expr(node.operand)
                return
            self._err(node, "Unsupported unary operator in DSL expression")
        if isinstance(node, ast.BinOp):
            self._binop(node.op)
            self._expr(node.left)
            self._expr(node.right)
            return
        if isinstance(node, ast.BoolOp):
            for value in node.values:
                self._expr(value)
            return
        if isinstance(node, ast.Compare):
            if len(node.ops) != 1 or len(node.comparators) != 1:
                self._err(node, "Chained comparisons are not supported in DSL")
            self._cmpop(node.ops[0])
            self._expr(node.left)
            self._expr(node.comparators[0])
            return
        if isinstance(node, ast.Call):
            self._call_name(node.func)
            if node.keywords:
                self._err(node, "Keyword arguments are not supported in DSL calls")
            for arg in node.args:
                self._expr(arg)
            return
        if isinstance(node, ast.IfExp):
            seg = ast.get_source_segment(self._source, node)
            col = getattr(node, "col_offset", 0) + 1
            if seg is not None:
                rel = seg.find(" if ")
                if rel >= 0:
                    col += rel + 1
            self._err(
                node,
                "Ternary expressions are not supported in DSL; use where(cond, a, b)",
                col=col,
            )
        self._err(node, f"Unsupported DSL expression: {type(node).__name__}")

    def _call_name(self, node: ast.AST) -> str:
        if isinstance(node, ast.Name):
            return node.id
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in {"np", "numpy", "math"}
        ):
            return node.attr
        self._err(node, "Unsupported call target in DSL")
        raise AssertionError("unreachable")

    def _binop(self, op: ast.operator):
        for k in self._binop_map:
            if isinstance(op, k):
                return
        self._err(op, "Unsupported binary operator in DSL")

    def _cmpop(self, op: ast.cmpop):
        for k in self._cmp_map:
            if isinstance(op, k):
                return
        self._err(op, "Unsupported comparison in DSL")


class DSLKernel:
    """Wrap a Python function and optionally extract a miniexpr DSL kernel from it."""

    def __init__(self, func):
        self.func = func
        self.__name__ = getattr(func, "__name__", self.__class__.__name__)
        self.__qualname__ = getattr(func, "__qualname__", self.__name__)
        self.__doc__ = getattr(func, "__doc__", None)
        try:
            sig = inspect.signature(func)
        except (TypeError, ValueError):
            sig = None
        self._sig = sig
        self._sig_has_varargs = False
        self._sig_npositional = None
        self._legacy_udf_signature = False
        if sig is not None:
            params = list(sig.parameters.values())
            positional_params = [p for p in params if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
            self._sig_has_varargs = any(p.kind == p.VAR_POSITIONAL for p in params)
            self._sig_npositional = len(positional_params)
            # Preserve support for classic lazyudf signature: (inputs_tuple, output, offset)
            if not self._sig_has_varargs and len(positional_params) == 3:
                p2 = positional_params[1].name.lower()
                p3 = positional_params[2].name.lower()
                self._legacy_udf_signature = p2 in {"output", "out"} and p3 == "offset"
        self.dsl_source = None
        self.input_names = None
        self.dsl_error = None
        try:
            dsl_source, input_names = self._extract_dsl(func)
        except DSLSyntaxError as e:
            # Preserve extracted source/signature for diagnostics even when DSL validation fails.
            try:
                dsl_source, input_names = self._extract_dsl(func, validate=False)
            except Exception:
                dsl_source = None
                input_names = None
            self.dsl_error = e
        except Exception:
            dsl_source = None
            input_names = None
        self.dsl_source = dsl_source
        self.input_names = input_names

    def _extract_dsl(self, func, validate: bool = True):
        source = inspect.getsource(func)
        source = textwrap.dedent(source)
        tree = ast.parse(source)
        func_node = None
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == func.__name__:
                func_node = node
                break
        if func_node is None:
            for node in tree.body:
                if isinstance(node, ast.FunctionDef):
                    func_node = node
                    break
        if func_node is None:
            raise ValueError("No function definition found for DSL extraction")

        dsl_source = self._slice_function_source(source, func_node)
        dsl_tree = ast.parse(dsl_source)
        dsl_func = next((node for node in dsl_tree.body if isinstance(node, ast.FunctionDef)), None)
        if dsl_func is None:
            raise ValueError("No function definition found in sliced DSL source")
        input_names = self._input_names_from_signature(dsl_func)
        if validate:
            DSLValidator(dsl_source, input_names=input_names).validate(dsl_func)
        if _PRINT_DSL_KERNEL:
            func_name = getattr(func, "__name__", "<dsl_kernel>")
            print(f"[DSLKernel:{func_name}] dsl_source (full):")
            print(dsl_source)
        return dsl_source, input_names

    @staticmethod
    def _slice_function_source(source: str, func_node: ast.FunctionDef) -> str:
        lines = source.splitlines()
        start = func_node.lineno - 1
        end_lineno = getattr(func_node, "end_lineno", None)
        if end_lineno is None:
            end = len(lines)
        else:
            end = end_lineno
        return "\n".join(lines[start:end])

    @staticmethod
    def _input_names_from_signature(func_node: ast.FunctionDef) -> list[str]:
        args = func_node.args
        if args.vararg or args.kwarg or args.kwonlyargs:
            raise ValueError("DSL kernel does not support *args/**kwargs/kwonly args")
        if args.defaults or args.kw_defaults:
            raise ValueError("DSL kernel does not support default arguments")
        return [a.arg for a in (args.posonlyargs + args.args)]

    def __call__(self, inputs_tuple, output, offset=None):
        if self.dsl_error is not None:
            raise self.dsl_error
        if self._legacy_udf_signature:
            return self.func(inputs_tuple, output, offset)

        n_inputs = len(inputs_tuple)
        if self._sig is not None and (
            self._sig_npositional in (n_inputs, n_inputs + 1) or self._sig_has_varargs
        ):
            if self._sig_npositional == n_inputs + 1:
                result = self.func(*inputs_tuple, offset)
            else:
                result = self.func(*inputs_tuple)
            output[...] = result
            return None

        try:
            return self.func(inputs_tuple, output, offset)
        except TypeError:
            result = self.func(*inputs_tuple)
            output[...] = result
            return None


def dsl_kernel(func):
    """Decorator to wrap a function in a DSLKernel."""

    return DSLKernel(func)


def validate_dsl(func):
    """Validate a DSL kernel function without executing it.

    Parameters
    ----------
    func
        A Python callable or :class:`DSLKernel`.

    Returns
    -------
    dict
        A dictionary with:
        - ``valid`` (bool): whether the DSL is valid
        - ``dsl_source`` (str | None): extracted DSL source when valid
        - ``input_names`` (list[str] | None): input signature names when valid
        - ``error`` (str | None): user-facing error message when invalid

    Examples
    --------
    >>> import blosc2
    >>> @blosc2.dsl_kernel
    ... def k(a, b):
    ...     return a * a + b * b
    >>> info = blosc2.validate_dsl(k)
    >>> info["valid"]
    True
    >>> info["input_names"]
    ['a', 'b']

    An unsupported construct is reported instead of raising:

    >>> @blosc2.dsl_kernel
    ... def bad(a):
    ...     return 1 if a > 0 else 0
    >>> blosc2.validate_dsl(bad)["valid"]
    False

    See Also
    --------
    validate_dsl_jit : Additionally probe whether the kernel JIT-compiles (vs falls
        back to the interpreter) for given operand/output dtypes.
    """

    kernel = func if isinstance(func, DSLKernel) else DSLKernel(func)
    err = kernel.dsl_error
    return {
        "valid": err is None,
        "dsl_source": kernel.dsl_source,
        "input_names": kernel.input_names,
        "error": None if err is None else str(err),
    }


def validate_dsl_jit(func, operands, out_dtype, *, shape=(64,), chunks=None, blocks=None):
    """Report whether a DSL kernel JIT-compiles (vs. interpreter fallback).

    Compiles a tiny probe of the kernel for the given operand/output dtypes and
    queries miniexpr for whether it produced a runtime JIT kernel.  The kernel is
    *not* run on real data, but compilation is dtype-specialized, so the operand and
    output dtypes must be supplied.

    Parameters
    ----------
    func
        A Python callable or :class:`DSLKernel`.
    operands
        One spec per kernel input.  A NumPy dtype (e.g. ``np.float64``, ``"f8"``)
        marks an *array* operand; a Python ``int``/``float``/``bool`` marks a *scalar*
        parameter inlined as a constant (as happens at run time).  Either a dict
        mapping input name -> spec, or a sequence matched positionally to the inputs.
    out_dtype
        Output dtype to specialize the codegen for.
    shape
        Probe shape; ``len(shape)`` sets the kernel dimensionality (kernels using
        ``_i1`` etc. need a matching ndim).  Defaults to ``(64,)``.

    Returns
    -------
    dict
        - ``valid`` (bool): DSL syntax is valid
        - ``jit`` (bool): a runtime JIT kernel was produced (vs interpreter fallback)
        - ``compiled`` (bool): miniexpr accepted the kernel
        - ``status`` (str | None): miniexpr compile-status name
        - ``error`` (str | None): syntax error message when ``valid`` is False

    Examples
    --------
    >>> import blosc2
    >>> import numpy as np
    >>> @blosc2.dsl_kernel
    ... def k(a, b):
    ...     return a * a + b * b
    >>> status = blosc2.validate_dsl_jit(k, [np.float64, np.float64], np.float64)
    >>> status["jit"]
    True
    >>> status["status"]
    'ME_COMPILE_SUCCESS'

    Pass array operands as dtypes and scalar parameters as values (the scalar is
    inlined as a constant, as at run time):

    >>> @blosc2.dsl_kernel
    ... def scaled(a, n):
    ...     return a * float(n)
    >>> blosc2.validate_dsl_jit(scaled, {"a": np.float64, "n": 3}, np.float64)["jit"]
    True

    See Also
    --------
    validate_dsl : Static DSL-syntax check only (no compilation; no dtypes needed).
    """
    import numpy as np

    import blosc2

    kernel = func if isinstance(func, DSLKernel) else DSLKernel(func)
    if kernel.dsl_error is not None:
        return {
            "valid": False,
            "jit": False,
            "compiled": False,
            "status": None,
            "error": str(kernel.dsl_error),
        }

    names = list(kernel.input_names or [])
    if isinstance(operands, dict):
        specs = dict(operands)
    else:
        operands = list(operands)
        if len(operands) != len(names):
            raise ValueError(f"expected {len(names)} operand specs for inputs {names}, got {len(operands)}")
        specs = dict(zip(names, operands, strict=True))

    scalars: dict = {}
    arrays: dict = {}
    for name in names:
        if name not in specs:
            raise ValueError(f"missing operand spec for input '{name}'")
        spec = specs[name]
        if isinstance(spec, bool | int | float):
            scalars[name] = spec
        else:
            arrays[name] = np.dtype(spec)

    source = kernel.dsl_source
    if scalars:
        source, _ = specialize_miniexpr_inputs(source, scalars)
    if not arrays:
        raise ValueError(
            "validate_dsl_jit needs at least one array operand (pass a dtype spec); "
            "all-scalar probes are not supported"
        )

    out = blosc2.empty(shape, dtype=np.dtype(out_dtype), chunks=chunks, blocks=blocks)
    probe_inputs = {name: np.empty(1, dtype=dt) for name, dt in arrays.items()}
    status = out._dsl_jit_status(source, probe_inputs)
    status["valid"] = True
    status["error"] = None
    return status


def kernel_from_source(source: str, name: str | None = None) -> DSLKernel:
    """Reconstruct a :class:`DSLKernel` from its stored source text.

    Executes *source* in a restricted namespace (builtins minus ``__import__``,
    plus ``np`` and ``blosc2``), extracts the defined function, and wraps it in
    a :class:`DSLKernel`.  This is the inverse of persisting
    :attr:`DSLKernel.dsl_source` and is shared by the persisted-``LazyUDF``
    decoder and the CTable DSL-column loaders.

    Parameters
    ----------
    source
        Complete, standalone function-definition source (as produced by
        :attr:`DSLKernel.dsl_source` or :func:`inspect.getsource`).
    name
        Name of the function to extract from *source*.  When omitted, the
        single top-level function definition in *source* is used.
    """
    import builtins
    import linecache

    import numpy as np

    import blosc2

    tree = ast.parse(source)
    func_defs = [n for n in tree.body if isinstance(n, ast.FunctionDef)]
    if name is None:
        if len(func_defs) != 1:
            raise ValueError(
                "kernel_from_source requires an explicit 'name' when 'source' does not "
                "contain exactly one top-level function definition"
            )
        name = func_defs[0].name

    # Security: reject sources that include top-level statements with
    # side effects (assignments, calls, loops, etc.) beyond the function
    # definition itself.  Imports and docstrings are permitted.
    allowed_nodes = (ast.FunctionDef, ast.Import, ast.ImportFrom)
    for node in tree.body:
        if isinstance(node, allowed_nodes):
            continue
        # A bare string literal at module level is a docstring — safe.
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            continue
        raise ValueError(
            f"kernel_from_source source must contain only a function definition "
            f"(optionally preceded by imports), but found: {ast.dump(node)[:120]}"
        )

    local_ns: dict = {}
    filename = f"<{name}>"
    safe_globals = {
        "__builtins__": {k: v for k, v in builtins.__dict__.items() if k != "__import__"},
        "np": np,
        "blosc2": blosc2,
    }
    linecache.cache[filename] = (len(source), None, source.splitlines(True), filename)
    exec(compile(source, filename, "exec"), safe_globals, local_ns)
    try:
        func = local_ns[name]
    except KeyError:
        raise ValueError(f"kernel_from_source: source did not define function {name!r}") from None
    if not isinstance(func, DSLKernel):
        func = DSLKernel(func)
    return func


class DSLBuilder:
    _binop_map: ClassVar[dict[type[ast.operator], str]] = {
        ast.Add: "+",
        ast.Sub: "-",
        ast.Mult: "*",
        ast.Div: "/",
        ast.FloorDiv: "//",
        ast.Mod: "%",
        ast.Pow: "**",
        ast.BitAnd: "&",
        ast.BitOr: "|",
        ast.BitXor: "^",
        ast.LShift: "<<",
        ast.RShift: ">>",
    }

    _cmp_map: ClassVar[dict[type[ast.cmpop], str]] = {
        ast.Eq: "==",
        ast.NotEq: "!=",
        ast.Lt: "<",
        ast.LtE: "<=",
        ast.Gt: ">",
        ast.GtE: ">=",
    }

    def __init__(self):
        self._lines = []

    def build(self, func_node: ast.FunctionDef):
        input_names = self._args(func_node.args)
        self._emit(f"def {func_node.name}({', '.join(input_names)}):", 0)
        if not func_node.body:
            raise ValueError("DSL kernel must have a body")
        for stmt in func_node.body:
            self._stmt(stmt, 4)
        return "\n".join(self._lines), input_names

    def _emit(self, line: str, indent: int):
        self._lines.append(" " * indent + line)

    def _args(self, args: ast.arguments):
        if args.vararg or args.kwarg or args.kwonlyargs:
            raise ValueError("DSL kernel does not support *args/**kwargs/kwonly args")
        if args.defaults or args.kw_defaults:
            raise ValueError("DSL kernel does not support default arguments")
        names = [a.arg for a in (args.posonlyargs + args.args)]
        if not names:
            raise ValueError("DSL kernel must accept at least one argument")
        return names

    def _stmt(self, node: ast.stmt, indent: int):
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                raise ValueError("Only simple assignments are supported in DSL kernels")
            target = node.targets[0].id
            value = self._expr(node.value)
            self._emit(f"{target} = {value}", indent)
            return
        if isinstance(node, ast.AugAssign):
            if not isinstance(node.target, ast.Name):
                raise ValueError("Only simple augmented assignments are supported")
            target = node.target.id
            op = self._binop(node.op)
            value = self._expr(node.value)
            self._emit(f"{target} = {target} {op} {value}", indent)
            return
        if isinstance(node, ast.Return):
            if node.value is None:
                raise ValueError("DSL kernel return must have a value")
            value = self._expr(node.value)
            self._emit(f"return {value}", indent)
            return
        if isinstance(node, ast.Expr):
            value = self._expr(node.value)
            self._emit(value, indent)
            return
        if isinstance(node, ast.If):
            self._if_stmt(node, indent)
            return
        if isinstance(node, ast.For):
            self._for_stmt(node, indent)
            return
        if isinstance(node, ast.While):
            self._while_stmt(node, indent)
            return
        if isinstance(node, ast.Break):
            self._emit("break", indent)
            return
        if isinstance(node, ast.Continue):
            self._emit("continue", indent)
            return
        raise ValueError(f"Unsupported DSL statement: {type(node).__name__}")

    def _stmt_block(self, body, indent: int):
        if not body:
            raise ValueError("Empty blocks are not supported in DSL kernels")
        i = 0
        while i < len(body):
            stmt = body[i]
            if (
                isinstance(stmt, ast.If)
                and not stmt.orelse
                and self._block_terminates(stmt.body)
                and i + 1 < len(body)
                and isinstance(body[i + 1], ast.If)
            ):
                merged = ast.If(test=stmt.test, body=stmt.body, orelse=[body[i + 1]])
                self._if_stmt(merged, indent)
                i += 2
                continue
            self._stmt(stmt, indent)
            i += 1

    def _block_terminates(self, body) -> bool:
        if not body:
            return False
        return self._stmt_terminates(body[-1])

    def _stmt_terminates(self, node: ast.stmt) -> bool:
        if isinstance(node, ast.Return | ast.Break | ast.Continue):
            return True
        if isinstance(node, ast.If) and node.orelse:
            return self._block_terminates(node.body) and self._block_terminates(node.orelse)
        return False

    def _if_stmt(self, node: ast.If, indent: int):
        current = node
        first = True
        while True:
            prefix = "if" if first else "elif"
            cond = self._expr(current.test)
            self._emit(f"{prefix} {cond}:", indent)
            self._stmt_block(current.body, indent + 4)
            first = False
            if current.orelse and len(current.orelse) == 1 and isinstance(current.orelse[0], ast.If):
                current = current.orelse[0]
                continue
            break
        if current.orelse:
            self._emit("else:", indent)
            self._stmt_block(current.orelse, indent + 4)

    def _for_stmt(self, node: ast.For, indent: int):
        if node.orelse:
            raise ValueError("for/else is not supported in DSL kernels")
        if not isinstance(node.target, ast.Name):
            raise ValueError("DSL for-loop target must be a simple name")
        if not isinstance(node.iter, ast.Call):
            raise ValueError("DSL for-loop must iterate over range()")
        func_name = self._call_name(node.iter.func)
        if func_name != "range":
            raise ValueError("DSL for-loop must iterate over range()")
        if node.iter.keywords or len(node.iter.args) != 1:
            raise ValueError("DSL range() must take a single argument")
        limit = self._expr(node.iter.args[0])
        self._emit(f"for {node.target.id} in range({limit}):", indent)
        self._stmt_block(node.body, indent + 4)

    def _while_stmt(self, node: ast.While, indent: int):
        if node.orelse:
            raise ValueError("while/else is not supported in DSL kernels")
        cond = self._expr(node.test)
        self._emit(f"while {cond}:", indent)
        self._stmt_block(node.body, indent + 4)

    def _expr(self, node: ast.AST) -> str:  # noqa: C901
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Constant):
            val = node.value
            if isinstance(val, bool):
                return "1" if val else "0"
            if isinstance(val, int | float):
                return repr(val)
            raise ValueError("Unsupported constant in DSL expression")
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.UAdd):
                return f"+{self._expr(node.operand)}"
            if isinstance(node.op, ast.USub):
                return f"-{self._expr(node.operand)}"
            if isinstance(node.op, ast.Not):
                return f"!{self._expr(node.operand)}"
            raise ValueError("Unsupported unary operator in DSL expression")
        if isinstance(node, ast.BinOp):
            left = self._expr(node.left)
            right = self._expr(node.right)
            op = self._binop(node.op)
            return f"({left} {op} {right})"
        if isinstance(node, ast.BoolOp):
            op = "&" if isinstance(node.op, ast.And) else "|"
            values = [self._expr(v) for v in node.values]
            expr = values[0]
            for val in values[1:]:
                expr = f"({expr} {op} {val})"
            return expr
        if isinstance(node, ast.Compare):
            if len(node.ops) != 1 or len(node.comparators) != 1:
                raise ValueError("Chained comparisons are not supported in DSL")
            left = self._expr(node.left)
            right = self._expr(node.comparators[0])
            op = self._cmpop(node.ops[0])
            return f"({left} {op} {right})"
        if isinstance(node, ast.Call):
            func_name = self._call_name(node.func)
            if node.keywords:
                raise ValueError("Keyword arguments are not supported in DSL calls")
            args = ", ".join(self._expr(a) for a in node.args)
            return f"{func_name}({args})"
        if isinstance(node, ast.IfExp):
            cond = self._expr(node.test)
            body = self._expr(node.body)
            orelse = self._expr(node.orelse)
            return f"where({cond}, {body}, {orelse})"
        raise ValueError(f"Unsupported DSL expression: {type(node).__name__}")

    def _call_name(self, node: ast.AST) -> str:
        if isinstance(node, ast.Name):
            return node.id
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in {"np", "numpy", "math"}
        ):
            return node.attr
        raise ValueError("Unsupported call target in DSL")

    def _binop(self, op: ast.operator) -> str:
        for k, v in self._binop_map.items():
            if isinstance(op, k):
                return v
        raise ValueError("Unsupported binary operator in DSL")

    def _cmpop(self, op: ast.cmpop) -> str:
        for k, v in self._cmp_map.items():
            if isinstance(op, k):
                return v
        raise ValueError("Unsupported comparison in DSL")


class DSLReducer:
    _binop_map: ClassVar[dict[type[ast.operator], str]] = DSLBuilder._binop_map
    _cmp_map: ClassVar[dict[type[ast.cmpop], str]] = DSLBuilder._cmp_map

    def __init__(self, max_unroll: int = 64):
        self._env: dict[str, str] = {}
        self._const_env: dict[str, object] = {}
        self._return_expr: str | None = None
        self._max_unroll = max_unroll

    def reduce(self, func_node: ast.FunctionDef):
        input_names = self._args(func_node.args)
        if not func_node.body:
            return None
        for stmt in func_node.body:
            if not self._stmt(stmt):
                return None
            if self._return_expr is not None:
                break
        if self._return_expr is None:
            return None
        return self._return_expr, input_names

    def _args(self, args: ast.arguments):
        if args.vararg or args.kwarg or args.kwonlyargs:
            raise ValueError("DSL kernel does not support *args/**kwargs/kwonly args")
        if args.defaults or args.kw_defaults:
            raise ValueError("DSL kernel does not support default arguments")
        names = [a.arg for a in (args.posonlyargs + args.args)]
        if not names:
            raise ValueError("DSL kernel must accept at least one argument")
        return names

    def _stmt(self, node: ast.stmt) -> bool:  # noqa: C901
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                return False
            target = node.targets[0].id
            value = self._expr(node.value)
            self._env[target] = value
            const_val = self._const_eval(node.value)
            if const_val is None:
                self._const_env.pop(target, None)
            else:
                self._const_env[target] = const_val
            return True
        if isinstance(node, ast.AugAssign):
            if not isinstance(node.target, ast.Name):
                return False
            target = node.target.id
            op = self._binop(node.op)
            value = self._expr(node.value)
            left = self._env.get(target, target)
            left_const = self._const_env.get(target)
            right_const = self._const_eval(node.value)
            simplified = self._simplify_binop_expr(op, left, value, left_const, right_const)
            self._env[target] = simplified
            if left_const is None or right_const is None:
                self._const_env.pop(target, None)
            else:
                self._const_env[target] = self._apply_binop(left_const, right_const, node.op)
            return True
        if isinstance(node, ast.Return):
            if node.value is None:
                return False
            self._return_expr = self._expr(node.value)
            return True
        if isinstance(node, ast.If):
            test_val = self._const_eval(node.test)
            if test_val is None:
                return False
            branch = node.body if bool(test_val) else node.orelse
            if not branch:
                return True
            for stmt in branch:
                if not self._stmt(stmt):
                    return False
                if self._return_expr is not None:
                    return True
            return True
        if isinstance(node, ast.For):
            if node.orelse:
                return False
            if not isinstance(node.target, ast.Name):
                return False
            if not isinstance(node.iter, ast.Call):
                return False
            func_name = self._call_name(node.iter.func)
            if func_name != "range":
                return False
            if node.iter.keywords or len(node.iter.args) != 1:
                return False
            limit_val = self._const_eval(node.iter.args[0])
            if limit_val is None or not isinstance(limit_val, int):
                return False
            if limit_val < 0 or limit_val > self._max_unroll:
                return False
            loop_var = node.target.id
            old_env = self._env.get(loop_var)
            old_const = self._const_env.get(loop_var)
            for i in range(limit_val):
                self._env[loop_var] = str(i)
                self._const_env[loop_var] = i
                for stmt in node.body:
                    if not self._stmt(stmt):
                        if old_env is None:
                            self._env.pop(loop_var, None)
                        else:
                            self._env[loop_var] = old_env
                        if old_const is None:
                            self._const_env.pop(loop_var, None)
                        else:
                            self._const_env[loop_var] = old_const
                        return False
                    if self._return_expr is not None:
                        break
                if self._return_expr is not None:
                    break
            if old_env is None:
                self._env.pop(loop_var, None)
            else:
                self._env[loop_var] = old_env
            if old_const is None:
                self._const_env.pop(loop_var, None)
            else:
                self._const_env[loop_var] = old_const
            return True
        return False

    def _expr(self, node: ast.AST) -> str:  # noqa: C901
        const_val = self._const_eval(node)
        if const_val is not None:
            if isinstance(const_val, bool):
                return "1" if const_val else "0"
            return repr(const_val)
        if isinstance(node, ast.Name):
            if node.id in self._env:
                val = self._env[node.id]
                # Avoid double-wrapping if already parenthesized or is a function call
                if (val.startswith("(") and val.endswith(")")) or "(" in val:
                    return val
                return f"({val})"
            return node.id
        if isinstance(node, ast.Constant):
            val = node.value
            if isinstance(val, bool):
                return "1" if val else "0"
            if isinstance(val, int | float):
                return repr(val)
            raise ValueError("Unsupported constant in DSL expression")
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.UAdd):
                return f"+{self._expr(node.operand)}"
            if isinstance(node.op, ast.USub):
                return f"-{self._expr(node.operand)}"
            if isinstance(node.op, ast.Not):
                return f"!{self._expr(node.operand)}"
            raise ValueError("Unsupported unary operator in DSL expression")
        if isinstance(node, ast.BinOp):
            left = self._expr(node.left)
            right = self._expr(node.right)
            op = self._binop(node.op)
            left_const = self._const_eval(node.left)
            right_const = self._const_eval(node.right)
            return self._simplify_binop_expr(op, left, right, left_const, right_const)
        if isinstance(node, ast.BoolOp):
            op = "&" if isinstance(node.op, ast.And) else "|"
            values = [self._expr(v) for v in node.values]
            expr = values[0]
            for val in values[1:]:
                expr = f"({expr} {op} {val})"
            return expr
        if isinstance(node, ast.Compare):
            if len(node.ops) != 1 or len(node.comparators) != 1:
                raise ValueError("Chained comparisons are not supported in DSL")
            left = self._expr(node.left)
            right = self._expr(node.comparators[0])
            op = self._cmpop(node.ops[0])
            return f"({left} {op} {right})"
        if isinstance(node, ast.Call):
            func_name = self._call_name(node.func)
            if node.keywords:
                raise ValueError("Keyword arguments are not supported in DSL calls")
            args = ", ".join(self._expr(a) for a in node.args)
            return f"{func_name}({args})"
        if isinstance(node, ast.IfExp):
            cond = self._expr(node.test)
            body = self._expr(node.body)
            orelse = self._expr(node.orelse)
            return f"where({cond}, {body}, {orelse})"
        raise ValueError(f"Unsupported DSL expression: {type(node).__name__}")

    def _call_name(self, node: ast.AST) -> str:
        if isinstance(node, ast.Name):
            return node.id
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in {"np", "numpy", "math"}
        ):
            return node.attr
        raise ValueError("Unsupported call target in DSL")

    def _binop(self, op: ast.operator) -> str:
        for k, v in self._binop_map.items():
            if isinstance(op, k):
                return v
        raise ValueError("Unsupported binary operator in DSL")

    def _cmpop(self, op: ast.cmpop) -> str:
        for k, v in self._cmp_map.items():
            if isinstance(op, k):
                return v
        raise ValueError("Unsupported comparison in DSL")

    def _const_eval(self, node: ast.AST):  # noqa: C901
        if isinstance(node, ast.Constant):
            if isinstance(node.value, int | float | bool):
                return node.value
            return None
        if isinstance(node, ast.Name):
            return self._const_env.get(node.id)
        if isinstance(node, ast.UnaryOp):
            val = self._const_eval(node.operand)
            if val is None:
                return None
            if isinstance(node.op, ast.UAdd):
                return +val
            if isinstance(node.op, ast.USub):
                return -val
            if isinstance(node.op, ast.Not):
                return not val
            return None
        if isinstance(node, ast.BinOp):
            left = self._const_eval(node.left)
            right = self._const_eval(node.right)
            if left is None or right is None:
                return None
            return self._apply_binop(left, right, node.op)
        if isinstance(node, ast.BoolOp):
            vals = [self._const_eval(v) for v in node.values]
            if any(v is None for v in vals):
                return None
            if isinstance(node.op, ast.And):
                return all(vals)
            if isinstance(node.op, ast.Or):
                return any(vals)
            return None
        if isinstance(node, ast.Compare):
            if len(node.ops) != 1 or len(node.comparators) != 1:
                return None
            left = self._const_eval(node.left)
            right = self._const_eval(node.comparators[0])
            if left is None or right is None:
                return None
            return self._apply_cmp(left, right, node.ops[0])
        return None

    def _apply_binop(self, left, right, op):
        if isinstance(op, ast.Add):
            return left + right
        if isinstance(op, ast.Sub):
            return left - right
        if isinstance(op, ast.Mult):
            return left * right
        if isinstance(op, ast.Div):
            return left / right
        if isinstance(op, ast.FloorDiv):
            return left // right
        if isinstance(op, ast.Mod):
            return left % right
        if isinstance(op, ast.Pow):
            return left**right
        if isinstance(op, ast.BitAnd):
            return left & right
        if isinstance(op, ast.BitOr):
            return left | right
        if isinstance(op, ast.BitXor):
            return left ^ right
        if isinstance(op, ast.LShift):
            return left << right
        if isinstance(op, ast.RShift):
            return left >> right
        return None

    def _apply_cmp(self, left, right, op):
        if isinstance(op, ast.Eq):
            return left == right
        if isinstance(op, ast.NotEq):
            return left != right
        if isinstance(op, ast.Lt):
            return left < right
        if isinstance(op, ast.LtE):
            return left <= right
        if isinstance(op, ast.Gt):
            return left > right
        if isinstance(op, ast.GtE):
            return left >= right
        return None

    def _simplify_binop_expr(self, op, left_expr, right_expr, left_const, right_const):
        if op == "+":
            if self._is_zero(left_const):
                return right_expr
            if self._is_zero(right_const):
                return left_expr
        if op == "-" and self._is_zero(right_const):
            return left_expr
        if op == "*":
            if self._is_one(left_const):
                return right_expr
            if self._is_one(right_const):
                return left_expr
        return f"({left_expr} {op} {right_expr})"

    def _is_zero(self, value):
        return isinstance(value, int | float | bool) and value == 0

    def _is_one(self, value):
        return isinstance(value, int | float | bool) and value == 1
