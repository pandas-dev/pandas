"""Transpile a blosc2 DSL kernel to JavaScript, and run it from a lazyudf callable.

Browser/Pyodide-only payoff: V8 JIT-compiles the emitted scalar loop to optimized native
code, which in the Newton-fractal demo beats blosc2's WASM JIT ~3.3x and the no-JIT
interpreter ~11x. See plans/dsl-js.md.

Public API:
    dsl_to_js(kernel)         -> (js_source, param_names)  # pure stdlib, runs anywhere
    build_js_module(k, ndim)  -> js_source string          # ndim needed for index symbols
    js_kernel(kernel, shape)  -> callable for lazyudf(...)  # needs Pyodide `js` to *run*

`kernel` may be a blosc2 DSLKernel (has .dsl_source), a plain function, or a source string.

Kernels may use index/shape symbols (`_i0`/`_i1`/.., `_n0`/.., `_flat_idx`); they become
trailing kernel parameters and the runtime driver reconstructs each element's global
coordinate per block. That requires the output rank, so `build_js_module`/`js_kernel` need
`ndim`/`shape`; without them, index-symbol kernels raise (and the caller falls back).
"""

from __future__ import annotations

import ast
import inspect
import json
import textwrap

# Wired into lazyexpr via jit_backend="js": a DSL kernel is transpiled here and run as a
# plain per-block callable. Browser/Pyodide only (js_kernel imports `js` at call time).

# Canonical signature order for index/shape symbols passed to the transpiled kernel.
_INDEX_SYMBOL_ORDER = ("_i0", "_i1", "_i2", "_n0", "_n1", "_n2", "_ndim", "_flat_idx")
_INDEX_SYMBOLS = set(_INDEX_SYMBOL_ORDER)

# numpy/math function name -> JS Math.* name (numpy aliases included).
_MATH = {
    "sin": "sin",
    "cos": "cos",
    "tan": "tan",
    "asin": "asin",
    "acos": "acos",
    "atan": "atan",
    "atan2": "atan2",
    "arcsin": "asin",
    "arccos": "acos",
    "arctan": "atan",
    "arctan2": "atan2",
    "sinh": "sinh",
    "cosh": "cosh",
    "tanh": "tanh",
    "exp": "exp",
    "log": "log",
    "log2": "log2",
    "log10": "log10",
    "sqrt": "sqrt",
    "cbrt": "cbrt",
    "pow": "pow",
    "power": "pow",
    "hypot": "hypot",
    "floor": "floor",
    "ceil": "ceil",
    "trunc": "trunc",
    "round": "round",
    "abs": "abs",
    "absolute": "abs",
    "fabs": "abs",
    "sign": "sign",
    "min": "min",
    "max": "max",
    "minimum": "min",
    "maximum": "max",
}

_BIN = {
    ast.Add: "+",
    ast.Sub: "-",
    ast.Mult: "*",
    ast.Div: "/",
    ast.BitAnd: "&",
    ast.BitOr: "|",
    ast.BitXor: "^",
    ast.LShift: "<<",
    ast.RShift: ">>",
}
_AUG = {
    ast.Add: "+=",
    ast.Sub: "-=",
    ast.Mult: "*=",
    ast.Div: "/=",
    ast.BitAnd: "&=",
    ast.BitOr: "|=",
    ast.BitXor: "^=",
    ast.LShift: "<<=",
    ast.RShift: ">>=",
}
_CMP = {
    ast.Eq: "===",
    ast.NotEq: "!==",
    ast.Lt: "<",
    ast.LtE: "<=",
    ast.Gt: ">",
    ast.GtE: ">=",
}

JS_PRELUDE = "const pymod = (a, b) => (((a % b) + b) % b);"


class _DSLToJSError(Exception):
    pass


def _get_source(obj) -> str:
    if hasattr(obj, "dsl_source"):
        src = obj.dsl_source
    elif isinstance(obj, str):
        src = obj
    elif callable(obj):
        src = inspect.getsource(obj)
    else:
        raise _DSLToJSError(f"cannot get DSL source from {obj!r}")
    return textwrap.dedent(src)


class _Transpiler:
    def transpile(self, func: ast.FunctionDef):
        self.params = [a.arg for a in func.args.args]
        used_index = self._collect_index_symbols(func)
        hoist = self._hoist_names(func)
        body = self._block(func.body, 1)
        # Index/shape symbols (`_i0`, `_n1`, `_flat_idx`, ...) become extra trailing
        # parameters; the runtime driver computes them per element (see _module_with_index).
        sig = self.params + used_index
        head = f"function {func.name}({', '.join(sig)}) {{\n"
        decl = f"  let {', '.join(sorted(hoist))};\n" if hoist else ""
        return head + decl + body + "}", list(self.params), used_index

    # -- scope analysis -------------------------------------------------
    def _collect_index_symbols(self, func):
        """Return the index/shape symbols the kernel references, in canonical order."""
        used = {
            node.id for node in ast.walk(func) if isinstance(node, ast.Name) and node.id in _INDEX_SYMBOLS
        }
        return [s for s in _INDEX_SYMBOL_ORDER if s in used]

    def _hoist_names(self, func):
        assigned, fortargets = set(), set()
        for node in ast.walk(func):
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        assigned.add(t.id)
            elif isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name):
                assigned.add(node.target.id)
            elif isinstance(node, ast.For) and isinstance(node.target, ast.Name):
                fortargets.add(node.target.id)
        return assigned - set(self.params) - fortargets

    # -- statements -----------------------------------------------------
    def _block(self, stmts, ind):
        return "".join(self._stmt(s, ind) for s in stmts)

    def _stmt(self, node, ind):
        pad = "  " * ind
        if isinstance(node, ast.Assign):
            return f"{pad}{node.targets[0].id} = {self._expr(node.value)};\n"
        if isinstance(node, ast.AugAssign):
            return pad + self._augassign(node) + "\n"
        if isinstance(node, ast.Return):
            return f"{pad}return {self._expr(node.value)};\n"
        if isinstance(node, ast.Expr):
            return f"{pad}{self._expr(node.value)};\n"
        if isinstance(node, ast.If):
            return self._if(node, ind)
        if isinstance(node, ast.For):
            return self._for(node, ind)
        if isinstance(node, ast.While):
            return f"{pad}while ({self._expr(node.test)}) {{\n{self._block(node.body, ind + 1)}{pad}}}\n"
        if isinstance(node, ast.Break):
            return f"{pad}break;\n"
        if isinstance(node, ast.Continue):
            return f"{pad}continue;\n"
        raise _DSLToJSError(f"unsupported statement: {type(node).__name__}")

    def _augassign(self, node):
        t, val, op = node.target.id, self._expr(node.value), type(node.op)
        if op in _AUG:
            return f"{t} {_AUG[op]} {val};"
        if op is ast.Pow:
            return f"{t} = Math.pow({t}, {val});"
        if op is ast.FloorDiv:
            return f"{t} = Math.floor({t} / {val});"
        if op is ast.Mod:
            return f"{t} = pymod({t}, {val});"
        raise _DSLToJSError(f"unsupported augmented op: {op.__name__}")

    def _if(self, node, ind):
        pad = "  " * ind
        s = f"{pad}if ({self._expr(node.test)}) {{\n{self._block(node.body, ind + 1)}{pad}}}"
        if node.orelse:
            if len(node.orelse) == 1 and isinstance(node.orelse[0], ast.If):
                s += " else " + self._if(node.orelse[0], ind).lstrip()
            else:
                s += f" else {{\n{self._block(node.orelse, ind + 1)}{pad}}}\n"
                return s
        return s + "\n"

    def _for(self, node, ind):
        pad = "  " * ind
        var = node.target.id
        args = node.iter.args
        if len(args) == 1:
            start, stop, step, stepnode = "0", self._expr(args[0]), "1", None
        elif len(args) == 2:
            start, stop, step, stepnode = self._expr(args[0]), self._expr(args[1]), "1", None
        else:
            start, stop, step, stepnode = (
                self._expr(args[0]),
                self._expr(args[1]),
                self._expr(args[2]),
                args[2],
            )
        cond = f"{var} > {stop}" if _neg_literal(stepnode) else f"{var} < {stop}"
        return (
            f"{pad}for (let {var} = {start}; {cond}; {var} += {step}) {{\n"
            f"{self._block(node.body, ind + 1)}{pad}}}\n"
        )

    # -- expressions ----------------------------------------------------
    def _expr(self, node):
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Constant):
            return _const(node.value)
        if isinstance(node, ast.UnaryOp):
            sym = {ast.Not: "!", ast.USub: "-", ast.UAdd: "+"}[type(node.op)]
            return f"({sym}{self._expr(node.operand)})"
        if isinstance(node, ast.BinOp):
            return self._binop(node)
        if isinstance(node, ast.BoolOp):
            sym = "&&" if isinstance(node.op, ast.And) else "||"
            return "(" + f" {sym} ".join(self._expr(v) for v in node.values) + ")"
        if isinstance(node, ast.Compare):
            op = _CMP[type(node.ops[0])]
            return f"({self._expr(node.left)} {op} {self._expr(node.comparators[0])})"
        if isinstance(node, ast.Call):
            return self._call(node)
        raise _DSLToJSError(f"unsupported expression: {type(node).__name__}")

    def _binop(self, node):
        left, right, op = self._expr(node.left), self._expr(node.right), type(node.op)
        if op is ast.Pow:
            return f"Math.pow({left}, {right})"
        if op is ast.FloorDiv:
            return f"Math.floor({left} / {right})"
        if op is ast.Mod:
            return f"pymod({left}, {right})"
        if op in _BIN:
            return f"({left} {_BIN[op]} {right})"
        raise _DSLToJSError(f"unsupported binary op: {op.__name__}")

    def _call(self, node):
        name = self._call_name(node.func)
        args = [self._expr(a) for a in node.args]
        if name == "where":
            if len(args) != 3:
                raise _DSLToJSError("where() needs 3 args: where(cond, a, b)")
            return f"({args[0]} ? {args[1]} : {args[2]})"
        if name == "int":
            return f"Math.trunc({args[0]})"
        if name == "float":
            return f"({args[0]})"
        if name == "bool":
            return f"(({args[0]}) != 0)"
        if name == "range":
            raise _DSLToJSError("range() is only valid as a for-loop iterator")
        if name in _MATH:
            return f"Math.{_MATH[name]}({', '.join(args)})"
        raise _DSLToJSError(f"unsupported call: {name}()")

    def _call_name(self, node):
        if isinstance(node, ast.Name):
            return node.id
        if (
            isinstance(node, ast.Attribute)
            and isinstance(node.value, ast.Name)
            and node.value.id in {"np", "numpy", "math"}
        ):
            return node.attr
        raise _DSLToJSError("unsupported call target")


def _neg_literal(node) -> bool:
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return node.value < 0
    return isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub)


def _const(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v != v:
            return "NaN"
        if v == float("inf"):
            return "Infinity"
        if v == float("-inf"):
            return "-Infinity"
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v)
    raise _DSLToJSError(f"unsupported constant: {v!r}")


# Transpiling is a pure function of the kernel source, so memoize it: the same kernel is
# typically re-transpiled on every lazyudf evaluation (e.g. an animation loop rebuilds the
# expression each frame), and ast.parse + codegen is non-trivial under WASM.
_TRANSPILE_CACHE: dict[str, tuple] = {}

# module string -> the V8-compiled `__run` function (a Pyodide JsProxy). Populated only when
# a kernel actually runs in-browser (see js_kernel's bridge); empty off-WASM.
_RUN_CACHE: dict = {}


def _transpile(kernel):
    """Transpile *kernel* to JS. Returns (js_source, params, index_symbols, func_name)."""
    src = _get_source(kernel)
    hit = _TRANSPILE_CACHE.get(src)
    if hit is not None:
        return hit
    tree = ast.parse(src)
    func = next((n for n in tree.body if isinstance(n, ast.FunctionDef)), None)
    if func is None:
        raise _DSLToJSError("no function definition found in DSL source")
    js_src, params, used_index = _Transpiler().transpile(func)
    result = (js_src, params, used_index, func.name)
    _TRANSPILE_CACHE[src] = result
    return result


def dsl_to_js(kernel):
    """Transpile a DSL kernel to a JS function string. Returns (js_source, param_names)."""
    js_src, params, _used, _name = _transpile(kernel)
    return js_src, params


def _max_index_dim(used_index) -> int:
    """Highest axis referenced by `_iK`/`_nK` symbols (-1 if only `_ndim`/`_flat_idx`/none)."""
    dims = [int(s[2:]) for s in used_index if s[:2] in ("_i", "_n") and s[2:].isdigit()]
    return max(dims) if dims else -1


def _op_call_args(nparams):
    return [f"(isarr[{k}] ? ops[{k}][i] : ops[{k}])" for k in range(nparams)]


def _module_with_index(kernel_js, fname, params, used_index) -> str:
    """Driver for kernels that use index/shape symbols.

    Signature ``__run(ops, isarr, out, n, off, gshape, cshape)``: `off` is the block's
    global start coord, `gshape` the whole-array shape, `cshape` the block shape (all JS
    arrays).  Each element's local coord is unravelled from its flat position (C-order),
    then the referenced symbols are derived and passed as trailing kernel args."""
    decls = []
    for s in used_index:
        if s in ("_flat_idx", "_ndim"):
            continue
        k = int(s[2:])
        rhs = f"off[{k}] + loc[{k}]" if s.startswith("_i") else f"gshape[{k}]"
        decls.append(f"    const {s} = {rhs};")
    if "_ndim" in used_index:
        decls.append("    const _ndim = d;")
    if "_flat_idx" in used_index:
        decls.append("    let _flat_idx = 0;")
        decls.append(
            "    for (let k = 0; k < d; k++) _flat_idx = _flat_idx * gshape[k] + (off[k] + loc[k]);"
        )
    call = f"{fname}({', '.join(_op_call_args(len(params)) + used_index)})"
    driver = "\n".join(
        [
            "function __run(ops, isarr, out, n, off, gshape, cshape) {",
            "  const d = cshape.length;",
            "  const loc = new Array(d);",
            "  for (let i = 0; i < n; i++) {",
            "    let rem = i;",
            "    for (let k = d - 1; k >= 0; k--) { loc[k] = rem % cshape[k]; rem = (rem - loc[k]) / cshape[k]; }",
            *decls,
            f"    out[i] = {call};",
            "  }",
            "}",
        ]
    )
    return f"{JS_PRELUDE}\n{kernel_js}\n{driver}\nreturn __run;"


def build_js_module(kernel, ndim: int | None = None) -> str:
    """Self-contained JS: prelude + kernel + a runtime element driver returning ``__run``.

    Kernels without index/shape symbols get a flat ``__run(ops, isarr, out, n)`` driver.
    Kernels that use `_i0`/`_n0`/`_flat_idx` get the index-aware driver (see
    :func:`_module_with_index`) and require *ndim* (the output rank) so the referenced
    axes can be validated; ``ndim=None`` raises for such kernels."""
    kernel_js, params, used_index, fname = _transpile(kernel)
    if used_index:
        if ndim is None:
            raise _DSLToJSError("kernel uses index/shape symbols; the output ndim must be supplied")
        max_dim = _max_index_dim(used_index)
        if max_dim >= ndim:
            raise _DSLToJSError(f"kernel references axis {max_dim} but the output is {ndim}-D")
        return _module_with_index(kernel_js, fname, params, used_index)
    call_args = ", ".join(_op_call_args(len(params)))
    driver = (
        f"function __run(ops, isarr, out, n) {{ "
        f"for (let i = 0; i < n; i++) out[i] = {fname}({call_args}); }}"
    )
    return f"{JS_PRELUDE}\n{kernel_js}\n{driver}\nreturn __run;"


def js_kernel(kernel, shape=None):
    """Return a lazyudf-compatible callable that runs the transpiled JS (Pyodide only).

    *shape* is the whole-array output shape; it is required for kernels that use index/shape
    symbols (so the driver knows the rank and global geometry) and ignored otherwise."""
    ndim = len(shape) if shape is not None else None
    _, _, used_index, _ = _transpile(kernel)  # cached
    module = build_js_module(kernel, ndim=ndim)
    uses_index = bool(used_index)
    gshape = tuple(int(s) for s in shape) if shape is not None else None
    run = None  # lazily created in-browser

    def bridge(inputs, output, offset=None):
        nonlocal run
        import numpy as np
        from js import Array, Float64Array, Uint8Array  # Pyodide

        if run is None:
            # Reuse the V8-compiled function across lazyudf evaluations of the same kernel:
            # the module is a pure function of (source, ndim), so js.eval need run only once
            # per distinct module instead of once per frame.
            run = _RUN_CACHE.get(module)
            if run is None:
                import js

                run = js.eval(f"(function() {{ {module} }})()")
                _RUN_CACHE[module] = run

        n = int(output.size)
        # Pass real JS Arrays, not Python lists: a Python list arrives in JS as a PyProxy,
        # so each ops[k][i] in the hot loop would cross the Python<->JS boundary (~10x slower).
        ops = Array.new()
        isarr = Array.new()
        for x in inputs:
            if isinstance(x, np.ndarray) and x.ndim > 0:
                ops.push(
                    _to_jsf64(
                        np.ascontiguousarray(x, dtype=np.float64).reshape(-1), Float64Array, Uint8Array
                    )
                )
                isarr.push(True)
            else:
                ops.push(float(x))
                isarr.push(False)
        out_js = Float64Array.new(n)
        if uses_index:
            off = offset if offset is not None else (0,) * output.ndim
            run(
                ops,
                isarr,
                out_js,
                n,
                _to_jsint(off, Array),
                _to_jsint(gshape, Array),
                _to_jsint(output.shape, Array),
            )
        else:
            run(ops, isarr, out_js, n)
        # ponytail: per-block to_js()/to_bytes() copies; swap to a zero-copy HEAPF64 view
        # onto WASM linear memory only if marshaling shows up as the bottleneck.
        res = np.frombuffer(bytes(out_js.to_bytes()), dtype=np.float64)
        output.reshape(-1)[:] = res
        return output

    bridge.js_source = module
    return bridge


def _to_jsf64(xf, Float64Array, Uint8Array):
    u8 = Uint8Array.new(xf.nbytes)
    u8.assign(xf.tobytes())  # Pyodide TypedArray.assign(buffer) copies bytes in
    return Float64Array.new(u8.buffer)


def _to_jsint(seq, Array):
    """Small geometry vector (offset/shape) -> a real JS Array of ints (avoids PyProxy)."""
    arr = Array.new()
    for v in seq:
        arr.push(int(v))
    return arr
