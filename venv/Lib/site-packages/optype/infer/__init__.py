"""Structurally infer the ``optype`` protocols required by a function."""

import re
from collections import defaultdict
from collections.abc import Callable, Collection, Iterator, Mapping
from inspect import Parameter, signature
from typing import Any, NamedTuple, cast, final

from ._spy import _Fork, _fork, _Spy, _SpyBytes, _SpyObject, _SpyStr, _TraceItem
from optype._core import _can, _has
from optype.inspect import get_protocol_members

__all__ = ("infer",)

type _AnyFunc = Callable[..., Any]

_VARIADIC = frozenset({Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD})

_TYPEVARS = "TUVWXYZ"

_ATTRIBUTE_DUNDERS = frozenset({
    "__delattr__",
    "__getattr__",
    "__getattribute__",
    "__setattr__",
})


def _dunder_protocols() -> dict[str, str]:
    return {
        dunder: name
        for name in _can.__all__
        if not name.endswith(("Self", "Same"))
        if (dunder := "__" + name.removeprefix("Can").lower() + "__")
        not in _ATTRIBUTE_DUNDERS
    }


def _attr_protocols() -> dict[str, str]:
    return {
        next(iter(members)): name
        for name in _has.__all__
        if len(members := get_protocol_members(getattr(_has, name))) == 1
    }


_DUNDER_PROTOCOL_MAP = _dunder_protocols()
_ATTR_PROTOCOL_MAP = _attr_protocols()

_COERCION_FALLBACK = {
    "__float__": ("__index__",),
    "__int__": ("__index__",),
    "__complex__": ("__float__", "__index__"),
}

_COERCION_UNION = {
    dunder: " | ".join(map(_DUNDER_PROTOCOL_MAP.__getitem__, (dunder, *fallback)))
    for dunder, fallback in _COERCION_FALLBACK.items()
}

_FORWARD_ARITH = frozenset(
    dunder
    for dunder, proto in _DUNDER_PROTOCOL_MAP.items()
    if "CanR" + proto.removeprefix("Can") in _DUNDER_PROTOCOL_MAP.values()
)


_DOC_SIGNATURE = re.compile(r"\b(\w+)\(([^)]*)\)")
_DOC_PARAM = re.compile(r"(?:^|,)\s*\**([a-zA-Z_]\w*)")


class _Op(NamedTuple):
    proto: str
    args: tuple[Any, ...]
    kwargs: dict[str, Any]
    ret: Any


type _Vars = dict[int, str]


def _resolve(trace: _TraceItem) -> _Op:
    if trace.attr in _ATTRIBUTE_DUNDERS:
        name = trace.args[0]
        if name not in _ATTR_PROTOCOL_MAP:
            raise NotImplementedError(name)
        return _Op(_ATTR_PROTOCOL_MAP[name], (), {}, trace.return_)
    if trace.attr in _COERCION_UNION:
        return _Op(_COERCION_UNION[trace.attr], (), {}, trace.return_)
    if trace.attr in _DUNDER_PROTOCOL_MAP:
        proto = _DUNDER_PROTOCOL_MAP[trace.attr]
        return _Op(proto, trace.args, trace.kwargs, trace.return_)
    raise NotImplementedError(trace.attr)


def _analyze(
    params: list[_SpyObject],
    results: list[object],
) -> tuple[list[_SpyObject], dict[int, int]]:
    appear: defaultdict[int, int] = defaultdict(int)
    for spy in params:
        appear[id(spy)] += 1
    for result in results:
        if isinstance(result, _SpyObject):
            appear[id(result)] += 1

    order: list[_SpyObject] = []
    seen: set[int] = set()
    stack = list(reversed(params))
    while stack:
        spy = stack.pop()
        if id(spy) in seen:
            continue
        seen.add(id(spy))
        order.append(spy)
        for op in spy.__optype_trace__:
            for value in (*op.args, *op.kwargs.values()):
                if isinstance(value, _SpyObject):
                    appear[id(value)] += 1
            if isinstance(op.return_, _SpyObject):
                appear[id(op.return_)] += 1
                stack.append(op.return_)
    return order, appear


def _select(params: tuple[str | int, ...], names: list[str]) -> list[str]:
    selected: list[str] = []
    for p in params:
        if isinstance(p, int):
            if not -len(names) <= p < len(names):
                msg = f"no parameter at position {p}"
                raise ValueError(msg)
            selected.append(names[p])
        elif p in names:
            selected.append(p)
        else:
            msg = f"unknown parameter {p!r}"
            raise ValueError(msg)
    return selected or names


def _return_spies(value: object) -> Iterator[_SpyObject]:
    match value:
        case _SpyObject():
            yield value
        case list() | set() | frozenset():
            for item in cast("Collection[object]", value):
                yield from _return_spies(item)
        case dict():
            mapping = cast("Mapping[object, object]", value)
            for key, val in mapping.items():
                yield from _return_spies(key)
                yield from _return_spies(val)
        case _:
            pass


@final
class _Renderer:
    """Render an inferred ``def`` signature from the recorded spy traces."""

    def __init__(
        self,
        names: list[str],
        selected: list[str],
        spies: dict[str, _SpyObject],
        results: list[object],
        optional: frozenset[str],
    ) -> None:
        self._selected = selected
        self._spies = spies
        self._results = results
        self._optional = optional

        param_spies = [spies[name] for name in names]
        order, appear = _analyze(param_spies, results)
        param_ids = {id(spy) for spy in param_spies}

        self._vars: _Vars = {}
        self._result_spies: list[_SpyObject] = []
        for result in results:
            for spy in _return_spies(result):
                if id(spy) in param_ids or id(spy) in self._vars:
                    continue
                n = len(self._result_spies)
                self._vars[id(spy)] = "R" if not n else f"R{n + 1}"
                self._result_spies.append(spy)
        self._pool = [spies[name] for name in names if appear[id(spies[name])] >= 2]
        self._pool += [
            spy
            for spy in order
            if appear[id(spy)] >= 2
            and id(spy) not in param_ids
            and id(spy) not in self._vars
        ]
        for i, spy in enumerate(self._pool):
            self._vars[id(spy)] = _TYPEVARS[i] if i < len(_TYPEVARS) else f"T{i}"

    def union(self, values: list[Any]) -> str | None:
        literals: list[Any] = []
        names: list[str] = []
        for value in values:
            if isinstance(value, int | str | bytes) and not isinstance(value, _Spy):
                literals.append(value)
            else:
                names.append(self.return_type(value))
        if "float" in names or "complex" in names:
            literals = [value for value in literals if not isinstance(value, int)]
        if "complex" in names:
            names = [name for name in names if name != "float"]

        parts: list[str] = []
        if literals:
            parts.append(f"Literal[{', '.join(map(repr, dict.fromkeys(literals)))}]")
        parts.extend(dict.fromkeys(names))
        return " | ".join(parts) or None

    def returns(self, members: list[_Op]) -> str | None:
        named: list[str] = []
        traces: list[_TraceItem] = []
        for m in members:
            if not isinstance(m.ret, _SpyObject):
                continue
            if (var := self._vars.get(id(m.ret))) is not None:
                named.append(var)
            else:
                traces.extend(m.ret.__optype_trace__)
        parts = list(dict.fromkeys(named))
        if traces and (formatted := self.traces(traces)):
            parts.append(formatted)
        return " & ".join(parts) or None

    def group(self, proto: str, members: list[_Op]) -> str:
        parts = [
            arg
            for i in range(len(members[0].args))
            if (arg := self.union([m.args[i] for m in members])) is not None
        ]
        parts += [
            f"{key}={kw}"
            for key in members[0].kwargs
            if (kw := self.union([m.kwargs[key] for m in members])) is not None
        ]
        if (ret := self.returns(members)) is not None:
            parts.append(ret)

        return f"{proto}[{', '.join(parts)}]" if parts else proto

    def traces(self, traces: list[_TraceItem]) -> str:
        groups: dict[tuple[str, int, tuple[str, ...]], list[_Op]] = {}
        for trace in traces:
            op = _resolve(trace)
            key = op.proto, len(op.args), tuple(sorted(op.kwargs))
            groups.setdefault(key, []).append(op)
        parts = [(key[0], self.group(key[0], group)) for key, group in groups.items()]
        wrap = len(parts) > 1
        return " & ".join(
            f"({s})" if wrap and " | " in proto else s for proto, s in parts
        )

    def spy(self, spy: _Spy) -> str:
        return self.traces(spy.__optype_trace__)

    def slot(self, spy: _SpyObject) -> str:
        return self._vars.get(id(spy)) or self.spy(spy) or "object"

    def typevar(self, spy: _SpyObject) -> str:
        var = self._vars[id(spy)]
        return f"{var}: {bound}" if (bound := self.spy(spy)) else var

    def return_type(self, result: object) -> str:
        match result:
            case _SpyObject():
                return self._vars.get(id(result), "object")
            case _SpyStr():
                return "str"
            case _SpyBytes():
                return "bytes"
            case None:
                return "None"
            case _:
                return self._container(result)

    def _container(self, result: object) -> str:
        name = type(result).__name__
        match result:
            case dict():
                mapping = cast("Mapping[object, object]", result)
                key = self.union([*mapping]) or "object"
                val = self.union([*mapping.values()]) or "object"
                return f"dict[{key}, {val}]" if mapping else name
            case list() | set() | frozenset():
                inner = self.union([*cast("Collection[object]", result)])
                return f"{name}[{inner}]" if inner else name
            case _:
                return name

    def return_types(self) -> str:
        return " | ".join(dict.fromkeys(map(self.return_type, self._results)))

    def render(self) -> str:
        typevars = [self.typevar(spy) for spy in self._pool + self._result_spies]
        generics = f"[{', '.join(typevars)}]" if typevars else ""
        params = ", ".join(
            f"{name}: {slot}"
            for name in self._selected
            if (slot := self.slot(self._spies[name])) != "object"
            or name not in self._optional
        )
        return f"{generics}({params}) -> {self.return_types()}"


def _doc_params(func: _AnyFunc) -> list[str] | None:
    name = getattr(func, "__name__", "")
    if not name:
        return None
    for match in _DOC_SIGNATURE.finditer(func.__doc__ or ""):
        if match[1] == name:
            params = match[2].replace("[", "").replace("]", "")
            return _DOC_PARAM.findall(params) or None
    return None


def _parameters(func: _AnyFunc) -> Mapping[str, Parameter]:
    try:
        return signature(func).parameters
    except ValueError as exc:
        if (names := _doc_params(func)) is None:
            raise NotImplementedError(str(exc)) from exc
        return {n: Parameter(n, Parameter.POSITIONAL_OR_KEYWORD) for n in names}


def _reflect(param_spies: list[_SpyObject], results: list[object]) -> None:
    order, _ = _analyze(param_spies, results)
    added: defaultdict[int, list[_TraceItem]] = defaultdict(list)
    for spy in order:
        keep: list[_TraceItem] = []
        for item in spy.__optype_trace__:
            rhs = item.args[0] if item.args else None
            if item.attr in _FORWARD_ARITH and isinstance(rhs, _SpyObject):
                reflected = _TraceItem("__r" + item.attr[2:], (spy,), {}, item.return_)
                added[id(rhs)].append(reflected)
            else:
                keep.append(item)
        spy.__optype_trace__ = keep

    for spy in order:
        spy.__optype_trace__ += added[id(spy)]


_FORK_LIMIT = 64  # bound fork depth so a forking loop can't diverge


def _explore[T](
    func: Callable[..., T],
    args: list[_SpyObject],
    kwds: dict[str, _SpyObject],
) -> list[T]:
    results: list[T] = []
    stack: list[list[bool]] = [[]]
    while stack:
        plan = stack.pop()
        token = _fork.set(iter(plan))
        try:
            results.append(func(*args, **kwds))
        except _Fork:
            if len(plan) < _FORK_LIMIT:
                stack.extend(([*plan, False], [*plan, True]))
        finally:
            _fork.reset(token)
    return results


def infer(func: _AnyFunc, /, *params: str | int) -> str:
    """Infer the ``optype`` protocol(s) required of ``func``'s parameters.

    Pass parameter names or positions to report only those parameters.

    >>> print(infer(lambda x: x + 1))
    [R](x: CanAdd[Literal[1], R]) -> R
    """
    parameters = _parameters(func)
    if any(p.kind in _VARIADIC for p in parameters.values()):
        raise NotImplementedError("variadic parameters")

    names = list(parameters)
    selected = _select(params, names)
    optional = frozenset(
        name for name, p in parameters.items() if p.default is not Parameter.empty
    )

    spies = {name: _SpyObject() for name in names}
    args: list[_SpyObject] = []
    kwds: dict[str, _SpyObject] = {}
    for name, param in parameters.items():
        if param.kind is Parameter.KEYWORD_ONLY:
            kwds[name] = spies[name]
        else:
            args.append(spies[name])
    results = _explore(func, args, kwds)

    sig1 = _Renderer(names, selected, spies, results, optional).render()
    _reflect(list(spies.values()), results)
    sig2 = _Renderer(names, selected, spies, results, optional).render()
    return "\n".join(dict.fromkeys((sig1, sig2)))
