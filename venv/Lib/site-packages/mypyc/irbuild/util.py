"""Various utilities that don't depend on other modules in mypyc.irbuild."""

from __future__ import annotations

from typing import Any, Final, Literal, TypedDict
from typing_extensions import NotRequired

from mypy.nodes import (
    ARG_NAMED,
    ARG_NAMED_OPT,
    ARG_OPT,
    ARG_POS,
    GDEF,
    ArgKind,
    BytesExpr,
    CallExpr,
    ClassDef,
    Decorator,
    Expression,
    FloatExpr,
    FuncDef,
    IntExpr,
    NameExpr,
    OverloadedFuncDef,
    RefExpr,
    StrExpr,
    TupleExpr,
    UnaryExpr,
    Var,
)
from mypy.semanal import refers_to_fullname
from mypy.types import FINAL_DECORATOR_NAMES
from mypyc.errors import Errors

MYPYC_ATTRS: Final[frozenset[MypycAttr]] = frozenset(
    ["native_class", "allow_interpreted_subclasses", "serializable", "free_list_len", "acyclic"]
)

DATACLASS_DECORATORS: Final = frozenset(["dataclasses.dataclass", "attr.s", "attr.attrs"])


MypycAttr = Literal[
    "native_class", "allow_interpreted_subclasses", "serializable", "free_list_len", "acyclic"
]


class MypycAttrs(TypedDict):
    native_class: NotRequired[bool]
    allow_interpreted_subclasses: NotRequired[bool]
    serializable: NotRequired[bool]
    free_list_len: NotRequired[int]
    acyclic: NotRequired[bool]


def is_final_decorator(d: Expression) -> bool:
    return refers_to_fullname(d, FINAL_DECORATOR_NAMES)


def is_trait_decorator(d: Expression) -> bool:
    return isinstance(d, RefExpr) and d.fullname == "mypy_extensions.trait"


def is_trait(cdef: ClassDef) -> bool:
    return any(is_trait_decorator(d) for d in cdef.decorators) or cdef.info.is_protocol


def dataclass_decorator_type(d: Expression) -> str | None:
    if isinstance(d, RefExpr) and d.fullname in DATACLASS_DECORATORS:
        return d.fullname.split(".")[0]
    elif (
        isinstance(d, CallExpr)
        and isinstance(d.callee, RefExpr)
        and d.callee.fullname in DATACLASS_DECORATORS
    ):
        name = d.callee.fullname.split(".")[0]
        if name == "attr" and "auto_attribs" in d.arg_names:
            # Note: the mypy attrs plugin checks that the value of auto_attribs is
            # not computed at runtime, so we don't need to perform that check here
            auto = d.args[d.arg_names.index("auto_attribs")]
            if isinstance(auto, NameExpr) and auto.name == "True":
                return "attr-auto"
        return name
    else:
        return None


def is_dataclass_decorator(d: Expression) -> bool:
    return dataclass_decorator_type(d) is not None


def is_dataclass(cdef: ClassDef) -> bool:
    return any(is_dataclass_decorator(d) for d in cdef.decorators)


# The string values returned by this function are inspected in
# mypyc/lib-rt/misc_ops.c:CPyDataclass_SleightOfHand(...).
def dataclass_type(cdef: ClassDef) -> str | None:
    for d in cdef.decorators:
        typ = dataclass_decorator_type(d)
        if typ is not None:
            return typ
    return None


def get_mypyc_attr_literal(e: Expression) -> Any:
    """Convert an expression from a mypyc_attr decorator to a value.

    Supports a pretty limited range."""
    if isinstance(e, (StrExpr, IntExpr, FloatExpr)):
        return e.value
    elif isinstance(e, RefExpr) and e.fullname == "builtins.True":
        return True
    elif isinstance(e, RefExpr) and e.fullname == "builtins.False":
        return False
    elif isinstance(e, RefExpr) and e.fullname == "builtins.None":
        return None
    elif isinstance(e, IntExpr):
        return e.value
    return NotImplemented


def get_mypyc_attr_call(d: Expression) -> CallExpr | None:
    """Check if an expression is a call to mypyc_attr and return it if so."""
    if (
        isinstance(d, CallExpr)
        and isinstance(d.callee, RefExpr)
        and d.callee.fullname == "mypy_extensions.mypyc_attr"
    ):
        return d
    return None


def get_mypyc_attrs(
    stmt: ClassDef | Decorator, path: str, errors: Errors
) -> tuple[MypycAttrs, dict[MypycAttr, int]]:
    """Collect all the mypyc_attr attributes on a class definition or a function."""
    attrs: MypycAttrs = {}
    lines: dict[MypycAttr, int] = {}

    def set_mypyc_attr(key: str, value: Any, line: int) -> None:
        if key in MYPYC_ATTRS:
            attrs[key] = value
            lines[key] = line
        else:
            errors.error(f'"{key}" is not a supported "mypyc_attr"', path, line)
            supported_keys = '", "'.join(sorted(MYPYC_ATTRS))
            errors.note(f'supported keys: "{supported_keys}"', path, line)

    for dec in stmt.decorators:
        if d := get_mypyc_attr_call(dec):
            line = d.line
            for name, arg in zip(d.arg_names, d.args):
                if name is None:
                    if isinstance(arg, StrExpr):
                        set_mypyc_attr(arg.value, True, line)
                    else:
                        errors.error(
                            'All "mypyc_attr" positional arguments must be string literals.',
                            path,
                            line,
                        )
                else:
                    arg_value = get_mypyc_attr_literal(arg)
                    set_mypyc_attr(name, arg_value, line)

    return attrs, lines


def is_extension_class(path: str, cdef: ClassDef, errors: Errors) -> bool:
    # Check for @mypyc_attr(native_class=True/False) decorator.
    explicit_native_class = get_explicit_native_class(path, cdef, errors)

    # Classes with native_class=False are explicitly marked as non extension.
    if explicit_native_class is False:
        return False

    implicit_extension_class, reason = is_implicit_extension_class(cdef)

    # Classes with native_class=True should be extension classes, but they might
    # not be able to be due to other reasons. Print an error in that case.
    if explicit_native_class is True and not implicit_extension_class:
        errors.error(
            f"Class is marked as native_class=True but it can't be a native class. {reason}",
            path,
            cdef.line,
        )

    return implicit_extension_class


def get_explicit_native_class(path: str, cdef: ClassDef, errors: Errors) -> bool | None:
    """Return value of @mypyc_attr(native_class=True/False) decorator.

    Look for a @mypyc_attr decorator with native_class=True/False and return
    the value assigned or None if it doesn't exist. Other values are an error.
    """

    for d in cdef.decorators:
        mypyc_attr_call = get_mypyc_attr_call(d)
        if not mypyc_attr_call:
            continue

        for i, name in enumerate(mypyc_attr_call.arg_names):
            if name != "native_class":
                continue

            arg = mypyc_attr_call.args[i]
            if not isinstance(arg, NameExpr):
                errors.error("native_class must be used with True or False only", path, cdef.line)
                return None

            if arg.name == "False":
                return False
            elif arg.name == "True":
                return True
            else:
                errors.error("native_class must be used with True or False only", path, cdef.line)
                return None
    return None


def is_implicit_extension_class(cdef: ClassDef) -> tuple[bool, str]:
    """Check if class can be extension class and return a user-friendly reason it can't be one."""

    for d in cdef.decorators:
        if (
            not is_trait_decorator(d)
            and not is_dataclass_decorator(d)
            and not get_mypyc_attr_call(d)
            and not is_final_decorator(d)
        ):
            return (
                False,
                "Classes that have decorators other than supported decorators"
                " can't be native classes.",
            )

    if cdef.info.typeddict_type:
        return False, "TypedDict classes can't be native classes."
    if cdef.info.is_named_tuple:
        return False, "NamedTuple classes can't be native classes."
    if cdef.info.metaclass_type and cdef.info.metaclass_type.type.fullname not in (
        "abc.ABCMeta",
        "typing.TypingMeta",
        "typing.GenericMeta",
    ):
        return (
            False,
            "Classes with a metaclass other than ABCMeta, TypingMeta or"
            " GenericMeta can't be native classes.",
        )
    return True, ""


def get_func_def(op: FuncDef | Decorator | OverloadedFuncDef) -> FuncDef:
    if isinstance(op, OverloadedFuncDef):
        assert op.impl
        op = op.impl
    if isinstance(op, Decorator):
        op = op.func
    return op


def concrete_arg_kind(kind: ArgKind) -> ArgKind:
    """Find the concrete version of an arg kind that is being passed."""
    if kind == ARG_OPT:
        return ARG_POS
    elif kind == ARG_NAMED_OPT:
        return ARG_NAMED
    else:
        return kind


def is_constant(e: Expression) -> bool:
    """Check whether we allow an expression to appear as a default value.

    We don't currently properly support storing the evaluated
    values for default arguments and default attribute values, so
    we restrict what expressions we allow.  We allow literals of
    primitives types, None, and references to Final global
    variables.
    """
    return (
        isinstance(e, (StrExpr, BytesExpr, IntExpr, FloatExpr))
        or (isinstance(e, UnaryExpr) and e.op == "-" and isinstance(e.expr, (IntExpr, FloatExpr)))
        or (isinstance(e, TupleExpr) and all(is_constant(e) for e in e.items))
        or (
            isinstance(e, RefExpr)
            and e.kind == GDEF
            and (
                e.fullname in ("builtins.True", "builtins.False", "builtins.None")
                or (isinstance(e.node, Var) and e.node.is_final)
            )
        )
    )


def bytes_from_str(value: str) -> bytes:
    """Convert a string representing bytes into actual bytes.

    This is needed because the literal characters of BytesExpr (the
    characters inside b'') are stored in BytesExpr.value, whose type is
    'str' not 'bytes'.
    """
    return bytes(value, "utf8").decode("unicode-escape").encode("raw-unicode-escape")
