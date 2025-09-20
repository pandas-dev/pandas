"""Shared definitions used by different parts of type checker."""

from __future__ import annotations

from abc import abstractmethod
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from typing import NamedTuple, overload

from mypy_extensions import trait

from mypy.errorcodes import ErrorCode
from mypy.errors import ErrorWatcher
from mypy.message_registry import ErrorMessage
from mypy.nodes import (
    ArgKind,
    Context,
    Expression,
    FuncItem,
    LambdaExpr,
    MypyFile,
    Node,
    RefExpr,
    SymbolNode,
    TypeInfo,
    Var,
)
from mypy.plugin import CheckerPluginInterface, Plugin
from mypy.types import (
    CallableType,
    Instance,
    LiteralValue,
    Overloaded,
    PartialType,
    TupleType,
    Type,
    TypedDictType,
    TypeType,
)
from mypy.typevars import fill_typevars


# An object that represents either a precise type or a type with an upper bound;
# it is important for correct type inference with isinstance.
class TypeRange(NamedTuple):
    item: Type
    is_upper_bound: bool  # False => precise type


@trait
class ExpressionCheckerSharedApi:
    @abstractmethod
    def accept(
        self,
        node: Expression,
        type_context: Type | None = None,
        allow_none_return: bool = False,
        always_allow_any: bool = False,
        is_callee: bool = False,
    ) -> Type:
        raise NotImplementedError

    @abstractmethod
    def analyze_ref_expr(self, e: RefExpr, lvalue: bool = False) -> Type:
        raise NotImplementedError

    @abstractmethod
    def check_call(
        self,
        callee: Type,
        args: list[Expression],
        arg_kinds: list[ArgKind],
        context: Context,
        arg_names: Sequence[str | None] | None = None,
        callable_node: Expression | None = None,
        callable_name: str | None = None,
        object_type: Type | None = None,
        original_type: Type | None = None,
    ) -> tuple[Type, Type]:
        raise NotImplementedError

    @abstractmethod
    def transform_callee_type(
        self,
        callable_name: str | None,
        callee: Type,
        args: list[Expression],
        arg_kinds: list[ArgKind],
        context: Context,
        arg_names: Sequence[str | None] | None = None,
        object_type: Type | None = None,
    ) -> Type:
        raise NotImplementedError

    @abstractmethod
    def method_fullname(self, object_type: Type, method_name: str) -> str | None:
        raise NotImplementedError

    @abstractmethod
    def check_method_call_by_name(
        self,
        method: str,
        base_type: Type,
        args: list[Expression],
        arg_kinds: list[ArgKind],
        context: Context,
        original_type: Type | None = None,
    ) -> tuple[Type, Type]:
        raise NotImplementedError

    @abstractmethod
    def visit_typeddict_index_expr(
        self, td_type: TypedDictType, index: Expression, setitem: bool = False
    ) -> tuple[Type, set[str]]:
        raise NotImplementedError

    @abstractmethod
    def infer_literal_expr_type(self, value: LiteralValue, fallback_name: str) -> Type:
        raise NotImplementedError

    @abstractmethod
    def analyze_static_reference(
        self,
        node: SymbolNode,
        ctx: Context,
        is_lvalue: bool,
        *,
        include_modules: bool = True,
        suppress_errors: bool = False,
    ) -> Type:
        raise NotImplementedError


@trait
class TypeCheckerSharedApi(CheckerPluginInterface):
    plugin: Plugin
    module_refs: set[str]
    scope: CheckerScope
    checking_missing_await: bool

    @property
    @abstractmethod
    def expr_checker(self) -> ExpressionCheckerSharedApi:
        raise NotImplementedError

    @abstractmethod
    def named_type(self, name: str) -> Instance:
        raise NotImplementedError

    @abstractmethod
    def lookup_typeinfo(self, fullname: str) -> TypeInfo:
        raise NotImplementedError

    @abstractmethod
    def lookup_type(self, node: Expression) -> Type:
        raise NotImplementedError

    @abstractmethod
    def handle_cannot_determine_type(self, name: str, context: Context) -> None:
        raise NotImplementedError

    @abstractmethod
    def handle_partial_var_type(
        self, typ: PartialType, is_lvalue: bool, node: Var, context: Context
    ) -> Type:
        raise NotImplementedError

    @overload
    @abstractmethod
    def check_subtype(
        self,
        subtype: Type,
        supertype: Type,
        context: Context,
        msg: str,
        subtype_label: str | None = None,
        supertype_label: str | None = None,
        *,
        notes: list[str] | None = None,
        code: ErrorCode | None = None,
        outer_context: Context | None = None,
    ) -> bool: ...

    @overload
    @abstractmethod
    def check_subtype(
        self,
        subtype: Type,
        supertype: Type,
        context: Context,
        msg: ErrorMessage,
        subtype_label: str | None = None,
        supertype_label: str | None = None,
        *,
        notes: list[str] | None = None,
        outer_context: Context | None = None,
    ) -> bool: ...

    # Unfortunately, mypyc doesn't support abstract overloads yet.
    @abstractmethod
    def check_subtype(
        self,
        subtype: Type,
        supertype: Type,
        context: Context,
        msg: str | ErrorMessage,
        subtype_label: str | None = None,
        supertype_label: str | None = None,
        *,
        notes: list[str] | None = None,
        code: ErrorCode | None = None,
        outer_context: Context | None = None,
    ) -> bool:
        raise NotImplementedError

    @abstractmethod
    def get_final_context(self) -> bool:
        raise NotImplementedError

    @overload
    @abstractmethod
    def conditional_types_with_intersection(
        self,
        expr_type: Type,
        type_ranges: list[TypeRange] | None,
        ctx: Context,
        default: None = None,
    ) -> tuple[Type | None, Type | None]: ...

    @overload
    @abstractmethod
    def conditional_types_with_intersection(
        self, expr_type: Type, type_ranges: list[TypeRange] | None, ctx: Context, default: Type
    ) -> tuple[Type, Type]: ...

    # Unfortunately, mypyc doesn't support abstract overloads yet.
    @abstractmethod
    def conditional_types_with_intersection(
        self,
        expr_type: Type,
        type_ranges: list[TypeRange] | None,
        ctx: Context,
        default: Type | None = None,
    ) -> tuple[Type | None, Type | None]:
        raise NotImplementedError

    @abstractmethod
    def check_deprecated(self, node: Node | None, context: Context) -> None:
        raise NotImplementedError

    @abstractmethod
    def warn_deprecated(self, node: Node | None, context: Context) -> None:
        raise NotImplementedError

    @abstractmethod
    def warn_deprecated_overload_item(
        self, node: Node | None, context: Context, *, target: Type, selftype: Type | None = None
    ) -> None:
        raise NotImplementedError

    @abstractmethod
    def type_is_iterable(self, type: Type) -> bool:
        raise NotImplementedError

    @abstractmethod
    def iterable_item_type(
        self, it: Instance | CallableType | TypeType | Overloaded, context: Context
    ) -> Type:
        raise NotImplementedError

    @abstractmethod
    @contextmanager
    def checking_await_set(self) -> Iterator[None]:
        raise NotImplementedError

    @abstractmethod
    def get_precise_awaitable_type(self, typ: Type, local_errors: ErrorWatcher) -> Type | None:
        raise NotImplementedError


class CheckerScope:
    # We keep two stacks combined, to maintain the relative order
    stack: list[TypeInfo | FuncItem | MypyFile]

    def __init__(self, module: MypyFile) -> None:
        self.stack = [module]

    def current_function(self) -> FuncItem | None:
        for e in reversed(self.stack):
            if isinstance(e, FuncItem):
                return e
        return None

    def top_level_function(self) -> FuncItem | None:
        """Return top-level non-lambda function."""
        for e in self.stack:
            if isinstance(e, FuncItem) and not isinstance(e, LambdaExpr):
                return e
        return None

    def active_class(self) -> TypeInfo | None:
        if isinstance(self.stack[-1], TypeInfo):
            return self.stack[-1]
        return None

    def enclosing_class(self, func: FuncItem | None = None) -> TypeInfo | None:
        """Is there a class *directly* enclosing this function?"""
        func = func or self.current_function()
        assert func, "This method must be called from inside a function"
        index = self.stack.index(func)
        assert index, "CheckerScope stack must always start with a module"
        enclosing = self.stack[index - 1]
        if isinstance(enclosing, TypeInfo):
            return enclosing
        return None

    def active_self_type(self) -> Instance | TupleType | None:
        """An instance or tuple type representing the current class.

        This returns None unless we are in class body or in a method.
        In particular, inside a function nested in method this returns None.
        """
        info = self.active_class()
        if not info and self.current_function():
            info = self.enclosing_class()
        if info:
            return fill_typevars(info)
        return None

    def current_self_type(self) -> Instance | TupleType | None:
        """Same as active_self_type() but handle functions nested in methods."""
        for item in reversed(self.stack):
            if isinstance(item, TypeInfo):
                return fill_typevars(item)
        return None

    @contextmanager
    def push_function(self, item: FuncItem) -> Iterator[None]:
        self.stack.append(item)
        yield
        self.stack.pop()

    @contextmanager
    def push_class(self, info: TypeInfo) -> Iterator[None]:
        self.stack.append(info)
        yield
        self.stack.pop()
