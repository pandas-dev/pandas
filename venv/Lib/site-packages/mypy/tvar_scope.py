from __future__ import annotations

from collections.abc import Callable
from typing import TypeAlias as _TypeAlias

from mypy.nodes import (
    Context,
    ParamSpecExpr,
    SymbolTableNode,
    TypeVarExpr,
    TypeVarLikeExpr,
    TypeVarTupleExpr,
)
from mypy.types import (
    AnyType,
    ParamSpecFlavor,
    ParamSpecType,
    TrivialSyntheticTypeTranslator,
    Type,
    TypeAliasType,
    TypeOfAny,
    TypeVarId,
    TypeVarLikeType,
    TypeVarTupleType,
    TypeVarType,
)

FailFunc: _TypeAlias = Callable[[str, Context], None]


class TypeVarLikeDefaultFixer(TrivialSyntheticTypeTranslator):
    """Set namespace for all TypeVarLikeTypes types."""

    def __init__(
        self,
        scope: TypeVarLikeScope,
        fail_func: FailFunc,
        source_tv: TypeVarLikeExpr,
        context: Context,
    ) -> None:
        self.scope = scope
        self.fail_func = fail_func
        self.source_tv = source_tv
        self.context = context
        super().__init__()

    def visit_type_var(self, t: TypeVarType) -> Type:
        existing = self.scope.get_binding(t.fullname)
        if existing is None:
            self._report_unbound_tvar(t)
            return AnyType(TypeOfAny.from_error)
        return existing

    def visit_param_spec(self, t: ParamSpecType) -> Type:
        existing = self.scope.get_binding(t.fullname)
        if existing is None:
            self._report_unbound_tvar(t)
            return AnyType(TypeOfAny.from_error)
        return existing

    def visit_type_var_tuple(self, t: TypeVarTupleType) -> Type:
        existing = self.scope.get_binding(t.fullname)
        if existing is None:
            self._report_unbound_tvar(t)
            return AnyType(TypeOfAny.from_error)
        return existing

    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        return t

    def _report_unbound_tvar(self, tvar: TypeVarLikeType) -> None:
        self.fail_func(
            f"Type variable {tvar.name} referenced in the default"
            f" of {self.source_tv.name} is unbound",
            self.context,
        )


class TypeVarLikeScope:
    """Scope that holds bindings for type variables and parameter specifications.

    Node fullname -> TypeVarLikeType.
    """

    def __init__(
        self,
        parent: TypeVarLikeScope | None = None,
        is_class_scope: bool = False,
        prohibited: TypeVarLikeScope | None = None,
        namespace: str = "",
    ) -> None:
        """Initializer for TypeVarLikeScope

        Parameters:
          parent: the outer scope for this scope
          is_class_scope: True if this represents a generic class
          prohibited: Type variables that aren't strictly in scope exactly,
                      but can't be bound because they're part of an outer class's scope.
        """
        self.scope: dict[str, TypeVarLikeType] = {}
        self.parent = parent
        self.func_id = 0
        self.class_id = 0
        self.is_class_scope = is_class_scope
        self.prohibited = prohibited
        self.namespace = namespace
        if parent is not None:
            self.func_id = parent.func_id
            self.class_id = parent.class_id

    def get_function_scope(self) -> TypeVarLikeScope | None:
        """Get the nearest parent that's a function scope, not a class scope"""
        it: TypeVarLikeScope | None = self
        while it is not None and it.is_class_scope:
            it = it.parent
        return it

    def allow_binding(self, fullname: str) -> bool:
        if fullname in self.scope:
            return False
        elif self.parent and not self.parent.allow_binding(fullname):
            return False
        elif self.prohibited and not self.prohibited.allow_binding(fullname):
            return False
        return True

    def method_frame(self, namespace: str) -> TypeVarLikeScope:
        """A new scope frame for binding a method"""
        return TypeVarLikeScope(self, False, None, namespace=namespace)

    def class_frame(self, namespace: str) -> TypeVarLikeScope:
        """A new scope frame for binding a class. Prohibits *this* class's tvars"""
        return TypeVarLikeScope(self.get_function_scope(), True, self, namespace=namespace)

    def new_unique_func_id(self) -> TypeVarId:
        """Used by plugin-like code that needs to make synthetic generic functions."""
        self.func_id -= 1
        return TypeVarId(self.func_id)

    def bind_new(
        self, name: str, tvar_expr: TypeVarLikeExpr, fail_func: FailFunc, context: Context
    ) -> TypeVarLikeType:
        if self.is_class_scope:
            self.class_id += 1
            i = self.class_id
        else:
            self.func_id -= 1
            i = self.func_id
        namespace = self.namespace

        # Defaults may reference other type variables. That is only valid when the
        # referenced variable is already in scope (textually precedes the definition we're
        # processing now).
        default = tvar_expr.default.accept(
            TypeVarLikeDefaultFixer(
                self, fail_func=fail_func, source_tv=tvar_expr, context=context
            )
        )

        if isinstance(tvar_expr, TypeVarExpr):
            tvar_def: TypeVarLikeType = TypeVarType(
                name=name,
                fullname=tvar_expr.fullname,
                id=TypeVarId(i, namespace=namespace),
                values=tvar_expr.values,
                upper_bound=tvar_expr.upper_bound,
                default=default,
                variance=tvar_expr.variance,
                line=tvar_expr.line,
                column=tvar_expr.column,
            )
        elif isinstance(tvar_expr, ParamSpecExpr):
            tvar_def = ParamSpecType(
                name=name,
                fullname=tvar_expr.fullname,
                id=TypeVarId(i, namespace=namespace),
                flavor=ParamSpecFlavor.BARE,
                upper_bound=tvar_expr.upper_bound,
                default=default,
                line=tvar_expr.line,
                column=tvar_expr.column,
            )
        elif isinstance(tvar_expr, TypeVarTupleExpr):
            tvar_def = TypeVarTupleType(
                name=name,
                fullname=tvar_expr.fullname,
                id=TypeVarId(i, namespace=namespace),
                upper_bound=tvar_expr.upper_bound,
                tuple_fallback=tvar_expr.tuple_fallback,
                default=default,
                line=tvar_expr.line,
                column=tvar_expr.column,
            )
        else:
            assert False
        self.scope[tvar_expr.fullname] = tvar_def
        return tvar_def

    def bind_existing(self, tvar_def: TypeVarLikeType) -> None:
        self.scope[tvar_def.fullname] = tvar_def

    def get_binding(self, item: str | SymbolTableNode) -> TypeVarLikeType | None:
        fullname = item.fullname if isinstance(item, SymbolTableNode) else item
        assert fullname
        if fullname in self.scope:
            return self.scope[fullname]
        elif self.parent is not None:
            return self.parent.get_binding(fullname)
        else:
            return None

    def __str__(self) -> str:
        me = ", ".join(f"{k}: {v.name}`{v.id}" for k, v in self.scope.items())
        if self.parent is None:
            return me
        return f"{self.parent} <- {me}"
