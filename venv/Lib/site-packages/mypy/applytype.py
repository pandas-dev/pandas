from __future__ import annotations

from typing import Callable, Iterable, Sequence

import mypy.subtypes
from mypy.erasetype import erase_typevars
from mypy.expandtype import expand_type
from mypy.nodes import Context, TypeInfo
from mypy.type_visitor import TypeTranslator
from mypy.typeops import get_all_type_vars
from mypy.types import (
    AnyType,
    CallableType,
    Instance,
    Parameters,
    ParamSpecFlavor,
    ParamSpecType,
    PartialType,
    ProperType,
    Type,
    TypeAliasType,
    TypeVarId,
    TypeVarLikeType,
    TypeVarTupleType,
    TypeVarType,
    UninhabitedType,
    UnpackType,
    get_proper_type,
    remove_dups,
)


def get_target_type(
    tvar: TypeVarLikeType,
    type: Type,
    callable: CallableType,
    report_incompatible_typevar_value: Callable[[CallableType, Type, str, Context], None],
    context: Context,
    skip_unsatisfied: bool,
) -> Type | None:
    p_type = get_proper_type(type)
    if isinstance(p_type, UninhabitedType) and tvar.has_default():
        return tvar.default
    if isinstance(tvar, ParamSpecType):
        return type
    if isinstance(tvar, TypeVarTupleType):
        return type
    assert isinstance(tvar, TypeVarType)
    values = tvar.values
    if values:
        if isinstance(p_type, AnyType):
            return type
        if isinstance(p_type, TypeVarType) and p_type.values:
            # Allow substituting T1 for T if every allowed value of T1
            # is also a legal value of T.
            if all(any(mypy.subtypes.is_same_type(v, v1) for v in values) for v1 in p_type.values):
                return type
        matching = []
        for value in values:
            if mypy.subtypes.is_subtype(type, value):
                matching.append(value)
        if matching:
            best = matching[0]
            # If there are more than one matching value, we select the narrowest
            for match in matching[1:]:
                if mypy.subtypes.is_subtype(match, best):
                    best = match
            return best
        if skip_unsatisfied:
            return None
        report_incompatible_typevar_value(callable, type, tvar.name, context)
    else:
        upper_bound = tvar.upper_bound
        if tvar.name == "Self":
            # Internally constructed Self-types contain class type variables in upper bound,
            # so we need to erase them to avoid false positives. This is safe because we do
            # not support type variables in upper bounds of user defined types.
            upper_bound = erase_typevars(upper_bound)
        if not mypy.subtypes.is_subtype(type, upper_bound):
            if skip_unsatisfied:
                return None
            report_incompatible_typevar_value(callable, type, tvar.name, context)
    return type


def apply_generic_arguments(
    callable: CallableType,
    orig_types: Sequence[Type | None],
    report_incompatible_typevar_value: Callable[[CallableType, Type, str, Context], None],
    context: Context,
    skip_unsatisfied: bool = False,
) -> CallableType:
    """Apply generic type arguments to a callable type.

    For example, applying [int] to 'def [T] (T) -> T' results in
    'def (int) -> int'.

    Note that each type can be None; in this case, it will not be applied.

    If `skip_unsatisfied` is True, then just skip the types that don't satisfy type variable
    bound or constraints, instead of giving an error.
    """
    tvars = callable.variables
    assert len(orig_types) <= len(tvars)
    # Check that inferred type variable values are compatible with allowed
    # values and bounds.  Also, promote subtype values to allowed values.
    # Create a map from type variable id to target type.
    id_to_type: dict[TypeVarId, Type] = {}

    for tvar, type in zip(tvars, orig_types):
        assert not isinstance(type, PartialType), "Internal error: must never apply partial type"
        if type is None:
            continue

        target_type = get_target_type(
            tvar, type, callable, report_incompatible_typevar_value, context, skip_unsatisfied
        )
        if target_type is not None:
            id_to_type[tvar.id] = target_type

    # TODO: validate arg_kinds/arg_names for ParamSpec and TypeVarTuple replacements,
    # not just type variable bounds above.
    param_spec = callable.param_spec()
    if param_spec is not None:
        nt = id_to_type.get(param_spec.id)
        if nt is not None:
            # ParamSpec expansion is special-cased, so we need to always expand callable
            # as a whole, not expanding arguments individually.
            callable = expand_type(callable, id_to_type)
            assert isinstance(callable, CallableType)
            return callable.copy_modified(
                variables=[tv for tv in tvars if tv.id not in id_to_type]
            )

    # Apply arguments to argument types.
    var_arg = callable.var_arg()
    if var_arg is not None and isinstance(var_arg.typ, UnpackType):
        # Same as for ParamSpec, callable with variadic types needs to be expanded as a whole.
        callable = expand_type(callable, id_to_type)
        assert isinstance(callable, CallableType)
        return callable.copy_modified(variables=[tv for tv in tvars if tv.id not in id_to_type])
    else:
        callable = callable.copy_modified(
            arg_types=[expand_type(at, id_to_type) for at in callable.arg_types]
        )

    # Apply arguments to TypeGuard and TypeIs if any.
    if callable.type_guard is not None:
        type_guard = expand_type(callable.type_guard, id_to_type)
    else:
        type_guard = None
    if callable.type_is is not None:
        type_is = expand_type(callable.type_is, id_to_type)
    else:
        type_is = None

    # The callable may retain some type vars if only some were applied.
    # TODO: move apply_poly() logic here when new inference
    # becomes universally used (i.e. in all passes + in unification).
    # With this new logic we can actually *add* some new free variables.
    remaining_tvars: list[TypeVarLikeType] = []
    for tv in tvars:
        if tv.id in id_to_type:
            continue
        if not tv.has_default():
            remaining_tvars.append(tv)
            continue
        # TypeVarLike isn't in id_to_type mapping.
        # Only expand the TypeVar default here.
        typ = expand_type(tv, id_to_type)
        assert isinstance(typ, TypeVarLikeType)
        remaining_tvars.append(typ)

    return callable.copy_modified(
        ret_type=expand_type(callable.ret_type, id_to_type),
        variables=remaining_tvars,
        type_guard=type_guard,
        type_is=type_is,
    )


def apply_poly(tp: CallableType, poly_tvars: Sequence[TypeVarLikeType]) -> CallableType | None:
    """Make free type variables generic in the type if possible.

    This will translate the type `tp` while trying to create valid bindings for
    type variables `poly_tvars` while traversing the type. This follows the same rules
    as we do during semantic analysis phase, examples:
      * Callable[Callable[[T], T], T] -> def [T] (def (T) -> T) -> T
      * Callable[[], Callable[[T], T]] -> def () -> def [T] (T -> T)
      * List[T] -> None (not possible)
    """
    try:
        return tp.copy_modified(
            arg_types=[t.accept(PolyTranslator(poly_tvars)) for t in tp.arg_types],
            ret_type=tp.ret_type.accept(PolyTranslator(poly_tvars)),
            variables=[],
        )
    except PolyTranslationError:
        return None


class PolyTranslationError(Exception):
    pass


class PolyTranslator(TypeTranslator):
    """Make free type variables generic in the type if possible.

    See docstring for apply_poly() for details.
    """

    def __init__(
        self,
        poly_tvars: Iterable[TypeVarLikeType],
        bound_tvars: frozenset[TypeVarLikeType] = frozenset(),
        seen_aliases: frozenset[TypeInfo] = frozenset(),
    ) -> None:
        super().__init__()
        self.poly_tvars = set(poly_tvars)
        # This is a simplified version of TypeVarScope used during semantic analysis.
        self.bound_tvars = bound_tvars
        self.seen_aliases = seen_aliases

    def collect_vars(self, t: CallableType | Parameters) -> list[TypeVarLikeType]:
        found_vars = []
        for arg in t.arg_types:
            for tv in get_all_type_vars(arg):
                if isinstance(tv, ParamSpecType):
                    normalized: TypeVarLikeType = tv.copy_modified(
                        flavor=ParamSpecFlavor.BARE, prefix=Parameters([], [], [])
                    )
                else:
                    normalized = tv
                if normalized in self.poly_tvars and normalized not in self.bound_tvars:
                    found_vars.append(normalized)
        return remove_dups(found_vars)

    def visit_callable_type(self, t: CallableType) -> Type:
        found_vars = self.collect_vars(t)
        self.bound_tvars |= set(found_vars)
        result = super().visit_callable_type(t)
        self.bound_tvars -= set(found_vars)

        assert isinstance(result, ProperType) and isinstance(result, CallableType)
        result.variables = list(result.variables) + found_vars
        return result

    def visit_type_var(self, t: TypeVarType) -> Type:
        if t in self.poly_tvars and t not in self.bound_tvars:
            raise PolyTranslationError()
        return super().visit_type_var(t)

    def visit_param_spec(self, t: ParamSpecType) -> Type:
        if t in self.poly_tvars and t not in self.bound_tvars:
            raise PolyTranslationError()
        return super().visit_param_spec(t)

    def visit_type_var_tuple(self, t: TypeVarTupleType) -> Type:
        if t in self.poly_tvars and t not in self.bound_tvars:
            raise PolyTranslationError()
        return super().visit_type_var_tuple(t)

    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        if not t.args:
            return t.copy_modified()
        if not t.is_recursive:
            return get_proper_type(t).accept(self)
        # We can't handle polymorphic application for recursive generic aliases
        # without risking an infinite recursion, just give up for now.
        raise PolyTranslationError()

    def visit_instance(self, t: Instance) -> Type:
        if t.type.has_param_spec_type:
            # We need this special-casing to preserve the possibility to store a
            # generic function in an instance type. Things like
            #     forall T . Foo[[x: T], T]
            # are not really expressible in current type system, but this looks like
            # a useful feature, so let's keep it.
            param_spec_index = next(
                i for (i, tv) in enumerate(t.type.defn.type_vars) if isinstance(tv, ParamSpecType)
            )
            p = get_proper_type(t.args[param_spec_index])
            if isinstance(p, Parameters):
                found_vars = self.collect_vars(p)
                self.bound_tvars |= set(found_vars)
                new_args = [a.accept(self) for a in t.args]
                self.bound_tvars -= set(found_vars)

                repl = new_args[param_spec_index]
                assert isinstance(repl, ProperType) and isinstance(repl, Parameters)
                repl.variables = list(repl.variables) + list(found_vars)
                return t.copy_modified(args=new_args)
        # There is the same problem with callback protocols as with aliases
        # (callback protocols are essentially more flexible aliases to callables).
        if t.args and t.type.is_protocol and t.type.protocol_members == ["__call__"]:
            if t.type in self.seen_aliases:
                raise PolyTranslationError()
            call = mypy.subtypes.find_member("__call__", t, t, is_operator=True)
            assert call is not None
            return call.accept(
                PolyTranslator(self.poly_tvars, self.bound_tvars, self.seen_aliases | {t.type})
            )
        return super().visit_instance(t)
