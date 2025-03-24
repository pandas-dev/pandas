"""Plugin for supporting the functools standard library module."""

from __future__ import annotations

from typing import Final, NamedTuple

import mypy.checker
import mypy.plugin
import mypy.semanal
from mypy.argmap import map_actuals_to_formals
from mypy.nodes import ARG_POS, ARG_STAR2, ArgKind, Argument, CallExpr, FuncItem, Var
from mypy.plugins.common import add_method_to_class
from mypy.typeops import get_all_type_vars
from mypy.types import (
    AnyType,
    CallableType,
    Instance,
    Overloaded,
    Type,
    TypeOfAny,
    TypeVarType,
    UnboundType,
    UnionType,
    get_proper_type,
)

functools_total_ordering_makers: Final = {"functools.total_ordering"}

_ORDERING_METHODS: Final = {"__lt__", "__le__", "__gt__", "__ge__"}

PARTIAL: Final = "functools.partial"


class _MethodInfo(NamedTuple):
    is_static: bool
    type: CallableType


def functools_total_ordering_maker_callback(
    ctx: mypy.plugin.ClassDefContext, auto_attribs_default: bool = False
) -> bool:
    """Add dunder methods to classes decorated with functools.total_ordering."""
    comparison_methods = _analyze_class(ctx)
    if not comparison_methods:
        ctx.api.fail(
            'No ordering operation defined when using "functools.total_ordering": < > <= >=',
            ctx.reason,
        )
        return True

    # prefer __lt__ to __le__ to __gt__ to __ge__
    root = max(comparison_methods, key=lambda k: (comparison_methods[k] is None, k))
    root_method = comparison_methods[root]
    if not root_method:
        # None of the defined comparison methods can be analysed
        return True

    other_type = _find_other_type(root_method)
    bool_type = ctx.api.named_type("builtins.bool")
    ret_type: Type = bool_type
    if root_method.type.ret_type != ctx.api.named_type("builtins.bool"):
        proper_ret_type = get_proper_type(root_method.type.ret_type)
        if not (
            isinstance(proper_ret_type, UnboundType)
            and proper_ret_type.name.split(".")[-1] == "bool"
        ):
            ret_type = AnyType(TypeOfAny.implementation_artifact)
    for additional_op in _ORDERING_METHODS:
        # Either the method is not implemented
        # or has an unknown signature that we can now extrapolate.
        if not comparison_methods.get(additional_op):
            args = [Argument(Var("other", other_type), other_type, None, ARG_POS)]
            add_method_to_class(ctx.api, ctx.cls, additional_op, args, ret_type)

    return True


def _find_other_type(method: _MethodInfo) -> Type:
    """Find the type of the ``other`` argument in a comparison method."""
    first_arg_pos = 0 if method.is_static else 1
    cur_pos_arg = 0
    other_arg = None
    for arg_kind, arg_type in zip(method.type.arg_kinds, method.type.arg_types):
        if arg_kind.is_positional():
            if cur_pos_arg == first_arg_pos:
                other_arg = arg_type
                break

            cur_pos_arg += 1
        elif arg_kind != ARG_STAR2:
            other_arg = arg_type
            break

    if other_arg is None:
        return AnyType(TypeOfAny.implementation_artifact)

    return other_arg


def _analyze_class(ctx: mypy.plugin.ClassDefContext) -> dict[str, _MethodInfo | None]:
    """Analyze the class body, its parents, and return the comparison methods found."""
    # Traverse the MRO and collect ordering methods.
    comparison_methods: dict[str, _MethodInfo | None] = {}
    # Skip object because total_ordering does not use methods from object
    for cls in ctx.cls.info.mro[:-1]:
        for name in _ORDERING_METHODS:
            if name in cls.names and name not in comparison_methods:
                node = cls.names[name].node
                if isinstance(node, FuncItem) and isinstance(node.type, CallableType):
                    comparison_methods[name] = _MethodInfo(node.is_static, node.type)
                    continue

                if isinstance(node, Var):
                    proper_type = get_proper_type(node.type)
                    if isinstance(proper_type, CallableType):
                        comparison_methods[name] = _MethodInfo(node.is_staticmethod, proper_type)
                        continue

                comparison_methods[name] = None

    return comparison_methods


def partial_new_callback(ctx: mypy.plugin.FunctionContext) -> Type:
    """Infer a more precise return type for functools.partial"""
    if not isinstance(ctx.api, mypy.checker.TypeChecker):  # use internals
        return ctx.default_return_type
    if len(ctx.arg_types) != 3:  # fn, *args, **kwargs
        return ctx.default_return_type
    if len(ctx.arg_types[0]) != 1:
        return ctx.default_return_type

    if isinstance(get_proper_type(ctx.arg_types[0][0]), Overloaded):
        # TODO: handle overloads, just fall back to whatever the non-plugin code does
        return ctx.default_return_type
    return handle_partial_with_callee(ctx, callee=ctx.arg_types[0][0])


def handle_partial_with_callee(ctx: mypy.plugin.FunctionContext, callee: Type) -> Type:
    if not isinstance(ctx.api, mypy.checker.TypeChecker):  # use internals
        return ctx.default_return_type

    if isinstance(callee_proper := get_proper_type(callee), UnionType):
        return UnionType.make_union(
            [handle_partial_with_callee(ctx, item) for item in callee_proper.items]
        )

    fn_type = ctx.api.extract_callable_type(callee, ctx=ctx.default_return_type)
    if fn_type is None:
        return ctx.default_return_type

    # We must normalize from the start to have coherent view together with TypeChecker.
    fn_type = fn_type.with_unpacked_kwargs().with_normalized_var_args()

    last_context = ctx.api.type_context[-1]
    if not fn_type.is_type_obj():
        # We wrap the return type to get use of a possible type context provided by caller.
        # We cannot do this in case of class objects, since otherwise the plugin may get
        # falsely triggered when evaluating the constructed call itself.
        ret_type: Type = ctx.api.named_generic_type(PARTIAL, [fn_type.ret_type])
        wrapped_return = True
    else:
        ret_type = fn_type.ret_type
        # Instead, for class objects we ignore any type context to avoid spurious errors,
        # since the type context will be partial[X] etc., not X.
        ctx.api.type_context[-1] = None
        wrapped_return = False

    # Flatten actual to formal mapping, since this is what check_call() expects.
    actual_args = []
    actual_arg_kinds = []
    actual_arg_names = []
    actual_types = []
    seen_args = set()
    for i, param in enumerate(ctx.args[1:], start=1):
        for j, a in enumerate(param):
            if a in seen_args:
                # Same actual arg can map to multiple formals, but we need to include
                # each one only once.
                continue
            # Here we rely on the fact that expressions are essentially immutable, so
            # they can be compared by identity.
            seen_args.add(a)
            actual_args.append(a)
            actual_arg_kinds.append(ctx.arg_kinds[i][j])
            actual_arg_names.append(ctx.arg_names[i][j])
            actual_types.append(ctx.arg_types[i][j])

    formal_to_actual = map_actuals_to_formals(
        actual_kinds=actual_arg_kinds,
        actual_names=actual_arg_names,
        formal_kinds=fn_type.arg_kinds,
        formal_names=fn_type.arg_names,
        actual_arg_type=lambda i: actual_types[i],
    )

    # We need to remove any type variables that appear only in formals that have
    # no actuals, to avoid eagerly binding them in check_call() below.
    can_infer_ids = set()
    for i, arg_type in enumerate(fn_type.arg_types):
        if not formal_to_actual[i]:
            continue
        can_infer_ids.update({tv.id for tv in get_all_type_vars(arg_type)})

    defaulted = fn_type.copy_modified(
        arg_kinds=[
            (
                ArgKind.ARG_OPT
                if k == ArgKind.ARG_POS
                else (ArgKind.ARG_NAMED_OPT if k == ArgKind.ARG_NAMED else k)
            )
            for k in fn_type.arg_kinds
        ],
        ret_type=ret_type,
        variables=[
            tv
            for tv in fn_type.variables
            # Keep TypeVarTuple/ParamSpec to avoid spurious errors on empty args.
            if tv.id in can_infer_ids or not isinstance(tv, TypeVarType)
        ],
    )
    if defaulted.line < 0:
        # Make up a line number if we don't have one
        defaulted.set_line(ctx.default_return_type)

    # Create a valid context for various ad-hoc inspections in check_call().
    call_expr = CallExpr(
        callee=ctx.args[0][0],
        args=actual_args,
        arg_kinds=actual_arg_kinds,
        arg_names=actual_arg_names,
        analyzed=ctx.context.analyzed if isinstance(ctx.context, CallExpr) else None,
    )
    call_expr.set_line(ctx.context)

    _, bound = ctx.api.expr_checker.check_call(
        callee=defaulted,
        args=actual_args,
        arg_kinds=actual_arg_kinds,
        arg_names=actual_arg_names,
        context=call_expr,
    )
    if not wrapped_return:
        # Restore previously ignored context.
        ctx.api.type_context[-1] = last_context

    bound = get_proper_type(bound)
    if not isinstance(bound, CallableType):
        return ctx.default_return_type

    if wrapped_return:
        # Reverse the wrapping we did above.
        ret_type = get_proper_type(bound.ret_type)
        if not isinstance(ret_type, Instance) or ret_type.type.fullname != PARTIAL:
            return ctx.default_return_type
        bound = bound.copy_modified(ret_type=ret_type.args[0])

    partial_kinds = []
    partial_types = []
    partial_names = []
    # We need to fully apply any positional arguments (they cannot be respecified)
    # However, keyword arguments can be respecified, so just give them a default
    for i, actuals in enumerate(formal_to_actual):
        if len(bound.arg_types) == len(fn_type.arg_types):
            arg_type = bound.arg_types[i]
            if not mypy.checker.is_valid_inferred_type(arg_type):
                arg_type = fn_type.arg_types[i]  # bit of a hack
        else:
            # TODO: I assume that bound and fn_type have the same arguments. It appears this isn't
            # true when PEP 646 things are happening. See testFunctoolsPartialTypeVarTuple
            arg_type = fn_type.arg_types[i]

        if not actuals or fn_type.arg_kinds[i] in (ArgKind.ARG_STAR, ArgKind.ARG_STAR2):
            partial_kinds.append(fn_type.arg_kinds[i])
            partial_types.append(arg_type)
            partial_names.append(fn_type.arg_names[i])
        else:
            assert actuals
            if any(actual_arg_kinds[j] in (ArgKind.ARG_POS, ArgKind.ARG_STAR) for j in actuals):
                # Don't add params for arguments passed positionally
                continue
            # Add defaulted params for arguments passed via keyword
            kind = actual_arg_kinds[actuals[0]]
            if kind == ArgKind.ARG_NAMED or kind == ArgKind.ARG_STAR2:
                kind = ArgKind.ARG_NAMED_OPT
            partial_kinds.append(kind)
            partial_types.append(arg_type)
            partial_names.append(fn_type.arg_names[i])

    ret_type = bound.ret_type
    if not mypy.checker.is_valid_inferred_type(ret_type):
        ret_type = fn_type.ret_type  # same kind of hack as above

    partially_applied = fn_type.copy_modified(
        arg_types=partial_types,
        arg_kinds=partial_kinds,
        arg_names=partial_names,
        ret_type=ret_type,
    )

    ret = ctx.api.named_generic_type(PARTIAL, [ret_type])
    ret = ret.copy_with_extra_attr("__mypy_partial", partially_applied)
    return ret


def partial_call_callback(ctx: mypy.plugin.MethodContext) -> Type:
    """Infer a more precise return type for functools.partial.__call__."""
    if (
        not isinstance(ctx.api, mypy.checker.TypeChecker)  # use internals
        or not isinstance(ctx.type, Instance)
        or ctx.type.type.fullname != PARTIAL
        or not ctx.type.extra_attrs
        or "__mypy_partial" not in ctx.type.extra_attrs.attrs
    ):
        return ctx.default_return_type

    partial_type = ctx.type.extra_attrs.attrs["__mypy_partial"]
    if len(ctx.arg_types) != 2:  # *args, **kwargs
        return ctx.default_return_type

    # See comments for similar actual to formal code above
    actual_args = []
    actual_arg_kinds = []
    actual_arg_names = []
    seen_args = set()
    for i, param in enumerate(ctx.args):
        for j, a in enumerate(param):
            if a in seen_args:
                continue
            seen_args.add(a)
            actual_args.append(a)
            actual_arg_kinds.append(ctx.arg_kinds[i][j])
            actual_arg_names.append(ctx.arg_names[i][j])

    result = ctx.api.expr_checker.check_call(
        callee=partial_type,
        args=actual_args,
        arg_kinds=actual_arg_kinds,
        arg_names=actual_arg_names,
        context=ctx.context,
    )
    return result[0]
