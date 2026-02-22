"""Type checking of attribute access"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Callable, TypeVar, cast

from mypy import message_registry, state, subtypes
from mypy.checker_shared import TypeCheckerSharedApi
from mypy.erasetype import erase_typevars
from mypy.expandtype import (
    expand_self_type,
    expand_type_by_instance,
    freshen_all_functions_type_vars,
)
from mypy.maptype import map_instance_to_supertype
from mypy.messages import MessageBuilder
from mypy.nodes import (
    ARG_POS,
    ARG_STAR,
    ARG_STAR2,
    EXCLUDED_ENUM_ATTRIBUTES,
    SYMBOL_FUNCBASE_TYPES,
    Context,
    Decorator,
    Expression,
    FuncBase,
    FuncDef,
    IndexExpr,
    MypyFile,
    NameExpr,
    OverloadedFuncDef,
    SymbolTable,
    TempNode,
    TypeAlias,
    TypeInfo,
    TypeVarLikeExpr,
    Var,
    is_final_node,
)
from mypy.plugin import AttributeContext
from mypy.typeops import (
    bind_self,
    erase_to_bound,
    freeze_all_type_vars,
    function_type,
    get_all_type_vars,
    get_type_vars,
    make_simplified_union,
    supported_self_type,
    tuple_fallback,
)
from mypy.types import (
    AnyType,
    CallableType,
    DeletedType,
    FunctionLike,
    Instance,
    LiteralType,
    NoneType,
    Overloaded,
    ParamSpecType,
    PartialType,
    ProperType,
    TupleType,
    Type,
    TypedDictType,
    TypeOfAny,
    TypeType,
    TypeVarLikeType,
    TypeVarTupleType,
    TypeVarType,
    UnionType,
    get_proper_type,
)


class MemberContext:
    """Information and objects needed to type check attribute access.

    Look at the docstring of analyze_member_access for more information.
    """

    def __init__(
        self,
        *,
        is_lvalue: bool,
        is_super: bool,
        is_operator: bool,
        original_type: Type,
        context: Context,
        chk: TypeCheckerSharedApi,
        self_type: Type | None = None,
        module_symbol_table: SymbolTable | None = None,
        no_deferral: bool = False,
        is_self: bool = False,
        rvalue: Expression | None = None,
        suppress_errors: bool = False,
        preserve_type_var_ids: bool = False,
    ) -> None:
        self.is_lvalue = is_lvalue
        self.is_super = is_super
        self.is_operator = is_operator
        self.original_type = original_type
        self.self_type = self_type or original_type
        self.context = context  # Error context
        self.chk = chk
        self.msg = chk.msg
        self.module_symbol_table = module_symbol_table
        self.no_deferral = no_deferral
        self.is_self = is_self
        if rvalue is not None:
            assert is_lvalue
        self.rvalue = rvalue
        self.suppress_errors = suppress_errors
        # This attribute is only used to preserve old protocol member access logic.
        # It is needed to avoid infinite recursion in cases involving self-referential
        # generic methods, see find_member() for details. Do not use for other purposes!
        self.preserve_type_var_ids = preserve_type_var_ids

    def named_type(self, name: str) -> Instance:
        return self.chk.named_type(name)

    def not_ready_callback(self, name: str, context: Context) -> None:
        self.chk.handle_cannot_determine_type(name, context)

    def fail(self, msg: str) -> None:
        if not self.suppress_errors:
            self.msg.fail(msg, self.context)

    def copy_modified(
        self,
        *,
        self_type: Type | None = None,
        is_lvalue: bool | None = None,
        original_type: Type | None = None,
    ) -> MemberContext:
        mx = MemberContext(
            is_lvalue=self.is_lvalue,
            is_super=self.is_super,
            is_operator=self.is_operator,
            original_type=self.original_type,
            context=self.context,
            chk=self.chk,
            self_type=self.self_type,
            module_symbol_table=self.module_symbol_table,
            no_deferral=self.no_deferral,
            rvalue=self.rvalue,
            suppress_errors=self.suppress_errors,
            preserve_type_var_ids=self.preserve_type_var_ids,
        )
        if self_type is not None:
            mx.self_type = self_type
        if is_lvalue is not None:
            mx.is_lvalue = is_lvalue
        if original_type is not None:
            mx.original_type = original_type
        return mx


def analyze_member_access(
    name: str,
    typ: Type,
    context: Context,
    *,
    is_lvalue: bool,
    is_super: bool,
    is_operator: bool,
    original_type: Type,
    chk: TypeCheckerSharedApi,
    override_info: TypeInfo | None = None,
    in_literal_context: bool = False,
    self_type: Type | None = None,
    module_symbol_table: SymbolTable | None = None,
    no_deferral: bool = False,
    is_self: bool = False,
    rvalue: Expression | None = None,
    suppress_errors: bool = False,
) -> Type:
    """Return the type of attribute 'name' of 'typ'.

    The actual implementation is in '_analyze_member_access' and this docstring
    also applies to it.

    This is a general operation that supports various different variations:

      1. lvalue or non-lvalue access (setter or getter access)
      2. supertype access when using super() (is_super == True and
         'override_info' should refer to the supertype)

    'original_type' is the most precise inferred or declared type of the base object
    that we have available. When looking for an attribute of 'typ', we may perform
    recursive calls targeting the fallback type, and 'typ' may become some supertype
    of 'original_type'. 'original_type' is always preserved as the 'typ' type used in
    the initial, non-recursive call. The 'self_type' is a component of 'original_type'
    to which generic self should be bound (a narrower type that has a fallback to instance).
    Currently, this is used only for union types.

    'module_symbol_table' is passed to this function if 'typ' is actually a module,
    and we want to keep track of the available attributes of the module (since they
    are not available via the type object directly)

    'rvalue' can be provided optionally to infer better setter type when is_lvalue is True,
    most notably this helps for descriptors with overloaded __set__() method.

    'suppress_errors' will skip any logic that is only needed to generate error messages.
    Note that this more of a performance optimization, one should not rely on this to not
    show any messages, as some may be show e.g. by callbacks called here,
    use msg.filter_errors(), if needed.
    """
    mx = MemberContext(
        is_lvalue=is_lvalue,
        is_super=is_super,
        is_operator=is_operator,
        original_type=original_type,
        context=context,
        chk=chk,
        self_type=self_type,
        module_symbol_table=module_symbol_table,
        no_deferral=no_deferral,
        is_self=is_self,
        rvalue=rvalue,
        suppress_errors=suppress_errors,
    )
    result = _analyze_member_access(name, typ, mx, override_info)
    possible_literal = get_proper_type(result)
    if (
        in_literal_context
        and isinstance(possible_literal, Instance)
        and possible_literal.last_known_value is not None
    ):
        return possible_literal.last_known_value
    else:
        return result


def _analyze_member_access(
    name: str, typ: Type, mx: MemberContext, override_info: TypeInfo | None = None
) -> Type:
    typ = get_proper_type(typ)
    if isinstance(typ, Instance):
        return analyze_instance_member_access(name, typ, mx, override_info)
    elif isinstance(typ, AnyType):
        # The base object has dynamic type.
        return AnyType(TypeOfAny.from_another_any, source_any=typ)
    elif isinstance(typ, UnionType):
        return analyze_union_member_access(name, typ, mx)
    elif isinstance(typ, FunctionLike) and typ.is_type_obj():
        return analyze_type_callable_member_access(name, typ, mx)
    elif isinstance(typ, TypeType):
        return analyze_type_type_member_access(name, typ, mx, override_info)
    elif isinstance(typ, TupleType):
        # Actually look up from the fallback instance type.
        return _analyze_member_access(name, tuple_fallback(typ), mx, override_info)
    elif isinstance(typ, (LiteralType, FunctionLike)):
        # Actually look up from the fallback instance type.
        return _analyze_member_access(name, typ.fallback, mx, override_info)
    elif isinstance(typ, TypedDictType):
        return analyze_typeddict_access(name, typ, mx, override_info)
    elif isinstance(typ, NoneType):
        return analyze_none_member_access(name, typ, mx)
    elif isinstance(typ, TypeVarLikeType):
        if isinstance(typ, TypeVarType) and typ.values:
            return _analyze_member_access(
                name, make_simplified_union(typ.values), mx, override_info
            )
        return _analyze_member_access(name, typ.upper_bound, mx, override_info)
    elif isinstance(typ, DeletedType):
        if not mx.suppress_errors:
            mx.msg.deleted_as_rvalue(typ, mx.context)
        return AnyType(TypeOfAny.from_error)
    return report_missing_attribute(mx.original_type, typ, name, mx)


def may_be_awaitable_attribute(
    name: str, typ: Type, mx: MemberContext, override_info: TypeInfo | None = None
) -> bool:
    """Check if the given type has the attribute when awaited."""
    if mx.chk.checking_missing_await:
        # Avoid infinite recursion.
        return False
    with mx.chk.checking_await_set(), mx.msg.filter_errors() as local_errors:
        aw_type = mx.chk.get_precise_awaitable_type(typ, local_errors)
        if aw_type is None:
            return False
        _ = _analyze_member_access(
            name, aw_type, mx.copy_modified(self_type=aw_type), override_info
        )
        return not local_errors.has_new_errors()


def report_missing_attribute(
    original_type: Type,
    typ: Type,
    name: str,
    mx: MemberContext,
    override_info: TypeInfo | None = None,
) -> Type:
    if mx.suppress_errors:
        return AnyType(TypeOfAny.from_error)
    error_code = mx.msg.has_no_attr(original_type, typ, name, mx.context, mx.module_symbol_table)
    if not mx.msg.prefer_simple_messages():
        if may_be_awaitable_attribute(name, typ, mx, override_info):
            mx.msg.possible_missing_await(mx.context, error_code)
    return AnyType(TypeOfAny.from_error)


# The several functions that follow implement analyze_member_access for various
# types and aren't documented individually.


def analyze_instance_member_access(
    name: str, typ: Instance, mx: MemberContext, override_info: TypeInfo | None
) -> Type:
    info = typ.type
    if override_info:
        info = override_info

    method = info.get_method(name)

    if name == "__init__" and not mx.is_super and not info.is_final:
        if not method or not method.is_final:
            # Accessing __init__ in statically typed code would compromise
            # type safety unless used via super() or the method/class is final.
            mx.fail(message_registry.CANNOT_ACCESS_INIT)
            return AnyType(TypeOfAny.from_error)

    # The base object has an instance type.

    if (
        state.find_occurrences
        and info.name == state.find_occurrences[0]
        and name == state.find_occurrences[1]
        and not mx.suppress_errors
    ):
        mx.msg.note("Occurrence of '{}.{}'".format(*state.find_occurrences), mx.context)

    # Look up the member. First look up the method dictionary.
    if method and not isinstance(method, Decorator):
        if mx.is_super and not mx.suppress_errors:
            validate_super_call(method, mx)

        if method.is_property:
            assert isinstance(method, OverloadedFuncDef)
            getter = method.items[0]
            assert isinstance(getter, Decorator)
            if mx.is_lvalue and getter.var.is_settable_property:
                mx.chk.warn_deprecated(method.setter, mx.context)
            return analyze_var(name, getter.var, typ, mx)

        if mx.is_lvalue and not mx.suppress_errors:
            mx.msg.cant_assign_to_method(mx.context)
        if not isinstance(method, OverloadedFuncDef):
            signature = function_type(method, mx.named_type("builtins.function"))
        else:
            if method.type is None:
                # Overloads may be not ready if they are decorated. Handle this in same
                # manner as we would handle a regular decorated function: defer if possible.
                if not mx.no_deferral and method.items:
                    mx.not_ready_callback(method.name, mx.context)
                return AnyType(TypeOfAny.special_form)
            assert isinstance(method.type, Overloaded)
            signature = method.type
        if not mx.preserve_type_var_ids:
            signature = freshen_all_functions_type_vars(signature)
        if not method.is_static:
            if isinstance(method, (FuncDef, OverloadedFuncDef)) and method.is_trivial_self:
                signature = bind_self_fast(signature, mx.self_type)
            else:
                signature = check_self_arg(
                    signature, mx.self_type, method.is_class, mx.context, name, mx.msg
                )
                signature = bind_self(signature, mx.self_type, is_classmethod=method.is_class)
        # TODO: should we skip these steps for static methods as well?
        # Since generic static methods should not be allowed.
        typ = map_instance_to_supertype(typ, method.info)
        member_type = expand_type_by_instance(signature, typ)
        freeze_all_type_vars(member_type)
        return member_type
    else:
        # Not a method.
        return analyze_member_var_access(name, typ, info, mx)


def validate_super_call(node: FuncBase, mx: MemberContext) -> None:
    unsafe_super = False
    if isinstance(node, FuncDef) and node.is_trivial_body:
        unsafe_super = True
    elif isinstance(node, OverloadedFuncDef):
        if node.impl:
            impl = node.impl if isinstance(node.impl, FuncDef) else node.impl.func
            unsafe_super = impl.is_trivial_body
        elif not node.is_property and node.items:
            assert isinstance(node.items[0], Decorator)
            unsafe_super = node.items[0].func.is_trivial_body
    if unsafe_super:
        mx.msg.unsafe_super(node.name, node.info.name, mx.context)


def analyze_type_callable_member_access(name: str, typ: FunctionLike, mx: MemberContext) -> Type:
    # Class attribute.
    # TODO super?
    ret_type = typ.items[0].ret_type
    assert isinstance(ret_type, ProperType)
    if isinstance(ret_type, TupleType):
        ret_type = tuple_fallback(ret_type)
    if isinstance(ret_type, TypedDictType):
        ret_type = ret_type.fallback
    if isinstance(ret_type, Instance):
        if not mx.is_operator:
            # When Python sees an operator (eg `3 == 4`), it automatically translates that
            # into something like `int.__eq__(3, 4)` instead of `(3).__eq__(4)` as an
            # optimization.
            #
            # While it normally it doesn't matter which of the two versions are used, it
            # does cause inconsistencies when working with classes. For example, translating
            # `int == int` to `int.__eq__(int)` would not work since `int.__eq__` is meant to
            # compare two int _instances_. What we really want is `type(int).__eq__`, which
            # is meant to compare two types or classes.
            #
            # This check makes sure that when we encounter an operator, we skip looking up
            # the corresponding method in the current instance to avoid this edge case.
            # See https://github.com/python/mypy/pull/1787 for more info.
            # TODO: do not rely on same type variables being present in all constructor overloads.
            result = analyze_class_attribute_access(
                ret_type, name, mx, original_vars=typ.items[0].variables, mcs_fallback=typ.fallback
            )
            if result:
                return result
        # Look up from the 'type' type.
        return _analyze_member_access(name, typ.fallback, mx)
    else:
        assert False, f"Unexpected type {ret_type!r}"


def analyze_type_type_member_access(
    name: str, typ: TypeType, mx: MemberContext, override_info: TypeInfo | None
) -> Type:
    # Similar to analyze_type_callable_attribute_access.
    item = None
    fallback = mx.named_type("builtins.type")
    if isinstance(typ.item, Instance):
        item = typ.item
    elif isinstance(typ.item, AnyType):
        with mx.msg.filter_errors():
            return _analyze_member_access(name, fallback, mx, override_info)
    elif isinstance(typ.item, TypeVarType):
        upper_bound = get_proper_type(typ.item.upper_bound)
        if isinstance(upper_bound, Instance):
            item = upper_bound
        elif isinstance(upper_bound, UnionType):
            return _analyze_member_access(
                name,
                TypeType.make_normalized(upper_bound, line=typ.line, column=typ.column),
                mx,
                override_info,
            )
        elif isinstance(upper_bound, TupleType):
            item = tuple_fallback(upper_bound)
        elif isinstance(upper_bound, AnyType):
            with mx.msg.filter_errors():
                return _analyze_member_access(name, fallback, mx, override_info)
    elif isinstance(typ.item, TupleType):
        item = tuple_fallback(typ.item)
    elif isinstance(typ.item, FunctionLike) and typ.item.is_type_obj():
        item = typ.item.fallback
    elif isinstance(typ.item, TypeType):
        # Access member on metaclass object via Type[Type[C]]
        if isinstance(typ.item.item, Instance):
            item = typ.item.item.type.metaclass_type
    ignore_messages = False

    if item is not None:
        fallback = item.type.metaclass_type or fallback

    if item and not mx.is_operator:
        # See comment above for why operators are skipped
        result = analyze_class_attribute_access(
            item, name, mx, mcs_fallback=fallback, override_info=override_info
        )
        if result:
            if not (isinstance(get_proper_type(result), AnyType) and item.type.fallback_to_any):
                return result
            else:
                # We don't want errors on metaclass lookup for classes with Any fallback
                ignore_messages = True

    with mx.msg.filter_errors(filter_errors=ignore_messages):
        return _analyze_member_access(name, fallback, mx, override_info)


def analyze_union_member_access(name: str, typ: UnionType, mx: MemberContext) -> Type:
    with mx.msg.disable_type_names():
        results = []
        for subtype in typ.relevant_items():
            # Self types should be bound to every individual item of a union.
            item_mx = mx.copy_modified(self_type=subtype)
            results.append(_analyze_member_access(name, subtype, item_mx))
    return make_simplified_union(results)


def analyze_none_member_access(name: str, typ: NoneType, mx: MemberContext) -> Type:
    if name == "__bool__":
        literal_false = LiteralType(False, fallback=mx.named_type("builtins.bool"))
        return CallableType(
            arg_types=[],
            arg_kinds=[],
            arg_names=[],
            ret_type=literal_false,
            fallback=mx.named_type("builtins.function"),
        )
    else:
        return _analyze_member_access(name, mx.named_type("builtins.object"), mx)


def analyze_member_var_access(
    name: str, itype: Instance, info: TypeInfo, mx: MemberContext
) -> Type:
    """Analyse attribute access that does not target a method.

    This is logically part of analyze_member_access and the arguments are similar.

    original_type is the type of E in the expression E.var
    """
    # It was not a method. Try looking up a variable.
    node = info.get(name)
    v = node.node if node else None

    mx.chk.warn_deprecated(v, mx.context)

    vv = v
    is_trivial_self = False
    if isinstance(vv, Decorator):
        # The associated Var node of a decorator contains the type.
        v = vv.var
        is_trivial_self = vv.func.is_trivial_self and not vv.decorators
        if mx.is_super and not mx.suppress_errors:
            validate_super_call(vv.func, mx)
    if isinstance(v, FuncDef):
        assert False, "Did not expect a function"
    if isinstance(v, MypyFile):
        mx.chk.module_refs.add(v.fullname)

    if isinstance(vv, (TypeInfo, TypeAlias, MypyFile, TypeVarLikeExpr)):
        # If the associated variable is a TypeInfo synthesize a Var node for
        # the purposes of type checking.  This enables us to type check things
        # like accessing class attributes on an inner class. Similar we allow
        # using qualified type aliases in runtime context. For example:
        #     class C:
        #         A = List[int]
        #     x = C.A() <- this is OK
        typ = mx.chk.expr_checker.analyze_static_reference(vv, mx.context, mx.is_lvalue)
        v = Var(name, type=typ)
        v.info = info

    if isinstance(v, Var):
        implicit = info[name].implicit

        # An assignment to final attribute is always an error,
        # independently of types.
        if mx.is_lvalue and not mx.chk.get_final_context():
            check_final_member(name, info, mx.msg, mx.context)

        return analyze_var(name, v, itype, mx, implicit=implicit, is_trivial_self=is_trivial_self)
    elif (
        not v
        and name not in ["__getattr__", "__setattr__", "__getattribute__"]
        and not mx.is_operator
        and mx.module_symbol_table is None
    ):
        # Above we skip ModuleType.__getattr__ etc. if we have a
        # module symbol table, since the symbol table allows precise
        # checking.
        if not mx.is_lvalue:
            for method_name in ("__getattribute__", "__getattr__"):
                method = info.get_method(method_name)

                # __getattribute__ is defined on builtins.object and returns Any, so without
                # the guard this search will always find object.__getattribute__ and conclude
                # that the attribute exists
                if method and method.info.fullname != "builtins.object":
                    bound_method = analyze_decorator_or_funcbase_access(
                        defn=method, itype=itype, name=method_name, mx=mx
                    )
                    typ = map_instance_to_supertype(itype, method.info)
                    getattr_type = get_proper_type(expand_type_by_instance(bound_method, typ))
                    if isinstance(getattr_type, CallableType):
                        result = getattr_type.ret_type
                    else:
                        result = getattr_type

                    # Call the attribute hook before returning.
                    fullname = f"{method.info.fullname}.{name}"
                    hook = mx.chk.plugin.get_attribute_hook(fullname)
                    if hook:
                        result = hook(
                            AttributeContext(
                                get_proper_type(mx.original_type),
                                result,
                                mx.is_lvalue,
                                mx.context,
                                mx.chk,
                            )
                        )
                    return result
        else:
            setattr_meth = info.get_method("__setattr__")
            if setattr_meth and setattr_meth.info.fullname != "builtins.object":
                bound_type = analyze_decorator_or_funcbase_access(
                    defn=setattr_meth,
                    itype=itype,
                    name="__setattr__",
                    mx=mx.copy_modified(is_lvalue=False),
                )
                typ = map_instance_to_supertype(itype, setattr_meth.info)
                setattr_type = get_proper_type(expand_type_by_instance(bound_type, typ))
                if isinstance(setattr_type, CallableType) and len(setattr_type.arg_types) > 0:
                    return setattr_type.arg_types[-1]

    if itype.type.fallback_to_any:
        return AnyType(TypeOfAny.special_form)

    # Could not find the member.
    if itype.extra_attrs and name in itype.extra_attrs.attrs:
        # For modules use direct symbol table lookup.
        if not itype.extra_attrs.mod_name:
            return itype.extra_attrs.attrs[name]

    if mx.is_super and not mx.suppress_errors:
        mx.msg.undefined_in_superclass(name, mx.context)
        return AnyType(TypeOfAny.from_error)
    else:
        ret = report_missing_attribute(mx.original_type, itype, name, mx)
        # Avoid paying double jeopardy if we can't find the member due to --no-implicit-reexport
        if (
            mx.module_symbol_table is not None
            and name in mx.module_symbol_table
            and not mx.module_symbol_table[name].module_public
        ):
            v = mx.module_symbol_table[name].node
            e = NameExpr(name)
            e.set_line(mx.context)
            e.node = v
            return mx.chk.expr_checker.analyze_ref_expr(e, lvalue=mx.is_lvalue)
        return ret


def check_final_member(name: str, info: TypeInfo, msg: MessageBuilder, ctx: Context) -> None:
    """Give an error if the name being assigned was declared as final."""
    for base in info.mro:
        sym = base.names.get(name)
        if sym and is_final_node(sym.node):
            msg.cant_assign_to_final(name, attr_assign=True, ctx=ctx)


def analyze_descriptor_access(descriptor_type: Type, mx: MemberContext) -> Type:
    """Type check descriptor access.

    Arguments:
        descriptor_type: The type of the descriptor attribute being accessed
            (the type of ``f`` in ``a.f`` when ``f`` is a descriptor).
        mx: The current member access context.
    Return:
        The return type of the appropriate ``__get__/__set__`` overload for the descriptor.
    """
    instance_type = get_proper_type(mx.self_type)
    orig_descriptor_type = descriptor_type
    descriptor_type = get_proper_type(descriptor_type)

    if isinstance(descriptor_type, UnionType):
        # Map the access over union types
        return make_simplified_union(
            [analyze_descriptor_access(typ, mx) for typ in descriptor_type.items]
        )
    elif not isinstance(descriptor_type, Instance):
        return orig_descriptor_type

    if not mx.is_lvalue and not descriptor_type.type.has_readable_member("__get__"):
        return orig_descriptor_type

    # We do this check first to accommodate for descriptors with only __set__ method.
    # If there is no __set__, we type-check that the assigned value matches
    # the return type of __get__. This doesn't match the python semantics,
    # (which allow you to override the descriptor with any value), but preserves
    # the type of accessing the attribute (even after the override).
    if mx.is_lvalue and descriptor_type.type.has_readable_member("__set__"):
        return analyze_descriptor_assign(descriptor_type, mx)

    if mx.is_lvalue and not descriptor_type.type.has_readable_member("__get__"):
        # This turned out to be not a descriptor after all.
        return orig_descriptor_type

    dunder_get = descriptor_type.type.get_method("__get__")
    if dunder_get is None:
        mx.fail(
            message_registry.DESCRIPTOR_GET_NOT_CALLABLE.format(
                descriptor_type.str_with_options(mx.msg.options)
            )
        )
        return AnyType(TypeOfAny.from_error)

    bound_method = analyze_decorator_or_funcbase_access(
        defn=dunder_get,
        itype=descriptor_type,
        name="__get__",
        mx=mx.copy_modified(self_type=descriptor_type),
    )

    typ = map_instance_to_supertype(descriptor_type, dunder_get.info)
    dunder_get_type = expand_type_by_instance(bound_method, typ)

    if isinstance(instance_type, FunctionLike) and instance_type.is_type_obj():
        owner_type = instance_type.items[0].ret_type
        instance_type = NoneType()
    elif isinstance(instance_type, TypeType):
        owner_type = instance_type.item
        instance_type = NoneType()
    else:
        owner_type = instance_type

    callable_name = mx.chk.expr_checker.method_fullname(descriptor_type, "__get__")
    dunder_get_type = mx.chk.expr_checker.transform_callee_type(
        callable_name,
        dunder_get_type,
        [
            TempNode(instance_type, context=mx.context),
            TempNode(TypeType.make_normalized(owner_type), context=mx.context),
        ],
        [ARG_POS, ARG_POS],
        mx.context,
        object_type=descriptor_type,
    )

    _, inferred_dunder_get_type = mx.chk.expr_checker.check_call(
        dunder_get_type,
        [
            TempNode(instance_type, context=mx.context),
            TempNode(TypeType.make_normalized(owner_type), context=mx.context),
        ],
        [ARG_POS, ARG_POS],
        mx.context,
        object_type=descriptor_type,
        callable_name=callable_name,
    )

    mx.chk.check_deprecated(dunder_get, mx.context)
    mx.chk.warn_deprecated_overload_item(
        dunder_get, mx.context, target=inferred_dunder_get_type, selftype=descriptor_type
    )

    inferred_dunder_get_type = get_proper_type(inferred_dunder_get_type)
    if isinstance(inferred_dunder_get_type, AnyType):
        # check_call failed, and will have reported an error
        return inferred_dunder_get_type

    if not isinstance(inferred_dunder_get_type, CallableType):
        mx.fail(
            message_registry.DESCRIPTOR_GET_NOT_CALLABLE.format(
                descriptor_type.str_with_options(mx.msg.options)
            )
        )
        return AnyType(TypeOfAny.from_error)

    return inferred_dunder_get_type.ret_type


def analyze_descriptor_assign(descriptor_type: Instance, mx: MemberContext) -> Type:
    instance_type = get_proper_type(mx.self_type)
    dunder_set = descriptor_type.type.get_method("__set__")
    if dunder_set is None:
        mx.fail(
            message_registry.DESCRIPTOR_SET_NOT_CALLABLE.format(
                descriptor_type.str_with_options(mx.msg.options)
            ).value
        )
        return AnyType(TypeOfAny.from_error)

    bound_method = analyze_decorator_or_funcbase_access(
        defn=dunder_set,
        itype=descriptor_type,
        name="__set__",
        mx=mx.copy_modified(is_lvalue=False, self_type=descriptor_type),
    )
    typ = map_instance_to_supertype(descriptor_type, dunder_set.info)
    dunder_set_type = expand_type_by_instance(bound_method, typ)

    callable_name = mx.chk.expr_checker.method_fullname(descriptor_type, "__set__")
    rvalue = mx.rvalue or TempNode(AnyType(TypeOfAny.special_form), context=mx.context)
    dunder_set_type = mx.chk.expr_checker.transform_callee_type(
        callable_name,
        dunder_set_type,
        [TempNode(instance_type, context=mx.context), rvalue],
        [ARG_POS, ARG_POS],
        mx.context,
        object_type=descriptor_type,
    )

    # For non-overloaded setters, the result should be type-checked like a regular assignment.
    # Hence, we first only try to infer the type by using the rvalue as type context.
    type_context = rvalue
    with mx.msg.filter_errors():
        _, inferred_dunder_set_type = mx.chk.expr_checker.check_call(
            dunder_set_type,
            [TempNode(instance_type, context=mx.context), type_context],
            [ARG_POS, ARG_POS],
            mx.context,
            object_type=descriptor_type,
            callable_name=callable_name,
        )

    # And now we in fact type check the call, to show errors related to wrong arguments
    # count, etc., replacing the type context for non-overloaded setters only.
    inferred_dunder_set_type = get_proper_type(inferred_dunder_set_type)
    if isinstance(inferred_dunder_set_type, CallableType):
        type_context = TempNode(AnyType(TypeOfAny.special_form), context=mx.context)
    mx.chk.expr_checker.check_call(
        dunder_set_type,
        [TempNode(instance_type, context=mx.context), type_context],
        [ARG_POS, ARG_POS],
        mx.context,
        object_type=descriptor_type,
        callable_name=callable_name,
    )

    # Search for possible deprecations:
    mx.chk.check_deprecated(dunder_set, mx.context)
    mx.chk.warn_deprecated_overload_item(
        dunder_set, mx.context, target=inferred_dunder_set_type, selftype=descriptor_type
    )

    # In the following cases, a message already will have been recorded in check_call.
    if (not isinstance(inferred_dunder_set_type, CallableType)) or (
        len(inferred_dunder_set_type.arg_types) < 2
    ):
        return AnyType(TypeOfAny.from_error)
    return inferred_dunder_set_type.arg_types[1]


def is_instance_var(var: Var) -> bool:
    """Return if var is an instance variable according to PEP 526."""
    return (
        # check the type_info node is the var (not a decorated function, etc.)
        var.name in var.info.names
        and var.info.names[var.name].node is var
        and not var.is_classvar
        # variables without annotations are treated as classvar
        and not var.is_inferred
    )


def analyze_var(
    name: str,
    var: Var,
    itype: Instance,
    mx: MemberContext,
    *,
    implicit: bool = False,
    is_trivial_self: bool = False,
) -> Type:
    """Analyze access to an attribute via a Var node.

    This is conceptually part of analyze_member_access and the arguments are similar.
    itype is the instance type in which attribute should be looked up
    original_type is the type of E in the expression E.var
    if implicit is True, the original Var was created as an assignment to self
    if is_trivial_self is True, we can use fast path for bind_self().
    """
    # Found a member variable.
    original_itype = itype
    itype = map_instance_to_supertype(itype, var.info)
    if var.is_settable_property and mx.is_lvalue:
        typ: Type | None = var.setter_type
        if typ is None and var.is_ready:
            # Existing synthetic properties may not set setter type. Fall back to getter.
            typ = var.type
    else:
        typ = var.type
    if typ:
        if isinstance(typ, PartialType):
            return mx.chk.handle_partial_var_type(typ, mx.is_lvalue, var, mx.context)
        if mx.is_lvalue and not mx.suppress_errors:
            if var.is_property and not var.is_settable_property:
                mx.msg.read_only_property(name, itype.type, mx.context)
            if var.is_classvar:
                mx.msg.cant_assign_to_classvar(name, mx.context)
        # This is the most common case for variables, so start with this.
        result = expand_without_binding(typ, var, itype, original_itype, mx)

        # A non-None value indicates that we should actually bind self for this variable.
        call_type: ProperType | None = None
        if var.is_initialized_in_class and (not is_instance_var(var) or mx.is_operator):
            typ = get_proper_type(typ)
            if isinstance(typ, FunctionLike) and not typ.is_type_obj():
                call_type = typ
            elif var.is_property:
                deco_mx = mx.copy_modified(original_type=typ, self_type=typ, is_lvalue=False)
                call_type = get_proper_type(_analyze_member_access("__call__", typ, deco_mx))
            else:
                call_type = typ

        # Bound variables with callable types are treated like methods
        # (these are usually method aliases like __rmul__ = __mul__).
        if isinstance(call_type, FunctionLike) and not call_type.is_type_obj():
            if mx.is_lvalue and not var.is_property and not mx.suppress_errors:
                mx.msg.cant_assign_to_method(mx.context)

        # Bind the self type for each callable component (when needed).
        if call_type and not var.is_staticmethod:
            bound_items = []
            for ct in call_type.items if isinstance(call_type, UnionType) else [call_type]:
                p_ct = get_proper_type(ct)
                if isinstance(p_ct, FunctionLike) and (not p_ct.bound() or var.is_property):
                    item = expand_and_bind_callable(p_ct, var, itype, name, mx, is_trivial_self)
                else:
                    item = expand_without_binding(ct, var, itype, original_itype, mx)
                bound_items.append(item)
            result = UnionType.make_union(bound_items)
    else:
        if not var.is_ready and not mx.no_deferral:
            mx.not_ready_callback(var.name, mx.context)
        # Implicit 'Any' type.
        result = AnyType(TypeOfAny.special_form)
    fullname = f"{var.info.fullname}.{name}"
    hook = mx.chk.plugin.get_attribute_hook(fullname)
    if result and not (implicit or var.info.is_protocol and is_instance_var(var)):
        result = analyze_descriptor_access(result, mx)
    if hook:
        result = hook(
            AttributeContext(
                get_proper_type(mx.original_type), result, mx.is_lvalue, mx.context, mx.chk
            )
        )
    return result


def expand_without_binding(
    typ: Type, var: Var, itype: Instance, original_itype: Instance, mx: MemberContext
) -> Type:
    if not mx.preserve_type_var_ids:
        typ = freshen_all_functions_type_vars(typ)
    typ = expand_self_type_if_needed(typ, mx, var, original_itype)
    expanded = expand_type_by_instance(typ, itype)
    freeze_all_type_vars(expanded)
    return expanded


def expand_and_bind_callable(
    functype: FunctionLike,
    var: Var,
    itype: Instance,
    name: str,
    mx: MemberContext,
    is_trivial_self: bool,
) -> Type:
    if not mx.preserve_type_var_ids:
        functype = freshen_all_functions_type_vars(functype)
    typ = get_proper_type(expand_self_type(var, functype, mx.original_type))
    assert isinstance(typ, FunctionLike)
    if is_trivial_self:
        typ = bind_self_fast(typ, mx.self_type)
    else:
        typ = check_self_arg(typ, mx.self_type, var.is_classmethod, mx.context, name, mx.msg)
        typ = bind_self(typ, mx.self_type, var.is_classmethod)
    expanded = expand_type_by_instance(typ, itype)
    freeze_all_type_vars(expanded)
    if not var.is_property:
        return expanded
    # TODO: a decorated property can result in Overloaded here.
    assert isinstance(expanded, CallableType)
    if var.is_settable_property and mx.is_lvalue and var.setter_type is not None:
        # TODO: use check_call() to infer better type, same as for __set__().
        if not expanded.arg_types:
            # This can happen when accessing invalid property from its own body,
            # error will be reported elsewhere.
            return AnyType(TypeOfAny.from_error)
        return expanded.arg_types[0]
    else:
        return expanded.ret_type


def expand_self_type_if_needed(
    t: Type, mx: MemberContext, var: Var, itype: Instance, is_class: bool = False
) -> Type:
    """Expand special Self type in a backwards compatible manner.

    This should ensure that mixing old-style and new-style self-types work
    seamlessly. Also, re-bind new style self-types in subclasses if needed.
    """
    original = get_proper_type(mx.self_type)
    if not (mx.is_self or mx.is_super):
        repl = mx.self_type
        if is_class:
            if isinstance(original, TypeType):
                repl = original.item
            elif isinstance(original, CallableType):
                # Problematic access errors should have been already reported.
                repl = erase_typevars(original.ret_type)
            else:
                repl = itype
        return expand_self_type(var, t, repl)
    elif supported_self_type(
        # Support compatibility with plain old style T -> T and Type[T] -> T only.
        get_proper_type(mx.self_type),
        allow_instances=False,
        allow_callable=False,
    ):
        repl = mx.self_type
        if is_class and isinstance(original, TypeType):
            repl = original.item
        return expand_self_type(var, t, repl)
    elif (
        mx.is_self
        and itype.type != var.info
        # If an attribute with Self-type was defined in a supertype, we need to
        # rebind the Self type variable to Self type variable of current class...
        and itype.type.self_type is not None
        # ...unless `self` has an explicit non-trivial annotation.
        and itype == mx.chk.scope.active_self_type()
    ):
        return expand_self_type(var, t, itype.type.self_type)
    else:
        return t


def check_self_arg(
    functype: FunctionLike,
    dispatched_arg_type: Type,
    is_classmethod: bool,
    context: Context,
    name: str,
    msg: MessageBuilder,
) -> FunctionLike:
    """Check that an instance has a valid type for a method with annotated 'self'.

    For example if the method is defined as:
        class A:
            def f(self: S) -> T: ...
    then for 'x.f' we check that type(x) <: S. If the method is overloaded, we select
    only overloads items that satisfy this requirement. If there are no matching
    overloads, an error is generated.
    """
    items = functype.items
    if not items:
        return functype
    new_items = []
    if is_classmethod:
        dispatched_arg_type = TypeType.make_normalized(dispatched_arg_type)

    for item in items:
        if not item.arg_types or item.arg_kinds[0] not in (ARG_POS, ARG_STAR):
            # No positional first (self) argument (*args is okay).
            msg.no_formal_self(name, item, context)
            # This is pretty bad, so just return the original signature if
            # there is at least one such error.
            return functype
        else:
            selfarg = get_proper_type(item.arg_types[0])
            # This matches similar special-casing in bind_self(), see more details there.
            self_callable = name == "__call__" and isinstance(selfarg, CallableType)
            if self_callable or subtypes.is_subtype(
                dispatched_arg_type,
                # This level of erasure matches the one in checker.check_func_def(),
                # better keep these two checks consistent.
                erase_typevars(erase_to_bound(selfarg)),
                # This is to work around the fact that erased ParamSpec and TypeVarTuple
                # callables are not always compatible with non-erased ones both ways.
                always_covariant=any(
                    not isinstance(tv, TypeVarType) for tv in get_all_type_vars(selfarg)
                ),
                ignore_pos_arg_names=True,
            ):
                new_items.append(item)
            elif isinstance(selfarg, ParamSpecType):
                # TODO: This is not always right. What's the most reasonable thing to do here?
                new_items.append(item)
            elif isinstance(selfarg, TypeVarTupleType):
                raise NotImplementedError
    if not new_items:
        # Choose first item for the message (it may be not very helpful for overloads).
        msg.incompatible_self_argument(
            name, dispatched_arg_type, items[0], is_classmethod, context
        )
        return functype
    if len(new_items) == 1:
        return new_items[0]
    return Overloaded(new_items)


def analyze_class_attribute_access(
    itype: Instance,
    name: str,
    mx: MemberContext,
    *,
    mcs_fallback: Instance,
    override_info: TypeInfo | None = None,
    original_vars: Sequence[TypeVarLikeType] | None = None,
) -> Type | None:
    """Analyze access to an attribute on a class object.

    itype is the return type of the class object callable, original_type is the type
    of E in the expression E.var, original_vars are type variables of the class callable
    (for generic classes).
    """
    info = itype.type
    if override_info:
        info = override_info

    fullname = f"{info.fullname}.{name}"
    hook = mx.chk.plugin.get_class_attribute_hook(fullname)

    node = info.get(name)
    if not node:
        if itype.extra_attrs and name in itype.extra_attrs.attrs:
            # For modules use direct symbol table lookup.
            if not itype.extra_attrs.mod_name:
                return itype.extra_attrs.attrs[name]
        if info.fallback_to_any or info.meta_fallback_to_any:
            return apply_class_attr_hook(mx, hook, AnyType(TypeOfAny.special_form))
        return None

    if (
        isinstance(node.node, Var)
        and not node.node.is_classvar
        and not hook
        and mcs_fallback.type.get(name)
    ):
        # If the same attribute is declared on the metaclass and the class but with different types,
        # and the attribute on the class is not a ClassVar,
        # the type of the attribute on the metaclass should take priority
        # over the type of the attribute on the class,
        # when the attribute is being accessed from the class object itself.
        #
        # Return `None` here to signify that the name should be looked up
        # on the class object itself rather than the instance.
        return None

    mx.chk.warn_deprecated(node.node, mx.context)

    is_decorated = isinstance(node.node, Decorator)
    is_method = is_decorated or isinstance(node.node, FuncBase)
    if mx.is_lvalue and not mx.suppress_errors:
        if is_method:
            mx.msg.cant_assign_to_method(mx.context)
        if isinstance(node.node, TypeInfo):
            mx.fail(message_registry.CANNOT_ASSIGN_TO_TYPE)

    # Refuse class attribute access if slot defined
    if info.slots and name in info.slots:
        mx.fail(message_registry.CLASS_VAR_CONFLICTS_SLOTS.format(name))

    # If a final attribute was declared on `self` in `__init__`, then it
    # can't be accessed on the class object.
    if node.implicit and isinstance(node.node, Var) and node.node.is_final:
        mx.fail(message_registry.CANNOT_ACCESS_FINAL_INSTANCE_ATTR.format(node.node.name))

    # An assignment to final attribute on class object is also always an error,
    # independently of types.
    if mx.is_lvalue and not mx.chk.get_final_context():
        check_final_member(name, info, mx.msg, mx.context)

    if info.is_enum and not (mx.is_lvalue or is_decorated or is_method):
        enum_class_attribute_type = analyze_enum_class_attribute_access(itype, name, mx)
        if enum_class_attribute_type:
            return apply_class_attr_hook(mx, hook, enum_class_attribute_type)

    t = node.type
    if t:
        if isinstance(t, PartialType):
            symnode = node.node
            assert isinstance(symnode, Var)
            return apply_class_attr_hook(
                mx, hook, mx.chk.handle_partial_var_type(t, mx.is_lvalue, symnode, mx.context)
            )

        # Find the class where method/variable was defined.
        if isinstance(node.node, Decorator):
            super_info: TypeInfo | None = node.node.var.info
        elif isinstance(node.node, (Var, SYMBOL_FUNCBASE_TYPES)):
            super_info = node.node.info
        else:
            super_info = None

        # Map the type to how it would look as a defining class. For example:
        #     class C(Generic[T]): ...
        #     class D(C[Tuple[T, S]]): ...
        #     D[int, str].method()
        # Here itype is D[int, str], isuper is C[Tuple[int, str]].
        if not super_info:
            isuper = None
        else:
            isuper = map_instance_to_supertype(itype, super_info)

        if isinstance(node.node, Var):
            assert isuper is not None
            # Check if original variable type has type variables. For example:
            #     class C(Generic[T]):
            #         x: T
            #     C.x  # Error, ambiguous access
            #     C[int].x  # Also an error, since C[int] is same as C at runtime
            # Exception is Self type wrapped in ClassVar, that is safe.
            def_vars = set(node.node.info.defn.type_vars)
            if not node.node.is_classvar and node.node.info.self_type:
                def_vars.add(node.node.info.self_type)
            # TODO: should we include ParamSpec etc. here (i.e. use get_all_type_vars)?
            typ_vars = set(get_type_vars(t))
            if def_vars & typ_vars:
                # Exception: access on Type[...], including first argument of class methods is OK.
                if not isinstance(get_proper_type(mx.original_type), TypeType) or node.implicit:
                    if node.node.is_classvar:
                        message = message_registry.GENERIC_CLASS_VAR_ACCESS
                    else:
                        message = message_registry.GENERIC_INSTANCE_VAR_CLASS_ACCESS
                    mx.fail(message)
            t = expand_self_type_if_needed(t, mx, node.node, itype, is_class=True)
            # Erase non-mapped variables, but keep mapped ones, even if there is an error.
            # In the above example this means that we infer following types:
            #     C.x -> Any
            #     C[int].x -> int
            t = erase_typevars(expand_type_by_instance(t, isuper), {tv.id for tv in def_vars})

        is_classmethod = (is_decorated and cast(Decorator, node.node).func.is_class) or (
            isinstance(node.node, SYMBOL_FUNCBASE_TYPES) and node.node.is_class
        )
        t = get_proper_type(t)
        is_trivial_self = False
        if isinstance(node.node, Decorator):
            # Use fast path if there are trivial decorators like @classmethod or @property
            is_trivial_self = node.node.func.is_trivial_self and not node.node.decorators
        elif isinstance(node.node, (FuncDef, OverloadedFuncDef)):
            is_trivial_self = node.node.is_trivial_self
        if isinstance(t, FunctionLike) and is_classmethod and not is_trivial_self:
            t = check_self_arg(t, mx.self_type, False, mx.context, name, mx.msg)
        result = add_class_tvars(
            t,
            isuper,
            is_classmethod,
            mx,
            original_vars=original_vars,
            is_trivial_self=is_trivial_self,
        )
        # __set__ is not called on class objects.
        if not mx.is_lvalue:
            result = analyze_descriptor_access(result, mx)

        return apply_class_attr_hook(mx, hook, result)
    elif isinstance(node.node, Var):
        mx.not_ready_callback(name, mx.context)
        return AnyType(TypeOfAny.special_form)

    if isinstance(node.node, (TypeInfo, TypeAlias, MypyFile, TypeVarLikeExpr)):
        # TODO: should we apply class plugin here (similar to instance access)?
        return mx.chk.expr_checker.analyze_static_reference(node.node, mx.context, mx.is_lvalue)

    if is_decorated:
        assert isinstance(node.node, Decorator)
        if node.node.type:
            return apply_class_attr_hook(mx, hook, node.node.type)
        else:
            mx.not_ready_callback(name, mx.context)
            return AnyType(TypeOfAny.from_error)
    else:
        assert isinstance(node.node, SYMBOL_FUNCBASE_TYPES)
        typ = function_type(node.node, mx.named_type("builtins.function"))
        # Note: if we are accessing class method on class object, the cls argument is bound.
        # Annotated and/or explicit class methods go through other code paths above, for
        # unannotated implicit class methods we do this here.
        if node.node.is_class:
            typ = bind_self_fast(typ)
        return apply_class_attr_hook(mx, hook, typ)


def apply_class_attr_hook(
    mx: MemberContext, hook: Callable[[AttributeContext], Type] | None, result: Type
) -> Type | None:
    if hook:
        result = hook(
            AttributeContext(
                get_proper_type(mx.original_type), result, mx.is_lvalue, mx.context, mx.chk
            )
        )
    return result


def analyze_enum_class_attribute_access(
    itype: Instance, name: str, mx: MemberContext
) -> Type | None:
    # Skip these since Enum will remove it
    if name in EXCLUDED_ENUM_ATTRIBUTES:
        return report_missing_attribute(mx.original_type, itype, name, mx)
    # Dunders and private names are not Enum members
    if name.startswith("__") and name.replace("_", "") != "":
        return None

    node = itype.type.get(name)
    if node and node.type:
        proper = get_proper_type(node.type)
        # Support `A = nonmember(1)` function call and decorator.
        if (
            isinstance(proper, Instance)
            and proper.type.fullname == "enum.nonmember"
            and proper.args
        ):
            return proper.args[0]

    enum_literal = LiteralType(name, fallback=itype)
    return itype.copy_modified(last_known_value=enum_literal)


def analyze_typeddict_access(
    name: str, typ: TypedDictType, mx: MemberContext, override_info: TypeInfo | None
) -> Type:
    if name == "__setitem__":
        if isinstance(mx.context, IndexExpr):
            # Since we can get this during `a['key'] = ...`
            # it is safe to assume that the context is `IndexExpr`.
            item_type, key_names = mx.chk.expr_checker.visit_typeddict_index_expr(
                typ, mx.context.index, setitem=True
            )
            assigned_readonly_keys = typ.readonly_keys & key_names
            if assigned_readonly_keys and not mx.suppress_errors:
                mx.msg.readonly_keys_mutated(assigned_readonly_keys, context=mx.context)
        else:
            # It can also be `a.__setitem__(...)` direct call.
            # In this case `item_type` can be `Any`,
            # because we don't have args available yet.
            # TODO: check in `default` plugin that `__setitem__` is correct.
            item_type = AnyType(TypeOfAny.implementation_artifact)
        return CallableType(
            arg_types=[mx.chk.named_type("builtins.str"), item_type],
            arg_kinds=[ARG_POS, ARG_POS],
            arg_names=[None, None],
            ret_type=NoneType(),
            fallback=mx.chk.named_type("builtins.function"),
            name=name,
        )
    elif name == "__delitem__":
        return CallableType(
            arg_types=[mx.chk.named_type("builtins.str")],
            arg_kinds=[ARG_POS],
            arg_names=[None],
            ret_type=NoneType(),
            fallback=mx.chk.named_type("builtins.function"),
            name=name,
        )
    return _analyze_member_access(name, typ.fallback, mx, override_info)


def add_class_tvars(
    t: ProperType,
    isuper: Instance | None,
    is_classmethod: bool,
    mx: MemberContext,
    original_vars: Sequence[TypeVarLikeType] | None = None,
    is_trivial_self: bool = False,
) -> Type:
    """Instantiate type variables during analyze_class_attribute_access,
    e.g T and Q in the following:

    class A(Generic[T]):
        @classmethod
        def foo(cls: Type[Q]) -> Tuple[T, Q]: ...

    class B(A[str]): pass
    B.foo()

    Args:
        t: Declared type of the method (or property)
        isuper: Current instance mapped to the superclass where method was defined, this
            is usually done by map_instance_to_supertype()
        is_classmethod: True if this method is decorated with @classmethod
        original_vars: Type variables of the class callable on which the method was accessed
        is_trivial_self: if True, we can use fast path for bind_self().
    Returns:
        Expanded method type with added type variables (when needed).
    """
    # TODO: verify consistency between Q and T

    # We add class type variables if the class method is accessed on class object
    # without applied type arguments, this matches the behavior of __init__().
    # For example (continuing the example in docstring):
    #     A       # The type of callable is def [T] () -> A[T], _not_ def () -> A[Any]
    #     A[int]  # The type of callable is def () -> A[int]
    # and
    #     A.foo       # The type is generic def [T] () -> Tuple[T, A[T]]
    #     A[int].foo  # The type is non-generic def () -> Tuple[int, A[int]]
    #
    # This behaviour is useful for defining alternative constructors for generic classes.
    # To achieve such behaviour, we add the class type variables that are still free
    # (i.e. appear in the return type of the class object on which the method was accessed).
    if isinstance(t, CallableType):
        tvars = original_vars if original_vars is not None else []
        if not mx.preserve_type_var_ids:
            t = freshen_all_functions_type_vars(t)
        if is_classmethod:
            if is_trivial_self:
                t = bind_self_fast(t, mx.self_type)
            else:
                t = bind_self(t, mx.self_type, is_classmethod=True)
        if isuper is not None:
            t = expand_type_by_instance(t, isuper)
        freeze_all_type_vars(t)
        return t.copy_modified(variables=list(tvars) + list(t.variables))
    elif isinstance(t, Overloaded):
        return Overloaded(
            [
                cast(
                    CallableType,
                    add_class_tvars(item, isuper, is_classmethod, mx, original_vars=original_vars),
                )
                for item in t.items
            ]
        )
    if isuper is not None:
        t = expand_type_by_instance(t, isuper)
    return t


def analyze_decorator_or_funcbase_access(
    defn: Decorator | FuncBase, itype: Instance, name: str, mx: MemberContext
) -> Type:
    """Analyzes the type behind method access.

    The function itself can possibly be decorated.
    See: https://github.com/python/mypy/issues/10409
    """
    if isinstance(defn, Decorator):
        return analyze_var(name, defn.var, itype, mx)
    typ = function_type(defn, mx.chk.named_type("builtins.function"))
    is_trivial_self = False
    if isinstance(defn, Decorator):
        # Use fast path if there are trivial decorators like @classmethod or @property
        is_trivial_self = defn.func.is_trivial_self and not defn.decorators
    elif isinstance(defn, (FuncDef, OverloadedFuncDef)):
        is_trivial_self = defn.is_trivial_self
    if is_trivial_self:
        return bind_self_fast(typ, mx.self_type)
    typ = check_self_arg(typ, mx.self_type, defn.is_class, mx.context, name, mx.msg)
    return bind_self(typ, original_type=mx.self_type, is_classmethod=defn.is_class)


F = TypeVar("F", bound=FunctionLike)


def bind_self_fast(method: F, original_type: Type | None = None) -> F:
    """Return a copy of `method`, with the type of its first parameter (usually
    self or cls) bound to original_type.

    This is a faster version of mypy.typeops.bind_self() that can be used for methods
    with trivial self/cls annotations.
    """
    if isinstance(method, Overloaded):
        items = [bind_self_fast(c, original_type) for c in method.items]
        return cast(F, Overloaded(items))
    assert isinstance(method, CallableType)
    func: CallableType = method
    if not func.arg_types:
        # Invalid method, return something.
        return method
    if func.arg_kinds[0] in (ARG_STAR, ARG_STAR2):
        # See typeops.py for details.
        return method
    original_type = get_proper_type(original_type)
    if isinstance(original_type, CallableType) and original_type.is_type_obj():
        original_type = TypeType.make_normalized(original_type.ret_type)
    res = func.copy_modified(
        arg_types=func.arg_types[1:],
        arg_kinds=func.arg_kinds[1:],
        arg_names=func.arg_names[1:],
        is_bound=True,
    )
    return cast(F, res)


def has_operator(typ: Type, op_method: str, named_type: Callable[[str], Instance]) -> bool:
    """Does type have operator with the given name?

    Note: this follows the rules for operator access, in particular:
    * __getattr__ is not considered
    * for class objects we only look in metaclass
    * instance level attributes (i.e. extra_attrs) are not considered
    """
    # This is much faster than analyze_member_access, and so using
    # it first as a filter is important for performance. This is mostly relevant
    # in situations where we can't expect that method is likely present,
    # e.g. for __OP__ vs __rOP__.
    typ = get_proper_type(typ)

    if isinstance(typ, TypeVarLikeType):
        typ = typ.values_or_bound()
    if isinstance(typ, AnyType):
        return True
    if isinstance(typ, UnionType):
        return all(has_operator(x, op_method, named_type) for x in typ.relevant_items())
    if isinstance(typ, FunctionLike) and typ.is_type_obj():
        return typ.fallback.type.has_readable_member(op_method)
    if isinstance(typ, TypeType):
        # Type[Union[X, ...]] is always normalized to Union[Type[X], ...],
        # so we don't need to care about unions here, but we need to care about
        # Type[T], where upper bound of T is a union.
        item = typ.item
        if isinstance(item, TypeVarType):
            item = item.values_or_bound()
        if isinstance(item, UnionType):
            return all(meta_has_operator(x, op_method, named_type) for x in item.relevant_items())
        return meta_has_operator(item, op_method, named_type)
    return instance_fallback(typ, named_type).type.has_readable_member(op_method)


def instance_fallback(typ: ProperType, named_type: Callable[[str], Instance]) -> Instance:
    if isinstance(typ, Instance):
        return typ
    if isinstance(typ, TupleType):
        return tuple_fallback(typ)
    if isinstance(typ, (LiteralType, TypedDictType)):
        return typ.fallback
    return named_type("builtins.object")


def meta_has_operator(item: Type, op_method: str, named_type: Callable[[str], Instance]) -> bool:
    item = get_proper_type(item)
    if isinstance(item, AnyType):
        return True
    item = instance_fallback(item, named_type)
    meta = item.type.metaclass_type or named_type("builtins.type")
    return meta.type.has_readable_member(op_method)
