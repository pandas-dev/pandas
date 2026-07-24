"""Transform mypy AST functions to IR (and related things).

Normal functions are translated into a list of basic blocks
containing various IR ops (defined in mypyc.ir.ops).

This also deals with generators, async functions and nested
functions. All of these are transformed into callable classes. These
have a custom __call__ method that implements the call, and state, such
as an environment containing non-local variables, is stored in the
instance of the callable class.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from typing import NamedTuple

from mypy.nodes import (
    ArgKind,
    ClassDef,
    Decorator,
    FuncBase,
    FuncDef,
    FuncItem,
    LambdaExpr,
    OverloadedFuncDef,
    TypeInfo,
    Var,
)
from mypy.types import CallableType, Type, UnboundType, get_proper_type
from mypyc.common import FAST_PREFIX, LAMBDA_NAME, PROPSET_PREFIX, SELF_NAME
from mypyc.ir.class_ir import ClassIR, NonExtClassInfo
from mypyc.ir.func_ir import (
    FUNC_CLASSMETHOD,
    FUNC_NORMAL,
    FUNC_STATICMETHOD,
    FuncDecl,
    FuncIR,
    FuncSignature,
    RuntimeArg,
)
from mypyc.ir.ops import (
    BasicBlock,
    ComparisonOp,
    GetAttr,
    Integer,
    LoadLiteral,
    Register,
    Return,
    SetAttr,
    Unbox,
    Unreachable,
    Value,
)
from mypyc.ir.rtypes import (
    RInstance,
    bool_rprimitive,
    c_int_rprimitive,
    dict_rprimitive,
    int_rprimitive,
    object_rprimitive,
)
from mypyc.irbuild.builder import IRBuilder, calculate_arg_defaults, gen_arg_defaults
from mypyc.irbuild.callable_class import (
    add_call_to_callable_class,
    add_get_to_callable_class,
    instantiate_callable_class,
    setup_callable_class,
)
from mypyc.irbuild.context import FuncInfo, GeneratorClass
from mypyc.irbuild.env_class import (
    add_vars_to_env,
    finalize_env_class,
    load_env_registers,
    setup_env_class,
)
from mypyc.irbuild.generator import gen_generator_func, gen_generator_func_body
from mypyc.irbuild.targets import AssignmentTarget
from mypyc.primitives.dict_ops import (
    dict_get_method_with_none,
    dict_new_op,
    exact_dict_set_item_op,
)
from mypyc.primitives.generic_ops import generic_getattr, generic_setattr, py_setattr_op
from mypyc.primitives.misc_ops import register_function
from mypyc.sametype import is_same_method_signature, is_same_type

# Top-level transform functions


def transform_func_def(builder: IRBuilder, fdef: FuncDef) -> None:
    sig = builder.mapper.fdef_to_sig(fdef, builder.options.strict_dunders_typing)
    func_ir, func_reg = gen_func_item(builder, fdef, fdef.name, sig)

    # If the function that was visited was a nested function, then either look it up in our
    # current environment or define it if it was not already defined.
    if func_reg:
        builder.assign(get_func_target(builder, fdef), func_reg, fdef.line)
    maybe_insert_into_registry_dict(builder, fdef)
    builder.add_function(func_ir, fdef.line)


def transform_overloaded_func_def(builder: IRBuilder, o: OverloadedFuncDef) -> None:
    # Handle regular overload case
    assert o.impl
    builder.accept(o.impl)


def transform_decorator(builder: IRBuilder, dec: Decorator) -> None:
    sig = builder.mapper.fdef_to_sig(dec.func, builder.options.strict_dunders_typing)
    func_ir, func_reg = gen_func_item(builder, dec.func, dec.func.name, sig)
    decorated_func: Value | None = None
    if func_reg:
        decorated_func = load_decorated_func(builder, dec.func, func_reg)
        builder.assign(get_func_target(builder, dec.func), decorated_func, dec.func.line)
    # If the prebuild pass didn't put this function in the function to decorators map (for example
    # if this is a registered singledispatch implementation with no other decorators), we should
    # treat this function as a regular function, not a decorated function
    elif dec.func in builder.fdefs_to_decorators:
        # Obtain the function name in order to construct the name of the helper function.
        name = dec.func.fullname.split(".")[-1]

        # Load the callable object representing the non-decorated function, and decorate it.
        orig_func = builder.load_global_str(name, dec.line)
        decorated_func = load_decorated_func(builder, dec.func, orig_func)

    if decorated_func is not None:
        # Set the callable object representing the decorated function as a global.
        builder.call_c(
            exact_dict_set_item_op,
            [builder.load_globals_dict(), builder.load_str(dec.func.name), decorated_func],
            decorated_func.line,
        )

    maybe_insert_into_registry_dict(builder, dec.func)

    builder.functions.append(func_ir)


def transform_lambda_expr(builder: IRBuilder, expr: LambdaExpr) -> Value:
    typ = get_proper_type(builder.types[expr])
    assert isinstance(typ, CallableType), typ

    runtime_args = []
    for arg, arg_type in zip(expr.arguments, typ.arg_types):
        arg.variable.type = arg_type
        runtime_args.append(
            RuntimeArg(arg.variable.name, builder.type_to_rtype(arg_type), arg.kind)
        )
    ret_type = builder.type_to_rtype(typ.ret_type)

    fsig = FuncSignature(runtime_args, ret_type)

    fname = f"{LAMBDA_NAME}{builder.lambda_counter}"
    builder.lambda_counter += 1
    func_ir, func_reg = gen_func_item(builder, expr, fname, fsig)
    assert func_reg is not None

    builder.functions.append(func_ir)
    return func_reg


# Internal functions


def gen_func_item(
    builder: IRBuilder,
    fitem: FuncItem,
    name: str,
    sig: FuncSignature,
    cdef: ClassDef | None = None,
    make_ext_method: bool = False,
) -> tuple[FuncIR, Value | None]:
    """Generate and return the FuncIR for a given FuncDef.

    If the given FuncItem is a nested function, then we generate a
    callable class representing the function and use that instead of
    the actual function. if the given FuncItem contains a nested
    function, then we generate an environment class so that inner
    nested functions can access the environment of the given FuncDef.

    Consider the following nested function:

        def a() -> None:
            def b() -> None:
                def c() -> None:
                    return None
                return None
            return None

    The classes generated would look something like the following.

                has pointer to        +-------+
        +-------------------------->  | a_env |
        |                             +-------+
        |                                 ^
        |                                 | has pointer to
    +-------+     associated with     +-------+
    | b_obj |   ------------------->  | b_env |
    +-------+                         +-------+
                                          ^
                                          |
    +-------+         has pointer to      |
    | c_obj |   --------------------------+
    +-------+
    """

    # TODO: do something about abstract methods.

    func_reg: Value | None = None

    # We treat lambdas as always being nested because we always generate
    # a class for lambdas, no matter where they are. (It would probably also
    # work to special case toplevel lambdas and generate a non-class function.)
    is_nested = fitem in builder.nested_fitems or isinstance(fitem, LambdaExpr)
    contains_nested = fitem in builder.encapsulating_funcs.keys()
    is_decorated = fitem in builder.fdefs_to_decorators
    is_singledispatch = fitem in builder.singledispatch_impls
    in_non_ext = False
    add_nested_funcs_to_env = has_nested_func_self_reference(builder, fitem)
    class_name = None
    if cdef:
        ir = builder.mapper.type_to_ir[cdef.info]
        in_non_ext = not ir.is_ext_class and not make_ext_method
        class_name = cdef.name

    if is_singledispatch:
        func_name = singledispatch_main_func_name(name)
    else:
        func_name = name

    fn_info = FuncInfo(
        fitem=fitem,
        name=func_name,
        class_name=class_name,
        namespace=gen_func_ns(builder),
        is_nested=is_nested,
        contains_nested=contains_nested,
        is_decorated=is_decorated,
        in_non_ext=in_non_ext,
        add_nested_funcs_to_env=add_nested_funcs_to_env,
    )
    is_generator = fn_info.is_generator
    builder.enter(fn_info, ret_type=sig.ret_type)

    if is_generator:
        fitem = builder.fn_info.fitem
        assert isinstance(fitem, FuncDef), fitem
        generator_class_ir = builder.mapper.fdef_to_generator[fitem]
        builder.fn_info.generator_class = GeneratorClass(generator_class_ir)

    # Functions that contain nested functions need an environment class to store variables that
    # are free in their nested functions. Generator functions need an environment class to
    # store a variable denoting the next instruction to be executed when the __next__ function
    # is called, along with all the variables inside the function itself.
    if contains_nested or (
        is_generator and not builder.fn_info.can_merge_generator_and_env_classes()
    ):
        setup_env_class(builder)

    if is_nested or in_non_ext:
        setup_callable_class(builder)

    if is_generator:
        # First generate a function that just constructs and returns a generator object.
        func_ir, func_reg = gen_generator_func(
            builder,
            lambda args, blocks, fn_info: gen_func_ir(
                builder, args, blocks, sig, fn_info, cdef, is_singledispatch
            ),
        )

        # Re-enter the FuncItem and visit the body of the function this time.
        gen_generator_func_body(builder, fn_info, func_reg)
    else:
        func_ir, func_reg = gen_func_body(builder, sig, cdef, is_singledispatch)

    if is_singledispatch:
        # add the generated main singledispatch function
        builder.functions.append(func_ir)
        # create the dispatch function
        assert isinstance(fitem, FuncDef), fitem
        return gen_dispatch_func_ir(builder, fitem, fn_info.name, name, sig)

    return func_ir, func_reg


def gen_func_body(
    builder: IRBuilder, sig: FuncSignature, cdef: ClassDef | None, is_singledispatch: bool
) -> tuple[FuncIR, Value | None]:
    load_env_registers(builder)
    gen_arg_defaults(builder)
    if builder.fn_info.contains_nested:
        finalize_env_class(builder)
    add_vars_to_env(builder)
    builder.accept(builder.fn_info.fitem.body)
    builder.maybe_add_implicit_return()

    # Hang on to the local symbol table for a while, since we use it
    # to calculate argument defaults below.
    symtable = builder.symtables[-1]

    args, _, blocks, ret_type, fn_info = builder.leave()

    func_ir, func_reg = gen_func_ir(builder, args, blocks, sig, fn_info, cdef, is_singledispatch)

    # Evaluate argument defaults in the surrounding scope, since we
    # calculate them *once* when the function definition is evaluated.
    calculate_arg_defaults(builder, fn_info, func_reg, symtable)
    return func_ir, func_reg


def has_nested_func_self_reference(builder: IRBuilder, fitem: FuncItem) -> bool:
    """Does a nested function contain a self-reference in its body?

    If a nested function only has references in the surrounding function,
    we don't need to add it to the environment.
    """
    if any(isinstance(sym, FuncBase) for sym in builder.free_variables.get(fitem, set())):
        return True
    return any(
        has_nested_func_self_reference(builder, nested)
        for nested in builder.encapsulating_funcs.get(fitem, [])
    )


def gen_func_ir(
    builder: IRBuilder,
    args: list[Register],
    blocks: list[BasicBlock],
    sig: FuncSignature,
    fn_info: FuncInfo,
    cdef: ClassDef | None,
    is_singledispatch_main_func: bool = False,
) -> tuple[FuncIR, Value | None]:
    """Generate the FuncIR for a function.

    This takes the basic blocks and function info of a particular
    function and returns the IR. If the function is nested,
    also returns the register containing the instance of the
    corresponding callable class.
    """
    func_reg: Value | None = None
    if fn_info.is_nested or fn_info.in_non_ext:
        func_ir = add_call_to_callable_class(builder, args, blocks, sig, fn_info)
        add_get_to_callable_class(builder, fn_info)
        func_reg = instantiate_callable_class(builder, fn_info)
    else:
        fitem = fn_info.fitem
        assert isinstance(fitem, FuncDef), fitem
        func_decl = builder.mapper.func_to_decl[fitem]
        if cdef and fn_info.name == FAST_PREFIX + func_decl.name:
            # Special-cased version of a method has a separate FuncDecl, use that one.
            func_decl = builder.mapper.type_to_ir[cdef.info].method_decls[fn_info.name]
        if fn_info.is_decorated or is_singledispatch_main_func:
            class_name = None if cdef is None else cdef.name
            func_decl = FuncDecl(
                fn_info.name,
                class_name,
                builder.module_name,
                sig,
                func_decl.kind,
                is_prop_getter=func_decl.is_prop_getter,
                is_prop_setter=func_decl.is_prop_setter,
                is_generator=func_decl.is_generator,
                is_coroutine=func_decl.is_coroutine,
                implicit=func_decl.implicit,
                internal=func_decl.internal,
            )
            func_ir = FuncIR(func_decl, args, blocks, fitem.line, traceback_name=fitem.name)
        else:
            func_ir = FuncIR(func_decl, args, blocks, fitem.line, traceback_name=fitem.name)
    return (func_ir, func_reg)


def generate_getattr_wrapper(builder: IRBuilder, cdef: ClassDef, getattr: FuncDef) -> None:
    """
    Generate a wrapper function for __getattr__ that can be put into the tp_getattro slot.
    The wrapper takes one argument besides self which is the attribute name.
    It first checks if the name matches any of the attributes of this class.
    If it does, it returns that attribute. If none match, it calls __getattr__.

    __getattr__ is not supported in classes that allow interpreted subclasses because the
    tp_getattro slot is inherited by subclasses and if the subclass overrides __getattr__,
    the override would be ignored in our wrapper. TODO: To support this, the wrapper would
    have to check type of self and if it's not the compiled class, resolve "__getattr__" against
    the type at runtime and call the returned method, like _Py_slot_tp_getattr_hook in cpython.

    __getattr__ is not supported in classes which inherit from non-native classes because those
    have __dict__ which currently has some strange interactions when class attributes and
    variables are assigned through __dict__ vs. through regular attribute access. Allowing
    __getattr__ on top of that could be problematic.
    """
    name = getattr.name + "__wrapper"
    ir = builder.mapper.type_to_ir[cdef.info]
    line = getattr.line

    error_base = f'"__getattr__" not supported in class "{cdef.name}" because '
    if ir.allow_interpreted_subclasses:
        builder.error(error_base + "it allows interpreted subclasses", line)
    if ir.inherits_python:
        builder.error(error_base + "it inherits from a non-native class", line)

    with builder.enter_method(ir, name, object_rprimitive, internal=True):
        attr_arg = builder.add_argument("attr", object_rprimitive)
        generic_getattr_result = builder.call_c(generic_getattr, [builder.self(), attr_arg], line)

        return_generic, call_getattr = BasicBlock(), BasicBlock()
        null = Integer(0, object_rprimitive, line)
        got_generic = builder.add(
            ComparisonOp(generic_getattr_result, null, ComparisonOp.NEQ, line)
        )
        builder.add_bool_branch(got_generic, return_generic, call_getattr)

        builder.activate_block(return_generic)
        builder.add(Return(generic_getattr_result, line))

        builder.activate_block(call_getattr)
        # No attribute matched so call user-provided __getattr__.
        getattr_result = builder.gen_method_call(
            builder.self(), getattr.name, [attr_arg], object_rprimitive, line
        )
        builder.add(Return(getattr_result, line))


def generate_setattr_wrapper(builder: IRBuilder, cdef: ClassDef, setattr: FuncDef) -> None:
    """
    Generate a wrapper function for __setattr__ that can be put into the tp_setattro slot.
    The wrapper takes two arguments besides self - attribute name and the new value.
    Returns 0 on success and -1 on failure. Restrictions are similar to the __getattr__
    wrapper above.

    The wrapper calls the user-defined __setattr__ when the value to set is not NULL.
    When it's NULL, this means that the call to tp_setattro comes from a del statement,
    so it calls __delattr__ instead. If __delattr__ is not overridden in the native class,
    this will call the base implementation in object which doesn't work without __dict__.
    """
    name = setattr.name + "__wrapper"
    ir = builder.mapper.type_to_ir[cdef.info]
    line = setattr.line

    error_base = f'"__setattr__" not supported in class "{cdef.name}" because '
    if ir.allow_interpreted_subclasses:
        builder.error(error_base + "it allows interpreted subclasses", line)
    if ir.inherits_python:
        builder.error(error_base + "it inherits from a non-native class", line)

    with builder.enter_method(ir, name, c_int_rprimitive, internal=True):
        attr_arg = builder.add_argument("attr", object_rprimitive)
        value_arg = builder.add_argument("value", object_rprimitive)

        call_delattr, call_setattr = BasicBlock(), BasicBlock()
        null = Integer(0, object_rprimitive, line)
        is_delattr = builder.add(ComparisonOp(value_arg, null, ComparisonOp.EQ, line))
        builder.add_bool_branch(is_delattr, call_delattr, call_setattr)

        builder.activate_block(call_delattr)
        delattr_symbol = cdef.info.get("__delattr__")
        delattr = delattr_symbol.node if delattr_symbol else None
        delattr_override = delattr is not None and not delattr.fullname.startswith("builtins.")
        if delattr_override:
            builder.gen_method_call(builder.self(), "__delattr__", [attr_arg], None, line)
        else:
            # Call internal function that cpython normally calls when deleting an attribute.
            # Cannot call object.__delattr__ here because it calls PyObject_SetAttr internally
            # which in turn calls our wrapper and recurses infinitely.
            # Note that since native classes don't have __dict__, this will raise AttributeError
            # for dynamic attributes.
            builder.call_c(generic_setattr, [builder.self(), attr_arg, null], line)
        builder.add(Return(Integer(0, c_int_rprimitive), line))

        builder.activate_block(call_setattr)
        builder.gen_method_call(builder.self(), setattr.name, [attr_arg, value_arg], None, line)
        builder.add(Return(Integer(0, c_int_rprimitive), line))


def handle_ext_method(builder: IRBuilder, cdef: ClassDef, fdef: FuncDef) -> None:
    # Perform the function of visit_method for methods inside extension classes.
    name = fdef.name
    class_ir = builder.mapper.type_to_ir[cdef.info]
    sig = builder.mapper.fdef_to_sig(fdef, builder.options.strict_dunders_typing)
    func_ir, func_reg = gen_func_item(builder, fdef, name, sig, cdef)
    builder.functions.append(func_ir)

    if is_decorated(builder, fdef):
        # Obtain the function name in order to construct the name of the helper function.
        _, _, name = fdef.fullname.rpartition(".")
        # Read the PyTypeObject representing the class, get the callable object
        # representing the non-decorated method
        typ = builder.load_native_type_object(cdef.fullname)
        orig_func = builder.py_get_attr(typ, name, fdef.line)

        # Decorate the non-decorated method
        decorated_func = load_decorated_func(builder, fdef, orig_func)

        # Set the callable object representing the decorated method as an attribute of the
        # extension class.
        builder.primitive_op(
            py_setattr_op, [typ, builder.load_str(name), decorated_func], fdef.line
        )

    if fdef.is_property:
        # If there is a property setter, it will be processed after the getter,
        # We populate the optional setter field with none for now.
        assert name not in class_ir.properties
        class_ir.properties[name] = (func_ir, None)

    elif fdef in builder.prop_setters:
        # The respective property getter must have been processed already
        assert name in class_ir.properties
        getter_ir, _ = class_ir.properties[name]
        class_ir.properties[name] = (getter_ir, func_ir)

    class_ir.methods[func_ir.decl.name] = func_ir

    # If this overrides a parent class method with a different type, we need
    # to generate a glue method to mediate between them.
    for base in class_ir.mro[1:]:
        if (
            name in base.method_decls
            and name != "__init__"
            and not is_same_method_signature(
                class_ir.method_decls[name].sig, base.method_decls[name].sig
            )
        ):
            # TODO: Support contravariant subtyping in the input argument for
            # property setters. Need to make a special glue method for handling this,
            # similar to gen_glue_property.

            f = gen_glue(builder, base.method_decls[name].sig, func_ir, class_ir, base, fdef)
            class_ir.glue_methods[(base, name)] = f
            builder.functions.append(f)

    # If the class allows interpreted children, create glue
    # methods that dispatch via the Python API. These will go in a
    # "shadow vtable" that will be assigned to interpreted
    # children.
    if class_ir.allow_interpreted_subclasses:
        f = gen_glue(builder, func_ir.sig, func_ir, class_ir, class_ir, fdef, do_py_ops=True)
        # Use func_ir.decl.name (unique) rather than fdef.name, because for properties
        # the getter and setter share the same fdef.name but have distinct decl names
        # (e.g. "prop" vs "__mypyc_setter__prop"). Using fdef.name would cause the
        # setter's glue to overwrite the getter's glue in the shadow vtable.
        class_ir.glue_methods[(class_ir, func_ir.decl.name)] = f
        builder.functions.append(f)

    if fdef.name == "__getattr__":
        generate_getattr_wrapper(builder, cdef, fdef)
    elif fdef.name == "__setattr__":
        generate_setattr_wrapper(builder, cdef, fdef)
    elif fdef.name == "__delattr__":
        setattr = cdef.info.get("__setattr__")
        if not setattr or not setattr.node or setattr.node.fullname.startswith("builtins."):
            builder.error(
                '"__delattr__" supported only in classes that also override "__setattr__", '
                + "or inherit from a native class that overrides it.",
                fdef.line,
            )


def handle_non_ext_method(
    builder: IRBuilder, non_ext: NonExtClassInfo, cdef: ClassDef, fdef: FuncDef
) -> None:
    # Perform the function of visit_method for methods inside non-extension classes.
    name = fdef.name
    sig = builder.mapper.fdef_to_sig(fdef, builder.options.strict_dunders_typing)
    func_ir, func_reg = gen_func_item(builder, fdef, name, sig, cdef)
    assert func_reg is not None
    builder.functions.append(func_ir)

    if is_decorated(builder, fdef):
        # The undecorated method is a generated callable class
        orig_func = func_reg
        func_reg = load_decorated_func(builder, fdef, orig_func)

    # TODO: Support property setters in non-extension classes
    if fdef.is_property:
        prop = builder.load_module_attr_by_fullname("builtins.property", fdef.line)
        func_reg = builder.py_call(prop, [func_reg], fdef.line)

    elif builder.mapper.func_to_decl[fdef].kind == FUNC_CLASSMETHOD:
        cls_meth = builder.load_module_attr_by_fullname("builtins.classmethod", fdef.line)
        func_reg = builder.py_call(cls_meth, [func_reg], fdef.line)

    elif builder.mapper.func_to_decl[fdef].kind == FUNC_STATICMETHOD:
        stat_meth = builder.load_module_attr_by_fullname("builtins.staticmethod", fdef.line)
        func_reg = builder.py_call(stat_meth, [func_reg], fdef.line)

    builder.add_to_non_ext_dict(non_ext, name, func_reg, fdef.line)

    # If we identified that this non-extension class method can be special-cased for
    # direct access during prepare phase, generate a "static" version of it.
    class_ir = builder.mapper.type_to_ir[cdef.info]
    name = FAST_PREFIX + fdef.name
    if name in class_ir.method_decls:
        func_ir, func_reg = gen_func_item(builder, fdef, name, sig, cdef, make_ext_method=True)
        class_ir.methods[name] = func_ir
        builder.functions.append(func_ir)


def gen_func_ns(builder: IRBuilder) -> str:
    """Generate a namespace for a nested function using its outer function names."""
    return "_".join(
        info.name + ("" if not info.class_name else "_" + info.class_name)
        for info in builder.fn_infos
        if info.name and info.name != "<module>"
    )


def load_decorated_func(builder: IRBuilder, fdef: FuncDef, orig_func_reg: Value) -> Value:
    """Apply decorators to a function.

    Given a decorated FuncDef and an instance of the callable class
    representing that FuncDef, apply the corresponding decorator
    functions on that decorated FuncDef and return the decorated
    function.
    """
    if not is_decorated(builder, fdef):
        # If there are no decorators associated with the function, then just return the
        # original function.
        return orig_func_reg

    decorators = builder.fdefs_to_decorators[fdef]
    func_reg = orig_func_reg
    for d in reversed(decorators):
        decorator = d.accept(builder.visitor)
        assert isinstance(decorator, Value), decorator
        func_reg = builder.py_call(decorator, [func_reg], func_reg.line)
    return func_reg


def is_decorated(builder: IRBuilder, fdef: FuncDef) -> bool:
    return fdef in builder.fdefs_to_decorators


def gen_glue(
    builder: IRBuilder,
    base_sig: FuncSignature,
    target: FuncIR,
    cls: ClassIR,
    base: ClassIR,
    fdef: FuncItem,
    *,
    do_py_ops: bool = False,
) -> FuncIR:
    """Generate glue methods that mediate between different method types in subclasses.

    Works on both properties and methods. See gen_glue_methods below
    for more details.

    If do_py_ops is True, then the glue methods should use generic
    C API operations instead of direct calls, to enable generating
    "shadow" glue methods that work with interpreted subclasses.
    """
    if fdef.is_property:
        return gen_glue_property(builder, base_sig, target, cls, base, fdef.line, do_py_ops)
    if do_py_ops and target.name.startswith(PROPSET_PREFIX):
        return gen_glue_property_setter(builder, base_sig, target, cls, base, fdef.line)
    return gen_glue_method(builder, base_sig, target, cls, base, fdef.line, do_py_ops)


class ArgInfo(NamedTuple):
    args: list[Value]
    arg_names: list[str | None]
    arg_kinds: list[ArgKind]


def get_args(builder: IRBuilder, rt_args: Sequence[RuntimeArg], line: int) -> ArgInfo:
    # The environment operates on Vars, so we make some up
    fake_vars = [(Var(arg.name), arg.type) for arg in rt_args]
    args = [
        builder.read(builder.add_local_reg(var, type, is_arg=True), line)
        for var, type in fake_vars
    ]
    arg_names = [
        arg.name if arg.kind.is_named() or (arg.kind.is_optional() and not arg.pos_only) else None
        for arg in rt_args
    ]
    arg_kinds = [arg.kind for arg in rt_args]
    return ArgInfo(args, arg_names, arg_kinds)


def gen_glue_method(
    builder: IRBuilder,
    base_sig: FuncSignature,
    target: FuncIR,
    cls: ClassIR,
    base: ClassIR,
    line: int,
    do_pycall: bool,
) -> FuncIR:
    """Generate glue methods that mediate between different method types in subclasses.

    For example, if we have:

    class A:
        def f(builder: IRBuilder, x: int) -> object: ...

    then it is totally permissible to have a subclass

    class B(A):
        def f(builder: IRBuilder, x: object) -> int: ...

    since '(object) -> int' is a subtype of '(int) -> object' by the usual
    contra/co-variant function subtyping rules.

    The trickiness here is that int and object have different
    runtime representations in mypyc, so A.f and B.f have
    different signatures at the native C level. To deal with this,
    we need to generate glue methods that mediate between the
    different versions by coercing the arguments and return
    values.

    If do_pycall is True, then make the call using the C API
    instead of a native call.
    """
    check_native_override(builder, base_sig, target.decl.sig, line)

    builder.enter()
    builder.ret_types[-1] = base_sig.ret_type

    rt_args = list(base_sig.args)
    if target.decl.kind == FUNC_NORMAL:
        rt_args[0] = RuntimeArg(base_sig.args[0].name, RInstance(cls))

    arg_info = get_args(builder, rt_args, line)
    args, arg_kinds, arg_names = arg_info.args, arg_info.arg_kinds, arg_info.arg_names

    bitmap_args = None
    if base_sig.num_bitmap_args:
        args = args[: -base_sig.num_bitmap_args]
        arg_kinds = arg_kinds[: -base_sig.num_bitmap_args]
        arg_names = arg_names[: -base_sig.num_bitmap_args]
        bitmap_args = list(builder.builder.args[-base_sig.num_bitmap_args :])

    # We can do a passthrough *args/**kwargs with a native call, but if the
    # args need to get distributed out to arguments, we just let python handle it
    if any(kind.is_star() for kind in arg_kinds) and any(
        not arg.kind.is_star() for arg in target.decl.sig.args
    ):
        do_pycall = True

    if do_pycall:
        if target.decl.kind == FUNC_STATICMETHOD:
            # FIXME: this won't work if we can do interpreted subclasses
            first = builder.builder.get_native_type(cls)
            st = 0
        else:
            first = args[0]
            st = 1
        retval = builder.builder.py_method_call(
            first, target.name, args[st:], line, arg_kinds[st:], arg_names[st:]
        )
    else:
        retval = builder.builder.call(
            target.decl, args, arg_kinds, arg_names, line, bitmap_args=bitmap_args
        )
    retval = builder.coerce(retval, base_sig.ret_type, line)
    builder.add(Return(retval))

    arg_regs, _, blocks, ret_type, _ = builder.leave()
    if base_sig.num_bitmap_args:
        rt_args = rt_args[: -base_sig.num_bitmap_args]
    return FuncIR(
        FuncDecl(
            target.name + "__" + base.name + "_glue",
            cls.name,
            builder.module_name,
            FuncSignature(rt_args, ret_type),
            target.decl.kind,
            is_coroutine=target.decl.is_coroutine,
        ),
        arg_regs,
        blocks,
    )


def check_native_override(
    builder: IRBuilder, base_sig: FuncSignature, sub_sig: FuncSignature, line: int
) -> None:
    """Report an error if an override changes signature in unsupported ways.

    Glue methods can work around many signature changes but not all of them.
    """
    for base_arg, sub_arg in zip(base_sig.real_args(), sub_sig.real_args()):
        if base_arg.type.error_overlap:
            if not base_arg.optional and sub_arg.optional and base_sig.num_bitmap_args:
                # This would change the meanings of bits in the argument defaults
                # bitmap, which we don't support. We'd need to do tricky bit
                # manipulations to support this generally.
                builder.error(
                    "An argument with type "
                    + f'"{base_arg.type}" cannot be given a default value in a method override',
                    line,
                )
        if base_arg.type.error_overlap or sub_arg.type.error_overlap:
            if not is_same_type(base_arg.type, sub_arg.type):
                # This would change from signaling a default via an error value to
                # signaling a default via bitmap, which we don't support.
                builder.error(
                    "Incompatible argument type "
                    + f'"{sub_arg.type}" (base class has type "{base_arg.type}")',
                    line,
                )


def gen_glue_property(
    builder: IRBuilder,
    sig: FuncSignature,
    target: FuncIR,
    cls: ClassIR,
    base: ClassIR,
    line: int,
    do_pygetattr: bool,
) -> FuncIR:
    """Generate glue methods for properties that mediate between different subclass types.

    Similarly to methods, properties of derived types can be covariantly subtyped. Thus,
    properties also require glue. However, this only requires the return type to change.
    Further, instead of a method call, an attribute get is performed.

    If do_pygetattr is True, then get the attribute using the Python C
    API instead of a native call.
    """
    builder.enter()

    rt_arg = RuntimeArg(SELF_NAME, RInstance(cls))
    self_target = builder.add_self_to_env(cls)
    arg = builder.read(self_target, line)
    builder.ret_types[-1] = sig.ret_type
    if do_pygetattr:
        retval = builder.py_get_attr(arg, target.name, line)
    else:
        retval = builder.add(GetAttr(arg, target.name, line))
    retbox = builder.coerce(retval, sig.ret_type, line)
    builder.add(Return(retbox))

    args, _, blocks, return_type, _ = builder.leave()
    return FuncIR(
        FuncDecl(
            target.name + "__" + base.name + "_glue",
            cls.name,
            builder.module_name,
            FuncSignature([rt_arg], return_type),
        ),
        args,
        blocks,
    )


def gen_glue_property_setter(
    builder: IRBuilder, sig: FuncSignature, target: FuncIR, cls: ClassIR, base: ClassIR, line: int
) -> FuncIR:
    """Generate a shadow glue method for a property setter.

    For interpreted subclasses, property setters can't be called via the
    internal __mypyc_setter__<name> method. Instead, use Python's setattr
    to set the property via the standard descriptor protocol.
    """
    builder.enter()
    builder.ret_types[-1] = sig.ret_type

    rt_args = list(sig.args)
    rt_args[0] = RuntimeArg(sig.args[0].name, RInstance(cls))

    arg_info = get_args(builder, rt_args, line)
    args = arg_info.args

    self_arg = args[0]
    value_arg = args[1]

    # Extract the property name from "__mypyc_setter__<name>"
    assert target.name.startswith(PROPSET_PREFIX)
    prop_name = target.name[len(PROPSET_PREFIX) :]

    builder.primitive_op(
        py_setattr_op,
        [
            self_arg,
            builder.load_str(prop_name),
            builder.coerce(value_arg, object_rprimitive, line),
        ],
        line,
    )
    retval = builder.coerce(builder.none(), sig.ret_type, line)
    builder.add(Return(retval))

    arg_regs, _, blocks, return_type, _ = builder.leave()
    return FuncIR(
        FuncDecl(
            target.name + "__" + base.name + "_glue",
            cls.name,
            builder.module_name,
            FuncSignature(rt_args, return_type),
        ),
        arg_regs,
        blocks,
    )


def get_func_target(builder: IRBuilder, fdef: FuncDef) -> AssignmentTarget:
    """Given a FuncDef, return the target for the instance of its callable class.

    If the function was not already defined somewhere, then define it
    and add it to the current environment.
    """
    if orig := fdef.original_def:
        if isinstance(orig, Decorator):
            orig = orig.func
        # Get the target associated with the previously defined FuncDef.
        return builder.lookup(orig)

    if builder.fn_info.is_generator or builder.fn_info.add_nested_funcs_to_env:
        return builder.lookup(fdef)

    return builder.add_local_reg(fdef, object_rprimitive)


def load_type(builder: IRBuilder, typ: TypeInfo, unbounded_type: Type | None, line: int) -> Value:
    # typ.fullname contains the module where the class object was defined. However, it is possible
    # that the class object's module was not imported in the file currently being compiled. So, we
    # use unbounded_type.name (if provided by caller) to load the class object through one of the
    # imported modules.
    # Example: for `json.JSONDecoder`, typ.fullname is `json.decoder.JSONDecoder` but the Python
    # file may import `json` not `json.decoder`.
    # Another corner case: The Python file being compiled imports mod1 and has a type hint
    # `mod1.OuterClass.InnerClass`. But, mod1/__init__.py might import OuterClass like this:
    # `from mod2.mod3 import OuterClass`. In this case, typ.fullname is
    # `mod2.mod3.OuterClass.InnerClass` and `unbounded_type.name` is `mod1.OuterClass.InnerClass`.
    # So, we must use unbounded_type.name to load the class object.
    # See issue mypyc/mypyc#1087.
    if typ in builder.mapper.type_to_ir:
        class_ir = builder.mapper.type_to_ir[typ]
        class_obj = builder.builder.get_native_type(class_ir)
    elif builtin := builder.load_builtin(typ.fullname, line):
        class_obj = builtin
    elif isinstance(unbounded_type, UnboundType):
        path_parts = unbounded_type.name.split(".")
        class_obj = builder.load_global_str(path_parts[0], line)
        for attr in path_parts[1:]:
            class_obj = builder.py_get_attr(class_obj, attr, line)
    else:
        class_obj = builder.load_global_str(typ.name, line)

    return class_obj


def load_func(builder: IRBuilder, func_name: str, fullname: str | None, line: int) -> Value:
    if fullname and not fullname.startswith(builder.current_module):
        # we're calling a function in a different module

        # We can't use load_module_attr_by_fullname here because we need to load the function using
        # func_name, not the name specified by fullname (which can be different for underscore
        # function)
        module = fullname.rsplit(".")[0]
        loaded_module = builder.load_module(module)

        func = builder.py_get_attr(loaded_module, func_name, line)
    else:
        func = builder.load_global_str(func_name, line)
    return func


def generate_singledispatch_dispatch_function(
    builder: IRBuilder, main_singledispatch_function_name: str, fitem: FuncDef
) -> None:
    line = fitem.line
    current_func_decl = builder.mapper.func_to_decl[fitem]
    arg_info = get_args(builder, current_func_decl.sig.args, line)

    dispatch_func_obj = builder.self()

    arg_type = builder.builder.get_type_of_obj(arg_info.args[0], line)
    dispatch_cache = builder.builder.get_attr(
        dispatch_func_obj, "dispatch_cache", dict_rprimitive, line
    )
    call_find_impl, use_cache, call_func = BasicBlock(), BasicBlock(), BasicBlock()
    get_result = builder.primitive_op(dict_get_method_with_none, [dispatch_cache, arg_type], line)
    is_not_none = builder.translate_is_op(get_result, builder.none_object(), "is not", line)
    impl_to_use = Register(object_rprimitive)
    builder.add_bool_branch(is_not_none, use_cache, call_find_impl)

    builder.activate_block(use_cache)
    builder.assign(impl_to_use, get_result, line)
    builder.goto(call_func)

    builder.activate_block(call_find_impl)
    find_impl = builder.load_module_attr_by_fullname("functools._find_impl", line)
    registry = load_singledispatch_registry(builder, dispatch_func_obj, line)
    uncached_impl = builder.py_call(find_impl, [arg_type, registry], line)
    builder.call_c(exact_dict_set_item_op, [dispatch_cache, arg_type, uncached_impl], line)
    builder.assign(impl_to_use, uncached_impl, line)
    builder.goto(call_func)

    builder.activate_block(call_func)
    gen_calls_to_correct_impl(builder, impl_to_use, arg_info, fitem, line)


def gen_calls_to_correct_impl(
    builder: IRBuilder, impl_to_use: Value, arg_info: ArgInfo, fitem: FuncDef, line: int
) -> None:
    current_func_decl = builder.mapper.func_to_decl[fitem]

    def gen_native_func_call_and_return(fdef: FuncDef) -> None:
        func_decl = builder.mapper.func_to_decl[fdef]
        ret_val = builder.builder.call(
            func_decl, arg_info.args, arg_info.arg_kinds, arg_info.arg_names, line
        )
        coerced = builder.coerce(ret_val, current_func_decl.sig.ret_type, line)
        builder.add(Return(coerced))

    int_type_obj = builder.load_builtin("builtins.int", line)
    assert int_type_obj
    is_int = builder.builder.type_is_op(impl_to_use, int_type_obj, line)

    native_call, non_native_call = BasicBlock(), BasicBlock()
    builder.add_bool_branch(is_int, native_call, non_native_call)
    builder.activate_block(native_call)

    passed_id = builder.add(Unbox(impl_to_use, int_rprimitive, line))

    native_ids = get_native_impl_ids(builder, fitem)
    for impl, i in native_ids.items():
        call_impl, next_impl = BasicBlock(), BasicBlock()

        current_id = builder.load_int(i)
        cond = builder.binary_op(passed_id, current_id, "==", line)
        builder.add_bool_branch(cond, call_impl, next_impl)

        # Call the registered implementation
        builder.activate_block(call_impl)

        gen_native_func_call_and_return(impl)
        builder.activate_block(next_impl)

    # We've already handled all the possible integer IDs, so we should never get here
    builder.add(Unreachable())

    builder.activate_block(non_native_call)
    ret_val = builder.py_call(
        impl_to_use, arg_info.args, line, arg_info.arg_kinds, arg_info.arg_names
    )
    coerced = builder.coerce(ret_val, current_func_decl.sig.ret_type, line)
    builder.add(Return(coerced))


def gen_dispatch_func_ir(
    builder: IRBuilder, fitem: FuncDef, main_func_name: str, dispatch_name: str, sig: FuncSignature
) -> tuple[FuncIR, Value]:
    """Create a dispatch function (a function that checks the first argument type and dispatches
    to the correct implementation)
    """
    builder.enter(FuncInfo(fitem, dispatch_name))
    setup_callable_class(builder)
    builder.fn_info.callable_class.ir.attributes["registry"] = dict_rprimitive
    builder.fn_info.callable_class.ir.attributes["dispatch_cache"] = dict_rprimitive
    builder.fn_info.callable_class.ir.has_dict = True
    builder.fn_info.callable_class.ir.needs_getseters = True
    generate_singledispatch_callable_class_ctor(builder)

    generate_singledispatch_dispatch_function(builder, main_func_name, fitem)
    args, _, blocks, _, fn_info = builder.leave()
    dispatch_callable_class = add_call_to_callable_class(builder, args, blocks, sig, fn_info)
    builder.functions.append(dispatch_callable_class)
    add_get_to_callable_class(builder, fn_info)
    add_register_method_to_callable_class(builder, fn_info)
    func_reg = instantiate_callable_class(builder, fn_info)
    dispatch_func_ir = generate_dispatch_glue_native_function(
        builder, fitem, dispatch_callable_class.decl, dispatch_name
    )

    return dispatch_func_ir, func_reg


def generate_dispatch_glue_native_function(
    builder: IRBuilder, fitem: FuncDef, callable_class_decl: FuncDecl, dispatch_name: str
) -> FuncIR:
    line = fitem.line
    builder.enter()
    # We store the callable class in the globals dict for this function
    callable_class = builder.load_global_str(dispatch_name, line)
    decl = builder.mapper.func_to_decl[fitem]
    arg_info = get_args(builder, decl.sig.args, line)
    args = [callable_class] + arg_info.args
    arg_kinds = [ArgKind.ARG_POS] + arg_info.arg_kinds
    arg_names = arg_info.arg_names
    arg_names.insert(0, "self")
    ret_val = builder.builder.call(callable_class_decl, args, arg_kinds, arg_names, line)
    builder.add(Return(ret_val))
    arg_regs, _, blocks, _, fn_info = builder.leave()
    return FuncIR(decl, arg_regs, blocks)


def generate_singledispatch_callable_class_ctor(builder: IRBuilder) -> None:
    """Create an __init__ that sets registry and dispatch_cache to empty dicts"""
    line = builder.fn_info.fitem.line
    class_ir = builder.fn_info.callable_class.ir
    with builder.enter_method(class_ir, "__init__", bool_rprimitive):
        empty_dict = builder.call_c(dict_new_op, [], line)
        builder.add(SetAttr(builder.self(), "registry", empty_dict, line))
        cache_dict = builder.call_c(dict_new_op, [], line)
        dispatch_cache_str = builder.load_str("dispatch_cache")
        # use the py_setattr_op instead of SetAttr so that it also gets added to our __dict__
        builder.primitive_op(py_setattr_op, [builder.self(), dispatch_cache_str, cache_dict], line)
        # the generated C code seems to expect that __init__ returns a char, so just return 1
        builder.add(Return(Integer(1, bool_rprimitive, line), line))


def add_register_method_to_callable_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    line = fn_info.fitem.line
    with builder.enter_method(fn_info.callable_class.ir, "register", object_rprimitive):
        cls_arg = builder.add_argument("cls", object_rprimitive)
        func_arg = builder.add_argument("func", object_rprimitive, ArgKind.ARG_OPT)
        ret_val = builder.call_c(register_function, [builder.self(), cls_arg, func_arg], line)
        builder.add(Return(ret_val, line))


def load_singledispatch_registry(builder: IRBuilder, dispatch_func_obj: Value, line: int) -> Value:
    return builder.builder.get_attr(dispatch_func_obj, "registry", dict_rprimitive, line)


def singledispatch_main_func_name(orig_name: str) -> str:
    return f"__mypyc_singledispatch_main_function_{orig_name}__"


def maybe_insert_into_registry_dict(builder: IRBuilder, fitem: FuncDef) -> None:
    line = fitem.line
    is_singledispatch_main_func = fitem in builder.singledispatch_impls
    # dict of singledispatch_func to list of register_types (fitem is the function to register)
    to_register: defaultdict[FuncDef, list[TypeInfo]] = defaultdict(list)
    for main_func, impls in builder.singledispatch_impls.items():
        for dispatch_type, impl in impls:
            if fitem == impl:
                to_register[main_func].append(dispatch_type)

    if not to_register and not is_singledispatch_main_func:
        return

    if is_singledispatch_main_func:
        main_func_name = singledispatch_main_func_name(fitem.name)
        main_func_obj = load_func(builder, main_func_name, fitem.fullname, line)

        loaded_object_type = builder.load_module_attr_by_fullname("builtins.object", line)
        registry_dict = builder.builder.make_dict([(loaded_object_type, main_func_obj)], line)

        dispatch_func_obj = builder.load_global_str(fitem.name, line)
        builder.primitive_op(
            py_setattr_op, [dispatch_func_obj, builder.load_str("registry"), registry_dict], line
        )

    for singledispatch_func, types in to_register.items():
        # TODO: avoid recomputing the native IDs for all the functions every time we find a new
        # function
        native_ids = get_native_impl_ids(builder, singledispatch_func)
        if fitem not in native_ids:
            to_insert = load_func(builder, fitem.name, fitem.fullname, line)
        else:
            current_id = native_ids[fitem]
            load_literal = LoadLiteral(current_id, object_rprimitive)
            to_insert = builder.add(load_literal)
        # TODO: avoid reloading the registry here if we just created it
        dispatch_func_obj = load_func(
            builder, singledispatch_func.name, singledispatch_func.fullname, line
        )
        registry = load_singledispatch_registry(builder, dispatch_func_obj, line)
        for typ in types:
            loaded_type = load_type(builder, typ, None, line)
            builder.call_c(exact_dict_set_item_op, [registry, loaded_type, to_insert], line)
        dispatch_cache = builder.builder.get_attr(
            dispatch_func_obj, "dispatch_cache", dict_rprimitive, line
        )
        builder.gen_method_call(dispatch_cache, "clear", [], None, line)


def get_native_impl_ids(builder: IRBuilder, singledispatch_func: FuncDef) -> dict[FuncDef, int]:
    """Return a dict of registered implementation to native implementation ID for all
    implementations
    """
    impls = builder.singledispatch_impls[singledispatch_func]
    return {impl: i for i, (typ, impl) in enumerate(impls) if not is_decorated(builder, impl)}


def gen_property_getter_ir(
    builder: IRBuilder, func_decl: FuncDecl, cdef: ClassDef, is_trait: bool
) -> FuncIR:
    """Generate an implicit trivial property getter for an attribute.

    These are used if an attribute can also be accessed as a property.
    """
    name = func_decl.name
    builder.enter(name)
    self_reg = builder.add_argument("self", func_decl.sig.args[0].type)
    line = func_decl._line or -1
    if not is_trait:
        value = builder.builder.get_attr(self_reg, name, func_decl.sig.ret_type, line)
        builder.add(Return(value, line))
    else:
        builder.add(Unreachable())
    args, _, blocks, ret_type, fn_info = builder.leave()
    return FuncIR(func_decl, args, blocks)


def gen_property_setter_ir(
    builder: IRBuilder, func_decl: FuncDecl, cdef: ClassDef, is_trait: bool
) -> FuncIR:
    """Generate an implicit trivial property setter for an attribute.

    These are used if an attribute can also be accessed as a property.
    """
    name = func_decl.name
    builder.enter(name)
    self_reg = builder.add_argument("self", func_decl.sig.args[0].type)
    value_reg = builder.add_argument("value", func_decl.sig.args[1].type)
    assert name.startswith(PROPSET_PREFIX)
    attr_name = name[len(PROPSET_PREFIX) :]
    line = func_decl._line or -1
    if not is_trait:
        builder.add(SetAttr(self_reg, attr_name, value_reg, line))
    builder.add(Return(builder.none(), line))
    args, _, blocks, ret_type, fn_info = builder.leave()
    return FuncIR(func_decl, args, blocks, line)
