"""Transform mypy AST functions to IR (and related things).

Normal functions are translated into a list of basic blocks
containing various IR ops (defined in mypyc.ir.ops).

This also deals with generators, async functions and nested
functions. All of these are transformed into callable classes. These
have a custom __call__ method that implements the call, and state, such
as an environment containing non-local variables, is stored in the
instance of the callable class.
"""

from typing import Optional, List, Tuple, Union, Dict

from mypy.nodes import (
    ClassDef, FuncDef, OverloadedFuncDef, Decorator, Var, YieldFromExpr, AwaitExpr, YieldExpr,
    FuncItem, LambdaExpr, SymbolNode, ARG_NAMED, ARG_NAMED_OPT
)
from mypy.types import CallableType, get_proper_type

from mypyc.ir.ops import (
    BasicBlock, Value,  Register, Return, SetAttr, Integer, GetAttr, Branch, InitStatic,
    LoadAddress
)
from mypyc.ir.rtypes import object_rprimitive, RInstance, object_pointer_rprimitive
from mypyc.ir.func_ir import (
    FuncIR, FuncSignature, RuntimeArg, FuncDecl, FUNC_CLASSMETHOD, FUNC_STATICMETHOD, FUNC_NORMAL
)
from mypyc.ir.class_ir import ClassIR, NonExtClassInfo
from mypyc.primitives.generic_ops import py_setattr_op, next_raw_op, iter_op
from mypyc.primitives.misc_ops import check_stop_op, yield_from_except_op, coro_op, send_op
from mypyc.primitives.dict_ops import dict_set_item_op
from mypyc.common import SELF_NAME, LAMBDA_NAME, decorator_helper_name
from mypyc.sametype import is_same_method_signature
from mypyc.irbuild.util import concrete_arg_kind, is_constant
from mypyc.irbuild.context import FuncInfo, ImplicitClass
from mypyc.irbuild.targets import AssignmentTarget
from mypyc.irbuild.statement import transform_try_except
from mypyc.irbuild.builder import IRBuilder, SymbolTarget, gen_arg_defaults
from mypyc.irbuild.callable_class import (
    setup_callable_class, add_call_to_callable_class, add_get_to_callable_class,
    instantiate_callable_class
)
from mypyc.irbuild.generator import (
    gen_generator_func, setup_env_for_generator_class, create_switch_for_generator_class,
    add_raise_exception_blocks_to_generator_class, populate_switch_for_generator_class,
    add_methods_to_generator_class
)
from mypyc.irbuild.env_class import (
    setup_env_class, load_outer_envs, load_env_registers, finalize_env_class,
    setup_func_for_recursive_call
)


# Top-level transform functions


def transform_func_def(builder: IRBuilder, fdef: FuncDef) -> None:
    func_ir, func_reg = gen_func_item(builder, fdef, fdef.name, builder.mapper.fdef_to_sig(fdef))

    # If the function that was visited was a nested function, then either look it up in our
    # current environment or define it if it was not already defined.
    if func_reg:
        builder.assign(get_func_target(builder, fdef), func_reg, fdef.line)
    builder.functions.append(func_ir)


def transform_overloaded_func_def(builder: IRBuilder, o: OverloadedFuncDef) -> None:
    # Handle regular overload case
    assert o.impl
    builder.accept(o.impl)


def transform_decorator(builder: IRBuilder, dec: Decorator) -> None:
    func_ir, func_reg = gen_func_item(
        builder,
        dec.func,
        dec.func.name,
        builder.mapper.fdef_to_sig(dec.func)
    )

    if dec.func in builder.nested_fitems:
        assert func_reg is not None
        decorated_func = load_decorated_func(builder, dec.func, func_reg)
        builder.assign(get_func_target(builder, dec.func), decorated_func, dec.func.line)
        func_reg = decorated_func
    else:
        # Obtain the the function name in order to construct the name of the helper function.
        name = dec.func.fullname.split('.')[-1]
        helper_name = decorator_helper_name(name)

        # Load the callable object representing the non-decorated function, and decorate it.
        orig_func = builder.load_global_str(helper_name, dec.line)
        decorated_func = load_decorated_func(builder, dec.func, orig_func)

        # Set the callable object representing the decorated function as a global.
        builder.call_c(dict_set_item_op,
                    [builder.load_globals_dict(),
                    builder.load_str(dec.func.name), decorated_func],
                    decorated_func.line)

    builder.functions.append(func_ir)


def transform_method(builder: IRBuilder,
                     cdef: ClassDef,
                     non_ext: Optional[NonExtClassInfo],
                     fdef: FuncDef) -> None:
    if non_ext:
        handle_non_ext_method(builder, non_ext, cdef, fdef)
    else:
        handle_ext_method(builder, cdef, fdef)


def transform_lambda_expr(builder: IRBuilder, expr: LambdaExpr) -> Value:
    typ = get_proper_type(builder.types[expr])
    assert isinstance(typ, CallableType)

    runtime_args = []
    for arg, arg_type in zip(expr.arguments, typ.arg_types):
        arg.variable.type = arg_type
        runtime_args.append(
            RuntimeArg(arg.variable.name, builder.type_to_rtype(arg_type), arg.kind))
    ret_type = builder.type_to_rtype(typ.ret_type)

    fsig = FuncSignature(runtime_args, ret_type)

    fname = '{}{}'.format(LAMBDA_NAME, builder.lambda_counter)
    builder.lambda_counter += 1
    func_ir, func_reg = gen_func_item(builder, expr, fname, fsig)
    assert func_reg is not None

    builder.functions.append(func_ir)
    return func_reg


def transform_yield_expr(builder: IRBuilder, expr: YieldExpr) -> Value:
    if expr.expr:
        retval = builder.accept(expr.expr)
    else:
        retval = builder.builder.none()
    return emit_yield(builder, retval, expr.line)


def transform_yield_from_expr(builder: IRBuilder, o: YieldFromExpr) -> Value:
    return handle_yield_from_and_await(builder, o)


def transform_await_expr(builder: IRBuilder, o: AwaitExpr) -> Value:
    return handle_yield_from_and_await(builder, o)


# Internal functions


def gen_func_item(builder: IRBuilder,
                  fitem: FuncItem,
                  name: str,
                  sig: FuncSignature,
                  cdef: Optional[ClassDef] = None,
                  ) -> Tuple[FuncIR, Optional[Value]]:
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

    func_reg = None  # type: Optional[Value]

    # We treat lambdas as always being nested because we always generate
    # a class for lambdas, no matter where they are. (It would probably also
    # work to special case toplevel lambdas and generate a non-class function.)
    is_nested = fitem in builder.nested_fitems or isinstance(fitem, LambdaExpr)
    contains_nested = fitem in builder.encapsulating_funcs.keys()
    is_decorated = fitem in builder.fdefs_to_decorators
    in_non_ext = False
    class_name = None
    if cdef:
        ir = builder.mapper.type_to_ir[cdef.info]
        in_non_ext = not ir.is_ext_class
        class_name = cdef.name

    builder.enter(FuncInfo(fitem, name, class_name, gen_func_ns(builder),
                           is_nested, contains_nested, is_decorated, in_non_ext))

    # Functions that contain nested functions need an environment class to store variables that
    # are free in their nested functions. Generator functions need an environment class to
    # store a variable denoting the next instruction to be executed when the __next__ function
    # is called, along with all the variables inside the function itself.
    if builder.fn_info.contains_nested or builder.fn_info.is_generator:
        setup_env_class(builder)

    if builder.fn_info.is_nested or builder.fn_info.in_non_ext:
        setup_callable_class(builder)

    if builder.fn_info.is_generator:
        # Do a first-pass and generate a function that just returns a generator object.
        gen_generator_func(builder)
        args, _, blocks, ret_type, fn_info = builder.leave()
        func_ir, func_reg = gen_func_ir(builder, args, blocks, sig, fn_info, cdef)

        # Re-enter the FuncItem and visit the body of the function this time.
        builder.enter(fn_info)
        setup_env_for_generator_class(builder)
        load_outer_envs(builder, builder.fn_info.generator_class)
        if builder.fn_info.is_nested and isinstance(fitem, FuncDef):
            setup_func_for_recursive_call(builder, fitem, builder.fn_info.generator_class)
        create_switch_for_generator_class(builder)
        add_raise_exception_blocks_to_generator_class(builder, fitem.line)
    else:
        load_env_registers(builder)
        gen_arg_defaults(builder)

    if builder.fn_info.contains_nested and not builder.fn_info.is_generator:
        finalize_env_class(builder)

    builder.ret_types[-1] = sig.ret_type

    # Add all variables and functions that are declared/defined within this
    # function and are referenced in functions nested within this one to this
    # function's environment class so the nested functions can reference
    # them even if they are declared after the nested function's definition.
    # Note that this is done before visiting the body of this function.

    env_for_func = builder.fn_info  # type: Union[FuncInfo, ImplicitClass]
    if builder.fn_info.is_generator:
        env_for_func = builder.fn_info.generator_class
    elif builder.fn_info.is_nested or builder.fn_info.in_non_ext:
        env_for_func = builder.fn_info.callable_class

    if builder.fn_info.fitem in builder.free_variables:
        # Sort the variables to keep things deterministic
        for var in sorted(builder.free_variables[builder.fn_info.fitem],
                          key=lambda x: x.name):
            if isinstance(var, Var):
                rtype = builder.type_to_rtype(var.type)
                builder.add_var_to_env_class(var, rtype, env_for_func, reassign=False)

    if builder.fn_info.fitem in builder.encapsulating_funcs:
        for nested_fn in builder.encapsulating_funcs[builder.fn_info.fitem]:
            if isinstance(nested_fn, FuncDef):
                # The return type is 'object' instead of an RInstance of the
                # callable class because differently defined functions with
                # the same name and signature across conditional blocks
                # will generate different callable classes, so the callable
                # class that gets instantiated must be generic.
                builder.add_var_to_env_class(
                    nested_fn, object_rprimitive, env_for_func, reassign=False
                )

    builder.accept(fitem.body)
    builder.maybe_add_implicit_return()

    if builder.fn_info.is_generator:
        populate_switch_for_generator_class(builder)

    # Hang on to the local symbol table for a while, since we use it
    # to calculate argument defaults below.
    symtable = builder.symtables[-1]

    args, _, blocks, ret_type, fn_info = builder.leave()

    if fn_info.is_generator:
        add_methods_to_generator_class(
            builder, fn_info, sig, args, blocks, fitem.is_coroutine)
    else:
        func_ir, func_reg = gen_func_ir(builder, args, blocks, sig, fn_info, cdef)

    # Evaluate argument defaults in the surrounding scope, since we
    # calculate them *once* when the function definition is evaluated.
    calculate_arg_defaults(builder, fn_info, func_reg, symtable)

    return (func_ir, func_reg)


def gen_func_ir(builder: IRBuilder,
                args: List[Register],
                blocks: List[BasicBlock],
                sig: FuncSignature,
                fn_info: FuncInfo,
                cdef: Optional[ClassDef]) -> Tuple[FuncIR, Optional[Value]]:
    """Generate the FuncIR for a function.

    This takes the basic blocks and function info of a particular
    function and returns the IR. If the function is nested,
    also returns the register containing the instance of the
    corresponding callable class.
    """
    func_reg = None  # type: Optional[Value]
    if fn_info.is_nested or fn_info.in_non_ext:
        func_ir = add_call_to_callable_class(builder, args, blocks, sig, fn_info)
        add_get_to_callable_class(builder, fn_info)
        func_reg = instantiate_callable_class(builder, fn_info)
    else:
        assert isinstance(fn_info.fitem, FuncDef)
        func_decl = builder.mapper.func_to_decl[fn_info.fitem]
        if fn_info.is_decorated:
            class_name = None if cdef is None else cdef.name
            func_decl = FuncDecl(fn_info.name, class_name, builder.module_name, sig,
                                 func_decl.kind,
                                 func_decl.is_prop_getter, func_decl.is_prop_setter)
            func_ir = FuncIR(func_decl, args, blocks, fn_info.fitem.line,
                             traceback_name=fn_info.fitem.name)
        else:
            func_ir = FuncIR(func_decl, args, blocks,
                             fn_info.fitem.line, traceback_name=fn_info.fitem.name)
    return (func_ir, func_reg)


def handle_ext_method(builder: IRBuilder, cdef: ClassDef, fdef: FuncDef) -> None:
    # Perform the function of visit_method for methods inside extension classes.
    name = fdef.name
    class_ir = builder.mapper.type_to_ir[cdef.info]
    func_ir, func_reg = gen_func_item(builder, fdef, name, builder.mapper.fdef_to_sig(fdef), cdef)
    builder.functions.append(func_ir)

    if is_decorated(builder, fdef):
        # Obtain the the function name in order to construct the name of the helper function.
        _, _, name = fdef.fullname.rpartition('.')
        helper_name = decorator_helper_name(name)
        # Read the PyTypeObject representing the class, get the callable object
        # representing the non-decorated method
        typ = builder.load_native_type_object(cdef.fullname)
        orig_func = builder.py_get_attr(typ, helper_name, fdef.line)

        # Decorate the non-decorated method
        decorated_func = load_decorated_func(builder, fdef, orig_func)

        # Set the callable object representing the decorated method as an attribute of the
        # extension class.
        builder.call_c(py_setattr_op,
                    [
                        typ,
                        builder.load_str(name),
                        decorated_func
                    ],
                    fdef.line)

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
        if (name in base.method_decls and name != '__init__'
                and not is_same_method_signature(class_ir.method_decls[name].sig,
                                                 base.method_decls[name].sig)):

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
        class_ir.glue_methods[(class_ir, name)] = f
        builder.functions.append(f)


def handle_non_ext_method(
        builder: IRBuilder, non_ext: NonExtClassInfo, cdef: ClassDef, fdef: FuncDef) -> None:
    # Perform the function of visit_method for methods inside non-extension classes.
    name = fdef.name
    func_ir, func_reg = gen_func_item(builder, fdef, name, builder.mapper.fdef_to_sig(fdef), cdef)
    assert func_reg is not None
    builder.functions.append(func_ir)

    if is_decorated(builder, fdef):
        # The undecorated method is a generated callable class
        orig_func = func_reg
        func_reg = load_decorated_func(builder, fdef, orig_func)

    # TODO: Support property setters in non-extension classes
    if fdef.is_property:
        prop = builder.load_module_attr_by_fullname('builtins.property', fdef.line)
        func_reg = builder.py_call(prop, [func_reg], fdef.line)

    elif builder.mapper.func_to_decl[fdef].kind == FUNC_CLASSMETHOD:
        cls_meth = builder.load_module_attr_by_fullname('builtins.classmethod', fdef.line)
        func_reg = builder.py_call(cls_meth, [func_reg], fdef.line)

    elif builder.mapper.func_to_decl[fdef].kind == FUNC_STATICMETHOD:
        stat_meth = builder.load_module_attr_by_fullname(
            'builtins.staticmethod', fdef.line
        )
        func_reg = builder.py_call(stat_meth, [func_reg], fdef.line)

    builder.add_to_non_ext_dict(non_ext, name, func_reg, fdef.line)


def calculate_arg_defaults(builder: IRBuilder,
                           fn_info: FuncInfo,
                           func_reg: Optional[Value],
                           symtable: Dict[SymbolNode, SymbolTarget]) -> None:
    """Calculate default argument values and store them.

    They are stored in statics for top level functions and in
    the function objects for nested functions (while constants are
    still stored computed on demand).
    """
    fitem = fn_info.fitem
    for arg in fitem.arguments:
        # Constant values don't get stored but just recomputed
        if arg.initializer and not is_constant(arg.initializer):
            value = builder.coerce(
                builder.accept(arg.initializer),
                symtable[arg.variable].type,
                arg.line
            )
            if not fn_info.is_nested:
                name = fitem.fullname + '.' + arg.variable.name
                builder.add(InitStatic(value, name, builder.module_name))
            else:
                assert func_reg is not None
                builder.add(SetAttr(func_reg, arg.variable.name, value, arg.line))


def gen_func_ns(builder: IRBuilder) -> str:
    """Generate a namespace for a nested function using its outer function names."""
    return '_'.join(info.name + ('' if not info.class_name else '_' + info.class_name)
                    for info in builder.fn_infos
                    if info.name and info.name != '<top level>')


def emit_yield(builder: IRBuilder, val: Value, line: int) -> Value:
    retval = builder.coerce(val, builder.ret_types[-1], line)

    cls = builder.fn_info.generator_class
    # Create a new block for the instructions immediately following the yield expression, and
    # set the next label so that the next time '__next__' is called on the generator object,
    # the function continues at the new block.
    next_block = BasicBlock()
    next_label = len(cls.continuation_blocks)
    cls.continuation_blocks.append(next_block)
    builder.assign(cls.next_label_target, Integer(next_label), line)
    builder.add(Return(retval))
    builder.activate_block(next_block)

    add_raise_exception_blocks_to_generator_class(builder, line)

    assert cls.send_arg_reg is not None
    return cls.send_arg_reg


def handle_yield_from_and_await(builder: IRBuilder, o: Union[YieldFromExpr, AwaitExpr]) -> Value:
    # This is basically an implementation of the code in PEP 380.

    # TODO: do we want to use the right types here?
    result = Register(object_rprimitive)
    to_yield_reg = Register(object_rprimitive)
    received_reg = Register(object_rprimitive)

    if isinstance(o, YieldFromExpr):
        iter_val = builder.call_c(iter_op, [builder.accept(o.expr)], o.line)
    else:
        iter_val = builder.call_c(coro_op, [builder.accept(o.expr)], o.line)

    iter_reg = builder.maybe_spill_assignable(iter_val)

    stop_block, main_block, done_block = BasicBlock(), BasicBlock(), BasicBlock()
    _y_init = builder.call_c(next_raw_op, [builder.read(iter_reg)], o.line)
    builder.add(Branch(_y_init, stop_block, main_block, Branch.IS_ERROR))

    # Try extracting a return value from a StopIteration and return it.
    # If it wasn't, this reraises the exception.
    builder.activate_block(stop_block)
    builder.assign(result, builder.call_c(check_stop_op, [], o.line), o.line)
    builder.goto(done_block)

    builder.activate_block(main_block)
    builder.assign(to_yield_reg, _y_init, o.line)

    # OK Now the main loop!
    loop_block = BasicBlock()
    builder.goto_and_activate(loop_block)

    def try_body() -> None:
        builder.assign(
            received_reg, emit_yield(builder, builder.read(to_yield_reg), o.line), o.line
        )

    def except_body() -> None:
        # The body of the except is all implemented in a C function to
        # reduce how much code we need to generate. It returns a value
        # indicating whether to break or yield (or raise an exception).
        val = Register(object_rprimitive)
        val_address = builder.add(LoadAddress(object_pointer_rprimitive, val))
        to_stop = builder.call_c(yield_from_except_op,
                                 [builder.read(iter_reg), val_address], o.line)

        ok, stop = BasicBlock(), BasicBlock()
        builder.add(Branch(to_stop, stop, ok, Branch.BOOL))

        # The exception got swallowed. Continue, yielding the returned value
        builder.activate_block(ok)
        builder.assign(to_yield_reg, val, o.line)
        builder.nonlocal_control[-1].gen_continue(builder, o.line)

        # The exception was a StopIteration. Stop iterating.
        builder.activate_block(stop)
        builder.assign(result, val, o.line)
        builder.nonlocal_control[-1].gen_break(builder, o.line)

    def else_body() -> None:
        # Do a next() or a .send(). It will return NULL on exception
        # but it won't automatically propagate.
        _y = builder.call_c(
            send_op, [builder.read(iter_reg), builder.read(received_reg)], o.line
        )
        ok, stop = BasicBlock(), BasicBlock()
        builder.add(Branch(_y, stop, ok, Branch.IS_ERROR))

        # Everything's fine. Yield it.
        builder.activate_block(ok)
        builder.assign(to_yield_reg, _y, o.line)
        builder.nonlocal_control[-1].gen_continue(builder, o.line)

        # Try extracting a return value from a StopIteration and return it.
        # If it wasn't, this rereaises the exception.
        builder.activate_block(stop)
        builder.assign(result, builder.call_c(check_stop_op, [], o.line), o.line)
        builder.nonlocal_control[-1].gen_break(builder, o.line)

    builder.push_loop_stack(loop_block, done_block)
    transform_try_except(
        builder, try_body, [(None, None, except_body)], else_body, o.line
    )
    builder.pop_loop_stack()

    builder.goto_and_activate(done_block)
    return builder.read(result)


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
        assert isinstance(decorator, Value)
        func_reg = builder.py_call(decorator, [func_reg], func_reg.line)
    return func_reg


def is_decorated(builder: IRBuilder, fdef: FuncDef) -> bool:
    return fdef in builder.fdefs_to_decorators


def gen_glue(builder: IRBuilder, sig: FuncSignature, target: FuncIR,
             cls: ClassIR, base: ClassIR, fdef: FuncItem,
             *,
             do_py_ops: bool = False
             ) -> FuncIR:
    """Generate glue methods that mediate between different method types in subclasses.

    Works on both properties and methods. See gen_glue_methods below
    for more details.

    If do_py_ops is True, then the glue methods should use generic
    C API operations instead of direct calls, to enable generating
    "shadow" glue methods that work with interpreted subclasses.
    """
    if fdef.is_property:
        return gen_glue_property(builder, sig, target, cls, base, fdef.line, do_py_ops)
    else:
        return gen_glue_method(builder, sig, target, cls, base, fdef.line, do_py_ops)


def gen_glue_method(builder: IRBuilder, sig: FuncSignature, target: FuncIR,
                    cls: ClassIR, base: ClassIR, line: int,
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
    builder.enter()
    builder.ret_types[-1] = sig.ret_type

    rt_args = list(sig.args)
    if target.decl.kind == FUNC_NORMAL:
        rt_args[0] = RuntimeArg(sig.args[0].name, RInstance(cls))

    # The environment operates on Vars, so we make some up
    fake_vars = [(Var(arg.name), arg.type) for arg in rt_args]
    args = [builder.read(builder.add_local_reg(var, type, is_arg=True), line)
            for var, type in fake_vars]
    arg_names = [arg.name if arg.kind in (ARG_NAMED, ARG_NAMED_OPT) else None
                 for arg in rt_args]
    arg_kinds = [concrete_arg_kind(arg.kind) for arg in rt_args]

    if do_pycall:
        retval = builder.builder.py_method_call(
            args[0], target.name, args[1:], line, arg_kinds[1:], arg_names[1:])
    else:
        retval = builder.builder.call(target.decl, args, arg_kinds, arg_names, line)
    retval = builder.coerce(retval, sig.ret_type, line)
    builder.add(Return(retval))

    arg_regs, _, blocks, ret_type, _ = builder.leave()
    return FuncIR(
        FuncDecl(target.name + '__' + base.name + '_glue',
                 cls.name, builder.module_name,
                 FuncSignature(rt_args, ret_type),
                 target.decl.kind),
        arg_regs, blocks)


def gen_glue_property(builder: IRBuilder,
                      sig: FuncSignature,
                      target: FuncIR,
                      cls: ClassIR,
                      base: ClassIR,
                      line: int,
                      do_pygetattr: bool) -> FuncIR:
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
        FuncDecl(target.name + '__' + base.name + '_glue',
                 cls.name, builder.module_name, FuncSignature([rt_arg], return_type)),
        args, blocks)


def get_func_target(builder: IRBuilder, fdef: FuncDef) -> AssignmentTarget:
    """Given a FuncDef, return the target for the instance of its callable class.

    If the function was not already defined somewhere, then define it
    and add it to the current environment.
    """
    if fdef.original_def:
        # Get the target associated with the previously defined FuncDef.
        return builder.lookup(fdef.original_def)

    if builder.fn_info.is_generator or builder.fn_info.contains_nested:
        return builder.lookup(fdef)

    return builder.add_local_reg(fdef, object_rprimitive)
