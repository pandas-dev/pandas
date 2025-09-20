"""Generate IR for generator functions.

A generator function is represented by a class that implements the
generator protocol and keeps track of the generator state, including
local variables.

The top-level logic for dealing with generator functions is in
mypyc.irbuild.function.
"""

from __future__ import annotations

from typing import Callable

from mypy.nodes import ARG_OPT, FuncDef, Var
from mypyc.common import ENV_ATTR_NAME, NEXT_LABEL_ATTR_NAME, SELF_NAME
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.func_ir import FuncDecl, FuncIR, FuncSignature, RuntimeArg
from mypyc.ir.ops import (
    NO_TRACEBACK_LINE_NO,
    BasicBlock,
    Branch,
    Call,
    Goto,
    Integer,
    MethodCall,
    RaiseStandardError,
    Register,
    Return,
    SetAttr,
    TupleSet,
    Unreachable,
    Value,
)
from mypyc.ir.rtypes import RInstance, int32_rprimitive, object_rprimitive
from mypyc.irbuild.builder import IRBuilder, calculate_arg_defaults, gen_arg_defaults
from mypyc.irbuild.context import FuncInfo, GeneratorClass
from mypyc.irbuild.env_class import (
    add_args_to_env,
    add_vars_to_env,
    finalize_env_class,
    load_env_registers,
    load_outer_env,
    load_outer_envs,
    setup_func_for_recursive_call,
)
from mypyc.irbuild.nonlocalcontrol import ExceptNonlocalControl
from mypyc.primitives.exc_ops import (
    error_catch_op,
    exc_matches_op,
    raise_exception_with_tb_op,
    reraise_exception_op,
    restore_exc_info_op,
)


def gen_generator_func(
    builder: IRBuilder,
    gen_func_ir: Callable[
        [list[Register], list[BasicBlock], FuncInfo], tuple[FuncIR, Value | None]
    ],
) -> tuple[FuncIR, Value | None]:
    """Generate IR for generator function that returns generator object."""
    setup_generator_class(builder)
    load_env_registers(builder)
    gen_arg_defaults(builder)
    if builder.fn_info.can_merge_generator_and_env_classes():
        gen = instantiate_generator_class(builder)
        builder.fn_info._curr_env_reg = gen
        finalize_env_class(builder)
    else:
        finalize_env_class(builder)
        gen = instantiate_generator_class(builder)
    builder.add(Return(gen))

    args, _, blocks, ret_type, fn_info = builder.leave()
    func_ir, func_reg = gen_func_ir(args, blocks, fn_info)
    return func_ir, func_reg


def gen_generator_func_body(
    builder: IRBuilder, fn_info: FuncInfo, sig: FuncSignature, func_reg: Value | None
) -> None:
    """Generate IR based on the body of a generator function.

    Add "__next__", "__iter__" and other generator methods to the generator
    class that implements the function (each function gets a separate class).

    Return the symbol table for the body.
    """
    builder.enter(fn_info, ret_type=sig.ret_type)
    setup_env_for_generator_class(builder)

    load_outer_envs(builder, builder.fn_info.generator_class)
    top_level = builder.top_level_fn_info()
    fitem = fn_info.fitem
    if (
        builder.fn_info.is_nested
        and isinstance(fitem, FuncDef)
        and top_level
        and top_level.add_nested_funcs_to_env
    ):
        setup_func_for_recursive_call(builder, fitem, builder.fn_info.generator_class)
    create_switch_for_generator_class(builder)
    add_raise_exception_blocks_to_generator_class(builder, fitem.line)

    add_vars_to_env(builder)

    builder.accept(fitem.body)
    builder.maybe_add_implicit_return()

    populate_switch_for_generator_class(builder)

    # Hang on to the local symbol table, since the caller will use it
    # to calculate argument defaults.
    symtable = builder.symtables[-1]

    args, _, blocks, ret_type, fn_info = builder.leave()

    add_methods_to_generator_class(builder, fn_info, sig, args, blocks, fitem.is_coroutine)

    # Evaluate argument defaults in the surrounding scope, since we
    # calculate them *once* when the function definition is evaluated.
    calculate_arg_defaults(builder, fn_info, func_reg, symtable)


def instantiate_generator_class(builder: IRBuilder) -> Value:
    fitem = builder.fn_info.fitem
    generator_reg = builder.add(Call(builder.fn_info.generator_class.ir.ctor, [], fitem.line))

    if builder.fn_info.can_merge_generator_and_env_classes():
        # Set the generator instance to the initial state (zero).
        zero = Integer(0)
        builder.add(SetAttr(generator_reg, NEXT_LABEL_ATTR_NAME, zero, fitem.line))
    else:
        # Get the current environment register. If the current function is nested, then the
        # generator class gets instantiated from the callable class' '__call__' method, and hence
        # we use the callable class' environment register. Otherwise, we use the original
        # function's environment register.
        if builder.fn_info.is_nested:
            curr_env_reg = builder.fn_info.callable_class.curr_env_reg
        else:
            curr_env_reg = builder.fn_info.curr_env_reg

        # Set the generator class' environment attribute to point at the environment class
        # defined in the current scope.
        builder.add(SetAttr(generator_reg, ENV_ATTR_NAME, curr_env_reg, fitem.line))

        # Set the generator instance's environment to the initial state (zero).
        zero = Integer(0)
        builder.add(SetAttr(curr_env_reg, NEXT_LABEL_ATTR_NAME, zero, fitem.line))
    return generator_reg


def setup_generator_class(builder: IRBuilder) -> ClassIR:
    name = f"{builder.fn_info.namespaced_name()}_gen"

    generator_class_ir = ClassIR(name, builder.module_name, is_generated=True, is_final_class=True)
    if builder.fn_info.can_merge_generator_and_env_classes():
        builder.fn_info.env_class = generator_class_ir
    else:
        generator_class_ir.attributes[ENV_ATTR_NAME] = RInstance(builder.fn_info.env_class)
    generator_class_ir.mro = [generator_class_ir]

    builder.classes.append(generator_class_ir)
    builder.fn_info.generator_class = GeneratorClass(generator_class_ir)
    return generator_class_ir


def create_switch_for_generator_class(builder: IRBuilder) -> None:
    builder.add(Goto(builder.fn_info.generator_class.switch_block))
    block = BasicBlock()
    builder.fn_info.generator_class.continuation_blocks.append(block)
    builder.activate_block(block)


def populate_switch_for_generator_class(builder: IRBuilder) -> None:
    cls = builder.fn_info.generator_class
    line = builder.fn_info.fitem.line

    builder.activate_block(cls.switch_block)
    for label, true_block in enumerate(cls.continuation_blocks):
        false_block = BasicBlock()
        comparison = builder.binary_op(cls.next_label_reg, Integer(label), "==", line)
        builder.add_bool_branch(comparison, true_block, false_block)
        builder.activate_block(false_block)

    builder.add(RaiseStandardError(RaiseStandardError.STOP_ITERATION, None, line))
    builder.add(Unreachable())


def add_raise_exception_blocks_to_generator_class(builder: IRBuilder, line: int) -> None:
    """Add error handling blocks to a generator class.

    Generates blocks to check if error flags are set while calling the
    helper method for generator functions, and raises an exception if
    those flags are set.
    """
    cls = builder.fn_info.generator_class
    assert cls.exc_regs is not None
    exc_type, exc_val, exc_tb = cls.exc_regs

    # Check to see if an exception was raised.
    error_block = BasicBlock()
    ok_block = BasicBlock()
    comparison = builder.translate_is_op(exc_type, builder.none_object(), "is not", line)
    builder.add_bool_branch(comparison, error_block, ok_block)

    builder.activate_block(error_block)
    builder.call_c(raise_exception_with_tb_op, [exc_type, exc_val, exc_tb], line)
    builder.add(Unreachable())
    builder.goto_and_activate(ok_block)


def add_methods_to_generator_class(
    builder: IRBuilder,
    fn_info: FuncInfo,
    sig: FuncSignature,
    arg_regs: list[Register],
    blocks: list[BasicBlock],
    is_coroutine: bool,
) -> None:
    helper_fn_decl = add_helper_to_generator_class(builder, arg_regs, blocks, sig, fn_info)
    add_next_to_generator_class(builder, fn_info, helper_fn_decl, sig)
    add_send_to_generator_class(builder, fn_info, helper_fn_decl, sig)
    add_iter_to_generator_class(builder, fn_info)
    add_throw_to_generator_class(builder, fn_info, helper_fn_decl, sig)
    add_close_to_generator_class(builder, fn_info)
    if is_coroutine:
        add_await_to_generator_class(builder, fn_info)


def add_helper_to_generator_class(
    builder: IRBuilder,
    arg_regs: list[Register],
    blocks: list[BasicBlock],
    sig: FuncSignature,
    fn_info: FuncInfo,
) -> FuncDecl:
    """Generates a helper method for a generator class, called by '__next__' and 'throw'."""
    sig = FuncSignature(
        (
            RuntimeArg(SELF_NAME, object_rprimitive),
            RuntimeArg("type", object_rprimitive),
            RuntimeArg("value", object_rprimitive),
            RuntimeArg("traceback", object_rprimitive),
            RuntimeArg("arg", object_rprimitive),
        ),
        sig.ret_type,
    )
    helper_fn_decl = FuncDecl(
        "__mypyc_generator_helper__",
        fn_info.generator_class.ir.name,
        builder.module_name,
        sig,
        internal=True,
    )
    helper_fn_ir = FuncIR(
        helper_fn_decl, arg_regs, blocks, fn_info.fitem.line, traceback_name=fn_info.fitem.name
    )
    fn_info.generator_class.ir.methods["__mypyc_generator_helper__"] = helper_fn_ir
    builder.functions.append(helper_fn_ir)
    fn_info.env_class.env_user_function = helper_fn_ir

    return helper_fn_decl


def add_iter_to_generator_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generates the '__iter__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, "__iter__", object_rprimitive, fn_info):
        builder.add(Return(builder.self()))


def add_next_to_generator_class(
    builder: IRBuilder, fn_info: FuncInfo, fn_decl: FuncDecl, sig: FuncSignature
) -> None:
    """Generates the '__next__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, "__next__", object_rprimitive, fn_info):
        none_reg = builder.none_object()
        # Call the helper function with error flags set to Py_None, and return that result.
        result = builder.add(
            Call(
                fn_decl,
                [builder.self(), none_reg, none_reg, none_reg, none_reg],
                fn_info.fitem.line,
            )
        )
        builder.add(Return(result))


def add_send_to_generator_class(
    builder: IRBuilder, fn_info: FuncInfo, fn_decl: FuncDecl, sig: FuncSignature
) -> None:
    """Generates the 'send' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, "send", object_rprimitive, fn_info):
        arg = builder.add_argument("arg", object_rprimitive)
        none_reg = builder.none_object()
        # Call the helper function with error flags set to Py_None, and return that result.
        result = builder.add(
            Call(
                fn_decl,
                [builder.self(), none_reg, none_reg, none_reg, builder.read(arg)],
                fn_info.fitem.line,
            )
        )
        builder.add(Return(result))


def add_throw_to_generator_class(
    builder: IRBuilder, fn_info: FuncInfo, fn_decl: FuncDecl, sig: FuncSignature
) -> None:
    """Generates the 'throw' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, "throw", object_rprimitive, fn_info):
        typ = builder.add_argument("type", object_rprimitive)
        val = builder.add_argument("value", object_rprimitive, ARG_OPT)
        tb = builder.add_argument("traceback", object_rprimitive, ARG_OPT)

        # Because the value and traceback arguments are optional and hence
        # can be NULL if not passed in, we have to assign them Py_None if
        # they are not passed in.
        none_reg = builder.none_object()
        builder.assign_if_null(val, lambda: none_reg, builder.fn_info.fitem.line)
        builder.assign_if_null(tb, lambda: none_reg, builder.fn_info.fitem.line)

        # Call the helper function using the arguments passed in, and return that result.
        result = builder.add(
            Call(
                fn_decl,
                [builder.self(), builder.read(typ), builder.read(val), builder.read(tb), none_reg],
                fn_info.fitem.line,
            )
        )
        builder.add(Return(result))


def add_close_to_generator_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generates the '__close__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, "close", object_rprimitive, fn_info):
        except_block, else_block = BasicBlock(), BasicBlock()
        builder.builder.push_error_handler(except_block)
        builder.goto_and_activate(BasicBlock())
        generator_exit = builder.load_module_attr_by_fullname(
            "builtins.GeneratorExit", fn_info.fitem.line
        )
        builder.add(
            MethodCall(
                builder.self(),
                "throw",
                [generator_exit, builder.none_object(), builder.none_object()],
            )
        )
        builder.goto(else_block)
        builder.builder.pop_error_handler()

        builder.activate_block(except_block)
        old_exc = builder.call_c(error_catch_op, [], fn_info.fitem.line)
        builder.nonlocal_control.append(
            ExceptNonlocalControl(builder.nonlocal_control[-1], old_exc)
        )
        stop_iteration = builder.load_module_attr_by_fullname(
            "builtins.StopIteration", fn_info.fitem.line
        )
        exceptions = builder.add(TupleSet([generator_exit, stop_iteration], fn_info.fitem.line))
        matches = builder.call_c(exc_matches_op, [exceptions], fn_info.fitem.line)

        match_block, non_match_block = BasicBlock(), BasicBlock()
        builder.add(Branch(matches, match_block, non_match_block, Branch.BOOL))

        builder.activate_block(match_block)
        builder.call_c(restore_exc_info_op, [builder.read(old_exc)], fn_info.fitem.line)
        builder.add(Return(builder.none_object()))

        builder.activate_block(non_match_block)
        builder.call_c(reraise_exception_op, [], NO_TRACEBACK_LINE_NO)
        builder.add(Unreachable())

        builder.nonlocal_control.pop()

        builder.activate_block(else_block)
        builder.add(
            RaiseStandardError(
                RaiseStandardError.RUNTIME_ERROR,
                "generator ignored GeneratorExit",
                fn_info.fitem.line,
            )
        )
        builder.add(Unreachable())


def add_await_to_generator_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generates the '__await__' method for a generator class."""
    with builder.enter_method(fn_info.generator_class.ir, "__await__", object_rprimitive, fn_info):
        builder.add(Return(builder.self()))


def setup_env_for_generator_class(builder: IRBuilder) -> None:
    """Populates the environment for a generator class."""
    fitem = builder.fn_info.fitem
    cls = builder.fn_info.generator_class
    self_target = builder.add_self_to_env(cls.ir)

    # Add the type, value, and traceback variables to the environment.
    exc_type = builder.add_local(Var("type"), object_rprimitive, is_arg=True)
    exc_val = builder.add_local(Var("value"), object_rprimitive, is_arg=True)
    exc_tb = builder.add_local(Var("traceback"), object_rprimitive, is_arg=True)
    # TODO: Use the right type here instead of object?
    exc_arg = builder.add_local(Var("arg"), object_rprimitive, is_arg=True)

    cls.exc_regs = (exc_type, exc_val, exc_tb)
    cls.send_arg_reg = exc_arg

    cls.self_reg = builder.read(self_target, fitem.line)
    if builder.fn_info.can_merge_generator_and_env_classes():
        cls.curr_env_reg = cls.self_reg
    else:
        cls.curr_env_reg = load_outer_env(builder, cls.self_reg, builder.symtables[-1])

    # Define a variable representing the label to go to the next time
    # the '__next__' function of the generator is called, and add it
    # as an attribute to the environment class.
    cls.next_label_target = builder.add_var_to_env_class(
        Var(NEXT_LABEL_ATTR_NAME), int32_rprimitive, cls, reassign=False, always_defined=True
    )

    # Add arguments from the original generator function to the
    # environment of the generator class.
    add_args_to_env(builder, local=False, base=cls, reassign=False)

    # Set the next label register for the generator class.
    cls.next_label_reg = builder.read(cls.next_label_target, fitem.line)
