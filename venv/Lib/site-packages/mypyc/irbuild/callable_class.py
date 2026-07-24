"""Generate a class that represents a nested function.

The class defines __call__ for calling the function and allows access to
non-local variables defined in outer scopes.
"""

from __future__ import annotations

from mypyc.common import CPYFUNCTION_NAME, ENV_ATTR_NAME, PROPSET_PREFIX, SELF_NAME
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.func_ir import FuncDecl, FuncIR, FuncSignature, RuntimeArg
from mypyc.ir.ops import BasicBlock, Call, GetAttr, Integer, Register, Return, SetAttr, Value
from mypyc.ir.rtypes import RInstance, c_pointer_rprimitive, int_rprimitive, object_rprimitive
from mypyc.irbuild.builder import IRBuilder
from mypyc.irbuild.context import FuncInfo, ImplicitClass
from mypyc.primitives.misc_ops import (
    cpyfunction_get_annotations,
    cpyfunction_get_code,
    cpyfunction_get_defaults,
    cpyfunction_get_kwdefaults,
    cpyfunction_get_name,
    cpyfunction_set_annotations,
    cpyfunction_set_name,
    method_new_op,
)


def setup_callable_class(builder: IRBuilder) -> None:
    """Generate an (incomplete) callable class representing a function.

    This can be a nested function or a function within a non-extension
    class.  Also set up the 'self' variable for that class.

    This takes the most recently visited function and returns a
    ClassIR to represent that function. Each callable class contains
    an environment attribute which points to another ClassIR
    representing the environment class where some of its variables can
    be accessed.

    Note that some methods, such as '__call__', are not yet
    created here. Use additional functions, such as
    add_call_to_callable_class(), to add them.

    Return a newly constructed ClassIR representing the callable
    class for the nested function.
    """
    # Check to see that the name has not already been taken. If so,
    # rename the class. We allow multiple uses of the same function
    # name because this is valid in if-else blocks. Example:
    #
    #     if True:
    #         def foo():          ---->    foo_obj()
    #             return True
    #     else:
    #         def foo():          ---->    foo_obj_0()
    #             return False
    name = base_name = f"{builder.fn_info.namespaced_name()}_obj"
    count = 0
    while name in builder.callable_class_names:
        name = base_name + "_" + str(count)
        count += 1
    builder.callable_class_names.add(name)

    # Define the actual callable class ClassIR, and set its
    # environment to point at the previously defined environment
    # class.
    callable_class_ir = ClassIR(name, builder.module_name, is_generated=True, is_final_class=True)
    callable_class_ir.reuse_freed_instance = True

    # The functools @wraps decorator attempts to call setattr on
    # nested functions, so we create a dict for these nested
    # functions.
    # https://github.com/python/cpython/blob/3.7/Lib/functools.py#L58
    if builder.fn_info.is_nested:
        callable_class_ir.has_dict = True

    # If the enclosing class doesn't contain nested (which will happen if
    # this is a toplevel lambda), don't set up an environment.
    if builder.fn_infos[-2].contains_nested:
        callable_class_ir.attributes[ENV_ATTR_NAME] = RInstance(builder.fn_infos[-2].env_class)
    callable_class_ir.mro = [callable_class_ir]
    builder.fn_info.callable_class = ImplicitClass(callable_class_ir)
    builder.classes.append(callable_class_ir)

    # Add a 'self' variable to the environment of the callable class,
    # and store that variable in a register to be accessed later.
    self_target = builder.add_self_to_env(callable_class_ir)
    builder.fn_info.callable_class.self_reg = builder.read(self_target, builder.fn_info.fitem.line)

    if not builder.fn_info.in_non_ext and builder.fn_info.is_coroutine:
        add_coroutine_properties(builder, callable_class_ir, builder.fn_info.name)


def add_coroutine_properties(
    builder: IRBuilder, callable_class_ir: ClassIR, coroutine_name: str
) -> None:
    """Adds properties to the class to make it look like a regular python function.
    Needed to make introspection functions like inspect.iscoroutinefunction work.
    """
    callable_class_ir.coroutine_name = coroutine_name
    callable_class_ir.attributes[CPYFUNCTION_NAME] = object_rprimitive

    properties = {
        "__name__": cpyfunction_get_name,
        "__code__": cpyfunction_get_code,
        "__annotations__": cpyfunction_get_annotations,
        "__defaults__": cpyfunction_get_defaults,
        "__kwdefaults__": cpyfunction_get_kwdefaults,
    }

    writable_props = {
        "__name__": cpyfunction_set_name,
        "__annotations__": cpyfunction_set_annotations,
    }

    line = builder.fn_info.fitem.line

    def get_func_wrapper() -> Value:
        return builder.add(GetAttr(builder.self(), CPYFUNCTION_NAME, line))

    for name, primitive in properties.items():
        with builder.enter_method(callable_class_ir, name, object_rprimitive, internal=True):
            func = get_func_wrapper()
            val = builder.primitive_op(primitive, [func, Integer(0, c_pointer_rprimitive)], line)
            builder.add(Return(val))

    for name, primitive in writable_props.items():
        with builder.enter_method(
            callable_class_ir, f"{PROPSET_PREFIX}{name}", int_rprimitive, internal=True
        ):
            value = builder.add_argument("value", object_rprimitive)
            func = get_func_wrapper()
            rv = builder.primitive_op(
                primitive, [func, value, Integer(0, c_pointer_rprimitive)], line
            )
            builder.add(Return(rv))

    for name in properties:
        getter = callable_class_ir.get_method(name)
        assert getter
        setter = callable_class_ir.get_method(f"{PROPSET_PREFIX}{name}")
        callable_class_ir.properties[name] = (getter, setter)


def add_call_to_callable_class(
    builder: IRBuilder,
    args: list[Register],
    blocks: list[BasicBlock],
    sig: FuncSignature,
    fn_info: FuncInfo,
) -> FuncIR:
    """Generate a '__call__' method for a callable class representing a nested function.

    This takes the blocks and signature associated with a function
    definition and uses those to build the '__call__' method of a
    given callable class, used to represent that function.
    """
    # Since we create a method, we also add a 'self' parameter.
    nargs = len(sig.args) - sig.num_bitmap_args
    sig = FuncSignature(
        (RuntimeArg(SELF_NAME, object_rprimitive),) + sig.args[:nargs], sig.ret_type
    )
    call_fn_decl = FuncDecl("__call__", fn_info.callable_class.ir.name, builder.module_name, sig)
    call_fn_ir = FuncIR(
        call_fn_decl, args, blocks, fn_info.fitem.line, traceback_name=fn_info.fitem.name
    )
    fn_info.callable_class.ir.methods["__call__"] = call_fn_ir
    fn_info.callable_class.ir.method_decls["__call__"] = call_fn_decl
    return call_fn_ir


def add_get_to_callable_class(builder: IRBuilder, fn_info: FuncInfo) -> None:
    """Generate the '__get__' method for a callable class."""
    line = fn_info.fitem.line
    with builder.enter_method(
        fn_info.callable_class.ir,
        "__get__",
        object_rprimitive,
        fn_info,
        self_type=object_rprimitive,
    ):
        instance = builder.add_argument("instance", object_rprimitive)
        builder.add_argument("owner", object_rprimitive)

        # If accessed through the class, just return the callable
        # object. If accessed through an object, create a new bound
        # instance method object.
        instance_block, class_block = BasicBlock(), BasicBlock()
        comparison = builder.translate_is_op(
            builder.read(instance), builder.none_object(line), "is", line
        )
        builder.add_bool_branch(comparison, class_block, instance_block)

        builder.activate_block(class_block)
        builder.add(Return(builder.self()))

        builder.activate_block(instance_block)
        builder.add(
            Return(builder.call_c(method_new_op, [builder.self(), builder.read(instance)], line))
        )


def instantiate_callable_class(builder: IRBuilder, fn_info: FuncInfo) -> Value:
    """Create an instance of a callable class for a function.

    Calls to the function will actually call this instance.

    Note that fn_info refers to the function being assigned, whereas
    builder.fn_info refers to the function encapsulating the function
    being turned into a callable class.
    """
    fitem = fn_info.fitem
    func_reg = builder.add(Call(fn_info.callable_class.ir.ctor, [], fitem.line))

    # Set the environment attribute of the callable class to point at
    # the environment class defined in the callable class' immediate
    # outer scope. Note that there are three possible environment
    # class registers we may use. This depends on what the encapsulating
    # (parent) function is:
    #
    # - A nested function: the callable class is instantiated
    #   from the current callable class' '__call__' function, and hence
    #   the callable class' environment register is used.
    # - A generator function: the callable class is instantiated
    #   from the '__next__' method of the generator class, and hence the
    #   environment of the generator class is used.
    # - Regular function or comprehension scope: we use the environment
    #   of the original function. Comprehension scopes are inlined (no
    #   callable class), so they fall into this case despite is_nested.
    curr_env_reg = None
    if builder.fn_info.is_generator:
        curr_env_reg = builder.fn_info.generator_class.curr_env_reg
    elif builder.fn_info.is_nested and not builder.fn_info.is_comprehension_scope:
        curr_env_reg = builder.fn_info.callable_class.curr_env_reg
    elif builder.fn_info.contains_nested:
        curr_env_reg = builder.fn_info.curr_env_reg
    if curr_env_reg:
        builder.add(SetAttr(func_reg, ENV_ATTR_NAME, curr_env_reg, fitem.line))
    # Initialize function wrapper for callable classes. As opposed to regular functions,
    # each instance of a callable class needs its own wrapper because they might be instantiated
    # inside other functions.
    if not fn_info.in_non_ext and fn_info.is_coroutine:
        builder.add_coroutine_setup_call(fn_info.callable_class.ir.name, func_reg)
    return func_reg
