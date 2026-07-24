"""Generate classes representing function environments (+ related operations).

If we have a nested function that has non-local (free) variables, access to the
non-locals is via an instance of an environment class. Example:

    def f() -> int:
        x = 0  # Make 'x' an attribute of an environment class instance

        def g() -> int:
            # We have access to the environment class instance to
            # allow accessing 'x'
            return x + 2

        x = x + 1  # Modify the attribute
        return g()
"""

from __future__ import annotations

from mypy.nodes import Argument, FuncDef, SymbolNode, Var
from mypyc.common import (
    BITMAP_BITS,
    ENV_ATTR_NAME,
    GENERATOR_ATTRIBUTE_PREFIX,
    SELF_NAME,
    bitmap_name,
)
from mypyc.ir.class_ir import ClassIR
from mypyc.ir.ops import Call, GetAttr, SetAttr, Value
from mypyc.ir.rtypes import RInstance, bitmap_rprimitive, object_rprimitive
from mypyc.irbuild.builder import IRBuilder, SymbolTarget
from mypyc.irbuild.context import FuncInfo, GeneratorClass, ImplicitClass
from mypyc.irbuild.targets import AssignmentTargetAttr


def setup_env_class(builder: IRBuilder) -> ClassIR:
    """Generate a class representing a function environment.

    Note that the variables in the function environment are not
    actually populated here. This is because when the environment
    class is generated, the function environment has not yet been
    visited. This behavior is allowed so that when the compiler visits
    nested functions, it can use the returned ClassIR instance to
    figure out free variables it needs to access.  The remaining
    attributes of the environment class are populated when the
    environment registers are loaded.

    Return a ClassIR representing an environment for a function
    containing a nested function.
    """
    env_class = ClassIR(
        f"{builder.fn_info.namespaced_name()}_env",
        builder.module_name,
        is_generated=True,
        is_final_class=True,
    )
    env_class.reuse_freed_instance = True
    env_class.attributes[SELF_NAME] = RInstance(env_class)
    if builder.fn_info.is_nested and builder.fn_infos[-2]._env_class is not None:
        # If the function is nested, its environment class must contain an environment
        # attribute pointing to its encapsulating functions' environment class.
        env_class.attributes[ENV_ATTR_NAME] = RInstance(builder.fn_infos[-2].env_class)
    env_class.mro = [env_class]
    builder.fn_info.env_class = env_class
    builder.classes.append(env_class)
    return env_class


def finalize_env_class(builder: IRBuilder, prefix: str = "") -> None:
    """Generate, instantiate, and set up the environment of an environment class."""
    if not builder.fn_info.can_merge_generator_and_env_classes():
        instantiate_env_class(builder)

    # Iterate through the function arguments and replace local definitions (using registers)
    # that were previously added to the environment with references to the function's
    # environment class. Comprehension scopes have no arguments to add.
    if not builder.fn_info.is_comprehension_scope:
        if builder.fn_info.is_nested:
            add_args_to_env(
                builder, local=False, base=builder.fn_info.callable_class, prefix=prefix
            )
        else:
            add_args_to_env(builder, local=False, base=builder.fn_info, prefix=prefix)


def instantiate_env_class(builder: IRBuilder) -> Value:
    """Assign an environment class to a register named after the given function definition."""
    curr_env_reg = builder.add(
        Call(builder.fn_info.env_class.ctor, [], builder.fn_info.fitem.line)
    )

    if builder.fn_info.is_nested and not builder.fn_info.is_comprehension_scope:
        builder.fn_info.callable_class._curr_env_reg = curr_env_reg
        builder.add(
            SetAttr(
                curr_env_reg,
                ENV_ATTR_NAME,
                builder.fn_info.callable_class.prev_env_reg,
                builder.fn_info.fitem.line,
            )
        )
    else:
        # Top-level functions and comprehension scopes store env reg directly.
        builder.fn_info._curr_env_reg = curr_env_reg
        # Comprehension scopes link to parent env if it exists.
        if (
            builder.fn_info.is_nested
            and builder.fn_infos[-2]._env_class is not None
            and builder.fn_infos[-2]._curr_env_reg is not None
        ):
            builder.add(
                SetAttr(
                    curr_env_reg,
                    ENV_ATTR_NAME,
                    builder.fn_infos[-2].curr_env_reg,
                    builder.fn_info.fitem.line,
                )
            )

    return curr_env_reg


def load_env_registers(builder: IRBuilder, prefix: str = "") -> None:
    """Load the registers for the current FuncItem being visited.

    Adds the arguments of the FuncItem to the environment. If the
    FuncItem is nested inside of another function, then this also
    loads all of the outer environments of the FuncItem into registers
    so that they can be used when accessing free variables.
    """
    add_args_to_env(builder, local=True, prefix=prefix)

    fn_info = builder.fn_info
    fitem = fn_info.fitem
    if fn_info.is_nested and builder.fn_infos[-2]._env_class is not None:
        load_outer_envs(builder, fn_info.callable_class)
        # If this is a FuncDef, then make sure to load the FuncDef into its own environment
        # class so that the function can be called recursively.
        if isinstance(fitem, FuncDef) and fn_info.add_nested_funcs_to_env:
            setup_func_for_recursive_call(builder, fitem, fn_info.callable_class, prefix=prefix)


def load_outer_env(
    builder: IRBuilder, base: Value, outer_env: dict[SymbolNode, SymbolTarget]
) -> Value:
    """Load the environment class for a given base into a register.

    Additionally, iterates through all of the SymbolNode and
    AssignmentTarget instances of the environment at the given index's
    symtable, and adds those instances to the environment of the
    current environment. This is done so that the current environment
    can access outer environment variables without having to reload
    all of the environment registers.

    Returns the register where the environment class was loaded.
    """
    env = builder.add(GetAttr(base, ENV_ATTR_NAME, builder.fn_info.fitem.line))
    assert isinstance(env.type, RInstance), f"{env} must be of type RInstance"

    for symbol, target in outer_env.items():
        attr_name = symbol.name
        if isinstance(target, AssignmentTargetAttr):
            attr_name = target.attr
        env.type.class_ir.attributes[attr_name] = target.type
        symbol_target = AssignmentTargetAttr(env, attr_name)
        builder.add_target(symbol, symbol_target)

    return env


def load_outer_envs(builder: IRBuilder, base: ImplicitClass) -> None:
    index = len(builder.builders) - 2

    # Load the first outer environment. This one is special because it gets saved in the
    # FuncInfo instance's prev_env_reg field.
    has_outer = index > 1 or (index == 1 and builder.fn_infos[1].contains_nested)
    if has_outer and builder.fn_infos[index]._env_class is not None:
        # outer_env = builder.fn_infos[index].environment
        outer_env = builder.symtables[index]
        if isinstance(base, GeneratorClass):
            base.prev_env_reg = load_outer_env(builder, base.curr_env_reg, outer_env)
        else:
            base.prev_env_reg = load_outer_env(builder, base.self_reg, outer_env)
        env_reg = base.prev_env_reg
        index -= 1

    # Load the remaining outer environments into registers.
    while index > 1:
        if builder.fn_infos[index]._env_class is None:
            break
        # outer_env = builder.fn_infos[index].environment
        outer_env = builder.symtables[index]
        env_reg = load_outer_env(builder, env_reg, outer_env)
        index -= 1


def num_bitmap_args(builder: IRBuilder, args: list[Argument]) -> int:
    n = 0
    for arg in args:
        t = builder.type_to_rtype(arg.variable.type)
        if t.error_overlap and arg.kind.is_optional():
            n += 1
    return (n + (BITMAP_BITS - 1)) // BITMAP_BITS


def add_args_to_env(
    builder: IRBuilder,
    local: bool = True,
    base: FuncInfo | ImplicitClass | None = None,
    reassign: bool = True,
    prefix: str = "",
) -> None:
    fn_info = builder.fn_info
    args = fn_info.fitem.arguments
    nb = num_bitmap_args(builder, args)
    if local:
        for arg in args:
            rtype = builder.type_to_rtype(arg.variable.type)
            builder.add_local_reg(arg.variable, rtype, is_arg=True)
        for i in reversed(range(nb)):
            builder.add_local_reg(Var(bitmap_name(i)), bitmap_rprimitive, is_arg=True)
    else:
        for arg in args:
            if (
                is_free_variable(builder, arg.variable)
                or fn_info.is_generator
                or fn_info.is_coroutine
            ):
                rtype = builder.type_to_rtype(arg.variable.type)
                assert base is not None, "base cannot be None for adding nonlocal args"
                builder.add_var_to_env_class(
                    arg.variable, rtype, base, reassign=reassign, prefix=prefix
                )


def add_vars_to_env(builder: IRBuilder, prefix: str = "") -> None:
    """Add relevant local variables and nested functions to the environment class.

    Add all variables and functions that are declared/defined within current
    function and are referenced in functions nested within this one to this
    function's environment class so the nested functions can reference
    them even if they are declared after the nested function's definition.
    Note that this is done before visiting the body of the function.
    """
    env_for_func: FuncInfo | ImplicitClass = builder.fn_info
    if builder.fn_info.is_generator:
        env_for_func = builder.fn_info.generator_class
    elif (
        builder.fn_info.is_nested or builder.fn_info.in_non_ext
    ) and not builder.fn_info.is_comprehension_scope:
        env_for_func = builder.fn_info.callable_class

    if builder.fn_info.fitem in builder.free_variables:
        # Sort the variables to keep things deterministic
        for var in sorted(builder.free_variables[builder.fn_info.fitem], key=lambda x: x.name):
            if isinstance(var, Var):
                rtype = builder.type_to_rtype(var.type)
                builder.add_var_to_env_class(
                    var, rtype, env_for_func, reassign=False, prefix=prefix
                )

    if builder.fn_info.fitem in builder.encapsulating_funcs:
        for nested_fn in builder.encapsulating_funcs[builder.fn_info.fitem]:
            if isinstance(nested_fn, FuncDef):
                # The return type is 'object' instead of an RInstance of the
                # callable class because differently defined functions with
                # the same name and signature across conditional blocks
                # will generate different callable classes, so the callable
                # class that gets instantiated must be generic.
                if nested_fn.is_generator or nested_fn.is_coroutine:
                    prefix = GENERATOR_ATTRIBUTE_PREFIX
                builder.add_var_to_env_class(
                    nested_fn, object_rprimitive, env_for_func, reassign=False, prefix=prefix
                )


def setup_func_for_recursive_call(
    builder: IRBuilder, fdef: FuncDef, base: ImplicitClass, prefix: str = ""
) -> None:
    """Enable calling a nested function (with a callable class) recursively.

    Adds the instance of the callable class representing the given
    FuncDef to a register in the environment so that the function can
    be called recursively. Note that this needs to be done only for
    nested functions.
    """
    # First, set the attribute of the environment class so that GetAttr can be called on it.
    prev_env = builder.fn_infos[-2].env_class
    attr_name = prefix + fdef.name
    prev_env.attributes[attr_name] = builder.type_to_rtype(fdef.type)
    line = fdef.line

    if isinstance(base, GeneratorClass):
        # If we are dealing with a generator class, then we need to first get the register
        # holding the current environment class, and load the previous environment class from
        # there.
        prev_env_reg = builder.add(GetAttr(base.curr_env_reg, ENV_ATTR_NAME, line))
    else:
        prev_env_reg = base.prev_env_reg

    # Obtain the instance of the callable class representing the FuncDef, and add it to the
    # current environment.
    val = builder.add(GetAttr(prev_env_reg, attr_name, line))
    target = builder.add_local_reg(fdef, object_rprimitive)
    builder.assign(target, val, line)


def is_free_variable(builder: IRBuilder, symbol: SymbolNode) -> bool:
    fitem = builder.fn_info.fitem
    return fitem in builder.free_variables and symbol in builder.free_variables[fitem]
