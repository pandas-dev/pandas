"""Special case IR generation of calls to specific builtin functions.

Most special cases should be handled using the data driven "primitive
ops" system, but certain operations require special handling that has
access to the AST/IR directly and can make decisions/optimizations
based on it. These special cases can be implemented here.

For example, we use specializers to statically emit the length of a
fixed length tuple and to emit optimized code for any()/all() calls with
generator comprehensions as the argument.

See comment below for more documentation.
"""

from typing import Callable, Optional, Dict, Tuple, List

from mypy.nodes import CallExpr, RefExpr, MemberExpr, TupleExpr, GeneratorExpr, ARG_POS
from mypy.types import AnyType, TypeOfAny

from mypyc.ir.ops import (
    Value, Register, BasicBlock, Integer, RaiseStandardError, Unreachable
)
from mypyc.ir.rtypes import (
    RType, RTuple, str_rprimitive, list_rprimitive, dict_rprimitive, set_rprimitive,
    bool_rprimitive, is_dict_rprimitive
)
from mypyc.primitives.dict_ops import dict_keys_op, dict_values_op, dict_items_op
from mypyc.primitives.list_ops import new_list_set_item_op
from mypyc.primitives.tuple_ops import new_tuple_set_item_op
from mypyc.irbuild.builder import IRBuilder
from mypyc.irbuild.for_helpers import (
    translate_list_comprehension, translate_set_comprehension,
    comprehension_helper, sequence_from_generator_preallocate_helper
)


# Specializers are attempted before compiling the arguments to the
# function.  Specializers can return None to indicate that they failed
# and the call should be compiled normally. Otherwise they should emit
# code for the call and return a Value containing the result.
#
# Specializers take three arguments: the IRBuilder, the CallExpr being
# compiled, and the RefExpr that is the left hand side of the call.
Specializer = Callable[['IRBuilder', CallExpr, RefExpr], Optional[Value]]

# Dictionary containing all configured specializers.
#
# Specializers can operate on methods as well, and are keyed on the
# name and RType in that case.
specializers = {}  # type: Dict[Tuple[str, Optional[RType]], List[Specializer]]


def specialize_function(
        name: str, typ: Optional[RType] = None) -> Callable[[Specializer], Specializer]:
    """Decorator to register a function as being a specializer.

    There may exist multiple specializers for one function. When translating method
    calls, the earlier appended specializer has higher priority.
    """
    def wrapper(f: Specializer) -> Specializer:
        specializers.setdefault((name, typ), []).append(f)
        return f
    return wrapper


@specialize_function('builtins.globals')
def translate_globals(builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case builtins.globals
    if len(expr.args) == 0:
        return builder.load_globals_dict()
    return None


@specialize_function('builtins.len')
def translate_len(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case builtins.len
    if (len(expr.args) == 1
            and expr.arg_kinds == [ARG_POS]):
        expr_rtype = builder.node_type(expr.args[0])
        if isinstance(expr_rtype, RTuple):
            # len() of fixed-length tuple can be trivially determined statically,
            # though we still need to evaluate it.
            builder.accept(expr.args[0])
            return Integer(len(expr_rtype.types))
        else:
            obj = builder.accept(expr.args[0])
            return builder.builtin_len(obj, expr.line)
    return None


@specialize_function('builtins.list')
def dict_methods_fast_path(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Specialize a common case when list() is called on a dictionary view
    # method call, for example foo = list(bar.keys()).
    if not (len(expr.args) == 1 and expr.arg_kinds == [ARG_POS]):
        return None
    arg = expr.args[0]
    if not (isinstance(arg, CallExpr) and not arg.args
            and isinstance(arg.callee, MemberExpr)):
        return None
    base = arg.callee.expr
    attr = arg.callee.name
    rtype = builder.node_type(base)
    if not (is_dict_rprimitive(rtype) and attr in ('keys', 'values', 'items')):
        return None

    obj = builder.accept(base)
    # Note that it is not safe to use fast methods on dict subclasses, so
    # the corresponding helpers in CPy.h fallback to (inlined) generic logic.
    if attr == 'keys':
        return builder.call_c(dict_keys_op, [obj], expr.line)
    elif attr == 'values':
        return builder.call_c(dict_values_op, [obj], expr.line)
    else:
        return builder.call_c(dict_items_op, [obj], expr.line)


@specialize_function('builtins.list')
def translate_list_from_generator_call(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case for simplest list comprehension, for example
    #     list(f(x) for x in other_list/other_tuple)
    # translate_list_comprehension() would take care of other cases if this fails.
    if (len(expr.args) == 1
            and expr.arg_kinds[0] == ARG_POS
            and isinstance(expr.args[0], GeneratorExpr)):
        return sequence_from_generator_preallocate_helper(
            builder, expr.args[0],
            empty_op_llbuilder=builder.builder.new_list_op_with_length,
            set_item_op=new_list_set_item_op)
    return None


@specialize_function('builtins.tuple')
def translate_tuple_from_generator_call(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case for simplest tuple creation from a generator, for example
    #     tuple(f(x) for x in other_list/other_tuple)
    # translate_safe_generator_call() would take care of other cases if this fails.
    if (len(expr.args) == 1
            and expr.arg_kinds[0] == ARG_POS
            and isinstance(expr.args[0], GeneratorExpr)):
        return sequence_from_generator_preallocate_helper(
            builder, expr.args[0],
            empty_op_llbuilder=builder.builder.new_tuple_with_length,
            set_item_op=new_tuple_set_item_op)
    return None


@specialize_function('builtins.set')
def translate_set_from_generator_call(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case for set creation from a generator：
    #     set(f(...) for ... in iterator/nested_generators...)
    if (len(expr.args) == 1
            and expr.arg_kinds[0] == ARG_POS
            and isinstance(expr.args[0], GeneratorExpr)):
        return translate_set_comprehension(builder, expr.args[0])
    return None


@specialize_function('builtins.tuple')
@specialize_function('builtins.frozenset')
@specialize_function('builtins.dict')
@specialize_function('builtins.sum')
@specialize_function('builtins.min')
@specialize_function('builtins.max')
@specialize_function('builtins.sorted')
@specialize_function('collections.OrderedDict')
@specialize_function('join', str_rprimitive)
@specialize_function('extend', list_rprimitive)
@specialize_function('update', dict_rprimitive)
@specialize_function('update', set_rprimitive)
def translate_safe_generator_call(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special cases for things that consume iterators where we know we
    # can safely compile a generator into a list.
    if (len(expr.args) > 0
            and expr.arg_kinds[0] == ARG_POS
            and isinstance(expr.args[0], GeneratorExpr)):
        if isinstance(callee, MemberExpr):
            return builder.gen_method_call(
                builder.accept(callee.expr), callee.name,
                ([translate_list_comprehension(builder, expr.args[0])]
                    + [builder.accept(arg) for arg in expr.args[1:]]),
                builder.node_type(expr), expr.line, expr.arg_kinds, expr.arg_names)
        else:
            return builder.call_refexpr_with_args(
                expr, callee,
                ([translate_list_comprehension(builder, expr.args[0])]
                    + [builder.accept(arg) for arg in expr.args[1:]]))
    return None


@specialize_function('builtins.any')
def translate_any_call(builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    if (len(expr.args) == 1
            and expr.arg_kinds == [ARG_POS]
            and isinstance(expr.args[0], GeneratorExpr)):
        return any_all_helper(builder, expr.args[0], builder.false, lambda x: x, builder.true)
    return None


@specialize_function('builtins.all')
def translate_all_call(builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    if (len(expr.args) == 1
            and expr.arg_kinds == [ARG_POS]
            and isinstance(expr.args[0], GeneratorExpr)):
        return any_all_helper(
            builder, expr.args[0],
            builder.true,
            lambda x: builder.unary_op(x, 'not', expr.line),
            builder.false
        )
    return None


def any_all_helper(builder: IRBuilder,
                   gen: GeneratorExpr,
                   initial_value: Callable[[], Value],
                   modify: Callable[[Value], Value],
                   new_value: Callable[[], Value]) -> Value:
    retval = Register(bool_rprimitive)
    builder.assign(retval, initial_value(), -1)
    loop_params = list(zip(gen.indices, gen.sequences, gen.condlists))
    true_block, false_block, exit_block = BasicBlock(), BasicBlock(), BasicBlock()

    def gen_inner_stmts() -> None:
        comparison = modify(builder.accept(gen.left_expr))
        builder.add_bool_branch(comparison, true_block, false_block)
        builder.activate_block(true_block)
        builder.assign(retval, new_value(), -1)
        builder.goto(exit_block)
        builder.activate_block(false_block)

    comprehension_helper(builder, loop_params, gen_inner_stmts, gen.line)
    builder.goto_and_activate(exit_block)

    return retval


@specialize_function('dataclasses.field')
@specialize_function('attr.Factory')
def translate_dataclasses_field_call(
        builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case for 'dataclasses.field' and 'attr.Factory' function calls
    # because the results of such calls are typechecked by mypy using the types
    # of the arguments to their respective functions, resulting in attempted
    # coercions by mypyc that throw a runtime error.
    builder.types[expr] = AnyType(TypeOfAny.from_error)
    return None


@specialize_function('builtins.next')
def translate_next_call(builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    # Special case for calling next() on a generator expression, an
    # idiom that shows up some in mypy.
    #
    # For example, next(x for x in l if x.id == 12, None) will
    # generate code that searches l for an element where x.id == 12
    # and produce the first such object, or None if no such element
    # exists.
    if not (expr.arg_kinds in ([ARG_POS], [ARG_POS, ARG_POS])
            and isinstance(expr.args[0], GeneratorExpr)):
        return None

    gen = expr.args[0]

    retval = Register(builder.node_type(expr))
    default_val = None
    if len(expr.args) > 1:
        default_val = builder.accept(expr.args[1])

    exit_block = BasicBlock()

    def gen_inner_stmts() -> None:
        # next takes the first element of the generator, so if
        # something gets produced, we are done.
        builder.assign(retval, builder.accept(gen.left_expr), gen.left_expr.line)
        builder.goto(exit_block)

    loop_params = list(zip(gen.indices, gen.sequences, gen.condlists))
    comprehension_helper(builder, loop_params, gen_inner_stmts, gen.line)

    # Now we need the case for when nothing got hit. If there was
    # a default value, we produce it, and otherwise we raise
    # StopIteration.
    if default_val:
        builder.assign(retval, default_val, gen.left_expr.line)
        builder.goto(exit_block)
    else:
        builder.add(RaiseStandardError(RaiseStandardError.STOP_ITERATION, None, expr.line))
        builder.add(Unreachable())

    builder.activate_block(exit_block)
    return retval


@specialize_function('builtins.isinstance')
def translate_isinstance(builder: IRBuilder, expr: CallExpr, callee: RefExpr) -> Optional[Value]:
    if (len(expr.args) == 2
            and expr.arg_kinds == [ARG_POS, ARG_POS]
            and isinstance(expr.args[1], (RefExpr, TupleExpr))):
        # Special case for builtins.isinstance
        # Prevent coercions on the thing we are checking the instance of - there is no need to
        # coerce something to a new type before checking what type it is, and the coercion could
        # lead to bugs.
        builder.types[expr.args[0]] = AnyType(TypeOfAny.from_error)

        irs = builder.flatten_classes(expr.args[1])
        if irs is not None:
            return builder.builder.isinstance_helper(builder.accept(expr.args[0]), irs, expr.line)
    return None
