import builtins
import collections
import dis
import operator
import logging
import textwrap

from numba.core import errors, ir, config
from numba.core.errors import (
    NotDefinedError,
    UnsupportedBytecodeError,
    error_extras,
)
from numba.core.ir_utils import get_definition, guard
from numba.core.utils import (PYVERSION, BINOPS_TO_OPERATORS,
                              INPLACE_BINOPS_TO_OPERATORS, _lazy_pformat,)
from numba.core.byteflow import Flow, AdaptDFA, AdaptCFA, BlockKind
from numba.core.unsafe import eh
from numba.cpython.unsafe.tuple import unpack_single_tuple


if PYVERSION in ((3, 12), (3, 13)):
    # Operands for CALL_INTRINSIC_1
    from numba.core.byteflow import CALL_INTRINSIC_1_Operand as ci1op
elif PYVERSION in ((3, 10), (3, 11)):
    pass
else:
    raise NotImplementedError(PYVERSION)


class _UNKNOWN_VALUE(object):
    """Represents an unknown value, this is for ease of debugging purposes only.
    """

    def __init__(self, varname):
        self._varname = varname

    def __repr__(self):
        return "_UNKNOWN_VALUE({})".format(self._varname)


_logger = logging.getLogger(__name__)


class Assigner(object):
    """
    This object keeps track of potential assignment simplifications
    inside a code block.
    For example `$O.1 = x` followed by `y = $0.1` can be simplified
    into `y = x`, but it's not possible anymore if we have `x = z`
    in-between those two instructions.

    NOTE: this is not only an optimization, but is actually necessary
    due to certain limitations of Numba - such as only accepting the
    returning of an array passed as function argument.
    """

    def __init__(self):
        # { destination variable name -> source Var object }
        self.dest_to_src = {}
        # Basically a reverse mapping of dest_to_src:
        # { source variable name -> all destination names in dest_to_src }
        self.src_invalidate = collections.defaultdict(list)
        self.unused_dests = set()

    def assign(self, srcvar, destvar):
        """
        Assign *srcvar* to *destvar*. Return either *srcvar* or a possible
        simplified assignment source (earlier assigned to *srcvar*).
        """
        srcname = srcvar.name
        destname = destvar.name
        if destname in self.src_invalidate:
            # destvar will change, invalidate all previously known
            # simplifications
            for d in self.src_invalidate.pop(destname):
                self.dest_to_src.pop(d)
        if srcname in self.dest_to_src:
            srcvar = self.dest_to_src[srcname]
        if destvar.is_temp:
            self.dest_to_src[destname] = srcvar
            self.src_invalidate[srcname].append(destname)
            self.unused_dests.add(destname)
        return srcvar

    def get_assignment_source(self, destname):
        """
        Get a possible assignment source (a ir.Var instance) to replace
        *destname*, otherwise None.
        """
        if destname in self.dest_to_src:
            return self.dest_to_src[destname]
        self.unused_dests.discard(destname)
        return None


def _remove_assignment_definition(old_body, idx, func_ir, already_deleted_defs):
    """
    Deletes the definition defined for old_body at index idx
    from func_ir. We assume this stmt will be deleted from
    new_body.

    In some optimizations we may update the same variable multiple times.
    In this situation, we only need to delete a particular definition once,
    this is tracked in already_deleted_def, which is a map from
    assignment name to the set of values that have already been
    deleted.
    """
    lhs = old_body[idx].target.name
    rhs = old_body[idx].value
    if rhs in func_ir._definitions[lhs]:
        func_ir._definitions[lhs].remove(rhs)
        already_deleted_defs[lhs].add(rhs)
    elif rhs not in already_deleted_defs[lhs]:
        raise UnsupportedBytecodeError(
            "Inconsistency found in the definitions while executing"
            " a peephole optimization. This suggests an internal"
            " error or inconsistency elsewhere in the compiler."
        )


def _call_function_ex_replace_kws_small(
    old_body,
    keyword_expr,
    new_body,
    buildmap_idx,
    func_ir,
    already_deleted_defs
):
    """
    Extracts the kws args passed as varkwarg
    for CALL_FUNCTION_EX. This pass is taken when
    n_kws <= 15 and the bytecode looks like:

        # Start for each argument
        LOAD_FAST  # Load each argument.
        # End for each argument
        ...
        BUILD_CONST_KEY_MAP # Build a map

    In the generated IR, the varkwarg refers
    to a single build_map that contains all of the
    kws. In addition to returning the kws, this
    function updates new_body to remove all usage
    of the map.
    """
    kws = keyword_expr.items.copy()
    # kws are required to have constant keys.
    # We update these with the value_indexes
    value_indexes = keyword_expr.value_indexes
    for key, index in value_indexes.items():
        kws[index] = (key, kws[index][1])
    # Remove the build_map by setting the list
    # index to None. Nones will be removed later.
    new_body[buildmap_idx] = None
    # Remove the definition.
    _remove_assignment_definition(
        old_body, buildmap_idx, func_ir, already_deleted_defs
    )
    return kws


def _call_function_ex_replace_kws_large(
    old_body,
    buildmap_name,
    buildmap_idx,
    search_end,
    new_body,
    func_ir,
    errmsg,
    already_deleted_defs
):
    """
    Extracts the kws args passed as varkwarg
    for CALL_FUNCTION_EX. This pass is taken when
    n_kws > 15 and the bytecode looks like:

        BUILD_MAP # Construct the map
        # Start for each argument
        LOAD_CONST # Load a constant for the name of the argument
        LOAD_FAST  # Load each argument.
        MAP_ADD # Append the (key, value) pair to the map
        # End for each argument

    In the IR generated, the initial build map is empty and a series
    of setitems are applied afterwards. THE IR looks like:

        $build_map_var = build_map(items=[])
        $constvar = const(str, ...) # create the const key
        # CREATE THE ARGUMENT, This may take multiple lines.
        $created_arg = ...
        $var = getattr(
            value=$build_map_var,
            attr=__setitem__,
        )
        $unused_var = call $var($constvar, $created_arg)

    We iterate through the IR, deleting all usages of the buildmap
    from the new_body, and adds the kws to a new kws list.
    """
    # Remove the build_map from the body.
    new_body[buildmap_idx] = None
    # Remove the definition.
    _remove_assignment_definition(
        old_body, buildmap_idx, func_ir, already_deleted_defs
    )
    kws = []
    search_start = buildmap_idx + 1
    while search_start <= search_end:
        # The first value must be a constant.
        const_stmt = old_body[search_start]
        if not (
            isinstance(const_stmt, ir.Assign)
            and isinstance(const_stmt.value, ir.Const)
        ):
            # We cannot handle this format so raise the
            # original error message.
            raise UnsupportedBytecodeError(errmsg)
        key_var_name = const_stmt.target.name
        key_val = const_stmt.value.value
        search_start += 1
        # Now we need to search for a getattr with setitem
        found_getattr = False
        while (
            search_start <= search_end
            and not found_getattr
        ):
            getattr_stmt = old_body[search_start]
            if (
                isinstance(getattr_stmt, ir.Assign)
                and isinstance(getattr_stmt.value, ir.Expr)
                and getattr_stmt.value.op == "getattr"
                and (
                    getattr_stmt.value.value.name
                    == buildmap_name
                )
                and getattr_stmt.value.attr == "__setitem__"
            ):
                found_getattr = True
            else:
                # If the argument is "created" in JIT, then there
                # will be intermediate operations in between setitems.
                # For example we have arg5=pow(arg5, 2),
                # then the IR would look like:
                #
                #   # Creation of the constant key.
                #   $const44.26 = const(str, arg5)
                #
                #   # Argument creation. This is the section we are skipping
                #   $46load_global.27 = global(pow: <built-in function pow>)
                #   $const50.29 = const(int, 2)
                #   $call.30 = call $46load_global.27(arg5, $const50.29)
                #
                #   # Setitem with arg5
                #   $54map_add.31 = getattr(value=$map.2, attr=__setitem__)
                #   $54map_add.32 = call $54map_add.31($const44.26, $call.30)
                search_start += 1
        if (
            not found_getattr
            or search_start == search_end
        ):
            # We cannot handle this format so raise the
            # original error message.
            raise UnsupportedBytecodeError(errmsg)
        setitem_stmt = old_body[search_start + 1]
        if not (
            isinstance(setitem_stmt, ir.Assign)
            and isinstance(setitem_stmt.value, ir.Expr)
            and setitem_stmt.value.op == "call"
            and (
                setitem_stmt.value.func.name
                == getattr_stmt.target.name
            )
            and len(setitem_stmt.value.args) == 2
            and (
                setitem_stmt.value.args[0].name
                == key_var_name
            )
        ):
            # A call statement should always immediately follow the
            # getattr. If for some reason this doesn't match the code
            # format, we raise the original error message. This check
            # is meant as a precaution.
            raise UnsupportedBytecodeError(errmsg)
        arg_var = setitem_stmt.value.args[1]
        # Append the (key, value) pair.
        kws.append((key_val, arg_var))
        # Remove the __setitem__ getattr and call
        new_body[search_start] = None
        new_body[search_start + 1] = None
        # Remove the definitions.
        _remove_assignment_definition(
            old_body, search_start, func_ir, already_deleted_defs
        )
        _remove_assignment_definition(
            old_body, search_start + 1, func_ir, already_deleted_defs
        )
        search_start += 2
    return kws


def _call_function_ex_replace_args_small(
    old_body,
    tuple_expr,
    new_body,
    buildtuple_idx,
    func_ir,
    already_deleted_defs
):
    """
    Extracts the args passed as vararg
    for CALL_FUNCTION_EX. This pass is taken when
    n_args <= 30 and the bytecode looks like:

        # Start for each argument
        LOAD_FAST  # Load each argument.
        # End for each argument
        ...
        BUILD_TUPLE # Create a tuple of the arguments

    In the IR generated, the vararg refer
    to a single build_tuple that contains all of the
    args. In addition to returning the args, this
    function updates new_body to remove all usage
    of the tuple.
    """
    # Delete the build tuple
    new_body[buildtuple_idx] = None
    # Remove the definition.
    _remove_assignment_definition(
        old_body, buildtuple_idx, func_ir, already_deleted_defs
    )
    # Return the args.
    return tuple_expr.items


def _call_function_ex_replace_args_large(
    old_body,
    vararg_stmt,
    new_body,
    search_end,
    func_ir,
    errmsg,
    already_deleted_defs
):
    """
    Extracts the args passed as vararg
    for CALL_FUNCTION_EX. This pass is taken when
    n_args > 30 and the bytecode looks like:

        BUILD_TUPLE # Create a list to append to
        # Start for each argument
        LOAD_FAST  # Load each argument.
        LIST_APPEND # Add the argument to the list
        # End for each argument
        ...
        LIST_TO_TUPLE # Convert the args to a tuple.

    In the IR generated, the tuple is created by concatenating
    together several 1 element tuples to an initial empty tuple.
    We traverse backwards in the IR, collecting args, until we
    find the original empty tuple. For example, the IR might
    look like:

        $orig_tuple = build_tuple(items=[])
        $first_var = build_tuple(items=[Var(arg0, test.py:6)])
        $next_tuple = $orig_tuple + $first_var
        ...
        $final_var = build_tuple(items=[Var(argn, test.py:6)])
        $final_tuple = $prev_tuple + $final_var
        $varargs_var = $final_tuple
    """
    # We traverse to the front of the block to look for the original
    # tuple.
    search_start = 0
    total_args = []
    if (
        isinstance(vararg_stmt, ir.Assign)
        and isinstance(vararg_stmt.value, ir.Var)
    ):
        target_name = vararg_stmt.value.name
        # If there is an initial assignment, delete it
        new_body[search_end] = None
        # Remove the definition.
        _remove_assignment_definition(
            old_body, search_end, func_ir, already_deleted_defs
        )
        search_end -= 1
    else:
        # There must always be an initial assignment
        # https://github.com/numba/numba/blob/59fa2e335be68148b3bd72a29de3ff011430038d/numba/core/interpreter.py#L259-L260
        # If this changes we may need to support this branch.
        raise AssertionError("unreachable")
    # Traverse backwards to find all concatenations
    # until eventually reaching the original empty tuple.
    while search_end >= search_start:
        concat_stmt = old_body[search_end]
        if (
            isinstance(concat_stmt, ir.Assign)
            and concat_stmt.target.name == target_name
            and isinstance(concat_stmt.value, ir.Expr)
            and concat_stmt.value.op == "build_tuple"
            and not concat_stmt.value.items
        ):
            new_body[search_end] = None
            # Remove the definition.
            _remove_assignment_definition(
                old_body, search_end, func_ir, already_deleted_defs
            )
            # If we have reached the build_tuple we exit.
            break
        else:
            # We expect to find another arg to append.
            # The first stmt must be a binop "add"
            if (search_end == search_start) or not (
                isinstance(concat_stmt, ir.Assign)
                and (
                    concat_stmt.target.name
                    == target_name
                )
                and isinstance(
                    concat_stmt.value, ir.Expr
                )
                and concat_stmt.value.op == "binop"
                and concat_stmt.value.fn == operator.add
            ):
                # We cannot handle this format.
                raise UnsupportedBytecodeError(errmsg)
            lhs_name = concat_stmt.value.lhs.name
            rhs_name = concat_stmt.value.rhs.name
            # The previous statement should be a
            # build_tuple containing the arg.
            arg_tuple_stmt = old_body[search_end - 1]
            if not (
                isinstance(arg_tuple_stmt, ir.Assign)
                and isinstance(
                    arg_tuple_stmt.value, ir.Expr
                )
                and (
                    arg_tuple_stmt.value.op
                    == "build_tuple"
                )
                and len(arg_tuple_stmt.value.items) == 1
            ):
                # We cannot handle this format.
                raise UnsupportedBytecodeError(errmsg)
            if arg_tuple_stmt.target.name == lhs_name:
                # The tuple should always be generated on the RHS.
                raise AssertionError("unreachable")
            elif arg_tuple_stmt.target.name == rhs_name:
                target_name = lhs_name
            else:
                # We cannot handle this format.
                raise UnsupportedBytecodeError(errmsg)
            total_args.append(
                arg_tuple_stmt.value.items[0]
            )
            new_body[search_end] = None
            new_body[search_end - 1] = None
            # Remove the definitions.
            _remove_assignment_definition(
                old_body, search_end, func_ir, already_deleted_defs
            )
            _remove_assignment_definition(
                old_body, search_end - 1, func_ir, already_deleted_defs
            )
            search_end -= 2
            # Avoid any space between appends
            keep_looking = True
            while search_end >= search_start and keep_looking:
                next_stmt = old_body[search_end]
                if (
                    isinstance(next_stmt, ir.Assign)
                    and (
                        next_stmt.target.name
                        == target_name
                    )
                ):
                    keep_looking = False
                else:
                    # If the argument is "created" in JIT, then there
                    # will be intermediate operations in between appends.
                    # For example if the next arg after arg4 is pow(arg5, 2),
                    # then the IR would look like:
                    #
                    #   # Appending arg4
                    #   $arg4_tup = build_tuple(items=[arg4])
                    #   $append_var.5 = $append_var.4 + $arg4_tup
                    #
                    #   # Creation of arg5.
                    #   # This is the section that we are skipping.
                    #   $32load_global.20 = global(pow: <built-in function pow>)
                    #   $const36.22 = const(int, 2)
                    #   $call.23 = call $32load_global.20(arg5, $const36.22)
                    #
                    #   # Appending arg5
                    #   $arg5_tup = build_tuple(items=[$call.23])
                    #   $append_var.6 = $append_var.5 + $arg5_tup
                    search_end -= 1
    if search_end == search_start:
        # If we reached the start we never found the build_tuple.
        # We cannot handle this format so raise the
        # original error message.
        raise UnsupportedBytecodeError(errmsg)
    # Reverse the arguments so we get the correct order.
    return total_args[::-1]


def peep_hole_call_function_ex_to_call_function_kw(func_ir):
    """
    This peephole rewrites a bytecode sequence unique to Python 3.10
    where CALL_FUNCTION_EX is used instead of CALL_FUNCTION_KW because of
    stack limitations set by CPython. This limitation is imposed whenever
    a function call has too many arguments or keyword arguments.

    https://github.com/python/cpython/blob/a58ebcc701dd6c43630df941481475ff0f615a81/Python/compile.c#L55
    https://github.com/python/cpython/blob/a58ebcc701dd6c43630df941481475ff0f615a81/Python/compile.c#L4442

    In particular, this change is imposed whenever (n_args / 2) + n_kws > 15.

    Different bytecode is generated for args depending on if n_args > 30
    or n_args <= 30 and similarly if n_kws > 15 or n_kws <= 15.

    This function unwraps the *args and **kwargs in the function call
    and places these values directly into the args and kwargs of the call.
    """
    # All changes are local to the a single block
    # so it can be traversed in any order.
    errmsg = textwrap.dedent("""
        CALL_FUNCTION_EX with **kwargs not supported.
        If you are not using **kwargs this may indicate that
        you have a large number of kwargs and are using inlined control
        flow. You can resolve this issue by moving the control flow out of
        the function call. For example, if you have

            f(a=1 if flag else 0, ...)

        Replace that with:

            a_val = 1 if flag else 0
            f(a=a_val, ...)""")

    # Track which definitions have already been deleted
    already_deleted_defs = collections.defaultdict(set)
    for blk in func_ir.blocks.values():
        blk_changed = False
        new_body = []
        for i, stmt in enumerate(blk.body):
            if (
                isinstance(stmt, ir.Assign)
                and isinstance(stmt.value, ir.Expr)
                and stmt.value.op == "call"
                and stmt.value.varkwarg is not None
            ):
                blk_changed = True
                call = stmt.value
                args = call.args
                kws = call.kws
                # We need to check the call expression contents if
                # it contains either vararg or varkwarg. If it contains
                # varkwarg we need to update the IR. If it just contains
                # vararg we don't need to update the IR, but we need to
                # check if peep_hole_list_to_tuple failed to replace the
                # vararg list with a tuple. If so, we output an error
                # message with suggested code changes.
                vararg = call.vararg
                varkwarg = call.varkwarg
                start_search = i - 1
                # varkwarg should be defined second so we start there.
                varkwarg_loc = start_search
                keyword_def = None
                found = False
                while varkwarg_loc >= 0 and not found:
                    keyword_def = blk.body[varkwarg_loc]
                    if (
                        isinstance(keyword_def, ir.Assign)
                        and keyword_def.target.name == varkwarg.name
                    ):
                        found = True
                    else:
                        varkwarg_loc -= 1
                if (
                    kws
                    or not found
                    or not (
                        isinstance(keyword_def.value, ir.Expr)
                        and keyword_def.value.op == "build_map"
                    )
                ):
                    # If we couldn't find where the kwargs are created
                    # then it should be a normal **kwargs call
                    # so we produce an unsupported message.
                    raise UnsupportedBytecodeError(errmsg)
                # Determine the kws
                if keyword_def.value.items:
                    # n_kws <= 15 case.
                    # Here the IR looks like a series of
                    # constants, then the arguments and finally
                    # a build_map that contains all of the pairs.
                    # For Example:
                    #
                    #   $const_n = const("arg_name")
                    #   $arg_n = ...
                    #   $kwargs_var = build_map(items=[
                    #              ($const_0, $arg_0),
                    #              ...,
                    #              ($const_n, $arg_n),])
                    kws = _call_function_ex_replace_kws_small(
                        blk.body,
                        keyword_def.value,
                        new_body,
                        varkwarg_loc,
                        func_ir,
                        already_deleted_defs,
                    )
                else:
                    # n_kws > 15 case.
                    # Here the IR is an initial empty build_map
                    # followed by a series of setitems with a constant
                    # key and then the argument.
                    # For example:
                    #
                    #   $kwargs_var = build_map(items=[])
                    #   $const_0 = const("arg_name")
                    #   $arg_0 = ...
                    #   $my_attr = getattr(const_0, attr=__setitem__)
                    #   $unused_var = call $my_attr($const_0, $arg_0)
                    #   ...
                    kws = _call_function_ex_replace_kws_large(
                        blk.body,
                        varkwarg.name,
                        varkwarg_loc,
                        i - 1,
                        new_body,
                        func_ir,
                        errmsg,
                        already_deleted_defs,
                    )
                start_search = varkwarg_loc
                # Vararg isn't required to be provided.
                if vararg is not None:
                    if args:
                        # If we have vararg then args is expected to
                        # be an empty list.
                        raise UnsupportedBytecodeError(errmsg)
                    vararg_loc = start_search
                    args_def = None
                    found = False
                    while vararg_loc >= 0 and not found:
                        args_def = blk.body[vararg_loc]
                        if (
                            isinstance(args_def, ir.Assign)
                            and args_def.target.name == vararg.name
                        ):
                            found = True
                        else:
                            vararg_loc -= 1
                    if not found:
                        # If we couldn't find where the args are created
                        # then we can't handle this format.
                        raise UnsupportedBytecodeError(errmsg)
                    if (
                        isinstance(args_def.value, ir.Expr)
                        and args_def.value.op == "build_tuple"
                    ):
                        # n_args <= 30 case.
                        # Here the IR is a simple build_tuple containing
                        # all of the args.
                        # For example:
                        #
                        #  $arg_n = ...
                        #  $varargs = build_tuple(
                        #   items=[$arg_0, ..., $arg_n]
                        #  )
                        args = _call_function_ex_replace_args_small(
                            blk.body,
                            args_def.value,
                            new_body,
                            vararg_loc,
                            func_ir,
                            already_deleted_defs,
                        )
                    elif (
                        isinstance(args_def.value, ir.Expr)
                        and args_def.value.op == "list_to_tuple"
                    ):
                        # If there is a call with vararg we need to check
                        # if the list -> tuple conversion failed and if so
                        # throw an error.
                        raise UnsupportedBytecodeError(errmsg)
                    else:
                        # Here the IR is an initial empty build_tuple.
                        # Then for each arg, a new tuple with a single
                        # element is created and one by one these are
                        # added to a growing tuple.
                        # For example:
                        #
                        #  $combo_tup_0 = build_tuple(items=[])
                        #  $arg0 = ...
                        #  $arg0_tup = build_tuple(items=[$arg0])
                        #  $combo_tup_1 = $combo_tup_0 + $arg0_tup
                        #  $arg1 = ...
                        #  $arg1_tup = build_tuple(items=[$arg1])
                        #  $combo_tup_2 = $combo_tup_1 + $arg1_tup
                        #  ...
                        #  $combo_tup_n = $combo_tup_{n-1} + $argn_tup
                        #
                        # In addition, the IR contains a final
                        # assignment for the varargs that looks like:
                        #
                        #  $varargs_var = $combo_tup_n
                        #
                        # Here args_def is expected to be a simple assignment.
                        args = _call_function_ex_replace_args_large(
                            blk.body,
                            args_def,
                            new_body,
                            vararg_loc,
                            func_ir,
                            errmsg,
                            already_deleted_defs,
                        )
                # Create a new call updating the args and kws
                new_call = ir.Expr.call(
                    call.func, args, kws, call.loc, target=call.target
                )
                # Drop the existing definition for this stmt.
                _remove_assignment_definition(
                    blk.body, i, func_ir, already_deleted_defs
                )
                # Update the statement
                stmt = ir.Assign(new_call, stmt.target, stmt.loc)
                # Update the definition
                func_ir._definitions[stmt.target.name].append(new_call)
            elif (
                isinstance(stmt, ir.Assign)
                and isinstance(stmt.value, ir.Expr)
                and stmt.value.op == "call"
                and stmt.value.vararg is not None
            ):
                # If there is a call with vararg we need to check
                # if the list -> tuple conversion failed and if so
                # throw an error.
                call = stmt.value
                vararg_name = call.vararg.name
                if (
                    vararg_name in func_ir._definitions
                    and len(func_ir._definitions[vararg_name]) == 1
                ):
                    # If this value is still a list to tuple raise the
                    # exception.
                    expr = func_ir._definitions[vararg_name][0]
                    if isinstance(expr, ir.Expr) and expr.op == "list_to_tuple":
                        raise UnsupportedBytecodeError(errmsg)

            new_body.append(stmt)
        # Replace the block body if we changed the IR
        if blk_changed:
            blk.body.clear()
            blk.body.extend([x for x in new_body if x is not None])
    return func_ir


def peep_hole_list_to_tuple(func_ir):
    """
    This peephole rewrites a bytecode sequence new to Python 3.9 that looks
    like e.g.:

    def foo(a):
        return (*a,)

    41          0 BUILD_LIST               0
                2 LOAD_FAST                0 (a)
                4 LIST_EXTEND              1
                6 LIST_TO_TUPLE
                8 RETURN_VAL

    essentially, the unpacking of tuples is written as a list which is appended
    to/extended and then "magicked" into a tuple by the new LIST_TO_TUPLE
    opcode.

    This peephole repeatedly analyses the bytecode in a block looking for a
    window between a `LIST_TO_TUPLE` and `BUILD_LIST` and...

    1. Turns the BUILD_LIST into a BUILD_TUPLE
    2. Sets an accumulator's initial value as the target of the BUILD_TUPLE
    3. Searches for 'extend' on the original list and turns these into binary
       additions on the accumulator.
    4. Searches for 'append' on the original list and turns these into a
       `BUILD_TUPLE` which is then appended via binary addition to the
       accumulator.
    5. Assigns the accumulator to the variable that exits the peephole and the
       rest of the block/code refers to as the result of the unpack operation.
    6. Patches up
    """
    _DEBUG = False

    # For all blocks
    for offset, blk in func_ir.blocks.items():
        # keep doing the peephole rewrite until nothing is left that matches
        while True:
            # first try and find a matching region
            # i.e. BUILD_LIST...<stuff>...LIST_TO_TUPLE
            def find_postive_region():
                found = False
                for idx in reversed(range(len(blk.body))):
                    stmt = blk.body[idx]
                    if isinstance(stmt, ir.Assign):
                        value = stmt.value
                        if (isinstance(value, ir.Expr) and
                                value.op == 'list_to_tuple'):
                            target_list = value.info[0]
                            found = True
                            bt = (idx, stmt)
                    if found:
                        if isinstance(stmt, ir.Assign):
                            if stmt.target.name == target_list:
                                region = (bt, (idx, stmt))
                                return region

            region = find_postive_region()
            # if there's a peep hole region then do something with it
            if region is not None:
                peep_hole = blk.body[region[1][0] : region[0][0]]
                if _DEBUG:
                    print("\nWINDOW:")
                    for x in peep_hole:
                        print(x)
                    print("")

                appends = []
                extends = []
                init = region[1][1]
                const_list = init.target.name
                # Walk through the peep_hole and find things that are being
                # "extend"ed and "append"ed to the BUILD_LIST
                for x in peep_hole:
                    if isinstance(x, ir.Assign):
                        if isinstance(x.value, ir.Expr):
                            expr = x.value
                            if (expr.op == 'getattr' and
                                    expr.value.name == const_list):
                                # it's not strictly necessary to split out
                                # extends and appends, but it helps with
                                # debugging to do so!
                                if expr.attr == 'extend':
                                    extends.append(x.target.name)
                                elif expr.attr == 'append':
                                    appends.append(x.target.name)
                                else:
                                    assert 0
                # go back through the peep hole build new IR based on it.
                new_hole = []

                def append_and_fix(x):
                    """ Adds to the new_hole and fixes up definitions"""
                    new_hole.append(x)
                    if x.target.name in func_ir._definitions:
                        # if there's already a definition, drop it, should only
                        # be 1 as the way cpython emits the sequence for
                        # `list_to_tuple` should ensure this.
                        assert len(func_ir._definitions[x.target.name]) == 1
                        func_ir._definitions[x.target.name].clear()
                    func_ir._definitions[x.target.name].append(x.value)

                the_build_list = init.target

                # Do the transform on the peep hole
                if _DEBUG:
                    print("\nBLOCK:")
                    blk.dump()

                # This section basically accumulates list appends and extends
                # as binop(+) on tuples, it drops all the getattr() for extend
                # and append as they are now dead and replaced with binop(+).
                # It also switches out the build_list for a build_tuple and then
                # ensures everything is wired up and defined ok.
                t2l_agn = region[0][1]
                acc = the_build_list
                for x in peep_hole:
                    if isinstance(x, ir.Assign):
                        if isinstance(x.value, ir.Expr):
                            expr = x.value
                            if expr.op == 'getattr':
                                if (x.target.name in extends or
                                        x.target.name in appends):
                                    # drop definition, it's being wholesale
                                    # replaced.
                                    func_ir._definitions.pop(x.target.name)
                                    continue
                                else:
                                    # a getattr on something we're not
                                    # interested in
                                    new_hole.append(x)
                            elif expr.op == 'call':
                                fname = expr.func.name
                                if fname in extends or fname in appends:
                                    arg = expr.args[0]
                                    if isinstance(arg, ir.Var):
                                        tmp_name = "%s_var_%s" % (fname,
                                                                  arg.name)
                                        if fname in appends:
                                            bt = ir.Expr.build_tuple([arg,],
                                                                     expr.loc)
                                        else:
                                            # Extend as tuple
                                            gv_tuple = ir.Global(
                                                name="tuple", value=tuple,
                                                loc=expr.loc,
                                            )
                                            tuple_var = arg.scope.redefine(
                                                "$_list_extend_gv_tuple",
                                                loc=expr.loc,
                                            )
                                            new_hole.append(
                                                ir.Assign(
                                                    target=tuple_var,
                                                    value=gv_tuple,
                                                    loc=expr.loc,
                                                ),
                                            )
                                            bt = ir.Expr.call(
                                                tuple_var, (arg,), (),
                                                loc=expr.loc,
                                            )
                                        var = ir.Var(arg.scope, tmp_name,
                                                     expr.loc)
                                        asgn = ir.Assign(bt, var, expr.loc)
                                        append_and_fix(asgn)
                                        arg = var

                                    # this needs to be a binary add
                                    new = ir.Expr.binop(fn=operator.add,
                                                        lhs=acc,
                                                        rhs=arg,
                                                        loc=x.loc)
                                    asgn = ir.Assign(new, x.target, expr.loc)
                                    append_and_fix(asgn)
                                    acc = asgn.target
                                else:
                                    # there could be a call in the unpack, like
                                    # *(a, x.append(y))
                                    new_hole.append(x)
                            elif (expr.op == 'build_list' and
                                    x.target.name == const_list):
                                new = ir.Expr.build_tuple(expr.items, expr.loc)
                                asgn = ir.Assign(new, x.target, expr.loc)
                                # Not a temporary any more
                                append_and_fix(asgn)
                            else:
                                new_hole.append(x)
                        else:
                            new_hole.append(x)

                    else:
                        # stick everything else in as-is
                        new_hole.append(x)
                # Finally write the result back into the original build list as
                # everything refers to it.
                append_and_fix(ir.Assign(acc, t2l_agn.target,
                                         the_build_list.loc))
                if _DEBUG:
                    print("\nNEW HOLE:")
                    for x in new_hole:
                        print(x)

                # and then update the block body with the modified region
                cpy = blk.body[:]
                head = cpy[:region[1][0]]
                tail = blk.body[region[0][0] + 1:]
                tmp = head + new_hole + tail
                blk.body.clear()
                blk.body.extend(tmp)

                if _DEBUG:
                    print("\nDUMP post hole:")
                    blk.dump()

            else:
                # else escape
                break

    return func_ir


def peep_hole_delete_with_exit(func_ir):
    """
    This rewrite removes variables used to store the `__exit__` function
    loaded by SETUP_WITH.
    """
    dead_vars = set()

    for blk in func_ir.blocks.values():
        for stmt in blk.body:
            # Any statement that uses a variable with the '$setup_with_exitfn'
            # prefix is considered dead.
            used = set(stmt.list_vars())
            for v in used:
                if v.name.startswith('$setup_with_exitfn'):
                    dead_vars.add(v)
            # Any assignment that uses any of the dead variable is considered
            # dead.
            if used & dead_vars:
                if isinstance(stmt, ir.Assign):
                    dead_vars.add(stmt.target)

        new_body = []
        for stmt in blk.body:
            # Skip any statements that uses anyone of the dead variable.
            if not (set(stmt.list_vars()) & dead_vars):
                new_body.append(stmt)
        blk.body.clear()
        blk.body.extend(new_body)

    return func_ir


def peep_hole_fuse_dict_add_updates(func_ir):
    """
    This rewrite removes d1._update_from_bytecode(d2)
    calls that are between two dictionaries, d1 and d2,
    in the same basic block. This pattern can appear as a
    result of Python 3.10 bytecode emission changes, which
    prevent large constant literal dictionaries
    (> 15 elements) from being constant. If both dictionaries
    are constant dictionaries defined in the same block and
    neither is used between the update call, then we replace d1
    with a new definition that combines the two dictionaries. At
    the bytecode translation stage we convert DICT_UPDATE into
    _update_from_bytecode, so we know that _update_from_bytecode
    always comes from the bytecode change and not user code.

    Python 3.10 may also rewrite the individual dictionaries
    as an empty build_map + many map_add. Here we again look
    for an _update_from_bytecode, and if so we replace these
    with a single constant dictionary.

    When running this algorithm we can always safely remove d2.

    This is the relevant section of the CPython 3.10 that causes
    this bytecode change:
    https://github.com/python/cpython/blob/3.10/Python/compile.c#L4048
    """

    # This algorithm fuses build_map expressions into the largest
    # possible build map before use. For example, if we have an
    # IR that looks like this:
    #
    #   $d1 = build_map([])
    #   $key = const("a")
    #   $value = const(2)
    #   $setitem_func = getattr($d1, "__setitem__")
    #   $unused1 = call (setitem_func, ($key, $value))
    #   $key2 = const("b")
    #   $value2 = const(3)
    #   $d2 = build_map([($key2, $value2)])
    #   $update_func = getattr($d1, "_update_from_bytecode")
    #   $unused2 = call ($update_func, ($d2,))
    #   $othervar = None
    #   $retvar = cast($othervar)
    #   return $retvar
    #
    # Then the IR is rewritten such that any __setitem__ and
    # _update_from_bytecode operations are fused into the original buildmap.
    # The new buildmap is then added to the
    # last location where it had previously had encountered a __setitem__,
    # _update_from_bytecode, or build_map before any other uses.
    # The new IR would look like:
    #
    #   $key = const("a")
    #   $value = const(2)
    #   $key2 = const("b")
    #   $value2 = const(3)
    #   $d1 = build_map([($key, $value), ($key2, $value2)])
    #   $othervar = None
    #   $retvar = cast($othervar)
    #   return $retvar
    #
    # Note that we don't push $d1 to the bottom of the block. This is because
    # some values may be found below this block (e.g pop_block) that are pattern
    # matched in other locations, such as objmode handling. It should be safe to
    # move a map to the last location at which there was _update_from_bytecode.

    errmsg = textwrap.dedent("""
        A DICT_UPDATE op-code was encountered that could not be replaced.
        If you have created a large constant dictionary, this may
        be an an indication that you are using inlined control
        flow. You can resolve this issue by moving the control flow out of
        the dicitonary constructor. For example, if you have

            d = {a: 1 if flag else 0, ...)

        Replace that with:

            a_val = 1 if flag else 0
            d = {a: a_val, ...)""")

    already_deleted_defs = collections.defaultdict(set)
    for blk in func_ir.blocks.values():
        new_body = []
        # literal map var name -> block idx of the original build_map
        lit_map_def_idx = {}
        # literal map var name -> list(map_uses)
        # This is the index of every build_map or __setitem__
        # in the IR that will need to be removed if the map
        # is updated.
        lit_map_use_idx = collections.defaultdict(list)
        # literal map var name -> list of key/value items for build map
        map_updates = {}
        blk_changed = False

        for i, stmt in enumerate(blk.body):
            # What instruction should we append
            new_inst = stmt
            # Name that should be skipped when tracking used
            # vars in statement. This is always the lhs with
            # a build_map.
            stmt_build_map_out = None
            if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Expr):
                if stmt.value.op == "build_map":
                    # Skip the output build_map when looking for used vars.
                    stmt_build_map_out = stmt.target.name
                    # If we encounter a build map add it to the
                    # tracked maps.
                    lit_map_def_idx[stmt.target.name] = i
                    lit_map_use_idx[stmt.target.name].append(i)
                    map_updates[stmt.target.name] = stmt.value.items.copy()
                elif stmt.value.op == "call" and i > 0:
                    # If we encounter a call we may need to replace
                    # the body
                    func_name = stmt.value.func.name
                    # If we have an update or a setitem
                    # it will be the previous expression.
                    getattr_stmt = blk.body[i - 1]
                    args = stmt.value.args
                    if (
                        isinstance(getattr_stmt, ir.Assign)
                        and getattr_stmt.target.name == func_name
                        and isinstance(getattr_stmt.value, ir.Expr)
                        and getattr_stmt.value.op == "getattr"
                        and getattr_stmt.value.attr in (
                            "__setitem__", "_update_from_bytecode"
                        )
                    ):
                        update_map_name = getattr_stmt.value.value.name
                        attr = getattr_stmt.value.attr
                        if (attr == "__setitem__"
                           and update_map_name in lit_map_use_idx):
                            # If we have a setitem, update the lists
                            map_updates[update_map_name].append(args)
                            # Update the list of instructions that would
                            # need to be removed to include the setitem
                            # and the the getattr
                            lit_map_use_idx[update_map_name].extend([i - 1, i])
                        elif attr == "_update_from_bytecode":
                            d2_map_name = args[0].name
                            if (update_map_name in lit_map_use_idx
                               and d2_map_name in lit_map_use_idx):
                                # If we have an update and the arg is also
                                # a literal dictionary, fuse the lists.
                                map_updates[update_map_name].extend(
                                    map_updates[d2_map_name]
                                )
                                # Delete the old IR for d1 and d2
                                lit_map_use_idx[update_map_name].extend(
                                    lit_map_use_idx[d2_map_name]
                                )
                                lit_map_use_idx[update_map_name].append(i - 1)
                                for linenum in lit_map_use_idx[update_map_name]:
                                    # Drop the existing definition.
                                    _remove_assignment_definition(
                                        blk.body,
                                        linenum,
                                        func_ir,
                                        already_deleted_defs,
                                    )
                                    # Delete it from the new block
                                    new_body[linenum] = None
                                # Delete the maps from dicts
                                del lit_map_def_idx[d2_map_name]
                                del lit_map_use_idx[d2_map_name]
                                del map_updates[d2_map_name]
                                # Add d1 as the new instruction, removing the
                                # old definition.
                                _remove_assignment_definition(
                                    blk.body, i, func_ir, already_deleted_defs
                                )
                                new_inst = _build_new_build_map(
                                    func_ir,
                                    update_map_name,
                                    blk.body,
                                    lit_map_def_idx[update_map_name],
                                    map_updates[update_map_name],
                                )
                                # Update d1 in lit_map_use_idx to just the new
                                # definition and clear the previous list.
                                lit_map_use_idx[update_map_name].clear()
                                lit_map_use_idx[update_map_name].append(i)
                                # Mark that this block has been modified
                                blk_changed = True
                            else:
                                # If we cannot remove _update_from_bytecode
                                # Then raise an error for the user.
                                raise UnsupportedBytecodeError(errmsg)

            # Check if we need to drop any maps from being tracked.
            # Skip the setitem/_update_from_bytecode getattr that
            # will be removed when handling their call in the next
            # iteration.
            if not (
                isinstance(stmt, ir.Assign)
                and isinstance(stmt.value, ir.Expr)
                and stmt.value.op == "getattr"
                and stmt.value.value.name in lit_map_use_idx
                and stmt.value.attr in ("__setitem__", "_update_from_bytecode")
            ):
                for var in stmt.list_vars():
                    # If a map is used it cannot be fused later in
                    # the block. As a result we delete it from
                    # the dicitonaries
                    if (
                        var.name in lit_map_use_idx
                        and var.name != stmt_build_map_out
                    ):
                        del lit_map_def_idx[var.name]
                        del lit_map_use_idx[var.name]
                        del map_updates[var.name]

            # Append the instruction to the new block
            new_body.append(new_inst)

        if blk_changed:
            # If the block is changed replace the block body.
            blk.body.clear()
            blk.body.extend([x for x in new_body if x is not None])

    return func_ir


def peep_hole_split_at_pop_block(func_ir):
    """
    Split blocks that contain ir.PopBlock.

    This rewrite restores the IR structure to pre 3.11 so that withlifting
    can work correctly.
    """
    new_block_map = {}
    sorted_blocks = sorted(func_ir.blocks.items())
    for blk_idx, (label, blk) in enumerate(sorted_blocks):
        # Gather locations of PopBlock
        pop_block_locs = []
        for i, inst in enumerate(blk.body):
            if isinstance(inst, ir.PopBlock):
                pop_block_locs.append(i)
        # Rewrite block with PopBlock
        if pop_block_locs:
            new_blocks = []
            for i in pop_block_locs:
                before_blk = ir.Block(blk.scope, loc=blk.loc)
                before_blk.body.extend(blk.body[:i])
                new_blocks.append(before_blk)

                popblk_blk = ir.Block(blk.scope, loc=blk.loc)
                popblk_blk.body.append(blk.body[i])
                new_blocks.append(popblk_blk)
            # Add jump instructions
            prev_label = label
            for newblk in new_blocks:
                new_block_map[prev_label] = newblk
                next_label = prev_label + 1
                newblk.body.append(ir.Jump(next_label, loc=blk.loc))
                prev_label = next_label
            # Check prev_label does not exceed current new block label
            if blk_idx + 1 < len(sorted_blocks):
                if prev_label >= sorted_blocks[blk_idx + 1][0]:
                    # Panic! Due to heuristic in with-lifting, block labels
                    # must be monotonically increasing. We cannot continue if we
                    # run out of usable label between the two blocks.
                    raise errors.InternalError("POP_BLOCK peephole failed")
            # Add tail block, which will get the original terminator
            tail_blk = ir.Block(blk.scope, loc=blk.loc)
            tail_blk.body.extend(blk.body[pop_block_locs[-1] + 1:])
            new_block_map[prev_label] = tail_blk

    func_ir.blocks.update(new_block_map)
    return func_ir


def _build_new_build_map(func_ir, name, old_body, old_lineno, new_items):
    """
    Create a new build_map with a new set of key/value items
    but all the other info the same.
    """
    old_assign = old_body[old_lineno]
    old_target = old_assign.target
    old_bm = old_assign.value
    # Build the literals
    literal_keys = []
    # Track the constant key/values to set the literal_value
    # field of build_map properly
    values = []
    for pair in new_items:
        k, v = pair
        key_def = guard(get_definition, func_ir, k)
        if isinstance(key_def, (ir.Const, ir.Global, ir.FreeVar)):
            literal_keys.append(key_def.value)
        value_def = guard(get_definition, func_ir, v)
        if isinstance(value_def, (ir.Const, ir.Global, ir.FreeVar)):
            values.append(value_def.value)
        else:
            # Append unknown value if not a literal.
            values.append(_UNKNOWN_VALUE(v.name))

    value_indexes = {}
    if len(literal_keys) == len(new_items):
        # All keys must be literals to have any literal values.
        literal_value = {x: y for x, y in zip(literal_keys, values)}
        for i, k in enumerate(literal_keys):
            value_indexes[k] = i
    else:
        literal_value = None

    # Construct a new build map.
    new_bm = ir.Expr.build_map(
        items=new_items,
        size=len(new_items),
        literal_value=literal_value,
        value_indexes=value_indexes,
        loc=old_bm.loc,
    )

    # The previous definition has already been removed
    # when updating the IR in peep_hole_fuse_dict_add_updates
    func_ir._definitions[name].append(new_bm)

    # Return a new assign.
    return ir.Assign(
        new_bm, ir.Var(old_target.scope, name, old_target.loc), new_bm.loc
    )


class Interpreter(object):
    """A bytecode interpreter that builds up the IR.
    """

    _DEBUG_PRINT = False

    def __init__(self, func_id):
        self.func_id = func_id
        if self._DEBUG_PRINT:
            print(func_id.func)
        self.arg_count = func_id.arg_count
        self.arg_names = func_id.arg_names
        self.loc = self.first_loc = ir.Loc.from_function_id(func_id)
        self.is_generator = func_id.is_generator

        # { inst offset : ir.Block }
        self.blocks = {}
        # { name: [definitions] } of local variables
        self.definitions = collections.defaultdict(list)
        # A set to keep track of all exception variables.
        # To be used in _legalize_exception_vars()
        self._exception_vars = set()

    def interpret(self, bytecode):
        """
        Generate IR for this bytecode.
        """
        self.bytecode = bytecode

        self.scopes = []
        global_scope = ir.Scope(parent=None, loc=self.loc)
        self.scopes.append(global_scope)

        flow = Flow(bytecode)
        flow.run()
        self.dfa = AdaptDFA(flow)
        self.cfa = AdaptCFA(flow)
        if config.DUMP_CFG:
            self.cfa.dump()

        # Temp states during interpretation
        self.current_block = None
        self.current_block_offset = None
        last_active_offset = 0
        for _, inst_blocks in self.cfa.blocks.items():
            if inst_blocks.body:
                last_active_offset = max(last_active_offset,
                                         max(inst_blocks.body))
        self.last_active_offset = last_active_offset

        if PYVERSION in ((3, 12), (3, 13)):
            self.active_exception_entries = tuple(
                [entry for entry in self.bytecode.exception_entries
                 if entry.start < self.last_active_offset])
        elif PYVERSION in ((3, 10), (3, 11)):
            pass
        else:
            raise NotImplementedError(PYVERSION)
        self.syntax_blocks = []
        self.dfainfo = None

        self.scopes.append(ir.Scope(parent=self.current_scope, loc=self.loc))

        # Interpret loop
        for inst, kws in self._iter_inst():
            self._dispatch(inst, kws)
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            # Insert end of try markers
            self._end_try_blocks()
        elif PYVERSION in ((3, 10),):
            pass
        else:
            raise NotImplementedError(PYVERSION)
        self._legalize_exception_vars()
        # Prepare FunctionIR
        func_ir = ir.FunctionIR(self.blocks, self.is_generator, self.func_id,
                                self.first_loc, self.definitions,
                                self.arg_count, self.arg_names)
        _logger.debug(_lazy_pformat(func_ir,
                                    lazy_func=lambda x: x.dump_to_string()))

        # post process the IR to rewrite opcodes/byte sequences that are too
        # involved to risk handling as part of direct interpretation
        peepholes = []
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            peepholes.append(peep_hole_split_at_pop_block)
        if PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            peepholes.append(peep_hole_list_to_tuple)
        peepholes.append(peep_hole_delete_with_exit)
        if PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            # peep_hole_call_function_ex_to_call_function_kw
            # depends on peep_hole_list_to_tuple converting
            # any large number of arguments from a list to a
            # tuple.
            peepholes.append(peep_hole_call_function_ex_to_call_function_kw)
            peepholes.append(peep_hole_fuse_dict_add_updates)

        post_processed_ir = self.post_process(peepholes, func_ir)

        return post_processed_ir

    def post_process(self, peepholes, func_ir):
        for peep in peepholes:
            func_ir = peep(func_ir)
        return func_ir

    def _end_try_blocks(self):
        """Closes all try blocks by inserting the required marker at the
        exception handler

        This is only needed for py3.11 because of the changes in exception
        handling. This merely maps the new py3.11 semantics back to the old way.

        What the code does:

        - For each block, compute the difference of blockstack to its incoming
          blocks' blockstack.
        - If the incoming blockstack has an extra TRY, the current block must
          be the EXCEPT block and we need to insert a marker.

        See also: _insert_try_block_end
        """
        assert PYVERSION in ((3, 11), (3, 12), (3, 13))
        graph = self.cfa.graph
        for offset, block in self.blocks.items():
            # Get current blockstack
            cur_bs = self.dfa.infos[offset].blockstack
            # Check blockstack of the incoming blocks
            for inc, _ in graph.predecessors(offset):
                inc_bs = self.dfa.infos[inc].blockstack

                # find first diff in the blockstack
                for i, (x, y) in enumerate(zip(cur_bs, inc_bs)):
                    if x != y:
                        break
                else:
                    i = min(len(cur_bs), len(inc_bs))

                def do_change(remain):
                    while remain:
                        ent = remain.pop()
                        if ent['kind'] == BlockKind('TRY'):
                            # Extend block with marker for end of try
                            self.current_block = block
                            oldbody = list(block.body)
                            block.body.clear()
                            self._insert_try_block_end()
                            block.body.extend(oldbody)
                            return True

                if do_change(list(inc_bs[i:])):
                    break

    def _legalize_exception_vars(self):
        """Search for unsupported use of exception variables.
        Note, they cannot be stored into user variable.
        """
        # Build a set of exception variables
        excvars = self._exception_vars.copy()
        # Propagate the exception variables to LHS of assignment
        for varname, defnvars in self.definitions.items():
            for v in defnvars:
                if isinstance(v, ir.Var):
                    k = v.name
                    if k in excvars:
                        excvars.add(varname)
        # Filter out the user variables.
        uservar = list(filter(lambda x: not x.startswith('$'), excvars))
        if uservar:
            # Complain about the first user-variable storing an exception
            first = uservar[0]
            loc = self.current_scope.get(first).loc
            msg = "Exception object cannot be stored into variable ({})."
            raise errors.UnsupportedBytecodeError(msg.format(first), loc=loc)

    def init_first_block(self):
        # Define variables receiving the function arguments
        for index, name in enumerate(self.arg_names):
            val = ir.Arg(index=index, name=name, loc=self.loc)
            self.store(val, name)

    def _iter_inst(self):
        for blkct, block in enumerate(self.cfa.iterliveblocks()):
            firstinst = self.bytecode[block.offset]
            # If its an END_FOR instruction, the start location of block
            # is set to start of the FOR loop, so take the location of
            # next instruction. This only affects the source location
            # marking and has no impact to semantic.
            if firstinst.opname == 'END_FOR':
                firstinst = self.bytecode[firstinst.next]
            self.loc = self.loc.with_lineno(firstinst.lineno)
            self._start_new_block(block.offset)
            if blkct == 0:
                # Is first block
                self.init_first_block()
            for offset, kws in self.dfainfo.insts:
                inst = self.bytecode[offset]
                self.loc = self.loc.with_lineno(inst.lineno)
                yield inst, kws
            self._end_current_block()

    def _start_new_block(self, offset):
        oldblock = self.current_block
        self.insert_block(offset)

        tryblk = self.dfainfo.active_try_block if self.dfainfo else None
        # Ensure the last block is terminated
        if oldblock is not None and not oldblock.is_terminated:
            # Handle ending try block.
            # If there's an active try-block and the handler block is live.
            if tryblk is not None and tryblk['end'] in self.cfa.graph.nodes():
                # We are in a try-block, insert a branch to except-block.
                # This logic cannot be in self._end_current_block()
                # because we don't know the non-raising next block-offset.
                branch = ir.Branch(
                    cond=self.get('$exception_check'),
                    truebr=tryblk['end'],
                    falsebr=offset,
                    loc=self.loc,
                )
                oldblock.append(branch)
            # Handle normal case
            else:
                jmp = ir.Jump(offset, loc=self.loc)
                oldblock.append(jmp)

        # Get DFA block info
        self.dfainfo = self.dfa.infos[self.current_block_offset]
        self.assigner = Assigner()
        # Check out-of-scope syntactic-block
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            # This is recreating pre-3.11 code structure
            while self.syntax_blocks:
                if offset >= self.syntax_blocks[-1].exit:
                    synblk = self.syntax_blocks.pop()
                    if isinstance(synblk, ir.With):
                        self.current_block.append(ir.PopBlock(self.loc))
                else:
                    break
            # inject try block:
            newtryblk = self.dfainfo.active_try_block
            if newtryblk is not None:
                if newtryblk is not tryblk:
                    self._insert_try_block_begin()
        elif PYVERSION in ((3, 10),):
            while self.syntax_blocks:
                if offset >= self.syntax_blocks[-1].exit:
                    self.syntax_blocks.pop()
                else:
                    break
        else:
            raise NotImplementedError(PYVERSION)

    def _end_current_block(self):
        # Handle try block
        if not self.current_block.is_terminated:
            tryblk = self.dfainfo.active_try_block
            if tryblk is not None:
                self._insert_exception_check()
        # Handle normal block cleanup
        self._remove_unused_temporaries()
        self._insert_outgoing_phis()

    def _inject_call(self, func, gv_name, res_name=None):
        """A helper function to inject a call to *func* which is a python
        function.
        Parameters
        ----------
        func : callable
            The function object to be called.
        gv_name : str
            The variable name to be used to store the function object.
        res_name : str; optional
            The variable name to be used to store the call result.
            If ``None``, a name is created automatically.
        """
        gv_fn = ir.Global(gv_name, func, loc=self.loc)
        self.store(value=gv_fn, name=gv_name, redefine=True)
        callres = ir.Expr.call(self.get(gv_name), (), (), loc=self.loc)
        res_name = res_name or '$callres_{}'.format(gv_name)
        self.store(value=callres, name=res_name, redefine=True)

    def _insert_try_block_begin(self):
        """Insert IR-nodes to mark the start of a `try` block.
        """
        self._inject_call(eh.mark_try_block, 'mark_try_block')

    def _insert_try_block_end(self):
        """Insert IR-nodes to mark the end of a `try` block.
        """
        self._inject_call(eh.end_try_block, 'end_try_block')

    def _insert_exception_variables(self):
        """Insert IR-nodes to initialize the exception variables.
        """
        tryblk = self.dfainfo.active_try_block
        # Get exception variables
        endblk = tryblk['end']
        edgepushed = self.dfainfo.outgoing_edgepushed.get(endblk)
        # Note: the last value on the stack is the exception value
        # Note: due to the current limitation, all exception variables are None
        if edgepushed:
            const_none = ir.Const(value=None, loc=self.loc)
            # For each variable going to the handler block.
            for var in edgepushed:
                if var in self.definitions:
                    raise AssertionError(
                        "exception variable CANNOT be defined by other code",
                    )
                self.store(value=const_none, name=var)
                self._exception_vars.add(var)

    def _insert_exception_check(self):
        """Called before the end of a block to inject checks if raised.
        """
        self._insert_exception_variables()
        # Do exception check
        self._inject_call(eh.exception_check, 'exception_check',
                          '$exception_check')

    def _remove_unused_temporaries(self):
        """
        Remove assignments to unused temporary variables from the
        current block.
        """
        new_body = []
        replaced_var = {}
        for inst in self.current_block.body:
            # the same temporary is assigned to multiple variables in cases
            # like a = b[i] = 1, so need to handle replaced temporaries in
            # later setitem/setattr nodes
            if (isinstance(inst, (ir.SetItem, ir.SetAttr))
                    and inst.value.name in replaced_var):
                inst.value = replaced_var[inst.value.name]
            elif isinstance(inst, ir.Assign):
                if (inst.target.is_temp
                        and inst.target.name in self.assigner.unused_dests):
                    continue
                # the same temporary is assigned to multiple variables in cases
                # like a = b = 1, so need to handle replaced temporaries in
                # later assignments
                if (isinstance(inst.value, ir.Var)
                        and inst.value.name in replaced_var):
                    inst.value = replaced_var[inst.value.name]
                    new_body.append(inst)
                    continue
                # chained unpack cases may reuse temporary
                # e.g. a = (b, c) = (x, y)
                if (isinstance(inst.value, ir.Expr)
                        and inst.value.op == "exhaust_iter"
                        and inst.value.value.name in replaced_var):
                    inst.value.value = replaced_var[inst.value.value.name]
                    new_body.append(inst)
                    continue
                # eliminate temporary variables that are assigned to user
                # variables right after creation. E.g.:
                # $1 = f(); a = $1 -> a = f()
                # the temporary variable is not reused elsewhere since CPython
                # bytecode is stack-based and this pattern corresponds to a pop
                if (isinstance(inst.value, ir.Var) and inst.value.is_temp
                        and new_body and isinstance(new_body[-1], ir.Assign)):
                    prev_assign = new_body[-1]
                    # _var_used_in_binop check makes sure we don't create a new
                    # inplace binop operation which can fail
                    # (see TestFunctionType.test_in_iter_func_call)
                    if (prev_assign.target.name == inst.value.name
                            and not self._var_used_in_binop(
                                inst.target.name, prev_assign.value)):
                        replaced_var[inst.value.name] = inst.target
                        prev_assign.target = inst.target
                        # replace temp var definition in target with proper defs
                        self.definitions[inst.target.name].remove(inst.value)
                        self.definitions[inst.target.name].extend(
                            self.definitions.pop(inst.value.name)
                        )
                        continue

            new_body.append(inst)

        self.current_block.body = new_body

    def _var_used_in_binop(self, varname, expr):
        """return True if 'expr' is a binary expression and 'varname' is used
        in it as an argument
        """
        return (isinstance(expr, ir.Expr)
                and expr.op in ("binop", "inplace_binop")
                and (varname == expr.lhs.name or varname == expr.rhs.name))

    def _insert_outgoing_phis(self):
        """
        Add assignments to forward requested outgoing values
        to subsequent blocks.
        """
        for phiname, varname in self.dfainfo.outgoing_phis.items():
            target = self.current_scope.get_or_define(phiname,
                                                      loc=self.loc)
            try:
                val = self.get(varname)
            except ir.NotDefinedError:
                # Hack to make sure exception variables are defined
                assert PYVERSION in ((3, 11), (3, 12), (3, 13)), \
                       "unexpected missing definition"
                val = ir.Const(value=None, loc=self.loc)
            stmt = ir.Assign(value=val, target=target,
                             loc=self.loc)
            self.definitions[target.name].append(stmt.value)
            if not self.current_block.is_terminated:
                self.current_block.append(stmt)
            else:
                self.current_block.insert_before_terminator(stmt)

    def get_global_value(self, name):
        """
        Get a global value from the func_global (first) or
        as a builtins (second).  If both failed, return a ir.UNDEFINED.
        """
        try:
            return self.func_id.func.__globals__[name]
        except KeyError:
            return getattr(builtins, name, ir.UNDEFINED)

    def get_closure_value(self, index):
        """
        Get a value from the cell contained in this function's closure.
        If not set, return a ir.UNDEFINED.
        """
        cell = self.func_id.func.__closure__[index]
        try:
            return cell.cell_contents
        except ValueError:
            return ir.UNDEFINED

    @property
    def current_scope(self):
        return self.scopes[-1]

    @property
    def code_consts(self):
        return self.bytecode.co_consts

    @property
    def code_locals(self):
        return self.bytecode.co_varnames

    @property
    def code_names(self):
        return self.bytecode.co_names

    @property
    def code_cellvars(self):
        return self.bytecode.co_cellvars

    @property
    def code_freevars(self):
        return self.bytecode.co_freevars

    def _dispatch(self, inst, kws):
        if self._DEBUG_PRINT:
            print(inst)
        assert self.current_block is not None
        if PYVERSION in ((3, 11), (3, 12), (3, 13)):
            if self.syntax_blocks:
                top = self.syntax_blocks[-1]
                if isinstance(top, ir.With) :
                    if inst.offset >= top.exit:
                        self.current_block.append(ir.PopBlock(loc=self.loc))
                        self.syntax_blocks.pop()
        elif PYVERSION in ((3, 10),):
            pass
        else:
            raise NotImplementedError(PYVERSION)

        fname = "op_%s" % inst.opname.replace('+', '_')
        try:
            fn = getattr(self, fname)
        except AttributeError:
            raise NotImplementedError(inst)
        else:
            try:
                return fn(inst, **kws)
            except errors.NotDefinedError as e:
                if e.loc is None:
                    loc = self.loc
                else:
                    loc = e.loc

                err = errors.NotDefinedError(e.name, loc=loc)
                if not config.FULL_TRACEBACKS:
                    raise err from None
                else:
                    m = f"handling op: {inst} | offset: {inst.offset}"
                    err.add_context(m)
                    err.add_context(self.bytecode.dump())
                    raise err

    # --- Scope operations ---

    def store(self, value, name, redefine=False):
        """
        Store *value* (a Expr or Var instance) into the variable named *name*
        (a str object). Returns the target variable.
        """
        if redefine or self.current_block_offset in self.cfa.backbone:
            rename = not (name in self.code_cellvars)
            target = self.current_scope.redefine(name, loc=self.loc,
                                                 rename=rename)
        else:
            target = self.current_scope.get_or_define(name, loc=self.loc)
        if isinstance(value, ir.Var):
            value = self.assigner.assign(value, target)
        stmt = ir.Assign(value=value, target=target, loc=self.loc)
        self.current_block.append(stmt)
        self.definitions[target.name].append(value)
        return target

    def get(self, name):
        """
        Get the variable (a Var instance) with the given *name*.
        """
        # Implicit argument for comprehension starts with '.'
        # See Parameter class in inspect.py (from Python source)
        if name[0] == '.' and name[1:].isdigit():
            name = 'implicit{}'.format(name[1:])

        # Try to simplify the variable lookup by returning an earlier
        # variable assigned to *name*.
        var = self.assigner.get_assignment_source(name)
        if var is None:
            var = self.current_scope.get(name)
        return var

    # --- Block operations ---

    def insert_block(self, offset, scope=None, loc=None):
        scope = scope or self.current_scope
        loc = loc or self.loc
        blk = ir.Block(scope=scope, loc=loc)
        self.blocks[offset] = blk
        self.current_block = blk
        self.current_block_offset = offset
        return blk

    # --- Bytecode handlers ---

    def op_NOP(self, inst):
        pass

    def op_RESUME(self, inst):
        pass

    def op_CACHE(self, inst):
        pass

    def op_PRECALL(self, inst):
        pass

    def op_PUSH_NULL(self, inst):
        pass

    def op_RETURN_GENERATOR(self, inst):
        pass

    def op_PRINT_ITEM(self, inst, item, printvar, res):
        item = self.get(item)
        printgv = ir.Global("print", print, loc=self.loc)
        self.store(value=printgv, name=printvar)
        call = ir.Expr.call(self.get(printvar), (item,), (), loc=self.loc)
        self.store(value=call, name=res)

    def op_PRINT_NEWLINE(self, inst, printvar, res):
        printgv = ir.Global("print", print, loc=self.loc)
        self.store(value=printgv, name=printvar)
        call = ir.Expr.call(self.get(printvar), (), (), loc=self.loc)
        self.store(value=call, name=res)

    def op_UNPACK_SEQUENCE(self, inst, iterable, stores, tupleobj):
        count = len(stores)
        # Exhaust the iterable into a tuple-like object
        tup = ir.Expr.exhaust_iter(value=self.get(iterable), loc=self.loc,
                                   count=count)
        self.store(name=tupleobj, value=tup)

        # then index the tuple-like object to extract the values
        for i, st in enumerate(stores):
            expr = ir.Expr.static_getitem(self.get(tupleobj),
                                          index=i, index_var=None,
                                          loc=self.loc)
            self.store(expr, st)

    def op_FORMAT_SIMPLE(self, inst, value, res, strvar):
        # Same as FORMAT_VALUE
        return self.op_FORMAT_VALUE(inst, value, res, strvar)

    def op_FORMAT_VALUE(self, inst, value, res, strvar):
        """
        FORMAT_VALUE(flags): flags argument specifies format spec which is not
        supported yet. Currently, str() is simply called on the value.
        https://docs.python.org/3/library/dis.html#opcode-FORMAT_VALUE
        """
        value = self.get(value)
        strgv = ir.Global("str", str, loc=self.loc)
        self.store(value=strgv, name=strvar)
        call = ir.Expr.call(self.get(strvar), (value,), (), loc=self.loc)
        self.store(value=call, name=res)

    def op_BUILD_STRING(self, inst, strings, tmps):
        """
        BUILD_STRING(count): Concatenates count strings.
        Required for supporting f-strings.
        https://docs.python.org/3/library/dis.html#opcode-BUILD_STRING
        """
        count = inst.arg
        # corner case: f""
        if count == 0:
            const = ir.Const("", loc=self.loc)
            self.store(const, tmps[-1])
            return

        prev = self.get(strings[0])
        for other, tmp in zip(strings[1:], tmps):
            other = self.get(other)
            expr = ir.Expr.binop(
                operator.add, lhs=prev, rhs=other, loc=self.loc
            )
            self.store(expr, tmp)
            prev = self.get(tmp)

    def op_BUILD_SLICE(self, inst, start, stop, step, res, slicevar):
        start = self.get(start)
        stop = self.get(stop)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        if step is None:
            sliceinst = ir.Expr.call(self.get(slicevar), (start, stop), (),
                                     loc=self.loc)
        else:
            step = self.get(step)
            sliceinst = ir.Expr.call(self.get(slicevar), (start, stop, step),
                                     (), loc=self.loc)
        self.store(value=sliceinst, name=res)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_BINARY_SLICE(self, inst, start, end, container, res, slicevar,
                            temp_res):
            start = self.get(start)
            end = self.get(end)
            slicegv = ir.Global("slice", slice, loc=self.loc)
            self.store(value=slicegv, name=slicevar)
            sliceinst = ir.Expr.call(self.get(slicevar), (start, end), (),
                                     loc=self.loc)
            self.store(value=sliceinst, name=temp_res)
            index = self.get(temp_res)
            target = self.get(container)
            expr = ir.Expr.getitem(target, index=index, loc=self.loc)
            self.store(expr, res)
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_STORE_SLICE(self, inst, start, end, container, value, res,
                           slicevar):
            start = self.get(start)
            end = self.get(end)
            slicegv = ir.Global("slice", slice, loc=self.loc)
            self.store(value=slicegv, name=slicevar)
            sliceinst = ir.Expr.call(self.get(slicevar), (start, end), (),
                                     loc=self.loc)
            self.store(value=sliceinst, name=res)
            index = self.get(res)
            target = self.get(container)
            value = self.get(value)

            stmt = ir.SetItem(target=target, index=index, value=value,
                              loc=self.loc)
            self.current_block.append(stmt)
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_SLICE_0(self, inst, base, res, slicevar, indexvar, nonevar):
        base = self.get(base)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        index = ir.Expr.call(self.get(slicevar), (none, none), (), loc=self.loc)
        self.store(value=index, name=indexvar)

        expr = ir.Expr.getitem(base, self.get(indexvar), loc=self.loc)
        self.store(value=expr, name=res)

    def op_SLICE_1(self, inst, base, start, nonevar, res, slicevar, indexvar):
        base = self.get(base)
        start = self.get(start)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (start, none), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        expr = ir.Expr.getitem(base, self.get(indexvar), loc=self.loc)
        self.store(value=expr, name=res)

    def op_SLICE_2(self, inst, base, nonevar, stop, res, slicevar, indexvar):
        base = self.get(base)
        stop = self.get(stop)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (none, stop,), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        expr = ir.Expr.getitem(base, self.get(indexvar), loc=self.loc)
        self.store(value=expr, name=res)

    def op_SLICE_3(self, inst, base, start, stop, res, slicevar, indexvar):
        base = self.get(base)
        start = self.get(start)
        stop = self.get(stop)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (start, stop), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        expr = ir.Expr.getitem(base, self.get(indexvar), loc=self.loc)
        self.store(value=expr, name=res)

    def op_STORE_SLICE_0(self, inst, base, value, slicevar, indexvar, nonevar):
        base = self.get(base)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        index = ir.Expr.call(self.get(slicevar), (none, none), (), loc=self.loc)
        self.store(value=index, name=indexvar)

        stmt = ir.SetItem(base, self.get(indexvar), self.get(value),
                          loc=self.loc)
        self.current_block.append(stmt)

    def op_STORE_SLICE_1(self, inst, base, start, nonevar, value, slicevar,
                         indexvar):
        base = self.get(base)
        start = self.get(start)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (start, none), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        stmt = ir.SetItem(base, self.get(indexvar), self.get(value),
                          loc=self.loc)
        self.current_block.append(stmt)

    def op_STORE_SLICE_2(self, inst, base, nonevar, stop, value, slicevar,
                         indexvar):
        base = self.get(base)
        stop = self.get(stop)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (none, stop,), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        stmt = ir.SetItem(base, self.get(indexvar), self.get(value),
                          loc=self.loc)
        self.current_block.append(stmt)

    def op_STORE_SLICE_3(self, inst, base, start, stop, value, slicevar,
                         indexvar):
        base = self.get(base)
        start = self.get(start)
        stop = self.get(stop)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (start, stop), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)
        stmt = ir.SetItem(base, self.get(indexvar), self.get(value),
                          loc=self.loc)
        self.current_block.append(stmt)

    def op_DELETE_SLICE_0(self, inst, base, slicevar, indexvar, nonevar):
        base = self.get(base)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        index = ir.Expr.call(self.get(slicevar), (none, none), (), loc=self.loc)
        self.store(value=index, name=indexvar)

        stmt = ir.DelItem(base, self.get(indexvar), loc=self.loc)
        self.current_block.append(stmt)

    def op_DELETE_SLICE_1(self, inst, base, start, nonevar, slicevar, indexvar):
        base = self.get(base)
        start = self.get(start)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (start, none), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        stmt = ir.DelItem(base, self.get(indexvar), loc=self.loc)
        self.current_block.append(stmt)

    def op_DELETE_SLICE_2(self, inst, base, nonevar, stop, slicevar, indexvar):
        base = self.get(base)
        stop = self.get(stop)

        nonegv = ir.Const(None, loc=self.loc)
        self.store(value=nonegv, name=nonevar)
        none = self.get(nonevar)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (none, stop,), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)

        stmt = ir.DelItem(base, self.get(indexvar), loc=self.loc)
        self.current_block.append(stmt)

    def op_DELETE_SLICE_3(self, inst, base, start, stop, slicevar, indexvar):
        base = self.get(base)
        start = self.get(start)
        stop = self.get(stop)

        slicegv = ir.Global("slice", slice, loc=self.loc)
        self.store(value=slicegv, name=slicevar)

        index = ir.Expr.call(self.get(slicevar), (start, stop), (),
                             loc=self.loc)
        self.store(value=index, name=indexvar)
        stmt = ir.DelItem(base, self.get(indexvar), loc=self.loc)
        self.current_block.append(stmt)

    def _op_LOAD_FAST(self, inst, res):
        srcname = self.code_locals[inst.arg]
        self.store(value=self.get(srcname), name=res)

    if PYVERSION in ((3, 13), ):
        def op_LOAD_FAST(self, inst, res, as_load_deref=False):
            if as_load_deref:
                self.op_LOAD_DEREF(inst, res)
            else:
                self._op_LOAD_FAST(inst, res)

    else:
        op_LOAD_FAST = _op_LOAD_FAST

    if PYVERSION in ((3, 13),):
        def op_LOAD_FAST_LOAD_FAST(self, inst, res1, res2):
            oparg = inst.arg
            oparg1 = oparg >> 4
            oparg2 = oparg & 15
            src1 = self.get(self.code_locals[oparg1])
            src2 = self.get(self.code_locals[oparg2])
            self.store(value=src1, name=res1)
            self.store(value=src2, name=res2)

        def op_STORE_FAST_LOAD_FAST(self, inst, store_value, load_res):
            oparg = inst.arg
            oparg1 = oparg >> 4
            oparg2 = oparg & 15

            dstname = self.code_locals[oparg1]
            dst_value = self.get(store_value)
            self.store(value=dst_value, name=dstname)

            src_value = self.get(self.code_locals[oparg2])
            self.store(value=src_value, name=load_res)

        def op_STORE_FAST_STORE_FAST(self, inst, value1, value2):
            oparg = inst.arg
            oparg1 = oparg >> 4
            oparg2 = oparg & 15

            dstname = self.code_locals[oparg1]
            self.store(value=self.get(value1), name=dstname)
            dstname = self.code_locals[oparg2]
            self.store(value=self.get(value2), name=dstname)

    elif PYVERSION in ((3, 10), (3, 11), (3, 12)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 12), (3, 13)):
        op_LOAD_FAST_CHECK = op_LOAD_FAST

        def op_LOAD_FAST_AND_CLEAR(self, inst, res):
            try:
                # try the regular LOAD_FAST logic
                srcname = self.code_locals[inst.arg]
                self.store(value=self.get(srcname), name=res)
            except NotDefinedError:
                # If the variable is not in the scope, set it to `undef`
                undef = ir.Expr.undef(loc=self.loc)
                self.store(undef, name=res)

    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_STORE_FAST(self, inst, value):
        dstname = self.code_locals[inst.arg]
        value = self.get(value)
        self.store(value=value, name=dstname)

    def op_DELETE_FAST(self, inst):
        dstname = self.code_locals[inst.arg]
        self.current_block.append(ir.Del(dstname, loc=self.loc))

    def op_DUP_TOPX(self, inst, orig, duped):
        for src, dst in zip(orig, duped):
            self.store(value=self.get(src), name=dst)

    op_DUP_TOP = op_DUP_TOPX
    op_DUP_TOP_TWO = op_DUP_TOPX

    def op_STORE_ATTR(self, inst, target, value):
        attr = self.code_names[inst.arg]
        sa = ir.SetAttr(target=self.get(target), value=self.get(value),
                        attr=attr, loc=self.loc)
        self.current_block.append(sa)

    def op_DELETE_ATTR(self, inst, target):
        attr = self.code_names[inst.arg]
        sa = ir.DelAttr(target=self.get(target), attr=attr, loc=self.loc)
        self.current_block.append(sa)

    def op_LOAD_ATTR(self, inst, item, res):
        item = self.get(item)
        if PYVERSION in ((3, 12), (3, 13)):
            attr = self.code_names[inst.arg >> 1]
        elif PYVERSION in ((3, 10), (3, 11)):
            attr = self.code_names[inst.arg]
        else:
            raise NotImplementedError(PYVERSION)
        getattr = ir.Expr.getattr(item, attr, loc=self.loc)
        self.store(getattr, res)

    def op_LOAD_CONST(self, inst, res):
        value = self.code_consts[inst.arg]
        if isinstance(value, tuple):
            st = []
            for x in value:
                nm = '$const_%s' % str(x)
                val_const = ir.Const(x, loc=self.loc)
                target = self.store(val_const, name=nm, redefine=True)
                st.append(target)
            const = ir.Expr.build_tuple(st, loc=self.loc)
        elif isinstance(value, frozenset):
            st = []
            for x in value:
                nm = '$const_%s' % str(x)
                val_const = ir.Const(x, loc=self.loc)
                target = self.store(val_const, name=nm, redefine=True)
                st.append(target)
            const = ir.Expr.build_set(st, loc=self.loc)
        else:
            const = ir.Const(value, loc=self.loc)
        self.store(const, res)

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_LOAD_GLOBAL(self, inst, idx, res):
            name = self.code_names[idx]
            value = self.get_global_value(name)
            gl = ir.Global(name, value, loc=self.loc)
            self.store(gl, res)
    elif PYVERSION in ((3, 10),):
        def op_LOAD_GLOBAL(self, inst, res):
            name = self.code_names[inst.arg]
            value = self.get_global_value(name)
            gl = ir.Global(name, value, loc=self.loc)
            self.store(gl, res)
    else:
        raise NotImplementedError(PYVERSION)

    def op_COPY_FREE_VARS(self, inst):
        pass

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_LOAD_DEREF(self, inst, res):
            name = self.func_id.func.__code__._varname_from_oparg(inst.arg)
            if name in self.code_cellvars:
                try:
                    gl = self.get(name)
                except NotDefinedError:
                    msg = "Unsupported use of cell variable encountered"
                    raise NotImplementedError(msg)
            elif name in self.code_freevars:
                idx = self.code_freevars.index(name)
                value = self.get_closure_value(idx)
                gl = ir.FreeVar(idx, name, value, loc=self.loc)
            self.store(gl, res)
    elif PYVERSION in ((3, 10),):
        def op_LOAD_DEREF(self, inst, res):
            n_cellvars = len(self.code_cellvars)
            if inst.arg < n_cellvars:
                name = self.code_cellvars[inst.arg]
                gl = self.get(name)
            else:
                idx = inst.arg - n_cellvars
                name = self.code_freevars[idx]
                value = self.get_closure_value(idx)
                gl = ir.FreeVar(idx, name, value, loc=self.loc)
            self.store(gl, res)
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_MAKE_CELL(self, inst):
            pass  # ignored bytecode

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_STORE_DEREF(self, inst, value):
            name = self.func_id.func.__code__._varname_from_oparg(inst.arg)
            value = self.get(value)
            self.store(value=value, name=name)
    elif PYVERSION in ((3, 10),):
        def op_STORE_DEREF(self, inst, value):
            n_cellvars = len(self.code_cellvars)
            if inst.arg < n_cellvars:
                dstname = self.code_cellvars[inst.arg]
            else:
                dstname = self.code_freevars[inst.arg - n_cellvars]
            value = self.get(value)
            self.store(value=value, name=dstname)
    else:
        raise NotImplementedError(PYVERSION)

    def op_SETUP_LOOP(self, inst):
        assert self.blocks[inst.offset] is self.current_block
        loop = ir.Loop(inst.offset, exit=(inst.next + inst.arg))
        self.syntax_blocks.append(loop)

    def op_SETUP_WITH(self, inst, contextmanager, exitfn=None):
        assert self.blocks[inst.offset] is self.current_block
        # Handle with
        exitpt = inst.next + inst.arg

        wth = ir.With(inst.offset, exit=exitpt)
        self.syntax_blocks.append(wth)
        ctxmgr = self.get(contextmanager)
        self.current_block.append(ir.EnterWith(contextmanager=ctxmgr,
                                               begin=inst.offset,
                                               end=exitpt, loc=self.loc,))

        # Store exit fn
        exit_fn_obj = ir.Const(None, loc=self.loc)
        self.store(value=exit_fn_obj, name=exitfn)

    def op_BEFORE_WITH(self, inst, contextmanager, exitfn, end):
        assert self.blocks[inst.offset] is self.current_block
        if PYVERSION in ((3, 12), (3, 13)):
            # Python 3.12 hack for handling nested with blocks
            if end > self.last_active_offset:
                # Use exception entries to figure out end of syntax block
                end = max([ex.end for ex in self.active_exception_entries
                           if ex.target == end])
        elif PYVERSION in ((3, 10), (3, 11)):
            pass
        else:
            raise NotImplementedError(PYVERSION)
        # Handle with
        wth = ir.With(inst.offset, exit=end)
        self.syntax_blocks.append(wth)
        ctxmgr = self.get(contextmanager)
        self.current_block.append(ir.EnterWith(contextmanager=ctxmgr,
                                               begin=inst.offset,
                                               end=end, loc=self.loc,))

        # Store exit function
        exit_fn_obj = ir.Const(None, loc=self.loc)
        self.store(value=exit_fn_obj, name=exitfn)

    def op_SETUP_FINALLY(self, inst):
        # Removed since python3.11
        self._insert_try_block_begin()

    def op_WITH_CLEANUP(self, inst):
        "no-op"

    def op_WITH_CLEANUP_START(self, inst):
        "no-op"

    def op_WITH_CLEANUP_FINISH(self, inst):
        "no-op"

    def op_END_FINALLY(self, inst):
        "no-op"

    def op_BEGIN_FINALLY(self, inst, temps):
        # The *temps* are the exception variables
        const_none = ir.Const(None, loc=self.loc)
        for tmp in temps:
            # Set to None for now
            self.store(const_none, name=tmp)
            self._exception_vars.add(tmp)

    def op_CALL(self, inst, func, args, kw_names, res):
        func = self.get(func)
        args = [self.get(x) for x in args]
        if kw_names is not None:
            assert PYVERSION < (3, 13)
            names = self.code_consts[kw_names]
            kwargs = list(zip(names, args[-len(names):]))
            args = args[:-len(names)]
        else:
            kwargs = ()
        expr = ir.Expr.call(func, args, kwargs, loc=self.loc)
        self.store(expr, res)

    if PYVERSION in ((3, 13),):
        def op_CALL_KW(self, inst, func, args, kw_names, res):
            func = self.get(func)
            args = [self.get(x) for x in args]
            consti = int(kw_names.rsplit('.', 2)[-1])
            names = self.code_consts[consti]
            kwargs = list(zip(names, args[-len(names):]))
            args = args[:-len(names)]
            expr = ir.Expr.call(func, args, kwargs, loc=self.loc)
            self.store(expr, res)
    else:
        assert PYVERSION < (3, 13)

    def op_CALL_FUNCTION(self, inst, func, args, res):
        func = self.get(func)
        args = [self.get(x) for x in args]
        expr = ir.Expr.call(func, args, (), loc=self.loc)
        self.store(expr, res)

    def op_CALL_FUNCTION_KW(self, inst, func, args, names, res):
        func = self.get(func)
        args = [self.get(x) for x in args]
        # Find names const
        names = self.get(names)
        for inst in self.current_block.body:
            if isinstance(inst, ir.Assign) and inst.target is names:
                self.current_block.remove(inst)
                # scan up the block looking for the values, remove them
                # and find their name strings
                named_items = []
                for x in inst.value.items:
                    for y in self.current_block.body[::-1]:
                        if x == y.target:
                            self.current_block.remove(y)
                            named_items.append(y.value.value)
                            break
                keys = named_items
                break

        nkeys = len(keys)
        posvals = args[:-nkeys]
        kwvals = args[-nkeys:]
        keyvalues = list(zip(keys, kwvals))

        expr = ir.Expr.call(func, posvals, keyvalues, loc=self.loc)
        self.store(expr, res)

    def op_CALL_FUNCTION_EX(self, inst, func, vararg, varkwarg, res):
        func = self.get(func)
        vararg = self.get(vararg)
        if varkwarg is not None:
            varkwarg = self.get(varkwarg)
        expr = ir.Expr.call(
            func, [], [], loc=self.loc, vararg=vararg, varkwarg=varkwarg
        )
        self.store(expr, res)

    def _build_tuple_unpack(self, inst, tuples, temps, is_assign):
        first = self.get(tuples[0])
        if is_assign:
            # it's assign-like, defer handling to an intrinsic that will have
            # type information.
            # Can deal with tuples only, i.e. y = (*x,). where x = <tuple>
            gv_name = "unpack_single_tuple"
            gv_fn = ir.Global(gv_name, unpack_single_tuple, loc=self.loc,)
            self.store(value=gv_fn, name=gv_name, redefine=True)
            exc = ir.Expr.call(self.get(gv_name), args=(first,), kws=(),
                               loc=self.loc,)
            self.store(exc, temps[0])
        else:
            loc = self.loc
            for other, tmp in zip(map(self.get, tuples[1:]), temps):
                # Emit as `first + tuple(other)`
                gv_tuple = ir.Global(
                    name="tuple", value=tuple,
                    loc=loc,
                )
                tuple_var = self.store(
                    gv_tuple, "$_list_extend_gv_tuple", redefine=True,
                )
                tuplify_val = ir.Expr.call(
                    tuple_var, (other,), (),
                    loc=loc,
                )
                tuplify_var = self.store(tuplify_val, "$_tuplify",
                                         redefine=True)
                out = ir.Expr.binop(
                    fn=operator.add, lhs=first, rhs=self.get(tuplify_var.name),
                    loc=self.loc,
                )
                self.store(out, tmp)
                first = self.get(tmp)

    def op_BUILD_TUPLE_UNPACK_WITH_CALL(self, inst, tuples, temps, is_assign):
        # just unpack the input tuple, call inst will be handled afterwards
        self._build_tuple_unpack(inst, tuples, temps, is_assign)

    def op_BUILD_TUPLE_UNPACK(self, inst, tuples, temps, is_assign):
        self._build_tuple_unpack(inst, tuples, temps, is_assign)

    def op_LIST_TO_TUPLE(self, inst, const_list, res):
        expr = ir.Expr.dummy('list_to_tuple', (const_list,), loc=self.loc)
        self.store(expr, res)

    def op_BUILD_CONST_KEY_MAP(self, inst, keys, keytmps, values, res):
        # Unpack the constant key-tuple and reused build_map which takes
        # a sequence of (key, value) pair.
        keyvar = self.get(keys)
        # TODO: refactor this pattern. occurred several times.
        for inst in self.current_block.body:
            if isinstance(inst, ir.Assign) and inst.target is keyvar:
                self.current_block.remove(inst)
                # scan up the block looking for the values, remove them
                # and find their name strings
                named_items = []
                for x in inst.value.items:
                    for y in self.current_block.body[::-1]:
                        if x == y.target:
                            self.current_block.remove(y)
                            named_items.append(y.value.value)
                            break
                keytup = named_items
                break
        assert len(keytup) == len(values)
        keyconsts = [ir.Const(value=x, loc=self.loc) for x in keytup]
        for kval, tmp in zip(keyconsts, keytmps):
            self.store(kval, tmp)
        items = list(zip(map(self.get, keytmps), map(self.get, values)))

        # sort out literal values
        literal_items = []
        for v in values:
            defns = self.definitions[v]
            if len(defns) != 1:
                break
            defn = defns[0]
            if not isinstance(defn, ir.Const):
                break
            literal_items.append(defn.value)

        def resolve_const(v):
            defns = self.definitions[v]
            if len(defns) != 1:
                return _UNKNOWN_VALUE(self.get(v).name)
            defn = defns[0]
            if not isinstance(defn, ir.Const):
                return _UNKNOWN_VALUE(self.get(v).name)
            return defn.value

        if len(literal_items) != len(values):
            literal_dict = {x: resolve_const(y) for x, y in
                            zip(keytup, values)}
        else:
            literal_dict = {x:y for x, y in zip(keytup, literal_items)}

        # to deal with things like {'a': 1, 'a': 'cat', 'b': 2, 'a': 2j}
        # store the index of the actual used value for a given key, this is
        # used when lowering to pull the right value out into the tuple repr
        # of a mixed value type dictionary.
        value_indexes = {}
        for i, k in enumerate(keytup):
            value_indexes[k] = i

        expr = ir.Expr.build_map(items=items,
                                 size=2,
                                 literal_value=literal_dict,
                                 value_indexes=value_indexes,
                                 loc=self.loc)

        self.store(expr, res)

    def op_GET_ITER(self, inst, value, res):
        expr = ir.Expr.getiter(value=self.get(value), loc=self.loc)
        self.store(expr, res)

    def op_FOR_ITER(self, inst, iterator, pair, indval, pred):
        """
        Assign new block other this instruction.
        """
        assert inst.offset in self.blocks, "FOR_ITER must be block head"

        # Emit code
        val = self.get(iterator)

        pairval = ir.Expr.iternext(value=val, loc=self.loc)
        self.store(pairval, pair)

        iternext = ir.Expr.pair_first(value=self.get(pair), loc=self.loc)
        self.store(iternext, indval)

        isvalid = ir.Expr.pair_second(value=self.get(pair), loc=self.loc)
        self.store(isvalid, pred)

        # Conditional jump
        br = ir.Branch(cond=self.get(pred), truebr=inst.next,
                       falsebr=inst.get_jump_target(),
                       loc=self.loc)
        self.current_block.append(br)

    def op_BINARY_SUBSCR(self, inst, target, index, res):
        index = self.get(index)
        target = self.get(target)
        expr = ir.Expr.getitem(target, index=index, loc=self.loc)
        self.store(expr, res)

    def op_STORE_SUBSCR(self, inst, target, index, value):
        index = self.get(index)
        target = self.get(target)
        value = self.get(value)
        stmt = ir.SetItem(target=target, index=index, value=value,
                          loc=self.loc)
        self.current_block.append(stmt)

    def op_DELETE_SUBSCR(self, inst, target, index):
        index = self.get(index)
        target = self.get(target)
        stmt = ir.DelItem(target=target, index=index, loc=self.loc)
        self.current_block.append(stmt)

    def op_BUILD_TUPLE(self, inst, items, res):
        expr = ir.Expr.build_tuple(items=[self.get(x) for x in items],
                                   loc=self.loc)
        self.store(expr, res)

    def op_BUILD_LIST(self, inst, items, res):
        expr = ir.Expr.build_list(items=[self.get(x) for x in items],
                                  loc=self.loc)
        self.store(expr, res)

    def op_BUILD_SET(self, inst, items, res):
        expr = ir.Expr.build_set(items=[self.get(x) for x in items],
                                 loc=self.loc)
        self.store(expr, res)

    def op_SET_UPDATE(self, inst, target, value, updatevar, res):
        target = self.get(target)
        value = self.get(value)
        updateattr = ir.Expr.getattr(target, 'update', loc=self.loc)
        self.store(value=updateattr, name=updatevar)
        updateinst = ir.Expr.call(self.get(updatevar), (value,), (),
                                  loc=self.loc)
        self.store(value=updateinst, name=res)

    def op_DICT_UPDATE(self, inst, target, value, updatevar, res):
        target = self.get(target)
        value = self.get(value)
        # We generate _update_from_bytecode instead of update so we can
        # differentiate between user .update() calls and those from the
        # bytecode. This is then used to recombine dictionaries in peephole
        # optimizations. See the dicussion in this PR about why:
        # https://github.com/numba/numba/pull/7964/files#r868229306
        updateattr = ir.Expr.getattr(
            target, '_update_from_bytecode', loc=self.loc
        )
        self.store(value=updateattr, name=updatevar)
        updateinst = ir.Expr.call(self.get(updatevar), (value,), (),
                                  loc=self.loc)
        self.store(value=updateinst, name=res)

    def op_BUILD_MAP(self, inst, items, size, res):
        got_items = [(self.get(k), self.get(v)) for k, v in items]

        # sort out literal values, this is a bit contrived but is to handle
        # situations like `{1: 10, 1: 10}` where the size of the literal dict
        # is smaller than the definition
        def get_literals(target):
            literal_items = []
            values = [self.get(v.name) for v in target]
            for v in values:
                defns = self.definitions[v.name]
                if len(defns) != 1:
                    break
                defn = defns[0]
                if not isinstance(defn, ir.Const):
                    break
                literal_items.append(defn.value)
            return literal_items

        literal_keys = get_literals(x[0] for x in got_items)
        literal_values = get_literals(x[1] for x in got_items)

        has_literal_keys = len(literal_keys) == len(got_items)
        has_literal_values = len(literal_values) == len(got_items)

        value_indexes = {}
        if not has_literal_keys and not has_literal_values:
            literal_dict = None
        elif has_literal_keys and not has_literal_values:
            literal_dict = {x: _UNKNOWN_VALUE(y[1]) for x, y in
                            zip(literal_keys, got_items)}
            for i, k in enumerate(literal_keys):
                value_indexes[k] = i
        else:
            literal_dict = {x: y for x, y in zip(literal_keys, literal_values)}
            for i, k in enumerate(literal_keys):
                value_indexes[k] = i

        expr = ir.Expr.build_map(items=got_items, size=size,
                                 literal_value=literal_dict,
                                 value_indexes=value_indexes,
                                 loc=self.loc)
        self.store(expr, res)

    def op_STORE_MAP(self, inst, dct, key, value):
        stmt = ir.StoreMap(dct=self.get(dct), key=self.get(key),
                           value=self.get(value), loc=self.loc)
        self.current_block.append(stmt)

    def op_UNARY_NEGATIVE(self, inst, value, res):
        value = self.get(value)
        expr = ir.Expr.unary('-', value=value, loc=self.loc)
        return self.store(expr, res)

    def op_UNARY_POSITIVE(self, inst, value, res):
        value = self.get(value)
        expr = ir.Expr.unary('+', value=value, loc=self.loc)
        return self.store(expr, res)

    def op_UNARY_INVERT(self, inst, value, res):
        value = self.get(value)
        expr = ir.Expr.unary('~', value=value, loc=self.loc)
        return self.store(expr, res)

    def op_UNARY_NOT(self, inst, value, res):
        value = self.get(value)
        expr = ir.Expr.unary('not', value=value, loc=self.loc)
        return self.store(expr, res)

    def _binop(self, op, lhs, rhs, res):
        op = BINOPS_TO_OPERATORS[op]
        lhs = self.get(lhs)
        rhs = self.get(rhs)
        expr = ir.Expr.binop(op, lhs=lhs, rhs=rhs, loc=self.loc)
        self.store(expr, res)

    def _inplace_binop(self, op, lhs, rhs, res):
        immuop = BINOPS_TO_OPERATORS[op]
        op = INPLACE_BINOPS_TO_OPERATORS[op + '=']
        lhs = self.get(lhs)
        rhs = self.get(rhs)
        expr = ir.Expr.inplace_binop(op, immuop, lhs=lhs, rhs=rhs,
                                     loc=self.loc)
        self.store(expr, res)

    def op_BINARY_OP(self, inst, op, lhs, rhs, res):
        if "=" in op:
            self._inplace_binop(op[:-1], lhs, rhs, res)
        else:
            self._binop(op, lhs, rhs, res)

    def op_BINARY_ADD(self, inst, lhs, rhs, res):
        self._binop('+', lhs, rhs, res)

    def op_BINARY_SUBTRACT(self, inst, lhs, rhs, res):
        self._binop('-', lhs, rhs, res)

    def op_BINARY_MULTIPLY(self, inst, lhs, rhs, res):
        self._binop('*', lhs, rhs, res)

    def op_BINARY_DIVIDE(self, inst, lhs, rhs, res):
        self._binop('/?', lhs, rhs, res)

    def op_BINARY_TRUE_DIVIDE(self, inst, lhs, rhs, res):
        self._binop('/', lhs, rhs, res)

    def op_BINARY_FLOOR_DIVIDE(self, inst, lhs, rhs, res):
        self._binop('//', lhs, rhs, res)

    def op_BINARY_MODULO(self, inst, lhs, rhs, res):
        self._binop('%', lhs, rhs, res)

    def op_BINARY_POWER(self, inst, lhs, rhs, res):
        self._binop('**', lhs, rhs, res)

    def op_BINARY_MATRIX_MULTIPLY(self, inst, lhs, rhs, res):
        self._binop('@', lhs, rhs, res)

    def op_BINARY_LSHIFT(self, inst, lhs, rhs, res):
        self._binop('<<', lhs, rhs, res)

    def op_BINARY_RSHIFT(self, inst, lhs, rhs, res):
        self._binop('>>', lhs, rhs, res)

    def op_BINARY_AND(self, inst, lhs, rhs, res):
        self._binop('&', lhs, rhs, res)

    def op_BINARY_OR(self, inst, lhs, rhs, res):
        self._binop('|', lhs, rhs, res)

    def op_BINARY_XOR(self, inst, lhs, rhs, res):
        self._binop('^', lhs, rhs, res)

    def op_INPLACE_ADD(self, inst, lhs, rhs, res):
        self._inplace_binop('+', lhs, rhs, res)

    def op_INPLACE_SUBTRACT(self, inst, lhs, rhs, res):
        self._inplace_binop('-', lhs, rhs, res)

    def op_INPLACE_MULTIPLY(self, inst, lhs, rhs, res):
        self._inplace_binop('*', lhs, rhs, res)

    def op_INPLACE_DIVIDE(self, inst, lhs, rhs, res):
        self._inplace_binop('/?', lhs, rhs, res)

    def op_INPLACE_TRUE_DIVIDE(self, inst, lhs, rhs, res):
        self._inplace_binop('/', lhs, rhs, res)

    def op_INPLACE_FLOOR_DIVIDE(self, inst, lhs, rhs, res):
        self._inplace_binop('//', lhs, rhs, res)

    def op_INPLACE_MODULO(self, inst, lhs, rhs, res):
        self._inplace_binop('%', lhs, rhs, res)

    def op_INPLACE_POWER(self, inst, lhs, rhs, res):
        self._inplace_binop('**', lhs, rhs, res)

    def op_INPLACE_MATRIX_MULTIPLY(self, inst, lhs, rhs, res):
        self._inplace_binop('@', lhs, rhs, res)

    def op_INPLACE_LSHIFT(self, inst, lhs, rhs, res):
        self._inplace_binop('<<', lhs, rhs, res)

    def op_INPLACE_RSHIFT(self, inst, lhs, rhs, res):
        self._inplace_binop('>>', lhs, rhs, res)

    def op_INPLACE_AND(self, inst, lhs, rhs, res):
        self._inplace_binop('&', lhs, rhs, res)

    def op_INPLACE_OR(self, inst, lhs, rhs, res):
        self._inplace_binop('|', lhs, rhs, res)

    def op_INPLACE_XOR(self, inst, lhs, rhs, res):
        self._inplace_binop('^', lhs, rhs, res)

    def op_JUMP_ABSOLUTE(self, inst):
        jmp = ir.Jump(inst.get_jump_target(), loc=self.loc)
        self.current_block.append(jmp)

    def op_JUMP_FORWARD(self, inst):
        jmp = ir.Jump(inst.get_jump_target(), loc=self.loc)
        self.current_block.append(jmp)

    def op_JUMP_BACKWARD(self, inst):
        jmp = ir.Jump(inst.get_jump_target(), loc=self.loc)
        self.current_block.append(jmp)

    op_JUMP_BACKWARD_NO_INTERRUPT = op_JUMP_BACKWARD

    def op_POP_BLOCK(self, inst, kind=None):
        if kind is None:
            self.syntax_blocks.pop()
        elif kind == 'with':
            d = ir.PopBlock(loc=self.loc)
            self.current_block.append(d)
        elif kind == 'try':
            self._insert_try_block_end()

    def op_RETURN_VALUE(self, inst, retval, castval):
        self.store(ir.Expr.cast(self.get(retval), loc=self.loc), castval)
        ret = ir.Return(self.get(castval), loc=self.loc)
        self.current_block.append(ret)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_RETURN_CONST(self, inst, retval, castval):
            value = self.code_consts[inst.arg]
            const = ir.Const(value, loc=self.loc)
            self.store(const, retval)
            self.store(ir.Expr.cast(self.get(retval), loc=self.loc), castval)
            ret = ir.Return(self.get(castval), loc=self.loc)
            self.current_block.append(ret)
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    if PYVERSION in ((3, 13),):
        def op_TO_BOOL(self, inst, val, res):
            self.store(self.get(val), res) # TODO: just a lazy hack

    elif PYVERSION in ((3, 10), (3, 11), (3, 12)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_COMPARE_OP(self, inst, lhs, rhs, res):
        if PYVERSION in ((3, 13),):
            op = dis.cmp_op[inst.arg >> 5]
            # TODO: fifth lowest bit now indicates a forced version to bool.
        elif PYVERSION in ((3, 12),):
            op = dis.cmp_op[inst.arg >> 4]
        elif PYVERSION in ((3, 10), (3, 11)):
            op = dis.cmp_op[inst.arg]
        else:
            raise NotImplementedError(PYVERSION)
        if op == 'in' or op == 'not in':
            lhs, rhs = rhs, lhs

        if op == 'not in':
            self._binop('in', lhs, rhs, res)
            tmp = self.get(res)
            out = ir.Expr.unary('not', value=tmp, loc=self.loc)
            self.store(out, res)
        elif op == 'exception match':
            gv_fn = ir.Global(
                "exception_match", eh.exception_match, loc=self.loc,
            )
            exc_match_name = '$exc_match'
            self.store(value=gv_fn, name=exc_match_name, redefine=True)
            lhs = self.get(lhs)
            rhs = self.get(rhs)
            exc = ir.Expr.call(
                self.get(exc_match_name), args=(lhs, rhs), kws=(), loc=self.loc,
            )
            self.store(exc, res)
        else:
            self._binop(op, lhs, rhs, res)

    def op_IS_OP(self, inst, lhs, rhs, res):
        # invert if op case is 1
        op = 'is not' if inst.arg == 1 else 'is'
        self._binop(op, lhs, rhs, res)

    def op_CONTAINS_OP(self, inst, lhs, rhs, res):
        lhs, rhs = rhs, lhs
        self._binop('in', lhs, rhs, res)
        # invert if op case is 1
        if inst.arg == 1:
            tmp = self.get(res)
            out = ir.Expr.unary('not', value=tmp, loc=self.loc)
            self.store(out, res)

    def op_BREAK_LOOP(self, inst, end=None):
        if end is None:
            loop = self.syntax_blocks[-1]
            assert isinstance(loop, ir.Loop)
            end = loop.exit
        jmp = ir.Jump(target=end, loc=self.loc)
        self.current_block.append(jmp)

    def _op_JUMP_IF(self, inst, pred, iftrue):
        brs = {
            True: inst.get_jump_target(),
            False: inst.next,
        }
        truebr = brs[iftrue]
        falsebr = brs[not iftrue]

        name = "$bool%s" % (inst.offset)
        gv_fn = ir.Global("bool", bool, loc=self.loc)
        self.store(value=gv_fn, name=name)

        callres = ir.Expr.call(self.get(name), (self.get(pred),), (),
                               loc=self.loc)

        pname = "$%spred" % (inst.offset)
        predicate = self.store(value=callres, name=pname)
        bra = ir.Branch(cond=predicate, truebr=truebr, falsebr=falsebr,
                        loc=self.loc)
        self.current_block.append(bra)

    def op_JUMP_IF_FALSE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=False)

    def op_JUMP_IF_TRUE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=True)

    def _jump_if_none(self, inst, pred, iftrue):
        # branch pruning assumes true falls through and false is jump
        truebr = inst.next
        falsebr = inst.get_jump_target()

        # this seems strange
        if not iftrue:
            op = BINOPS_TO_OPERATORS["is"]
        else:
            op = BINOPS_TO_OPERATORS["is not"]

        rhs = self.store(value=ir.Const(None, loc=self.loc),
                         name=f"$constNone{inst.offset}")
        lhs = self.get(pred)
        isnone = ir.Expr.binop(op, lhs=lhs, rhs=rhs, loc=self.loc)

        maybeNone = f"$maybeNone{inst.offset}"
        self.store(value=isnone, name=maybeNone)

        name = f"$bool{inst.offset}"
        gv_fn = ir.Global("bool", bool, loc=self.loc)
        self.store(value=gv_fn, name=name)

        callres = ir.Expr.call(self.get(name), (self.get(maybeNone),), (),
                               loc=self.loc)

        pname = f"$pred{inst.offset}"
        predicate = self.store(value=callres, name=pname)
        branch = ir.Branch(cond=predicate,
                           truebr=truebr,
                           falsebr=falsebr,
                           loc=self.loc)
        self.current_block.append(branch)

    def op_POP_JUMP_FORWARD_IF_NONE(self, inst, pred):
        self._jump_if_none(inst, pred, True)

    def op_POP_JUMP_FORWARD_IF_NOT_NONE(self, inst, pred):
        self._jump_if_none(inst, pred, False)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_POP_JUMP_IF_NONE(self, inst, pred):
            self._jump_if_none(inst, pred, True)

        def op_POP_JUMP_IF_NOT_NONE(self, inst, pred):
            self._jump_if_none(inst, pred, False)
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)

    def op_POP_JUMP_BACKWARD_IF_NONE(self, inst, pred):
        self._jump_if_none(inst, pred, True)

    def op_POP_JUMP_BACKWARD_IF_NOT_NONE(self, inst, pred):
        self._jump_if_none(inst, pred, False)

    def op_POP_JUMP_FORWARD_IF_FALSE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=False)

    def op_POP_JUMP_FORWARD_IF_TRUE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=True)

    def op_POP_JUMP_BACKWARD_IF_FALSE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=False)

    def op_POP_JUMP_BACKWARD_IF_TRUE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=True)

    def op_POP_JUMP_IF_FALSE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=False)

    def op_POP_JUMP_IF_TRUE(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=True)

    def op_JUMP_IF_FALSE_OR_POP(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=False)

    def op_JUMP_IF_TRUE_OR_POP(self, inst, pred):
        self._op_JUMP_IF(inst, pred=pred, iftrue=True)

    def op_CHECK_EXC_MATCH(self, inst, pred, tos, tos1):
        gv_fn = ir.Global(
            "exception_match", eh.exception_match, loc=self.loc,
        )
        exc_match_name = '$exc_match'
        self.store(value=gv_fn, name=exc_match_name, redefine=True)
        lhs = self.get(tos1)
        rhs = self.get(tos)
        exc = ir.Expr.call(
            self.get(exc_match_name), args=(lhs, rhs), kws=(), loc=self.loc,
        )
        self.store(exc, pred)

    def op_JUMP_IF_NOT_EXC_MATCH(self, inst, pred, tos, tos1):
        truebr = inst.next
        falsebr = inst.get_jump_target()
        gv_fn = ir.Global(
            "exception_match", eh.exception_match, loc=self.loc,
        )
        exc_match_name = '$exc_match'
        self.store(value=gv_fn, name=exc_match_name, redefine=True)
        lhs = self.get(tos1)
        rhs = self.get(tos)
        exc = ir.Expr.call(
            self.get(exc_match_name), args=(lhs, rhs), kws=(), loc=self.loc,
        )
        predicate = self.store(exc, pred)
        bra = ir.Branch(cond=predicate, truebr=truebr, falsebr=falsebr,
                        loc=self.loc)
        self.current_block.append(bra)

    def op_RERAISE(self, inst, exc):
        tryblk = self.dfainfo.active_try_block
        if tryblk is not None:
            stmt = ir.TryRaise(exception=None, loc=self.loc)
            self.current_block.append(stmt)
            self._insert_try_block_end()
            self.current_block.append(ir.Jump(tryblk['end'], loc=self.loc))
        else:
            # Numba can't handle this case and it's caught else where, this is a
            # runtime guard in case this is reached by unknown means.
            msg = (f"Unreachable condition reached (op code RERAISE executed)"
                   f"{error_extras['reportable']}")
            stmt = ir.StaticRaise(AssertionError, (msg,), self.loc)
            self.current_block.append(stmt)

    def op_RAISE_VARARGS(self, inst, exc):
        if exc is not None:
            exc = self.get(exc)
        tryblk = self.dfainfo.active_try_block
        if tryblk is not None:
            # In a try block
            stmt = ir.TryRaise(exception=exc, loc=self.loc)
            self.current_block.append(stmt)
            self._insert_try_block_end()
            self.current_block.append(ir.Jump(tryblk['end'], loc=self.loc))
        else:
            # Not in a try block
            stmt = ir.Raise(exception=exc, loc=self.loc)
            self.current_block.append(stmt)

    def op_YIELD_VALUE(self, inst, value, res):
        # initialize index to None.  it's being set later in post-processing
        index = None
        inst = ir.Yield(value=self.get(value), index=index, loc=self.loc)
        return self.store(inst, res)

    def op_MAKE_FUNCTION(self, inst, name, code, closure, annotations,
                         kwdefaults, defaults, res):
        # annotations are ignored by numba but useful for static analysis
        # re. https://github.com/numba/numba/issues/7269
        if kwdefaults is not None:
            msg = "op_MAKE_FUNCTION with kwdefaults is not implemented"
            raise NotImplementedError(msg)
        if defaults:
            if isinstance(defaults, tuple):
                defaults = tuple([self.get(name) for name in defaults])
            else:
                defaults = self.get(defaults)

        assume_code_const = self.definitions[code][0]
        if not isinstance(assume_code_const, ir.Const):
            msg = (
                "Unsupported use of closure. "
                "Probably caused by complex control-flow constructs; "
                "e.g. try-except"
            )
            raise errors.UnsupportedBytecodeError(msg, loc=self.loc)
        fcode = assume_code_const.value
        if name:
            name = self.get(name)
        if closure:
            closure = self.get(closure)
        expr = ir.Expr.make_function(name, fcode, closure, defaults, self.loc)
        self.store(expr, res)

    def op_MAKE_CLOSURE(self, inst, name, code, closure, annotations,
                        kwdefaults, defaults, res):
        self.op_MAKE_FUNCTION(inst, name, code, closure, annotations,
                              kwdefaults, defaults, res)

    if PYVERSION in ((3, 11), (3, 12), (3, 13)):
        def op_LOAD_CLOSURE(self, inst, res):
            name = self.func_id.func.__code__._varname_from_oparg(inst.arg)
            if name in self.code_cellvars:
                try:
                    gl = self.get(name)
                except NotDefinedError:
                    msg = "Unsupported use of cell variable encountered"
                    raise NotImplementedError(msg)
            elif name in self.code_freevars:
                idx = self.code_freevars.index(name)
                value = self.get_closure_value(idx)
                gl = ir.FreeVar(idx, name, value, loc=self.loc)
            else:
                assert 0, "unreachable"
            self.store(gl, res)

    elif PYVERSION in ((3, 10),):
        def op_LOAD_CLOSURE(self, inst, res):
            n_cellvars = len(self.code_cellvars)
            if inst.arg < n_cellvars:
                name = self.code_cellvars[inst.arg]
                try:
                    gl = self.get(name)
                except NotDefinedError:
                    msg = "Unsupported use of cell variable encountered"
                    raise NotImplementedError(msg)
            else:
                idx = inst.arg - n_cellvars
                name = self.code_freevars[idx]
                value = self.get_closure_value(idx)
                gl = ir.FreeVar(idx, name, value, loc=self.loc)
            self.store(gl, res)
    else:
        raise NotImplementedError(PYVERSION)

    def op_LIST_APPEND(self, inst, target, value, appendvar, res):
        target = self.get(target)
        value = self.get(value)
        appendattr = ir.Expr.getattr(target, 'append', loc=self.loc)
        self.store(value=appendattr, name=appendvar)
        appendinst = ir.Expr.call(self.get(appendvar), (value,), (),
                                  loc=self.loc)
        self.store(value=appendinst, name=res)

    def op_LIST_EXTEND(self, inst, target, value, extendvar, res):
        target = self.get(target)
        value = self.get(value)
        # If the statements between the current instruction and the target
        # are N * consts followed by build_tuple AND the target has no items,
        # it's a situation where a list is being statically initialised, rewrite
        # the build_tuple as a build_list, drop the extend, and wire up the
        # target as the result from the build_tuple that's been rewritten.

        # See if this is the first statement in a block, if so its probably from
        # control flow in a tuple unpack like:
        # `(*(1, (2,) if predicate else (3,)))`
        # this cannot be handled as present so raise
        msg = ("An unsupported bytecode sequence has been encountered: "
               "op_LIST_EXTEND at the start of a block.\n\nThis could be "
               "due to the use of a branch in a tuple unpacking statement.")
        if not self.current_block.body:
            raise errors.UnsupportedBytecodeError(msg)

        # is last emitted statement a build_tuple?
        stmt = self.current_block.body[-1]
        ok = isinstance(stmt.value, ir.Expr) and stmt.value.op == "build_tuple"
        # check statements from self.current_block.body[-1] through to target,
        # make sure they are consts
        build_empty_list = None
        if ok:
            for stmt in reversed(self.current_block.body[:-1]):
                if not isinstance(stmt, ir.Assign):
                    ok = False
                    break
                # if its not a const, it needs to be the `build_list` for the
                # target, else it's something else we don't know about so just
                # bail
                if isinstance(stmt.value, ir.Const):
                    continue

                # it's not a const, check for target
                elif isinstance(stmt.value, ir.Expr) and stmt.target == target:
                    build_empty_list = stmt
                    # it's only ok to do this if the target has no initializer
                    # already
                    ok = not stmt.value.items
                    break
                else:
                    ok = False
                    break
        if ok and build_empty_list is None:
            raise errors.UnsupportedBytecodeError(msg)
        if ok:
            stmts = self.current_block.body
            build_tuple_asgn = self.current_block.body[-1]
            # move build list to last issued statement
            stmts.append(stmts.pop(stmts.index(build_empty_list)))
            # fix the build list
            build_tuple = build_tuple_asgn.value
            build_list = build_empty_list.value
            build_list.items = build_tuple.items
        else:
            # it's just a list extend with no static init, let it be
            extendattr = ir.Expr.getattr(target, 'extend', loc=self.loc)
            self.store(value=extendattr, name=extendvar)
            extendinst = ir.Expr.call(self.get(extendvar), (value,), (),
                                      loc=self.loc)
            self.store(value=extendinst, name=res)

    def op_MAP_ADD(self, inst, target, key, value, setitemvar, res):
        target = self.get(target)
        key = self.get(key)
        value = self.get(value)
        setitemattr = ir.Expr.getattr(target, '__setitem__', loc=self.loc)
        self.store(value=setitemattr, name=setitemvar)
        appendinst = ir.Expr.call(self.get(setitemvar), (key, value,), (),
                                  loc=self.loc)
        self.store(value=appendinst, name=res)

    def op_LOAD_ASSERTION_ERROR(self, inst, res):
        gv_fn = ir.Global("AssertionError", AssertionError, loc=self.loc)
        self.store(value=gv_fn, name=res)

    # NOTE: The LOAD_METHOD opcode is implemented as a LOAD_ATTR for ease,
    # however this means a new object (the bound-method instance) could be
    # created. Conversely, using a pure LOAD_METHOD no intermediary is present
    # and it is essentially like a pointer grab and forward to CALL_METHOD. The
    # net outcome is that the implementation in Numba produces the same result,
    # but in object mode it may be that it runs more slowly than it would if
    # run in CPython.

    def op_LOAD_METHOD(self, *args, **kws):
        self.op_LOAD_ATTR(*args, **kws)

    def op_CALL_METHOD(self, *args, **kws):
        self.op_CALL_FUNCTION(*args, **kws)

    if PYVERSION in ((3, 12), (3, 13)):
        def op_CALL_INTRINSIC_1(self, inst, operand, **kwargs):
            if operand == ci1op.INTRINSIC_STOPITERATION_ERROR:
                stmt = ir.StaticRaise(INTRINSIC_STOPITERATION_ERROR, (),
                                      self.loc)
                self.current_block.append(stmt)
                return
            elif operand == ci1op.UNARY_POSITIVE:
                self.op_UNARY_POSITIVE(inst, **kwargs)
                return
            elif operand == ci1op.INTRINSIC_LIST_TO_TUPLE:
                self.op_LIST_TO_TUPLE(inst, **kwargs)
                return
            else:
                raise NotImplementedError(operand)
    elif PYVERSION in ((3, 10), (3, 11)):
        pass
    else:
        raise NotImplementedError(PYVERSION)


if PYVERSION in ((3, 12), (3, 13)):
    class INTRINSIC_STOPITERATION_ERROR(AssertionError):
        pass
elif PYVERSION in ((3, 10), (3, 11)):
    pass
else:
    raise NotImplementedError(PYVERSION)
