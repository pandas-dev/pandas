import copy
import operator

import types as pytypes
import operator
import warnings
from dataclasses import make_dataclass

import llvmlite.ir
import numpy as np

import numba
from numba.parfors import parfor
from numba.core import types, ir, config, compiler, sigutils, cgutils
from numba.core.ir_utils import (
    add_offset_to_labels,
    replace_var_names,
    remove_dels,
    legalize_names,
    rename_labels,
    get_name_var_table,
    visit_vars_inner,
    get_definition,
    guard,
    get_call_table,
    is_pure,
    get_np_ufunc_typ,
    get_unused_var_name,
    is_const_call,
    fixup_var_define_in_scope,
    transfer_scope,
    find_max_label,
    get_global_func_typ,
    find_topo_order,
)
from numba.core.typing import signature
from numba.core import lowering
from numba.parfors.parfor import ensure_parallel_support
from numba.core.errors import (
    NumbaParallelSafetyWarning, NotDefinedError, CompilerError, InternalError,
)
from numba.parfors.parfor_lowering_utils import ParforLoweringBuilder


class ParforLower(lowering.Lower):
    """This is a custom lowering class that extends standard lowering so as
    to accommodate parfor.Parfor nodes."""

    # custom instruction lowering to handle parfor nodes
    def lower_inst(self, inst):
        if isinstance(inst, parfor.Parfor):
            _lower_parfor_parallel(self, inst)
        else:
            super().lower_inst(inst)

    @property
    def _disable_sroa_like_opt(self):
        """
        Force disable this because Parfor use-defs is incompatible---it only
        considers use-defs in blocks that must be executing.
        See https://github.com/numba/numba/commit/017e2ff9db87fc34149b49dd5367ecbf0bb45268
        """
        return True


def _lower_parfor_parallel(lowerer, parfor):
    if parfor.lowerer is None:
        return _lower_parfor_parallel_std(lowerer, parfor)
    else:
        return parfor.lowerer(lowerer, parfor)


def _lower_parfor_parallel_std(lowerer, parfor):
    """Lowerer that handles LLVM code generation for parfor.
    This function lowers a parfor IR node to LLVM.
    The general approach is as follows:
    1) The code from the parfor's init block is lowered normally
       in the context of the current function.
    2) The body of the parfor is transformed into a gufunc function.
    3) Code is inserted into the main function that calls do_scheduling
       to divide the iteration space for each thread, allocates
       reduction arrays, calls the gufunc function, and then invokes
       the reduction function across the reduction arrays to produce
       the final reduction values.
    """
    from numba.np.ufunc.parallel import get_thread_count

    ensure_parallel_support()
    typingctx = lowerer.context.typing_context
    targetctx = lowerer.context
    builder = lowerer.builder
    # We copy the typemap here because for race condition variable we'll
    # update their type to array so they can be updated by the gufunc.
    orig_typemap = lowerer.fndesc.typemap
    # replace original typemap with copy and restore the original at the end.
    lowerer.fndesc.typemap = copy.copy(orig_typemap)
    if config.DEBUG_ARRAY_OPT:
        print("lowerer.fndesc", lowerer.fndesc, type(lowerer.fndesc))
    typemap = lowerer.fndesc.typemap
    varmap = lowerer.varmap

    if config.DEBUG_ARRAY_OPT:
        print("_lower_parfor_parallel")
        parfor.dump()

    loc = parfor.init_block.loc
    scope = parfor.init_block.scope

    # produce instructions for init_block
    if config.DEBUG_ARRAY_OPT:
        print("init_block = ", parfor.init_block, " ", type(parfor.init_block))
    for instr in parfor.init_block.body:
        if config.DEBUG_ARRAY_OPT:
            print("lower init_block instr = ", instr)
        lowerer.lower_inst(instr)

    for racevar in parfor.races:
        if racevar not in varmap:
            rvtyp = typemap[racevar]
            rv = ir.Var(scope, racevar, loc)
            lowerer._alloca_var(rv.name, rvtyp)

    alias_map = {}
    arg_aliases = {}
    numba.parfors.parfor.find_potential_aliases_parfor(parfor, parfor.params, typemap,
                                        lowerer.func_ir, alias_map, arg_aliases)
    if config.DEBUG_ARRAY_OPT:
        print("alias_map", alias_map)
        print("arg_aliases", arg_aliases)

    # run get_parfor_outputs() and get_parfor_reductions() before gufunc creation
    # since Jumps are modified so CFG of loop_body dict will become invalid
    assert parfor.params is not None

    parfor_output_arrays = numba.parfors.parfor.get_parfor_outputs(
        parfor, parfor.params)
    parfor_redvars, parfor_reddict = parfor.redvars, parfor.reddict
    if config.DEBUG_ARRAY_OPT:
        print("parfor_redvars:", parfor_redvars)
        print("parfor_reddict:", parfor_reddict)

    # init reduction array allocation here.
    nredvars = len(parfor_redvars)
    redarrs = {}
    to_cleanup = []
    if nredvars > 0:
        # reduction arrays outer dimension equal to thread count
        scope = parfor.init_block.scope
        loc = parfor.init_block.loc
        pfbdr = ParforLoweringBuilder(lowerer=lowerer, scope=scope, loc=loc)

        # Get the Numba internal function to call to get the thread count.
        get_num_threads = pfbdr.bind_global_function(
            fobj=numba.np.ufunc.parallel._iget_num_threads,
            ftype=get_global_func_typ(numba.np.ufunc.parallel._iget_num_threads),
            args=()
        )

        # Insert the call to assign the thread count to a variable.
        num_threads_var = pfbdr.assign(
            rhs=pfbdr.call(get_num_threads, args=[]),
            typ=types.intp,
            name="num_threads_var")

        # For each reduction variable...
        for i in range(nredvars):
            red_name = parfor_redvars[i]
            # Get the type of the reduction variable.
            redvar_typ = lowerer.fndesc.typemap[red_name]
            # Get the ir.Var for the reduction variable.
            redvar = ir.Var(scope, red_name, loc)
            # Get the type of the array that holds the per-thread
            # reduction variables.
            redarrvar_typ = redtyp_to_redarraytype(redvar_typ)
            reddtype = redarrvar_typ.dtype
            if config.DEBUG_ARRAY_OPT:
                print(
                    "reduction_info",
                    red_name,
                    redvar_typ,
                    redarrvar_typ,
                    reddtype,
                    types.DType(reddtype),
                    num_threads_var,
                    type(num_threads_var)
                )

            # If this is reduction over an array,
            # the reduction array has just one added per-worker dimension.
            if isinstance(redvar_typ, types.npytypes.Array):
                redarrdim = redvar_typ.ndim + 1
            else:
                redarrdim = 1

            # Reduction array is created and initialized to the initial reduction value.

            # First create a var for the numpy empty ufunc.
            glbl_np_empty = pfbdr.bind_global_function(
                fobj=np.empty,
                ftype=get_np_ufunc_typ(np.empty),
                args=(
                    types.UniTuple(types.intp, redarrdim),
                ),
                kws={'dtype': types.DType(reddtype)}
            )

            size_var_list = [num_threads_var]

            # If this is a reduction over an array...
            if isinstance(redvar_typ, types.npytypes.Array):
                # Add code to get the shape of the array being reduced over.
                redshape_var = pfbdr.assign(
                    rhs=ir.Expr.getattr(redvar, "shape", loc),
                    typ=types.UniTuple(types.intp, redvar_typ.ndim),
                    name="redarr_shape",
                )

                # Add the dimension sizes of the array being reduced over to the tuple of sizes pass to empty.
                for j in range(redvar_typ.ndim):
                    onedimvar = pfbdr.assign(
                        rhs=ir.Expr.static_getitem(redshape_var, j, None, loc),
                        typ=types.intp,
                        name="redshapeonedim",
                    )
                    size_var_list.append(onedimvar)

            # Empty call takes tuple of sizes.  Create here and fill in outer dimension (num threads).
            size_var = pfbdr.make_tuple_variable(
                size_var_list, name='tuple_size_var',
            )

            # Resolve dtype
            cval = pfbdr._typingctx.resolve_value_type(reddtype)
            dt = pfbdr.make_const_variable(cval=cval, typ=types.DType(reddtype))
            # Add call to empty passing the size var tuple.
            empty_call = pfbdr.call(glbl_np_empty, args=[size_var, dt])

            redarr_var = pfbdr.assign(
                rhs=empty_call, typ=redarrvar_typ, name="redarr",
            )

            # Remember mapping of original reduction array to the newly created per-worker reduction array.
            redarrs[redvar.name] = redarr_var
            to_cleanup.append(redarr_var)

            init_val = parfor_reddict[red_name].init_val

            if init_val is not None:
                if isinstance(redvar_typ, types.npytypes.Array):
                    # Create an array of identity values for the reduction.
                    # First, create a variable for np.full.
                    full_func_node = pfbdr.bind_global_function(
                        fobj=np.full,
                        ftype=get_np_ufunc_typ(np.full),
                        args=(
                            types.UniTuple(types.intp, redvar_typ.ndim),
                            reddtype,
                        ),
                        kws={'dtype': types.DType(reddtype)},
                    )

                    # Then create a var with the identify value.
                    init_val_var = pfbdr.make_const_variable(
                        cval=init_val,
                        typ=reddtype,
                        name="init_val",
                    )

                    # Then, call np.full with the shape of the reduction array and the identity value.
                    full_call = pfbdr.call(
                        full_func_node, args=[redshape_var, init_val_var, dt],
                    )

                    redtoset = pfbdr.assign(
                        rhs=full_call,
                        typ=redvar_typ,
                        name="redtoset",
                    )
                    # rettoset is an array from np.full() and must be released
                    to_cleanup.append(redtoset)
                else:
                    redtoset = pfbdr.make_const_variable(
                        cval=init_val,
                        typ=reddtype,
                        name="redtoset",
                    )
            else:
                redtoset = redvar

                if config.DEBUG_ARRAY_OPT_RUNTIME:
                    res_print_str = "res_print1 for redvar " + str(redvar) + ":"
                    strconsttyp = types.StringLiteral(res_print_str)

                    lhs = pfbdr.make_const_variable(
                        cval=res_print_str,
                        typ=strconsttyp,
                        name="str_const",
                    )

                    res_print = ir.Print(args=[lhs, redvar],
                                         vararg=None, loc=loc)
                    lowerer.fndesc.calltypes[res_print] = signature(types.none,
                                                             typemap[lhs.name],
                                                             typemap[redvar.name])
                    print("res_print_redvar", res_print)
                    lowerer.lower_inst(res_print)


            # For each thread, initialize the per-worker reduction array to
            # the current reduction array value.

            # Get the Numba type of the variable that holds the thread count.
            num_thread_type = typemap[num_threads_var.name]
            # Get the LLVM type of the thread count variable.
            ntllvm_type = targetctx.get_value_type(num_thread_type)
            # Create a LLVM variable to hold the loop index.
            alloc_loop_var = cgutils.alloca_once(builder, ntllvm_type)
            # Associate this LLVM variable to a Numba IR variable so that
            # we can use setitem IR builder.
            # Create a Numba IR variable.
            numba_ir_loop_index_var = scope.redefine("$loop_index", loc)
            # Give that variable the right type.
            typemap[numba_ir_loop_index_var.name] = num_thread_type
            # Associate this Numba variable to the LLVM variable in the
            # lowerer's varmap.
            lowerer.varmap[numba_ir_loop_index_var.name] = alloc_loop_var
            # Insert a loop into the outputed LLVM that goes from 0 to
            # the current thread count.
            with cgutils.for_range(builder, lowerer.loadvar(num_threads_var.name), intp=ntllvm_type) as loop:
                # Store the loop index into the alloca'd LLVM loop index variable.
                builder.store(loop.index, alloc_loop_var)
                # Initialize one element of the reduction array using the Numba
                # IR variable associated with this loop's index.
                pfbdr.setitem(obj=redarr_var, index=numba_ir_loop_index_var, val=redtoset)

    # compile parfor body as a separate function to be used with GUFuncWrapper
    flags = parfor.flags.copy()
    flags.error_model = "numpy"
    # Can't get here unless  flags.auto_parallel == ParallelOptions(True)
    index_var_typ = typemap[parfor.loop_nests[0].index_variable.name]
    # index variables should have the same type, check rest of indices
    for l in parfor.loop_nests[1:]:
        assert typemap[l.index_variable.name] == index_var_typ
    numba.parfors.parfor.sequential_parfor_lowering = True
    try:
        (func,
         func_args,
         func_sig,
         func_arg_types,
         exp_name_to_tuple_var) = _create_gufunc_for_parfor_body(
            lowerer, parfor, typemap, typingctx, targetctx, flags, {},
            bool(alias_map), index_var_typ, parfor.races)
    finally:
        numba.parfors.parfor.sequential_parfor_lowering = False

    # get the shape signature
    func_args = ['sched'] + func_args
    num_reductions = len(parfor_redvars)
    num_inputs = len(func_args) - len(parfor_output_arrays) - num_reductions
    if config.DEBUG_ARRAY_OPT:
        print("func_args = ", func_args)
        print("num_inputs = ", num_inputs)
        print("parfor_outputs = ", parfor_output_arrays)
        print("parfor_redvars = ", parfor_redvars)
        print("num_reductions = ", num_reductions)
    gu_signature = _create_shape_signature(
        parfor.get_shape_classes,
        num_inputs,
        num_reductions,
        func_args,
        func_sig,
        parfor.races,
        typemap)
    if config.DEBUG_ARRAY_OPT:
        print("gu_signature = ", gu_signature)

    # call the func in parallel by wrapping it with ParallelGUFuncBuilder
    loop_ranges = [(l.start, l.stop, l.step) for l in parfor.loop_nests]
    if config.DEBUG_ARRAY_OPT:
        print("loop_nests = ", parfor.loop_nests)
        print("loop_ranges = ", loop_ranges)
    call_parallel_gufunc(
        lowerer,
        func,
        gu_signature,
        func_sig,
        func_args,
        func_arg_types,
        loop_ranges,
        parfor_redvars,
        parfor_reddict,
        redarrs,
        parfor.init_block,
        index_var_typ,
        parfor.races,
        exp_name_to_tuple_var)

    if nredvars > 0:
        _parfor_lowering_finalize_reduction(
            parfor, redarrs, lowerer, parfor_reddict, num_threads_var,
        )

    # Cleanup reduction variable
    for v in to_cleanup:
        lowerer.lower_inst(ir.Del(v.name, loc=loc))
    # Restore the original typemap of the function that was replaced temporarily at the
    # Beginning of this function.
    lowerer.fndesc.typemap = orig_typemap

    if config.DEBUG_ARRAY_OPT:
        print("_lower_parfor_parallel done")


_ReductionInfo = make_dataclass(
    "_ReductionInfo",
    [
        "redvar_info",
        "redvar_name",
        "redvar_typ",
        "redarr_var",
        "redarr_typ",
        "init_val",
    ],
    frozen=True,
)


def _parfor_lowering_finalize_reduction(
        parfor,
        redarrs,
        lowerer,
        parfor_reddict,
        thread_count_var,
    ):
    """Emit code to finalize the reduction from the intermediate values of
    each thread.
    """
    # For each reduction variable
    for redvar_name, redarr_var in redarrs.items():
        # Pseudo-code for this loop body:
        #     tmp = redarr[0]
        #     for i in range(1, thread_count):
        #         tmp = reduce_op(redarr[i], tmp)
        #     reduction_result = tmp
        redvar_typ = lowerer.fndesc.typemap[redvar_name]
        redarr_typ = lowerer.fndesc.typemap[redarr_var.name]
        init_val = lowerer.loadvar(redvar_name)

        reduce_info = _ReductionInfo(
            redvar_info = parfor_reddict[redvar_name],
            redvar_name=redvar_name,
            redvar_typ=redvar_typ,
            redarr_var=redarr_var,
            redarr_typ=redarr_typ,
            init_val=init_val,
        )
        # generate code for combining reduction variable with thread output
        handler = (_lower_trivial_inplace_binops
                   if reduce_info.redvar_info.redop is not None
                   else _lower_non_trivial_reduce)
        handler(parfor, lowerer, thread_count_var, reduce_info)


class ParforsUnexpectedReduceNodeError(InternalError):
    def __init__(self, inst):
        super().__init__(f"Unknown reduce instruction node: {inst}")


def _lower_trivial_inplace_binops(parfor, lowerer, thread_count_var, reduce_info):
    """Lower trivial inplace-binop reduction.
    """
    for inst in reduce_info.redvar_info.reduce_nodes:
        # Var assigns to Var?
        if _lower_var_to_var_assign(lowerer, inst):
            pass
        # Is inplace-binop for the reduction?
        elif _is_right_op_and_rhs_is_init(inst, reduce_info.redvar_name, "inplace_binop"):
            fn = inst.value.fn
            redvar_result = _emit_binop_reduce_call(
                fn, lowerer, thread_count_var, reduce_info,
            )
            lowerer.storevar(redvar_result, name=inst.target.name)
        # Is binop for the reduction?
        elif _is_right_op_and_rhs_is_init(inst, reduce_info.redvar_name, "binop"):
            fn = inst.value.fn
            redvar_result = _emit_binop_reduce_call(
                fn, lowerer, thread_count_var, reduce_info,
            )
            lowerer.storevar(redvar_result, name=inst.target.name)
        # Otherwise?
        else:
            raise ParforsUnexpectedReduceNodeError(inst)

        # XXX: This seems like a hack to stop the loop with this condition.
        if _fix_redvar_name_ssa_mismatch(parfor, lowerer, inst,
                                   reduce_info.redvar_name):
            break
    if config.DEBUG_ARRAY_OPT_RUNTIME:
        varname = reduce_info.redvar_name
        lowerer.print_variable(
            f"{parfor.loc}: parfor {fn.__name__} reduction {varname} =",
            varname,
        )


def _lower_non_trivial_reduce(parfor, lowerer, thread_count_var, reduce_info):
    """Lower non-trivial reduction such as call to `functools.reduce()`.
    """
    init_name = f"{reduce_info.redvar_name}#init"
    # The init_name variable is not defined at this point.
    lowerer.fndesc.typemap.setdefault(init_name, reduce_info.redvar_typ)
    # Emit a sequence of the reduction operation for each intermediate result
    # of each thread.
    num_thread_llval = lowerer.loadvar(thread_count_var.name)
    with cgutils.for_range(lowerer.builder, num_thread_llval) as loop:
        tid = loop.index
        for inst in reduce_info.redvar_info.reduce_nodes:
            # Var assigns to Var?
            if _lower_var_to_var_assign(lowerer, inst):
                pass
            # The reduction operation?
            elif (isinstance(inst, ir.Assign)
                    and any(var.name == init_name for var in inst.list_vars())):
                elem = _emit_getitem_call(tid, lowerer, reduce_info)
                lowerer.storevar(elem, init_name)
                lowerer.lower_inst(inst)

            # Otherwise?
            else:
                raise ParforsUnexpectedReduceNodeError(inst)

            # XXX: This seems like a hack to stop the loop with this condition.
            if _fix_redvar_name_ssa_mismatch(parfor, lowerer, inst,
                                       reduce_info.redvar_name):
                break

    if config.DEBUG_ARRAY_OPT_RUNTIME:
        varname = reduce_info.redvar_name
        lowerer.print_variable(
            f"{parfor.loc}: parfor non-trivial reduction {varname} =",
            varname,
        )

def _lower_var_to_var_assign(lowerer, inst):
    """Lower Var->Var assignment.

    Returns True if-and-only-if `inst` is a Var->Var assignment.
    """
    if isinstance(inst, ir.Assign) and isinstance(inst.value, ir.Var):
        loaded = lowerer.loadvar(inst.value.name)
        lowerer.storevar(loaded, name=inst.target.name)
        return True
    return False

def _emit_getitem_call(idx, lowerer, reduce_info):
    """Emit call to ``redarr_var[idx]``
    """
    def reducer_getitem(redarr, index):
        return redarr[index]

    builder = lowerer.builder
    ctx = lowerer.context
    redarr_typ = reduce_info.redarr_typ
    arg_arr = lowerer.loadvar(reduce_info.redarr_var.name)
    args = (arg_arr, idx)
    sig = signature(reduce_info.redvar_typ, redarr_typ, types.intp)
    elem = ctx.compile_internal(builder, reducer_getitem, sig, args)
    return elem


def _emit_binop_reduce_call(binop, lowerer, thread_count_var, reduce_info):
    """Emit call to the ``binop`` for the reduction variable.
    """

    def reduction_add(thread_count, redarr, init):
        c = init
        for i in range(thread_count):
            c += redarr[i]
        return c

    def reduction_mul(thread_count, redarr, init):
        c = init
        for i in range(thread_count):
            c *= redarr[i]
        return c

    kernel = {
        operator.iadd: reduction_add,
        operator.isub: reduction_add,
        operator.add: reduction_add,
        operator.sub: reduction_add,
        operator.imul: reduction_mul,
        operator.ifloordiv: reduction_mul,
        operator.itruediv: reduction_mul,
        operator.mul: reduction_mul,
        operator.floordiv: reduction_mul,
        operator.truediv: reduction_mul,
    }[binop]

    ctx = lowerer.context
    builder = lowerer.builder
    redarr_typ = reduce_info.redarr_typ
    arg_arr = lowerer.loadvar(reduce_info.redarr_var.name)

    if config.DEBUG_ARRAY_OPT_RUNTIME:
        init_var = reduce_info.redarr_var.scope.get(reduce_info.redvar_name)
        res_print = ir.Print(
            args=[reduce_info.redarr_var, init_var], vararg=None,
            loc=lowerer.loc,
        )
        typemap = lowerer.fndesc.typemap
        lowerer.fndesc.calltypes[res_print] = signature(
            types.none, typemap[reduce_info.redarr_var.name],
            typemap[init_var.name],
        )
        lowerer.lower_inst(res_print)

    arg_thread_count = lowerer.loadvar(thread_count_var.name)
    args = (arg_thread_count, arg_arr, reduce_info.init_val)
    sig = signature(
        reduce_info.redvar_typ, types.uintp, redarr_typ, reduce_info.redvar_typ,
    )

    redvar_result = ctx.compile_internal(builder, kernel, sig, args)
    return redvar_result


def _is_right_op_and_rhs_is_init(inst, redvar_name, op):
    """Is ``inst`` an inplace-binop and the RHS is the reduction init?
    """
    if not isinstance(inst, ir.Assign):
        return False
    rhs = inst.value
    if not isinstance(rhs, ir.Expr):
        return False
    if rhs.op != op:
        return False
    if rhs.rhs.name != f"{redvar_name}#init":
        return False
    return True


def _fix_redvar_name_ssa_mismatch(parfor, lowerer, inst, redvar_name):
    """Fix reduction variable name mismatch due to SSA.
    """
    # Only process reduction statements post-gufunc execution
    # until we see an assignment with a left-hand side to the
    # reduction variable's name.  This fixes problems with
    # cases where there are multiple assignments to the
    # reduction variable in the parfor.
    scope = parfor.init_block.scope
    if isinstance(inst, ir.Assign):
        try:
            reduction_var = scope.get_exact(redvar_name)
        except NotDefinedError:
            # Ideally, this shouldn't happen. The redvar name
            # missing from scope indicates an error from
            # other rewrite passes.
            is_same_source_var = redvar_name == inst.target.name
        else:
            # Because of SSA, the redvar and target var of
            # the current assignment would be different even
            # though they refer to the same source-level var.
            redvar_unver_name = reduction_var.unversioned_name
            target_unver_name = inst.target.unversioned_name
            is_same_source_var = redvar_unver_name == target_unver_name

        if is_same_source_var:
            # If redvar is different from target var, add an
            # assignment to put target var into redvar.
            if redvar_name != inst.target.name:
                val = lowerer.loadvar(inst.target.name)
                lowerer.storevar(val, name=redvar_name)
                return True

    return False

def _create_shape_signature(
        get_shape_classes,
        num_inputs,
        num_reductions,
        args,
        func_sig,
        races,
        typemap):
    '''Create shape signature for GUFunc
    '''
    if config.DEBUG_ARRAY_OPT:
        print("_create_shape_signature", num_inputs, num_reductions, args, races)
        for i in args[1:]:
            print("argument", i, type(i), get_shape_classes(i, typemap=typemap))

    num_inouts = len(args) - num_reductions
    # maximum class number for array shapes
    classes = [get_shape_classes(var, typemap=typemap) if var not in races else (-1,) for var in args[1:]]
    class_set = set()
    for _class in classes:
        if _class:
            for i in _class:
                class_set.add(i)
    max_class = max(class_set) + 1 if class_set else 0
    classes.insert(0, (max_class,)) # force set the class of 'sched' argument
    class_set.add(max_class)
    thread_num_class = max_class + 1
    class_set.add(thread_num_class)
    class_map = {}
    # TODO: use prefix + class number instead of single char
    alphabet = ord('a')
    for n in class_set:
       if n >= 0:
           class_map[n] = chr(alphabet)
           alphabet += 1
    threadcount_ordinal = chr(alphabet)

    alpha_dict = {'latest_alpha' : alphabet}

    def bump_alpha(c, class_map):
        if c >= 0:
            return class_map[c]
        else:
            alpha_dict['latest_alpha'] += 1
            return chr(alpha_dict['latest_alpha'])

    gu_sin = []
    gu_sout = []
    count = 0
    syms_sin = ()
    if config.DEBUG_ARRAY_OPT:
        print("args", args)
        print("classes", classes)
        print("threadcount_ordinal", threadcount_ordinal)
    for cls, arg in zip(classes, args):
        count = count + 1
        if cls:
            dim_syms = tuple(bump_alpha(c, class_map) for c in cls)
        else:
            dim_syms = ()
        if (count > num_inouts):
            # Add the threadcount_ordinal to represent the thread count
            # to the start of the reduction array.
            gu_sin.append(tuple([threadcount_ordinal] + list(dim_syms[1:])))
        else:
            gu_sin.append(dim_syms)
            syms_sin += dim_syms
    return (gu_sin, gu_sout)

def _print_block(block):
    for i, inst in enumerate(block.body):
        print("    ", i, " ", inst)

def _print_body(body_dict):
    '''Pretty-print a set of IR blocks.
    '''
    topo_order = wrap_find_topo(body_dict)
    for label in topo_order:
        block = body_dict[label]
        print("label: ", label)
        _print_block(block)


def wrap_loop_body(loop_body):
    blocks = loop_body.copy()  # shallow copy is enough
    first_label = min(blocks.keys())
    last_label = max(blocks.keys())
    loc = blocks[last_label].loc
    blocks[last_label].body.append(ir.Jump(first_label, loc))
    return blocks

def unwrap_loop_body(loop_body):
    last_label = max(loop_body.keys())
    loop_body[last_label].body = loop_body[last_label].body[:-1]

def add_to_def_once_sets(a_def, def_once, def_more):
    '''If the variable is already defined more than once, do nothing.
       Else if defined exactly once previously then transition this
       variable to the defined more than once set (remove it from
       def_once set and add to def_more set).
       Else this must be the first time we've seen this variable defined
       so add to def_once set.
    '''
    if a_def in def_more:
        pass
    elif a_def in def_once:
        def_more.add(a_def)
        def_once.remove(a_def)
    else:
        def_once.add(a_def)

def compute_def_once_block(block, def_once, def_more, getattr_taken, typemap, module_assigns):
    '''Effect changes to the set of variables defined once or more than once
       for a single block.
       block - the block to process
       def_once - set of variable names known to be defined exactly once
       def_more - set of variable names known to be defined more than once
       getattr_taken - dict mapping variable name to tuple of object and attribute taken
       module_assigns - dict mapping variable name to the Global that they came from
    '''
    # The only "defs" occur in assignments, so find such instructions.
    assignments = block.find_insts(ir.Assign)
    # For each assignment...
    for one_assign in assignments:
        # Get the LHS/target of the assignment.
        a_def = one_assign.target.name
        # Add variable to def sets.
        add_to_def_once_sets(a_def, def_once, def_more)

        rhs = one_assign.value
        if isinstance(rhs, ir.Global):
            # Remember assignments of the form "a = Global(...)"
            # Is this a module?
            if isinstance(rhs.value, pytypes.ModuleType):
                module_assigns[a_def] = rhs.value.__name__
        if isinstance(rhs, ir.Expr) and rhs.op == 'getattr' and rhs.value.name in def_once:
            # Remember assignments of the form "a = b.c"
            getattr_taken[a_def] = (rhs.value.name, rhs.attr)
        if isinstance(rhs, ir.Expr) and rhs.op == 'call' and rhs.func.name in getattr_taken:
            # If "a" is being called then lookup the getattr definition of "a"
            # as above, getting the module variable "b" (base_obj)
            # and the attribute "c" (base_attr).
            base_obj, base_attr = getattr_taken[rhs.func.name]
            if base_obj in module_assigns:
                # If we know the definition of the module variable then get the module
                # name from module_assigns.
                base_mod_name = module_assigns[base_obj]
                if not is_const_call(base_mod_name, base_attr):
                    # Calling a method on an object could modify the object and is thus
                    # like a def of that object.  We call is_const_call to see if this module/attribute
                    # combination is known to not modify the module state.  If we don't know that
                    # the combination is safe then we have to assume there could be a modification to
                    # the module and thus add the module variable as defined more than once.
                    add_to_def_once_sets(base_obj, def_once, def_more)
            else:
                # Assume the worst and say that base_obj could be modified by the call.
                add_to_def_once_sets(base_obj, def_once, def_more)
        if isinstance(rhs, ir.Expr) and rhs.op == 'call':
            # If a mutable object is passed to a function, then it may be changed and
            # therefore can't be hoisted.
            # For each argument to the function...
            for argvar in rhs.args:
                # Get the argument's type.
                if isinstance(argvar, ir.Var):
                    argvar = argvar.name
                avtype = typemap[argvar]
                # If that type doesn't have a mutable attribute or it does and it's set to
                # not mutable then this usage is safe for hoisting.
                if getattr(avtype, 'mutable', False):
                    # Here we have a mutable variable passed to a function so add this variable
                    # to the def lists.
                    add_to_def_once_sets(argvar, def_once, def_more)

def wrap_find_topo(loop_body):
    blocks = wrap_loop_body(loop_body)
    topo_order = find_topo_order(blocks)
    unwrap_loop_body(loop_body)
    return topo_order

def compute_def_once_internal(loop_body, def_once, def_more, getattr_taken, typemap, module_assigns):
    '''Compute the set of variables defined exactly once in the given set of blocks
       and use the given sets for storing which variables are defined once, more than
       once and which have had a getattr call on them.
    '''
    # For each block in topological order...
    topo_order = wrap_find_topo(loop_body)
    for label in topo_order:
        block = loop_body[label]
        # Scan this block and effect changes to def_once, def_more, and getattr_taken
        # based on the instructions in that block.
        compute_def_once_block(block, def_once, def_more, getattr_taken, typemap, module_assigns)
        # Have to recursively process parfors manually here.
        for inst in block.body:
            if isinstance(inst, parfor.Parfor):
                # Recursively compute for the parfor's init block.
                compute_def_once_block(inst.init_block, def_once, def_more, getattr_taken, typemap, module_assigns)
                # Recursively compute for the parfor's loop body.
                compute_def_once_internal(inst.loop_body, def_once, def_more, getattr_taken, typemap, module_assigns)

def compute_def_once(loop_body, typemap):
    '''Compute the set of variables defined exactly once in the given set of blocks.
    '''
    def_once = set()   # set to hold variables defined exactly once
    def_more = set()   # set to hold variables defined more than once
    getattr_taken = {}
    module_assigns = {}
    compute_def_once_internal(loop_body, def_once, def_more, getattr_taken, typemap, module_assigns)
    return def_once, def_more

def find_vars(var, varset):
    assert isinstance(var, ir.Var)
    varset.add(var.name)
    return var

def _hoist_internal(inst, dep_on_param, call_table, hoisted, not_hoisted,
                    typemap, stored_arrays):
    if inst.target.name in stored_arrays:
        not_hoisted.append((inst, "stored array"))
        if config.DEBUG_ARRAY_OPT >= 1:
            print("Instruction", inst, "could not be hoisted because the created array is stored.")
        return False

    target_type = typemap[inst.target.name]

    uses = set()
    # Get vars used by this statement.
    visit_vars_inner(inst.value, find_vars, uses)
    # Filter out input parameters from the set of variable usages.
    unhoistable = {assgn.target.name for assgn, _ in not_hoisted}
    use_unhoist = uses & unhoistable
    diff = uses.difference(dep_on_param)
    diff |= use_unhoist
    if config.DEBUG_ARRAY_OPT >= 1:
        print("_hoist_internal:", inst, "uses:", uses, "diff:", diff)
    if len(diff) == 0 and is_pure(inst.value, None, call_table):
        if config.DEBUG_ARRAY_OPT >= 1:
            print("Will hoist instruction", inst, target_type)
        hoisted.append(inst)
        if not isinstance(target_type, types.npytypes.Array):
            dep_on_param += [inst.target.name]
        return True
    else:
        if len(diff) > 0:
            not_hoisted.append((inst, "dependency"))
            if config.DEBUG_ARRAY_OPT >= 1:
                print("Instruction", inst, "could not be hoisted because of a dependency.")
        else:
            not_hoisted.append((inst, "not pure"))
            if config.DEBUG_ARRAY_OPT >= 1:
                print("Instruction", inst, "could not be hoisted because it isn't pure.")
    return False

def find_setitems_block(setitems, itemsset, block, typemap):
    for inst in block.body:
        if isinstance(inst, (ir.StaticSetItem, ir.SetItem)):
            setitems.add(inst.target.name)
            # If we store a non-mutable object into an array then that is safe to hoist.
            # If the stored object is mutable and you hoist then multiple entries in the
            # outer array could reference the same object and changing one index would then
            # change other indices.
            if getattr(typemap[inst.value.name], "mutable", False):
                itemsset.add(inst.value.name)
        elif isinstance(inst, parfor.Parfor):
            find_setitems_block(setitems, itemsset, inst.init_block, typemap)
            find_setitems_body(setitems, itemsset, inst.loop_body, typemap)
        elif isinstance(inst, ir.Assign):
            # If something of mutable type is given to a build_tuple or
            # used in a call then consider it unanalyzable and so
            # unavailable for hoisting.
            rhs = inst.value
            def add_to_itemset(item):
                assert isinstance(item, ir.Var), rhs
                if getattr(typemap[item.name], "mutable", False):
                    itemsset.add(item.name)

            if isinstance(rhs, ir.Expr):
                if rhs.op in ["build_tuple", "build_list", "build_set"]:
                    for item in rhs.items:
                        add_to_itemset(item)
                elif rhs.op == "build_map":
                    for pair in rhs.items:
                        for item in pair:
                            add_to_itemset(item)
                elif rhs.op == "call":
                    for item in list(rhs.args) + [x[1] for x in rhs.kws]:
                        add_to_itemset(item)

def find_setitems_body(setitems, itemsset, loop_body, typemap):
    """
      Find the arrays that are written into (goes into setitems) and the
      mutable objects (mostly arrays) that are written into other arrays
      (goes into itemsset).
    """
    for label, block in loop_body.items():
        find_setitems_block(setitems, itemsset, block, typemap)

def empty_container_allocator_hoist(inst, dep_on_param, call_table, hoisted,
                                    not_hoisted, typemap, stored_arrays):
    if (isinstance(inst, ir.Assign) and
        isinstance(inst.value, ir.Expr) and
        inst.value.op == 'call' and
        inst.value.func.name in call_table):
        call_list = call_table[inst.value.func.name]
        if call_list == ['empty', np]:
            return _hoist_internal(inst, dep_on_param, call_table, hoisted,
                                   not_hoisted, typemap, stored_arrays)
    return False

def hoist(parfor_params, loop_body, typemap, wrapped_blocks):
    dep_on_param = copy.copy(parfor_params)
    hoisted = []
    not_hoisted = []

    # Compute the set of variable defined exactly once in the loop body.
    def_once, def_more = compute_def_once(loop_body, typemap)
    (call_table, reverse_call_table) = get_call_table(wrapped_blocks)

    setitems = set()
    itemsset = set()
    find_setitems_body(setitems, itemsset, loop_body, typemap)
    dep_on_param = list(set(dep_on_param).difference(setitems))
    if config.DEBUG_ARRAY_OPT >= 1:
        print("hoist - def_once:", def_once, "setitems:", setitems, "itemsset:", itemsset, "dep_on_param:", dep_on_param, "parfor_params:", parfor_params)
    for si in setitems:
        add_to_def_once_sets(si, def_once, def_more)

    for label, block in loop_body.items():
        new_block = []
        for inst in block.body:
            if empty_container_allocator_hoist(inst, dep_on_param, call_table,
                                   hoisted, not_hoisted, typemap, itemsset):
                continue
            elif isinstance(inst, ir.Assign) and inst.target.name in def_once:
                if _hoist_internal(inst, dep_on_param, call_table,
                                   hoisted, not_hoisted, typemap, itemsset):
                    # don't add this instruction to the block since it is
                    # hoisted
                    continue
            elif isinstance(inst, parfor.Parfor):
                new_init_block = []
                if config.DEBUG_ARRAY_OPT >= 1:
                    print("parfor")
                    inst.dump()
                for ib_inst in inst.init_block.body:
                    if empty_container_allocator_hoist(ib_inst, dep_on_param,
                        call_table, hoisted, not_hoisted, typemap, itemsset):
                        continue
                    elif (isinstance(ib_inst, ir.Assign) and
                        ib_inst.target.name in def_once):
                        if _hoist_internal(ib_inst, dep_on_param, call_table,
                                           hoisted, not_hoisted, typemap,
                                           itemsset):
                            # don't add this instruction to the block since it is hoisted
                            continue
                    new_init_block.append(ib_inst)
                inst.init_block.body = new_init_block

            new_block.append(inst)
        block.body = new_block
    return hoisted, not_hoisted

def redtyp_is_scalar(redtype):
    return not isinstance(redtype, types.npytypes.Array)

def redtyp_to_redarraytype(redtyp):
    """Go from a reducation variable type to a reduction array type used to hold
       per-worker results.
    """
    redarrdim = 1
    # If the reduction type is an array then allocate reduction array with ndim+1 dimensions.
    if isinstance(redtyp, types.npytypes.Array):
        redarrdim += redtyp.ndim
        # We don't create array of array but multi-dimensional reduction array with same dtype.
        redtyp = redtyp.dtype
    return types.npytypes.Array(redtyp, redarrdim, "C")

def redarraytype_to_sig(redarraytyp):
    """Given a reduction array type, find the type of the reduction argument to the gufunc.
    """
    assert isinstance(redarraytyp, types.npytypes.Array)
    return types.npytypes.Array(redarraytyp.dtype, redarraytyp.ndim, redarraytyp.layout)

def legalize_names_with_typemap(names, typemap):
    """ We use ir_utils.legalize_names to replace internal IR variable names
        containing illegal characters (e.g. period) with a legal character
        (underscore) so as to create legal variable names.
        The original variable names are in the typemap so we also
        need to add the legalized name to the typemap as well.
    """
    outdict = legalize_names(names)
    # For each pair in the dict of legalized names...
    for x, y in outdict.items():
        # If the name had some legalization change to it...
        if x != y:
            # Set the type of the new name the same as the type of the old name.
            typemap[y] = typemap[x]
    return outdict

def to_scalar_from_0d(x):
    if isinstance(x, types.ArrayCompatible):
        if x.ndim == 0:
            return x.dtype
    return x

def _create_gufunc_for_parfor_body(
        lowerer,
        parfor,
        typemap,
        typingctx,
        targetctx,
        flags,
        locals,
        has_aliases,
        index_var_typ,
        races):
    '''
    Takes a parfor and creates a gufunc function for its body.
    There are two parts to this function.
    1) Code to iterate across the iteration space as defined by the schedule.
    2) The parfor body that does the work for a single point in the iteration space.
    Part 1 is created as Python text for simplicity with a sentinel assignment to mark the point
    in the IR where the parfor body should be added.
    This Python text is 'exec'ed into existence and its IR retrieved with run_frontend.
    The IR is scanned for the sentinel assignment where that basic block is split and the IR
    for the parfor body inserted.
    '''
    if config.DEBUG_ARRAY_OPT >= 1:
        print("starting _create_gufunc_for_parfor_body")

    loc = parfor.init_block.loc

    # The parfor body and the main function body share ir.Var nodes.
    # We have to do some replacements of Var names in the parfor body to make them
    # legal parameter names.  If we don't copy then the Vars in the main function also
    # would incorrectly change their name.
    loop_body = copy.copy(parfor.loop_body)
    remove_dels(loop_body)

    parfor_dim = len(parfor.loop_nests)
    loop_indices = [l.index_variable.name for l in parfor.loop_nests]

    # Get all the parfor params.
    parfor_params = parfor.params
    # Get just the outputs of the parfor.
    parfor_outputs = numba.parfors.parfor.get_parfor_outputs(parfor, parfor_params)
    # Get all parfor reduction vars, and operators.
    typemap = lowerer.fndesc.typemap
    parfor_redvars, parfor_reddict = numba.parfors.parfor.get_parfor_reductions(
        lowerer.func_ir, parfor, parfor_params, lowerer.fndesc.calltypes)
    # Compute just the parfor inputs as a set difference.
    parfor_inputs = sorted(
        list(
            set(parfor_params) -
            set(parfor_outputs) -
            set(parfor_redvars)))

    if config.DEBUG_ARRAY_OPT >= 1:
        print("parfor_params = ", parfor_params, " ", type(parfor_params))
        print("parfor_outputs = ", parfor_outputs, " ", type(parfor_outputs))
        print("parfor_inputs = ", parfor_inputs, " ", type(parfor_inputs))
        print("parfor_redvars = ", parfor_redvars, " ", type(parfor_redvars))

    # -------------------------------------------------------------------------
    # Convert tuples to individual parameters.
    tuple_expanded_parfor_inputs = []
    tuple_var_to_expanded_names = {}
    expanded_name_to_tuple_var = {}
    next_expanded_tuple_var = 0
    parfor_tuple_params = []
    # For each input to the parfor.
    for pi in parfor_inputs:
        # Get the type of the input.
        pi_type = typemap[pi]
        # If it is a UniTuple or Tuple we will do the conversion.
        if isinstance(pi_type, types.UniTuple) or isinstance(pi_type, types.NamedUniTuple):
            # Get the size and dtype of the tuple.
            tuple_count = pi_type.count
            tuple_dtype = pi_type.dtype
            # Only do tuples up to config.PARFOR_MAX_TUPLE_SIZE length.
            assert(tuple_count <= config.PARFOR_MAX_TUPLE_SIZE)
            this_var_expansion = []
            for i in range(tuple_count):
                # Generate a new name for the individual part of the tuple var.
                expanded_name = "expanded_tuple_var_" + str(next_expanded_tuple_var)
                # Add that name to the new list of inputs to the gufunc.
                tuple_expanded_parfor_inputs.append(expanded_name)
                this_var_expansion.append(expanded_name)
                # Remember a mapping from new param name to original tuple
                # var and the index within the tuple.
                expanded_name_to_tuple_var[expanded_name] = (pi, i)
                next_expanded_tuple_var += 1
                # Set the type of the new parameter.
                typemap[expanded_name] = tuple_dtype
            # Remember a mapping from the original tuple var to the
            # individual parts.
            tuple_var_to_expanded_names[pi] = this_var_expansion
            parfor_tuple_params.append(pi)
        elif isinstance(pi_type, types.Tuple) or isinstance(pi_type, types.NamedTuple):
            # This is the same as above for UniTuple except that each part of
            # the tuple can have a different type and we fetch that type with
            # pi_type.types[offset].
            tuple_count = pi_type.count
            tuple_types = pi_type.types
            # Only do tuples up to config.PARFOR_MAX_TUPLE_SIZE length.
            assert(tuple_count <= config.PARFOR_MAX_TUPLE_SIZE)
            this_var_expansion = []
            for i in range(tuple_count):
                expanded_name = "expanded_tuple_var_" + str(next_expanded_tuple_var)
                tuple_expanded_parfor_inputs.append(expanded_name)
                this_var_expansion.append(expanded_name)
                expanded_name_to_tuple_var[expanded_name] = (pi, i)
                next_expanded_tuple_var += 1
                typemap[expanded_name] = tuple_types[i]
            tuple_var_to_expanded_names[pi] = this_var_expansion
            parfor_tuple_params.append(pi)
        else:
            tuple_expanded_parfor_inputs.append(pi)
    parfor_inputs = tuple_expanded_parfor_inputs
    if config.DEBUG_ARRAY_OPT >= 1:
        print("parfor_inputs post tuple handling = ", parfor_inputs, " ", type(parfor_inputs))
    # -------------------------------------------------------------------------

    races = races.difference(set(parfor_redvars))
    for race in races:
        msg = ("Variable %s used in parallel loop may be written "
               "to simultaneously by multiple workers and may result "
               "in non-deterministic or unintended results." % race)
        warnings.warn(NumbaParallelSafetyWarning(msg, loc))
    replace_var_with_array(races, loop_body, typemap, lowerer.fndesc.calltypes)

    # Reduction variables are represented as arrays, so they go under
    # different names.
    parfor_redarrs = []
    parfor_red_arg_types = []
    for var in parfor_redvars:
        arr = var + "_arr"
        parfor_redarrs.append(arr)
        redarraytype = redtyp_to_redarraytype(typemap[var])
        parfor_red_arg_types.append(redarraytype)
        redarrsig = redarraytype_to_sig(redarraytype)
        if arr in typemap:
            assert(typemap[arr] == redarrsig)
        else:
            typemap[arr] = redarrsig

    # Reorder all the params so that inputs go first then outputs.
    parfor_params = parfor_inputs + parfor_outputs + parfor_redarrs

    if config.DEBUG_ARRAY_OPT >= 1:
        print("parfor_params = ", parfor_params, " ", type(parfor_params))
        print("loop_indices = ", loop_indices, " ", type(loop_indices))
        print("loop_body = ", loop_body, " ", type(loop_body))
        _print_body(loop_body)

    # Some Var are not legal parameter names so create a dict of potentially illegal
    # param name to guaranteed legal name.
    param_dict = legalize_names_with_typemap(parfor_params + parfor_redvars + parfor_tuple_params, typemap)
    if config.DEBUG_ARRAY_OPT >= 1:
        print(
            "param_dict = ",
            sorted(
                param_dict.items()),
            " ",
            type(param_dict))

    # Some loop_indices are not legal parameter names so create a dict of potentially illegal
    # loop index to guaranteed legal name.
    ind_dict = legalize_names_with_typemap(loop_indices, typemap)
    # Compute a new list of legal loop index names.
    legal_loop_indices = [ind_dict[v] for v in loop_indices]
    if config.DEBUG_ARRAY_OPT >= 1:
        print("ind_dict = ", sorted(ind_dict.items()), " ", type(ind_dict))
        print(
            "legal_loop_indices = ",
            legal_loop_indices,
            " ",
            type(legal_loop_indices))
        for pd in parfor_params:
            print("pd = ", pd)
            print("pd type = ", typemap[pd], " ", type(typemap[pd]))

    # Get the types of each parameter.
    param_types = [to_scalar_from_0d(typemap[v]) for v in parfor_params]
    # Calculate types of args passed to gufunc.
    func_arg_types = [typemap[v] for v in (parfor_inputs + parfor_outputs)] + parfor_red_arg_types
    if config.DEBUG_ARRAY_OPT >= 1:
        print("new param_types:", param_types)
        print("new func_arg_types:", func_arg_types)

    # Replace illegal parameter names in the loop body with legal ones.
    replace_var_names(loop_body, param_dict)
    # remember the name before legalizing as the actual arguments
    parfor_args = parfor_params
    # Change parfor_params to be legal names.
    parfor_params = [param_dict[v] for v in parfor_params]
    parfor_params_orig = parfor_params

    parfor_params = []
    ascontig = False
    for pindex in range(len(parfor_params_orig)):
        if (ascontig and
            pindex < len(parfor_inputs) and
            isinstance(param_types[pindex], types.npytypes.Array)):
            parfor_params.append(parfor_params_orig[pindex]+"param")
        else:
            parfor_params.append(parfor_params_orig[pindex])

    # Change parfor body to replace illegal loop index vars with legal ones.
    replace_var_names(loop_body, ind_dict)
    loop_body_var_table = get_name_var_table(loop_body)
    sentinel_name = get_unused_var_name("__sentinel__", loop_body_var_table)

    if config.DEBUG_ARRAY_OPT >= 1:
        print(
            "legal parfor_params = ",
            parfor_params,
            " ",
            type(parfor_params))

    # Determine the unique names of the scheduling and gufunc functions.
    # sched_func_name = "__numba_parfor_sched_%s" % (hex(hash(parfor)).replace("-", "_"))
    gufunc_name = "__numba_parfor_gufunc_%s" % (
        hex(hash(parfor)).replace("-", "_"))
    if config.DEBUG_ARRAY_OPT:
        # print("sched_func_name ", type(sched_func_name), " ", sched_func_name)
        print("gufunc_name ", type(gufunc_name), " ", gufunc_name)

    gufunc_txt = ""

    # Create the gufunc function.
    gufunc_txt += "def " + gufunc_name + \
        "(sched, " + (", ".join(parfor_params)) + "):\n"

    globls = {"np": np, "numba": numba}

    # First thing in the gufunc, we reconstruct tuples from their
    # individual parts, e.g., orig_tup_name = (part1, part2,).
    # The rest of the code of the function will use the original tuple name.
    for tup_var, exp_names in tuple_var_to_expanded_names.items():
        tup_type = typemap[tup_var]
        gufunc_txt += "    " + param_dict[tup_var]
        # Determine if the tuple is a named tuple.
        if (isinstance(tup_type, types.NamedTuple) or
            isinstance(tup_type, types.NamedUniTuple)):
            named_tup = True
        else:
            named_tup = False

        if named_tup:
            # It is a named tuple so try to find the global that defines the
            # named tuple.
            func_def = guard(get_definition, lowerer.func_ir, tup_var)
            named_tuple_def = None
            if config.DEBUG_ARRAY_OPT:
                print("func_def:", func_def, type(func_def))
            if func_def is not None:
                if (isinstance(func_def, ir.Expr) and
                    func_def.op == 'call'):
                    named_tuple_def = guard(get_definition, lowerer.func_ir, func_def.func)
                    if config.DEBUG_ARRAY_OPT:
                        print("named_tuple_def:", named_tuple_def, type(named_tuple_def))
                elif isinstance(func_def, ir.Arg):
                    named_tuple_def = typemap[func_def.name]
                    if config.DEBUG_ARRAY_OPT:
                        print("named_tuple_def:", named_tuple_def,
                              type(named_tuple_def), named_tuple_def.name)
            if named_tuple_def is not None:
                if (isinstance(named_tuple_def, ir.Global) or
                    isinstance(named_tuple_def, ir.FreeVar)):
                    gval = named_tuple_def.value
                    if config.DEBUG_ARRAY_OPT:
                        print("gval:", gval, type(gval))
                    globls[named_tuple_def.name] = gval
                elif isinstance(named_tuple_def, types.containers.BaseNamedTuple):
                    named_tuple_name = named_tuple_def.name.split('(')[0]
                    if config.DEBUG_ARRAY_OPT:
                        print("name:", named_tuple_name,
                              named_tuple_def.instance_class,
                              type(named_tuple_def.instance_class))
                    globls[named_tuple_name] = named_tuple_def.instance_class
            else:
                if config.DEBUG_ARRAY_OPT:
                    print("Didn't find definition of namedtuple for globls.")
                raise CompilerError("Could not find definition of " + str(tup_var),
                                     tup_var.loc)
            gufunc_txt += " = " + tup_type.instance_class.__name__ + "("
            for name, field_name in zip(exp_names, tup_type.fields):
                gufunc_txt += field_name + "=" + param_dict[name] + ","
        else:
            # Just a regular tuple so use (part0, part1, ...)
            gufunc_txt += " = (" + ", ".join([param_dict[x] for x in exp_names])
            if len(exp_names) == 1:
                # Add comma for tuples with singular values.  We can't unilaterally
                # add a comma always because (,) isn't valid.
                gufunc_txt += ","

        gufunc_txt += ")\n"

    for pindex in range(len(parfor_inputs)):
        if ascontig and isinstance(param_types[pindex], types.npytypes.Array):
            gufunc_txt += ("    " + parfor_params_orig[pindex]
                + " = np.ascontiguousarray(" + parfor_params[pindex] + ")\n")

    gufunc_thread_id_var = "ParallelAcceleratorGufuncThreadId"

    if len(parfor_redarrs) > 0:
        gufunc_txt += "    " + gufunc_thread_id_var + " = "
        gufunc_txt += "numba.np.ufunc.parallel._iget_thread_id()\n"

    # Add initialization of reduction variables
    for arr, var in zip(parfor_redarrs, parfor_redvars):
        gufunc_txt += "    " + param_dict[var] + \
             "=" + param_dict[arr] + "[" + gufunc_thread_id_var + "]\n"
        if config.DEBUG_ARRAY_OPT_RUNTIME:
            gufunc_txt += "    print(\"thread id =\", ParallelAcceleratorGufuncThreadId)\n"
            gufunc_txt += "    print(\"initial reduction value\",ParallelAcceleratorGufuncThreadId," + param_dict[var] + "," + param_dict[var] + ".shape)\n"
            gufunc_txt += "    print(\"reduction array\",ParallelAcceleratorGufuncThreadId," + param_dict[arr] + "," + param_dict[arr] + ".shape)\n"

    # For each dimension of the parfor, create a for loop in the generated gufunc function.
    # Iterate across the proper values extracted from the schedule.
    # The form of the schedule is start_dim0, start_dim1, ..., start_dimN, end_dim0,
    # end_dim1, ..., end_dimN
    for eachdim in range(parfor_dim):
        for indent in range(eachdim + 1):
            gufunc_txt += "    "
        sched_dim = eachdim
        gufunc_txt += ("for " +
                       legal_loop_indices[eachdim] +
                       " in range(sched[" +
                       str(sched_dim) +
                       "], sched[" +
                       str(sched_dim +
                           parfor_dim) +
                       "] + np.uint8(1)):\n")

    if config.DEBUG_ARRAY_OPT_RUNTIME:
        for indent in range(parfor_dim + 1):
            gufunc_txt += "    "
        gufunc_txt += "print("
        for eachdim in range(parfor_dim):
            gufunc_txt += "\"" + legal_loop_indices[eachdim] + "\"," + legal_loop_indices[eachdim] + ","
        gufunc_txt += ")\n"

    # Add the sentinel assignment so that we can find the loop body position
    # in the IR.
    for indent in range(parfor_dim + 1):
        gufunc_txt += "    "
    gufunc_txt += sentinel_name + " = 0\n"
    # Add assignments of reduction variables (for returning the value)
    for arr, var in zip(parfor_redarrs, parfor_redvars):
        if config.DEBUG_ARRAY_OPT_RUNTIME:
            gufunc_txt += "    print(\"final reduction value\",ParallelAcceleratorGufuncThreadId," + param_dict[var] + ")\n"
            gufunc_txt += "    print(\"final reduction array\",ParallelAcceleratorGufuncThreadId," + param_dict[arr] + ")\n"
        # After the gufunc loops, copy the accumulated temp value back to reduction array.
        gufunc_txt += "    " + param_dict[arr] + \
            "[" + gufunc_thread_id_var + "] = " + param_dict[var] + "\n"
    gufunc_txt += "    return None\n"

    if config.DEBUG_ARRAY_OPT:
        print("gufunc_txt = ", type(gufunc_txt), "\n", gufunc_txt)
        print("globls:", globls, type(globls))
    # Force gufunc outline into existence.
    locls = {}
    exec(gufunc_txt, globls, locls)
    gufunc_func = locls[gufunc_name]

    if config.DEBUG_ARRAY_OPT:
        print("gufunc_func = ", type(gufunc_func), "\n", gufunc_func)
    # Get the IR for the gufunc outline.
    gufunc_ir = compiler.run_frontend(gufunc_func)
    if config.DEBUG_ARRAY_OPT:
        print("gufunc_ir dump ", type(gufunc_ir))
        gufunc_ir.dump()
        print("loop_body dump ", type(loop_body))
        _print_body(loop_body)

    # rename all variables in gufunc_ir afresh
    var_table = get_name_var_table(gufunc_ir.blocks)
    new_var_dict = {}
    reserved_names = [sentinel_name] + \
        list(param_dict.values()) + legal_loop_indices
    for name, var in var_table.items():
        if not (name in reserved_names):
            new_var_dict[name] = parfor.init_block.scope.redefine(name, loc).name
    replace_var_names(gufunc_ir.blocks, new_var_dict)
    if config.DEBUG_ARRAY_OPT:
        print("gufunc_ir dump after renaming ")
        gufunc_ir.dump()
    gufunc_param_types = [types.npytypes.Array(
            index_var_typ, 1, "C")] + param_types
    if config.DEBUG_ARRAY_OPT:
        print(
            "gufunc_param_types = ",
            type(gufunc_param_types),
            "\n",
            gufunc_param_types)

    gufunc_stub_last_label = find_max_label(gufunc_ir.blocks) + 1

    # Add gufunc stub last label to each parfor.loop_body label to prevent
    # label conflicts.
    loop_body = add_offset_to_labels(loop_body, gufunc_stub_last_label)
    # new label for splitting sentinel block
    new_label = find_max_label(loop_body) + 1

    # If enabled, add a print statement after every assignment.
    if config.DEBUG_ARRAY_OPT_RUNTIME:
        for label, block in loop_body.items():
            new_block = block.copy()
            new_block.clear()
            loc = block.loc
            scope = block.scope
            for inst in block.body:
                new_block.append(inst)
                # Append print after assignment
                if isinstance(inst, ir.Assign):
                    # Only apply to numbers
                    if typemap[inst.target.name] not in types.number_domain:
                        continue

                    # Make constant string
                    strval = "{} =".format(inst.target.name)
                    strconsttyp = types.StringLiteral(strval)

                    lhs = scope.redefine("str_const", loc)
                    # lhs = ir.Var(scope, mk_unique_var("str_const"), loc)
                    assign_lhs = ir.Assign(value=ir.Const(value=strval, loc=loc),
                                           target=lhs, loc=loc)
                    typemap[lhs.name] = strconsttyp
                    new_block.append(assign_lhs)

                    # Make print node
                    print_node = ir.Print(args=[lhs, inst.target], vararg=None, loc=loc)
                    new_block.append(print_node)
                    sig = numba.core.typing.signature(types.none,
                                           typemap[lhs.name],
                                           typemap[inst.target.name])
                    lowerer.fndesc.calltypes[print_node] = sig
            loop_body[label] = new_block

    if config.DEBUG_ARRAY_OPT:
        print("parfor loop body")
        _print_body(loop_body)

    wrapped_blocks = wrap_loop_body(loop_body)
    hoisted, not_hoisted = hoist(parfor_params, loop_body, typemap, wrapped_blocks)
    start_block = gufunc_ir.blocks[min(gufunc_ir.blocks.keys())]
    start_block.body = start_block.body[:-1] + hoisted + [start_block.body[-1]]
    unwrap_loop_body(loop_body)

    # store hoisted into diagnostics
    diagnostics = lowerer.metadata['parfor_diagnostics']
    diagnostics.hoist_info[parfor.id] = {'hoisted': hoisted,
                                         'not_hoisted': not_hoisted}

    if config.DEBUG_ARRAY_OPT:
        print("After hoisting")
        _print_body(loop_body)

    # Search all the block in the gufunc outline for the sentinel assignment.
    for label, block in gufunc_ir.blocks.items():
        for i, inst in enumerate(block.body):
            if isinstance(
                    inst,
                    ir.Assign) and inst.target.name == sentinel_name:
                # We found the sentinel assignment.
                loc = inst.loc
                scope = block.scope
                # split block across __sentinel__
                # A new block is allocated for the statements prior to the sentinel
                # but the new block maintains the current block label.
                prev_block = ir.Block(scope, loc)
                prev_block.body = block.body[:i]
                # The current block is used for statements after the sentinel.
                block.body = block.body[i + 1:]
                # But the current block gets a new label.
                body_first_label = min(loop_body.keys())

                # The previous block jumps to the minimum labelled block of the
                # parfor body.
                prev_block.append(ir.Jump(body_first_label, loc))
                # Add all the parfor loop body blocks to the gufunc function's
                # IR.
                for (l, b) in loop_body.items():
                    gufunc_ir.blocks[l] = transfer_scope(b, scope)
                body_last_label = max(loop_body.keys())
                gufunc_ir.blocks[new_label] = block
                gufunc_ir.blocks[label] = prev_block
                # Add a jump from the last parfor body block to the block containing
                # statements after the sentinel.
                gufunc_ir.blocks[body_last_label].append(
                    ir.Jump(new_label, loc))
                break
        else:
            continue
        break

    if config.DEBUG_ARRAY_OPT:
        print("gufunc_ir last dump before renaming")
        gufunc_ir.dump()

    gufunc_ir.blocks = rename_labels(gufunc_ir.blocks)
    remove_dels(gufunc_ir.blocks)

    if config.DEBUG_ARRAY_OPT:
        print("gufunc_ir last dump")
        gufunc_ir.dump()
        print("flags", flags)
        print("typemap", typemap)

    old_alias = flags.noalias
    if not has_aliases:
        if config.DEBUG_ARRAY_OPT:
            print("No aliases found so adding noalias flag.")
        flags.noalias = True

    fixup_var_define_in_scope(gufunc_ir.blocks)

    class ParforGufuncCompiler(compiler.CompilerBase):
        def define_pipelines(self):
            from numba.core.compiler_machinery import PassManager
            dpb = compiler.DefaultPassBuilder
            pm = PassManager("full_parfor_gufunc")
            parfor_gufunc_passes = dpb.define_parfor_gufunc_pipeline(self.state)
            pm.passes.extend(parfor_gufunc_passes.passes)
            lowering_passes = dpb.define_parfor_gufunc_nopython_lowering_pipeline(self.state)
            pm.passes.extend(lowering_passes.passes)

            pm.finalize()
            return [pm]

    kernel_func = compiler.compile_ir(
        typingctx,
        targetctx,
        gufunc_ir,
        gufunc_param_types,
        types.none,
        flags,
        locals,
        pipeline_class=ParforGufuncCompiler)

    flags.noalias = old_alias

    kernel_sig = signature(types.none, *gufunc_param_types)
    if config.DEBUG_ARRAY_OPT:
        print("finished create_gufunc_for_parfor_body. kernel_sig = ", kernel_sig)

    return kernel_func, parfor_args, kernel_sig, func_arg_types, expanded_name_to_tuple_var

def replace_var_with_array_in_block(vars, block, typemap, calltypes):
    new_block = []
    for inst in block.body:
        if isinstance(inst, ir.Assign) and inst.target.name in vars:
            loc = inst.loc
            scope = inst.target.scope

            const_node = ir.Const(0, loc)
            const_var = scope.redefine("$const_ind_0", loc)
            typemap[const_var.name] = types.uintp
            const_assign = ir.Assign(const_node, const_var, loc)
            new_block.append(const_assign)

            val_var = scope.redefine("$val", loc)
            typemap[val_var.name] = typemap[inst.target.name]
            new_block.append(ir.Assign(inst.value, val_var, loc))
            setitem_node = ir.SetItem(inst.target, const_var, val_var, loc)
            calltypes[setitem_node] = signature(
                types.none, types.npytypes.Array(typemap[inst.target.name], 1, "C"), types.intp, typemap[inst.target.name])
            new_block.append(setitem_node)
            continue
        elif isinstance(inst, parfor.Parfor):
            replace_var_with_array_internal(vars, {0: inst.init_block}, typemap, calltypes)
            replace_var_with_array_internal(vars, inst.loop_body, typemap, calltypes)

        new_block.append(inst)
    return new_block

def replace_var_with_array_internal(vars, loop_body, typemap, calltypes):
    for label, block in loop_body.items():
        block.body = replace_var_with_array_in_block(vars, block, typemap, calltypes)

def replace_var_with_array(vars, loop_body, typemap, calltypes):
    replace_var_with_array_internal(vars, loop_body, typemap, calltypes)
    for v in vars:
        el_typ = typemap[v]
        typemap.pop(v, None)
        typemap[v] = types.npytypes.Array(el_typ, 1, "C")

def call_parallel_gufunc(lowerer, cres, gu_signature, outer_sig, expr_args, expr_arg_types,
                         loop_ranges, redvars, reddict, redarrdict, init_block, index_var_typ, races,
                         exp_name_to_tuple_var):
    '''
    Adds the call to the gufunc function from the main function.
    '''
    context = lowerer.context
    builder = lowerer.builder

    from numba.np.ufunc.parallel import (build_gufunc_wrapper,
                           _launch_threads)

    if config.DEBUG_ARRAY_OPT:
        print("make_parallel_loop")
        print("outer_sig = ", outer_sig.args, outer_sig.return_type,
              outer_sig.recvr, outer_sig.pysig)
        print("loop_ranges = ", loop_ranges)
        print("expr_args", expr_args)
        print("expr_arg_types", expr_arg_types)
        print("gu_signature", gu_signature)

    # Build the wrapper for GUFunc
    args, return_type = sigutils.normalize_signature(outer_sig)
    llvm_func = cres.library.get_function(cres.fndesc.llvm_func_name)
    sin, sout = gu_signature

    # These are necessary for build_gufunc_wrapper to find external symbols
    _launch_threads()

    info = build_gufunc_wrapper(llvm_func, cres, sin, sout,
                                cache=False, is_parfors=True)
    wrapper_name = info.name
    cres.library._ensure_finalized()

    if config.DEBUG_ARRAY_OPT:
        print("parallel function = ", wrapper_name, cres)

    # loadvars for loop_ranges
    def load_range(v):
        if isinstance(v, ir.Var):
            return lowerer.loadvar(v.name)
        else:
            return context.get_constant(types.uintp, v)

    num_dim = len(loop_ranges)
    for i in range(num_dim):
        start, stop, step = loop_ranges[i]
        start = load_range(start)
        stop = load_range(stop)
        assert(step == 1)  # We do not support loop steps other than 1
        step = load_range(step)
        loop_ranges[i] = (start, stop, step)

        if config.DEBUG_ARRAY_OPT:
            print("call_parallel_gufunc loop_ranges[{}] = ".format(i), start,
                  stop, step)
            cgutils.printf(builder, "loop range[{}]: %d %d (%d)\n".format(i),
                           start, stop, step)

    # Commonly used LLVM types and constants
    byte_t = llvmlite.ir.IntType(8)
    byte_ptr_t = llvmlite.ir.PointerType(byte_t)
    byte_ptr_ptr_t = llvmlite.ir.PointerType(byte_ptr_t)
    intp_t = context.get_value_type(types.intp)
    uintp_t = context.get_value_type(types.uintp)
    intp_ptr_t = llvmlite.ir.PointerType(intp_t)
    intp_ptr_ptr_t = llvmlite.ir.PointerType(intp_ptr_t)
    uintp_ptr_t = llvmlite.ir.PointerType(uintp_t)
    uintp_ptr_ptr_t = llvmlite.ir.PointerType(uintp_ptr_t)
    zero = context.get_constant(types.uintp, 0)
    one = context.get_constant(types.uintp, 1)
    one_type = one.type
    sizeof_intp = context.get_abi_sizeof(intp_t)

    # Prepare sched, first pop it out of expr_args, outer_sig, and gu_signature
    expr_args.pop(0)
    sched_sig = sin.pop(0)

    if config.DEBUG_ARRAY_OPT:
        print("Parfor has potentially negative start", index_var_typ.signed)

    if index_var_typ.signed:
        sched_type = intp_t
        sched_ptr_type = intp_ptr_t
        sched_ptr_ptr_type = intp_ptr_ptr_t
    else:
        sched_type = uintp_t
        sched_ptr_type = uintp_ptr_t
        sched_ptr_ptr_type = uintp_ptr_ptr_t

    # Call do_scheduling with appropriate arguments
    dim_starts = cgutils.alloca_once(
        builder, sched_type, size=context.get_constant(
            types.uintp, num_dim), name="dim_starts")
    dim_stops = cgutils.alloca_once(
        builder, sched_type, size=context.get_constant(
            types.uintp, num_dim), name="dim_stops")
    for i in range(num_dim):
        start, stop, step = loop_ranges[i]
        if start.type != one_type:
            start = builder.sext(start, one_type)
        if stop.type != one_type:
            stop = builder.sext(stop, one_type)
        if step.type != one_type:
            step = builder.sext(step, one_type)
        # substract 1 because do-scheduling takes inclusive ranges
        stop = builder.sub(stop, one)
        builder.store(
            start, builder.gep(
                dim_starts, [
                    context.get_constant(
                        types.uintp, i)]))
        builder.store(stop, builder.gep(dim_stops,
                                        [context.get_constant(types.uintp, i)]))

    # Prepare to call get/set parallel_chunksize and get the number of threads.
    get_chunksize = cgutils.get_or_insert_function(
        builder.module,
        llvmlite.ir.FunctionType(uintp_t, []),
        name="get_parallel_chunksize")

    set_chunksize = cgutils.get_or_insert_function(
        builder.module,
        llvmlite.ir.FunctionType(llvmlite.ir.VoidType(), [uintp_t]),
        name="set_parallel_chunksize")

    get_num_threads = cgutils.get_or_insert_function(
        builder.module,
        llvmlite.ir.FunctionType(llvmlite.ir.IntType(types.intp.bitwidth), []),
        "get_num_threads")

    # Get the current number of threads.
    num_threads = builder.call(get_num_threads, [])
    # Get the current chunksize so we can use it and restore the value later.
    current_chunksize = builder.call(get_chunksize, [])

    with cgutils.if_unlikely(builder, builder.icmp_signed('<=', num_threads,
                                                  num_threads.type(0))):
        cgutils.printf(builder, "num_threads: %d\n", num_threads)
        context.call_conv.return_user_exc(builder, RuntimeError,
                                                  ("Invalid number of threads. "
                                                   "This likely indicates a bug in Numba.",))

    # Call get_sched_size from gufunc_scheduler.cpp that incorporates the size of the work,
    # the number of threads and the selected chunk size.  This will tell us how many entries
    # in the schedule we will need.
    get_sched_size_fnty = llvmlite.ir.FunctionType(uintp_t, [uintp_t, uintp_t, intp_ptr_t, intp_ptr_t])
    get_sched_size = cgutils.get_or_insert_function(
        builder.module,
        get_sched_size_fnty,
        name="get_sched_size")
    num_divisions = builder.call(get_sched_size, [num_threads,
                                                  context.get_constant(types.uintp, num_dim),
                                                  dim_starts,
                                                  dim_stops])
    # Set the chunksize to zero so that any nested calls get the default chunk size behavior.
    builder.call(set_chunksize, [zero])

    # Each entry in the schedule is 2 times the number of dimensions long.
    multiplier = context.get_constant(types.uintp, num_dim * 2)
    # Compute the total number of entries in the schedule.
    sched_size = builder.mul(num_divisions, multiplier)

    # Prepare to dynamically allocate memory to hold the schedule.
    alloc_sched_fnty = llvmlite.ir.FunctionType(sched_ptr_type, [uintp_t])
    alloc_sched_func = cgutils.get_or_insert_function(
        builder.module,
        alloc_sched_fnty,
        name="allocate_sched")
    # Call gufunc_scheduler.cpp to allocate the schedule.
    # This may or may not do pooling.
    alloc_space = builder.call(alloc_sched_func, [sched_size])
    # Allocate a slot in the entry block to store the schedule pointer.
    sched = cgutils.alloca_once(builder, sched_ptr_type)
    # Store the schedule pointer into that slot.
    builder.store(alloc_space, sched)

    debug_flag = 1 if config.DEBUG_ARRAY_OPT else 0
    scheduling_fnty = llvmlite.ir.FunctionType(
        intp_ptr_t, [uintp_t, intp_ptr_t, intp_ptr_t, uintp_t, sched_ptr_type, intp_t])
    if index_var_typ.signed:
        do_scheduling = cgutils.get_or_insert_function(builder.module,
                                                       scheduling_fnty,
                                                       name="do_scheduling_signed")
    else:
        do_scheduling = cgutils.get_or_insert_function(builder.module,
                                                       scheduling_fnty,
                                                       name="do_scheduling_unsigned")

    # Call the scheduling routine that decides how to break up the work.
    builder.call(
        do_scheduling, [
            context.get_constant(
                types.uintp, num_dim), dim_starts, dim_stops, num_divisions,
            builder.load(sched), context.get_constant(
                    types.intp, debug_flag)])

    # Get the LLVM vars for the Numba IR reduction array vars.
    redarrs = [lowerer.loadvar(redarrdict[x].name) for x in redvars]

    nredvars = len(redvars)
    ninouts = len(expr_args) - nredvars

    def load_potential_tuple_var(x):
        """Given a variable name, if that variable is not a new name
           introduced as the extracted part of a tuple then just return
           the variable loaded from its name.  However, if the variable
           does represent part of a tuple, as recognized by the name of
           the variable being present in the exp_name_to_tuple_var dict,
           then we load the original tuple var instead that we get from
           the dict and then extract the corresponding element of the
           tuple, also stored and returned to use in the dict (i.e., offset).
        """
        if x in exp_name_to_tuple_var:
            orig_tup, offset = exp_name_to_tuple_var[x]
            tup_var = lowerer.loadvar(orig_tup)
            res = builder.extract_value(tup_var, offset)
            return res
        else:
            return lowerer.loadvar(x)

    # ----------------------------------------------------------------------------
    # Prepare arguments: args, shapes, steps, data
    all_args = [load_potential_tuple_var(x) for x in expr_args[:ninouts]] + redarrs
    num_args = len(all_args)
    num_inps = len(sin) + 1
    args = cgutils.alloca_once(
        builder,
        byte_ptr_t,
        size=context.get_constant(
            types.intp,
            1 + num_args),
        name="pargs")
    array_strides = []
    # sched goes first
    builder.store(builder.bitcast(builder.load(sched), byte_ptr_t), args)
    array_strides.append(context.get_constant(types.intp, sizeof_intp))
    rv_to_arg_dict = {}
    # followed by other arguments
    for i in range(num_args):
        arg = all_args[i]
        var = expr_args[i]
        aty = expr_arg_types[i]
        dst = builder.gep(args, [context.get_constant(types.intp, i + 1)])
        if i >= ninouts:  # reduction variables
            ary = context.make_array(aty)(context, builder, arg)
            strides = cgutils.unpack_tuple(builder, ary.strides, aty.ndim)
            # Start from 1 because we skip the first dimension of length num_threads just like sched.
            for j in range(len(strides)):
                array_strides.append(strides[j])
            builder.store(builder.bitcast(ary.data, byte_ptr_t), dst)
        elif isinstance(aty, types.ArrayCompatible):
            if var in races:
                typ = (context.get_data_type(aty.dtype)
                       if aty.dtype != types.boolean
                       else llvmlite.ir.IntType(1))

                rv_arg = cgutils.alloca_once(builder, typ)
                builder.store(arg, rv_arg)
                builder.store(builder.bitcast(rv_arg, byte_ptr_t), dst)
                rv_to_arg_dict[var] = (arg, rv_arg)

                array_strides.append(context.get_constant(types.intp, context.get_abi_sizeof(typ)))
            else:
                ary = context.make_array(aty)(context, builder, arg)
                strides = cgutils.unpack_tuple(builder, ary.strides, aty.ndim)
                for j in range(len(strides)):
                    array_strides.append(strides[j])
                builder.store(builder.bitcast(ary.data, byte_ptr_t), dst)
        else:
            if i < num_inps:
                # Scalar input, need to store the value in an array of size 1
                if isinstance(aty, types.Optional):
                    # Unpack optional type
                    unpacked_aty = aty.type
                    arg = context.cast(builder, arg, aty, unpacked_aty)
                else:
                    unpacked_aty = aty
                typ = (context.get_data_type(unpacked_aty)
                       if not isinstance(unpacked_aty, types.Boolean)
                       else llvmlite.ir.IntType(1))
                ptr = cgutils.alloca_once(builder, typ)
                builder.store(arg, ptr)
            else:
                # Scalar output, must allocate
                typ = (context.get_data_type(aty)
                       if not isinstance(aty, types.Boolean)
                       else llvmlite.ir.IntType(1))
                ptr = cgutils.alloca_once(builder, typ)
            builder.store(builder.bitcast(ptr, byte_ptr_t), dst)

    # ----------------------------------------------------------------------------
    # Next, we prepare the individual dimension info recorded in gu_signature
    sig_dim_dict = {}
    occurrences = []
    occurrences = [sched_sig[0]]
    sig_dim_dict[sched_sig[0]] = context.get_constant(types.intp, 2 * num_dim)
    assert len(expr_args) == len(all_args)
    assert len(expr_args) == len(expr_arg_types)
    assert len(expr_args) == len(sin + sout)
    assert len(expr_args) == len(outer_sig.args[1:])
    for var, arg, aty, gu_sig in zip(expr_args, all_args,
                                     expr_arg_types, sin + sout):
        if isinstance(aty, types.npytypes.Array):
            i = aty.ndim - len(gu_sig)
        else:
            i = 0
        if config.DEBUG_ARRAY_OPT:
            print("var =", var, "gu_sig =", gu_sig, "type =", aty, "i =", i)

        for dim_sym in gu_sig:
            if config.DEBUG_ARRAY_OPT:
                print("var = ", var, " type = ", aty)
            if var in races:
                sig_dim_dict[dim_sym] = context.get_constant(types.intp, 1)
            else:
                ary = context.make_array(aty)(context, builder, arg)
                shapes = cgutils.unpack_tuple(builder, ary.shape, aty.ndim)
                sig_dim_dict[dim_sym] = shapes[i]

            if not (dim_sym in occurrences):
                if config.DEBUG_ARRAY_OPT:
                    print("dim_sym = ", dim_sym, ", i = ", i)
                    cgutils.printf(builder, dim_sym + " = %d\n", sig_dim_dict[dim_sym])
                occurrences.append(dim_sym)
            i = i + 1

    # ----------------------------------------------------------------------------
    # Prepare shapes, which is a single number (outer loop size), followed by
    # the size of individual shape variables.
    nshapes = len(sig_dim_dict) + 1
    shapes = cgutils.alloca_once(builder, intp_t, size=nshapes, name="pshape")
    # For now, outer loop size is the same as number of threads
    builder.store(num_divisions, shapes)
    # Individual shape variables go next
    i = 1
    for dim_sym in occurrences:
        if config.DEBUG_ARRAY_OPT:
            cgutils.printf(builder, dim_sym + " = %d\n", sig_dim_dict[dim_sym])
        builder.store(
            sig_dim_dict[dim_sym], builder.gep(
                shapes, [
                    context.get_constant(
                        types.intp, i)]))
        i = i + 1

    # ----------------------------------------------------------------------------
    # Prepare steps for each argument. Note that all steps are counted in
    # bytes.
    num_steps = num_args + 1 + len(array_strides)
    steps = cgutils.alloca_once(
        builder, intp_t, size=context.get_constant(
            types.intp, num_steps), name="psteps")
    # First goes the step size for sched, which is 2 * num_dim
    builder.store(context.get_constant(types.intp, 2 * num_dim * sizeof_intp),
                  steps)
    # The steps for all others are 0, except for reduction results.
    for i in range(num_args):
        # steps are strides from one thread to the next
        stepsize = zero

        dst = builder.gep(steps, [context.get_constant(types.intp, 1 + i)])
        builder.store(stepsize, dst)
    for j in range(len(array_strides)):
        dst = builder.gep(
            steps, [
                context.get_constant(
                    types.intp, 1 + num_args + j)])
        builder.store(array_strides[j], dst)

    # ----------------------------------------------------------------------------
    # prepare data
    data = cgutils.get_null_value(byte_ptr_t)

    fnty = llvmlite.ir.FunctionType(llvmlite.ir.VoidType(),
                                    [byte_ptr_ptr_t, intp_ptr_t,
                                     intp_ptr_t, byte_ptr_t])

    fn = cgutils.get_or_insert_function(builder.module, fnty, wrapper_name)
    context.active_code_library.add_linking_library(info.library)

    if config.DEBUG_ARRAY_OPT:
        cgutils.printf(builder, "before calling kernel %p\n", fn)
    builder.call(fn, [args, shapes, steps, data])
    if config.DEBUG_ARRAY_OPT:
        cgutils.printf(builder, "after calling kernel %p\n", fn)

    builder.call(set_chunksize, [current_chunksize])

    # Deallocate the schedule's memory.
    dealloc_sched_fnty = llvmlite.ir.FunctionType(llvmlite.ir.VoidType(), [sched_ptr_type])
    dealloc_sched_func = cgutils.get_or_insert_function(
        builder.module,
        dealloc_sched_fnty,
        name="deallocate_sched")
    builder.call(dealloc_sched_func, [builder.load(sched)])

    for k, v in rv_to_arg_dict.items():
        arg, rv_arg = v
        only_elem_ptr = builder.gep(rv_arg, [context.get_constant(types.intp, 0)])
        builder.store(builder.load(only_elem_ptr), lowerer.getvar(k))

    context.active_code_library.add_linking_library(cres.library)
