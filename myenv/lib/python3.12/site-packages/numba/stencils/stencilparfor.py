#
# Copyright (c) 2017 Intel Corporation
# SPDX-License-Identifier: BSD-2-Clause
#

import numbers
import copy
import types as pytypes
from operator import add
import operator

import numpy as np

import numba.parfors.parfor
from numba.core import types, ir, rewrites, config, ir_utils
from numba.core.typing.templates import infer_global, AbstractTemplate
from numba.core.typing import signature
from numba.core import  utils, typing
from numba.core.ir_utils import (get_call_table, mk_unique_var,
                            compile_to_numba_ir, replace_arg_nodes, guard,
                            find_callname, require, find_const, GuardException)
from numba.core.errors import NumbaValueError
from numba.core.utils import OPERATORS_TO_BUILTINS
from numba.np import numpy_support


def _compute_last_ind(dim_size, index_const):
    if index_const > 0:
        return dim_size - index_const
    else:
        return dim_size

class StencilPass(object):
    def __init__(self, func_ir, typemap, calltypes, array_analysis, typingctx,
                 targetctx, flags):
        self.func_ir = func_ir
        self.typemap = typemap
        self.calltypes = calltypes
        self.array_analysis = array_analysis
        self.typingctx = typingctx
        self.targetctx = targetctx
        self.flags = flags

    def run(self):
        """ Finds all calls to StencilFuncs in the IR and converts them to parfor.
        """
        from numba.stencils.stencil import StencilFunc

        # Get all the calls in the function IR.
        call_table, _ = get_call_table(self.func_ir.blocks)
        stencil_calls = []
        stencil_dict = {}
        for call_varname, call_list in call_table.items():
            for one_call in call_list:
                if isinstance(one_call, StencilFunc):
                    # Remember all calls to StencilFuncs.
                    stencil_calls.append(call_varname)
                    stencil_dict[call_varname] = one_call
        if not stencil_calls:
            return  # return early if no stencil calls found

        # find and transform stencil calls
        for label, block in self.func_ir.blocks.items():
            for i, stmt in reversed(list(enumerate(block.body))):
                # Found a call to a StencilFunc.
                if (isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Expr)
                        and stmt.value.op == 'call'
                        and stmt.value.func.name in stencil_calls):
                    kws = dict(stmt.value.kws)
                    # Create dictionary of input argument number to
                    # the argument itself.
                    input_dict = {i: stmt.value.args[i] for i in
                                                    range(len(stmt.value.args))}
                    in_args = stmt.value.args
                    arg_typemap = tuple(self.typemap[i.name] for i in in_args)
                    for arg_type in arg_typemap:
                        if isinstance(arg_type, types.BaseTuple):
                            raise ValueError("Tuple parameters not supported " \
                                "for stencil kernels in parallel=True mode.")

                    out_arr = kws.get('out')

                    # Get the StencilFunc object corresponding to this call.
                    sf = stencil_dict[stmt.value.func.name]
                    stencil_ir, rt, arg_to_arr_dict = get_stencil_ir(sf,
                            self.typingctx, arg_typemap,
                            block.scope, block.loc, input_dict,
                            self.typemap, self.calltypes)
                    index_offsets = sf.options.get('index_offsets', None)
                    gen_nodes = self._mk_stencil_parfor(label, in_args, out_arr,
                            stencil_ir, index_offsets, stmt.target, rt, sf,
                            arg_to_arr_dict)
                    block.body = block.body[:i] + gen_nodes + block.body[i+1:]
                # Found a call to a stencil via numba.stencil().
                elif (isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Expr)
                        and stmt.value.op == 'call'
                        and guard(find_callname, self.func_ir, stmt.value)
                                    == ('stencil', 'numba')):
                    # remove dummy stencil() call
                    stmt.value = ir.Const(0, stmt.loc)

    def replace_return_with_setitem(self, blocks, exit_value_var,
                                    parfor_body_exit_label):
        """
        Find return statements in the IR and replace them with a SetItem
        call of the value "returned" by the kernel into the result array.
        Returns the block labels that contained return statements.
        """
        for label, block in blocks.items():
            scope = block.scope
            loc = block.loc
            new_body = []
            for stmt in block.body:
                if isinstance(stmt, ir.Return):
                    # previous stmt should have been a cast
                    prev_stmt = new_body.pop()
                    assert (isinstance(prev_stmt, ir.Assign)
                        and isinstance(prev_stmt.value, ir.Expr)
                        and prev_stmt.value.op == 'cast')

                    new_body.append(ir.Assign(prev_stmt.value.value, exit_value_var, loc))
                    new_body.append(ir.Jump(parfor_body_exit_label, loc))
                else:
                    new_body.append(stmt)
            block.body = new_body

    def _mk_stencil_parfor(self, label, in_args, out_arr, stencil_ir,
                           index_offsets, target, return_type, stencil_func,
                           arg_to_arr_dict):
        """ Converts a set of stencil kernel blocks to a parfor.
        """
        gen_nodes = []
        stencil_blocks = stencil_ir.blocks

        if config.DEBUG_ARRAY_OPT >= 1:
            print("_mk_stencil_parfor", label, in_args, out_arr, index_offsets,
                   return_type, stencil_func, stencil_blocks)
            ir_utils.dump_blocks(stencil_blocks)

        in_arr = in_args[0]
        # run copy propagate to replace in_args copies (e.g. a = A)
        in_arr_typ = self.typemap[in_arr.name]
        in_cps, out_cps = ir_utils.copy_propagate(stencil_blocks, self.typemap)
        name_var_table = ir_utils.get_name_var_table(stencil_blocks)

        ir_utils.apply_copy_propagate(
            stencil_blocks,
            in_cps,
            name_var_table,
            self.typemap,
            self.calltypes)
        if config.DEBUG_ARRAY_OPT >= 1:
            print("stencil_blocks after copy_propagate")
            ir_utils.dump_blocks(stencil_blocks)
        ir_utils.remove_dead(stencil_blocks, self.func_ir.arg_names, stencil_ir,
                             self.typemap)
        if config.DEBUG_ARRAY_OPT >= 1:
            print("stencil_blocks after removing dead code")
            ir_utils.dump_blocks(stencil_blocks)

        # create parfor vars
        ndims = self.typemap[in_arr.name].ndim
        scope = in_arr.scope
        loc = in_arr.loc
        parfor_vars = []
        for i in range(ndims):
            parfor_var = ir.Var(scope, mk_unique_var(
                "$parfor_index_var"), loc)
            self.typemap[parfor_var.name] = types.intp
            parfor_vars.append(parfor_var)

        start_lengths, end_lengths = self._replace_stencil_accesses(
             stencil_ir, parfor_vars, in_args, index_offsets, stencil_func,
             arg_to_arr_dict)

        if config.DEBUG_ARRAY_OPT >= 1:
            print("stencil_blocks after replace stencil accesses")
            print("start_lengths:", start_lengths)
            print("end_lengths:", end_lengths)
            ir_utils.dump_blocks(stencil_blocks)

        # create parfor loop nests
        loopnests = []
        equiv_set = self.array_analysis.get_equiv_set(label)
        in_arr_dim_sizes = equiv_set.get_shape(in_arr)

        assert ndims == len(in_arr_dim_sizes)
        start_inds = []
        last_inds = []
        for i in range(ndims):
            last_ind = self._get_stencil_last_ind(in_arr_dim_sizes[i],
                                        end_lengths[i], gen_nodes, scope, loc)
            start_ind = self._get_stencil_start_ind(
                                        start_lengths[i], gen_nodes, scope, loc)
            start_inds.append(start_ind)
            last_inds.append(last_ind)
            # start from stencil size to avoid invalid array access
            loopnests.append(numba.parfors.parfor.LoopNest(parfor_vars[i],
                                start_ind, last_ind, 1))

        # We have to guarantee that the exit block has maximum label and that
        # there's only one exit block for the parfor body.
        # So, all return statements will change to jump to the parfor exit block.
        parfor_body_exit_label = max(stencil_blocks.keys()) + 1
        stencil_blocks[parfor_body_exit_label] = ir.Block(scope, loc)
        exit_value_var = ir.Var(scope, mk_unique_var("$parfor_exit_value"), loc)
        self.typemap[exit_value_var.name] = return_type.dtype

        # create parfor index var
        for_replacing_ret = []
        if ndims == 1:
            parfor_ind_var = parfor_vars[0]
        else:
            parfor_ind_var = ir.Var(scope, mk_unique_var(
                "$parfor_index_tuple_var"), loc)
            self.typemap[parfor_ind_var.name] = types.containers.UniTuple(
                types.intp, ndims)
            tuple_call = ir.Expr.build_tuple(parfor_vars, loc)
            tuple_assign = ir.Assign(tuple_call, parfor_ind_var, loc)
            for_replacing_ret.append(tuple_assign)

        if config.DEBUG_ARRAY_OPT >= 1:
            print("stencil_blocks after creating parfor index var")
            ir_utils.dump_blocks(stencil_blocks)

        # empty init block
        init_block = ir.Block(scope, loc)
        if out_arr is None:
            in_arr_typ = self.typemap[in_arr.name]

            shape_name = ir_utils.mk_unique_var("in_arr_shape")
            shape_var = ir.Var(scope, shape_name, loc)
            shape_getattr = ir.Expr.getattr(in_arr, "shape", loc)
            self.typemap[shape_name] = types.containers.UniTuple(types.intp,
                                                               in_arr_typ.ndim)
            init_block.body.extend([ir.Assign(shape_getattr, shape_var, loc)])

            zero_name = ir_utils.mk_unique_var("zero_val")
            zero_var = ir.Var(scope, zero_name, loc)
            if "cval" in stencil_func.options:
                cval = stencil_func.options["cval"]
                # TODO: Loosen this restriction to adhere to casting rules.
                cval_ty = typing.typeof.typeof(cval)
                if not self.typingctx.can_convert(cval_ty, return_type.dtype):
                    raise ValueError("cval type does not match stencil return type.")

                temp2 = return_type.dtype(cval)
            else:
                temp2 = return_type.dtype(0)
            full_const = ir.Const(temp2, loc)
            self.typemap[zero_name] = return_type.dtype
            init_block.body.extend([ir.Assign(full_const, zero_var, loc)])

            so_name = ir_utils.mk_unique_var("stencil_output")
            out_arr = ir.Var(scope, so_name, loc)
            self.typemap[out_arr.name] = numba.core.types.npytypes.Array(
                                                           return_type.dtype,
                                                           in_arr_typ.ndim,
                                                           in_arr_typ.layout)
            dtype_g_np_var = ir.Var(scope, mk_unique_var("$np_g_var"), loc)
            self.typemap[dtype_g_np_var.name] = types.misc.Module(np)
            dtype_g_np = ir.Global('np', np, loc)
            dtype_g_np_assign = ir.Assign(dtype_g_np, dtype_g_np_var, loc)
            init_block.body.append(dtype_g_np_assign)

            return_type_name = numpy_support.as_dtype(
                               return_type.dtype).type.__name__
            if return_type_name == 'bool':
                return_type_name = 'bool_'
            dtype_np_attr_call = ir.Expr.getattr(dtype_g_np_var, return_type_name, loc)
            dtype_attr_var = ir.Var(scope, mk_unique_var("$np_attr_attr"), loc)
            self.typemap[dtype_attr_var.name] = types.functions.NumberClass(return_type.dtype)
            dtype_attr_assign = ir.Assign(dtype_np_attr_call, dtype_attr_var, loc)
            init_block.body.append(dtype_attr_assign)

            stmts = ir_utils.gen_np_call("empty",
                                       np.empty,
                                       out_arr,
                                       [shape_var, dtype_attr_var],
                                       self.typingctx,
                                       self.typemap,
                                       self.calltypes)
            # ------------------
            # Generate the code to fill just the border with zero_var.

            # Generate a none var to use in slicing.
            none_var = ir.Var(scope, mk_unique_var("$none_var"), loc)
            none_assign = ir.Assign(ir.Const(None, loc), none_var, loc)
            stmts.append(none_assign)
            self.typemap[none_var.name] = types.none
            # Generate a zero var to use in slicing.
            zero_index_var = ir.Var(scope, mk_unique_var("$zero_index_var"), loc)
            zero_index_assign = ir.Assign(ir.Const(0, loc), zero_index_var, loc)
            stmts.append(zero_index_assign)
            self.typemap[zero_index_var.name] = types.intp
            # Generate generic ":" slice.
            # ---- Generate var to hold slice func var.
            slice_func_var = ir.Var(scope, mk_unique_var("$slice_func_var"), loc)
            slice_fn_ty = self.typingctx.resolve_value_type(slice)
            self.typemap[slice_func_var.name] = slice_fn_ty
            slice_g = ir.Global('slice', slice, loc)
            slice_assign = ir.Assign(slice_g, slice_func_var, loc)
            stmts.append(slice_assign)
            # ---- Generate call to slice func.
            sig = self.typingctx.resolve_function_type(slice_fn_ty,
                                                       (types.none,) * 2,
                                                       {})
            slice_callexpr = ir.Expr.call(func=slice_func_var,
                                          args=(none_var, none_var),
                                          kws=(),
                                          loc=loc)
            self.calltypes[slice_callexpr] = sig
            # ---- Generate slice var
            slice_var = ir.Var(scope, mk_unique_var("$slice"), loc)
            self.typemap[slice_var.name] = types.slice2_type
            slice_assign = ir.Assign(slice_callexpr, slice_var, loc)
            stmts.append(slice_assign)

            def handle_border(slice_fn_ty,
                              dim,
                              scope,
                              loc,
                              slice_func_var,
                              stmts,
                              border_inds,
                              border_tuple_items,
                              other_arg,
                              other_first):
                # Handle the border for start or end of the index range.
                # ---- Generate call to slice func.
                sig = self.typingctx.resolve_function_type(
                    slice_fn_ty,
                    (types.intp,) * 2,
                    {})
                si = border_inds[dim]
                assert(isinstance(si, (int, ir.Var)))
                si_var = ir.Var(scope, mk_unique_var("$border_ind"), loc)
                self.typemap[si_var.name] = types.intp
                if isinstance(si, int):
                    si_assign = ir.Assign(ir.Const(si, loc), si_var, loc)
                else:
                    si_assign = ir.Assign(si, si_var, loc)
                stmts.append(si_assign)

                slice_callexpr = ir.Expr.call(
                    func=slice_func_var,
                    args=(other_arg, si_var) if other_first else (si_var, other_arg),
                    kws=(),
                    loc=loc)
                self.calltypes[slice_callexpr] = sig
                # ---- Generate slice var
                border_slice_var = ir.Var(scope, mk_unique_var("$slice"), loc)
                self.typemap[border_slice_var.name] = types.slice2_type
                slice_assign = ir.Assign(slice_callexpr, border_slice_var, loc)
                stmts.append(slice_assign)

                border_tuple_items[dim] = border_slice_var
                border_ind_var = ir.Var(scope, mk_unique_var(
                    "$border_index_tuple_var"), loc)
                self.typemap[border_ind_var.name] = types.containers.UniTuple(
                    types.slice2_type, ndims)
                tuple_call = ir.Expr.build_tuple(border_tuple_items, loc)
                tuple_assign = ir.Assign(tuple_call, border_ind_var, loc)
                stmts.append(tuple_assign)

                setitem_call = ir.SetItem(out_arr, border_ind_var, zero_var, loc)
                self.calltypes[setitem_call] = signature(
                                                types.none, self.typemap[out_arr.name],
                                                self.typemap[border_ind_var.name],
                                                self.typemap[out_arr.name].dtype
                                                )
                stmts.append(setitem_call)

            # For each dimension, add setitem to set border values.
            for dim in range(in_arr_typ.ndim):
                # First, fill all entries with ":".
                start_tuple_items = [slice_var] * in_arr_typ.ndim
                last_tuple_items = [slice_var] * in_arr_typ.ndim

                handle_border(slice_fn_ty,
                              dim,
                              scope,
                              loc,
                              slice_func_var,
                              stmts,
                              start_inds,
                              start_tuple_items,
                              zero_index_var,
                              True)
                handle_border(slice_fn_ty,
                              dim,
                              scope,
                              loc,
                              slice_func_var,
                              stmts,
                              last_inds,
                              last_tuple_items,
                              in_arr_dim_sizes[dim],
                              False)

            # ------------------

            equiv_set.insert_equiv(out_arr, in_arr_dim_sizes)
            init_block.body.extend(stmts)
        else: # out is present
            if "cval" in stencil_func.options: # do out[:] = cval
                cval = stencil_func.options["cval"]
                # TODO: Loosen this restriction to adhere to casting rules.
                cval_ty = typing.typeof.typeof(cval)
                if not self.typingctx.can_convert(cval_ty, return_type.dtype):
                    msg = "cval type does not match stencil return type."
                    raise NumbaValueError(msg)

                # get slice ref
                slice_var = ir.Var(scope, mk_unique_var("$py_g_var"), loc)
                slice_fn_ty = self.typingctx.resolve_value_type(slice)
                self.typemap[slice_var.name] = slice_fn_ty
                slice_g = ir.Global('slice', slice, loc)
                slice_assigned = ir.Assign(slice_g, slice_var, loc)
                init_block.body.append(slice_assigned)

                sig = self.typingctx.resolve_function_type(slice_fn_ty,
                                                           (types.none,) * 2,
                                                           {})

                callexpr = ir.Expr.call(func=slice_var, args=(), kws=(),
                                        loc=loc)

                self.calltypes[callexpr] = sig
                slice_inst_var = ir.Var(scope, mk_unique_var("$slice_inst"),
                                        loc)
                self.typemap[slice_inst_var.name] = types.slice2_type
                slice_assign = ir.Assign(callexpr, slice_inst_var, loc)
                init_block.body.append(slice_assign)

                # get const val for cval
                cval_const_val = ir.Const(return_type.dtype(cval), loc)
                cval_const_var = ir.Var(scope, mk_unique_var("$cval_const"),
                                            loc)
                self.typemap[cval_const_var.name] = return_type.dtype
                cval_const_assign = ir.Assign(cval_const_val,
                                              cval_const_var, loc)
                init_block.body.append(cval_const_assign)

                # do setitem on `out` array
                setitemexpr = ir.StaticSetItem(out_arr, slice(None, None),
                                               slice_inst_var, cval_const_var,
                                               loc)
                init_block.body.append(setitemexpr)
                sig = signature(types.none, self.typemap[out_arr.name],
                                self.typemap[slice_inst_var.name],
                                self.typemap[out_arr.name].dtype)
                self.calltypes[setitemexpr] = sig


        self.replace_return_with_setitem(stencil_blocks, exit_value_var,
                                         parfor_body_exit_label)

        if config.DEBUG_ARRAY_OPT >= 1:
            print("stencil_blocks after replacing return")
            ir_utils.dump_blocks(stencil_blocks)

        setitem_call = ir.SetItem(out_arr, parfor_ind_var, exit_value_var, loc)
        self.calltypes[setitem_call] = signature(
                                        types.none, self.typemap[out_arr.name],
                                        self.typemap[parfor_ind_var.name],
                                        self.typemap[out_arr.name].dtype
                                        )
        stencil_blocks[parfor_body_exit_label].body.extend(for_replacing_ret)
        stencil_blocks[parfor_body_exit_label].body.append(setitem_call)

        # simplify CFG of parfor body (exit block could be simplified often)
        # add dummy return to enable CFG
        dummy_loc = ir.Loc("stencilparfor_dummy", -1)
        ret_const_var = ir.Var(scope, mk_unique_var("$cval_const"), dummy_loc)
        cval_const_assign = ir.Assign(ir.Const(0, loc=dummy_loc), ret_const_var, dummy_loc)
        stencil_blocks[parfor_body_exit_label].body.append(cval_const_assign)

        stencil_blocks[parfor_body_exit_label].body.append(
            ir.Return(ret_const_var, dummy_loc),
        )
        stencil_blocks = ir_utils.simplify_CFG(stencil_blocks)
        stencil_blocks[max(stencil_blocks.keys())].body.pop()

        if config.DEBUG_ARRAY_OPT >= 1:
            print("stencil_blocks after adding SetItem")
            ir_utils.dump_blocks(stencil_blocks)

        pattern = ('stencil', [start_lengths, end_lengths])
        parfor = numba.parfors.parfor.Parfor(loopnests, init_block, stencil_blocks,
                                     loc, parfor_ind_var, equiv_set, pattern, self.flags)
        gen_nodes.append(parfor)
        gen_nodes.append(ir.Assign(out_arr, target, loc))
        return gen_nodes

    def _get_stencil_last_ind(self, dim_size, end_length, gen_nodes, scope,
                                                                        loc):
        last_ind = dim_size
        if end_length != 0:
            # set last index to size minus stencil size to avoid invalid
            # memory access
            index_const = ir.Var(scope, mk_unique_var("stencil_const_var"),
                                                                        loc)
            self.typemap[index_const.name] = types.intp
            if isinstance(end_length, numbers.Number):
                const_assign = ir.Assign(ir.Const(end_length, loc),
                                                        index_const, loc)
            else:
                const_assign = ir.Assign(end_length, index_const, loc)

            gen_nodes.append(const_assign)
            last_ind = ir.Var(scope, mk_unique_var("last_ind"), loc)
            self.typemap[last_ind.name] = types.intp

            g_var = ir.Var(scope, mk_unique_var("compute_last_ind_var"), loc)
            check_func = numba.njit(_compute_last_ind)
            func_typ = types.functions.Dispatcher(check_func)
            self.typemap[g_var.name] = func_typ
            g_obj = ir.Global("_compute_last_ind", check_func, loc)
            g_assign = ir.Assign(g_obj, g_var, loc)
            gen_nodes.append(g_assign)
            index_call = ir.Expr.call(g_var, [dim_size, index_const], (), loc)
            self.calltypes[index_call] = func_typ.get_call_type(
                self.typingctx, [types.intp, types.intp], {})
            index_assign = ir.Assign(index_call, last_ind, loc)
            gen_nodes.append(index_assign)

        return last_ind

    def _get_stencil_start_ind(self, start_length, gen_nodes, scope, loc):
        if isinstance(start_length, int):
            return abs(min(start_length, 0))
        def get_start_ind(s_length):
            return abs(min(s_length, 0))
        f_ir = compile_to_numba_ir(get_start_ind, {}, self.typingctx,
                                   self.targetctx, (types.intp,), self.typemap,
                                   self.calltypes)
        assert len(f_ir.blocks) == 1
        block = f_ir.blocks.popitem()[1]
        replace_arg_nodes(block, [start_length])
        gen_nodes += block.body[:-2]
        ret_var = block.body[-2].value.value
        return ret_var

    def _replace_stencil_accesses(self, stencil_ir, parfor_vars, in_args,
                                  index_offsets, stencil_func, arg_to_arr_dict):
        """ Convert relative indexing in the stencil kernel to standard indexing
            by adding the loop index variables to the corresponding dimensions
            of the array index tuples.
        """
        stencil_blocks = stencil_ir.blocks
        in_arr = in_args[0]
        in_arg_names = [x.name for x in in_args]

        if "standard_indexing" in stencil_func.options:
            for x in stencil_func.options["standard_indexing"]:
                if x not in arg_to_arr_dict:
                    raise ValueError("Standard indexing requested for an array " \
                        "name not present in the stencil kernel definition.")
            standard_indexed = [arg_to_arr_dict[x] for x in
                                     stencil_func.options["standard_indexing"]]
        else:
            standard_indexed = []

        if in_arr.name in standard_indexed:
            raise ValueError("The first argument to a stencil kernel must use " \
                "relative indexing, not standard indexing.")

        ndims = self.typemap[in_arr.name].ndim
        scope = in_arr.scope
        loc = in_arr.loc
        # replace access indices, find access lengths in each dimension
        need_to_calc_kernel = stencil_func.neighborhood is None

        # If we need to infer the kernel size then initialize the minimum and
        # maximum seen indices for each dimension to 0.  If we already have
        # the neighborhood calculated then just convert from neighborhood format
        # to the separate start and end lengths format used here.
        if need_to_calc_kernel:
            start_lengths = ndims*[0]
            end_lengths = ndims*[0]
        else:
            start_lengths = [x[0] for x in stencil_func.neighborhood]
            end_lengths   = [x[1] for x in stencil_func.neighborhood]

        # Get all the tuples defined in the stencil blocks.
        tuple_table = ir_utils.get_tuple_table(stencil_blocks)

        found_relative_index = False

        # For all blocks in the stencil kernel...
        for label, block in stencil_blocks.items():
            new_body = []
            # For all statements in those blocks...
            for stmt in block.body:
                # Reject assignments to input arrays.
                if ((isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Expr)
                        and stmt.value.op in ['setitem', 'static_setitem']
                        and stmt.value.value.name in in_arg_names) or
                   ((isinstance(stmt, ir.SetItem) or
                     isinstance(stmt, ir.StaticSetItem))
                        and stmt.target.name in in_arg_names)):
                    raise ValueError("Assignments to arrays passed to stencil kernels is not allowed.")
                # We found a getitem for some array.  If that array is an input
                # array and isn't in the list of standard indexed arrays then
                # update min and max seen indices if we are inferring the
                # kernel size and create a new tuple where the relative offsets
                # are added to loop index vars to get standard indexing.
                if (isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Expr)
                        and stmt.value.op in ['static_getitem', 'getitem']
                        and stmt.value.value.name in in_arg_names
                        and stmt.value.value.name not in standard_indexed):
                    index_list = stmt.value.index
                    # handle 1D case
                    if ndims == 1:
                        index_list = [index_list]
                    else:
                        if hasattr(index_list, 'name') and index_list.name in tuple_table:
                            index_list = tuple_table[index_list.name]
                    # indices can be inferred as constant in simple expressions
                    # like -c where c is constant
                    # handled here since this is a common stencil index pattern
                    stencil_ir._definitions = ir_utils.build_definitions(stencil_blocks)
                    index_list = [_get_const_index_expr(
                        stencil_ir, self.func_ir, v) for v in index_list]
                    if index_offsets:
                        index_list = self._add_index_offsets(index_list,
                                    list(index_offsets), new_body, scope, loc)

                    # update min and max indices
                    if need_to_calc_kernel:
                        # all indices should be integer to be able to calculate
                        # neighborhood automatically
                        if (isinstance(index_list, ir.Var) or
                            any([not isinstance(v, int) for v in index_list])):
                            raise ValueError("Variable stencil index only "
                                "possible with known neighborhood")
                        start_lengths = list(map(min, start_lengths,
                                                                    index_list))
                        end_lengths = list(map(max, end_lengths, index_list))
                        found_relative_index = True

                    # update access indices
                    index_vars = self._add_index_offsets(parfor_vars,
                                list(index_list), new_body, scope, loc)

                    # new access index tuple
                    if ndims == 1:
                        ind_var = index_vars[0]
                    else:
                        ind_var = ir.Var(scope, mk_unique_var(
                            "$parfor_index_ind_var"), loc)
                        self.typemap[ind_var.name] = types.containers.UniTuple(
                            types.intp, ndims)
                        tuple_call = ir.Expr.build_tuple(index_vars, loc)
                        tuple_assign = ir.Assign(tuple_call, ind_var, loc)
                        new_body.append(tuple_assign)

                    # getitem return type is scalar if all indices are integer
                    if all([self.typemap[v.name] == types.intp
                                                        for v in index_vars]):
                        getitem_return_typ = self.typemap[
                                                    stmt.value.value.name].dtype
                    else:
                        # getitem returns an array
                        getitem_return_typ = self.typemap[stmt.value.value.name]
                    # new getitem with the new index var
                    getitem_call = ir.Expr.getitem(stmt.value.value, ind_var,
                                                                            loc)
                    self.calltypes[getitem_call] = signature(
                        getitem_return_typ,
                        self.typemap[stmt.value.value.name],
                        self.typemap[ind_var.name])
                    stmt.value = getitem_call

                new_body.append(stmt)
            block.body = new_body
        if need_to_calc_kernel and not found_relative_index:
            raise ValueError("Stencil kernel with no accesses to " \
                "relatively indexed arrays.")

        return start_lengths, end_lengths

    def _add_index_offsets(self, index_list, index_offsets, new_body,
                           scope, loc):
        """ Does the actual work of adding loop index variables to the
            relative index constants or variables.
        """
        assert len(index_list) == len(index_offsets)

        # shortcut if all values are integer
        if all([isinstance(v, int) for v in index_list+index_offsets]):
            # add offsets in all dimensions
            return list(map(add, index_list, index_offsets))

        out_nodes = []
        index_vars = []
        for i in range(len(index_list)):
            # new_index = old_index + offset
            old_index_var = index_list[i]
            if isinstance(old_index_var, int):
                old_index_var = ir.Var(scope,
                                mk_unique_var("old_index_var"), loc)
                self.typemap[old_index_var.name] = types.intp
                const_assign = ir.Assign(ir.Const(index_list[i], loc),
                                                    old_index_var, loc)
                out_nodes.append(const_assign)

            offset_var = index_offsets[i]
            if isinstance(offset_var, int):
                offset_var = ir.Var(scope,
                                mk_unique_var("offset_var"), loc)
                self.typemap[offset_var.name] = types.intp
                const_assign = ir.Assign(ir.Const(index_offsets[i], loc),
                                                offset_var, loc)
                out_nodes.append(const_assign)

            if (isinstance(old_index_var, slice)
                    or isinstance(self.typemap[old_index_var.name],
                                    types.misc.SliceType)):
                # only one arg can be slice
                assert self.typemap[offset_var.name] == types.intp
                index_var = self._add_offset_to_slice(old_index_var, offset_var,
                                                        out_nodes, scope, loc)
                index_vars.append(index_var)
                continue

            if (isinstance(offset_var, slice)
                    or isinstance(self.typemap[offset_var.name],
                                    types.misc.SliceType)):
                # only one arg can be slice
                assert self.typemap[old_index_var.name] == types.intp
                index_var = self._add_offset_to_slice(offset_var, old_index_var,
                                                        out_nodes, scope, loc)
                index_vars.append(index_var)
                continue

            index_var = ir.Var(scope,
                            mk_unique_var("offset_stencil_index"), loc)
            self.typemap[index_var.name] = types.intp
            index_call = ir.Expr.binop(operator.add, old_index_var,
                                                offset_var, loc)
            self.calltypes[index_call] = self.typingctx.resolve_function_type(
                                         operator.add, (types.intp, types.intp), {})
            index_assign = ir.Assign(index_call, index_var, loc)
            out_nodes.append(index_assign)
            index_vars.append(index_var)

        new_body.extend(out_nodes)
        return index_vars

    def _add_offset_to_slice(self, slice_var, offset_var, out_nodes, scope,
                                loc):
        if isinstance(slice_var, slice):
            f_text = """def f(offset):
                return slice({} + offset, {} + offset)
            """.format(slice_var.start, slice_var.stop)
            loc = {}
            exec(f_text, {}, loc)
            f = loc['f']
            args = [offset_var]
            arg_typs = (types.intp,)
        else:
            def f(old_slice, offset):
                return slice(old_slice.start + offset, old_slice.stop + offset)
            args = [slice_var, offset_var]
            slice_type = self.typemap[slice_var.name]
            arg_typs = (slice_type, types.intp,)
        _globals = self.func_ir.func_id.func.__globals__
        f_ir = compile_to_numba_ir(f, _globals, self.typingctx, self.targetctx,
                                   arg_typs, self.typemap, self.calltypes)
        _, block = f_ir.blocks.popitem()
        replace_arg_nodes(block, args)
        new_index = block.body[-2].value.value
        out_nodes.extend(block.body[:-2])  # ignore return nodes
        return new_index

def get_stencil_ir(sf, typingctx, args, scope, loc, input_dict, typemap,
                                                                    calltypes):
    """get typed IR from stencil bytecode
    """
    from numba.core.cpu import CPUContext
    from numba.core.registry import cpu_target
    from numba.core.annotations import type_annotations
    from numba.core.typed_passes import type_inference_stage

    # get untyped IR
    stencil_func_ir = sf.kernel_ir.copy()
    # copy the IR nodes to avoid changing IR in the StencilFunc object
    stencil_blocks = copy.deepcopy(stencil_func_ir.blocks)
    stencil_func_ir.blocks = stencil_blocks

    name_var_table = ir_utils.get_name_var_table(stencil_func_ir.blocks)
    if "out" in name_var_table:
        raise ValueError("Cannot use the reserved word 'out' in stencil kernels.")

    # get typed IR with a dummy pipeline (similar to test_parfors.py)
    from numba.core.registry import cpu_target
    targetctx = cpu_target.target_context

    tp = DummyPipeline(typingctx, targetctx, args, stencil_func_ir)

    rewrites.rewrite_registry.apply('before-inference', tp.state)

    tp.state.typemap, tp.state.return_type, tp.state.calltypes, _ = type_inference_stage(
        tp.state.typingctx, tp.state.targetctx, tp.state.func_ir,
        tp.state.args, None)

    type_annotations.TypeAnnotation(
        func_ir=tp.state.func_ir,
        typemap=tp.state.typemap,
        calltypes=tp.state.calltypes,
        lifted=(),
        lifted_from=None,
        args=tp.state.args,
        return_type=tp.state.return_type,
        html_output=config.HTML)

    # make block labels unique
    stencil_blocks = ir_utils.add_offset_to_labels(stencil_blocks,
                                                        ir_utils.next_label())
    min_label = min(stencil_blocks.keys())
    max_label = max(stencil_blocks.keys())
    ir_utils._the_max_label.update(max_label)

    if config.DEBUG_ARRAY_OPT >= 1:
        print("Initial stencil_blocks")
        ir_utils.dump_blocks(stencil_blocks)

    # rename variables,
    var_dict = {}
    for v, typ in tp.state.typemap.items():
        new_var = ir.Var(scope, mk_unique_var(v), loc)
        var_dict[v] = new_var
        typemap[new_var.name] = typ  # add new var type for overall function
    ir_utils.replace_vars(stencil_blocks, var_dict)

    if config.DEBUG_ARRAY_OPT >= 1:
        print("After replace_vars")
        ir_utils.dump_blocks(stencil_blocks)

    # add call types to overall function
    for call, call_typ in tp.state.calltypes.items():
        calltypes[call] = call_typ

    arg_to_arr_dict = {}
    # replace arg with arr
    for block in stencil_blocks.values():
        for stmt in block.body:
            if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Arg):
                if config.DEBUG_ARRAY_OPT >= 1:
                    print("input_dict", input_dict, stmt.value.index,
                               stmt.value.name, stmt.value.index in input_dict)
                arg_to_arr_dict[stmt.value.name] = input_dict[stmt.value.index].name
                stmt.value = input_dict[stmt.value.index]

    if config.DEBUG_ARRAY_OPT >= 1:
        print("arg_to_arr_dict", arg_to_arr_dict)
        print("After replace arg with arr")
        ir_utils.dump_blocks(stencil_blocks)

    ir_utils.remove_dels(stencil_blocks)
    stencil_func_ir.blocks = stencil_blocks
    return stencil_func_ir, sf.get_return_type(args)[0], arg_to_arr_dict

class DummyPipeline(object):
    def __init__(self, typingctx, targetctx, args, f_ir):
        from numba.core.compiler import StateDict
        self.state = StateDict()
        self.state.typingctx = typingctx
        self.state.targetctx = targetctx
        self.state.args = args
        self.state.func_ir = f_ir
        self.state.typemap = None
        self.state.return_type = None
        self.state.calltypes = None


def _get_const_index_expr(stencil_ir, func_ir, index_var):
    """
    infer index_var as constant if it is of a expression form like c-1 where c
    is a constant in the outer function.
    index_var is assumed to be inside stencil kernel
    """
    const_val = guard(
        _get_const_index_expr_inner, stencil_ir, func_ir, index_var)
    if const_val is not None:
        return const_val
    return index_var

def _get_const_index_expr_inner(stencil_ir, func_ir, index_var):
    """inner constant inference function that calls constant, unary and binary
    cases.
    """
    require(isinstance(index_var, ir.Var))
    # case where the index is a const itself in outer function
    var_const =  guard(_get_const_two_irs, stencil_ir, func_ir, index_var)
    if var_const is not None:
        return var_const
    # get index definition
    index_def = ir_utils.get_definition(stencil_ir, index_var)
    # match inner_var = unary(index_var)
    var_const = guard(
        _get_const_unary_expr, stencil_ir, func_ir, index_def)
    if var_const is not None:
        return var_const
    # match inner_var = arg1 + arg2
    var_const = guard(
        _get_const_binary_expr, stencil_ir, func_ir, index_def)
    if var_const is not None:
        return var_const
    raise GuardException

def _get_const_two_irs(ir1, ir2, var):
    """get constant in either of two IRs if available
    otherwise, throw GuardException
    """
    var_const = guard(find_const, ir1, var)
    if var_const is not None:
        return var_const
    var_const = guard(find_const, ir2, var)
    if var_const is not None:
        return var_const
    raise GuardException

def _get_const_unary_expr(stencil_ir, func_ir, index_def):
    """evaluate constant unary expr if possible
    otherwise, raise GuardException
    """
    require(isinstance(index_def, ir.Expr) and index_def.op == 'unary')
    inner_var = index_def.value
    # return -c as constant
    const_val = _get_const_index_expr_inner(stencil_ir, func_ir, inner_var)
    op = OPERATORS_TO_BUILTINS[index_def.fn]
    return eval("{}{}".format(op, const_val))

def _get_const_binary_expr(stencil_ir, func_ir, index_def):
    """evaluate constant binary expr if possible
    otherwise, raise GuardException
    """
    require(isinstance(index_def, ir.Expr) and index_def.op == 'binop')
    arg1 = _get_const_index_expr_inner(stencil_ir, func_ir, index_def.lhs)
    arg2 = _get_const_index_expr_inner(stencil_ir, func_ir, index_def.rhs)
    op = OPERATORS_TO_BUILTINS[index_def.fn]
    return eval("{}{}{}".format(arg1, op, arg2))
