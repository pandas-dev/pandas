#
# Copyright (c) 2017 Intel Corporation
# SPDX-License-Identifier: BSD-2-Clause
#

import copy

import numpy as np
from llvmlite import ir as lir

from numba.core import types, typing, utils, ir, config, ir_utils, registry
from numba.core.typing.templates import (CallableTemplate, signature,
                                         infer_global, AbstractTemplate)
from numba.core.imputils import lower_builtin
from numba.core.extending import register_jitable
from numba.core.errors import NumbaValueError
from numba.misc.special import literal_unroll
import numba

import operator
from numba.np import numpy_support

class StencilFuncLowerer(object):
    '''Callable class responsible for lowering calls to a specific StencilFunc.
    '''
    def __init__(self, sf):
        self.stencilFunc = sf

    def __call__(self, context, builder, sig, args):
        cres = self.stencilFunc.compile_for_argtys(sig.args, {},
                    sig.return_type, None)
        res = context.call_internal(builder, cres.fndesc, sig, args)
        context.add_linking_libs([cres.library])
        return res

@register_jitable
def raise_if_incompatible_array_sizes(a, *args):
    ashape = a.shape

    # We need literal_unroll here because the stencil might take
    # multiple input arrays with different types that are not compatible
    # (e.g. values as float[:] and flags as bool[:])
    # When more than three total arrays are given, the second and third
    # are iterated over in the loop below. Without literal_unroll, their
    # types have to match.
    # An example failing signature without literal_unroll might be
    # (float[:], float[:], bool[:]) (Just (float[:], bool[:]) wouldn't fail)
    for arg in literal_unroll(args):
        if a.ndim != arg.ndim:
            raise ValueError("Secondary stencil array does not have same number "
                             " of dimensions as the first stencil input.")
        argshape = arg.shape
        for i in range(len(ashape)):
            if ashape[i] > argshape[i]:
                raise ValueError("Secondary stencil array has some dimension "
                                 "smaller the same dimension in the first "
                                 "stencil input.")

def slice_addition(the_slice, addend):
    """ Called by stencil in Python mode to add the loop index to a
        user-specified slice.
    """
    return slice(the_slice.start + addend, the_slice.stop + addend)

class StencilFunc(object):
    """
    A special type to hold stencil information for the IR.
    """

    id_counter = 0

    def __init__(self, kernel_ir, mode, options):
        self.id = type(self).id_counter
        type(self).id_counter += 1
        self.kernel_ir = kernel_ir
        self.mode = mode
        self.options = options
        self.kws = []       # remember original kws arguments

        # stencils only supported for CPU context currently
        self._typingctx = registry.cpu_target.typing_context
        self._targetctx = registry.cpu_target.target_context
        self._install_type(self._typingctx)
        self.neighborhood = self.options.get("neighborhood")
        self._type_cache = {}
        self._lower_me = StencilFuncLowerer(self)

    def replace_return_with_setitem(self, blocks, index_vars, out_name):
        """
        Find return statements in the IR and replace them with a SetItem
        call of the value "returned" by the kernel into the result array.
        Returns the block labels that contained return statements.
        """
        ret_blocks = []

        for label, block in blocks.items():
            scope = block.scope
            loc = block.loc
            new_body = []
            for stmt in block.body:
                if isinstance(stmt, ir.Return):
                    ret_blocks.append(label)
                    # If 1D array then avoid the tuple construction.
                    if len(index_vars) == 1:
                        rvar = ir.Var(scope, out_name, loc)
                        ivar = ir.Var(scope, index_vars[0], loc)
                        new_body.append(ir.SetItem(rvar, ivar, stmt.value, loc))
                    else:
                        # Convert the string names of the index variables into
                        # ir.Var's.
                        var_index_vars = []
                        for one_var in index_vars:
                            index_var = ir.Var(scope, one_var, loc)
                            var_index_vars += [index_var]

                        s_index_var = scope.redefine("stencil_index", loc)
                        # Build a tuple from the index ir.Var's.
                        tuple_call = ir.Expr.build_tuple(var_index_vars, loc)
                        new_body.append(ir.Assign(tuple_call, s_index_var, loc))
                        rvar = ir.Var(scope, out_name, loc)
                        # Write the return statements original value into
                        # the array using the tuple index.
                        si = ir.SetItem(rvar, s_index_var, stmt.value, loc)
                        new_body.append(si)
                else:
                    new_body.append(stmt)
            block.body = new_body
        return ret_blocks

    def add_indices_to_kernel(self, kernel, index_names, ndim,
                              neighborhood, standard_indexed, typemap, calltypes):
        """
        Transforms the stencil kernel as specified by the user into one
        that includes each dimension's index variable as part of the getitem
        calls.  So, in effect array[-1] becomes array[index0-1].
        """
        const_dict = {}
        kernel_consts = []

        if config.DEBUG_ARRAY_OPT >= 1:
            print("add_indices_to_kernel", ndim, neighborhood)
            ir_utils.dump_blocks(kernel.blocks)

        if neighborhood is None:
            need_to_calc_kernel = True
        else:
            need_to_calc_kernel = False
            if len(neighborhood) != ndim:
                raise NumbaValueError("%d dimensional neighborhood specified "
                                      "for %d dimensional input array" %
                                      (len(neighborhood), ndim))

        tuple_table = ir_utils.get_tuple_table(kernel.blocks)

        relatively_indexed = set()

        for block in kernel.blocks.values():
            scope = block.scope
            loc = block.loc
            new_body = []
            for stmt in block.body:
                if (isinstance(stmt, ir.Assign) and
                    isinstance(stmt.value, ir.Const)):
                    if config.DEBUG_ARRAY_OPT >= 1:
                        print("remembering in const_dict", stmt.target.name,
                              stmt.value.value)
                    # Remember consts for use later.
                    const_dict[stmt.target.name] = stmt.value.value
                if ((isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Expr)
                        and stmt.value.op in ['setitem', 'static_setitem']
                        and stmt.value.value.name in kernel.arg_names) or
                   (isinstance(stmt, ir.SetItem)
                        and stmt.target.name in kernel.arg_names)):
                    raise NumbaValueError("Assignments to arrays passed to " \
                                          "stencil kernels is not allowed.")
                if (isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Expr)
                        and stmt.value.op in ['getitem', 'static_getitem']
                        and stmt.value.value.name in kernel.arg_names
                        and stmt.value.value.name not in standard_indexed):
                    # We found a getitem from the input array.
                    if stmt.value.op == 'getitem':
                        stmt_index_var = stmt.value.index
                    else:
                        stmt_index_var = stmt.value.index_var
                        # allow static_getitem since rewrite passes are applied
                        #raise ValueError("Unexpected static_getitem in add_indices_to_kernel.")

                    relatively_indexed.add(stmt.value.value.name)

                    # Store the index used after looking up the variable in
                    # the const dictionary.
                    if need_to_calc_kernel:
                        assert hasattr(stmt_index_var, 'name')

                        if stmt_index_var.name in tuple_table:
                            kernel_consts += [tuple_table[stmt_index_var.name]]
                        elif stmt_index_var.name in const_dict:
                            kernel_consts += [const_dict[stmt_index_var.name]]
                        else:
                            raise NumbaValueError("stencil kernel index is not "
                                "constant, 'neighborhood' option required")

                    if ndim == 1:
                        # Single dimension always has index variable 'index0'.
                        # tmpvar will hold the real index and is computed by
                        # adding the relative offset in stmt.value.index to
                        # the current absolute location in index0.
                        index_var = ir.Var(scope, index_names[0], loc)
                        tmpvar = scope.redefine("stencil_index", loc)
                        stmt_index_var_typ = typemap[stmt_index_var.name]
                        # If the array is indexed with a slice then we
                        # have to add the index value with a call to
                        # slice_addition.
                        if isinstance(stmt_index_var_typ, types.misc.SliceType):
                            sa_var = scope.redefine("slice_addition", loc)
                            sa_func = numba.njit(slice_addition)
                            sa_func_typ = types.functions.Dispatcher(sa_func)
                            typemap[sa_var.name] = sa_func_typ
                            g_sa = ir.Global("slice_addition", sa_func, loc)
                            new_body.append(ir.Assign(g_sa, sa_var, loc))
                            slice_addition_call = ir.Expr.call(sa_var, [stmt_index_var, index_var], (), loc)
                            calltypes[slice_addition_call] = sa_func_typ.get_call_type(self._typingctx, [stmt_index_var_typ, types.intp], {})
                            new_body.append(ir.Assign(slice_addition_call, tmpvar, loc))
                            new_body.append(ir.Assign(
                                           ir.Expr.getitem(stmt.value.value, tmpvar, loc),
                                           stmt.target, loc))
                        else:
                            acc_call = ir.Expr.binop(operator.add, stmt_index_var,
                                                     index_var, loc)
                            new_body.append(ir.Assign(acc_call, tmpvar, loc))
                            new_body.append(ir.Assign(
                                           ir.Expr.getitem(stmt.value.value, tmpvar, loc),
                                           stmt.target, loc))
                    else:
                        index_vars = []
                        sum_results = []
                        s_index_var = scope.redefine("stencil_index", loc)
                        const_index_vars = []
                        ind_stencils = []

                        stmt_index_var_typ = typemap[stmt_index_var.name]
                        # Same idea as above but you have to extract
                        # individual elements out of the tuple indexing
                        # expression and add the corresponding index variable
                        # to them and then reconstitute as a tuple that can
                        # index the array.
                        for dim in range(ndim):
                            tmpvar = scope.redefine("const_index", loc)
                            new_body.append(ir.Assign(ir.Const(dim, loc),
                                                      tmpvar, loc))
                            const_index_vars += [tmpvar]
                            index_var = ir.Var(scope, index_names[dim], loc)
                            index_vars += [index_var]

                            tmpvar = scope.redefine("ind_stencil_index", loc)
                            ind_stencils += [tmpvar]
                            getitemvar = scope.redefine("getitem", loc)
                            getitemcall = ir.Expr.getitem(stmt_index_var,
                                                       const_index_vars[dim], loc)
                            new_body.append(ir.Assign(getitemcall, getitemvar, loc))
                            # Get the type of this particular part of the index tuple.
                            if isinstance(stmt_index_var_typ, types.ConstSized):
                                one_index_typ = stmt_index_var_typ[dim]
                            else:
                                one_index_typ = stmt_index_var_typ[:]
                            # If the array is indexed with a slice then we
                            # have to add the index value with a call to
                            # slice_addition.
                            if isinstance(one_index_typ, types.misc.SliceType):
                                sa_var = scope.redefine("slice_addition", loc)
                                sa_func = numba.njit(slice_addition)
                                sa_func_typ = types.functions.Dispatcher(sa_func)
                                typemap[sa_var.name] = sa_func_typ
                                g_sa = ir.Global("slice_addition", sa_func, loc)
                                new_body.append(ir.Assign(g_sa, sa_var, loc))
                                slice_addition_call = ir.Expr.call(sa_var, [getitemvar, index_vars[dim]], (), loc)
                                calltypes[slice_addition_call] = sa_func_typ.get_call_type(self._typingctx, [one_index_typ, types.intp], {})
                                new_body.append(ir.Assign(slice_addition_call, tmpvar, loc))
                            else:
                                acc_call = ir.Expr.binop(operator.add, getitemvar,
                                                         index_vars[dim], loc)
                                new_body.append(ir.Assign(acc_call, tmpvar, loc))

                        tuple_call = ir.Expr.build_tuple(ind_stencils, loc)
                        new_body.append(ir.Assign(tuple_call, s_index_var, loc))
                        new_body.append(ir.Assign(
                                  ir.Expr.getitem(stmt.value.value,s_index_var,loc),
                                  stmt.target,loc))
                else:
                    new_body.append(stmt)
            block.body = new_body

        if need_to_calc_kernel:
            # Find the size of the kernel by finding the maximum absolute value
            # index used in the kernel specification.
            neighborhood = [[0,0] for _ in range(ndim)]
            if len(kernel_consts) == 0:
                raise NumbaValueError("Stencil kernel with no accesses to "
                                      "relatively indexed arrays.")

            for index in kernel_consts:
                if isinstance(index, tuple) or isinstance(index, list):
                    for i in range(len(index)):
                        te = index[i]
                        if isinstance(te, ir.Var) and te.name in const_dict:
                            te = const_dict[te.name]
                        if isinstance(te, int):
                            neighborhood[i][0] = min(neighborhood[i][0], te)
                            neighborhood[i][1] = max(neighborhood[i][1], te)
                        else:
                            raise NumbaValueError(
                                "stencil kernel index is not constant,"
                                "'neighborhood' option required")
                    index_len = len(index)
                elif isinstance(index, int):
                    neighborhood[0][0] = min(neighborhood[0][0], index)
                    neighborhood[0][1] = max(neighborhood[0][1], index)
                    index_len = 1
                else:
                    raise NumbaValueError(
                        "Non-tuple or non-integer used as stencil index.")
                if index_len != ndim:
                    raise NumbaValueError(
                        "Stencil index does not match array dimensionality.")

        return (neighborhood, relatively_indexed)


    def get_return_type(self, argtys):
        if config.DEBUG_ARRAY_OPT >= 1:
            print("get_return_type", argtys)
            ir_utils.dump_blocks(self.kernel_ir.blocks)

        if not isinstance(argtys[0], types.npytypes.Array):
            raise NumbaValueError("The first argument to a stencil kernel must "
                                  "be the primary input array.")

        from numba.core import typed_passes
        typemap, return_type, calltypes, _ = typed_passes.type_inference_stage(
                self._typingctx,
                self._targetctx,
                self.kernel_ir,
                argtys,
                None,
                {})
        if isinstance(return_type, types.npytypes.Array):
            raise NumbaValueError(
                "Stencil kernel must return a scalar and not a numpy array.")

        real_ret = types.npytypes.Array(return_type, argtys[0].ndim,
                                                     argtys[0].layout)
        return (real_ret, typemap, calltypes)

    def _install_type(self, typingctx):
        """Constructs and installs a typing class for a StencilFunc object in
        the input typing context.
        """
        _ty_cls = type('StencilFuncTyping_' +
                       str(self.id),
                       (AbstractTemplate,),
                       dict(key=self, generic=self._type_me))
        typingctx.insert_user_function(self, _ty_cls)

    def compile_for_argtys(self, argtys, kwtys, return_type, sigret):
        # look in the type cache to find if result array is passed
        (_, result, typemap, calltypes) = self._type_cache[argtys]
        new_func = self._stencil_wrapper(result, sigret, return_type,
                                         typemap, calltypes, *argtys)
        return new_func

    def _type_me(self, argtys, kwtys):
        """
        Implement AbstractTemplate.generic() for the typing class
        built by StencilFunc._install_type().
        Return the call-site signature.
        """
        if (self.neighborhood is not None and
            len(self.neighborhood) != argtys[0].ndim):
            raise NumbaValueError("%d dimensional neighborhood specified "
                                  "for %d dimensional input array" %
                                  (len(self.neighborhood), argtys[0].ndim))

        argtys_extra = argtys
        sig_extra = ""
        result = None
        if 'out' in kwtys:
            argtys_extra += (kwtys['out'],)
            sig_extra += ", out=None"
            result = kwtys['out']

        if 'neighborhood' in kwtys:
            argtys_extra += (kwtys['neighborhood'],)
            sig_extra += ", neighborhood=None"

        # look in the type cache first
        if argtys_extra in self._type_cache:
            (_sig, _, _, _) = self._type_cache[argtys_extra]
            return _sig

        (real_ret, typemap, calltypes) = self.get_return_type(argtys)
        sig = signature(real_ret, *argtys_extra)
        dummy_text = ("def __numba_dummy_stencil({}{}):\n    pass\n".format(
                        ",".join(self.kernel_ir.arg_names), sig_extra))
        dct = {}
        exec(dummy_text, dct)
        dummy_func = dct["__numba_dummy_stencil"]
        sig = sig.replace(pysig=utils.pysignature(dummy_func))
        self._targetctx.insert_func_defn([(self._lower_me, self, argtys_extra)])
        self._type_cache[argtys_extra] = (sig, result, typemap, calltypes)
        return sig

    def copy_ir_with_calltypes(self, ir, calltypes):
        """
        Create a copy of a given IR along with its calltype information.
        We need a copy of the calltypes because copy propagation applied
        to the copied IR will change the calltypes and make subsequent
        uses of the original IR invalid.
        """
        copy_calltypes = {}
        kernel_copy = ir.copy()
        kernel_copy.blocks = {}
        # For each block...
        for (block_label, block) in ir.blocks.items():
            new_block = copy.deepcopy(ir.blocks[block_label])
            new_block.body = []
            # For each statement in each block...
            for stmt in ir.blocks[block_label].body:
                # Copy the statement to the new copy of the kernel
                # and if the original statement is in the original
                # calltypes then add the type associated with this
                # statement to the calltypes copy.
                scopy = copy.deepcopy(stmt)
                new_block.body.append(scopy)
                if stmt in calltypes:
                    copy_calltypes[scopy] = calltypes[stmt]
            kernel_copy.blocks[block_label] = new_block
        return (kernel_copy, copy_calltypes)

    def _stencil_wrapper(self, result, sigret, return_type, typemap, calltypes, *args):
        # Overall approach:
        # 1) Construct a string containing a function definition for the stencil function
        #    that will execute the stencil kernel.  This function definition includes a
        #    unique stencil function name, the parameters to the stencil kernel, loop
        #    nests across the dimensions of the input array.  Those loop nests use the
        #    computed stencil kernel size so as not to try to compute elements where
        #    elements outside the bounds of the input array would be needed.
        # 2) The but of the loop nest in this new function is a special sentinel
        #    assignment.
        # 3) Get the IR of this new function.
        # 4) Split the block containing the sentinel assignment and remove the sentinel
        #    assignment.  Insert the stencil kernel IR into the stencil function IR
        #    after label and variable renaming of the stencil kernel IR to prevent
        #    conflicts with the stencil function IR.
        # 5) Compile the combined stencil function IR + stencil kernel IR into existence.

        # Copy the kernel so that our changes for this callsite
        # won't effect other callsites.
        (kernel_copy, copy_calltypes) = self.copy_ir_with_calltypes(
                                            self.kernel_ir, calltypes)
        # The stencil kernel body becomes the body of a loop, for which args aren't needed.
        ir_utils.remove_args(kernel_copy.blocks)
        first_arg = kernel_copy.arg_names[0]

        in_cps, out_cps = ir_utils.copy_propagate(kernel_copy.blocks, typemap)
        name_var_table = ir_utils.get_name_var_table(kernel_copy.blocks)
        ir_utils.apply_copy_propagate(
            kernel_copy.blocks,
            in_cps,
            name_var_table,
            typemap,
            copy_calltypes)

        if "out" in name_var_table:
            raise NumbaValueError("Cannot use the reserved word 'out' in stencil kernels.")

        sentinel_name = ir_utils.get_unused_var_name("__sentinel__", name_var_table)
        if config.DEBUG_ARRAY_OPT >= 1:
            print("name_var_table", name_var_table, sentinel_name)

        the_array = args[0]

        if config.DEBUG_ARRAY_OPT >= 1:
            print("_stencil_wrapper", return_type, return_type.dtype,
                                      type(return_type.dtype), args)
            ir_utils.dump_blocks(kernel_copy.blocks)

        # We generate a Numba function to execute this stencil and here
        # create the unique name of this function.
        stencil_func_name = "__numba_stencil_%s_%s" % (
                                        hex(id(the_array)).replace("-", "_"),
                                        self.id)

        # We will put a loop nest in the generated function for each
        # dimension in the input array.  Here we create the name for
        # the index variable for each dimension.  index0, index1, ...
        index_vars = []
        for i in range(the_array.ndim):
            index_var_name = ir_utils.get_unused_var_name("index" + str(i),
                                                          name_var_table)
            index_vars += [index_var_name]

        # Create extra signature for out and neighborhood.
        out_name = ir_utils.get_unused_var_name("out", name_var_table)
        neighborhood_name = ir_utils.get_unused_var_name("neighborhood",
                                                         name_var_table)
        sig_extra = ""
        if result is not None:
            sig_extra += ", {}=None".format(out_name)
        if "neighborhood" in dict(self.kws):
            sig_extra += ", {}=None".format(neighborhood_name)

        # Get a list of the standard indexed array names.
        standard_indexed = self.options.get("standard_indexing", [])

        if first_arg in standard_indexed:
            raise NumbaValueError("The first argument to a stencil kernel must "
                                  "use relative indexing, not standard indexing.")

        if len(set(standard_indexed) - set(kernel_copy.arg_names)) != 0:
            raise NumbaValueError("Standard indexing requested for an array name "
                                  "not present in the stencil kernel definition.")

        # Add index variables to getitems in the IR to transition the accesses
        # in the kernel from relative to regular Python indexing.  Returns the
        # computed size of the stencil kernel and a list of the relatively indexed
        # arrays.
        kernel_size, relatively_indexed = self.add_indices_to_kernel(
                kernel_copy, index_vars, the_array.ndim,
                self.neighborhood, standard_indexed, typemap, copy_calltypes)
        if self.neighborhood is None:
            self.neighborhood = kernel_size

        if config.DEBUG_ARRAY_OPT >= 1:
            print("After add_indices_to_kernel")
            ir_utils.dump_blocks(kernel_copy.blocks)

        # The return in the stencil kernel becomes a setitem for that
        # particular point in the iteration space.
        ret_blocks = self.replace_return_with_setitem(kernel_copy.blocks,
                                                      index_vars, out_name)

        if config.DEBUG_ARRAY_OPT >= 1:
            print("After replace_return_with_setitem", ret_blocks)
            ir_utils.dump_blocks(kernel_copy.blocks)

        # Start to form the new function to execute the stencil kernel.
        func_text = "def {}({}{}):\n".format(stencil_func_name,
                        ",".join(kernel_copy.arg_names), sig_extra)

        # Get loop ranges for each dimension, which could be either int
        # or variable. In the latter case we'll use the extra neighborhood
        # argument to the function.
        ranges = []
        for i in range(the_array.ndim):
            if isinstance(kernel_size[i][0], int):
                lo = kernel_size[i][0]
                hi = kernel_size[i][1]
            else:
                lo = "{}[{}][0]".format(neighborhood_name, i)
                hi = "{}[{}][1]".format(neighborhood_name, i)
            ranges.append((lo, hi))

        # If there are more than one relatively indexed arrays, add a call to
        # a function that will raise an error if any of the relatively indexed
        # arrays are of different size than the first input array.
        if len(relatively_indexed) > 1:
            func_text += "    raise_if_incompatible_array_sizes(" + first_arg
            for other_array in relatively_indexed:
                if other_array != first_arg:
                    func_text += "," + other_array
            func_text += ")\n"

        # Get the shape of the first input array.
        shape_name = ir_utils.get_unused_var_name("full_shape", name_var_table)
        func_text += "    {} = {}.shape\n".format(shape_name, first_arg)

        # Converts cval to a string constant
        def cval_as_str(cval):
            if not np.isfinite(cval):
                # See if this is a string-repr numerical const, issue #7286
                if np.isnan(cval):
                    return "np.nan"
                elif np.isinf(cval):
                    if cval < 0:
                        return "-np.inf"
                    else:
                        return "np.inf"
            else:
                return str(cval)

        # If we have to allocate the output array (the out argument was not used)
        # then us numpy.full if the user specified a cval stencil decorator option
        # or np.zeros if they didn't to allocate the array.
        if result is None:
            return_type_name = numpy_support.as_dtype(
                               return_type.dtype).type.__name__
            out_init ="{} = np.empty({}, dtype=np.{})\n".format(
                        out_name, shape_name, return_type_name)

            if "cval" in self.options:
                cval = self.options["cval"]
                cval_ty = typing.typeof.typeof(cval)
                if not self._typingctx.can_convert(cval_ty, return_type.dtype):
                    msg = "cval type does not match stencil return type."
                    raise NumbaValueError(msg)
            else:
                 cval = 0
            func_text += "    " + out_init
            for dim in range(the_array.ndim):
                start_items = [":"] * the_array.ndim
                end_items = [":"] * the_array.ndim
                start_items[dim] = ":-{}".format(self.neighborhood[dim][0])
                end_items[dim] = "-{}:".format(self.neighborhood[dim][1])
                func_text += "    " + "{}[{}] = {}\n".format(out_name, ",".join(start_items), cval_as_str(cval))
                func_text += "    " + "{}[{}] = {}\n".format(out_name, ",".join(end_items), cval_as_str(cval))
        else: # result is present, if cval is set then use it
            if "cval" in self.options:
                cval = self.options["cval"]
                cval_ty = typing.typeof.typeof(cval)
                if not self._typingctx.can_convert(cval_ty, return_type.dtype):
                    msg = "cval type does not match stencil return type."
                    raise NumbaValueError(msg)
                out_init = "{}[:] = {}\n".format(out_name, cval_as_str(cval))
                func_text += "    " + out_init

        offset = 1
        # Add the loop nests to the new function.
        for i in range(the_array.ndim):
            for j in range(offset):
                func_text += "    "
            # ranges[i][0] is the minimum index used in the i'th dimension
            # but minimum's greater than 0 don't preclude any entry in the array.
            # So, take the minimum of 0 and the minimum index found in the kernel
            # and this will be a negative number (potentially -0).  Then, we do
            # unary - on that to get the positive offset in this dimension whose
            # use is precluded.
            # ranges[i][1] is the maximum of 0 and the observed maximum index
            # in this dimension because negative maximums would not cause us to
            # preclude any entry in the array from being used.
            func_text += ("for {} in range(-min(0,{}),"
                          "{}[{}]-max(0,{})):\n").format(
                            index_vars[i],
                            ranges[i][0],
                            shape_name,
                            i,
                            ranges[i][1])
            offset += 1

        for j in range(offset):
            func_text += "    "
        # Put a sentinel in the code so we can locate it in the IR.  We will
        # remove this sentinel assignment and replace it with the IR for the
        # stencil kernel body.
        func_text += "{} = 0\n".format(sentinel_name)
        func_text += "    return {}\n".format(out_name)

        if config.DEBUG_ARRAY_OPT >= 1:
            print("new stencil func text")
            print(func_text)

        # Force the new stencil function into existence.
        dct = {}
        dct.update(globals())
        exec(func_text, dct)
        stencil_func = dct[stencil_func_name]
        if sigret is not None:
            pysig = utils.pysignature(stencil_func)
            sigret.pysig = pysig
        # Get the IR for the newly created stencil function.
        from numba.core import compiler
        stencil_ir = compiler.run_frontend(stencil_func)
        ir_utils.remove_dels(stencil_ir.blocks)

        # rename all variables in stencil_ir afresh
        var_table = ir_utils.get_name_var_table(stencil_ir.blocks)
        new_var_dict = {}
        reserved_names = ([sentinel_name, out_name, neighborhood_name,
                           shape_name] + kernel_copy.arg_names + index_vars)
        for name, var in var_table.items():
            if not name in reserved_names:
                assert isinstance(var, ir.Var)
                new_var = var.scope.redefine(var.name, var.loc)
                new_var_dict[name] = new_var.name
        ir_utils.replace_var_names(stencil_ir.blocks, new_var_dict)

        stencil_stub_last_label = max(stencil_ir.blocks.keys()) + 1

        # Shift labels in the kernel copy so they are guaranteed unique
        # and don't conflict with any labels in the stencil_ir.
        kernel_copy.blocks = ir_utils.add_offset_to_labels(
                                kernel_copy.blocks, stencil_stub_last_label)
        new_label = max(kernel_copy.blocks.keys()) + 1
        # Adjust ret_blocks to account for addition of the offset.
        ret_blocks = [x + stencil_stub_last_label for x in ret_blocks]

        if config.DEBUG_ARRAY_OPT >= 1:
            print("ret_blocks w/ offsets", ret_blocks, stencil_stub_last_label)
            print("before replace sentinel stencil_ir")
            ir_utils.dump_blocks(stencil_ir.blocks)
            print("before replace sentinel kernel_copy")
            ir_utils.dump_blocks(kernel_copy.blocks)

        # Search all the block in the stencil outline for the sentinel.
        for label, block in stencil_ir.blocks.items():
            for i, inst in enumerate(block.body):
                if (isinstance( inst, ir.Assign) and
                    inst.target.name == sentinel_name):
                    # We found the sentinel assignment.
                    loc = inst.loc
                    scope = block.scope
                    # split block across __sentinel__
                    # A new block is allocated for the statements prior to the
                    # sentinel but the new block maintains the current block
                    # label.
                    prev_block = ir.Block(scope, loc)
                    prev_block.body = block.body[:i]
                    # The current block is used for statements after sentinel.
                    block.body = block.body[i + 1:]
                    # But the current block gets a new label.
                    body_first_label = min(kernel_copy.blocks.keys())

                    # The previous block jumps to the minimum labelled block of
                    # the parfor body.
                    prev_block.append(ir.Jump(body_first_label, loc))
                    # Add all the parfor loop body blocks to the gufunc
                    # function's IR.
                    for (l, b) in kernel_copy.blocks.items():
                        stencil_ir.blocks[l] = b

                    stencil_ir.blocks[new_label] = block
                    stencil_ir.blocks[label] = prev_block
                    # Add a jump from all the blocks that previously contained
                    # a return in the stencil kernel to the block
                    # containing statements after the sentinel.
                    for ret_block in ret_blocks:
                        stencil_ir.blocks[ret_block].append(
                            ir.Jump(new_label, loc))
                    break
            else:
                continue
            break

        stencil_ir.blocks = ir_utils.rename_labels(stencil_ir.blocks)
        ir_utils.remove_dels(stencil_ir.blocks)

        assert(isinstance(the_array, types.Type))
        array_types = args

        new_stencil_param_types = list(array_types)

        if config.DEBUG_ARRAY_OPT >= 1:
            print("new_stencil_param_types", new_stencil_param_types)
            ir_utils.dump_blocks(stencil_ir.blocks)

        # Compile the combined stencil function with the replaced loop
        # body in it.
        ir_utils.fixup_var_define_in_scope(stencil_ir.blocks)
        new_func = compiler.compile_ir(
            self._typingctx,
            self._targetctx,
            stencil_ir,
            new_stencil_param_types,
            None,
            compiler.DEFAULT_FLAGS,
            {})
        return new_func

    def __call__(self, *args, **kwargs):
        self._typingctx.refresh()
        if (self.neighborhood is not None and
            len(self.neighborhood) != args[0].ndim):
            raise NumbaValueError("{} dimensional neighborhood specified for "
                                  "{} dimensional input array".format(
                                  len(self.neighborhood), args[0].ndim))

        if 'out' in kwargs:
            result = kwargs['out']
            rdtype = result.dtype
            rttype = numpy_support.from_dtype(rdtype)
            result_type = types.npytypes.Array(rttype, result.ndim,
                                               numpy_support.map_layout(result))
            array_types = tuple([typing.typeof.typeof(x) for x in args])
            array_types_full = tuple([typing.typeof.typeof(x) for x in args] +
                                     [result_type])
        else:
            result = None
            array_types = tuple([typing.typeof.typeof(x) for x in args])
            array_types_full = array_types

        if config.DEBUG_ARRAY_OPT >= 1:
            print("__call__", array_types, args, kwargs)

        (real_ret, typemap, calltypes) = self.get_return_type(array_types)
        new_func = self._stencil_wrapper(result, None, real_ret, typemap,
                                         calltypes, *array_types_full)

        if result is None:
            return new_func.entry_point(*args)
        else:
            return new_func.entry_point(*(args+(result,)))

def stencil(func_or_mode='constant', **options):
    # called on function without specifying mode style
    if not isinstance(func_or_mode, str):
        mode = 'constant'  # default style
        func = func_or_mode
    else:
        mode = func_or_mode
        func = None

    for option in options:
        if option not in ["cval", "standard_indexing", "neighborhood"]:
            raise NumbaValueError("Unknown stencil option " + option)

    wrapper = _stencil(mode, options)
    if func is not None:
        return wrapper(func)
    return wrapper

def _stencil(mode, options):
    if mode != 'constant':
        raise NumbaValueError("Unsupported mode style " + mode)

    def decorated(func):
        from numba.core import compiler
        kernel_ir = compiler.run_frontend(func)
        return StencilFunc(kernel_ir, mode, options)

    return decorated

@lower_builtin(stencil)
def stencil_dummy_lower(context, builder, sig, args):
    "lowering for dummy stencil calls"
    return lir.Constant(lir.IntType(types.intp.bitwidth), 0)
