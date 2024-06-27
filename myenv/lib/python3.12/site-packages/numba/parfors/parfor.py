#
# Copyright (c) 2017 Intel Corporation
# SPDX-License-Identifier: BSD-2-Clause
#

"""
This module transforms data-parallel operations such as Numpy calls into
'Parfor' nodes, which are nested loops that can be parallelized.
It also implements optimizations such as loop fusion, and extends the rest of
compiler analysis and optimizations to support Parfors.
This is similar to ParallelAccelerator package in Julia:
https://github.com/IntelLabs/ParallelAccelerator.jl
'Parallelizing Julia with a Non-invasive DSL', T. Anderson et al., ECOOP'17.
"""
import types as pytypes  # avoid confusion with numba.types
import sys, math
import os
import textwrap
import copy
import inspect
import linecache
from functools import reduce
from collections import defaultdict, OrderedDict, namedtuple
from contextlib import contextmanager
import operator
from dataclasses import make_dataclass
import warnings

from llvmlite import ir as lir
from numba.core.imputils import impl_ret_untracked
import numba.core.ir
from numba.core import types, typing, utils, errors, ir, analysis, postproc, rewrites, typeinfer, config, ir_utils
from numba import prange, pndindex
from numba.np.npdatetime_helpers import datetime_minimum, datetime_maximum
from numba.np.numpy_support import as_dtype, numpy_version
from numba.core.typing.templates import infer_global, AbstractTemplate
from numba.stencils.stencilparfor import StencilPass
from numba.core.extending import register_jitable, lower_builtin

from numba.core.ir_utils import (
    mk_unique_var,
    next_label,
    mk_alloc,
    get_np_ufunc_typ,
    mk_range_block,
    mk_loop_header,
    get_name_var_table,
    replace_vars,
    replace_vars_inner,
    visit_vars,
    visit_vars_inner,
    remove_dead,
    copy_propagate,
    get_block_copies,
    apply_copy_propagate,
    dprint_func_ir,
    find_topo_order,
    get_stmt_writes,
    rename_labels,
    get_call_table,
    simplify,
    simplify_CFG,
    has_no_side_effect,
    canonicalize_array_math,
    add_offset_to_labels,
    find_callname,
    find_build_sequence,
    guard,
    require,
    GuardException,
    compile_to_numba_ir,
    get_definition,
    build_definitions,
    replace_arg_nodes,
    replace_returns,
    is_getitem,
    is_setitem,
    is_get_setitem,
    index_var_of_get_setitem,
    set_index_var_of_get_setitem,
    find_potential_aliases,
    replace_var_names,
    transfer_scope,
)

from numba.core.analysis import (compute_use_defs, compute_live_map,
                            compute_dead_maps, compute_cfg_from_blocks)
from numba.core.controlflow import CFGraph
from numba.core.typing import npydecl, signature
from numba.core.types.functions import Function
from numba.parfors.array_analysis import (random_int_args, random_1arg_size,
                                  random_2arg_sizelast, random_3arg_sizelast,
                                  random_calls, assert_equiv)
from numba.core.extending import overload
import copy
import numpy
import numpy as np
from numba.parfors import array_analysis
import numba.cpython.builtins
from numba.stencils import stencilparfor
# circular dependency: import numba.npyufunc.dufunc.DUFunc

# wrapped pretty print
_termwidth = 80
_txtwrapper = textwrap.TextWrapper(width=_termwidth, drop_whitespace=False)
def print_wrapped(x):
    for l in x.splitlines():
        [print(y) for y in _txtwrapper.wrap(l)]

sequential_parfor_lowering = False

# init_prange is a sentinel call that specifies the start of the initialization
# code for the computation in the upcoming prange call
# This lets the prange pass to put the code in the generated parfor's init_block
def init_prange():
    return

@overload(init_prange)
def init_prange_overload():
    def no_op():
        return
    return no_op

class internal_prange(object):

    def __new__(cls, *args):
        return range(*args)


def min_parallel_impl(return_type, arg):
    # XXX: use prange for 1D arrays since pndindex returns a 1-tuple instead of
    # integer. This causes type and fusion issues.
    if arg.ndim == 0:
        def min_1(in_arr):
            return in_arr[()]
    elif arg.ndim == 1:
        if isinstance(arg.dtype, (types.NPDatetime, types.NPTimedelta)):
            # NaT is always returned if it is in the array
            def min_1(in_arr):
                numba.parfors.parfor.init_prange()
                min_checker(len(in_arr))
                val = numba.cpython.builtins.get_type_max_value(in_arr.dtype)
                for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                    val = datetime_minimum(val, in_arr[i])
                return val
        else:
            def min_1(in_arr):
                numba.parfors.parfor.init_prange()
                min_checker(len(in_arr))
                val = numba.cpython.builtins.get_type_max_value(in_arr.dtype)
                for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                    val = min(val, in_arr[i])
                return val
    else:
        def min_1(in_arr):
            numba.parfors.parfor.init_prange()
            min_checker(len(in_arr))
            val = numba.cpython.builtins.get_type_max_value(in_arr.dtype)
            for i in numba.pndindex(in_arr.shape):
                val = min(val, in_arr[i])
            return val
    return min_1

def max_parallel_impl(return_type, arg):
    if arg.ndim == 0:
        def max_1(in_arr):
            return in_arr[()]
    elif arg.ndim == 1:
        if isinstance(arg.dtype, (types.NPDatetime, types.NPTimedelta)):
            # NaT is always returned if it is in the array
            def max_1(in_arr):
                numba.parfors.parfor.init_prange()
                max_checker(len(in_arr))
                val = numba.cpython.builtins.get_type_min_value(in_arr.dtype)
                for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                    val = datetime_maximum(val, in_arr[i])
                return val
        else:
            def max_1(in_arr):
                numba.parfors.parfor.init_prange()
                max_checker(len(in_arr))
                val = numba.cpython.builtins.get_type_min_value(in_arr.dtype)
                for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                    val = max(val, in_arr[i])
                return val
    else:
        def max_1(in_arr):
            numba.parfors.parfor.init_prange()
            max_checker(len(in_arr))
            val = numba.cpython.builtins.get_type_min_value(in_arr.dtype)
            for i in numba.pndindex(in_arr.shape):
                val = max(val, in_arr[i])
            return val
    return max_1

def argmin_parallel_impl(in_arr):
    numba.parfors.parfor.init_prange()
    argmin_checker(len(in_arr))
    A = in_arr.ravel()
    init_val = numba.cpython.builtins.get_type_max_value(A.dtype)
    ival = typing.builtins.IndexValue(0, init_val)
    for i in numba.parfors.parfor.internal_prange(len(A)):
        curr_ival = typing.builtins.IndexValue(i, A[i])
        ival = min(ival, curr_ival)
    return ival.index

def argmax_parallel_impl(in_arr):
    numba.parfors.parfor.init_prange()
    argmax_checker(len(in_arr))
    A = in_arr.ravel()
    init_val = numba.cpython.builtins.get_type_min_value(A.dtype)
    ival = typing.builtins.IndexValue(0, init_val)
    for i in numba.parfors.parfor.internal_prange(len(A)):
        curr_ival = typing.builtins.IndexValue(i, A[i])
        ival = max(ival, curr_ival)
    return ival.index

def dotvv_parallel_impl(a, b):
    numba.parfors.parfor.init_prange()
    l = a.shape[0]
    m = b.shape[0]
    # TODO: investigate assert_equiv
    #assert_equiv("sizes of l, m do not match", l, m)
    s = 0
    for i in numba.parfors.parfor.internal_prange(l):
        s += a[i] * b[i]
    return s

def dotvm_parallel_impl(a, b):
    numba.parfors.parfor.init_prange()
    l = a.shape
    m, n = b.shape
    # TODO: investigate assert_equiv
    #assert_equiv("Sizes of l, m do not match", l, m)
    c = np.zeros(n, a.dtype)
    # TODO: evaluate dotvm implementation options
    #for i in prange(n):
    #    s = 0
    #    for j in range(m):
    #        s += a[j] * b[j, i]
    #    c[i] = s
    for i in numba.parfors.parfor.internal_prange(m):
        c += a[i] * b[i, :]
    return c

def dotmv_parallel_impl(a, b):
    numba.parfors.parfor.init_prange()
    m, n = a.shape
    l = b.shape
    # TODO: investigate assert_equiv
    #assert_equiv("sizes of n, l do not match", n, l)
    c = np.empty(m, a.dtype)
    for i in numba.parfors.parfor.internal_prange(m):
        s = 0
        for j in range(n):
            s += a[i, j] * b[j]
        c[i] = s
    return c

def dot_parallel_impl(return_type, atyp, btyp):
    # Note that matrix matrix multiply is not translated.
    if (isinstance(atyp, types.npytypes.Array) and
        isinstance(btyp, types.npytypes.Array)):
        if atyp.ndim == btyp.ndim == 1:
            return dotvv_parallel_impl
        # TODO: evaluate support for dotvm and enable
        #elif atyp.ndim == 1 and btyp.ndim == 2:
        #    return dotvm_parallel_impl
        elif atyp.ndim == 2 and btyp.ndim == 1:
            return dotmv_parallel_impl

def sum_parallel_impl(return_type, arg):
    zero = return_type(0)

    if arg.ndim == 0:
        def sum_1(in_arr):
            return in_arr[()]
    elif arg.ndim == 1:
        def sum_1(in_arr):
            numba.parfors.parfor.init_prange()
            val = zero
            for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                val += in_arr[i]
            return val
    else:
        def sum_1(in_arr):
            numba.parfors.parfor.init_prange()
            val = zero
            for i in numba.pndindex(in_arr.shape):
                val += in_arr[i]
            return val
    return sum_1

def prod_parallel_impl(return_type, arg):
    one = return_type(1)

    if arg.ndim == 0:
        def prod_1(in_arr):
            return in_arr[()]
    elif arg.ndim == 1:
        def prod_1(in_arr):
            numba.parfors.parfor.init_prange()
            val = one
            for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                val *= in_arr[i]
            return val
    else:
        def prod_1(in_arr):
            numba.parfors.parfor.init_prange()
            val = one
            for i in numba.pndindex(in_arr.shape):
                val *= in_arr[i]
            return val
    return prod_1


def mean_parallel_impl(return_type, arg):
    # can't reuse sum since output type is different
    zero = return_type(0)

    if arg.ndim == 0:
        def mean_1(in_arr):
            return in_arr[()]
    elif arg.ndim == 1:
        def mean_1(in_arr):
            numba.parfors.parfor.init_prange()
            val = zero
            for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                val += in_arr[i]
            return val/len(in_arr)
    else:
        def mean_1(in_arr):
            numba.parfors.parfor.init_prange()
            val = zero
            for i in numba.pndindex(in_arr.shape):
                val += in_arr[i]
            return val/in_arr.size
    return mean_1

def var_parallel_impl(return_type, arg):
    if arg.ndim == 0:
        def var_1(in_arr):
            return 0
    elif arg.ndim == 1:
        def var_1(in_arr):
            # Compute the mean
            m = in_arr.mean()
            # Compute the sum of square diffs
            numba.parfors.parfor.init_prange()
            ssd = 0
            for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                val = in_arr[i] - m
                ssd += np.real(val * np.conj(val))
            return ssd / len(in_arr)
    else:
        def var_1(in_arr):
            # Compute the mean
            m = in_arr.mean()
            # Compute the sum of square diffs
            numba.parfors.parfor.init_prange()
            ssd = 0
            for i in numba.pndindex(in_arr.shape):
                val = in_arr[i] - m
                ssd += np.real(val * np.conj(val))
            return ssd / in_arr.size
    return var_1

def std_parallel_impl(return_type, arg):
    def std_1(in_arr):
        return in_arr.var() ** 0.5
    return std_1

def arange_parallel_impl(return_type, *args, dtype=None):
    inferred_dtype = as_dtype(return_type.dtype)

    def arange_1(stop):
        return np.arange(0, stop, 1, inferred_dtype)

    def arange_1_dtype(stop, dtype):
        return np.arange(0, stop, 1, dtype)

    def arange_2(start, stop):
        return np.arange(start, stop, 1, inferred_dtype)

    def arange_2_dtype(start, stop, dtype):
        return np.arange(start, stop, 1, dtype)

    def arange_3(start, stop, step):
        return np.arange(start, stop, step, inferred_dtype)

    def arange_3_dtype(start, stop, step, dtype):
        return np.arange(start, stop, step, dtype)

    if any(isinstance(a, types.Complex) for a in args):
        def arange_4(start, stop, step, dtype):
            numba.parfors.parfor.init_prange()
            nitems_c = (stop - start) / step
            nitems_r = math.ceil(nitems_c.real)
            nitems_i = math.ceil(nitems_c.imag)
            nitems = int(max(min(nitems_i, nitems_r), 0))
            arr = np.empty(nitems, dtype)
            for i in numba.parfors.parfor.internal_prange(nitems):
                arr[i] = start + i * step
            return arr
    else:
        def arange_4(start, stop, step, dtype):
            numba.parfors.parfor.init_prange()
            nitems_r = math.ceil((stop - start) / step)
            nitems = int(max(nitems_r, 0))
            arr = np.empty(nitems, dtype)
            val = start
            for i in numba.parfors.parfor.internal_prange(nitems):
                arr[i] = start + i * step
            return arr

    if len(args) == 1:
        return arange_1 if dtype is None else arange_1_dtype
    elif len(args) == 2:
        return arange_2  if dtype is None else arange_2_dtype
    elif len(args) == 3:
        return arange_3  if dtype is None else arange_3_dtype
    elif len(args) == 4:
        return arange_4
    else:
        raise ValueError("parallel arange with types {}".format(args))

def linspace_parallel_impl(return_type, *args):
    dtype = as_dtype(return_type.dtype)

    def linspace_2(start, stop):
        return np.linspace(start, stop, 50)

    def linspace_3(start, stop, num):
        numba.parfors.parfor.init_prange()
        arr = np.empty(num, dtype)
        div = num - 1
        delta = stop - start
        arr[0] = start
        for i in numba.parfors.parfor.internal_prange(num):
            arr[i] = start + delta * (i / div)
        return arr

    if len(args) == 2:
        return linspace_2
    elif len(args) == 3:
        return linspace_3
    else:
        raise ValueError("parallel linspace with types {}".format(args))

swap_functions_map = {
    ('argmin', 'numpy'): lambda r,a: argmin_parallel_impl,
    ('argmax', 'numpy'): lambda r,a: argmax_parallel_impl,
    ('min', 'numpy'): min_parallel_impl,
    ('max', 'numpy'): max_parallel_impl,
    ('amin', 'numpy'): min_parallel_impl,
    ('amax', 'numpy'): max_parallel_impl,
    ('sum', 'numpy'): sum_parallel_impl,
    ('prod', 'numpy'): prod_parallel_impl,
    ('mean', 'numpy'): mean_parallel_impl,
    ('var', 'numpy'): var_parallel_impl,
    ('std', 'numpy'): std_parallel_impl,
    ('dot', 'numpy'): dot_parallel_impl,
    ('arange', 'numpy'): arange_parallel_impl,
    ('linspace', 'numpy'): linspace_parallel_impl,
}

def fill_parallel_impl(return_type, arr, val):
    """Parallel implementation of ndarray.fill.  The array on
       which to operate is retrieved from get_call_name and
       is passed along with the value to fill.
    """
    if arr.ndim == 1:
        def fill_1(in_arr, val):
            numba.parfors.parfor.init_prange()
            for i in numba.parfors.parfor.internal_prange(len(in_arr)):
                in_arr[i] = val
            return None
    else:
        def fill_1(in_arr, val):
            numba.parfors.parfor.init_prange()
            for i in numba.pndindex(in_arr.shape):
                in_arr[i] = val
            return None
    return fill_1

replace_functions_ndarray = {
    'fill': fill_parallel_impl,
}

@register_jitable
def max_checker(arr_size):
    if arr_size == 0:
        raise ValueError(("zero-size array to reduction operation "
                            "maximum which has no identity"))

@register_jitable
def min_checker(arr_size):
    if arr_size == 0:
        raise ValueError(("zero-size array to reduction operation "
                            "minimum which has no identity"))

@register_jitable
def argmin_checker(arr_size):
    if arr_size == 0:
        raise ValueError("attempt to get argmin of an empty sequence")

@register_jitable
def argmax_checker(arr_size):
    if arr_size == 0:
        raise ValueError("attempt to get argmax of an empty sequence")

checker_impl = namedtuple('checker_impl', ['name', 'func'])

replace_functions_checkers_map = {
    ('argmin', 'numpy') : checker_impl('argmin_checker', argmin_checker),
    ('argmax', 'numpy') : checker_impl('argmax_checker', argmax_checker),
    ('min', 'numpy') : checker_impl('min_checker', min_checker),
    ('max', 'numpy') : checker_impl('max_checker', max_checker),
    ('amin', 'numpy') : checker_impl('min_checker', min_checker),
    ('amax', 'numpy') : checker_impl('max_checker', max_checker),
}


class LoopNest(object):

    '''The LoopNest class holds information of a single loop including
    the index variable (of a non-negative integer value), and the
    range variable, e.g. range(r) is 0 to r-1 with step size 1.
    '''

    def __init__(self, index_variable, start, stop, step):
        self.index_variable = index_variable
        self.start = start
        self.stop = stop
        self.step = step


    def __repr__(self):
        return ("LoopNest(index_variable = {}, range = ({}, {}, {}))".
                format(self.index_variable, self.start, self.stop, self.step))

    def list_vars(self):
       all_uses = []
       all_uses.append(self.index_variable)
       if isinstance(self.start, ir.Var):
           all_uses.append(self.start)
       if isinstance(self.stop, ir.Var):
           all_uses.append(self.stop)
       if isinstance(self.step, ir.Var):
           all_uses.append(self.step)
       return all_uses

class Parfor(ir.Expr, ir.Stmt):

    id_counter = 0

    def __init__(
            self,
            loop_nests,
            init_block,
            loop_body,
            loc,
            index_var,
            equiv_set,
            pattern,
            flags,
            *,  #only specify the options below by keyword
            no_sequential_lowering=False,
            races=set()):
        super(Parfor, self).__init__(
            op='parfor',
            loc=loc
        )

        self.id = type(self).id_counter
        type(self).id_counter += 1
        #self.input_info  = input_info
        #self.output_info = output_info
        self.loop_nests = loop_nests
        self.init_block = init_block
        self.loop_body = loop_body
        self.index_var = index_var
        self.params = None  # filled right before parallel lowering
        self.equiv_set = equiv_set
        # The parallel patterns this parfor was generated from and their options
        # for example, a parfor could be from the stencil pattern with
        # the neighborhood option
        assert len(pattern) > 1
        self.patterns = [pattern]
        self.flags = flags
        # if True, this parfor shouldn't be lowered sequentially even with the
        # sequential lowering option
        self.no_sequential_lowering = no_sequential_lowering
        self.races = races
        self.redvars = []
        self.reddict = {}
        # If the lowerer is None then the standard lowerer will be used.
        # This can be set to a function to have that function act as the lowerer
        # for this parfor.  This lowerer field will also prevent parfors from
        # being fused unless they have use the same lowerer.
        self.lowerer = None
        if config.DEBUG_ARRAY_OPT_STATS:
            fmt = 'Parallel for-loop #{} is produced from pattern \'{}\' at {}'
            print(fmt.format(
                  self.id, pattern, loc))

    def __repr__(self):
        return "id=" + str(self.id) + repr(self.loop_nests) + \
            repr(self.loop_body) + repr(self.index_var)

    def get_loop_nest_vars(self):
        return [x.index_variable for x in self.loop_nests]

    def list_vars(self):
        """list variables used (read/written) in this parfor by
        traversing the body and combining block uses.
        """
        all_uses = []
        for l, b in self.loop_body.items():
            for stmt in b.body:
                all_uses += stmt.list_vars()

        for loop in self.loop_nests:
            all_uses += loop.list_vars()

        for stmt in self.init_block.body:
            all_uses += stmt.list_vars()

        return all_uses

    def get_shape_classes(self, var, typemap=None):
        """get the shape classes for a given variable.
        If a typemap is specified then use it for type resolution
        """
        # We get shape classes from the equivalence set but that
        # keeps its own typemap at a time prior to lowering.  So
        # if something is added during lowering then we can pass
        # in a type map to use.  We temporarily replace the
        # equivalence set typemap, do the work and then restore
        # the original on the way out.
        if typemap is not None:
            save_typemap = self.equiv_set.typemap
            self.equiv_set.typemap = typemap
        res = self.equiv_set.get_shape_classes(var)
        if typemap is not None:
            self.equiv_set.typemap = save_typemap
        return res

    def dump(self, file=None):
        file = file or sys.stdout
        print(("begin parfor {}".format(self.id)).center(20, '-'), file=file)
        print("index_var = ", self.index_var, file=file)
        print("params = ", self.params, file=file)
        print("races = ", self.races, file=file)
        for loopnest in self.loop_nests:
            print(loopnest, file=file)
        print("init block:", file=file)
        self.init_block.dump(file)
        for offset, block in sorted(self.loop_body.items()):
            print('label %s:' % (offset,), file=file)
            block.dump(file)
        print(("end parfor {}".format(self.id)).center(20, '-'), file=file)

    def validate_params(self, typemap):
        """
        Check that Parfors params are of valid types.
        """
        if self.params is None:
            msg = ("Cannot run parameter validation on a Parfor with params "
                   "not set")
            raise ValueError(msg)
        for p in self.params:
            ty = typemap.get(p)
            if ty is None:
                msg = ("Cannot validate parameter %s, there is no type "
                       "information available")
                raise ValueError(msg)
            if isinstance(ty, types.BaseTuple):
                if ty.count > config.PARFOR_MAX_TUPLE_SIZE:
                    msg = ("Use of a tuple (%s) of length %d in a parallel region "
                           "exceeds the maximum supported tuple size.  Since "
                           "Generalized Universal Functions back parallel regions "
                           "and those do not support tuples, tuples passed to "
                           "parallel regions are unpacked if their size is below "
                           "a certain threshold, currently configured to be %d. "
                           "This threshold can be modified using the Numba "
                           "environment variable NUMBA_PARFOR_MAX_TUPLE_SIZE.")
                    raise errors.UnsupportedParforsError(msg %
                              (p, ty.count, config.PARFOR_MAX_TUPLE_SIZE),
                              self.loc)


def _analyze_parfor(parfor, equiv_set, typemap, array_analysis):
    """Recursive array analysis for parfor nodes.
    """
    func_ir = array_analysis.func_ir
    parfor_blocks = wrap_parfor_blocks(parfor)
    # Since init_block get label 0 after wrap, we need to save
    # the equivset for the real block label 0.
    backup_equivset = array_analysis.equiv_sets.get(0, None)
    array_analysis.run(parfor_blocks, equiv_set)
    unwrap_parfor_blocks(parfor, parfor_blocks)
    parfor.equiv_set = array_analysis.equiv_sets[0]
    # Restore equivset for block 0 after parfor is unwrapped
    if backup_equivset:
        array_analysis.equiv_sets[0] = backup_equivset
    return [], []

array_analysis.array_analysis_extensions[Parfor] = _analyze_parfor

class ParforDiagnostics(object):
    """Holds parfor diagnostic info, this is accumulated throughout the
    PreParforPass and ParforPass, also in the closure inlining!
    """
    def __init__(self):
        # holds ref to the function for which this is providing diagnostics
        self.func = None
        # holds a map of the replaced functions
        self.replaced_fns = dict()
        # used to identify "internal" parfor functions
        self.internal_name = '__numba_parfor_gufunc'
        self.fusion_info = defaultdict(list)
        self.nested_fusion_info = defaultdict(list)
        self.fusion_reports = []
        self.hoist_info = {}
        self.has_setup = False

    def setup(self, func_ir, fusion_enabled):
        self.func_ir = func_ir
        self.name = self.func_ir.func_id.func_qualname
        self.line = self.func_ir.loc
        self.fusion_enabled = fusion_enabled
        if self.internal_name in self.name:
            self.purpose = 'Internal parallel function'
        else:
            self.purpose = 'Function %s, %s' % (self.name, self.line)
        # we store a reference to the parfors prior to fusion etc, the parfors
        # do get mangled in the fusion process but in a predetermined manner
        # and by holding a reference here the "before" state can be printed
        self.initial_parfors = self.get_parfors()
        self.has_setup = True

    @property
    def has_setup(self):
        return self._has_setup

    @has_setup.setter
    def has_setup(self, state):
        self._has_setup = state

    def count_parfors(self, blocks=None):
        return len(self.get_parfors())

    def _get_nested_parfors(self, parfor, parfors_list):
        blocks = wrap_parfor_blocks(parfor)
        self._get_parfors(blocks, parfors_list)
        unwrap_parfor_blocks(parfor)

    def _get_parfors(self, blocks, parfors_list):
        for label, blk in blocks.items():
            for stmt in blk.body:
                if isinstance(stmt, Parfor):
                    parfors_list.append(stmt)
                    self._get_nested_parfors(stmt, parfors_list)

    def get_parfors(self):
        parfors_list = []
        self._get_parfors(self.func_ir.blocks, parfors_list)
        return parfors_list

    def hoisted_allocations(self):
        allocs = []
        for pf_id, data in self.hoist_info.items():
            stmt = data.get('hoisted', [])
            for inst in stmt:
                if isinstance(inst.value, ir.Expr):
                    if inst.value.op == 'call':
                        call = guard(find_callname, self.func_ir, inst.value)
                        if call is not None and call == ('empty', 'numpy'):
                            allocs.append(inst)
        return allocs

    def compute_graph_info(self, _a):
        """
        compute adjacency list of the fused loops
        and find the roots in of the lists
        """
        a = copy.deepcopy(_a)
        if a == {}:
            return [], set()

        vtx = set()
        for v in a.values():
            for x in v:
                vtx.add(x)

        # find roots
        potential_roots = set(a.keys())
        roots = potential_roots - vtx
        if roots is None:
            roots = set()

        # populate rest of adjacency list
        not_roots = set()
        for x in range(max(set(a.keys()).union(vtx)) + 1):
            val = a.get(x)
            if val is not None:
                a[x] = val
            elif val == []:
                not_roots.add(x) # debug only
            else:
                a[x] = []


        # fold adjacency list into an actual list ordered
        # by vtx
        l = []
        for x in sorted(a.keys()):
            l.append(a[x])

        return l, roots #, not_roots

    def get_stats(self, fadj, nadj, root):
        """
        Computes the number of fused and serialized loops
        based on a fusion adjacency list `fadj` and a nested
        parfors adjacency list `nadj` for the root, `root`
        """
        def count_root(fadj, nadj, root, nfused, nserial):
            for k in nadj[root]:
                nserial += 1
                if nadj[k] == []:
                    nfused += len(fadj[k])
                else:
                    nf, ns = count_root(fadj, nadj, k, nfused, nserial)
                    nfused += nf
                    nserial = ns
            return nfused, nserial
        nfused, nserial = count_root(fadj, nadj, root, 0, 0)
        return nfused, nserial

    def reachable_nodes(self, adj, root):
        """
        returns a list of nodes reachable in an adjacency list from a
        specified root
        """
        fusers = []
        fusers.extend(adj[root])
        for k in adj[root]:
            if adj[k] != []:
                fusers.extend(self.reachable_nodes(adj, k))
        return fusers

    def sort_pf_by_line(self, pf_id, parfors_simple):
        """
        pd_id - the parfors id
        parfors_simple - the simple parfors map
        """
        # this sorts parfors by source line number
        pf = parfors_simple[pf_id][0]
        pattern = pf.patterns[0]
        line = max(0, pf.loc.line - 1) # why are these out by 1 ?!
        filename = self.func_ir.loc.filename
        nadj, nroots = self.compute_graph_info(self.nested_fusion_info)
        fadj, froots = self.compute_graph_info(self.fusion_info)
        graphs = [nadj, fadj]

        # If the parfor is internal, like internal prange, then the
        # default line number is from its location in the numba source
        # To get a more accurate line number, this first checks the
        # adjacency graph for fused parfors that might not be internal
        # and uses the minimum line number from there. If that fails
        # (case where there's just a single internal parfor) the IR
        # is walked backwards from the parfor location and the first non
        # parfor statement line number is used.
        if isinstance(pattern, tuple):
            if pattern[1] == 'internal':
                reported_loc = pattern[2][1]
                if reported_loc.filename == filename:
                    return max(0, reported_loc.line - 1)
                else:
                    # first recurse and check the adjacency list for
                    # something that is not an in internal parfor
                    tmp = []
                    for adj in graphs:
                        if adj: # graph may be empty, e.g. no nesting
                            for k in adj[pf_id]:
                                tmp.append(self.sort_pf_by_line(k, parfors_simple))
                            if tmp:
                                return max(0, min(tmp) - 1)
                    # second run through the parfor block to see if there's
                    # and reference to a line number in the user source
                    for blk in pf.loop_body.values():
                        for stmt in blk.body:
                            if stmt.loc.filename == filename:
                                return max(0, stmt.loc.line - 1)
                    # finally run through the func_ir and look for the
                    # first non-parfor statement prior to this one and
                    # grab the line from that
                    for blk in self.func_ir.blocks.values():
                        try:
                            idx = blk.body.index(pf)
                            for i in range(idx - 1, 0, -1):
                                stmt = blk.body[i]
                                if not isinstance(stmt, Parfor):
                                    line = max(0, stmt.loc.line - 1)
                                    break
                        except ValueError:
                            pass
        return line

    def get_parfors_simple(self, print_loop_search):
        parfors_simple = dict()

        # print in line order, parfors loop id is based on discovery order
        for pf in sorted(self.initial_parfors, key=lambda x: x.loc.line):
            # use 0 here, the parfors are mutated by the time this routine
            # is called, however, fusion appends the patterns so we can just
            # pull in the first as a "before fusion" emulation
            r_pattern = pf.patterns[0]
            pattern = pf.patterns[0]
            loc = pf.loc
            if isinstance(pattern, tuple):
                if pattern[0] == 'prange':
                    if pattern[1] == 'internal':
                        replfn = '.'.join(reversed(list(pattern[2][0])))
                        loc = pattern[2][1]
                        r_pattern = '%s %s' % (replfn, '(internal parallel version)')
                    elif pattern[1] == 'user':
                        r_pattern = "user defined prange"
                    elif pattern[1] == 'pndindex':
                        r_pattern = "internal pndindex" #FIXME: trace this!
                    else:
                        assert 0
            fmt = 'Parallel for-loop #%s: is produced from %s:\n    %s\n \n'
            if print_loop_search:
                print_wrapped(fmt % (pf.id, loc, r_pattern))
            parfors_simple[pf.id] = (pf, loc, r_pattern)
        return parfors_simple

    def get_all_lines(self, parfors_simple):
        # ensure adjacency lists are the same size for both sets of info
        # (nests and fusion may not traverse the same space, for
        # convenience [] is used as a condition to halt recursion)
        fadj, froots = self.compute_graph_info(self.fusion_info)
        nadj, _nroots = self.compute_graph_info(self.nested_fusion_info)

        if len(fadj) > len(nadj):
            lim = len(fadj)
            tmp = nadj
        else:
            lim = len(nadj)
            tmp = fadj
        for x in range(len(tmp), lim):
            tmp.append([])

        # This computes the roots of true loop nests (i.e. loops containing
        # loops opposed to just a loop that's a root).
        nroots = set()
        if _nroots:
            for r in _nroots:
                if nadj[r] != []:
                    nroots.add(r)
        all_roots = froots ^ nroots

        # This computes all the parfors at the top level that are either:
        # - roots of loop fusion
        # - roots of true loop nests
        # it then combines these based on source line number for ease of
        # producing output ordered in a manner similar to the code structure
        froots_lines = {}
        for x in froots:
            line = self.sort_pf_by_line(x, parfors_simple)
            froots_lines[line] = 'fuse', x, fadj

        nroots_lines = {}
        for x in nroots:
            line = self.sort_pf_by_line(x, parfors_simple)
            nroots_lines[line] = 'nest', x, nadj

        all_lines = froots_lines.copy()
        all_lines.update(nroots_lines)
        return all_lines

    def source_listing(self, parfors_simple, purpose_str):
        filename = self.func_ir.loc.filename
        count = self.count_parfors()
        func_name = self.func_ir.func_id.func
        try:
            lines = inspect.getsource(func_name).splitlines()
        except OSError: # generated function
            lines = None
        if lines and parfors_simple:
            src_width = max([len(x) for x in lines])
            map_line_to_pf = defaultdict(list) # parfors can alias lines
            for k, v in parfors_simple.items():
                # TODO: do a better job of tracking parfors that are not in
                # this file but are referred to, e.g. np.arange()
                if parfors_simple[k][1].filename == filename:
                    match_line = self.sort_pf_by_line(k, parfors_simple)
                    map_line_to_pf[match_line].append(str(k))

            max_pf_per_line = max([1] + [len(x) for x in map_line_to_pf.values()])
            width = src_width + (1 + max_pf_per_line * (len(str(count)) + 2))
            newlines = []
            newlines.append('\n')
            newlines.append('Parallel loop listing for %s' % purpose_str)
            newlines.append(width * '-' + '|loop #ID')
            fmt = '{0:{1}}| {2}'
            # why are these off by 1?
            lstart = max(0, self.func_ir.loc.line - 1)
            for no, line in enumerate(lines, lstart):
                pf_ids = map_line_to_pf.get(no, None)
                if pf_ids is not None:
                    pfstr = '#' + ', '.join(pf_ids)
                else:
                    pfstr = ''
                stripped = line.strip('\n')
                srclen = len(stripped)
                if pf_ids:
                    l = fmt.format(width * '-', width, pfstr)
                else:
                    l = fmt.format(width * ' ', width, pfstr)
                newlines.append(stripped + l[srclen:])
            print('\n'.join(newlines))
        else:
            print("No source available")

    def print_unoptimised(self, lines):
        # This prints the unoptimised parfors state
        sword = '+--'
        fac = len(sword)
        fadj, froots = self.compute_graph_info(self.fusion_info)
        nadj, _nroots = self.compute_graph_info(self.nested_fusion_info)

        if len(fadj) > len(nadj):
            lim = len(fadj)
            tmp = nadj
        else:
            lim = len(nadj)
            tmp = fadj
        for x in range(len(tmp), lim):
            tmp.append([])

        def print_nest(fadj_, nadj_, theroot, reported, region_id):
            def print_g(fadj_, nadj_, nroot, depth):
                print_wrapped(fac * depth * ' ' + '%s%s %s' % (sword, nroot, '(parallel)'))
                for k in nadj_[nroot]:
                    if nadj_[k] == []:
                        msg = []
                        msg.append(fac * (depth + 1) * ' ' + '%s%s %s' % (sword, k, '(parallel)'))
                        if fadj_[k] != [] and k not in reported:
                            fused = self.reachable_nodes(fadj_, k)
                            for i in fused:
                                msg.append(fac * (depth + 1) * ' ' + '%s%s %s' % (sword, i, '(parallel)'))
                        reported.append(k)
                        print_wrapped('\n'.join(msg))
                    else:
                        print_g(fadj_, nadj_, k, depth + 1)

            if nadj_[theroot] != []:
                print_wrapped("Parallel region %s:" % region_id)
                print_g(fadj_, nadj_, theroot, 0)
                print("\n")
                region_id = region_id + 1
            return region_id

        def print_fuse(ty, pf_id, adj, depth, region_id):
            msg = []
            print_wrapped("Parallel region %s:" % region_id)
            msg.append(fac * depth * ' ' + '%s%s %s' % (sword, pf_id, '(parallel)'))
            if adj[pf_id] != []:
                fused = sorted(self.reachable_nodes(adj, pf_id))
                for k in fused:
                    msg.append(fac * depth * ' ' + '%s%s %s' % (sword, k, '(parallel)'))
            region_id = region_id + 1
            print_wrapped('\n'.join(msg))
            print("\n")
            return region_id

        # Walk the parfors by src line and print optimised structure
        region_id = 0
        reported = []
        for line, info in sorted(lines.items()):
            opt_ty, pf_id, adj = info
            if opt_ty == 'fuse':
                if pf_id not in reported:
                    region_id = print_fuse('f', pf_id, adj, 0, region_id)
            elif opt_ty == 'nest':
                region_id = print_nest(fadj, nadj, pf_id, reported, region_id)
            else:
                assert 0

    def print_optimised(self, lines):
        # This prints the optimised output based on the transforms that
        # occurred during loop fusion and rewriting of loop nests
        sword = '+--'
        fac = len(sword)
        fadj, froots = self.compute_graph_info(self.fusion_info)
        nadj, _nroots = self.compute_graph_info(self.nested_fusion_info)

        if len(fadj) > len(nadj):
            lim = len(fadj)
            tmp = nadj
        else:
            lim = len(nadj)
            tmp = fadj
        for x in range(len(tmp), lim):
            tmp.append([])

        summary = dict()
        # region : {fused, serialized}

        def print_nest(fadj_, nadj_, theroot, reported, region_id):
            def print_g(fadj_, nadj_, nroot, depth):
                for k in nadj_[nroot]:
                    msg = fac * depth * ' ' + '%s%s %s' % (sword, k, '(serial')
                    if nadj_[k] == []:
                        fused = []
                        if fadj_[k] != [] and k not in reported:
                            fused = sorted(self.reachable_nodes(fadj_, k))
                            msg += ", fused with loop(s): "
                            msg += ', '.join([str(x) for x in fused])
                        msg += ')'
                        reported.append(k)
                        print_wrapped(msg)
                        summary[region_id]['fused'] += len(fused)
                    else:
                        print_wrapped(msg + ')')
                        print_g(fadj_, nadj_, k, depth + 1)
                    summary[region_id]['serialized'] += 1

            if nadj_[theroot] != []:
                print_wrapped("Parallel region %s:" % region_id)
                print_wrapped('%s%s %s' % (sword, theroot, '(parallel)'))
                summary[region_id] = {'root': theroot, 'fused': 0, 'serialized': 0}
                print_g(fadj_, nadj_, theroot, 1)
                print("\n")
                region_id = region_id + 1
            return region_id

        def print_fuse(ty, pf_id, adj, depth, region_id):
            print_wrapped("Parallel region %s:" % region_id)
            msg = fac * depth * ' ' + '%s%s %s' % (sword, pf_id, '(parallel')
            fused = []
            if adj[pf_id] != []:
                fused = sorted(self.reachable_nodes(adj, pf_id))
                msg += ", fused with loop(s): "
                msg += ', '.join([str(x) for x in fused])

            summary[region_id] = {'root': pf_id, 'fused': len(fused), 'serialized': 0}
            msg += ')'
            print_wrapped(msg)
            print("\n")
            region_id = region_id + 1
            return region_id

        # Walk the parfors by src line and print optimised structure
        region_id = 0
        reported = []
        for line, info in sorted(lines.items()):
            opt_ty, pf_id, adj = info
            if opt_ty == 'fuse':
                if pf_id not in reported:
                    region_id = print_fuse('f', pf_id, adj, 0, region_id)
            elif opt_ty == 'nest':
                region_id = print_nest(fadj, nadj, pf_id, reported, region_id)
            else:
                assert 0

        # print the summary of the fuse/serialize rewrite
        if summary:
            for k, v in sorted(summary.items()):
                msg = ('\n \nParallel region %s (loop #%s) had %s '
                    'loop(s) fused')
                root = v['root']
                fused = v['fused']
                serialized = v['serialized']
                if serialized != 0:
                    msg += (' and %s loop(s) '
                    'serialized as part of the larger '
                    'parallel loop (#%s).')
                    print_wrapped(msg % (k, root, fused, serialized, root))
                else:
                    msg += '.'
                    print_wrapped(msg % (k, root, fused))
        else:
            print_wrapped("Parallel structure is already optimal.")

    def allocation_hoist(self):
        found = False
        print('Allocation hoisting:')
        for pf_id, data in self.hoist_info.items():
            stmt = data.get('hoisted', [])
            for inst in stmt:
                if isinstance(inst.value, ir.Expr):
                    try:
                        attr = inst.value.attr
                        if attr == 'empty':
                            msg = ("The memory allocation derived from the "
                                "instruction at %s is hoisted out of the "
                                "parallel loop labelled #%s (it will be "
                                "performed before the loop is executed and "
                                "reused inside the loop):")
                            loc = inst.loc
                            print_wrapped(msg % (loc, pf_id))
                            try:
                                path = os.path.relpath(loc.filename)
                            except ValueError:
                                path = os.path.abspath(loc.filename)
                            lines = linecache.getlines(path)
                            if lines and loc.line:
                                print_wrapped("   Allocation:: " + lines[0 if loc.line < 2 else loc.line - 1].strip())
                            print_wrapped("    - numpy.empty() is used for the allocation.\n")
                            found = True
                    except (KeyError, AttributeError):
                        pass
        if not found:
            print_wrapped('No allocation hoisting found')

    def instruction_hoist(self):
        print("")
        print('Instruction hoisting:')
        hoist_info_printed = False
        if self.hoist_info:
            for pf_id, data in self.hoist_info.items():
                hoisted = data.get('hoisted', None)
                not_hoisted = data.get('not_hoisted', None)
                if not hoisted and not not_hoisted:
                    print("loop #%s has nothing to hoist." % pf_id)
                    continue

                print("loop #%s:" % pf_id)
                if hoisted:
                    print("  Has the following hoisted:")
                    [print("    %s" % y) for y in hoisted]
                    hoist_info_printed = True
                if not_hoisted:
                    print("  Failed to hoist the following:")
                    [print("    %s: %s" % (y, x)) for x, y in not_hoisted]
                    hoist_info_printed = True
        if not hoist_info_printed:
            print_wrapped('No instruction hoisting found')
        print_wrapped(80 * '-')

    def dump(self, level=1):
        if not self.has_setup:
            raise RuntimeError("self.setup has not been called")
        name = self.func_ir.func_id.func_qualname
        line = self.func_ir.loc
        if self.internal_name in name:
            purpose_str = 'Internal parallel functions '
            purpose = 'internal'
        else:
            purpose_str = ' Function %s, %s ' % (name, line)
            purpose = 'user'

        print_loop_search = False
        print_source_listing = False
        print_fusion_search = False
        print_fusion_summary = False
        print_loopnest_rewrite = False
        print_pre_optimised = False
        print_post_optimised = False
        print_allocation_hoist = False
        print_instruction_hoist = False
        print_internal = False

        # each level switches on progressively more output
        if level in (1, 2, 3, 4):
            print_source_listing = True
            print_post_optimised = True
        else:
            raise ValueError("Report level unknown, should be one of 1, 2, 3, 4")

        if level in (2, 3, 4):
            print_pre_optimised = True

        if level in (3, 4):
            print_allocation_hoist = True

        if level == 3:
            print_fusion_summary = True
            print_loopnest_rewrite = True

        if level == 4:
            print_fusion_search = True
            print_instruction_hoist = True
            print_internal = True

        if purpose == 'internal' and not print_internal:
            return

        print_wrapped('\n ')
        print_wrapped(_termwidth * "=")
        print_wrapped((" Parallel Accelerator Optimizing: %s " % purpose_str).center(_termwidth, '='))
        print_wrapped(_termwidth * "=")
        print_wrapped("")

#----------- search section
        if print_loop_search:
            print_wrapped('Looking for parallel loops'.center(_termwidth, '-'))
        parfors_simple = self.get_parfors_simple(print_loop_search)

        count = self.count_parfors()
        if print_loop_search:
            print_wrapped("\nFound %s parallel loops." % count)
            print_wrapped('-' * _termwidth)

#----------- augmented source section
        filename = self.func_ir.loc.filename
        try:
            # Try to get a relative path
            # ipython/jupyter input just returns as filename
            path = os.path.relpath(filename)
        except ValueError:
            # Fallback to absolute path if error occurred in getting the
            # relative path.
            # This may happen on windows if the drive is different
            path = os.path.abspath(filename)

        if print_source_listing:
            self.source_listing(parfors_simple, purpose_str)

#---------- these are used a lot here on in
        sword = '+--'
        parfors = self.get_parfors() # this is the mutated parfors
        parfor_ids = [x.id for x in parfors]
        n_parfors = len(parfor_ids)

#----------- loop fusion section
        if print_fusion_search or print_fusion_summary:
            if not sequential_parfor_lowering:
                print_wrapped(' Fusing loops '.center(_termwidth, '-'))
                msg = ("Attempting fusion of parallel loops (combines loops "
                        "with similar properties)...\n")
                print_wrapped(msg)
            else:
                msg = "Performing sequential lowering of loops...\n"
                print_wrapped(msg)
                print_wrapped(_termwidth * '-')
        # if there are some parfors, print information about them!
        if n_parfors > -1:
            def dump_graph_indented(a, root_msg, node_msg):
                fac = len(sword)
                def print_graph(adj, roots):
                    def print_g(adj, root, depth):
                        for k in adj[root]:
                            print_wrapped(fac * depth * ' ' + '%s%s %s' % (sword, k, node_msg))
                            if adj[k] != []:
                                print_g(adj, k, depth + 1)
                    for r in roots:
                        print_wrapped('%s%s %s' % (sword, r, root_msg))
                        print_g(l, r, 1)
                        print_wrapped("")
                l, roots = self.compute_graph_info(a)
                print_graph(l, roots)

            if print_fusion_search:
                for report in self.fusion_reports:
                    l1, l2, msg = report
                    print_wrapped("  Trying to fuse loops #%s and #%s:" % (l1, l2))
                    print_wrapped("    %s" % msg)

            if self.fusion_info != {}:
                if print_fusion_summary:
                    print_wrapped("\n \nFused loop summary:\n")

                    dump_graph_indented(self.fusion_info, 'has the following loops fused into it:', '(fused)')

            if print_fusion_summary:
                if self.fusion_enabled:
                    after_fusion = "Following the attempted fusion of parallel for-loops"
                else:
                    after_fusion = "With fusion disabled"

                print_wrapped(('\n{} there are {} parallel for-loop(s) (originating from loops labelled: {}).').format(
                        after_fusion, n_parfors, ', '.join(['#%s' % x for x in parfor_ids])))
                print_wrapped(_termwidth * '-')
                print_wrapped("")

#----------- loop nest section
            if print_loopnest_rewrite:
                if self.nested_fusion_info != {}:
                    print_wrapped((" Optimising loop nests ").center(_termwidth, '-'))
                    print_wrapped("Attempting loop nest rewrites (optimising for the largest parallel loops)...\n ")
                    root_msg = 'is a parallel loop'
                    node_msg = '--> rewritten as a serial loop'
                    dump_graph_indented(self.nested_fusion_info, root_msg, node_msg)
                    print_wrapped(_termwidth * '-')
                    print_wrapped("")

#---------- compute various properties and orderings in the data for subsequent use
            all_lines = self.get_all_lines(parfors_simple)

            if print_pre_optimised:
                print(' Before Optimisation '.center(_termwidth,'-'))
                self.print_unoptimised(all_lines)
                print(_termwidth * '-')

            if print_post_optimised:
                print(' After Optimisation '.center(_termwidth,'-'))
                self.print_optimised(all_lines)
                print(_termwidth * '-')
            print_wrapped("")
            print_wrapped(_termwidth * '-')
            print_wrapped("\n ")

#----------- LICM section
            if print_allocation_hoist or print_instruction_hoist:
                print_wrapped('Loop invariant code motion'.center(80, '-'))

            if print_allocation_hoist:
                self.allocation_hoist()

            if print_instruction_hoist:
                self.instruction_hoist()

        else: # there are no parfors
            print_wrapped('Function %s, %s, has no parallel for-loops.'.format(name, line))

    def __str__(self):
        r  = "ParforDiagnostics:\n"
        r += repr(self.replaced_fns)
        return r

    def __repr__(self):
        r  = "ParforDiagnostics"
        return r


class PreParforPass(object):
    """Preprocessing for the Parfor pass. It mostly inlines parallel
    implementations of numpy functions if available.
    """
    def __init__(self, func_ir, typemap, calltypes, typingctx, targetctx,
                 options, swapped={}, replace_functions_map=None):
        self.func_ir = func_ir
        self.typemap = typemap
        self.calltypes = calltypes
        self.typingctx = typingctx
        self.targetctx = targetctx
        self.options = options
        # diagnostics
        self.swapped = swapped
        if replace_functions_map is None:
            replace_functions_map = swap_functions_map
        self.replace_functions_map = replace_functions_map
        self.stats = {
            'replaced_func': 0,
            'replaced_dtype': 0,
        }

    def run(self):
        """Run pre-parfor processing pass.
        """
        # e.g. convert A.sum() to np.sum(A) for easier match and optimization
        canonicalize_array_math(self.func_ir, self.typemap,
                                self.calltypes, self.typingctx)
        if self.options.numpy:
            self._replace_parallel_functions(self.func_ir.blocks)
        self.func_ir.blocks = simplify_CFG(self.func_ir.blocks)

    def _replace_parallel_functions(self, blocks):
        """
        Replace functions with their parallel implementation in
        replace_functions_map if available.
        The implementation code is inlined to enable more optimization.
        """
        swapped = self.swapped
        from numba.core.inline_closurecall import inline_closure_call
        work_list = list(blocks.items())
        while work_list:
            label, block = work_list.pop()
            for i, instr in enumerate(block.body):
                if isinstance(instr, ir.Assign):
                    lhs = instr.target
                    lhs_typ = self.typemap[lhs.name]
                    expr = instr.value
                    if isinstance(expr, ir.Expr) and expr.op == 'call':
                        # Try and inline known calls with their parallel implementations
                        def replace_func():
                            func_def = get_definition(self.func_ir, expr.func)
                            callname = find_callname(self.func_ir, expr)
                            repl_func = self.replace_functions_map.get(callname, None)
                            # Handle method on array type
                            if (repl_func is None and
                                len(callname) == 2 and
                                isinstance(callname[1], ir.Var) and
                                isinstance(self.typemap[callname[1].name],
                                           types.npytypes.Array)):
                                repl_func = replace_functions_ndarray.get(callname[0], None)
                                if repl_func is not None:
                                    # Add the array that the method is on to the arg list.
                                    expr.args.insert(0, callname[1])

                            require(repl_func is not None)
                            typs = tuple(self.typemap[x.name] for x in expr.args)
                            kws_typs = {k: self.typemap[x.name] for k, x in expr.kws}
                            try:
                                new_func =  repl_func(lhs_typ, *typs, **kws_typs)
                            except:
                                new_func = None
                            require(new_func is not None)
                            # bind arguments to the new_func
                            typs = utils.pysignature(new_func).bind(*typs, **kws_typs).args

                            g = copy.copy(self.func_ir.func_id.func.__globals__)
                            g['numba'] = numba
                            g['np'] = numpy
                            g['math'] = math
                            # if the function being inlined has a function
                            # checking the inputs, find it and add it to globals
                            check = replace_functions_checkers_map.get(callname,
                                                                       None)
                            if check is not None:
                                g[check.name] = check.func
                            # inline the parallel implementation
                            new_blocks, _ = inline_closure_call(self.func_ir, g,
                                            block, i, new_func, self.typingctx, self.targetctx,
                                            typs, self.typemap, self.calltypes, work_list)
                            call_table = get_call_table(new_blocks, topological_ordering=False)

                            # find the prange in the new blocks and record it for use in diagnostics
                            for call in call_table:
                                for k, v in call.items():
                                    if v[0] == 'internal_prange':
                                        swapped[k] = [callname, repl_func.__name__, func_def, block.body[i].loc]
                                        break
                            return True
                        if guard(replace_func):
                            self.stats['replaced_func'] += 1
                            break
                    elif (isinstance(expr, ir.Expr) and expr.op == 'getattr' and
                          expr.attr == 'dtype'):
                        # Replace getattr call "A.dtype" with numpy.dtype(<actual type>).
                        # This helps remove superfluous dependencies from parfor.
                        typ = self.typemap[expr.value.name]
                        if isinstance(typ, types.npytypes.Array):
                            # Convert A.dtype to four statements.
                            # 1) Get numpy global.
                            # 2) Create var for known type of array as string
                            #    constant. e.g. 'float64'
                            # 3) Get dtype function from numpy module.
                            # 4) Create var for numpy.dtype(var from #2).

                            # Create var for numpy module.
                            dtype = typ.dtype
                            scope = block.scope
                            loc = instr.loc
                            g_np_var = ir.Var(scope, mk_unique_var("$np_g_var"), loc)
                            self.typemap[g_np_var.name] = types.misc.Module(numpy)
                            g_np = ir.Global('np', numpy, loc)
                            g_np_assign = ir.Assign(g_np, g_np_var, loc)

                            # Create var for the inferred type of the array
                            # e.g., 'float64'
                            dtype_str = str(dtype)
                            if dtype_str == 'bool':
                                dtype_str = 'bool_'
                            typ_var = ir.Var(
                                scope, mk_unique_var("$np_typ_var"), loc)
                            self.typemap[typ_var.name] = types.StringLiteral(
                                dtype_str)
                            typ_var_assign = ir.Assign(
                                ir.Const(dtype_str, loc), typ_var, loc)

                            # Get the dtype function from the numpy module.
                            dtype_attr_var = ir.Var(scope, mk_unique_var("$dtype_attr_var"), loc)
                            temp = find_template(numpy.dtype)
                            tfunc = numba.core.types.Function(temp)
                            tfunc.get_call_type(self.typingctx, (self.typemap[typ_var.name],), {})
                            self.typemap[dtype_attr_var.name] = types.functions.Function(temp)
                            dtype_attr_getattr = ir.Expr.getattr(g_np_var, 'dtype', loc)
                            dtype_attr_assign = ir.Assign(dtype_attr_getattr, dtype_attr_var, loc)

                            # Call numpy.dtype on the statically coded type two steps above.
                            dtype_var = ir.Var(scope, mk_unique_var("$dtype_var"), loc)
                            self.typemap[dtype_var.name] = types.npytypes.DType(dtype)
                            dtype_getattr = ir.Expr.call(dtype_attr_var, [typ_var], (), loc)
                            dtype_assign = ir.Assign(dtype_getattr, dtype_var, loc)
                            self.calltypes[dtype_getattr] = signature(
                                self.typemap[dtype_var.name], self.typemap[typ_var.name])

                            # The original A.dtype rhs is replaced with result of this call.
                            instr.value = dtype_var
                            # Add statements to body of the code.
                            block.body.insert(0, dtype_assign)
                            block.body.insert(0, dtype_attr_assign)
                            block.body.insert(0, typ_var_assign)
                            block.body.insert(0, g_np_assign)
                            self.stats['replaced_dtype'] += 1
                            break

def find_template(op):
    for ft in numba.core.typing.templates.builtin_registry.functions:
        if ft.key == op:
            return ft


class ParforPassStates:
    """This class encapsulates all internal states of the ParforPass.
    """

    def __init__(self, func_ir, typemap, calltypes, return_type, typingctx,
                 targetctx,  options, flags, metadata,
                 diagnostics=ParforDiagnostics()):
        self.func_ir = func_ir
        self.typemap = typemap
        self.calltypes = calltypes
        self.typingctx = typingctx
        self.targetctx = targetctx
        self.return_type = return_type
        self.options = options
        self.diagnostics = diagnostics
        self.swapped_fns = diagnostics.replaced_fns
        self.fusion_info = diagnostics.fusion_info
        self.nested_fusion_info = diagnostics.nested_fusion_info

        self.array_analysis = array_analysis.ArrayAnalysis(
            self.typingctx, self.func_ir, self.typemap, self.calltypes,
        )

        ir_utils._the_max_label.update(max(func_ir.blocks.keys()))
        self.flags = flags
        self.metadata = metadata
        if "parfors" not in metadata:
            metadata["parfors"] = {}


class ConvertInplaceBinop:
    """Parfor subpass to convert setitem on Arrays
    """
    def __init__(self, pass_states):
        """
        Parameters
        ----------
        pass_states : ParforPassStates
        """
        self.pass_states = pass_states
        self.rewritten = []

    def run(self, blocks):
        pass_states = self.pass_states
        # convert expressions like A += ... where A is an array.
        topo_order = find_topo_order(blocks)
        # variables available in the program so far (used for finding map
        # functions in array_expr lowering)
        for label in topo_order:
            block = blocks[label]
            new_body = []
            equiv_set = pass_states.array_analysis.get_equiv_set(label)
            for instr in block.body:
                if isinstance(instr, ir.Assign):
                    lhs = instr.target
                    expr = instr.value
                    if isinstance(expr, ir.Expr) and expr.op == 'inplace_binop':
                        loc = expr.loc
                        target = expr.lhs
                        value = expr.rhs
                        target_typ = pass_states.typemap[target.name]
                        value_typ = pass_states.typemap[value.name]
                        # Handle A op= ...
                        if isinstance(target_typ, types.npytypes.Array):
                            # RHS is an array
                            if isinstance(value_typ, types.npytypes.Array):
                                new_instr = self._inplace_binop_to_parfor(equiv_set,
                                        loc, expr.immutable_fn, target, value)
                                self.rewritten.append(
                                    dict(old=instr, new=new_instr,
                                         reason='inplace_binop'),
                                )
                                instr = [new_instr, ir.Assign(target, lhs, loc)]
                if isinstance(instr, list):
                    new_body.extend(instr)
                else:
                    new_body.append(instr)
            block.body = new_body

    def _inplace_binop_to_parfor(self, equiv_set, loc, op, target, value):
        """generate parfor from setitem node with a boolean or slice array indices.
        The value can be either a scalar or an array variable, and if a boolean index
        is used for the latter case, the same index must be used for the value too.
        """
        pass_states = self.pass_states
        scope = target.scope
        arr_typ = pass_states.typemap[target.name]
        el_typ = arr_typ.dtype
        init_block = ir.Block(scope, loc)
        value_typ = pass_states.typemap[value.name]

        size_vars = equiv_set.get_shape(target)

        # generate loopnests and size variables from target correlations
        index_vars, loopnests = _mk_parfor_loops(pass_states.typemap, size_vars, scope, loc)

        # generate body
        body_label = next_label()
        body_block = ir.Block(scope, loc)
        index_var, index_var_typ = _make_index_var(
                pass_states.typemap, scope, index_vars, body_block)

        # Read value.
        value_var = ir.Var(scope, mk_unique_var("$value_var"), loc)
        pass_states.typemap[value_var.name] = value_typ.dtype
        getitem_call = ir.Expr.getitem(value, index_var, loc)
        pass_states.calltypes[getitem_call] = signature(
            value_typ.dtype, value_typ, index_var_typ)
        body_block.body.append(ir.Assign(getitem_call, value_var, loc))

        # Read target
        target_var = ir.Var(scope, mk_unique_var("$target_var"), loc)
        pass_states.typemap[target_var.name] = el_typ
        getitem_call = ir.Expr.getitem(target, index_var, loc)
        pass_states.calltypes[getitem_call] = signature(
            el_typ, arr_typ, index_var_typ)
        body_block.body.append(ir.Assign(getitem_call, target_var, loc))

        # Create temp to hold result.
        expr_out_var = ir.Var(scope, mk_unique_var("$expr_out_var"), loc)
        pass_states.typemap[expr_out_var.name] = el_typ

        # Create binop and assign result to temporary.
        binop_expr = ir.Expr.binop(op, target_var, value_var, loc)
        body_block.body.append(ir.Assign(binop_expr, expr_out_var, loc))
        unified_type = self.pass_states.typingctx.unify_pairs(el_typ, value_typ.dtype)
        pass_states.calltypes[binop_expr] = signature(
            unified_type, unified_type, unified_type)

        # Write to target
        setitem_node = ir.SetItem(target, index_var, expr_out_var, loc)
        pass_states.calltypes[setitem_node] = signature(
            types.none, arr_typ, index_var_typ, el_typ)
        body_block.body.append(setitem_node)

        parfor = Parfor(loopnests, init_block, {}, loc, index_var, equiv_set,
                        ('inplace_binop', ''), pass_states.flags)
        parfor.loop_body = {body_label: body_block}
        if config.DEBUG_ARRAY_OPT >= 1:
            print("parfor from inplace_binop")
            parfor.dump()
        return parfor

    def _type_getitem(self, args):
        fnty = operator.getitem
        return self.pass_states.typingctx.resolve_function_type(fnty, tuple(args), {})


def get_index_var(x):
    return x.index if isinstance(x, ir.SetItem) else x.index_var


class ConvertSetItemPass:
    """Parfor subpass to convert setitem on Arrays
    """
    def __init__(self, pass_states):
        """
        Parameters
        ----------
        pass_states : ParforPassStates
        """
        self.pass_states = pass_states
        self.rewritten = []

    def run(self, blocks):
        pass_states = self.pass_states
        # convert setitem expressions like A[C] = c or A[C] = B[C] to parfor,
        # where C is a boolean array.
        topo_order = find_topo_order(blocks)
        # variables available in the program so far (used for finding map
        # functions in array_expr lowering)
        for label in topo_order:
            block = blocks[label]
            new_body = []
            equiv_set = pass_states.array_analysis.get_equiv_set(label)
            for instr in block.body:
                if isinstance(instr, (ir.StaticSetItem, ir.SetItem)):
                    loc = instr.loc
                    target = instr.target
                    index = get_index_var(instr)
                    value = instr.value
                    target_typ = pass_states.typemap[target.name]
                    index_typ = pass_states.typemap[index.name]
                    value_typ = pass_states.typemap[value.name]
                    # Handle A[boolean_array] = <scalar or array>
                    if isinstance(target_typ, types.npytypes.Array):
                        if (isinstance(index_typ, types.npytypes.Array) and
                            isinstance(index_typ.dtype, types.Boolean) and
                            target_typ.ndim == index_typ.ndim):
                            # RHS is a scalar number
                            if isinstance(value_typ, types.Number):
                                new_instr = self._setitem_to_parfor(equiv_set,
                                        loc, target, index, value)
                                self.rewritten.append(
                                    dict(old=instr, new=new_instr,
                                        reason='masked_assign_broadcast_scalar'),
                                )
                                instr = new_instr
                            # RHS is an array
                            elif isinstance(value_typ, types.npytypes.Array):
                                val_def = guard(get_definition, pass_states.func_ir,
                                                value.name)
                                if (isinstance(val_def, ir.Expr) and
                                    val_def.op == 'getitem' and
                                    val_def.index.name == index.name):
                                    new_instr = self._setitem_to_parfor(equiv_set,
                                            loc, target, index, val_def.value)
                                    self.rewritten.append(
                                        dict(old=instr, new=new_instr,
                                             reason='masked_assign_array'),
                                    )
                                    instr = new_instr
                        else:
                            # Handle A[:] = x
                            shape = equiv_set.get_shape(instr)
                            # Don't converted broadcasted setitems into parfors.
                            if isinstance(index_typ, types.BaseTuple):
                                # The sliced dims are those in the index that
                                # are made of slices.  Count the numbers of slices
                                # in the index tuple.
                                sliced_dims = len(list(filter(
                                    lambda x: isinstance(x, types.misc.SliceType),
                                    index_typ.types)))
                            elif isinstance(index_typ, types.misc.SliceType):
                                # For singular indices there can be a bare slice
                                # and if so there is one dimension being set.
                                sliced_dims = 1
                            else:
                                sliced_dims = 0

                            # Only create a parfor for this setitem if we know the
                            # shape of the output and number of dimensions set is
                            # equal to the number of dimensions on the right side.
                            if (shape is not None and
                                (not isinstance(value_typ, types.npytypes.Array) or
                                 sliced_dims == value_typ.ndim)):
                                new_instr = self._setitem_to_parfor(equiv_set,
                                        loc, target, index, value, shape=shape)
                                self.rewritten.append(
                                        dict(old=instr, new=new_instr,
                                             reason='slice'),
                                )
                                instr = new_instr
                new_body.append(instr)
            block.body = new_body

    def _setitem_to_parfor(self, equiv_set, loc, target, index, value, shape=None):
        """generate parfor from setitem node with a boolean or slice array indices.
        The value can be either a scalar or an array variable, and if a boolean index
        is used for the latter case, the same index must be used for the value too.
        """
        pass_states = self.pass_states
        scope = target.scope
        arr_typ = pass_states.typemap[target.name]
        el_typ = arr_typ.dtype
        index_typ = pass_states.typemap[index.name]
        init_block = ir.Block(scope, loc)

        if shape:
            # Slice index is being used on the target array, we'll have to create
            # a sub-array so that the target dimension matches the given shape.
            assert(isinstance(index_typ, types.BaseTuple) or
                   isinstance(index_typ, types.SliceType))
            # setitem has a custom target shape
            size_vars = shape
            # create a new target array via getitem
            subarr_var = ir.Var(scope, mk_unique_var("$subarr"), loc)
            getitem_call = ir.Expr.getitem(target, index, loc)
            subarr_typ = typing.arraydecl.get_array_index_type( arr_typ, index_typ).result
            pass_states.typemap[subarr_var.name] = subarr_typ
            pass_states.calltypes[getitem_call] = self._type_getitem((arr_typ, index_typ))
            init_block.append(ir.Assign(getitem_call, subarr_var, loc))
            target = subarr_var
        else:
            # Otherwise it is a boolean array that is used as index.
            assert(isinstance(index_typ, types.ArrayCompatible))
            size_vars = equiv_set.get_shape(target)
            bool_typ = index_typ.dtype

        # generate loopnests and size variables from lhs correlations
        loopnests = []
        index_vars = []
        for size_var in size_vars:
            index_var = ir.Var(scope, mk_unique_var("parfor_index"), loc)
            index_vars.append(index_var)
            pass_states.typemap[index_var.name] = types.uintp
            loopnests.append(LoopNest(index_var, 0, size_var, 1))

        # generate body
        body_label = next_label()
        body_block = ir.Block(scope, loc)
        index_var, index_var_typ = _make_index_var(
                pass_states.typemap, scope, index_vars, body_block)
        parfor = Parfor(loopnests, init_block, {}, loc, index_var, equiv_set,
                        ('setitem', ''), pass_states.flags)
        if shape:
            # slice subarray
            parfor.loop_body = {body_label: body_block}
            true_block = body_block
            end_label = None
        else:
            # boolean mask
            true_label = next_label()
            true_block = ir.Block(scope, loc)
            end_label = next_label()
            end_block = ir.Block(scope, loc)
            parfor.loop_body = {body_label: body_block,
                                true_label: true_block,
                                end_label:  end_block,
                                }
            mask_var = ir.Var(scope, mk_unique_var("$mask_var"), loc)
            pass_states.typemap[mask_var.name] = bool_typ
            mask_val = ir.Expr.getitem(index, index_var, loc)
            body_block.body.extend([
               ir.Assign(mask_val, mask_var, loc),
               ir.Branch(mask_var, true_label, end_label, loc)
            ])

        value_typ = pass_states.typemap[value.name]
        if isinstance(value_typ, types.npytypes.Array):
            value_var = ir.Var(scope, mk_unique_var("$value_var"), loc)
            pass_states.typemap[value_var.name] = value_typ.dtype
            getitem_call = ir.Expr.getitem(value, index_var, loc)
            pass_states.calltypes[getitem_call] = signature(
                value_typ.dtype, value_typ, index_var_typ)
            true_block.body.append(ir.Assign(getitem_call, value_var, loc))
        else:
            value_var = value
        setitem_node = ir.SetItem(target, index_var, value_var, loc)
        pass_states.calltypes[setitem_node] = signature(
            types.none, pass_states.typemap[target.name], index_var_typ, el_typ)
        true_block.body.append(setitem_node)
        if end_label:
            true_block.body.append(ir.Jump(end_label, loc))

        if config.DEBUG_ARRAY_OPT >= 1:
            print("parfor from setitem")
            parfor.dump()
        return parfor

    def _type_getitem(self, args):
        fnty = operator.getitem
        return self.pass_states.typingctx.resolve_function_type(fnty, tuple(args), {})


def _make_index_var(typemap, scope, index_vars, body_block, force_tuple=False):
    """ When generating a SetItem call to an array in a parfor, the general
    strategy is to generate a tuple if the array is more than 1 dimension.
    If it is 1 dimensional then you can use a simple variable.  This routine
    is also used when converting pndindex to parfor but pndindex requires a
    tuple even if the iteration space is 1 dimensional.  The pndindex use of
    this function will use force_tuple to make the output index a tuple even
    if it is one dimensional.
    """
    ndims = len(index_vars)
    loc = body_block.loc
    if ndims > 1 or force_tuple:
        tuple_var = ir.Var(scope, mk_unique_var(
            "$parfor_index_tuple_var"), loc)
        typemap[tuple_var.name] = types.containers.UniTuple(
            types.uintp, ndims)
        tuple_call = ir.Expr.build_tuple(list(index_vars), loc)
        tuple_assign = ir.Assign(tuple_call, tuple_var, loc)
        body_block.body.append(tuple_assign)
        return tuple_var, types.containers.UniTuple(types.uintp, ndims)
    elif ndims == 1:
        return index_vars[0], types.uintp
    else:
        raise errors.UnsupportedRewriteError(
            "Parfor does not handle arrays of dimension 0",
            loc=loc,
        )


def _mk_parfor_loops(typemap, size_vars, scope, loc):
    """
    Create loop index variables and build LoopNest objects for a parfor.
    """
    loopnests = []
    index_vars = []
    for size_var in size_vars:
        index_var = ir.Var(scope, mk_unique_var("parfor_index"), loc)
        index_vars.append(index_var)
        typemap[index_var.name] = types.uintp
        loopnests.append(LoopNest(index_var, 0, size_var, 1))
    return index_vars, loopnests

class ConvertNumpyPass:
    """
    Convert supported Numpy functions, as well as arrayexpr nodes, to
    parfor nodes.
    """
    def __init__(self, pass_states):
        self.pass_states = pass_states
        self.rewritten = []

    def run(self, blocks):
        pass_states = self.pass_states
        topo_order = find_topo_order(blocks)
        # variables available in the program so far (used for finding map
        # functions in array_expr lowering)
        avail_vars = []
        for label in topo_order:
            block = blocks[label]
            new_body = []
            equiv_set = pass_states.array_analysis.get_equiv_set(label)
            for instr in block.body:
                if isinstance(instr, ir.Assign):
                    expr = instr.value
                    lhs = instr.target
                    lhs_typ = self.pass_states.typemap[lhs.name]
                    if self._is_C_or_F_order(lhs_typ):
                        if guard(self._is_supported_npycall, expr):
                            new_instr = self._numpy_to_parfor(equiv_set, lhs, expr)
                            if new_instr is not None:
                                self.rewritten.append(dict(
                                    old=instr,
                                    new=new_instr,
                                    reason='numpy_allocator',
                                ))
                                instr = new_instr
                        elif isinstance(expr, ir.Expr) and expr.op == 'arrayexpr':
                            new_instr = self._arrayexpr_to_parfor(
                                equiv_set, lhs, expr, avail_vars)
                            self.rewritten.append(dict(
                                old=instr,
                                new=new_instr,
                                reason='arrayexpr',
                            ))
                            instr = new_instr
                    avail_vars.append(lhs.name)
                new_body.append(instr)
            block.body = new_body

    def _is_C_order(self, arr_name):
        if isinstance(arr_name, types.npytypes.Array):
            return arr_name.layout == 'C' and arr_name.ndim > 0
        elif arr_name is str:
            typ = self.pass_states.typemap[arr_name]
            return (isinstance(typ, types.npytypes.Array) and
                    typ.layout == 'C' and
                    typ.ndim > 0)
        else:
            return False

    def _is_C_or_F_order(self, arr_name):
        if isinstance(arr_name, types.npytypes.Array):
            return (arr_name.layout == 'C' or arr_name.layout == 'F') and arr_name.ndim > 0
        elif arr_name is str:
            typ = self.pass_states.typemap[arr_name]
            return (isinstance(typ, types.npytypes.Array) and
                    (typ.layout == 'C' or typ.layout == 'F') and
                    typ.ndim > 0)
        else:
            return False

    def _arrayexpr_to_parfor(self, equiv_set, lhs, arrayexpr, avail_vars):
        """generate parfor from arrayexpr node, which is essentially a
        map with recursive tree.
        """
        pass_states = self.pass_states
        scope = lhs.scope
        loc = lhs.loc
        expr = arrayexpr.expr
        arr_typ = pass_states.typemap[lhs.name]
        el_typ = arr_typ.dtype

        # generate loopnests and size variables from lhs correlations
        size_vars = equiv_set.get_shape(lhs)
        index_vars, loopnests = _mk_parfor_loops(pass_states.typemap, size_vars, scope, loc)

        # generate init block and body
        init_block = ir.Block(scope, loc)
        init_block.body = mk_alloc(
            pass_states.typingctx,
            pass_states.typemap, pass_states.calltypes, lhs,
            tuple(size_vars), el_typ, scope, loc,
            pass_states.typemap[lhs.name])
        body_label = next_label()
        body_block = ir.Block(scope, loc)
        expr_out_var = ir.Var(scope, mk_unique_var("$expr_out_var"), loc)
        pass_states.typemap[expr_out_var.name] = el_typ

        index_var, index_var_typ = _make_index_var(
            pass_states.typemap, scope, index_vars, body_block)

        body_block.body.extend(
            _arrayexpr_tree_to_ir(
                pass_states.func_ir,
                pass_states.typingctx,
                pass_states.typemap,
                pass_states.calltypes,
                equiv_set,
                init_block,
                expr_out_var,
                expr,
                index_var,
                index_vars,
                avail_vars))

        pat = ('array expression {}'.format(repr_arrayexpr(arrayexpr.expr)),)

        parfor = Parfor(loopnests, init_block, {}, loc, index_var, equiv_set, pat[0], pass_states.flags)

        setitem_node = ir.SetItem(lhs, index_var, expr_out_var, loc)
        pass_states.calltypes[setitem_node] = signature(
            types.none, pass_states.typemap[lhs.name], index_var_typ, el_typ)
        body_block.body.append(setitem_node)
        parfor.loop_body = {body_label: body_block}
        if config.DEBUG_ARRAY_OPT >= 1:
            print("parfor from arrayexpr")
            parfor.dump()
        return parfor

    def _is_supported_npycall(self, expr):
        """check if we support parfor translation for
        this Numpy call.
        """
        call_name, mod_name = find_callname(self.pass_states.func_ir, expr)
        if not (isinstance(mod_name, str) and mod_name.startswith('numpy')):
            return False
        if call_name in ['zeros', 'ones']:
            return True
        if mod_name == 'numpy.random' and call_name in random_calls:
            return True
        # TODO: add more calls
        return False

    def _numpy_to_parfor(self, equiv_set, lhs, expr):
        call_name, mod_name = find_callname(self.pass_states.func_ir, expr)
        args = expr.args
        kws = dict(expr.kws)
        if call_name in ['zeros', 'ones'] or mod_name == 'numpy.random':
            return self._numpy_map_to_parfor(equiv_set, call_name, lhs, args, kws, expr)
        # return error if we couldn't handle it (avoid rewrite infinite loop)
        raise errors.UnsupportedRewriteError(
            f"parfor translation failed for {expr}", loc=expr.loc,
        )

    def _numpy_map_to_parfor(self, equiv_set, call_name, lhs, args, kws, expr):
        """generate parfor from Numpy calls that are maps.
        """
        pass_states = self.pass_states
        scope = lhs.scope
        loc = lhs.loc
        arr_typ = pass_states.typemap[lhs.name]
        el_typ = arr_typ.dtype

        # generate loopnests and size variables from lhs correlations
        size_vars = equiv_set.get_shape(lhs)
        if size_vars is None:
            if config.DEBUG_ARRAY_OPT >= 1:
                print("Could not convert numpy map to parfor, unknown size")
            return None

        index_vars, loopnests = _mk_parfor_loops(pass_states.typemap, size_vars, scope, loc)

        # generate init block and body
        init_block = ir.Block(scope, loc)
        init_block.body = mk_alloc(
            pass_states.typingctx,
            pass_states.typemap, pass_states.calltypes, lhs,
            tuple(size_vars), el_typ, scope, loc,
            pass_states.typemap[lhs.name])
        body_label = next_label()
        body_block = ir.Block(scope, loc)
        expr_out_var = ir.Var(scope, mk_unique_var("$expr_out_var"), loc)
        pass_states.typemap[expr_out_var.name] = el_typ

        index_var, index_var_typ = _make_index_var(
            pass_states.typemap, scope, index_vars, body_block)

        if call_name == 'zeros':
            value = ir.Const(el_typ(0), loc)
        elif call_name == 'ones':
            value = ir.Const(el_typ(1), loc)
        elif call_name in random_calls:
            # remove size arg to reuse the call expr for single value
            _remove_size_arg(call_name, expr)
            # update expr type
            new_arg_typs, new_kw_types = _get_call_arg_types(
                expr, pass_states.typemap)
            pass_states.calltypes.pop(expr)
            pass_states.calltypes[expr] = pass_states.typemap[expr.func.name].get_call_type(
                typing.Context(), new_arg_typs, new_kw_types)
            value = expr
        else:
            raise NotImplementedError(
                "Map of numpy.{} to parfor is not implemented".format(call_name))

        value_assign = ir.Assign(value, expr_out_var, loc)
        body_block.body.append(value_assign)

        setitem_node = ir.SetItem(lhs, index_var, expr_out_var, loc)
        pass_states.calltypes[setitem_node] = signature(
            types.none, pass_states.typemap[lhs.name], index_var_typ, el_typ)
        body_block.body.append(setitem_node)

        parfor = Parfor(loopnests, init_block, {}, loc, index_var, equiv_set,
                        ('{} function'.format(call_name,), 'NumPy mapping'),
                        pass_states.flags)
        parfor.loop_body = {body_label: body_block}
        if config.DEBUG_ARRAY_OPT >= 1:
            print("generated parfor for numpy map:")
            parfor.dump()
        return parfor


class ConvertReducePass:
    """
    Find reduce() calls and convert them to parfors.
    """
    def __init__(self, pass_states):
        self.pass_states = pass_states
        self.rewritten = []

    def run(self, blocks):
        pass_states = self.pass_states

        topo_order = find_topo_order(blocks)
        for label in topo_order:
            block = blocks[label]
            new_body = []
            equiv_set = pass_states.array_analysis.get_equiv_set(label)
            for instr in block.body:
                parfor = None
                if isinstance(instr, ir.Assign):
                    loc = instr.loc
                    lhs = instr.target
                    expr = instr.value
                    callname = guard(find_callname, pass_states.func_ir, expr)
                    if (callname == ('reduce', 'builtins')
                        or callname == ('reduce', '_functools')):
                        # reduce function with generic function
                        parfor = guard(self._reduce_to_parfor, equiv_set, lhs,
                                       expr.args, loc)
                    if parfor:
                        self.rewritten.append(dict(
                            new=parfor,
                            old=instr,
                            reason='reduce',
                        ))
                        instr = parfor
                new_body.append(instr)
            block.body = new_body
        return

    def _reduce_to_parfor(self, equiv_set, lhs, args, loc):
        """
        Convert a reduce call to a parfor.
        The call arguments should be (call_name, array, init_value).
        """
        pass_states = self.pass_states

        scope = lhs.scope
        call_name = args[0]
        in_arr = args[1]
        arr_def = get_definition(pass_states.func_ir, in_arr.name)

        mask_var = None
        mask_indices = None

        # Search for array[boolean_mask]
        mask_query_result = guard(_find_mask, pass_states.typemap, pass_states.func_ir, arr_def)
        if mask_query_result:
            in_arr, mask_var, mask_typ, mask_indices = mask_query_result

        init_val = args[2]
        size_vars = equiv_set.get_shape(in_arr if mask_indices is None else mask_var)
        if size_vars is None:
            return None

        index_vars, loopnests = _mk_parfor_loops(pass_states.typemap, size_vars, scope, loc)
        mask_index = index_vars
        if mask_indices:
            # the following is never tested
            raise AssertionError("unreachable")
            index_vars = tuple(x if x else index_vars[0] for x in mask_indices)
        acc_var = lhs

        # init block has to init the reduction variable
        init_block = ir.Block(scope, loc)
        init_block.body.append(ir.Assign(init_val, acc_var, loc))

        # produce loop body
        body_label = next_label()
        index_var, loop_body = self._mk_reduction_body(call_name,
                                scope, loc, index_vars, in_arr, acc_var)
        if mask_indices:
            # the following is never tested
            raise AssertionError("unreachable")
            index_var = mask_index[0]

        if mask_var is not None:
            true_label = min(loop_body.keys())
            false_label = max(loop_body.keys())
            body_block = ir.Block(scope, loc)
            loop_body[body_label] = body_block
            mask = ir.Var(scope, mk_unique_var("$mask_val"), loc)
            pass_states.typemap[mask.name] = mask_typ
            mask_val = ir.Expr.getitem(mask_var, index_var, loc)
            body_block.body.extend([
               ir.Assign(mask_val, mask, loc),
               ir.Branch(mask, true_label, false_label, loc)
            ])

        parfor = Parfor(loopnests, init_block, loop_body, loc, index_var,
                        equiv_set, ('{} function'.format(call_name),
                                    'reduction'), pass_states.flags)
        if config.DEBUG_ARRAY_OPT >= 1:
            print("parfor from reduction")
            parfor.dump()
        return parfor

    def _mk_reduction_body(self, call_name, scope, loc,
                           index_vars, in_arr, acc_var):
        """
        Produce the body blocks for a reduction function indicated by call_name.
        """
        from numba.core.inline_closurecall import check_reduce_func

        pass_states = self.pass_states
        reduce_func = get_definition(pass_states.func_ir, call_name)
        fcode = check_reduce_func(pass_states.func_ir, reduce_func)

        arr_typ = pass_states.typemap[in_arr.name]
        in_typ = arr_typ.dtype
        body_block = ir.Block(scope, loc)
        index_var, index_var_type = _make_index_var(
            pass_states.typemap, scope, index_vars, body_block)

        tmp_var = ir.Var(scope, mk_unique_var("$val"), loc)
        pass_states.typemap[tmp_var.name] = in_typ
        getitem_call = ir.Expr.getitem(in_arr, index_var, loc)
        pass_states.calltypes[getitem_call] = signature(
            in_typ, arr_typ, index_var_type)
        body_block.append(ir.Assign(getitem_call, tmp_var, loc))

        reduce_f_ir = compile_to_numba_ir(fcode,
                                        pass_states.func_ir.func_id.func.__globals__,
                                        pass_states.typingctx,
                                        pass_states.targetctx,
                                        (in_typ, in_typ),
                                        pass_states.typemap,
                                        pass_states.calltypes)
        loop_body = reduce_f_ir.blocks
        end_label = next_label()
        end_block = ir.Block(scope, loc)
        loop_body[end_label] = end_block
        first_reduce_label = min(reduce_f_ir.blocks.keys())
        first_reduce_block = reduce_f_ir.blocks[first_reduce_label]
        body_block.body.extend(first_reduce_block.body)
        first_reduce_block.body = body_block.body
        replace_arg_nodes(first_reduce_block, [acc_var, tmp_var])
        replace_returns(loop_body, acc_var, end_label)
        return index_var, loop_body


class ConvertLoopPass:
    """Build Parfor nodes from prange loops.
    """
    def __init__(self, pass_states):
        self.pass_states = pass_states
        self.rewritten = []

    def run(self, blocks):
        pass_states = self.pass_states

        call_table, _ = get_call_table(blocks)
        cfg = compute_cfg_from_blocks(blocks)
        usedefs = compute_use_defs(blocks)
        live_map = compute_live_map(cfg, blocks, usedefs.usemap, usedefs.defmap)
        loops = cfg.loops()
        sized_loops = [(loops[k], len(loops[k].body)) for k in loops.keys()]
        moved_blocks = []
        # We go over all loops, smaller loops first (inner first)
        for loop, s in sorted(sized_loops, key=lambda tup: tup[1]):
            if len(loop.entries) != 1 or len(loop.exits) != 1:
                if not config.DISABLE_PERFORMANCE_WARNINGS:
                    for entry in loop.entries:
                        for inst in blocks[entry].body:
                            # if prange or pndindex call
                            if (
                                isinstance(inst, ir.Assign)
                                and isinstance(inst.value, ir.Expr)
                                and inst.value.op == "call"
                                and self._is_parallel_loop(
                                    inst.value.func.name, call_table)
                            ):
                                msg = "\nprange or pndindex loop " \
                                      "will not be executed in " \
                                      "parallel due to there being more than one " \
                                      "entry to or exit from the loop (e.g., an " \
                                      "assertion)."
                                warnings.warn(
                                    errors.NumbaPerformanceWarning(
                                        msg, inst.loc))
                continue

            entry = list(loop.entries)[0]
            for inst in blocks[entry].body:
                # if prange or pndindex call
                if (isinstance(inst, ir.Assign)
                        and isinstance(inst.value, ir.Expr)
                        and inst.value.op == 'call'
                        and self._is_parallel_loop(inst.value.func.name, call_table)):
                    # Here we've found a parallel loop, either prange or pndindex.
                    # We create a parfor from this loop and then overwrite the contents
                    # of the original loop header block to contain this parfor and then
                    # a jump to the original loop exit block.  Other blocks in the
                    # original loop are discarded.
                    body_labels = [ l for l in loop.body if
                                    l in blocks and l != loop.header ]
                    args = inst.value.args
                    loop_kind, loop_replacing = self._get_loop_kind(inst.value.func.name,
                                                                    call_table)
                    # Get the body of the header of the loops minus the branch terminator
                    # The general approach is to prepend the header block to the first
                    # body block and then let dead code removal handle removing unneeded
                    # statements.  Not all statements in the header block are unnecessary.
                    header_body = blocks[loop.header].body[:-1]
                    # find loop index variable (pair_first in header block)
                    loop_index = None
                    for hbi, stmt in enumerate(header_body):
                        if (isinstance(stmt, ir.Assign)
                                and isinstance(stmt.value, ir.Expr)
                                and stmt.value.op == 'pair_first'):
                            loop_index = stmt.target.name
                            li_index = hbi
                            break
                    assert(loop_index is not None)
                    # Remove pair_first from header.
                    # We have to remove the pair_first by hand since it causes problems
                    # for some code below if we don't.
                    header_body = header_body[:li_index] + header_body[li_index+1:]

                    # loop_index may be assigned to other vars
                    # get header copies to find all of them
                    cps, _ = get_block_copies({0: blocks[loop.header]},
                                              pass_states.typemap)
                    cps = cps[0]
                    loop_index_vars = set(t for t, v in cps if v == loop_index)
                    loop_index_vars.add(loop_index)

                    scope = blocks[entry].scope
                    loc = inst.loc
                    equiv_set = pass_states.array_analysis.get_equiv_set(loop.header)
                    init_block = ir.Block(scope, loc)
                    init_block.body = self._get_prange_init_block(blocks[entry],
                                                            call_table, args)
                    loop_body = {l: blocks[l] for l in body_labels}
                    # Add an empty block to the end of loop body
                    end_label = next_label()
                    loop_body[end_label] = ir.Block(scope, loc)

                    # Detect races in the prange.
                    # Races are defs in the parfor body that are live at the exit block.
                    bodydefs = set()
                    for bl in body_labels:
                        bodydefs = bodydefs.union(usedefs.defmap[bl])
                    exit_lives = set()
                    for bl in loop.exits:
                        exit_lives = exit_lives.union(live_map[bl])
                    races = bodydefs.intersection(exit_lives)
                    # It is possible for the result of an ir.Global to be flagged
                    # as a race if it is defined in this Parfor and then used in
                    # a subsequent Parfor.  push_call_vars() in the Parfor pass
                    # copies such ir.Global nodes into the Parfors in which they
                    # are used so no need to treat things of type Module as a race.
                    races = races.intersection({x for x in races
                              if not isinstance(pass_states.typemap[x], types.misc.Module)})

                    # replace jumps to header block with the end block
                    for l in body_labels:
                        last_inst = loop_body[l].body[-1]
                        if (isinstance(last_inst, ir.Jump) and
                            last_inst.target == loop.header):
                            last_inst.target = end_label

                    def find_indexed_arrays():
                        """find expressions that involve getitem using the
                        index variable. Return both the arrays and expressions.
                        """
                        indices = copy.copy(loop_index_vars)
                        for block in loop_body.values():
                            for inst in block.find_insts(ir.Assign):
                                if (isinstance(inst.value, ir.Var) and
                                    inst.value.name in indices):
                                    indices.add(inst.target.name)
                        arrs = []
                        exprs = []
                        for block in loop_body.values():
                            for inst in block.body:
                                lv = set(x.name for x in inst.list_vars())
                                if lv & indices:
                                    if lv.issubset(indices):
                                        continue
                                    require(isinstance(inst, ir.Assign))
                                    expr = inst.value
                                    require(isinstance(expr, ir.Expr) and
                                       expr.op in ['getitem', 'static_getitem'])
                                    arrs.append(expr.value.name)
                                    exprs.append(expr)
                        return arrs, exprs

                    mask_var = None
                    mask_indices = None
                    def find_mask_from_size(size_var):
                        """Find the case where size_var is defined by A[M].shape,
                        where M is a boolean array.
                        """
                        size_def = get_definition(pass_states.func_ir, size_var)
                        require(size_def and isinstance(size_def, ir.Expr) and
                                size_def.op == 'getattr' and size_def.attr == 'shape')
                        arr_var = size_def.value
                        live_vars = set.union(*[live_map[l] for l in loop.exits])
                        index_arrs, index_exprs = find_indexed_arrays()
                        require([arr_var.name] == list(index_arrs))
                        # input array has to be dead after loop
                        require(arr_var.name not in live_vars)
                        # loop for arr's definition, where size = arr.shape
                        arr_def = get_definition(pass_states.func_ir, size_def.value)
                        result = _find_mask(pass_states.typemap, pass_states.func_ir, arr_def)

                        # The following is never tested.
                        raise AssertionError("unreachable")
                        # Found the mask.
                        # Replace B[i] with A[i], where B = A[M]
                        for expr in index_exprs:
                            expr.value = result[0]
                        return result

                    # pndindex and prange are provably positive except when
                    # user provides negative start to prange()
                    unsigned_index = True
                    # TODO: support array mask optimization for prange
                    # TODO: refactor and simplify array mask optimization
                    if loop_kind == 'pndindex':
                        assert(equiv_set.has_shape(args[0]))
                        # see if input array to pndindex is output of array
                        # mask like B = A[M]
                        result = guard(find_mask_from_size, args[0])
                        if result:
                            in_arr, mask_var, mask_typ, mask_indices = result
                        else:
                            in_arr = args[0]
                        assert(isinstance(in_arr, ir.Var))
                        in_arr_typ = pass_states.typemap[in_arr.name]
                        if isinstance(in_arr_typ, types.Integer):
                            index_var = ir.Var(scope, mk_unique_var("parfor_index"), loc)
                            pass_states.typemap[index_var.name] = types.uintp
                            loops = [LoopNest(index_var, 0, in_arr, 1)]
                            index_vars = [index_var]
                        else:
                            size_vars = equiv_set.get_shape(in_arr
                                          if mask_indices is None else mask_var)
                            index_vars, loops = _mk_parfor_loops(
                                pass_states.typemap, size_vars, scope, loc,
                            )
                        assert(len(loops) > 0)
                        orig_index = index_vars
                        if mask_indices:
                            # replace mask indices if required;
                            # integer indices of original array should be used
                            # instead of parfor indices
                            index_vars = tuple(x if x else index_vars[0]
                                               for x in mask_indices)
                        first_body_block = loop_body[min(loop_body.keys())]
                        body_block = ir.Block(scope, loc)
                        index_var, index_var_typ = _make_index_var(
                            pass_states.typemap, scope, index_vars, body_block,
                            force_tuple=True
                        )
                        body = body_block.body + first_body_block.body
                        first_body_block.body = body
                        if mask_indices:
                            orig_index_var = orig_index[0]
                        else:
                            orig_index_var = index_var

                        # if masked array optimization is being applied, create
                        # the branch for array selection
                        if mask_var is not None:
                            # The following code are not tested
                            raise AssertionError("unreachable")
                            body_label = next_label()
                            # loop_body needs new labels greater than body_label
                            loop_body = add_offset_to_labels(loop_body,
                                            body_label - min(loop_body.keys()) + 1)
                            labels = loop_body.keys()
                            true_label = min(labels)
                            false_label = max(labels)
                            body_block = ir.Block(scope, loc)
                            loop_body[body_label] = body_block
                            mask = ir.Var(scope, mk_unique_var("$mask_val"), loc)
                            pass_states.typemap[mask.name] = mask_typ
                            mask_val = ir.Expr.getitem(mask_var, orig_index_var, loc)
                            body_block.body.extend([
                               ir.Assign(mask_val, mask, loc),
                               ir.Branch(mask, true_label, false_label, loc)
                            ])
                    else: # prange
                        start = 0
                        step = 1
                        size_var = args[0]
                        if len(args) == 2:
                            start = args[0]
                            size_var = args[1]
                        if len(args) == 3:
                            start = args[0]
                            size_var = args[1]
                            try:
                                step = pass_states.func_ir.get_definition(args[2])
                            except KeyError:
                                raise errors.UnsupportedRewriteError(
                                    "Only known step size is supported for prange",
                                    loc=inst.loc,
                                )
                            if not isinstance(step, ir.Const):
                                raise errors.UnsupportedRewriteError(
                                    "Only constant step size is supported for prange",
                                    loc=inst.loc,
                                )
                            step = step.value
                            if step != 1:
                                raise errors.UnsupportedRewriteError(
                                    "Only constant step size of 1 is supported for prange",
                                    loc=inst.loc,
                                )
                        index_var = ir.Var(scope, mk_unique_var("parfor_index"), loc)
                        # assume user-provided start to prange can be negative
                        # this is the only case parfor can have negative index
                        if isinstance(start, int) and start >= 0:
                            index_var_typ = types.uintp
                        else:
                            index_var_typ = types.intp
                            unsigned_index = False
                        loops = [LoopNest(index_var, start, size_var, step)]
                        pass_states.typemap[index_var.name] = index_var_typ

                        # We can't just drop the header block since there can be things
                        # in there other than the prange looping infrastructure.
                        # So we just add the header to the first loop body block (minus the
                        # branch) and let dead code elimination remove the unnecessary parts.
                        first_body_label = min(loop_body.keys())
                        loop_body[first_body_label].body = header_body + loop_body[first_body_label].body

                    index_var_map = {v: index_var for v in loop_index_vars}
                    replace_vars(loop_body, index_var_map)
                    if unsigned_index:
                        # need to replace signed array access indices to enable
                        # optimizations (see #2846)
                        self._replace_loop_access_indices(
                            loop_body, loop_index_vars, index_var)
                    parfor = Parfor(loops, init_block, loop_body, loc,
                                    orig_index_var if mask_indices else index_var,
                                    equiv_set,
                                    ("prange", loop_kind, loop_replacing),
                                    pass_states.flags, races=races)

                    blocks[loop.header].body = [parfor]
                    # We have to insert the header_body after the parfor because in
                    # a Numba loop this will be executed one more times before the
                    # branch and may contain instructions such as variable renamings
                    # that are relied upon later.
                    blocks[loop.header].body.extend(header_body)
                    blocks[loop.header].body.append(ir.Jump(list(loop.exits)[0], loc))
                    self.rewritten.append(dict(
                        old_loop=loop,
                        new=parfor,
                        reason='loop',
                    ))
                    # remove loop blocks from top level dict
                    for l in body_labels:
                        if l != loop.header:
                            blocks.pop(l)
                    if config.DEBUG_ARRAY_OPT >= 1:
                        print("parfor from loop")
                        parfor.dump()

    def _is_parallel_loop(self, func_var, call_table):
        # prange can be either getattr (numba.prange) or global (prange)
        if func_var not in call_table:
            return False
        call = call_table[func_var]
        return len(call) > 0 and (call[0] == 'prange' or call[0] == prange
                or call[0] == 'internal_prange' or call[0] == internal_prange
                or call[0] == 'pndindex' or call[0] == pndindex)

    def _get_loop_kind(self, func_var, call_table):
        """see if prange is user prange or internal"""
        pass_states = self.pass_states
        # prange can be either getattr (numba.prange) or global (prange)
        assert func_var in call_table
        call = call_table[func_var]
        assert len(call) > 0
        kind = 'user', ''
        if call[0] == 'internal_prange' or call[0] == internal_prange:
            try:
                kind = 'internal', (pass_states.swapped_fns[func_var][0], pass_states.swapped_fns[func_var][-1])
            except KeyError:
                # FIXME: Fix this issue... the code didn't manage to trace the
                # swapout for func_var so set the kind as internal so that the
                # transform can occur, it's just not tracked
                kind = 'internal', ('', '')
        elif call[0] == 'pndindex' or call[0] == pndindex:
            kind = 'pndindex', ''
        return kind

    def _get_prange_init_block(self, entry_block, call_table, prange_args):
        """
        If there is init_prange, find the code between init_prange and prange
        calls. Remove the code from entry_block and return it.
        """
        init_call_ind = -1
        prange_call_ind = -1
        init_body = []
        for i, inst in enumerate(entry_block.body):
            # if init_prange call
            if (isinstance(inst, ir.Assign) and isinstance(inst.value, ir.Expr)
                    and inst.value.op == 'call'
                    and self._is_prange_init(inst.value.func.name, call_table)):
                init_call_ind = i
            if (isinstance(inst, ir.Assign) and isinstance(inst.value, ir.Expr)
                    and inst.value.op == 'call'
                    and self._is_parallel_loop(inst.value.func.name, call_table)):
                prange_call_ind = i
        if init_call_ind != -1 and prange_call_ind != -1:
            # we save instructions that are used to calculate prange call args
            # in the entry block. The rest go to parfor init_block
            arg_related_vars = {v.name for v in prange_args}
            saved_nodes = []
            for i in reversed(range(init_call_ind+1, prange_call_ind)):
                inst = entry_block.body[i]
                inst_vars = {v.name for v in inst.list_vars()}
                if arg_related_vars & inst_vars:
                    arg_related_vars |= inst_vars
                    saved_nodes.append(inst)
                else:
                    init_body.append(inst)

            init_body.reverse()
            saved_nodes.reverse()
            entry_block.body = (entry_block.body[:init_call_ind]
                        + saved_nodes + entry_block.body[prange_call_ind+1:])

        return init_body

    def _is_prange_init(self, func_var, call_table):
        if func_var not in call_table:
            return False
        call = call_table[func_var]
        return len(call) > 0 and (call[0] == 'init_prange' or call[0] == init_prange)

    def _replace_loop_access_indices(self, loop_body, index_set, new_index):
        """
        Replace array access indices in a loop body with a new index.
        index_set has all the variables that are equivalent to loop index.
        """
        # treat new index like others since replacing it with itself is ok
        index_set.add(new_index.name)

        with dummy_return_in_loop_body(loop_body):
            labels = find_topo_order(loop_body)

        first_label = labels[0]
        added_indices = set()

        # traverse loop body and replace indices in getitem/setitem with
        # new_index if possible.
        # also, find equivalent indices defined in first block.
        for l in labels:
            block = loop_body[l]
            for stmt in block.body:
                if (isinstance(stmt, ir.Assign)
                        and isinstance(stmt.value, ir.Var)):
                    # the first block dominates others so we can use copies
                    # of indices safely
                    if (l == first_label and stmt.value.name in index_set
                            and stmt.target.name not in index_set):
                        index_set.add(stmt.target.name)
                        added_indices.add(stmt.target.name)
                    # make sure parallel index is not overwritten
                    else:
                        scope = block.scope

                        def unver(name):
                            from numba.core import errors
                            try:
                                return scope.get_exact(name).unversioned_name
                            except errors.NotDefinedError:
                                return name

                        if unver(stmt.target.name) in map(unver, index_set) and unver(stmt.target.name) != unver(stmt.value.name):
                            raise errors.UnsupportedRewriteError(
                                "Overwrite of parallel loop index",
                                loc=stmt.target.loc,
                            )



                if is_get_setitem(stmt):
                    index = index_var_of_get_setitem(stmt)
                    # statics can have none indices
                    if index is None:
                        continue
                    ind_def = guard(get_definition, self.pass_states.func_ir,
                                    index, lhs_only=True)
                    if (index.name in index_set
                            or (ind_def is not None
                                and ind_def.name in index_set)):
                        set_index_var_of_get_setitem(stmt, new_index)
                    # corner case where one dimension of a multi-dim access
                    # should be replaced
                    guard(self._replace_multi_dim_ind, ind_def, index_set,
                                                                     new_index)

                if isinstance(stmt, Parfor):
                    self._replace_loop_access_indices(stmt.loop_body, index_set, new_index)

        # remove added indices for correct recursive parfor handling
        index_set -= added_indices
        return

    def _replace_multi_dim_ind(self, ind_var, index_set, new_index):
        """
        replace individual indices in multi-dimensional access variable, which
        is a build_tuple
        """
        pass_states = self.pass_states
        require(ind_var is not None)
        # check for Tuple instead of UniTuple since some dims could be slices
        require(isinstance(pass_states.typemap[ind_var.name],
                (types.Tuple, types.UniTuple)))
        ind_def_node = get_definition(pass_states.func_ir, ind_var)
        require(isinstance(ind_def_node, ir.Expr)
                and ind_def_node.op == 'build_tuple')
        ind_def_node.items = [new_index if v.name in index_set else v
                              for v in ind_def_node.items]


def _find_mask(typemap, func_ir, arr_def):
    """check if an array is of B[...M...], where M is a
    boolean array, and other indices (if available) are ints.
    If found, return B, M, M's type, and a tuple representing mask indices.
    Otherwise, raise GuardException.
    """
    require(isinstance(arr_def, ir.Expr) and arr_def.op == 'getitem')
    value = arr_def.value
    index = arr_def.index
    value_typ = typemap[value.name]
    index_typ = typemap[index.name]
    ndim = value_typ.ndim
    require(isinstance(value_typ, types.npytypes.Array))
    if (isinstance(index_typ, types.npytypes.Array) and
        isinstance(index_typ.dtype, types.Boolean) and
        ndim == index_typ.ndim):
        return value, index, index_typ.dtype, None
    elif isinstance(index_typ, types.BaseTuple):
        # Handle multi-dimension differently by requiring
        # all indices to be constant except the one for mask.
        seq, op = find_build_sequence(func_ir, index)
        require(op == 'build_tuple' and len(seq) == ndim)
        count_consts = 0
        mask_indices = []
        mask_var = None
        for ind in seq:
            index_typ = typemap[ind.name]
            # Handle boolean mask
            if (isinstance(index_typ, types.npytypes.Array) and
                isinstance(index_typ.dtype, types.Boolean)):
                mask_var = ind
                mask_typ = index_typ.dtype
                mask_indices.append(None)
            # Handle integer array selector
            elif (isinstance(index_typ, types.npytypes.Array) and
                isinstance(index_typ.dtype, types.Integer)):
                mask_var = ind
                mask_typ = index_typ.dtype
                mask_indices.append(None)
            # Handle integer index
            elif isinstance(index_typ, types.Integer):
                count_consts += 1
                mask_indices.append(ind)

        require(mask_var and count_consts == ndim - 1)
        return value, mask_var, mask_typ, mask_indices
    raise GuardException


class ParforPass(ParforPassStates):

    """ParforPass class is responsible for converting NumPy
    calls in Numba intermediate representation to Parfors, which
    will lower into either sequential or parallel loops during lowering
    stage.
    """

    def _pre_run(self):
        # run array analysis, a pre-requisite for parfor translation
        self.array_analysis.run(self.func_ir.blocks)
        # NOTE: Prepare _the_max_label. See #6102
        ir_utils._the_max_label.update(
            ir_utils.find_max_label(self.func_ir.blocks))

    def run(self):
        """run parfor conversion pass: replace Numpy calls
        with Parfors when possible and optimize the IR."""
        self._pre_run()
        # run stencil translation to parfor
        if self.options.stencil:
            stencil_pass = StencilPass(self.func_ir, self.typemap,
                                       self.calltypes, self.array_analysis,
                                       self.typingctx, self.targetctx,
                                       self.flags)
            stencil_pass.run()
        if self.options.setitem:
            ConvertSetItemPass(self).run(self.func_ir.blocks)
        if self.options.numpy:
            ConvertNumpyPass(self).run(self.func_ir.blocks)
        if self.options.reduction:
            ConvertReducePass(self).run(self.func_ir.blocks)
        if self.options.prange:
            ConvertLoopPass(self).run(self.func_ir.blocks)
        if self.options.inplace_binop:
            ConvertInplaceBinop(self).run(self.func_ir.blocks)

        # setup diagnostics now parfors are found
        self.diagnostics.setup(self.func_ir, self.options.fusion)

        dprint_func_ir(self.func_ir, "after parfor pass")

    def _find_mask(self, arr_def):
        """check if an array is of B[...M...], where M is a
        boolean array, and other indices (if available) are ints.
        If found, return B, M, M's type, and a tuple representing mask indices.
        Otherwise, raise GuardException.
        """
        return _find_mask(self.typemap, self.func_ir, arr_def)

    def _mk_parfor_loops(self, size_vars, scope, loc):
        """
        Create loop index variables and build LoopNest objects for a parfor.
        """
        return _mk_parfor_loops(self.typemap, size_vars, scope, loc)


class ParforFusionPass(ParforPassStates):

    """ParforFusionPass class is responsible for fusing parfors
    """

    def run(self):
        """run parfor fusion pass"""

        # simplify CFG of parfor body loops since nested parfors with extra
        # jumps can be created with prange conversion
        n_parfors = simplify_parfor_body_CFG(self.func_ir.blocks)
        # simplify before fusion
        simplify(self.func_ir, self.typemap, self.calltypes, self.metadata["parfors"])
        # need two rounds of copy propagation to enable fusion of long sequences
        # of parfors like test_fuse_argmin (some PYTHONHASHSEED values since
        # apply_copies_parfor depends on set order for creating dummy assigns)
        simplify(self.func_ir, self.typemap, self.calltypes, self.metadata["parfors"])

        if self.options.fusion and n_parfors >= 2:
            self.func_ir._definitions = build_definitions(self.func_ir.blocks)
            self.array_analysis.equiv_sets = dict()
            self.array_analysis.run(self.func_ir.blocks)

            # Get parfor params to calculate reductions below.
            _, parfors = get_parfor_params(self.func_ir.blocks,
                                           self.options.fusion,
                                           self.nested_fusion_info)

            # Find reductions so that fusion can be disallowed if a
            # subsequent parfor read a reduction variable.
            for p in parfors:
                p.redvars, p.reddict = get_parfor_reductions(self.func_ir,
                                                             p,
                                                             p.params,
                                                             self.calltypes)

            # reorder statements to maximize fusion
            # push non-parfors down
            maximize_fusion(self.func_ir, self.func_ir.blocks, self.typemap,
                                                            up_direction=False)
            dprint_func_ir(self.func_ir, "after maximize fusion down")
            self.fuse_parfors(self.array_analysis,
                              self.func_ir.blocks,
                              self.func_ir,
                              self.typemap)
            dprint_func_ir(self.func_ir, "after first fuse")
            # push non-parfors up
            maximize_fusion(self.func_ir, self.func_ir.blocks, self.typemap)
            dprint_func_ir(self.func_ir, "after maximize fusion up")
            # try fuse again after maximize
            self.fuse_parfors(self.array_analysis,
                              self.func_ir.blocks,
                              self.func_ir,
                              self.typemap)
            dprint_func_ir(self.func_ir, "after fusion")
            # remove dead code after fusion to remove extra arrays and variables
            simplify(self.func_ir, self.typemap, self.calltypes, self.metadata["parfors"])

    def fuse_parfors(self, array_analysis, blocks, func_ir, typemap):
        for label, block in blocks.items():
            equiv_set = array_analysis.get_equiv_set(label)
            fusion_happened = True
            while fusion_happened:
                fusion_happened = False
                new_body = []
                i = 0
                while i < len(block.body) - 1:
                    stmt = block.body[i]
                    next_stmt = block.body[i + 1]
                    if isinstance(stmt, Parfor) and isinstance(next_stmt, Parfor):
                        # we have to update equiv_set since they have changed due to
                        # variables being renamed before fusion.
                        equiv_set = array_analysis.get_equiv_set(label)
                        stmt.equiv_set = equiv_set
                        next_stmt.equiv_set = equiv_set
                        fused_node, fuse_report = try_fuse(equiv_set, stmt, next_stmt,
                            self.metadata["parfors"], func_ir, typemap)
                        # accumulate fusion reports
                        self.diagnostics.fusion_reports.append(fuse_report)
                        if fused_node is not None:
                            fusion_happened = True
                            self.diagnostics.fusion_info[stmt.id].extend([next_stmt.id])
                            new_body.append(fused_node)
                            self.fuse_recursive_parfor(fused_node, equiv_set, func_ir, typemap)
                            i += 2
                            continue
                    new_body.append(stmt)
                    if isinstance(stmt, Parfor):
                        self.fuse_recursive_parfor(stmt, equiv_set, func_ir, typemap)
                    i += 1
                new_body.append(block.body[-1])
                block.body = new_body
        return

    def fuse_recursive_parfor(self, parfor, equiv_set, func_ir, typemap):
        blocks = wrap_parfor_blocks(parfor)
        maximize_fusion(self.func_ir, blocks, self.typemap)
        dprint_func_ir(self.func_ir, "after recursive maximize fusion down", blocks)
        arr_analysis = array_analysis.ArrayAnalysis(self.typingctx, self.func_ir,
                                                self.typemap, self.calltypes)
        arr_analysis.run(blocks, equiv_set)
        self.fuse_parfors(arr_analysis, blocks, func_ir, typemap)
        unwrap_parfor_blocks(parfor)


class ParforPreLoweringPass(ParforPassStates):

    """ParforPreLoweringPass class is responsible for preparing parfors for lowering.
    """

    def run(self):
        """run parfor prelowering pass"""

        # push function call variables inside parfors so gufunc function
        # wouldn't need function variables as argument
        push_call_vars(self.func_ir.blocks, {}, {}, self.typemap)
        dprint_func_ir(self.func_ir, "after push call vars")
        # simplify again
        simplify(self.func_ir, self.typemap, self.calltypes, self.metadata["parfors"])
        dprint_func_ir(self.func_ir, "after optimization")
        if config.DEBUG_ARRAY_OPT >= 1:
            print("variable types: ", sorted(self.typemap.items()))
            print("call types: ", self.calltypes)

        if config.DEBUG_ARRAY_OPT >= 3:
            for(block_label, block) in self.func_ir.blocks.items():
                new_block = []
                scope = block.scope
                for stmt in block.body:
                    new_block.append(stmt)
                    if isinstance(stmt, ir.Assign):
                        loc = stmt.loc
                        lhs = stmt.target
                        rhs = stmt.value
                        lhs_typ = self.typemap[lhs.name]
                        print("Adding print for assignment to ", lhs.name, lhs_typ, type(lhs_typ))
                        if lhs_typ in types.number_domain or isinstance(lhs_typ, types.Literal):
                            str_var = ir.Var(scope, mk_unique_var("str_var"), loc)
                            self.typemap[str_var.name] = types.StringLiteral(lhs.name)
                            lhs_const = ir.Const(lhs.name, loc)
                            str_assign = ir.Assign(lhs_const, str_var, loc)
                            new_block.append(str_assign)
                            str_print = ir.Print([str_var], None, loc)
                            self.calltypes[str_print] = signature(types.none, self.typemap[str_var.name])
                            new_block.append(str_print)
                            ir_print = ir.Print([lhs], None, loc)
                            self.calltypes[ir_print] = signature(types.none, lhs_typ)
                            new_block.append(ir_print)
                block.body = new_block

        if self.func_ir.is_generator:
            fix_generator_types(self.func_ir.generator_info, self.return_type,
                                self.typemap)
        if sequential_parfor_lowering:
            lower_parfor_sequential(
                self.typingctx, self.func_ir, self.typemap, self.calltypes, self.metadata)
        else:
            # prepare for parallel lowering
            # add parfor params to parfors here since lowering is destructive
            # changing the IR after this is not allowed
            parfor_ids, parfors = get_parfor_params(self.func_ir.blocks,
                                                    self.options.fusion,
                                                    self.nested_fusion_info)

            # Validate reduction in parfors.
            for p in parfors:
                p.redvars, p.reddict = get_parfor_reductions(self.func_ir,
                                                             p,
                                                             p.params,
                                                             self.calltypes)

            # Validate parameters:
            for p in parfors:
                p.validate_params(self.typemap)

            if config.DEBUG_ARRAY_OPT_STATS:
                name = self.func_ir.func_id.func_qualname
                n_parfors = len(parfor_ids)
                if n_parfors > 0:
                    after_fusion = ("After fusion" if self.options.fusion
                                    else "With fusion disabled")
                    print(('{}, function {} has '
                           '{} parallel for-loop(s) #{}.').format(
                           after_fusion, name, n_parfors, parfor_ids))
                else:
                    print('Function {} has no Parfor.'.format(name))


def _remove_size_arg(call_name, expr):
    "remove size argument from args or kws"
    # remove size kwarg
    kws = dict(expr.kws)
    kws.pop('size', '')
    expr.kws = tuple(kws.items())

    # remove size arg if available
    if call_name in random_1arg_size + random_int_args:
        # these calls have only a "size" argument or list of ints
        # so remove all args
        expr.args = []

    if call_name in random_3arg_sizelast:
        # normal, uniform, ... have 3 args, last one is size
        if len(expr.args) == 3:
            expr.args.pop()

    if call_name in random_2arg_sizelast:
        # have 2 args, last one is size
        if len(expr.args) == 2:
            expr.args.pop()

    if call_name == 'randint':
        # has 4 args, 3rd one is size
        if len(expr.args) == 3:
            expr.args.pop()
        if len(expr.args) == 4:
            dt_arg = expr.args.pop()
            expr.args.pop()  # remove size
            expr.args.append(dt_arg)

    if call_name == 'triangular':
        # has 4 args, last one is size
        if len(expr.args) == 4:
            expr.args.pop()

    return


def _get_call_arg_types(expr, typemap):
    new_arg_typs = []
    for arg in expr.args:
        new_arg_typs.append(typemap[arg.name])

    new_kw_types = {}
    for name, arg in expr.kws:
        new_kw_types[name] = typemap[arg.name]

    return tuple(new_arg_typs), new_kw_types


def _arrayexpr_tree_to_ir(
        func_ir,
        typingctx,
        typemap,
        calltypes,
        equiv_set,
        init_block,
        expr_out_var,
        expr,
        parfor_index_tuple_var,
        all_parfor_indices,
        avail_vars):
    """generate IR from array_expr's expr tree recursively. Assign output to
    expr_out_var and returns the whole IR as a list of Assign nodes.
    """
    el_typ = typemap[expr_out_var.name]
    scope = expr_out_var.scope
    loc = expr_out_var.loc
    out_ir = []

    if isinstance(expr, tuple):
        op, arr_expr_args = expr
        arg_vars = []
        for arg in arr_expr_args:
            arg_out_var = ir.Var(scope, mk_unique_var("$arg_out_var"), loc)
            typemap[arg_out_var.name] = el_typ
            out_ir += _arrayexpr_tree_to_ir(func_ir,
                                            typingctx,
                                            typemap,
                                            calltypes,
                                            equiv_set,
                                            init_block,
                                            arg_out_var,
                                            arg,
                                            parfor_index_tuple_var,
                                            all_parfor_indices,
                                            avail_vars)
            arg_vars.append(arg_out_var)
        if op in npydecl.supported_array_operators:
            el_typ1 = typemap[arg_vars[0].name]
            if len(arg_vars) == 2:
                el_typ2 = typemap[arg_vars[1].name]
                func_typ = typingctx.resolve_function_type(op, (el_typ1,
                                                                el_typ2), {})
                ir_expr = ir.Expr.binop(op, arg_vars[0], arg_vars[1], loc)
                if op == operator.truediv:
                    func_typ, ir_expr = _gen_np_divide(
                        arg_vars[0], arg_vars[1], out_ir, typemap)
            else:
                func_typ = typingctx.resolve_function_type(op, (el_typ1,), {})
                ir_expr = ir.Expr.unary(op, arg_vars[0], loc)
            calltypes[ir_expr] = func_typ
            el_typ = func_typ.return_type
            out_ir.append(ir.Assign(ir_expr, expr_out_var, loc))
        for T in array_analysis.MAP_TYPES:
            if isinstance(op, T):
                # elif isinstance(op, (np.ufunc, DUFunc)):
                # function calls are stored in variables which are not removed
                # op is typing_key to the variables type
                func_var_name = _find_func_var(typemap, op, avail_vars, loc=loc)
                func_var = ir.Var(scope, mk_unique_var(func_var_name), loc)
                typemap[func_var.name] = typemap[func_var_name]
                func_var_def = copy.deepcopy(func_ir.get_definition(func_var_name))
                if isinstance(func_var_def, ir.Expr) and func_var_def.op == 'getattr' and func_var_def.attr == 'sqrt':
                     g_math_var = ir.Var(scope, mk_unique_var("$math_g_var"), loc)
                     typemap[g_math_var.name] = types.misc.Module(math)
                     g_math = ir.Global('math', math, loc)
                     g_math_assign = ir.Assign(g_math, g_math_var, loc)
                     func_var_def = ir.Expr.getattr(g_math_var, 'sqrt', loc)
                     out_ir.append(g_math_assign)
#                     out_ir.append(func_var_def)
                ir_expr = ir.Expr.call(func_var, arg_vars, (), loc)
                call_typ = typemap[func_var.name].get_call_type(
                    typingctx, tuple(typemap[a.name] for a in arg_vars), {})
                calltypes[ir_expr] = call_typ
                el_typ = call_typ.return_type
                #signature(el_typ, el_typ)
                out_ir.append(ir.Assign(func_var_def, func_var, loc))
                out_ir.append(ir.Assign(ir_expr, expr_out_var, loc))
    elif isinstance(expr, ir.Var):
        var_typ = typemap[expr.name]
        if isinstance(var_typ, types.Array):
            el_typ = var_typ.dtype
            ir_expr = _gen_arrayexpr_getitem(
                equiv_set,
                expr,
                parfor_index_tuple_var,
                all_parfor_indices,
                el_typ,
                calltypes,
                typingctx,
                typemap,
                init_block,
                out_ir)
        else:
            # assert typemap[expr.name]==el_typ
            el_typ = var_typ
            ir_expr = expr
        out_ir.append(ir.Assign(ir_expr, expr_out_var, loc))
    elif isinstance(expr, ir.Const):
        el_typ = typing.Context().resolve_value_type(expr.value)
        out_ir.append(ir.Assign(expr, expr_out_var, loc))

    if len(out_ir) == 0:
        raise errors.UnsupportedRewriteError(
            f"Don't know how to translate array expression '{expr:r}'",
            loc=expr.loc,
        )
    typemap.pop(expr_out_var.name, None)
    typemap[expr_out_var.name] = el_typ
    return out_ir


def _gen_np_divide(arg1, arg2, out_ir, typemap):
    """generate np.divide() instead of / for array_expr to get numpy error model
    like inf for division by zero (test_division_by_zero).
    """
    scope = arg1.scope
    loc = arg1.loc
    # g_np_var = Global(numpy)
    g_np_var = ir.Var(scope, mk_unique_var("$np_g_var"), loc)
    typemap[g_np_var.name] = types.misc.Module(numpy)
    g_np = ir.Global('np', numpy, loc)
    g_np_assign = ir.Assign(g_np, g_np_var, loc)
    # attr call: div_attr = getattr(g_np_var, divide)
    div_attr_call = ir.Expr.getattr(g_np_var, "divide", loc)
    attr_var = ir.Var(scope, mk_unique_var("$div_attr"), loc)
    func_var_typ = get_np_ufunc_typ(numpy.divide)
    typemap[attr_var.name] = func_var_typ
    attr_assign = ir.Assign(div_attr_call, attr_var, loc)
    # divide call:  div_attr(arg1, arg2)
    div_call = ir.Expr.call(attr_var, [arg1, arg2], (), loc)
    func_typ = func_var_typ.get_call_type(
        typing.Context(), [typemap[arg1.name], typemap[arg2.name]], {})
    out_ir.extend([g_np_assign, attr_assign])
    return func_typ, div_call


def _gen_arrayexpr_getitem(
        equiv_set,
        var,
        parfor_index_tuple_var,
        all_parfor_indices,
        el_typ,
        calltypes,
        typingctx,
        typemap,
        init_block,
        out_ir):
    """if there is implicit dimension broadcast, generate proper access variable
    for getitem. For example, if indices are (i1,i2,i3) but shape is (c1,0,c3),
    generate a tuple with (i1,0,i3) for access.  Another example: for (i1,i2,i3)
    and (c1,c2) generate (i2,i3).
    """
    loc = var.loc
    index_var = parfor_index_tuple_var
    var_typ =  typemap[var.name]
    ndims = typemap[var.name].ndim
    num_indices = len(all_parfor_indices)
    size_vars = equiv_set.get_shape(var) or []
    size_consts = [equiv_set.get_equiv_const(x) for x in size_vars]
    # Handle array-scalar
    if ndims == 0:
        # call np.ravel
        ravel_var = ir.Var(var.scope, mk_unique_var("$ravel"), loc)
        ravel_typ = types.npytypes.Array(dtype=var_typ.dtype, ndim=1, layout='C')
        typemap[ravel_var.name] = ravel_typ
        stmts = ir_utils.gen_np_call('ravel', numpy.ravel, ravel_var, [var], typingctx, typemap, calltypes)
        init_block.body.extend(stmts)
        var = ravel_var
        # Const(0)
        const_node = ir.Const(0, var.loc)
        const_var = ir.Var(var.scope, mk_unique_var("$const_ind_0"), loc)
        typemap[const_var.name] = types.uintp
        const_assign = ir.Assign(const_node, const_var, loc)
        out_ir.append(const_assign)
        index_var = const_var
    # Handle 1d array
    elif ndims == 1:
        # Use last index for 1D arrays
        index_var = all_parfor_indices[-1]
    # Handle known constant size
    elif any([x is not None for x in size_consts]):
        # Need a tuple as index
        ind_offset = num_indices - ndims
        tuple_var = ir.Var(var.scope, mk_unique_var(
            "$parfor_index_tuple_var_bcast"), loc)
        typemap[tuple_var.name] = types.containers.UniTuple(types.uintp, ndims)
        # Just in case, const var for size 1 dim access index: $const0 =
        # Const(0)
        const_node = ir.Const(0, var.loc)
        const_var = ir.Var(var.scope, mk_unique_var("$const_ind_0"), loc)
        typemap[const_var.name] = types.uintp
        const_assign = ir.Assign(const_node, const_var, loc)
        out_ir.append(const_assign)
        index_vars = []
        for i in reversed(range(ndims)):
            size_var = size_vars[i]
            size_const = size_consts[i]
            if size_const == 1:
                index_vars.append(const_var)
            else:
                index_vars.append(all_parfor_indices[ind_offset + i])
        index_vars = list(reversed(index_vars))
        tuple_call = ir.Expr.build_tuple(index_vars, loc)
        tuple_assign = ir.Assign(tuple_call, tuple_var, loc)
        out_ir.append(tuple_assign)
        index_var = tuple_var

    ir_expr = ir.Expr.getitem(var, index_var, loc)
    calltypes[ir_expr] = signature(el_typ, typemap[var.name],
                                   typemap[index_var.name])
    return ir_expr


def _find_func_var(typemap, func, avail_vars, loc):
    """find variable in typemap which represents the function func.
    """
    for v in avail_vars:
        t = typemap[v]
        # Function types store actual functions in typing_key.
        if isinstance(t, Function) and t.typing_key == func:
            return v
    raise errors.UnsupportedRewriteError("ufunc call variable not found", loc=loc)


def lower_parfor_sequential(typingctx, func_ir, typemap, calltypes, metadata):
    ir_utils._the_max_label.update(ir_utils.find_max_label(func_ir.blocks))
    parfor_found = False
    new_blocks = {}
    scope = next(iter(func_ir.blocks.values())).scope
    for (block_label, block) in func_ir.blocks.items():
        block_label, parfor_found = _lower_parfor_sequential_block(
            block_label, block, new_blocks, typemap, calltypes, parfor_found,
            scope=scope)
        # old block stays either way
        new_blocks[block_label] = block
    func_ir.blocks = new_blocks
    # rename only if parfor found and replaced (avoid test_flow_control error)
    if parfor_found:
        func_ir.blocks = rename_labels(func_ir.blocks)
    dprint_func_ir(func_ir, "after parfor sequential lowering")
    simplify(func_ir, typemap, calltypes, metadata["parfors"])
    dprint_func_ir(func_ir, "after parfor sequential simplify")


def _lower_parfor_sequential_block(
        block_label,
        block,
        new_blocks,
        typemap,
        calltypes,
        parfor_found,
        scope):
    i = _find_first_parfor(block.body)
    while i != -1:
        parfor_found = True
        inst = block.body[i]
        loc = inst.init_block.loc
        # split block across parfor
        prev_block = ir.Block(scope, loc)
        prev_block.body = block.body[:i]
        block.body = block.body[i + 1:]
        # previous block jump to parfor init block
        init_label = next_label()
        prev_block.body.append(ir.Jump(init_label, loc))
        new_blocks[init_label] = transfer_scope(inst.init_block, scope)
        new_blocks[block_label] = prev_block
        block_label = next_label()

        ndims = len(inst.loop_nests)
        for i in range(ndims):
            loopnest = inst.loop_nests[i]
            # create range block for loop
            range_label = next_label()
            header_label = next_label()
            range_block = mk_range_block(
                typemap,
                loopnest.start,
                loopnest.stop,
                loopnest.step,
                calltypes,
                scope,
                loc)
            range_block.body[-1].target = header_label  # fix jump target
            phi_var = range_block.body[-2].target
            new_blocks[range_label] = range_block
            header_block = mk_loop_header(typemap, phi_var, calltypes,
                                          scope, loc)
            header_block.body[-2].target = loopnest.index_variable
            new_blocks[header_label] = header_block
            # jump to this new inner loop
            if i == 0:
                inst.init_block.body.append(ir.Jump(range_label, loc))
                header_block.body[-1].falsebr = block_label
            else:
                new_blocks[prev_header_label].body[-1].truebr = range_label
                header_block.body[-1].falsebr = prev_header_label
            prev_header_label = header_label  # to set truebr next loop

        # last body block jump to inner most header
        body_last_label = max(inst.loop_body.keys())
        inst.loop_body[body_last_label].body.append(
            ir.Jump(header_label, loc))
        # inner most header jumps to first body block
        body_first_label = min(inst.loop_body.keys())
        header_block.body[-1].truebr = body_first_label
        # add parfor body to blocks
        for (l, b) in inst.loop_body.items():
            l, parfor_found = _lower_parfor_sequential_block(
                l, b, new_blocks, typemap, calltypes, parfor_found,
                scope=scope)
            new_blocks[l] = transfer_scope(b, scope)
        i = _find_first_parfor(block.body)
    return block_label, parfor_found


def _find_first_parfor(body):
    for (i, inst) in enumerate(body):
        if isinstance(inst, Parfor) and not inst.no_sequential_lowering:
            return i
    return -1


def get_parfor_params(blocks, options_fusion, fusion_info):
    """find variables used in body of parfors from outside and save them.
    computed as live variables at entry of first block.
    """

    # since parfor wrap creates a back-edge to first non-init basic block,
    # live_map[first_non_init_block] contains variables defined in parfor body
    # that could be undefined before. So we only consider variables that are
    # actually defined before the parfor body in the program.
    parfor_ids = set()
    parfors = []
    pre_defs = set()
    _, all_defs = compute_use_defs(blocks)
    topo_order = find_topo_order(blocks)
    for label in topo_order:
        block = blocks[label]
        for i, parfor in _find_parfors(block.body):
            # find variable defs before the parfor in the same block
            dummy_block = ir.Block(block.scope, block.loc)
            dummy_block.body = block.body[:i]
            before_defs = compute_use_defs({0: dummy_block}).defmap[0]
            pre_defs |= before_defs
            params = get_parfor_params_inner(
                parfor, pre_defs, options_fusion, fusion_info,
            )
            parfor.params, parfor.races = _combine_params_races_for_ssa_names(
                block.scope, params, parfor.races,
            )
            parfor_ids.add(parfor.id)
            parfors.append(parfor)

        pre_defs |= all_defs[label]
    return parfor_ids, parfors


def _combine_params_races_for_ssa_names(scope, params, races):
    """Returns `(params|races1, races1)`, where `races1` contains all variables
    in `races` are NOT referring to the same unversioned (SSA) variables in
    `params`.
    """
    def unversion(k):
        try:
            return scope.get_exact(k).unversioned_name
        except ir.NotDefinedError:
            # XXX: it's a bug that something references an undefined name
            return k

    races1 = set(races)
    unver_params = list(map(unversion, params))

    for rv in races:
        if any(unversion(rv) == pv for pv in unver_params):
            races1.discard(rv)
        else:
            break

    return params | races1, races1


def get_parfor_params_inner(parfor, pre_defs, options_fusion, fusion_info):
    blocks = wrap_parfor_blocks(parfor)
    cfg = compute_cfg_from_blocks(blocks)
    usedefs = compute_use_defs(blocks)
    live_map = compute_live_map(cfg, blocks, usedefs.usemap, usedefs.defmap)
    parfor_ids, _ = get_parfor_params(blocks, options_fusion, fusion_info)
    n_parfors = len(parfor_ids)
    if n_parfors > 0:
        if config.DEBUG_ARRAY_OPT_STATS:
            after_fusion = ("After fusion" if options_fusion
                            else "With fusion disabled")
            print(('{}, parallel for-loop {} has '
                'nested Parfor(s) #{}.').format(
                after_fusion, parfor.id, n_parfors, parfor_ids))
        fusion_info[parfor.id] = list(parfor_ids)

    unwrap_parfor_blocks(parfor)
    keylist = sorted(live_map.keys())
    init_block = keylist[0]
    first_non_init_block = keylist[1]

    before_defs = usedefs.defmap[init_block] | pre_defs
    params = live_map[first_non_init_block] & before_defs
    return params


def _find_parfors(body):
    for i, inst in enumerate(body):
        if isinstance(inst, Parfor):
            yield i, inst


def _is_indirect_index(func_ir, index, nest_indices):
    index_def = guard(get_definition, func_ir, index.name)
    if isinstance(index_def, ir.Expr) and index_def.op == 'build_tuple':
        if [x.name for x in index_def.items] == [x.name for x in nest_indices]:
            return True

    return False


def get_array_indexed_with_parfor_index_internal(loop_body,
                                                 index,
                                                 ret_indexed,
                                                 ret_not_indexed,
                                                 nest_indices,
                                                 func_ir):
    for blk in loop_body:
        for stmt in blk.body:
            if isinstance(stmt, (ir.StaticSetItem, ir.SetItem)):
                setarray_index = get_index_var(stmt)
                if (isinstance(setarray_index, ir.Var) and
                    (setarray_index.name == index or
                     _is_indirect_index(
                         func_ir,
                         setarray_index,
                         nest_indices))):
                    ret_indexed.add(stmt.target.name)
                else:
                    ret_not_indexed.add(stmt.target.name)
            elif (isinstance(stmt, ir.Assign) and
                  isinstance(stmt.value, ir.Expr) and
                  stmt.value.op in ['getitem', 'static_getitem']):
                getarray_index = stmt.value.index
                getarray_name = stmt.value.value.name
                if (isinstance(getarray_index, ir.Var) and
                    (getarray_index.name == index or
                     _is_indirect_index(
                         func_ir,
                         getarray_index,
                         nest_indices))):
                    ret_indexed.add(getarray_name)
                else:
                    ret_not_indexed.add(getarray_name)
            elif isinstance(stmt, Parfor):
                get_array_indexed_with_parfor_index_internal(
                    stmt.loop_body.values(),
                    index,
                    ret_indexed,
                    ret_not_indexed,
                    nest_indices,
                    func_ir)


def get_array_indexed_with_parfor_index(loop_body,
                                        index,
                                        nest_indices,
                                        func_ir):
    ret_indexed = set()
    ret_not_indexed = set()
    get_array_indexed_with_parfor_index_internal(
        loop_body,
        index,
        ret_indexed,
        ret_not_indexed,
        nest_indices,
        func_ir)
    return ret_indexed, ret_not_indexed


def get_parfor_outputs(parfor, parfor_params):
    """get arrays that are written to inside the parfor and need to be passed
    as parameters to gufunc.
    """
    # FIXME: The following assumes the target of all SetItem are outputs,
    # which is wrong!
    last_label = max(parfor.loop_body.keys())
    outputs = []
    for blk in parfor.loop_body.values():
        for stmt in blk.body:
            if (isinstance(stmt, (ir.StaticSetItem, ir.SetItem)) and
                get_index_var(stmt).name == parfor.index_var.name):
                outputs.append(stmt.target.name)
    # make sure these written arrays are in parfor parameters (live coming in)
    outputs = list(set(outputs) & set(parfor_params))
    return sorted(outputs)


_RedVarInfo = make_dataclass(
    "_RedVarInfo",
    ["init_val", "reduce_nodes", "redop"],
    frozen=True,
)


def get_parfor_reductions(func_ir, parfor, parfor_params, calltypes, reductions=None,
        reduce_varnames=None, param_uses=None, param_nodes=None,
        var_to_param=None):
    """find variables that are updated using their previous values and an array
    item accessed with parfor index, e.g. s = s+A[i]
    """
    if reductions is None:
        reductions = {}
    if reduce_varnames is None:
        reduce_varnames = []

    # for each param variable, find what other variables are used to update it
    # also, keep the related nodes
    if param_uses is None:
        param_uses = defaultdict(list)
    if param_nodes is None:
        param_nodes = defaultdict(list)
    if var_to_param is None:
        var_to_param = {}

    blocks = wrap_parfor_blocks(parfor)
    topo_order = find_topo_order(blocks)
    topo_order = topo_order[1:]  # ignore init block
    unwrap_parfor_blocks(parfor)

    for label in reversed(topo_order):
        for stmt in reversed(parfor.loop_body[label].body):
            if (isinstance(stmt, ir.Assign)
                    and (stmt.target.name in parfor_params
                      or stmt.target.name in var_to_param)):
                lhs = stmt.target
                rhs = stmt.value
                cur_param = lhs if lhs.name in parfor_params else var_to_param[lhs.name]
                used_vars = []
                if isinstance(rhs, ir.Var):
                    used_vars = [rhs.name]
                elif isinstance(rhs, ir.Expr):
                    used_vars = [v.name for v in stmt.value.list_vars()]
                param_uses[cur_param].extend(used_vars)
                for v in used_vars:
                    var_to_param[v] = cur_param
                # save copy of dependent stmt
                stmt_cp = copy.deepcopy(stmt)
                if stmt.value in calltypes:
                    calltypes[stmt_cp.value] = calltypes[stmt.value]
                param_nodes[cur_param].append(stmt_cp)
            if isinstance(stmt, Parfor):
                # recursive parfors can have reductions like test_prange8
                get_parfor_reductions(func_ir, stmt, parfor_params, calltypes,
                    reductions, reduce_varnames, None, param_nodes, var_to_param)

    for param, used_vars in param_uses.items():
        # a parameter is a reduction variable if its value is used to update it
        # check reduce_varnames since recursive parfors might have processed
        # param already
        param_name = param.name
        if param_name in used_vars and param_name not in reduce_varnames:
            param_nodes[param].reverse()
            reduce_nodes = get_reduce_nodes(param, param_nodes[param], func_ir)
            # Certain kinds of ill-formed Python (like potentially undefined
            # variables) in combination with SSA can make things look like
            # reductions except that they don't have reduction operators.
            # If we get to this point but don't find a reduction operator
            # then assume it is this situation and just don't treat this
            # variable as a reduction.
            if reduce_nodes is not None:
                reduce_varnames.append(param_name)
                check_conflicting_reduction_operators(param, reduce_nodes)
                gri_out = guard(get_reduction_init, reduce_nodes)
                if gri_out is not None:
                    init_val, redop = gri_out
                else:
                    init_val = None
                    redop = None
                reductions[param_name] = _RedVarInfo(
                    init_val=init_val,
                    reduce_nodes=reduce_nodes,
                    redop=redop,
                )

    return reduce_varnames, reductions

def check_conflicting_reduction_operators(param, nodes):
    """In prange, a user could theoretically specify conflicting
       reduction operators.  For example, in one spot it is += and
       another spot *=.  Here, we raise an exception if multiple
       different reduction operators are used in one prange.
    """
    first_red_func = None
    for node in nodes:
        if (isinstance(node, ir.Assign) and
            isinstance(node.value, ir.Expr) and
            node.value.op=='inplace_binop'):
            if first_red_func is None:
                first_red_func = node.value.fn
            else:
                if first_red_func != node.value.fn:
                    msg = ("Reduction variable %s has multiple conflicting "
                           "reduction operators." % param.unversioned_name)
                    raise errors.UnsupportedRewriteError(msg, node.loc)

def get_reduction_init(nodes):
    """
    Get initial value for known reductions.
    Currently, only += and *= are supported.
    """
    require(len(nodes) >= 1)
    # there could be multiple extra assignments after the reduce node
    # See: test_reduction_var_reuse
    acc_expr = list(filter(lambda x: isinstance(x.value, ir.Expr), nodes))[-1].value
    require(isinstance(acc_expr, ir.Expr) and acc_expr.op in ['inplace_binop', 'binop'])
    acc_expr_fn = acc_expr.fn
    if acc_expr.op == 'binop':
        if acc_expr_fn == operator.add:
            acc_expr_fn = operator.iadd
        elif acc_expr_fn == operator.sub:
            acc_expr_fn = operator.isub
        elif acc_expr_fn == operator.mul:
            acc_expr_fn = operator.imul
        elif acc_expr_fn == operator.truediv:
            acc_expr_fn = operator.itruediv
    if acc_expr_fn == operator.iadd or acc_expr_fn == operator.isub:
        return 0, acc_expr_fn
    if (  acc_expr_fn == operator.imul
       or acc_expr_fn == operator.itruediv ):
        return 1, acc_expr_fn
    return None, None

def supported_reduction(x, func_ir):
    if x.op == 'inplace_binop' or x.op == 'binop':
        if x.fn == operator.ifloordiv or x.fn == operator.floordiv:
            raise errors.NumbaValueError(("Parallel floordiv reductions are not supported. "
                              "If all divisors are integers then a floordiv "
                              "reduction can in some cases be parallelized as "
                              "a multiply reduction followed by a floordiv of "
                              "the resulting product."), x.loc)
        supps = [operator.iadd,
                 operator.isub,
                 operator.imul,
                 operator.itruediv,
                 operator.add,
                 operator.sub,
                 operator.mul,
                 operator.truediv]
        return x.fn in supps
    if x.op == 'call':
        callname = guard(find_callname, func_ir, x)
        if callname in [
            ('max', 'builtins'),
            ('min', 'builtins'),
            ('datetime_minimum', 'numba.np.npdatetime_helpers'),
            ('datetime_maximum', 'numba.np.npdatetime_helpers'),
        ]:
            return True
    return False

def get_reduce_nodes(reduction_node, nodes, func_ir):
    """
    Get nodes that combine the reduction variable with a sentinel variable.
    Recognizes the first node that combines the reduction variable with another
    variable.
    """
    reduce_nodes = None
    defs = {}

    def cyclic_lookup(var, varonly=True, start=None):
        """Lookup definition of ``var``.
        Returns ``None`` if variable definition forms a cycle.
        """
        lookedup_var = defs.get(var.name, None)
        if isinstance(lookedup_var, ir.Var):
            if start is None:
                start = lookedup_var
            elif start == lookedup_var:
                # cycle detected
                return None
            return cyclic_lookup(lookedup_var, start=start)
        else:
            return var if (varonly or lookedup_var is None) else lookedup_var

    def noncyclic_lookup(*args, **kwargs):
        """Similar to cyclic_lookup but raise AssertionError if a cycle is
        detected.
        """
        res = cyclic_lookup(*args, **kwargs)
        if res is None:
            raise AssertionError("unexpected cycle in lookup()")
        return res

    name = reduction_node.name
    unversioned_name = reduction_node.unversioned_name
    for i, stmt in enumerate(nodes):
        lhs = stmt.target
        rhs = stmt.value
        defs[lhs.name] = rhs
        if isinstance(rhs, ir.Var) and rhs.name in defs:
            rhs = cyclic_lookup(rhs)
        if isinstance(rhs, ir.Expr):
            in_vars = set(noncyclic_lookup(v, True).name
                          for v in rhs.list_vars())
            if name in in_vars:
                # reductions like sum have an assignment afterwards
                # e.g. $2 = a + $1; a = $2
                # reductions that are functions calls like max() don't have an
                # extra assignment afterwards

                # This code was created when Numba had an IR generation strategy
                # where a binop for a reduction would be followed by an
                # assignment as follows:
                #$c.4.15 = inplace_binop(fn=<iadd>, ...>, lhs=c.3, rhs=$const20)
                #c.4 = $c.4.15

                # With Python 3.12 changes, Numba may separate that assignment
                # to a new basic block.  The code below looks and sees if an
                # assignment to the reduction var follows the reduction operator
                # and if not it searches the rest of the reduction nodes to find
                # the assignment that should follow the reduction operator
                # and then reorders the reduction nodes so that assignment
                # follows the reduction operator.
                if (i + 1 < len(nodes) and
                    ((not isinstance(nodes[i + 1], ir.Assign)) or
                     nodes[i + 1].target.unversioned_name != unversioned_name)):
                    foundj = None
                    # Iterate through the rest of the reduction nodes.
                    for j, jstmt in enumerate(nodes[i + 1:]):
                        # If this stmt is an assignment where the right-hand
                        # side of the assignment is the output of the reduction
                        # operator.
                        if isinstance(jstmt, ir.Assign) and jstmt.value == lhs:
                            # Remember the index of this node.  Because of
                            # nodes[i+1] above, we have to add i + 1 to j below
                            # to get the index in the original nodes list.
                            foundj = i + j + 1
                            break
                    if foundj is not None:
                        # If we found the correct assignment then move it to
                        # after the reduction operator.
                        nodes = (nodes[:i + 1] +   # nodes up to operator
                                 nodes[foundj:foundj + 1] + # assignment node
                                 nodes[i + 1:foundj] + # between op and assign
                                 nodes[foundj + 1:]) # after assignment node

                if (not (i+1 < len(nodes) and isinstance(nodes[i+1], ir.Assign)
                        and nodes[i+1].target.unversioned_name == unversioned_name)
                        and lhs.unversioned_name != unversioned_name):
                    raise ValueError(
                        f"Use of reduction variable {unversioned_name!r} other "
                        "than in a supported reduction function is not "
                        "permitted."
                    )

                if not supported_reduction(rhs, func_ir):
                    raise ValueError(("Use of reduction variable " + unversioned_name +
                                      " in an unsupported reduction function."))
                args = [(x.name, noncyclic_lookup(x, True))
                        for x in get_expr_args(rhs) ]
                non_red_args = [ x for (x, y) in args if y.name != name ]
                assert len(non_red_args) == 1
                args = [ (x, y) for (x, y) in args if x != y.name ]
                replace_dict = dict(args)
                replace_dict[non_red_args[0]] = ir.Var(lhs.scope, name+"#init", lhs.loc)
                replace_vars_inner(rhs, replace_dict)
                reduce_nodes = nodes[i:]
                break
    return reduce_nodes

def get_expr_args(expr):
    """
    Get arguments of an expression node
    """
    if expr.op in ['binop', 'inplace_binop']:
        return [expr.lhs, expr.rhs]
    if expr.op == 'call':
        return [v for v in expr.args]
    raise NotImplementedError("get arguments for expression {}".format(expr))

def visit_parfor_pattern_vars(parfor, callback, cbdata):
    # currently, only stencil pattern has variables
    for pattern in parfor.patterns:
        if pattern[0] == 'stencil':
            left_lengths = pattern[1][0]
            for i in range(len(left_lengths)):
                if isinstance(left_lengths[i], ir.Var):
                    left_lengths[i] = visit_vars_inner(left_lengths[i],
                                                            callback, cbdata)
            right_lengths = pattern[1][1]
            for i in range(len(right_lengths)):
                if isinstance(right_lengths[i], ir.Var):
                    right_lengths[i] = visit_vars_inner(right_lengths[i],
                                                            callback, cbdata)

def visit_vars_parfor(parfor, callback, cbdata):
    if config.DEBUG_ARRAY_OPT >= 1:
        print("visiting parfor vars for:", parfor)
        print("cbdata: ", sorted(cbdata.items()))
    for l in parfor.loop_nests:
        l.index_variable = visit_vars_inner(l.index_variable, callback, cbdata)
        if isinstance(l.start, ir.Var):
            l.start = visit_vars_inner(l.start, callback, cbdata)
        if isinstance(l.stop, ir.Var):
            l.stop = visit_vars_inner(l.stop, callback, cbdata)
        if isinstance(l.step, ir.Var):
            l.step = visit_vars_inner(l.step, callback, cbdata)
    visit_vars({-1: parfor.init_block}, callback, cbdata)
    visit_parfor_pattern_vars(parfor, callback, cbdata)
    visit_vars(parfor.loop_body, callback, cbdata)
    return


# add call to visit parfor variable
ir_utils.visit_vars_extensions[Parfor] = visit_vars_parfor


def parfor_defs(parfor, use_set=None, def_set=None):
    """list variables written in this parfor by recursively
    calling compute_use_defs() on body and combining block defs.
    """
    if use_set is None:
        use_set = set()
    if def_set is None:
        def_set = set()
    blocks = wrap_parfor_blocks(parfor)
    uses, defs = compute_use_defs(blocks)
    cfg = compute_cfg_from_blocks(blocks)
    last_label = max(blocks.keys())
    unwrap_parfor_blocks(parfor)

    # Conservatively, only add defs for blocks that are definitely executed
    # Go through blocks in order, as if they are statements of the block that
    # includes the parfor, and update uses/defs.

    # no need for topo order of ir_utils
    topo_order = cfg.topo_order()
    # blocks that dominate last block are definitely executed
    definitely_executed = cfg.dominators()[last_label]
    # except loop bodies that might not execute
    for loop in cfg.loops().values():
        definitely_executed -= loop.body
    for label in topo_order:
        if label in definitely_executed:
            # see compute_use_defs() in analysis.py
            # variables defined in the block that includes the parfor are not
            # uses of that block (are not potentially live in the beginning of
            # the block)
            use_set.update(uses[label] - def_set)
            def_set.update(defs[label])
        else:
            use_set.update(uses[label] - def_set)

    # treat loop variables and size variables as use
    loop_vars = {
        l.start.name for l in parfor.loop_nests if isinstance(
            l.start, ir.Var)}
    loop_vars |= {
        l.stop.name for l in parfor.loop_nests if isinstance(
            l.stop, ir.Var)}
    loop_vars |= {
        l.step.name for l in parfor.loop_nests if isinstance(
            l.step, ir.Var)}
    use_set.update(loop_vars - def_set)
    use_set |= get_parfor_pattern_vars(parfor)

    return analysis._use_defs_result(usemap=use_set, defmap=def_set)


analysis.ir_extension_usedefs[Parfor] = parfor_defs


def _parfor_use_alloca(parfor, alloca_set):
    """
    Reduction variables for parfors and the reduction variables within
    nested parfors must be stack allocated.
    """
    alloca_set |= set(parfor.redvars)

    blocks = wrap_parfor_blocks(parfor)
    alloca_set |= analysis.must_use_alloca(blocks)

    unwrap_parfor_blocks(parfor)


analysis.ir_extension_use_alloca[Parfor] = _parfor_use_alloca


def parfor_insert_dels(parfor, curr_dead_set):
    """insert dels in parfor. input: dead variable set right after parfor.
    returns the variables for which del was inserted.
    """
    blocks = wrap_parfor_blocks(parfor)
    cfg = compute_cfg_from_blocks(blocks)
    usedefs = compute_use_defs(blocks)
    live_map = compute_live_map(cfg, blocks, usedefs.usemap, usedefs.defmap)
    dead_map = compute_dead_maps(cfg, blocks, live_map, usedefs.defmap)

    # treat loop variables and size variables as live
    loop_vars = {
        l.start.name for l in parfor.loop_nests if isinstance(
            l.start, ir.Var)}
    loop_vars |= {
        l.stop.name for l in parfor.loop_nests if isinstance(
            l.stop, ir.Var)}
    loop_vars |= {
        l.step.name for l in parfor.loop_nests if isinstance(
            l.step, ir.Var)}
    loop_vars |= {l.index_variable.name for l in parfor.loop_nests}
    # for var_list in parfor.array_analysis.array_size_vars.values():
    #    loop_vars |= {v.name for v in var_list if isinstance(v, ir.Var)}

    dead_set = set()
    for label in blocks.keys():
        # only kill vars that are actually dead at the parfor's block
        dead_map.internal[label] &= curr_dead_set
        dead_map.internal[label] -= loop_vars
        dead_set |= dead_map.internal[label]
        dead_map.escaping[label] &= curr_dead_set
        dead_map.escaping[label] -= loop_vars
        dead_set |= dead_map.escaping[label]

    # dummy class to replace func_ir. _patch_var_dels only accesses blocks
    class DummyFuncIR(object):

        def __init__(self, blocks):
            self.blocks = blocks
    post_proc = postproc.PostProcessor(DummyFuncIR(blocks))
    post_proc._patch_var_dels(dead_map.internal, dead_map.escaping)
    unwrap_parfor_blocks(parfor)

    return dead_set | loop_vars


postproc.ir_extension_insert_dels[Parfor] = parfor_insert_dels


def maximize_fusion(func_ir, blocks, typemap, up_direction=True):
    """
    Reorder statements to maximize parfor fusion. Push all parfors up or down
    so they are adjacent.
    """
    call_table, _ = get_call_table(blocks)
    alias_map, arg_aliases = find_potential_aliases(
                                 blocks,
                                 func_ir.arg_names,
                                 typemap,
                                 func_ir
                             )
    for block in blocks.values():
        order_changed = True
        while order_changed:
            order_changed = maximize_fusion_inner(
                                func_ir,
                                block,
                                call_table,
                                alias_map,
                                arg_aliases,
                                up_direction
                            )

def maximize_fusion_inner(func_ir, block, call_table, alias_map,
                          arg_aliases, up_direction=True):
    order_changed = False
    i = 0
    # i goes to body[-3] (i+1 to body[-2]) since body[-1] is terminator and
    # shouldn't be reordered
    while i < len(block.body) - 2:
        stmt = block.body[i]
        next_stmt = block.body[i+1]
        can_reorder = (_can_reorder_stmts(stmt, next_stmt, func_ir,
                                          call_table, alias_map, arg_aliases)
                        if up_direction else _can_reorder_stmts(next_stmt, stmt,
                        func_ir, call_table, alias_map, arg_aliases))
        if can_reorder:
            block.body[i] = next_stmt
            block.body[i+1] = stmt
            order_changed = True
        i += 1
    return order_changed

def expand_aliases(the_set, alias_map, arg_aliases):
    ret = set()
    for i in the_set:
        if i in alias_map:
            ret = ret.union(alias_map[i])
        if i in arg_aliases:
            ret = ret.union(arg_aliases)
        ret.add(i)
    return ret

def _can_reorder_stmts(stmt, next_stmt, func_ir, call_table,
                       alias_map, arg_aliases):
    """
    Check dependencies to determine if a parfor can be reordered in the IR block
    with a non-parfor statement.
    """
    # swap only parfors with non-parfors
    # don't reorder calls with side effects (e.g. file close)
    # only read-read dependencies are OK
    # make sure there is no write-write, write-read dependencies
    if (isinstance(stmt, Parfor)
        and not isinstance(next_stmt, Parfor)
        and not isinstance(next_stmt, ir.Print)
        and (not isinstance(next_stmt, ir.Assign)
            or has_no_side_effect(next_stmt.value, set(), call_table)
            or guard(is_assert_equiv, func_ir, next_stmt.value))):
        stmt_accesses = expand_aliases({v.name for v in stmt.list_vars()},
                                       alias_map, arg_aliases)
        stmt_writes = expand_aliases(get_parfor_writes(stmt),
                                     alias_map, arg_aliases)
        next_accesses = expand_aliases({v.name for v in next_stmt.list_vars()},
                                       alias_map, arg_aliases)
        next_writes = expand_aliases(get_stmt_writes(next_stmt),
                                     alias_map, arg_aliases)
        if len((stmt_writes & next_accesses)
                | (next_writes & stmt_accesses)) == 0:
            return True
    return False

def is_assert_equiv(func_ir, expr):
    func_name, mod_name = find_callname(func_ir, expr)
    return func_name == 'assert_equiv'


def get_parfor_writes(parfor):
    assert isinstance(parfor, Parfor)
    writes = set()
    blocks = parfor.loop_body.copy()
    blocks[-1] = parfor.init_block
    for block in blocks.values():
        for stmt in block.body:
            writes.update(get_stmt_writes(stmt))
            if isinstance(stmt, Parfor):
                writes.update(get_parfor_writes(stmt))
    return writes

FusionReport = namedtuple('FusionReport', ['first', 'second', 'message'])

def try_fuse(equiv_set, parfor1, parfor2, metadata, func_ir, typemap):
    """try to fuse parfors and return a fused parfor, otherwise return None
    """
    dprint("try_fuse: trying to fuse \n", parfor1, "\n", parfor2)

    # default report is None
    report = None

    # fusion of parfors with different lowerers is not possible
    if parfor1.lowerer != parfor2.lowerer:
        dprint("try_fuse: parfors different lowerers")
        msg = "- fusion failed: lowerer mismatch"
        report = FusionReport(parfor1.id, parfor2.id, msg)
        return None, report

    # fusion of parfors with different dimensions not supported yet
    if len(parfor1.loop_nests) != len(parfor2.loop_nests):
        dprint("try_fuse: parfors number of dimensions mismatch")
        msg = "- fusion failed: number of loops mismatched, %s, %s."
        fmt = "parallel loop #%s has a nest of %s loops"
        l1 = fmt % (parfor1.id, len(parfor1.loop_nests))
        l2 = fmt % (parfor2.id, len(parfor2.loop_nests))
        report = FusionReport(parfor1.id, parfor2.id, msg % (l1, l2))
        return None, report

    ndims = len(parfor1.loop_nests)
    # all loops should be equal length

    def is_equiv(x, y):
        return x == y or equiv_set.is_equiv(x, y)

    def get_user_varname(v):
        """get original variable name by user if possible"""
        if not isinstance(v, ir.Var):
            return v
        v = v.name
        if "var_rename_map" in metadata and v in metadata["var_rename_map"]:
            user_varname = metadata["var_rename_map"][v]
            return user_varname
        return v

    for i in range(ndims):
        nest1 = parfor1.loop_nests[i]
        nest2 = parfor2.loop_nests[i]
        if not (is_equiv(nest1.start, nest2.start) and
                is_equiv(nest1.stop, nest2.stop) and
                is_equiv(nest1.step, nest2.step)):
            dprint("try_fuse: parfor dimension correlation mismatch", i)
            msg = "- fusion failed: loop dimension mismatched in axis %s. "
            msg += "slice(%s, %s, %s) != " % (get_user_varname(nest1.start),
                get_user_varname(nest1.stop), get_user_varname(nest1.step))
            msg += "slice(%s, %s, %s)" % (get_user_varname(nest2.start),
                get_user_varname(nest2.stop), get_user_varname(nest2.step))
            report = FusionReport(parfor1.id, parfor2.id, msg % i)
            return None, report

    func_ir._definitions = build_definitions(func_ir.blocks)
    p1_cross_dep, p1_ip, p1_ia, p1_non_ia = has_cross_iter_dep(parfor1, func_ir, typemap)
    if not p1_cross_dep:
        p2_cross_dep = has_cross_iter_dep(parfor2, func_ir, typemap, p1_ip, p1_ia, p1_non_ia)[0]
    else:
        p2_cross_dep = True

    if p1_cross_dep or p2_cross_dep:
        dprint("try_fuse: parfor cross iteration dependency found")
        msg = ("- fusion failed: cross iteration dependency found "
                "between loops #%s and #%s")
        report = FusionReport(parfor1.id, parfor2.id,
                                msg % (parfor1.id, parfor2.id))
        return None, report


    # find parfor1's defs, only body is considered since init_block will run
    # first after fusion as well
    p1_body_usedefs = compute_use_defs(parfor1.loop_body)
    p1_body_defs = set()
    for defs in p1_body_usedefs.defmap.values():
        p1_body_defs |= defs
    p1_body_defs |= get_parfor_writes(parfor1)
    # Add reduction variables from parfor1 to the set of body defs
    # so that if parfor2 reads the reduction variable it won't fuse.
    p1_body_defs |= set(parfor1.redvars)

    p2_usedefs = compute_use_defs(parfor2.loop_body)
    p2_uses = compute_use_defs({0: parfor2.init_block}).usemap[0]
    for uses in p2_usedefs.usemap.values():
        p2_uses |= uses

    overlap = p1_body_defs.intersection(p2_uses)
    # overlap are those variable defined in first parfor and used in the second
    if len(overlap) != 0:
        # Get all the arrays
        _, p2arraynotindexed = get_array_indexed_with_parfor_index(
            parfor2.loop_body.values(),
            parfor2.index_var.name,
            parfor2.get_loop_nest_vars(),
            func_ir)

        unsafe_var = (not isinstance(typemap[x], types.ArrayCompatible) or x in p2arraynotindexed for x in overlap)

        if any(unsafe_var):
            dprint("try_fuse: parfor2 depends on parfor1 body")
            msg = ("- fusion failed: parallel loop %s has a dependency on the "
                    "body of parallel loop %s. ")
            report = FusionReport(parfor1.id, parfor2.id,
                                    msg % (parfor1.id, parfor2.id))
            return None, report

    return fuse_parfors_inner(parfor1, parfor2)


def fuse_parfors_inner(parfor1, parfor2):
    # fuse parfor2 into parfor1
    # append parfor2's init block on parfor1's
    parfor1.init_block.body.extend(parfor2.init_block.body)

    # append parfor2's first block to parfor1's last block
    parfor2_first_label = min(parfor2.loop_body.keys())
    parfor2_first_block = parfor2.loop_body[parfor2_first_label].body
    parfor1_first_label = min(parfor1.loop_body.keys())
    parfor1_last_label = max(parfor1.loop_body.keys())
    parfor1.loop_body[parfor1_last_label].body.extend(parfor2_first_block)

    # add parfor2 body blocks to parfor1's except first
    parfor1.loop_body.update(parfor2.loop_body)
    parfor1.loop_body.pop(parfor2_first_label)

    # replace parfor2 indices with parfor1's
    ndims = len(parfor1.loop_nests)
    index_dict = {parfor2.index_var.name: parfor1.index_var}
    for i in range(ndims):
        index_dict[parfor2.loop_nests[i].index_variable.name] = parfor1.loop_nests[
            i].index_variable
    replace_vars(parfor1.loop_body, index_dict)

    # re-order labels from min to max
    blocks = wrap_parfor_blocks(parfor1, entry_label=parfor1_first_label)
    blocks = rename_labels(blocks)
    unwrap_parfor_blocks(parfor1, blocks)

    nameset = set(x.name for x in index_dict.values())
    remove_duplicate_definitions(parfor1.loop_body, nameset)
    parfor1.patterns.extend(parfor2.patterns)
    if config.DEBUG_ARRAY_OPT_STATS:
        print('Parallel for-loop #{} is fused into for-loop #{}.'.format(
              parfor2.id, parfor1.id))

    msg = '- fusion succeeded: parallel for-loop #{} is fused into for-loop #{}.'
    msg = msg.format(parfor2.id, parfor1.id)
    report = FusionReport(parfor1.id, parfor2.id, msg)

    return parfor1, report


def remove_duplicate_definitions(blocks, nameset):
    """Remove duplicated definition for variables in the given nameset, which
    is often a result of parfor fusion.
    """
    for label, block in blocks.items():
        body = block.body
        new_body = []
        defined = set()
        for inst in body:
            if isinstance(inst, ir.Assign):
                name = inst.target.name
                if name in nameset:
                    if name in defined:
                        continue
                    defined.add(name)
            new_body.append(inst)
        block.body = new_body
    return


def has_cross_iter_dep(
        parfor,
        func_ir,
        typemap,
        index_positions=None,
        indexed_arrays=None,
        non_indexed_arrays=None):
    # We should assume there is cross iteration dependency unless we can
    # prove otherwise.  Return True if there is a cross-iter dependency
    # that should prevent fusion, False if fusion is okay.
    # Also returns index_positions, indexed_arrays, and non_indexed_arrays
    # who purpose is described below so that subsequent additional
    # has_cross_iter_dep calls for other parfors can build on the same
    # data structures to make sure that the array accesses generate no
    # cross-iter dependencies both within a parfor but also across parfors.

    # TODO: make it more accurate using ud-chains

    # Get the index variable used by this parfor.
    # This will hold other variables with equivalent value, e.g., a = index_var
    indices = {l.index_variable.name for l in parfor.loop_nests}
    # This set will store variables that are (potentially recursively)
    # defined in relation to an index variable, e.g., a = index_var + 1.
    # A getitem that uses an index variable from this set will be considered
    # as potentially having a cross-iter dependency and so won't fuse.
    derived_from_indices = set()
    # For the first parfor considered for fusion, the latter 3 args will be None
    # and initialized to empty.  For the second parfor, the structures from the
    # previous parfor are passed in so that cross-parfor violations of the
    # below comments can prevent fusion.
    #
    # index_positions keeps track of which index positions have had an index
    # variable used for them and which ones haven't for each possible array
    # dimensionality.  After the first array access is seen, if subsequent
    # ones use a parfor index for a different dimension then we conservatively
    # say that we can't fuse.  For example, if a 2D array is accessed with
    # a[parfor_index, 0] then index_positions[2] will be (True, False) and
    # if a[0, parfor_index] happens later which is (False, True) then this
    # conflicts with the previous value and will prevent fusion.
    #
    # indexed_arrays records arrays that are accessed with at least one
    # parfor index.  If such an array is later accessed with indices that
    # don't include a parfor index then conservatively assume we can't fuse.
    #
    # non_indexed_arrays holds arrays that are indexed without a parfor index.
    # If an array first accessed without a parfor index is later indexed
    # with one then conservatively assume we can't fuse.
    if index_positions is None:
        index_positions = {}
    if indexed_arrays is None:
        indexed_arrays = set()
    if non_indexed_arrays is None:
        non_indexed_arrays = set()

    def add_check_position(new_position,
                           array_accessed,
                           index_positions,
                           indexed_arrays,
                           non_indexed_arrays):
        """Returns True if there is a reason to prevent fusion based
           on the rules described above.
           new_position will be a list or tuple of booleans that
           says whether the index in that spot is a parfor index
           or not.  array_accessed is the array on which the access
           is occurring."""

        # Convert list indices to tuple for generality.
        if isinstance(new_position, list):
            new_position = tuple(new_position)

        # If none of the indices are based on a parfor index.
        if True not in new_position:
            # See if this array has been accessed before with a
            # a parfor index and if so say that we can't fuse.
            if array_accessed in indexed_arrays:
                return True
            else:
                # Either array is already in non_indexed arrays or we
                # will add it.  Either way, this index usage can fuse.
                non_indexed_arrays.add(array_accessed)
                return False

        # Fallthrough for cases where one of the indices is a parfor index.
        # If this array was previously accessed without a parfor index then
        # conservatively say we can't fuse.
        if array_accessed in non_indexed_arrays:
            return True

        indexed_arrays.add(array_accessed)

        npsize = len(new_position)
        # See if we have not seen a npsize dimensioned array accessed before.
        if npsize not in index_positions:
            # If not then add current set of parfor/non-parfor indices and
            # indicate it is safe as it is the first usage.
            index_positions[npsize] = new_position
            return False

        # Here we have a subsequent access to a npsize-dimensioned array.
        # Make sure we see the same combination of parfor/non-parfor index
        # indices that we've seen before.  If not then return True saying
        # that we can't fuse.
        return index_positions[npsize] != new_position

    def check_index(stmt_index,
                    array_accessed,
                    index_positions,
                    indexed_arrays,
                    non_indexed_arrays,
                    derived_from_indices):
        """Looks at the indices of a getitem or setitem to see if there
           is a reason that they would prevent fusion.
           Returns True if fusion should be prohibited, False otherwise.
        """
        if isinstance(stmt_index, ir.Var):
            # If the array is 2+ dimensions then the index should be a tuple.
            if isinstance(typemap[stmt_index.name], types.BaseTuple):
                # Get how the index tuple is constructed.
                fbs_res = guard(find_build_sequence, func_ir, stmt_index)
                if fbs_res is not None:
                    ind_seq, _ = fbs_res
                    # If any indices are derived from an index is used then
                    # return True to say we can't fuse.
                    if (all([x.name in indices or
                        x.name not in derived_from_indices for x in ind_seq])):
                        # Get position in index tuple where parfor indices used.
                        new_index_positions = [x.name in indices for x in ind_seq]
                        # Make sure that we aren't accessing a given array with
                        # different indices in a different order.
                        return add_check_position(new_index_positions,
                                                  array_accessed,
                                                  index_positions,
                                                  indexed_arrays,
                                                  non_indexed_arrays)
                    else:
                        # index derived from a parfor index used so no fusion
                        return True
                else:
                    # Don't know how the index tuple is built so
                    # have to assume fusion can't happen.
                    return True
            else:
                # Should be for 1D arrays.
                if stmt_index.name in indices:
                    # Array indexed by a parfor index variable.
                    # Make sure this 1D access is consistent with prior ones.
                    return add_check_position((True,),
                                              array_accessed,
                                              index_positions,
                                              indexed_arrays,
                                              non_indexed_arrays)
                elif stmt_index.name in derived_from_indices:
                    # If we ever index an array with something calculated
                    # from an index then no fusion.
                    return True
                else:
                    # Some kind of index that isn't a parfor index or
                    # one derived from one, e.g., a constant.  Make sure
                    # this is consistent with prior accessed of this array.
                    return add_check_position((False,),
                                              array_accessed,
                                              index_positions,
                                              indexed_arrays,
                                              non_indexed_arrays)
        else:
            # We don't know how to handle non-Var indices so no fusion.
            return True

        # All branches above should cover all the cases and each should
        # return so we should never get here.
        raise errors.InternalError("Some code path in the parfor fusion "
                                   "cross-iteration dependency checker "
                                   "check_index didn't return a result.")

    # Iterate through all the statements in the parfor.
    for b in parfor.loop_body.values():
        for stmt in b.body:
            # Make sure SetItem accesses are fusion safe.
            if isinstance(stmt, (ir.SetItem, ir.StaticSetItem)):
                if isinstance(typemap[stmt.target.name], types.npytypes.Array):
                    # Check index safety with prior array accesses.
                    if check_index(stmt.index,
                                   stmt.target.name,
                                   index_positions,
                                   indexed_arrays,
                                   non_indexed_arrays,
                                   derived_from_indices):
                        return True, index_positions, indexed_arrays, non_indexed_arrays
                # Fusion safe so go to next statement.
                continue
            elif isinstance(stmt, ir.Assign):
                # If stmt of form a = parfor_index then add "a" to set of
                # parfor indices.
                if isinstance(stmt.value, ir.Var):
                    if stmt.value.name in indices:
                        indices.add(stmt.target.name)
                        continue
                elif isinstance(stmt.value, ir.Expr):
                    op = stmt.value.op
                    # Make sure getitem accesses are fusion safe.
                    if op in ['getitem', 'static_getitem']:
                        if isinstance(typemap[stmt.value.value.name], types.npytypes.Array):
                            # Check index safety with prior array accesses.
                            if check_index(stmt.value.index,
                                           stmt.value.value.name,
                                           index_positions,
                                           indexed_arrays,
                                           non_indexed_arrays,
                                           derived_from_indices):
                                return True, index_positions, indexed_arrays, non_indexed_arrays
                        # Fusion safe so go to next statement.
                        continue
                    elif op == 'call':
                        # If there is a call in the parfor body that takes some
                        # array parameter then we have no way to analyze what
                        # that call is doing so presume it is unsafe for fusion.
                        if (any([isinstance(typemap[x.name], types.npytypes.Array)
                            for x in stmt.value.list_vars()])):
                            return True, index_positions, indexed_arrays, non_indexed_arrays

                    # Get the vars used by this non-setitem/getitem statement.
                    rhs_vars = [x.name for x in stmt.value.list_vars()]
                    # If a parfor index is used as part of this statement or
                    # something previous determined to be derived from a parfor
                    # index then add the target variable to the set of
                    # variables that are derived from parfors and so should
                    # prevent fusion if used as an index.
                    if (not indices.isdisjoint(rhs_vars) or
                        not derived_from_indices.isdisjoint(rhs_vars)):
                        derived_from_indices.add(stmt.target.name)

    return False, index_positions, indexed_arrays, non_indexed_arrays


def dprint(*s):
    if config.DEBUG_ARRAY_OPT >= 1:
        print(*s)

def get_parfor_pattern_vars(parfor):
    """ get the variables used in parfor pattern information
    """
    out = set()
    # currently, only stencil pattern has variables
    for pattern in parfor.patterns:
        if pattern[0] == 'stencil':
            left_lengths = pattern[1][0]
            right_lengths = pattern[1][1]
            for v in left_lengths+right_lengths:
                if isinstance(v, ir.Var):
                    out.add(v.name)
    return out

def remove_dead_parfor(parfor, lives, lives_n_aliases, arg_aliases, alias_map, func_ir, typemap):
    """ remove dead code inside parfor including get/sets
    """

    with dummy_return_in_loop_body(parfor.loop_body):
        labels = find_topo_order(parfor.loop_body)

    # get/setitem replacement should ideally use dataflow to propagate setitem
    # saved values, but for simplicity we handle the common case of propagating
    # setitems in the first block (which is dominant) if the array is not
    # potentially changed in any way
    first_label = labels[0]
    first_block_saved_values = {}
    _update_parfor_get_setitems(
        parfor.loop_body[first_label].body,
        parfor.index_var, alias_map,
        first_block_saved_values,
        lives_n_aliases
        )

    # remove saved first block setitems if array potentially changed later
    saved_arrs = set(first_block_saved_values.keys())
    for l in labels:
        if l == first_label:
            continue
        for stmt in parfor.loop_body[l].body:
            if (isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Expr)
                    and stmt.value.op == 'getitem'
                    and stmt.value.index.name == parfor.index_var.name):
                continue
            varnames = set(v.name for v in stmt.list_vars())
            rm_arrs = varnames & saved_arrs
            for a in rm_arrs:
                first_block_saved_values.pop(a, None)


    # replace getitems with available value
    # e.g. A[i] = v; ... s = A[i]  ->  s = v
    for l in labels:
        if l == first_label:
            continue
        block = parfor.loop_body[l]
        saved_values = first_block_saved_values.copy()
        _update_parfor_get_setitems(block.body, parfor.index_var, alias_map,
                                        saved_values, lives_n_aliases)


    # after getitem replacement, remove extra setitems
    blocks = parfor.loop_body.copy()  # shallow copy is enough
    last_label = max(blocks.keys())
    return_label, tuple_var = _add_liveness_return_block(blocks, lives_n_aliases, typemap)
    # jump to return label
    jump = ir.Jump(return_label, ir.Loc("parfors_dummy", -1))
    blocks[last_label].body.append(jump)
    cfg = compute_cfg_from_blocks(blocks)
    usedefs = compute_use_defs(blocks)
    live_map = compute_live_map(cfg, blocks, usedefs.usemap, usedefs.defmap)
    alias_set = set(alias_map.keys())

    for label, block in blocks.items():
        new_body = []
        in_lives = {v.name for v in block.terminator.list_vars()}
        # find live variables at the end of block
        for out_blk, _data in cfg.successors(label):
            in_lives |= live_map[out_blk]
        for stmt in reversed(block.body):
            # aliases of lives are also live for setitems
            alias_lives = in_lives & alias_set
            for v in alias_lives:
                in_lives |= alias_map[v]
            if (isinstance(stmt, (ir.StaticSetItem, ir.SetItem)) and
                get_index_var(stmt).name == parfor.index_var.name and
                stmt.target.name not in in_lives and
                stmt.target.name not in arg_aliases):
                continue
            in_lives |= {v.name for v in stmt.list_vars()}
            new_body.append(stmt)
        new_body.reverse()
        block.body = new_body

    typemap.pop(tuple_var.name)  # remove dummy tuple type
    blocks[last_label].body.pop()  # remove jump

    """
      Process parfor body recursively.
      Note that this is the only place in this function that uses the
      argument lives instead of lives_n_aliases.  The former does not
      include the aliases of live variables but only the live variable
      names themselves.  See a comment in this function for how that
      is used.
    """
    remove_dead_parfor_recursive(
        parfor, lives, arg_aliases, alias_map, func_ir, typemap)

    # remove parfor if empty
    is_empty = len(parfor.init_block.body) == 0
    for block in parfor.loop_body.values():
        is_empty &= len(block.body) == 0
    if is_empty:
        return None
    return parfor

def _update_parfor_get_setitems(block_body, index_var, alias_map,
                                  saved_values, lives):
    """
    replace getitems of a previously set array in a block of parfor loop body
    """
    for stmt in block_body:
        if (isinstance(stmt, (ir.StaticSetItem, ir.SetItem)) and
            get_index_var(stmt).name == index_var.name and
            stmt.target.name not in lives):
            # saved values of aliases of SetItem target array are invalid
            for w in alias_map.get(stmt.target.name, []):
                saved_values.pop(w, None)
            # set saved value after invalidation since alias_map may
            # contain the array itself (e.g. pi example)
            saved_values[stmt.target.name] = stmt.value
            continue
        if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Expr):
            rhs = stmt.value
            if rhs.op == 'getitem' and isinstance(rhs.index, ir.Var):
                if rhs.index.name == index_var.name:
                    # replace getitem if value saved
                    stmt.value = saved_values.get(rhs.value.name, rhs)
                    continue
        # conservative assumption: array is modified if referenced
        # remove all referenced arrays
        for v in stmt.list_vars():
            saved_values.pop(v.name, None)
            # aliases are potentially modified as well
            for w in alias_map.get(v.name, []):
                saved_values.pop(w, None)

    return

ir_utils.remove_dead_extensions[Parfor] = remove_dead_parfor


def remove_dead_parfor_recursive(parfor, lives, arg_aliases, alias_map,
                                                             func_ir, typemap):
    """create a dummy function from parfor and call remove dead recursively
    """
    blocks = parfor.loop_body.copy()  # shallow copy is enough
    first_body_block = min(blocks.keys())
    assert first_body_block > 0  # we are using 0 for init block here
    last_label = max(blocks.keys())

    """
      Previously, this statement used lives_n_aliases.  That had the effect of
      keeping variables in the init_block alive if they aliased an array that
      was later written to.  By using just lives to indicate which variables
      names are live at exit of the parfor but then using alias_map for the
      actual recursive dead code removal, we keep any writes to aliased arrays
      alive but also allow aliasing assignments (i.e., a = b) to be eliminated
      so long as 'b' is not written to through the variable 'a' later on.
      This makes assignment handling of remove_dead_block work properly since
      it allows distinguishing between live variables and their aliases.
    """
    return_label, tuple_var = _add_liveness_return_block(blocks, lives, typemap)

    # branch back to first body label to simulate loop
    scope = blocks[last_label].scope

    branchcond = ir.Var(scope, mk_unique_var("$branchcond"), ir.Loc("parfors_dummy", -1))
    typemap[branchcond.name] = types.boolean

    branch = ir.Branch(branchcond, first_body_block, return_label, ir.Loc("parfors_dummy", -1))
    blocks[last_label].body.append(branch)

    # add dummy jump in init_block for CFG to work
    blocks[0] = parfor.init_block
    blocks[0].body.append(ir.Jump(first_body_block, ir.Loc("parfors_dummy", -1)))

    # args var including aliases is ok
    remove_dead(blocks, arg_aliases, func_ir, typemap, alias_map, arg_aliases)
    typemap.pop(tuple_var.name)  # remove dummy tuple type
    blocks[0].body.pop()  # remove dummy jump
    blocks[last_label].body.pop()  # remove branch
    return

def _add_liveness_return_block(blocks, lives, typemap):
    last_label = max(blocks.keys())
    return_label = last_label + 1

    loc = blocks[last_label].loc
    scope = blocks[last_label].scope
    blocks[return_label] = ir.Block(scope, loc)

    # add lives in a dummpy return to last block to avoid their removal
    tuple_var = ir.Var(scope, mk_unique_var("$tuple_var"), loc)
    # dummy type for tuple_var
    typemap[tuple_var.name] = types.containers.UniTuple(
        types.uintp, 2)
    live_vars = [ir.Var(scope, v, loc) for v in lives]
    tuple_call = ir.Expr.build_tuple(live_vars, loc)
    blocks[return_label].body.append(ir.Assign(tuple_call, tuple_var, loc))
    blocks[return_label].body.append(ir.Return(tuple_var, loc))
    return return_label, tuple_var


def find_potential_aliases_parfor(parfor, args, typemap, func_ir, alias_map, arg_aliases):
    blocks = wrap_parfor_blocks(parfor)
    ir_utils.find_potential_aliases(
        blocks, args, typemap, func_ir, alias_map, arg_aliases)
    unwrap_parfor_blocks(parfor)
    return

ir_utils.alias_analysis_extensions[Parfor] = find_potential_aliases_parfor

def simplify_parfor_body_CFG(blocks):
    """simplify CFG of body loops in parfors"""
    n_parfors = 0
    for block in blocks.values():
        for stmt in block.body:
            if isinstance(stmt, Parfor):
                n_parfors += 1
                parfor = stmt
                # add dummy return to enable CFG creation
                # can't use dummy_return_in_loop_body since body changes
                last_block = parfor.loop_body[max(parfor.loop_body.keys())]
                scope = last_block.scope
                loc = ir.Loc("parfors_dummy", -1)
                const = ir.Var(scope, mk_unique_var("$const"), loc)
                last_block.body.append(ir.Assign(ir.Const(0, loc), const, loc))
                last_block.body.append(ir.Return(const, loc))
                parfor.loop_body = simplify_CFG(parfor.loop_body)
                last_block = parfor.loop_body[max(parfor.loop_body.keys())]
                last_block.body.pop()
                # call on body recursively
                simplify_parfor_body_CFG(parfor.loop_body)
    return n_parfors


def wrap_parfor_blocks(parfor, entry_label = None):
    """wrap parfor blocks for analysis/optimization like CFG"""
    blocks = parfor.loop_body.copy()  # shallow copy is enough
    if entry_label is None:
        entry_label = min(blocks.keys())
    assert entry_label > 0  # we are using 0 for init block here

    # add dummy jump in init_block for CFG to work
    blocks[0] = parfor.init_block
    blocks[0].body.append(ir.Jump(entry_label, blocks[0].loc))
    for block in blocks.values():
        if len(block.body) == 0 or (not block.body[-1].is_terminator):
            block.body.append(ir.Jump(entry_label, block.loc))
    return blocks


def unwrap_parfor_blocks(parfor, blocks=None):
    """
    unwrap parfor blocks after analysis/optimization.
    Allows changes to the parfor loop.
    """
    if blocks is not None:
        # make sure init block isn't removed
        init_block_label = min(blocks.keys())
        # update loop body blocks
        blocks.pop(init_block_label)
        parfor.loop_body = blocks

    # make sure dummy jump to loop body isn't altered
    first_body_label = min(parfor.loop_body.keys())
    assert isinstance(parfor.init_block.body[-1], ir.Jump)

    # remove dummy jump to loop body
    parfor.init_block.body.pop()

    # make sure dummy jump back to loop body isn't altered
    for block in parfor.loop_body.values():
        if (isinstance(block.body[-1], ir.Jump) and
            block.body[-1].target == first_body_label):
            # remove dummy jump back to loop
            block.body.pop()
    return


def get_copies_parfor(parfor, typemap):
    """find copies generated/killed by parfor"""
    blocks = wrap_parfor_blocks(parfor)
    in_copies_parfor, out_copies_parfor = copy_propagate(blocks, typemap)
    in_gen_copies, in_extra_kill = get_block_copies(blocks, typemap)
    unwrap_parfor_blocks(parfor)

    # parfor's extra kill is kills of its init block,
    # and all possible gens and kills of it's body loop.
    # body doesn't gen and only kills since it may or may not run
    # TODO: save copies that are repeated in parfor
    kill_set = in_extra_kill[0]
    for label in parfor.loop_body.keys():
        kill_set |= {l for l, r in in_gen_copies[label]}
        kill_set |= in_extra_kill[label]

    # gen copies is copies generated by init that are not killed by body
    last_label = max(parfor.loop_body.keys())
    gens = out_copies_parfor[last_label] & in_gen_copies[0]

    if config.DEBUG_ARRAY_OPT >= 1:
        print("copy propagate parfor gens:", gens, "kill_set", kill_set)
    return gens, kill_set


ir_utils.copy_propagate_extensions[Parfor] = get_copies_parfor


def apply_copies_parfor(parfor, var_dict, name_var_table,
                        typemap, calltypes, save_copies):
    """apply copy propagate recursively in parfor"""
    # replace variables in pattern metadata like stencil neighborhood
    for i, pattern in enumerate(parfor.patterns):
        if pattern[0] == 'stencil':
            parfor.patterns[i] = ('stencil',
                replace_vars_inner(pattern[1], var_dict))

    # replace loop boundary variables
    for l in parfor.loop_nests:
        l.start = replace_vars_inner(l.start, var_dict)
        l.stop = replace_vars_inner(l.stop, var_dict)
        l.step = replace_vars_inner(l.step, var_dict)

    blocks = wrap_parfor_blocks(parfor)
    # add dummy assigns for each copy
    assign_list = []
    for lhs_name, rhs in var_dict.items():
        assign_list.append(ir.Assign(rhs, name_var_table[lhs_name],
                                     ir.Loc("dummy", -1)))
    blocks[0].body = assign_list + blocks[0].body
    in_copies_parfor, out_copies_parfor = copy_propagate(blocks, typemap)
    apply_copy_propagate(blocks, in_copies_parfor, name_var_table, typemap,
                         calltypes, save_copies)
    unwrap_parfor_blocks(parfor)
    # remove dummy assignments
    blocks[0].body = blocks[0].body[len(assign_list):]
    return


ir_utils.apply_copy_propagate_extensions[Parfor] = apply_copies_parfor


def push_call_vars(blocks, saved_globals, saved_getattrs, typemap, nested=False):
    """push call variables to right before their call site.
    assuming one global/getattr is created for each call site and control flow
    doesn't change it.
    """
    for block in blocks.values():
        new_body = []
        # global/attr variables that are defined in this block already,
        #   no need to reassign them
        block_defs = set()
        # Some definitions are copied right before the call but then we
        # need to rename that symbol in that block so that typing won't
        # generate an error trying to lock the save var twice.
        # In rename_dict, we collect the symbols that must be renamed in
        # this block.  We collect them then apply the renaming at the end.
        rename_dict = {}
        for stmt in block.body:
            def process_assign(stmt):
                if isinstance(stmt, ir.Assign):
                    rhs = stmt.value
                    lhs = stmt.target
                    if (isinstance(rhs, ir.Global)):
                        saved_globals[lhs.name] = stmt
                        block_defs.add(lhs.name)
                    elif isinstance(rhs, ir.Expr) and rhs.op == 'getattr':
                        if (rhs.value.name in saved_globals
                                or rhs.value.name in saved_getattrs):
                            saved_getattrs[lhs.name] = stmt
                            block_defs.add(lhs.name)

            if not nested and isinstance(stmt, Parfor):
                for s in stmt.init_block.body:
                    process_assign(s)
                pblocks = stmt.loop_body.copy()
                push_call_vars(pblocks, saved_globals, saved_getattrs, typemap, nested=True)
                new_body.append(stmt)
                continue
            else:
                process_assign(stmt)
            for v in stmt.list_vars():
                new_body += _get_saved_call_nodes(v.name, saved_globals,
                                                  saved_getattrs, block_defs, rename_dict)
            new_body.append(stmt)
        block.body = new_body
        # If there is anything to rename then apply the renaming here.
        if len(rename_dict) > 0:
            # Fix-up the typing for the renamed vars.
            for k, v in rename_dict.items():
                typemap[v] = typemap[k]
            # This is only to call replace_var_names which takes a dict.
            temp_blocks = {0: block}
            replace_var_names(temp_blocks, rename_dict)

    return


def _get_saved_call_nodes(fname, saved_globals, saved_getattrs, block_defs, rename_dict):
    """ Implement the copying of globals or getattrs for the purposes noted in
        push_call_vars.  We make a new var and assign to it a copy of the
        global or getattr.  We remember this new assignment node and add an
        entry in the renaming dictionary so that for this block the original
        var name is replaced by the new var name we created.
    """
    nodes = []
    while (fname not in block_defs and (fname in saved_globals
                                        or fname in saved_getattrs)):
        def rename_global_or_getattr(obj, var_base, nodes, block_defs, rename_dict):
            assert(isinstance(obj, ir.Assign))
            renamed_var = ir.Var(obj.target.scope,
                                 mk_unique_var(var_base),
                                 obj.target.loc)
            renamed_assign = ir.Assign(copy.deepcopy(obj.value),
                                       renamed_var,
                                       obj.loc)
            nodes.append(renamed_assign)
            block_defs.add(obj.target.name)
            rename_dict[obj.target.name] = renamed_assign.target.name

        if fname in saved_globals:
            rename_global_or_getattr(saved_globals[fname], "$push_global_to_block",
                                     nodes, block_defs, rename_dict)
            fname = '_PA_DONE'
        elif fname in saved_getattrs:
            rename_global_or_getattr(saved_getattrs[fname], "$push_getattr_to_block",
                                     nodes, block_defs, rename_dict)
            fname = saved_getattrs[fname].value.value.name
    nodes.reverse()
    return nodes

def repr_arrayexpr(arrayexpr):
    """Extract operators from arrayexpr to represent it abstractly as a string.
    """
    if isinstance(arrayexpr, tuple):
        opr = arrayexpr[0]
        # sometimes opr is not string like '+', but is a ufunc object
        if not isinstance(opr, str):
            if hasattr(opr, '__name__'):
                opr = opr.__name__
            else:
                opr = '_'  # can return dummy since repr is not critical
        args = arrayexpr[1]
        if len(args) == 1:
            return '({}({}))'.format(opr, repr_arrayexpr(args[0]))
        else:
            opr = ' ' + opr + ' '
            return '({})'.format(opr.join([ repr_arrayexpr(x) for x in args ]))
    elif isinstance(arrayexpr, numba.core.ir.Var):
        name = arrayexpr.name
        if name.startswith('$'):
            return '\'%s\' (temporary variable)' % name
        else:
            return name
    elif isinstance(arrayexpr, numba.core.ir.Const):
        return repr(arrayexpr.value)
    else:
        return '_'

def fix_generator_types(generator_info, return_type, typemap):
    """postproc updates generator_info with live variables after transformations
    but generator variables have types in return_type that are updated here.
    """
    new_state_types = []
    for v in generator_info.state_vars:
        new_state_types.append(typemap[v])
    return_type.state_types = tuple(new_state_types)
    return


def get_parfor_call_table(parfor, call_table=None, reverse_call_table=None):
    if call_table is None:
        call_table = {}
    if reverse_call_table is None:
        reverse_call_table = {}
    blocks = wrap_parfor_blocks(parfor)
    call_table, reverse_call_table = get_call_table(blocks, call_table,
                                                    reverse_call_table)
    unwrap_parfor_blocks(parfor)
    return call_table, reverse_call_table


ir_utils.call_table_extensions[Parfor] = get_parfor_call_table


def get_parfor_tuple_table(parfor, tuple_table=None):
    if tuple_table is None:
        tuple_table = {}
    blocks = wrap_parfor_blocks(parfor)
    tuple_table = ir_utils.get_tuple_table(blocks, tuple_table)
    unwrap_parfor_blocks(parfor)
    return tuple_table


ir_utils.tuple_table_extensions[Parfor] = get_parfor_tuple_table


def get_parfor_array_accesses(parfor, accesses=None):
    if accesses is None:
        accesses = set()
    blocks = wrap_parfor_blocks(parfor)
    accesses = ir_utils.get_array_accesses(blocks, accesses)
    unwrap_parfor_blocks(parfor)
    return accesses


# parfor handler is same as
ir_utils.array_accesses_extensions[Parfor] = get_parfor_array_accesses


def parfor_add_offset_to_labels(parfor, offset):
    blocks = wrap_parfor_blocks(parfor)
    blocks = add_offset_to_labels(blocks, offset)
    blocks[0] = blocks[offset]
    blocks.pop(offset)
    unwrap_parfor_blocks(parfor, blocks)
    return


ir_utils.add_offset_to_labels_extensions[Parfor] = parfor_add_offset_to_labels


def parfor_find_max_label(parfor):
    blocks = wrap_parfor_blocks(parfor)
    max_label = ir_utils.find_max_label(blocks)
    unwrap_parfor_blocks(parfor)
    return max_label

ir_utils.find_max_label_extensions[Parfor] = parfor_find_max_label


def parfor_typeinfer(parfor, typeinferer):
    save_blocks = typeinferer.blocks
    blocks = wrap_parfor_blocks(parfor)
    index_vars = [l.index_variable for l in parfor.loop_nests]
    # no need to handle parfor.index_var (tuple of variables), since it will be
    # assigned to a tuple from individual indices
    first_block = min(blocks.keys())
    loc = blocks[first_block].loc
    # XXX
    index_assigns = [ir.Assign(ir.Const(1, loc=loc, use_literal_type=False), v, loc) for v in index_vars]
    save_first_block_body = blocks[first_block].body
    blocks[first_block].body = index_assigns + blocks[first_block].body
    typeinferer.blocks = blocks
    typeinferer.build_constraint()
    typeinferer.blocks = save_blocks
    blocks[first_block].body = save_first_block_body
    unwrap_parfor_blocks(parfor)


typeinfer.typeinfer_extensions[Parfor] = parfor_typeinfer

def build_parfor_definitions(parfor, definitions=None):
    """get variable definition table for parfors"""
    if definitions is None:
        definitions = defaultdict(list)

    # avoid wrap_parfor_blocks() since build_definitions is called inside
    # find_potential_aliases_parfor where the parfor is already wrapped
    build_definitions(parfor.loop_body, definitions)
    build_definitions({0: parfor.init_block}, definitions)
    return definitions

ir_utils.build_defs_extensions[Parfor] = build_parfor_definitions

@contextmanager
def dummy_return_in_loop_body(loop_body):
    """adds dummy return to last block of parfor loop body for CFG computation
    """
    # max is last block since we add it manually for prange
    last_label = max(loop_body.keys())
    scope = loop_body[last_label].scope
    const = ir.Var(scope, mk_unique_var("$const"), ir.Loc("parfors_dummy", -1))
    loop_body[last_label].body.append(
        ir.Return(const, ir.Loc("parfors_dummy", -1)))
    yield
    # remove dummy return
    loop_body[last_label].body.pop()

@infer_global(reduce)
class ReduceInfer(AbstractTemplate):
    def generic(self, args, kws):
        assert not kws
        if len(args) != 3:
            raise errors.NumbaAssertionError("len(args) != 3")
        assert isinstance(args[1], types.Array)
        return signature(args[1].dtype, *args)


def ensure_parallel_support():
    """Check if the platform supports parallel=True and raise if it does not.
    """
    if config.IS_32BITS:
        msg = ("The 'parallel' target is not currently supported on 32 bit "
               "hardware.")
        raise errors.UnsupportedParforsError(msg)
