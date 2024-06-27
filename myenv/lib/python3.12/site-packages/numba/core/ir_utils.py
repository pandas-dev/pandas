#
# Copyright (c) 2017 Intel Corporation
# SPDX-License-Identifier: BSD-2-Clause
#

import numpy
import math

import types as pytypes
import collections
import warnings

import numba
from numba.core.extending import _Intrinsic
from numba.core import types, typing, ir, analysis, postproc, rewrites, config
from numba.core.typing.templates import signature
from numba.core.analysis import (compute_live_map, compute_use_defs,
                            compute_cfg_from_blocks)
from numba.core.errors import (TypingError, UnsupportedError,
                               NumbaPendingDeprecationWarning,
                               CompilerError)

import copy

_unique_var_count = 0


def mk_unique_var(prefix):
    global _unique_var_count
    var = prefix + "." + str(_unique_var_count)
    _unique_var_count = _unique_var_count + 1
    return var


class _MaxLabel:
    def __init__(self, value=0):
        self._value = value

    def next(self):
        self._value += 1
        return self._value

    def update(self, newval):
        self._value = max(newval, self._value)


_the_max_label = _MaxLabel()
del _MaxLabel


def get_unused_var_name(prefix, var_table):
    """ Get a new var name with a given prefix and
        make sure it is unused in the given variable table.
    """
    cur = 0
    while True:
        var = prefix + str(cur)
        if var not in var_table:
            return var
        cur += 1


def next_label():
    return _the_max_label.next()


def mk_alloc(typingctx, typemap, calltypes, lhs, size_var, dtype, scope, loc,
             lhs_typ):
    """generate an array allocation with np.empty() and return list of nodes.
    size_var can be an int variable or tuple of int variables.
    lhs_typ is the type of the array being allocated.
    """
    out = []
    ndims = 1
    size_typ = types.intp
    if isinstance(size_var, tuple):
        if len(size_var) == 1:
            size_var = size_var[0]
            size_var = convert_size_to_var(size_var, typemap, scope, loc, out)
        else:
            # tuple_var = build_tuple([size_var...])
            ndims = len(size_var)
            tuple_var = ir.Var(scope, mk_unique_var("$tuple_var"), loc)
            if typemap:
                typemap[tuple_var.name] = types.containers.UniTuple(
                    types.intp, ndims)
            # constant sizes need to be assigned to vars
            new_sizes = [convert_size_to_var(s, typemap, scope, loc, out)
                         for s in size_var]
            tuple_call = ir.Expr.build_tuple(new_sizes, loc)
            tuple_assign = ir.Assign(tuple_call, tuple_var, loc)
            out.append(tuple_assign)
            size_var = tuple_var
            size_typ = types.containers.UniTuple(types.intp, ndims)
    if hasattr(lhs_typ, "__allocate__"):
        return lhs_typ.__allocate__(
            typingctx,
            typemap,
            calltypes,
            lhs,
            size_var,
            dtype,
            scope,
            loc,
            lhs_typ,
            size_typ,
            out,
        )
    # g_np_var = Global(numpy)
    g_np_var = ir.Var(scope, mk_unique_var("$np_g_var"), loc)
    if typemap:
        typemap[g_np_var.name] = types.misc.Module(numpy)
    g_np = ir.Global('np', numpy, loc)
    g_np_assign = ir.Assign(g_np, g_np_var, loc)
    # attr call: empty_attr = getattr(g_np_var, empty)
    empty_attr_call = ir.Expr.getattr(g_np_var, "empty", loc)
    attr_var = ir.Var(scope, mk_unique_var("$empty_attr_attr"), loc)
    if typemap:
        typemap[attr_var.name] = get_np_ufunc_typ(numpy.empty)
    attr_assign = ir.Assign(empty_attr_call, attr_var, loc)
     # Assume str(dtype) returns a valid type
    dtype_str = str(dtype)
    # alloc call: lhs = empty_attr(size_var, typ_var)
    typ_var = ir.Var(scope, mk_unique_var("$np_typ_var"), loc)
    if typemap:
        typemap[typ_var.name] = types.functions.NumberClass(dtype)
    # If dtype is a datetime/timedelta with a unit,
    # then it won't return a valid type and instead can be created
    # with a string. i.e. "datetime64[ns]")
    if (isinstance(dtype, (types.NPDatetime, types.NPTimedelta)) and
        dtype.unit != ''):
            typename_const = ir.Const(dtype_str, loc)
            typ_var_assign = ir.Assign(typename_const, typ_var, loc)
    else:
        if dtype_str=='bool':
            # empty doesn't like 'bool' sometimes (e.g. kmeans example)
            dtype_str = 'bool_'
        np_typ_getattr = ir.Expr.getattr(g_np_var, dtype_str, loc)
        typ_var_assign = ir.Assign(np_typ_getattr, typ_var, loc)
    alloc_call = ir.Expr.call(attr_var, [size_var, typ_var], (), loc)

    if calltypes:
        cac = typemap[attr_var.name].get_call_type(
            typingctx, [size_typ, types.functions.NumberClass(dtype)], {})
        # By default, all calls to "empty" are typed as returning a standard
        # NumPy ndarray.  If we are allocating a ndarray subclass here then
        # just change the return type to be that of the subclass.
        cac._return_type = (lhs_typ.copy(layout='C')
                            if lhs_typ.layout == 'F'
                            else lhs_typ)
        calltypes[alloc_call] = cac
    if lhs_typ.layout == 'F':
        empty_c_typ = lhs_typ.copy(layout='C')
        empty_c_var = ir.Var(scope, mk_unique_var("$empty_c_var"), loc)
        if typemap:
            typemap[empty_c_var.name] = lhs_typ.copy(layout='C')
        empty_c_assign = ir.Assign(alloc_call, empty_c_var, loc)

        # attr call: asfortranarray = getattr(g_np_var, asfortranarray)
        asfortranarray_attr_call = ir.Expr.getattr(g_np_var, "asfortranarray", loc)
        afa_attr_var = ir.Var(scope, mk_unique_var("$asfortran_array_attr"), loc)
        if typemap:
            typemap[afa_attr_var.name] = get_np_ufunc_typ(numpy.asfortranarray)
        afa_attr_assign = ir.Assign(asfortranarray_attr_call, afa_attr_var, loc)
        # call asfortranarray
        asfortranarray_call = ir.Expr.call(afa_attr_var, [empty_c_var], (), loc)
        if calltypes:
            calltypes[asfortranarray_call] = typemap[afa_attr_var.name].get_call_type(
                typingctx, [empty_c_typ], {})

        asfortranarray_assign = ir.Assign(asfortranarray_call, lhs, loc)

        out.extend([g_np_assign, attr_assign, typ_var_assign, empty_c_assign,
                    afa_attr_assign, asfortranarray_assign])
    else:
        alloc_assign = ir.Assign(alloc_call, lhs, loc)
        out.extend([g_np_assign, attr_assign, typ_var_assign, alloc_assign])

    return out


def convert_size_to_var(size_var, typemap, scope, loc, nodes):
    if isinstance(size_var, int):
        new_size = ir.Var(scope, mk_unique_var("$alloc_size"), loc)
        if typemap:
            typemap[new_size.name] = types.intp
        size_assign = ir.Assign(ir.Const(size_var, loc), new_size, loc)
        nodes.append(size_assign)
        return new_size
    assert isinstance(size_var, ir.Var)
    return size_var


def get_np_ufunc_typ(func):
    """get type of the incoming function from builtin registry"""
    for (k, v) in typing.npydecl.registry.globals:
        if k == func:
            return v
    for (k, v) in typing.templates.builtin_registry.globals:
        if k == func:
            return v
    raise RuntimeError("type for func ", func, " not found")


def mk_range_block(typemap, start, stop, step, calltypes, scope, loc):
    """make a block that initializes loop range and iteration variables.
    target label in jump needs to be set.
    """
    # g_range_var = Global(range)
    g_range_var = ir.Var(scope, mk_unique_var("$range_g_var"), loc)
    typemap[g_range_var.name] = get_global_func_typ(range)
    g_range = ir.Global('range', range, loc)
    g_range_assign = ir.Assign(g_range, g_range_var, loc)
    arg_nodes, args = _mk_range_args(typemap, start, stop, step, scope, loc)
    # range_call_var = call g_range_var(start, stop, step)
    range_call = ir.Expr.call(g_range_var, args, (), loc)
    calltypes[range_call] = typemap[g_range_var.name].get_call_type(
        typing.Context(), [types.intp] * len(args), {})
    #signature(types.range_state64_type, types.intp)
    range_call_var = ir.Var(scope, mk_unique_var("$range_c_var"), loc)
    typemap[range_call_var.name] = types.iterators.RangeType(types.intp)
    range_call_assign = ir.Assign(range_call, range_call_var, loc)
    # iter_var = getiter(range_call_var)
    iter_call = ir.Expr.getiter(range_call_var, loc)
    calltypes[iter_call] = signature(types.range_iter64_type,
                                     types.range_state64_type)
    iter_var = ir.Var(scope, mk_unique_var("$iter_var"), loc)
    typemap[iter_var.name] = types.iterators.RangeIteratorType(types.intp)
    iter_call_assign = ir.Assign(iter_call, iter_var, loc)
    # $phi = iter_var
    phi_var = ir.Var(scope, mk_unique_var("$phi"), loc)
    typemap[phi_var.name] = types.iterators.RangeIteratorType(types.intp)
    phi_assign = ir.Assign(iter_var, phi_var, loc)
    # jump to header
    jump_header = ir.Jump(-1, loc)
    range_block = ir.Block(scope, loc)
    range_block.body = arg_nodes + [g_range_assign, range_call_assign,
                                    iter_call_assign, phi_assign, jump_header]
    return range_block


def _mk_range_args(typemap, start, stop, step, scope, loc):
    nodes = []
    if isinstance(stop, ir.Var):
        g_stop_var = stop
    else:
        assert isinstance(stop, int)
        g_stop_var = ir.Var(scope, mk_unique_var("$range_stop"), loc)
        if typemap:
            typemap[g_stop_var.name] = types.intp
        stop_assign = ir.Assign(ir.Const(stop, loc), g_stop_var, loc)
        nodes.append(stop_assign)
    if start == 0 and step == 1:
        return nodes, [g_stop_var]

    if isinstance(start, ir.Var):
        g_start_var = start
    else:
        assert isinstance(start, int)
        g_start_var = ir.Var(scope, mk_unique_var("$range_start"), loc)
        if typemap:
            typemap[g_start_var.name] = types.intp
        start_assign = ir.Assign(ir.Const(start, loc), g_start_var, loc)
        nodes.append(start_assign)
    if step == 1:
        return nodes, [g_start_var, g_stop_var]

    if isinstance(step, ir.Var):
        g_step_var = step
    else:
        assert isinstance(step, int)
        g_step_var = ir.Var(scope, mk_unique_var("$range_step"), loc)
        if typemap:
            typemap[g_step_var.name] = types.intp
        step_assign = ir.Assign(ir.Const(step, loc), g_step_var, loc)
        nodes.append(step_assign)

    return nodes, [g_start_var, g_stop_var, g_step_var]


def get_global_func_typ(func):
    """get type variable for func() from builtin registry"""
    for (k, v) in typing.templates.builtin_registry.globals:
        if k == func:
            return v
    raise RuntimeError("func type not found {}".format(func))


def mk_loop_header(typemap, phi_var, calltypes, scope, loc):
    """make a block that is a loop header updating iteration variables.
    target labels in branch need to be set.
    """
    # iternext_var = iternext(phi_var)
    iternext_var = ir.Var(scope, mk_unique_var("$iternext_var"), loc)
    typemap[iternext_var.name] = types.containers.Pair(
        types.intp, types.boolean)
    iternext_call = ir.Expr.iternext(phi_var, loc)
    calltypes[iternext_call] = signature(
        types.containers.Pair(
            types.intp,
            types.boolean),
        types.range_iter64_type)
    iternext_assign = ir.Assign(iternext_call, iternext_var, loc)
    # pair_first_var = pair_first(iternext_var)
    pair_first_var = ir.Var(scope, mk_unique_var("$pair_first_var"), loc)
    typemap[pair_first_var.name] = types.intp
    pair_first_call = ir.Expr.pair_first(iternext_var, loc)
    pair_first_assign = ir.Assign(pair_first_call, pair_first_var, loc)
    # pair_second_var = pair_second(iternext_var)
    pair_second_var = ir.Var(scope, mk_unique_var("$pair_second_var"), loc)
    typemap[pair_second_var.name] = types.boolean
    pair_second_call = ir.Expr.pair_second(iternext_var, loc)
    pair_second_assign = ir.Assign(pair_second_call, pair_second_var, loc)
    # phi_b_var = pair_first_var
    phi_b_var = ir.Var(scope, mk_unique_var("$phi"), loc)
    typemap[phi_b_var.name] = types.intp
    phi_b_assign = ir.Assign(pair_first_var, phi_b_var, loc)
    # branch pair_second_var body_block out_block
    branch = ir.Branch(pair_second_var, -1, -1, loc)
    header_block = ir.Block(scope, loc)
    header_block.body = [iternext_assign, pair_first_assign,
                         pair_second_assign, phi_b_assign, branch]
    return header_block


def legalize_names(varnames):
    """returns a dictionary for conversion of variable names to legal
    parameter names.
    """
    var_map = {}
    for var in varnames:
        new_name = var.replace("_", "__").replace("$", "_").replace(".", "_")
        assert new_name not in var_map
        var_map[var] = new_name
    return var_map


def get_name_var_table(blocks):
    """create a mapping from variable names to their ir.Var objects"""
    def get_name_var_visit(var, namevar):
        namevar[var.name] = var
        return var
    namevar = {}
    visit_vars(blocks, get_name_var_visit, namevar)
    return namevar


def replace_var_names(blocks, namedict):
    """replace variables (ir.Var to ir.Var) from dictionary (name -> name)"""
    # remove identity values to avoid infinite loop
    new_namedict = {}
    for l, r in namedict.items():
        if l != r:
            new_namedict[l] = r

    def replace_name(var, namedict):
        assert isinstance(var, ir.Var)
        while var.name in namedict:
            var = ir.Var(var.scope, namedict[var.name], var.loc)
        return var
    visit_vars(blocks, replace_name, new_namedict)


def replace_var_callback(var, vardict):
    assert isinstance(var, ir.Var)
    while var.name in vardict.keys():
        assert(vardict[var.name].name != var.name)
        new_var = vardict[var.name]
        var = ir.Var(new_var.scope, new_var.name, new_var.loc)
    return var


def replace_vars(blocks, vardict):
    """replace variables (ir.Var to ir.Var) from dictionary (name -> ir.Var)"""
    # remove identity values to avoid infinite loop
    new_vardict = {}
    for l, r in vardict.items():
        if l != r.name:
            new_vardict[l] = r
    visit_vars(blocks, replace_var_callback, new_vardict)


def replace_vars_stmt(stmt, vardict):
    visit_vars_stmt(stmt, replace_var_callback, vardict)


def replace_vars_inner(node, vardict):
    return visit_vars_inner(node, replace_var_callback, vardict)


# other packages that define new nodes add calls to visit variables in them
# format: {type:function}
visit_vars_extensions = {}


def visit_vars(blocks, callback, cbdata):
    """go over statements of block bodies and replace variable names with
    dictionary.
    """
    for block in blocks.values():
        for stmt in block.body:
            visit_vars_stmt(stmt, callback, cbdata)
    return


def visit_vars_stmt(stmt, callback, cbdata):
    # let external calls handle stmt if type matches
    for t, f in visit_vars_extensions.items():
        if isinstance(stmt, t):
            f(stmt, callback, cbdata)
            return
    if isinstance(stmt, ir.Assign):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
        stmt.value = visit_vars_inner(stmt.value, callback, cbdata)
    elif isinstance(stmt, ir.Arg):
        stmt.name = visit_vars_inner(stmt.name, callback, cbdata)
    elif isinstance(stmt, ir.Return):
        stmt.value = visit_vars_inner(stmt.value, callback, cbdata)
    elif isinstance(stmt, ir.Raise):
        stmt.exception = visit_vars_inner(stmt.exception, callback, cbdata)
    elif isinstance(stmt, ir.Branch):
        stmt.cond = visit_vars_inner(stmt.cond, callback, cbdata)
    elif isinstance(stmt, ir.Jump):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
    elif isinstance(stmt, ir.Del):
        # Because Del takes only a var name, we make up by
        # constructing a temporary variable.
        var = ir.Var(None, stmt.value, stmt.loc)
        var = visit_vars_inner(var, callback, cbdata)
        stmt.value = var.name
    elif isinstance(stmt, ir.DelAttr):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
        stmt.attr = visit_vars_inner(stmt.attr, callback, cbdata)
    elif isinstance(stmt, ir.SetAttr):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
        stmt.attr = visit_vars_inner(stmt.attr, callback, cbdata)
        stmt.value = visit_vars_inner(stmt.value, callback, cbdata)
    elif isinstance(stmt, ir.DelItem):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
        stmt.index = visit_vars_inner(stmt.index, callback, cbdata)
    elif isinstance(stmt, ir.StaticSetItem):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
        stmt.index_var = visit_vars_inner(stmt.index_var, callback, cbdata)
        stmt.value = visit_vars_inner(stmt.value, callback, cbdata)
    elif isinstance(stmt, ir.SetItem):
        stmt.target = visit_vars_inner(stmt.target, callback, cbdata)
        stmt.index = visit_vars_inner(stmt.index, callback, cbdata)
        stmt.value = visit_vars_inner(stmt.value, callback, cbdata)
    elif isinstance(stmt, ir.Print):
        stmt.args = [visit_vars_inner(x, callback, cbdata) for x in stmt.args]
    else:
        # TODO: raise NotImplementedError("no replacement for IR node: ", stmt)
        pass
    return


def visit_vars_inner(node, callback, cbdata):
    if isinstance(node, ir.Var):
        return callback(node, cbdata)
    elif isinstance(node, list):
        return [visit_vars_inner(n, callback, cbdata) for n in node]
    elif isinstance(node, tuple):
        return tuple([visit_vars_inner(n, callback, cbdata) for n in node])
    elif isinstance(node, ir.Expr):
        # if node.op in ['binop', 'inplace_binop']:
        #     lhs = node.lhs.name
        #     rhs = node.rhs.name
        #     node.lhs.name = callback, cbdata.get(lhs, lhs)
        #     node.rhs.name = callback, cbdata.get(rhs, rhs)
        for arg in node._kws.keys():
            node._kws[arg] = visit_vars_inner(node._kws[arg], callback, cbdata)
    elif isinstance(node, ir.Yield):
        node.value = visit_vars_inner(node.value, callback, cbdata)
    return node


add_offset_to_labels_extensions = {}


def add_offset_to_labels(blocks, offset):
    """add an offset to all block labels and jump/branch targets
    """
    new_blocks = {}
    for l, b in blocks.items():
        # some parfor last blocks might be empty
        term = None
        if b.body:
            term = b.body[-1]
            for inst in b.body:
                for T, f in add_offset_to_labels_extensions.items():
                    if isinstance(inst, T):
                        f_max = f(inst, offset)
        if isinstance(term, ir.Jump):
            b.body[-1] = ir.Jump(term.target + offset, term.loc)
        if isinstance(term, ir.Branch):
            b.body[-1] = ir.Branch(term.cond, term.truebr + offset,
                                   term.falsebr + offset, term.loc)
        new_blocks[l + offset] = b
    return new_blocks


find_max_label_extensions = {}


def find_max_label(blocks):
    max_label = 0
    for l, b in blocks.items():
        term = None
        if b.body:
            term = b.body[-1]
            for inst in b.body:
                for T, f in find_max_label_extensions.items():
                    if isinstance(inst, T):
                        f_max = f(inst)
                        if f_max > max_label:
                            max_label = f_max
        if l > max_label:
            max_label = l
    return max_label


def flatten_labels(blocks):
    """makes the labels in range(0, len(blocks)), useful to compare CFGs
    """
    # first bulk move the labels out of the rewrite range
    blocks = add_offset_to_labels(blocks, find_max_label(blocks) + 1)
    # order them in topo order because it's easier to read
    new_blocks = {}
    topo_order = find_topo_order(blocks)
    l_map = dict()
    idx = 0
    for x in topo_order:
        l_map[x] = idx
        idx += 1

    for t_node in topo_order:
        b = blocks[t_node]
        # some parfor last blocks might be empty
        term = None
        if b.body:
            term = b.body[-1]
        if isinstance(term, ir.Jump):
            b.body[-1] = ir.Jump(l_map[term.target], term.loc)
        if isinstance(term, ir.Branch):
            b.body[-1] = ir.Branch(term.cond, l_map[term.truebr],
                                   l_map[term.falsebr], term.loc)
        new_blocks[l_map[t_node]] = b
    return new_blocks


def remove_dels(blocks):
    """remove ir.Del nodes"""
    for block in blocks.values():
        new_body = []
        for stmt in block.body:
            if not isinstance(stmt, ir.Del):
                new_body.append(stmt)
        block.body = new_body
    return


def remove_args(blocks):
    """remove ir.Arg nodes"""
    for block in blocks.values():
        new_body = []
        for stmt in block.body:
            if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Arg):
                continue
            new_body.append(stmt)
        block.body = new_body
    return


def dead_code_elimination(func_ir, typemap=None, alias_map=None,
                          arg_aliases=None):
    """ Performs dead code elimination and leaves the IR in a valid state on
    exit
    """
    do_post_proc = False
    while (remove_dead(func_ir.blocks, func_ir.arg_names, func_ir, typemap,
                       alias_map, arg_aliases)):
        do_post_proc = True

    if do_post_proc:
        post_proc = postproc.PostProcessor(func_ir)
        post_proc.run()


def remove_dead(blocks, args, func_ir, typemap=None, alias_map=None, arg_aliases=None):
    """dead code elimination using liveness and CFG info.
    Returns True if something has been removed, or False if nothing is removed.
    """
    cfg = compute_cfg_from_blocks(blocks)
    usedefs = compute_use_defs(blocks)
    live_map = compute_live_map(cfg, blocks, usedefs.usemap, usedefs.defmap)
    call_table, _ = get_call_table(blocks)
    if alias_map is None or arg_aliases is None:
        alias_map, arg_aliases = find_potential_aliases(blocks, args, typemap,
                                                        func_ir)
    if config.DEBUG_ARRAY_OPT >= 1:
        print("args:", args)
        print("alias map:", alias_map)
        print("arg_aliases:", arg_aliases)
        print("live_map:", live_map)
        print("usemap:", usedefs.usemap)
        print("defmap:", usedefs.defmap)
    # keep set for easier search
    alias_set = set(alias_map.keys())

    removed = False
    for label, block in blocks.items():
        # find live variables at each statement to delete dead assignment
        lives = {v.name for v in block.terminator.list_vars()}
        if config.DEBUG_ARRAY_OPT >= 2:
            print("remove_dead processing block", label, lives)
        # find live variables at the end of block
        for out_blk, _data in cfg.successors(label):
            if config.DEBUG_ARRAY_OPT >= 2:
                print("succ live_map", out_blk, live_map[out_blk])
            lives |= live_map[out_blk]
        removed |= remove_dead_block(block, lives, call_table, arg_aliases,
                                     alias_map, alias_set, func_ir, typemap)

    return removed


# other packages that define new nodes add calls to remove dead code in them
# format: {type:function}
remove_dead_extensions = {}


def remove_dead_block(block, lives, call_table, arg_aliases, alias_map,
                                                  alias_set, func_ir, typemap):
    """remove dead code using liveness info.
    Mutable arguments (e.g. arrays) that are not definitely assigned are live
    after return of function.
    """
    # TODO: find mutable args that are not definitely assigned instead of
    # assuming all args are live after return
    removed = False

    # add statements in reverse order
    new_body = [block.terminator]
    # for each statement in reverse order, excluding terminator
    for stmt in reversed(block.body[:-1]):
        if config.DEBUG_ARRAY_OPT >= 2:
            print("remove_dead_block", stmt)
        # aliases of lives are also live
        alias_lives = set()
        init_alias_lives = lives & alias_set
        for v in init_alias_lives:
            alias_lives |= alias_map[v]
        lives_n_aliases = lives | alias_lives | arg_aliases

        # let external calls handle stmt if type matches
        if type(stmt) in remove_dead_extensions:
            f = remove_dead_extensions[type(stmt)]
            stmt = f(stmt, lives, lives_n_aliases, arg_aliases, alias_map, func_ir,
                     typemap)
            if stmt is None:
                if config.DEBUG_ARRAY_OPT >= 2:
                    print("Statement was removed.")
                removed = True
                continue

        # ignore assignments that their lhs is not live or lhs==rhs
        if isinstance(stmt, ir.Assign):
            lhs = stmt.target
            rhs = stmt.value
            if lhs.name not in lives and has_no_side_effect(
                    rhs, lives_n_aliases, call_table):
                if config.DEBUG_ARRAY_OPT >= 2:
                    print("Statement was removed.")
                removed = True
                continue
            if isinstance(rhs, ir.Var) and lhs.name == rhs.name:
                if config.DEBUG_ARRAY_OPT >= 2:
                    print("Statement was removed.")
                removed = True
                continue
            # TODO: remove other nodes like SetItem etc.

        if isinstance(stmt, ir.Del):
            if stmt.value not in lives:
                if config.DEBUG_ARRAY_OPT >= 2:
                    print("Statement was removed.")
                removed = True
                continue

        if isinstance(stmt, ir.SetItem):
            name = stmt.target.name
            if name not in lives_n_aliases:
                if config.DEBUG_ARRAY_OPT >= 2:
                    print("Statement was removed.")
                continue

        if type(stmt) in analysis.ir_extension_usedefs:
            def_func = analysis.ir_extension_usedefs[type(stmt)]
            uses, defs = def_func(stmt)
            lives -= defs
            lives |= uses
        else:
            lives |= {v.name for v in stmt.list_vars()}
            if isinstance(stmt, ir.Assign):
                # make sure lhs is not used in rhs, e.g. a = g(a)
                if isinstance(stmt.value, ir.Expr):
                    rhs_vars = {v.name for v in stmt.value.list_vars()}
                    if lhs.name not in rhs_vars:
                        lives.remove(lhs.name)
                else:
                    lives.remove(lhs.name)

        new_body.append(stmt)
    new_body.reverse()
    block.body = new_body
    return removed

# list of functions
remove_call_handlers = []

def remove_dead_random_call(rhs, lives, call_list):
    if len(call_list) == 3 and call_list[1:] == ['random', numpy]:
        return call_list[0] not in {'seed', 'shuffle'}
    return False

remove_call_handlers.append(remove_dead_random_call)

def has_no_side_effect(rhs, lives, call_table):
    """ Returns True if this expression has no side effects that
        would prevent re-ordering.
    """
    from numba.parfors import array_analysis, parfor
    from numba.misc.special import prange
    if isinstance(rhs, ir.Expr) and rhs.op == 'call':
        func_name = rhs.func.name
        if func_name not in call_table or call_table[func_name] == []:
            return False
        call_list = call_table[func_name]
        if (call_list == ['empty', numpy] or
            call_list == [slice] or
            call_list == ['stencil', numba] or
            call_list == ['log', numpy] or
            call_list == ['dtype', numpy] or
            call_list == [array_analysis.wrap_index] or
            call_list == [prange] or
            call_list == ['prange', numba] or
            call_list == ['pndindex', numba] or
            call_list == [parfor.internal_prange] or
            call_list == ['ceil', math] or
            call_list == [max] or
            call_list == [int]):
            return True
        elif (isinstance(call_list[0], _Intrinsic) and
              (call_list[0]._name == 'empty_inferred' or
               call_list[0]._name == 'unsafe_empty_inferred')):
            return True
        from numba.core.registry import CPUDispatcher
        from numba.np.linalg import dot_3_mv_check_args
        if isinstance(call_list[0], CPUDispatcher):
            py_func = call_list[0].py_func
            if py_func == dot_3_mv_check_args:
                return True
        for f in remove_call_handlers:
            if f(rhs, lives, call_list):
                return True
        return False
    if isinstance(rhs, ir.Expr) and rhs.op == 'inplace_binop':
        return rhs.lhs.name not in lives
    if isinstance(rhs, ir.Yield):
        return False
    if isinstance(rhs, ir.Expr) and rhs.op == 'pair_first':
        # don't remove pair_first since prange looks for it
        return False
    return True

is_pure_extensions = []

def is_pure(rhs, lives, call_table):
    """ Returns True if every time this expression is evaluated it
        returns the same result.  This is not the case for things
        like calls to numpy.random.
    """
    if isinstance(rhs, ir.Expr):
        if rhs.op == 'call':
            func_name = rhs.func.name
            if func_name not in call_table or call_table[func_name] == []:
                return False
            call_list = call_table[func_name]
            if (call_list == [slice] or
                call_list == ['log', numpy] or
                call_list == ['empty', numpy] or
                call_list == ['ceil', math] or
                call_list == [max] or
                call_list == [int]):
                return True
            for f in is_pure_extensions:
                if f(rhs, lives, call_list):
                    return True
            return False
        elif rhs.op == 'getiter' or rhs.op == 'iternext':
            return False
    if isinstance(rhs, ir.Yield):
        return False
    return True

def is_const_call(module_name, func_name):
    # Returns True if there is no state in the given module changed by the given function.
    if module_name == 'numpy':
        if func_name in ['empty']:
            return True
    return False

alias_analysis_extensions = {}
alias_func_extensions = {}

def get_canonical_alias(v, alias_map):
    if v not in alias_map:
        return v

    v_aliases = sorted(list(alias_map[v]))
    return v_aliases[0]

def find_potential_aliases(blocks, args, typemap, func_ir, alias_map=None,
                                                           arg_aliases=None):
    "find all array aliases and argument aliases to avoid remove as dead"
    if alias_map is None:
        alias_map = {}
    if arg_aliases is None:
        arg_aliases = set(a for a in args if not is_immutable_type(a, typemap))

    # update definitions since they are not guaranteed to be up-to-date
    # FIXME keep definitions up-to-date to avoid the need for rebuilding
    func_ir._definitions = build_definitions(func_ir.blocks)
    np_alias_funcs = ['ravel', 'transpose', 'reshape']

    for bl in blocks.values():
        for instr in bl.body:
            if type(instr) in alias_analysis_extensions:
                f = alias_analysis_extensions[type(instr)]
                f(instr, args, typemap, func_ir, alias_map, arg_aliases)
            if isinstance(instr, ir.Assign):
                expr = instr.value
                lhs = instr.target.name
                # only mutable types can alias
                if is_immutable_type(lhs, typemap):
                    continue
                if isinstance(expr, ir.Var) and lhs!=expr.name:
                    _add_alias(lhs, expr.name, alias_map, arg_aliases)
                # subarrays like A = B[0] for 2D B
                if (isinstance(expr, ir.Expr) and (expr.op == 'cast' or
                    expr.op in ['getitem', 'static_getitem'])):
                    _add_alias(lhs, expr.value.name, alias_map, arg_aliases)
                if isinstance(expr, ir.Expr) and expr.op == 'inplace_binop':
                    _add_alias(lhs, expr.lhs.name, alias_map, arg_aliases)
                # array attributes like A.T
                if (isinstance(expr, ir.Expr) and expr.op == 'getattr'
                        and expr.attr in ['T', 'ctypes', 'flat']):
                    _add_alias(lhs, expr.value.name, alias_map, arg_aliases)
                # a = b.c.  a should alias b
                if (isinstance(expr, ir.Expr) and expr.op == 'getattr'
                        and expr.attr not in ['shape']
                        and expr.value.name in arg_aliases):
                    _add_alias(lhs, expr.value.name, alias_map, arg_aliases)
                # calls that can create aliases such as B = A.ravel()
                if isinstance(expr, ir.Expr) and expr.op == 'call':
                    fdef = guard(find_callname, func_ir, expr, typemap)
                    # TODO: sometimes gufunc backend creates duplicate code
                    # causing find_callname to fail. Example: test_argmax
                    # ignored here since those cases don't create aliases
                    # but should be fixed in general
                    if fdef is None:
                        continue
                    fname, fmod = fdef
                    if fdef in alias_func_extensions:
                        alias_func = alias_func_extensions[fdef]
                        alias_func(lhs, expr.args, alias_map, arg_aliases)
                    if fmod == 'numpy' and fname in np_alias_funcs:
                        _add_alias(lhs, expr.args[0].name, alias_map, arg_aliases)
                    if isinstance(fmod, ir.Var) and fname in np_alias_funcs:
                        _add_alias(lhs, fmod.name, alias_map, arg_aliases)

    # copy to avoid changing size during iteration
    old_alias_map = copy.deepcopy(alias_map)
    # combine all aliases transitively
    for v in old_alias_map:
        for w in old_alias_map[v]:
            alias_map[v] |= alias_map[w]
        for w in old_alias_map[v]:
            alias_map[w] = alias_map[v]

    return alias_map, arg_aliases

def _add_alias(lhs, rhs, alias_map, arg_aliases):
    if rhs in arg_aliases:
        arg_aliases.add(lhs)
    else:
        if rhs not in alias_map:
            alias_map[rhs] = set()
        if lhs not in alias_map:
            alias_map[lhs] = set()
        alias_map[rhs].add(lhs)
        alias_map[lhs].add(rhs)
    return

def is_immutable_type(var, typemap):
    # Conservatively, assume mutable if type not available
    if typemap is None or var not in typemap:
        return False
    typ = typemap[var]
    # TODO: add more immutable types
    if isinstance(typ, (types.Number, types.scalars._NPDatetimeBase,
                        types.iterators.RangeType)):
        return True
    if typ==types.string:
        return True
    # conservatively, assume mutable
    return False

def copy_propagate(blocks, typemap):
    """compute copy propagation information for each block using fixed-point
     iteration on data flow equations:
     in_b = intersect(predec(B))
     out_b = gen_b | (in_b - kill_b)
    """
    cfg = compute_cfg_from_blocks(blocks)
    entry = cfg.entry_point()

    # format: dict of block labels to copies as tuples
    # label -> (l,r)
    c_data = init_copy_propagate_data(blocks, entry, typemap)
    (gen_copies, all_copies, kill_copies, in_copies, out_copies) = c_data

    old_point = None
    new_point = copy.deepcopy(out_copies)
    # comparison works since dictionary of built-in types
    while old_point != new_point:
        for label in blocks.keys():
            if label == entry:
                continue
            predecs = [i for i, _d in cfg.predecessors(label)]
            # in_b =  intersect(predec(B))
            in_copies[label] = out_copies[predecs[0]].copy()
            for p in predecs:
                in_copies[label] &= out_copies[p]

            # out_b = gen_b | (in_b - kill_b)
            out_copies[label] = (gen_copies[label]
                                 | (in_copies[label] - kill_copies[label]))
        old_point = new_point
        new_point = copy.deepcopy(out_copies)
    if config.DEBUG_ARRAY_OPT >= 1:
        print("copy propagate out_copies:", out_copies)
    return in_copies, out_copies


def init_copy_propagate_data(blocks, entry, typemap):
    """get initial condition of copy propagation data flow for each block.
    """
    # gen is all definite copies, extra_kill is additional ones that may hit
    # for example, parfors can have control flow so they may hit extra copies
    gen_copies, extra_kill = get_block_copies(blocks, typemap)
    # set of all program copies
    all_copies = set()
    for l, s in gen_copies.items():
        all_copies |= gen_copies[l]
    kill_copies = {}
    for label, gen_set in gen_copies.items():
        kill_copies[label] = set()
        for lhs, rhs in all_copies:
            if lhs in extra_kill[label] or rhs in extra_kill[label]:
                kill_copies[label].add((lhs, rhs))
            # a copy is killed if it is not in this block and lhs or rhs are
            # assigned in this block
            assigned = {lhs for lhs, rhs in gen_set}
            if ((lhs, rhs) not in gen_set
                    and (lhs in assigned or rhs in assigned)):
                kill_copies[label].add((lhs, rhs))
    # set initial values
    # all copies are in for all blocks except entry
    in_copies = {l: all_copies.copy() for l in blocks.keys()}
    in_copies[entry] = set()
    out_copies = {}
    for label in blocks.keys():
        # out_b = gen_b | (in_b - kill_b)
        out_copies[label] = (gen_copies[label]
                             | (in_copies[label] - kill_copies[label]))
    out_copies[entry] = gen_copies[entry]
    return (gen_copies, all_copies, kill_copies, in_copies, out_copies)


# other packages that define new nodes add calls to get copies in them
# format: {type:function}
copy_propagate_extensions = {}


def get_block_copies(blocks, typemap):
    """get copies generated and killed by each block
    """
    block_copies = {}
    extra_kill = {}
    for label, block in blocks.items():
        assign_dict = {}
        extra_kill[label] = set()
        # assignments as dict to replace with latest value
        for stmt in block.body:
            for T, f in copy_propagate_extensions.items():
                if isinstance(stmt, T):
                    gen_set, kill_set = f(stmt, typemap)
                    for lhs, rhs in gen_set:
                        assign_dict[lhs] = rhs
                    # if a=b is in dict and b is killed, a is also killed
                    new_assign_dict = {}
                    for l, r in assign_dict.items():
                        if l not in kill_set and r not in kill_set:
                            new_assign_dict[l] = r
                        if r in kill_set:
                            extra_kill[label].add(l)
                    assign_dict = new_assign_dict
                    extra_kill[label] |= kill_set
            if isinstance(stmt, ir.Assign):
                lhs = stmt.target.name
                if isinstance(stmt.value, ir.Var):
                    rhs = stmt.value.name
                    # copy is valid only if same type (see
                    # TestCFunc.test_locals)
                    # Some transformations can produce assignments of the
                    # form A = A.  We don't put these mapping in the
                    # copy propagation set because then you get cycles and
                    # infinite loops in the replacement phase.
                    if typemap[lhs] == typemap[rhs] and lhs != rhs:
                        assign_dict[lhs] = rhs
                        continue
                if isinstance(stmt.value,
                              ir.Expr) and stmt.value.op == 'inplace_binop':
                    in1_var = stmt.value.lhs.name
                    in1_typ = typemap[in1_var]
                    # inplace_binop assigns first operand if mutable
                    if not (isinstance(in1_typ, types.Number)
                            or in1_typ == types.string):
                        extra_kill[label].add(in1_var)
                        # if a=b is in dict and b is killed, a is also killed
                        new_assign_dict = {}
                        for l, r in assign_dict.items():
                            if l != in1_var and r != in1_var:
                                new_assign_dict[l] = r
                            if r == in1_var:
                                extra_kill[label].add(l)
                        assign_dict = new_assign_dict
                extra_kill[label].add(lhs)
        block_cps = set(assign_dict.items())
        block_copies[label] = block_cps
    return block_copies, extra_kill


# other packages that define new nodes add calls to apply copy propagate in them
# format: {type:function}
apply_copy_propagate_extensions = {}


def apply_copy_propagate(blocks, in_copies, name_var_table, typemap, calltypes,
                         save_copies=None):
    """apply copy propagation to IR: replace variables when copies available"""
    # save_copies keeps an approximation of the copies that were applied, so
    # that the variable names of removed user variables can be recovered to some
    # extent.
    if save_copies is None:
        save_copies = []

    for label, block in blocks.items():
        var_dict = {l: name_var_table[r] for l, r in in_copies[label]}
        # assignments as dict to replace with latest value
        for stmt in block.body:
            if type(stmt) in apply_copy_propagate_extensions:
                f = apply_copy_propagate_extensions[type(stmt)]
                f(stmt, var_dict, name_var_table,
                    typemap, calltypes, save_copies)
            # only rhs of assignments should be replaced
            # e.g. if x=y is available, x in x=z shouldn't be replaced
            elif isinstance(stmt, ir.Assign):
                stmt.value = replace_vars_inner(stmt.value, var_dict)
            else:
                replace_vars_stmt(stmt, var_dict)
            fix_setitem_type(stmt, typemap, calltypes)
            for T, f in copy_propagate_extensions.items():
                if isinstance(stmt, T):
                    gen_set, kill_set = f(stmt, typemap)
                    for lhs, rhs in gen_set:
                        if rhs in name_var_table:
                            var_dict[lhs] = name_var_table[rhs]
                    for l, r in var_dict.copy().items():
                        if l in kill_set or r.name in kill_set:
                            var_dict.pop(l)
            if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Var):
                lhs = stmt.target.name
                rhs = stmt.value.name
                # rhs could be replaced with lhs from previous copies
                if lhs != rhs:
                    # copy is valid only if same type (see
                    # TestCFunc.test_locals)
                    if typemap[lhs] == typemap[rhs] and rhs in name_var_table:
                        var_dict[lhs] = name_var_table[rhs]
                    else:
                        var_dict.pop(lhs, None)
                    # a=b kills previous t=a
                    lhs_kill = []
                    for k, v in var_dict.items():
                        if v.name == lhs:
                            lhs_kill.append(k)
                    for k in lhs_kill:
                        var_dict.pop(k, None)
            if (isinstance(stmt, ir.Assign)
                                        and not isinstance(stmt.value, ir.Var)):
                lhs = stmt.target.name
                var_dict.pop(lhs, None)
                # previous t=a is killed if a is killed
                lhs_kill = []
                for k, v in var_dict.items():
                    if v.name == lhs:
                        lhs_kill.append(k)
                for k in lhs_kill:
                    var_dict.pop(k, None)
        save_copies.extend(var_dict.items())

    return save_copies

def fix_setitem_type(stmt, typemap, calltypes):
    """Copy propagation can replace setitem target variable, which can be array
    with 'A' layout. The replaced variable can be 'C' or 'F', so we update
    setitem call type reflect this (from matrix power test)
    """
    if not isinstance(stmt, (ir.SetItem, ir.StaticSetItem)):
        return
    t_typ = typemap[stmt.target.name]
    s_typ = calltypes[stmt].args[0]
    # test_optional t_typ can be Optional with array
    if not isinstance(
            s_typ,
            types.npytypes.Array) or not isinstance(
            t_typ,
            types.npytypes.Array):
        return
    if s_typ.layout == 'A' and t_typ.layout != 'A':
        new_s_typ = s_typ.copy(layout=t_typ.layout)
        calltypes[stmt].args = (
            new_s_typ,
            calltypes[stmt].args[1],
            calltypes[stmt].args[2])
    return


def dprint_func_ir(func_ir, title, blocks=None):
    """Debug print function IR, with an optional blocks argument
    that may differ from the IR's original blocks.
    """
    if config.DEBUG_ARRAY_OPT >= 1:
        ir_blocks = func_ir.blocks
        func_ir.blocks = ir_blocks if blocks == None else blocks
        name = func_ir.func_id.func_qualname
        print(("IR %s: %s" % (title, name)).center(80, "-"))
        func_ir.dump()
        print("-" * 40)
        func_ir.blocks = ir_blocks


def find_topo_order(blocks, cfg = None):
    """find topological order of blocks such that true branches are visited
    first (e.g. for_break test in test_dataflow).
    """
    if cfg is None:
        cfg = compute_cfg_from_blocks(blocks)
    post_order = []
    seen = set()

    def _dfs_rec(node):
        if node not in seen:
            seen.add(node)
            succs = cfg._succs[node]
            last_inst = blocks[node].body[-1]
            if isinstance(last_inst, ir.Branch):
                succs = [last_inst.falsebr, last_inst.truebr]
            for dest in succs:
                if (node, dest) not in cfg._back_edges:
                    _dfs_rec(dest)
            post_order.append(node)

    _dfs_rec(cfg.entry_point())
    post_order.reverse()
    return post_order


# other packages that define new nodes add calls to get call table
# format: {type:function}
call_table_extensions = {}


def get_call_table(blocks, call_table=None, reverse_call_table=None, topological_ordering=True):
    """returns a dictionary of call variables and their references.
    """
    # call_table example: c = np.zeros becomes c:["zeroes", np]
    # reverse_call_table example: c = np.zeros becomes np_var:c
    if call_table is None:
        call_table = {}
    if reverse_call_table is None:
        reverse_call_table = {}

    if topological_ordering:
        order = find_topo_order(blocks)
    else:
        order = list(blocks.keys())

    for label in reversed(order):
        for inst in reversed(blocks[label].body):
            if isinstance(inst, ir.Assign):
                lhs = inst.target.name
                rhs = inst.value
                if isinstance(rhs, ir.Expr) and rhs.op == 'call':
                    call_table[rhs.func.name] = []
                if isinstance(rhs, ir.Expr) and rhs.op == 'getattr':
                    if lhs in call_table:
                        call_table[lhs].append(rhs.attr)
                        reverse_call_table[rhs.value.name] = lhs
                    if lhs in reverse_call_table:
                        call_var = reverse_call_table[lhs]
                        call_table[call_var].append(rhs.attr)
                        reverse_call_table[rhs.value.name] = call_var
                if isinstance(rhs, ir.Global):
                    if lhs in call_table:
                        call_table[lhs].append(rhs.value)
                    if lhs in reverse_call_table:
                        call_var = reverse_call_table[lhs]
                        call_table[call_var].append(rhs.value)
                if isinstance(rhs, ir.FreeVar):
                    if lhs in call_table:
                        call_table[lhs].append(rhs.value)
                    if lhs in reverse_call_table:
                        call_var = reverse_call_table[lhs]
                        call_table[call_var].append(rhs.value)
                if isinstance(rhs, ir.Var):
                    if lhs in call_table:
                        call_table[lhs].append(rhs.name)
                        reverse_call_table[rhs.name] = lhs
                    if lhs in reverse_call_table:
                        call_var = reverse_call_table[lhs]
                        call_table[call_var].append(rhs.name)
            for T, f in call_table_extensions.items():
                if isinstance(inst, T):
                    f(inst, call_table, reverse_call_table)
    return call_table, reverse_call_table


# other packages that define new nodes add calls to get tuple table
# format: {type:function}
tuple_table_extensions = {}


def get_tuple_table(blocks, tuple_table=None):
    """returns a dictionary of tuple variables and their values.
    """
    if tuple_table is None:
        tuple_table = {}

    for block in blocks.values():
        for inst in block.body:
            if isinstance(inst, ir.Assign):
                lhs = inst.target.name
                rhs = inst.value
                if isinstance(rhs, ir.Expr) and rhs.op == 'build_tuple':
                    tuple_table[lhs] = rhs.items
                if isinstance(rhs, ir.Const) and isinstance(rhs.value, tuple):
                    tuple_table[lhs] = rhs.value
            for T, f in tuple_table_extensions.items():
                if isinstance(inst, T):
                    f(inst, tuple_table)
    return tuple_table


def get_stmt_writes(stmt):
    writes = set()
    if isinstance(stmt, (ir.Assign, ir.SetItem, ir.StaticSetItem)):
        writes.add(stmt.target.name)
    return writes


def rename_labels(blocks):
    """rename labels of function body blocks according to topological sort.
    The set of labels of these blocks will remain unchanged.
    """
    topo_order = find_topo_order(blocks)

    # make a block with return last if available (just for readability)
    return_label = -1
    for l, b in blocks.items():
        if isinstance(b.body[-1], ir.Return):
            return_label = l
    # some cases like generators can have no return blocks
    if return_label != -1:
        topo_order.remove(return_label)
        topo_order.append(return_label)

    label_map = {}
    all_labels = sorted(topo_order, reverse=True)
    for label in topo_order:
        label_map[label] = all_labels.pop()
    # update target labels in jumps/branches
    for b in blocks.values():
        term = b.terminator
        if isinstance(term, ir.Jump):
            term.target = label_map[term.target]
        if isinstance(term, ir.Branch):
            term.truebr = label_map[term.truebr]
            term.falsebr = label_map[term.falsebr]
    # update blocks dictionary keys
    new_blocks = {}
    for k, b in blocks.items():
        new_label = label_map[k]
        new_blocks[new_label] = b

    return new_blocks


def simplify_CFG(blocks):
    """transform chains of blocks that have no loop into a single block"""
    # first, inline single-branch-block to its predecessors
    cfg = compute_cfg_from_blocks(blocks)
    def find_single_branch(label):
        block = blocks[label]
        return len(block.body) == 1 and isinstance(block.body[0], ir.Branch)
    single_branch_blocks = list(filter(find_single_branch, blocks.keys()))
    marked_for_del = set()
    for label in single_branch_blocks:
        inst = blocks[label].body[0]
        predecessors = cfg.predecessors(label)
        delete_block = True
        for (p, q) in predecessors:
            block = blocks[p]
            if isinstance(block.body[-1], ir.Jump):
                block.body[-1] = copy.copy(inst)
            else:
                delete_block = False
        if delete_block:
            marked_for_del.add(label)
    # Delete marked labels
    for label in marked_for_del:
        del blocks[label]
    merge_adjacent_blocks(blocks)
    return rename_labels(blocks)


arr_math = ['min', 'max', 'sum', 'prod', 'mean', 'var', 'std',
            'cumsum', 'cumprod', 'argmax', 'argmin', 'argsort',
            'nonzero', 'ravel']


def canonicalize_array_math(func_ir, typemap, calltypes, typingctx):
    # save array arg to call
    # call_varname -> array
    blocks = func_ir.blocks
    saved_arr_arg = {}
    topo_order = find_topo_order(blocks)
    for label in topo_order:
        block = blocks[label]
        new_body = []
        for stmt in block.body:
            if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Expr):
                lhs = stmt.target.name
                rhs = stmt.value
                # replace A.func with np.func, and save A in saved_arr_arg
                if (rhs.op == 'getattr' and rhs.attr in arr_math
                        and isinstance(
                            typemap[rhs.value.name], types.npytypes.Array)):
                    rhs = stmt.value
                    arr = rhs.value
                    saved_arr_arg[lhs] = arr
                    scope = arr.scope
                    loc = arr.loc
                    # g_np_var = Global(numpy)
                    g_np_var = ir.Var(scope, mk_unique_var("$np_g_var"), loc)
                    typemap[g_np_var.name] = types.misc.Module(numpy)
                    g_np = ir.Global('np', numpy, loc)
                    g_np_assign = ir.Assign(g_np, g_np_var, loc)
                    rhs.value = g_np_var
                    new_body.append(g_np_assign)
                    func_ir._definitions[g_np_var.name] = [g_np]
                    # update func var type
                    func = getattr(numpy, rhs.attr)
                    func_typ = get_np_ufunc_typ(func)
                    typemap.pop(lhs)
                    typemap[lhs] = func_typ
                if rhs.op == 'call' and rhs.func.name in saved_arr_arg:
                    # add array as first arg
                    arr = saved_arr_arg[rhs.func.name]
                    # update call type signature to include array arg
                    old_sig = calltypes.pop(rhs)
                    # argsort requires kws for typing so sig.args can't be used
                    # reusing sig.args since some types become Const in sig
                    argtyps = old_sig.args[:len(rhs.args)]
                    kwtyps = {name: typemap[v.name] for name, v in rhs.kws}
                    calltypes[rhs] = typemap[rhs.func.name].get_call_type(
                        typingctx, [typemap[arr.name]] + list(argtyps), kwtyps)
                    rhs.args = [arr] + rhs.args

            new_body.append(stmt)
        block.body = new_body
    return


# format: {type:function}
array_accesses_extensions = {}


def get_array_accesses(blocks, accesses=None):
    """returns a set of arrays accessed and their indices.
    """
    if accesses is None:
        accesses = set()

    for block in blocks.values():
        for inst in block.body:
            if isinstance(inst, ir.SetItem):
                accesses.add((inst.target.name, inst.index.name))
            if isinstance(inst, ir.StaticSetItem):
                accesses.add((inst.target.name, inst.index_var.name))
            if isinstance(inst, ir.Assign):
                lhs = inst.target.name
                rhs = inst.value
                if isinstance(rhs, ir.Expr) and rhs.op == 'getitem':
                    accesses.add((rhs.value.name, rhs.index.name))
                if isinstance(rhs, ir.Expr) and rhs.op == 'static_getitem':
                    index = rhs.index
                    # slice is unhashable, so just keep the variable
                    if index is None or is_slice_index(index):
                        index = rhs.index_var.name
                    accesses.add((rhs.value.name, index))
            for T, f in array_accesses_extensions.items():
                if isinstance(inst, T):
                    f(inst, accesses)
    return accesses

def is_slice_index(index):
    """see if index is a slice index or has slice in it"""
    if isinstance(index, slice):
        return True
    if isinstance(index, tuple):
        for i in index:
            if isinstance(i, slice):
                return True
    return False

def merge_adjacent_blocks(blocks):
    cfg = compute_cfg_from_blocks(blocks)
    # merge adjacent blocks
    removed = set()
    for label in list(blocks.keys()):
        if label in removed:
            continue
        block = blocks[label]
        succs = list(cfg.successors(label))
        while True:
            if len(succs) != 1:
                break
            next_label = succs[0][0]
            if next_label in removed:
                break
            preds = list(cfg.predecessors(next_label))
            succs = list(cfg.successors(next_label))
            if len(preds) != 1 or preds[0][0] != label:
                break
            next_block = blocks[next_label]
            # XXX: commented out since scope objects are not consistent
            # throughout the compiler. for example, pieces of code are compiled
            # and inlined on the fly without proper scope merge.
            # if block.scope != next_block.scope:
            #     break
            # merge
            block.body.pop()  # remove Jump
            block.body += next_block.body
            del blocks[next_label]
            removed.add(next_label)
            label = next_label


def restore_copy_var_names(blocks, save_copies, typemap):
    """
    restores variable names of user variables after applying copy propagation
    """
    if not save_copies:
        return {}

    rename_dict = {}
    var_rename_map = {}
    for (a, b) in save_copies:
        # a is string name, b is variable
        # if a is user variable and b is generated temporary and b is not
        # already renamed
        if (not a.startswith('$') and b.name.startswith('$')
                                                and b.name not in rename_dict):
            new_name = mk_unique_var('${}'.format(a));
            rename_dict[b.name] = new_name
            var_rename_map[new_name] = a
            typ = typemap.pop(b.name)
            typemap[new_name] = typ

    replace_var_names(blocks, rename_dict)
    return var_rename_map


def simplify(func_ir, typemap, calltypes, metadata):
    # get copies in to blocks and out from blocks
    in_cps, _ = copy_propagate(func_ir.blocks, typemap)
    # table mapping variable names to ir.Var objects to help replacement
    name_var_table = get_name_var_table(func_ir.blocks)
    save_copies = apply_copy_propagate(
        func_ir.blocks,
        in_cps,
        name_var_table,
        typemap,
        calltypes)
    var_rename_map = restore_copy_var_names(func_ir.blocks, save_copies, typemap)
    if "var_rename_map" not in metadata:
            metadata["var_rename_map"] = {}
    metadata["var_rename_map"].update(var_rename_map)
    # remove dead code to enable fusion
    if config.DEBUG_ARRAY_OPT >= 1:
        dprint_func_ir(func_ir, "after copy prop")
    remove_dead(func_ir.blocks, func_ir.arg_names, func_ir, typemap)
    func_ir.blocks = simplify_CFG(func_ir.blocks)
    if config.DEBUG_ARRAY_OPT >= 1:
        dprint_func_ir(func_ir, "after simplify")


class GuardException(Exception):
    pass


def require(cond):
    """
    Raise GuardException if the given condition is False.
    """
    if not cond:
       raise GuardException

def guard(func, *args, **kwargs):
    """
    Run a function with given set of arguments, and guard against
    any GuardException raised by the function by returning None,
    or the expected return results if no such exception was raised.
    """
    try:
        return func(*args, **kwargs)
    except GuardException:
        return None

def get_definition(func_ir, name, **kwargs):
    """
    Same as func_ir.get_definition(name), but raise GuardException if
    exception KeyError is caught.
    """
    try:
        return func_ir.get_definition(name, **kwargs)
    except KeyError:
        raise GuardException

def build_definitions(blocks, definitions=None):
    """Build the definitions table of the given blocks by scanning
    through all blocks and instructions, useful when the definitions
    table is out-of-sync.
    Will return a new definition table if one is not passed.
    """
    if definitions is None:
        definitions = collections.defaultdict(list)

    for block in blocks.values():
        for inst in block.body:
            if isinstance(inst, ir.Assign):
                name = inst.target.name
                definition = definitions.get(name, [])
                if definition == []:
                    definitions[name] = definition
                definition.append(inst.value)
            if type(inst) in build_defs_extensions:
                f = build_defs_extensions[type(inst)]
                f(inst, definitions)

    return definitions

build_defs_extensions = {}

def find_callname(func_ir, expr, typemap=None, definition_finder=get_definition):
    """Try to find a call expression's function and module names and return
    them as strings for unbounded calls. If the call is a bounded call, return
    the self object instead of module name. Raise GuardException if failed.

    Providing typemap can make the call matching more accurate in corner cases
    such as bounded call on an object which is inside another object.
    """
    require(isinstance(expr, ir.Expr) and expr.op == 'call')
    callee = expr.func
    callee_def = definition_finder(func_ir, callee)
    attrs = []
    obj = None
    while True:
        if isinstance(callee_def, (ir.Global, ir.FreeVar)):
            # require(callee_def.value == numpy)
            # these checks support modules like numpy, numpy.random as well as
            # calls like len() and intrinsics like assertEquiv
            keys = ['name', '_name', '__name__']
            value = None
            for key in keys:
                if hasattr(callee_def.value, key):
                    value = getattr(callee_def.value, key)
                    break
            if not value or not isinstance(value, str):
                raise GuardException
            attrs.append(value)
            def_val = callee_def.value
            # get the underlying definition of Intrinsic object to be able to
            # find the module effectively.
            # Otherwise, it will return numba.extending
            if isinstance(def_val, _Intrinsic):
                def_val = def_val._defn
            if hasattr(def_val, '__module__'):
                mod_name = def_val.__module__
                # The reason for first checking if the function is in NumPy's
                # top level name space by module is that some functions are
                # deprecated in NumPy but the functions' names are aliased with
                # other common names. This prevents deprecation warnings on
                # e.g. getattr(numpy, 'bool') were a bool the target.
                # For context see #6175, impacts NumPy>=1.20.
                mod_not_none = mod_name is not None
                numpy_toplevel = (mod_not_none and
                                  (mod_name == 'numpy'
                                   or mod_name.startswith('numpy.')))
                # it might be a numpy function imported directly
                if (numpy_toplevel and hasattr(numpy, value)
                        and def_val == getattr(numpy, value)):
                    attrs += ['numpy']
                # it might be a np.random function imported directly
                elif (hasattr(numpy.random, value)
                        and def_val == getattr(numpy.random, value)):
                    attrs += ['random', 'numpy']
                elif mod_not_none:
                    attrs.append(mod_name)
            else:
                class_name = def_val.__class__.__name__
                if class_name == 'builtin_function_or_method':
                    class_name = 'builtin'
                if class_name != 'module':
                    attrs.append(class_name)
            break
        elif isinstance(callee_def, ir.Expr) and callee_def.op == 'getattr':
            obj = callee_def.value
            attrs.append(callee_def.attr)
            if typemap and obj.name in typemap:
                typ = typemap[obj.name]
                if not isinstance(typ, types.Module):
                    return attrs[0], obj
            callee_def = definition_finder(func_ir, obj)
        else:
            # obj.func calls where obj is not np array
            if obj is not None:
                return '.'.join(reversed(attrs)), obj
            raise GuardException
    return attrs[0], '.'.join(reversed(attrs[1:]))

def find_build_sequence(func_ir, var):
    """Check if a variable is constructed via build_tuple or
    build_list or build_set, and return the sequence and the
    operator, or raise GuardException otherwise.
    Note: only build_tuple is immutable, so use with care.
    """
    require(isinstance(var, ir.Var))
    var_def = get_definition(func_ir, var)
    require(isinstance(var_def, ir.Expr))
    build_ops = ['build_tuple', 'build_list', 'build_set']
    require(var_def.op in build_ops)
    return var_def.items, var_def.op

def find_const(func_ir, var):
    """Check if a variable is defined as constant, and return
    the constant value, or raise GuardException otherwise.
    """
    require(isinstance(var, ir.Var))
    var_def = get_definition(func_ir, var)
    require(isinstance(var_def, (ir.Const, ir.Global, ir.FreeVar)))
    return var_def.value

def compile_to_numba_ir(mk_func, glbls, typingctx=None, targetctx=None,
                        arg_typs=None, typemap=None, calltypes=None):
    """
    Compile a function or a make_function node to Numba IR.

    Rename variables and
    labels to avoid conflict if inlined somewhere else. Perform type inference
    if typingctx and other typing inputs are available and update typemap and
    calltypes.
    """
    from numba.core import typed_passes
    # mk_func can be actual function or make_function node, or a njit function
    if hasattr(mk_func, 'code'):
        code = mk_func.code
    elif hasattr(mk_func, '__code__'):
        code = mk_func.__code__
    else:
        raise NotImplementedError("function type not recognized {}".format(mk_func))
    f_ir = get_ir_of_code(glbls, code)
    remove_dels(f_ir.blocks)

    # relabel by adding an offset
    f_ir.blocks = add_offset_to_labels(f_ir.blocks, _the_max_label.next())
    max_label = max(f_ir.blocks.keys())
    _the_max_label.update(max_label)

    # rename all variables to avoid conflict
    var_table = get_name_var_table(f_ir.blocks)
    new_var_dict = {}
    for name, var in var_table.items():
        new_var_dict[name] = mk_unique_var(name)
    replace_var_names(f_ir.blocks, new_var_dict)

    # perform type inference if typingctx is available and update type
    # data structures typemap and calltypes
    if typingctx:
        f_typemap, f_return_type, f_calltypes, _ = typed_passes.type_inference_stage(
                typingctx, targetctx, f_ir, arg_typs, None)
        # remove argument entries like arg.a from typemap
        arg_names = [vname for vname in f_typemap if vname.startswith("arg.")]
        for a in arg_names:
            f_typemap.pop(a)
        typemap.update(f_typemap)
        calltypes.update(f_calltypes)
    return f_ir

def _create_function_from_code_obj(fcode, func_env, func_arg, func_clo, glbls):
    """
    Creates a function from a code object. Args:
    * fcode - the code object
    * func_env - string for the freevar placeholders
    * func_arg - string for the function args (e.g. "a, b, c, d=None")
    * func_clo - string for the closure args
    * glbls - the function globals
    """
    sanitized_co_name = fcode.co_name.replace('<', '_').replace('>', '_')
    func_text = (f"def closure():\n{func_env}\n"
                 f"\tdef {sanitized_co_name}({func_arg}):\n"
                 f"\t\treturn ({func_clo})\n"
                 f"\treturn {sanitized_co_name}")
    loc = {}
    exec(func_text, glbls, loc)

    f = loc['closure']()
    # replace the code body
    f.__code__ = fcode
    f.__name__ = fcode.co_name
    return f

def get_ir_of_code(glbls, fcode):
    """
    Compile a code object to get its IR, ir.Del nodes are emitted
    """
    nfree = len(fcode.co_freevars)
    func_env = "\n".join(["\tc_%d = None" % i for i in range(nfree)])
    func_clo = ",".join(["c_%d" % i for i in range(nfree)])
    func_arg = ",".join(["x_%d" % i for i in range(fcode.co_argcount)])

    f = _create_function_from_code_obj(fcode, func_env, func_arg, func_clo,
                                       glbls)

    from numba.core import compiler
    ir = compiler.run_frontend(f)
    # we need to run the before inference rewrite pass to normalize the IR
    # XXX: check rewrite pass flag?
    # for example, Raise nodes need to become StaticRaise before type inference
    class DummyPipeline(object):
        def __init__(self, f_ir):
            self.state = compiler.StateDict()
            self.state.typingctx = None
            self.state.targetctx = None
            self.state.args = None
            self.state.func_ir = f_ir
            self.state.typemap = None
            self.state.return_type = None
            self.state.calltypes = None
    state = DummyPipeline(ir).state
    rewrites.rewrite_registry.apply('before-inference', state)
    # call inline pass to handle cases like stencils and comprehensions
    swapped = {} # TODO: get this from diagnostics store
    import numba.core.inline_closurecall
    inline_pass = numba.core.inline_closurecall.InlineClosureCallPass(
        ir, numba.core.cpu.ParallelOptions(False), swapped)
    inline_pass.run()

    # TODO: DO NOT ADD MORE THINGS HERE!
    # If adding more things here is being contemplated, it really is time to
    # retire this function and work on getting the InlineWorker class from
    # numba.core.inline_closurecall into sufficient shape as a replacement.
    # The issue with `get_ir_of_code` is that it doesn't run a full compilation
    # pipeline and as a result various additional things keep needing to be
    # added to create valid IR.

    # rebuild IR in SSA form
    from numba.core.untyped_passes import ReconstructSSA
    from numba.core.typed_passes import PreLowerStripPhis
    reconstruct_ssa = ReconstructSSA()
    phistrip = PreLowerStripPhis()
    reconstruct_ssa.run_pass(state)
    phistrip.run_pass(state)

    post_proc = postproc.PostProcessor(ir)
    post_proc.run(True)
    return ir

def replace_arg_nodes(block, args):
    """
    Replace ir.Arg(...) with variables
    """
    for stmt in block.body:
        if isinstance(stmt, ir.Assign) and isinstance(stmt.value, ir.Arg):
            idx = stmt.value.index
            assert(idx < len(args))
            stmt.value = args[idx]
    return


def replace_returns(blocks, target, return_label):
    """
    Return return statement by assigning directly to target, and a jump.
    """
    for block in blocks.values():
        # some blocks may be empty during transformations
        if not block.body:
            continue
        stmt = block.terminator
        if isinstance(stmt, ir.Return):
            block.body.pop()  # remove return
            cast_stmt = block.body.pop()
            assert (isinstance(cast_stmt, ir.Assign)
                and isinstance(cast_stmt.value, ir.Expr)
                and cast_stmt.value.op == 'cast'), "invalid return cast"
            block.body.append(ir.Assign(cast_stmt.value.value, target, stmt.loc))
            block.body.append(ir.Jump(return_label, stmt.loc))


def gen_np_call(func_as_str, func, lhs, args, typingctx, typemap, calltypes):
    scope = args[0].scope
    loc = args[0].loc

    # g_np_var = Global(numpy)
    g_np_var = ir.Var(scope, mk_unique_var("$np_g_var"), loc)
    typemap[g_np_var.name] = types.misc.Module(numpy)
    g_np = ir.Global('np', numpy, loc)
    g_np_assign = ir.Assign(g_np, g_np_var, loc)
    # attr call: <something>_attr = getattr(g_np_var, func_as_str)
    np_attr_call = ir.Expr.getattr(g_np_var, func_as_str, loc)
    attr_var = ir.Var(scope, mk_unique_var("$np_attr_attr"), loc)
    func_var_typ = get_np_ufunc_typ(func)
    typemap[attr_var.name] = func_var_typ
    attr_assign = ir.Assign(np_attr_call, attr_var, loc)
    # np call: lhs = np_attr(*args)
    np_call = ir.Expr.call(attr_var, args, (), loc)
    arg_types = [typemap[x.name] for x in args]
    func_typ = func_var_typ.get_call_type(typingctx, arg_types, {})
    calltypes[np_call] = func_typ
    np_assign = ir.Assign(np_call, lhs, loc)
    return [g_np_assign, attr_assign, np_assign]

def dump_block(label, block):
    print(label, ":")
    for stmt in block.body:
        print("    ", stmt)

def dump_blocks(blocks):
    for label, block in blocks.items():
        dump_block(label, block)

def is_operator_or_getitem(expr):
    """true if expr is unary or binary operator or getitem"""
    return (isinstance(expr, ir.Expr)
            and getattr(expr, 'op', False)
            and expr.op in ['unary', 'binop', 'inplace_binop', 'getitem', 'static_getitem'])

def is_get_setitem(stmt):
    """stmt is getitem assignment or setitem (and static cases)"""
    return is_getitem(stmt) or is_setitem(stmt)


def is_getitem(stmt):
    """true if stmt is a getitem or static_getitem assignment"""
    return (isinstance(stmt, ir.Assign)
            and isinstance(stmt.value, ir.Expr)
            and stmt.value.op in ['getitem', 'static_getitem'])

def is_setitem(stmt):
    """true if stmt is a SetItem or StaticSetItem node"""
    return isinstance(stmt, (ir.SetItem, ir.StaticSetItem))

def index_var_of_get_setitem(stmt):
    """get index variable for getitem/setitem nodes (and static cases)"""
    if is_getitem(stmt):
        if stmt.value.op == 'getitem':
            return stmt.value.index
        else:
            return stmt.value.index_var

    if is_setitem(stmt):
        if isinstance(stmt, ir.SetItem):
            return stmt.index
        else:
            return stmt.index_var

    return None

def set_index_var_of_get_setitem(stmt, new_index):
    if is_getitem(stmt):
        if stmt.value.op == 'getitem':
            stmt.value.index = new_index
        else:
            stmt.value.index_var = new_index
    elif is_setitem(stmt):
        if isinstance(stmt, ir.SetItem):
            stmt.index = new_index
        else:
            stmt.index_var = new_index
    else:
        raise ValueError("getitem or setitem node expected but received {}".format(
                     stmt))


def is_namedtuple_class(c):
    """check if c is a namedtuple class"""
    if not isinstance(c, type):
        return False
    # should have only tuple as superclass
    bases = c.__bases__
    if len(bases) != 1 or bases[0] != tuple:
        return False
    # should have _make method
    if not hasattr(c, '_make'):
        return False
    # should have _fields that is all string
    fields = getattr(c, '_fields', None)
    if not isinstance(fields, tuple):
        return False
    return all(isinstance(f, str) for f in fields)


def fill_block_with_call(newblock, callee, label_next, inputs, outputs):
    """Fill *newblock* to call *callee* with arguments listed in *inputs*.
    The returned values are unwrapped into variables in *outputs*.
    The block would then jump to *label_next*.
    """
    scope = newblock.scope
    loc = newblock.loc

    fn = ir.Const(value=callee, loc=loc)
    fnvar = scope.make_temp(loc=loc)
    newblock.append(ir.Assign(target=fnvar, value=fn, loc=loc))
    # call
    args = [scope.get_exact(name) for name in inputs]
    callexpr = ir.Expr.call(func=fnvar, args=args, kws=(), loc=loc)
    callres = scope.make_temp(loc=loc)
    newblock.append(ir.Assign(target=callres, value=callexpr, loc=loc))
    # unpack return value
    for i, out in enumerate(outputs):
        target = scope.get_exact(out)
        getitem = ir.Expr.static_getitem(value=callres, index=i,
                                         index_var=None, loc=loc)
        newblock.append(ir.Assign(target=target, value=getitem, loc=loc))
    # jump to next block
    newblock.append(ir.Jump(target=label_next, loc=loc))
    return newblock


def fill_callee_prologue(block, inputs, label_next):
    """
    Fill a new block *block* that unwraps arguments using names in *inputs* and
    then jumps to *label_next*.

    Expected to use with *fill_block_with_call()*
    """
    scope = block.scope
    loc = block.loc
    # load args
    args = [ir.Arg(name=k, index=i, loc=loc)
            for i, k in enumerate(inputs)]
    for aname, aval in zip(inputs, args):
        tmp = ir.Var(scope=scope, name=aname, loc=loc)
        block.append(ir.Assign(target=tmp, value=aval, loc=loc))
    # jump to loop entry
    block.append(ir.Jump(target=label_next, loc=loc))
    return block


def fill_callee_epilogue(block, outputs):
    """
    Fill a new block *block* to prepare the return values.
    This block is the last block of the function.

    Expected to use with *fill_block_with_call()*
    """
    scope = block.scope
    loc = block.loc
    # prepare tuples to return
    vals = [scope.get_exact(name=name) for name in outputs]
    tupexpr = ir.Expr.build_tuple(items=vals, loc=loc)
    tup = scope.make_temp(loc=loc)
    block.append(ir.Assign(target=tup, value=tupexpr, loc=loc))
    # return
    block.append(ir.Return(value=tup, loc=loc))
    return block


def find_outer_value(func_ir, var):
    """Check if a variable is a global value, and return the value,
    or raise GuardException otherwise.
    """
    dfn = get_definition(func_ir, var)
    if isinstance(dfn, (ir.Global, ir.FreeVar)):
        return dfn.value

    if isinstance(dfn, ir.Expr) and dfn.op == 'getattr':
        prev_val = find_outer_value(func_ir, dfn.value)
        try:
            val = getattr(prev_val, dfn.attr)
            return val
        except AttributeError:
            raise GuardException

    raise GuardException


def raise_on_unsupported_feature(func_ir, typemap):
    """
    Helper function to walk IR and raise if it finds op codes
    that are unsupported. Could be extended to cover IR sequences
    as well as op codes. Intended use is to call it as a pipeline
    stage just prior to lowering to prevent LoweringErrors for known
    unsupported features.
    """
    gdb_calls = [] # accumulate calls to gdb/gdb_init

    # issue 2195: check for excessively large tuples
    for arg_name in func_ir.arg_names:
        if arg_name in typemap and \
           isinstance(typemap[arg_name], types.containers.UniTuple) and \
           typemap[arg_name].count > 1000:
            # Raise an exception when len(tuple) > 1000. The choice of this number (1000)
            # was entirely arbitrary
            msg = ("Tuple '{}' length must be smaller than 1000.\n"
                   "Large tuples lead to the generation of a prohibitively large "
                   "LLVM IR which causes excessive memory pressure "
                   "and large compile times.\n"
                   "As an alternative, the use of a 'list' is recommended in "
                   "place of a 'tuple' as lists do not suffer from this problem.".format(arg_name))
            raise UnsupportedError(msg, func_ir.loc)

    for blk in func_ir.blocks.values():
        for stmt in blk.find_insts(ir.Assign):
            # This raises on finding `make_function`
            if isinstance(stmt.value, ir.Expr):
                if stmt.value.op == 'make_function':
                    val = stmt.value

                    # See if the construct name can be refined
                    code = getattr(val, 'code', None)
                    if code is not None:
                        # check if this is a closure, the co_name will
                        # be the captured function name which is not
                        # useful so be explicit
                        if getattr(val, 'closure', None) is not None:
                            use = '<creating a function from a closure>'
                            expr = ''
                        else:
                            use = code.co_name
                            expr = '(%s) ' % use
                    else:
                        use = '<could not ascertain use case>'
                        expr = ''

                    msg = ("Numba encountered the use of a language "
                            "feature it does not support in this context: "
                            "%s (op code: make_function not supported). If "
                            "the feature is explicitly supported it is "
                            "likely that the result of the expression %s"
                            "is being used in an unsupported manner.") % \
                            (use, expr)
                    raise UnsupportedError(msg, stmt.value.loc)

            # this checks for gdb initialization calls, only one is permitted
            if isinstance(stmt.value, (ir.Global, ir.FreeVar)):
                val = stmt.value
                val = getattr(val, 'value', None)
                if val is None:
                    continue

                # check global function
                found = False
                if isinstance(val, pytypes.FunctionType):
                    found = val in {numba.gdb, numba.gdb_init}
                if not found: # freevar bind to intrinsic
                    found = getattr(val, '_name', "") == "gdb_internal"
                if found:
                    gdb_calls.append(stmt.loc) # report last seen location

            # this checks that np.<type> was called if view is called
            if isinstance(stmt.value, ir.Expr):
                if stmt.value.op == 'getattr' and stmt.value.attr == 'view':
                    var = stmt.value.value.name
                    if isinstance(typemap[var], types.Array):
                        continue
                    df = func_ir.get_definition(var)
                    cn = guard(find_callname, func_ir, df)
                    if cn and cn[1] == 'numpy':
                        ty = getattr(numpy, cn[0])
                        if (numpy.issubdtype(ty, numpy.integer) or
                                numpy.issubdtype(ty, numpy.floating)):
                            continue

                    vardescr = '' if var.startswith('$') else "'{}' ".format(var)
                    raise TypingError(
                        "'view' can only be called on NumPy dtypes, "
                        "try wrapping the variable {}with 'np.<dtype>()'".
                        format(vardescr), loc=stmt.loc)

            # checks for globals that are also reflected
            if isinstance(stmt.value, ir.Global):
                ty = typemap[stmt.target.name]
                msg = ("The use of a %s type, assigned to variable '%s' in "
                       "globals, is not supported as globals are considered "
                       "compile-time constants and there is no known way to "
                       "compile a %s type as a constant.")
                if (getattr(ty, 'reflected', False) or
                    isinstance(ty, (types.DictType, types.ListType))):
                    raise TypingError(msg % (ty, stmt.value.name, ty), loc=stmt.loc)

            # checks for generator expressions (yield in use when func_ir has
            # not been identified as a generator).
            if isinstance(stmt.value, ir.Yield) and not func_ir.is_generator:
                msg = "The use of generator expressions is unsupported."
                raise UnsupportedError(msg, loc=stmt.loc)

    # There is more than one call to function gdb/gdb_init
    if len(gdb_calls) > 1:
        msg = ("Calling either numba.gdb() or numba.gdb_init() more than once "
               "in a function is unsupported (strange things happen!), use "
               "numba.gdb_breakpoint() to create additional breakpoints "
               "instead.\n\nRelevant documentation is available here:\n"
               "https://numba.readthedocs.io/en/stable/user/troubleshoot.html"
               "#using-numba-s-direct-gdb-bindings-in-nopython-mode\n\n"
               "Conflicting calls found at:\n %s")
        buf = '\n'.join([x.strformat() for x in gdb_calls])
        raise UnsupportedError(msg % buf)


def warn_deprecated(func_ir, typemap):
    # first pass, just walk the type map
    for name, ty in typemap.items():
        # the Type Metaclass has a reflected member
        if ty.reflected:
            # if its an arg, report function call
            if name.startswith('arg.'):
                loc = func_ir.loc
                arg = name.split('.')[1]
                fname = func_ir.func_id.func_qualname
                tyname = 'list' if isinstance(ty, types.List) else 'set'
                url = ("https://numba.readthedocs.io/en/stable/reference/"
                       "deprecation.html#deprecation-of-reflection-for-list-and"
                       "-set-types")
                msg = ("\nEncountered the use of a type that is scheduled for "
                        "deprecation: type 'reflected %s' found for argument "
                        "'%s' of function '%s'.\n\nFor more information visit "
                        "%s" % (tyname, arg, fname, url))
                warnings.warn(NumbaPendingDeprecationWarning(msg, loc=loc))


def resolve_func_from_module(func_ir, node):
    """
    This returns the python function that is being getattr'd from a module in
    some IR, it resolves import chains/submodules recursively. Should it not be
    possible to find the python function being called None will be returned.

    func_ir - the FunctionIR object
    node - the IR node from which to start resolving (should be a `getattr`).
    """
    getattr_chain = []
    def resolve_mod(mod):
        if getattr(mod, 'op', False) == 'getattr':
            getattr_chain.insert(0, mod.attr)
            try:
                mod = func_ir.get_definition(mod.value)
            except KeyError: # multiple definitions
                return None
            return resolve_mod(mod)
        elif isinstance(mod, (ir.Global, ir.FreeVar)):
            if isinstance(mod.value, pytypes.ModuleType):
                return mod
        return None

    mod = resolve_mod(node)
    if mod is not None:
        defn = mod.value
        for x in getattr_chain:
            defn = getattr(defn, x, False)
            if not defn:
                break
        else:
            return defn
    else:
        return None


def enforce_no_dels(func_ir):
    """
    Enforce there being no ir.Del nodes in the IR.
    """
    for blk in func_ir.blocks.values():
        dels = [x for x in blk.find_insts(ir.Del)]
        if dels:
            msg = "Illegal IR, del found at: %s" % dels[0]
            raise CompilerError(msg, loc=dels[0].loc)

def enforce_no_phis(func_ir):
    """
    Enforce there being no ir.Expr.phi nodes in the IR.
    """
    for blk in func_ir.blocks.values():
        phis = [x for x in blk.find_exprs(op='phi')]
        if phis:
            msg = "Illegal IR, phi found at: %s" % phis[0]
            raise CompilerError(msg, loc=phis[0].loc)


def legalize_single_scope(blocks):
    """Check the given mapping of ir.Block for containing a single scope.
    """
    return len({blk.scope for blk in blocks.values()}) == 1


def check_and_legalize_ir(func_ir, flags: "numba.core.compiler.Flags"):
    """
    This checks that the IR presented is legal
    """
    enforce_no_phis(func_ir)
    enforce_no_dels(func_ir)
    # postprocess and emit ir.Dels
    post_proc = postproc.PostProcessor(func_ir)
    post_proc.run(True, extend_lifetimes=flags.dbg_extend_lifetimes)

def convert_code_obj_to_function(code_obj, caller_ir):
    """
    Converts a code object from a `make_function.code` attr in the IR into a
    python function, caller_ir is the FunctionIR of the caller and is used for
    the resolution of freevars.
    """
    fcode = code_obj.code
    nfree = len(fcode.co_freevars)

    # try and resolve freevars if they are consts in the caller's IR
    # these can be baked into the new function
    freevars = []
    for x in fcode.co_freevars:
        # not using guard here to differentiate between multiple definition and
        # non-const variable
        try:
            freevar_def = caller_ir.get_definition(x)
        except KeyError:
            msg = ("Cannot capture a constant value for variable '%s' as there "
                   "are multiple definitions present." % x)
            raise TypingError(msg, loc=code_obj.loc)
        if isinstance(freevar_def, ir.Const):
            freevars.append(freevar_def.value)
        else:
            msg = ("Cannot capture the non-constant value associated with "
                   "variable '%s' in a function that may escape." % x)
            raise TypingError(msg, loc=code_obj.loc)

    func_env = "\n".join(["\tc_%d = %s" % (i, x) for i, x in enumerate(freevars)])
    func_clo = ",".join(["c_%d" % i for i in range(nfree)])
    co_varnames = list(fcode.co_varnames)

    # This is horrible. The code object knows about the number of args present
    # it also knows the name of the args but these are bundled in with other
    # vars in `co_varnames`. The make_function IR node knows what the defaults
    # are, they are defined in the IR as consts. The following finds the total
    # number of args (args + kwargs with defaults), finds the default values
    # and infers the number of "kwargs with defaults" from this and then infers
    # the number of actual arguments from that.
    n_kwargs = 0
    n_allargs = fcode.co_argcount
    kwarg_defaults = caller_ir.get_definition(code_obj.defaults)
    if kwarg_defaults is not None:
        if isinstance(kwarg_defaults, tuple):
            d = [caller_ir.get_definition(x).value for x in kwarg_defaults]
            kwarg_defaults_tup = tuple(d)
        else:
            d = [caller_ir.get_definition(x).value
                 for x in kwarg_defaults.items]
            kwarg_defaults_tup = tuple(d)
        n_kwargs = len(kwarg_defaults_tup)
    nargs = n_allargs - n_kwargs

    func_arg = ",".join(["%s" % (co_varnames[i]) for i in range(nargs)])
    if n_kwargs:
        kw_const = ["%s = %s" % (co_varnames[i + nargs], kwarg_defaults_tup[i])
                    for i in range(n_kwargs)]
        func_arg += ", "
        func_arg += ", ".join(kw_const)

    # globals are the same as those in the caller
    glbls = caller_ir.func_id.func.__globals__

    # create the function and return it
    return _create_function_from_code_obj(fcode, func_env, func_arg, func_clo,
                                          glbls)


def fixup_var_define_in_scope(blocks):
    """Fixes the mapping of ir.Block to ensure all referenced ir.Var are
    defined in every scope used by the function. Such that looking up a variable
    from any scope in this function will not fail.

    Note: This is a workaround. Ideally, all the blocks should refer to the
    same ir.Scope, but that property is not maintained by all the passes.
    """
    # Scan for all used variables
    used_var = {}
    for blk in blocks.values():
        scope = blk.scope
        for inst in blk.body:
            for var in inst.list_vars():
                used_var[var] = inst
    # Note: not all blocks share a single scope even though they should.
    # Ensure the scope of each block defines all used variables.
    for blk in blocks.values():
        scope = blk.scope
        for var, inst in used_var.items():
            # add this variable if it's not in scope
            if var.name not in scope.localvars:
                # Note: using a internal method to reuse the same
                scope.localvars.define(var.name, var)


def transfer_scope(block, scope):
    """Transfer the ir.Block to use the given ir.Scope.
    """
    old_scope = block.scope
    if old_scope is scope:
        # bypass if the block is already using the given scope
        return block
    # Ensure variables are defined in the new scope
    for var in old_scope.localvars._con.values():
        if var.name not in scope.localvars:
            scope.localvars.define(var.name, var)
    # replace scope
    block.scope = scope
    return block


def is_setup_with(stmt):
    return isinstance(stmt, ir.EnterWith)


def is_terminator(stmt):
    return isinstance(stmt, ir.Terminator)


def is_raise(stmt):
    return isinstance(stmt, ir.Raise)


def is_return(stmt):
    return isinstance(stmt, ir.Return)


def is_pop_block(stmt):
    return isinstance(stmt, ir.PopBlock)
