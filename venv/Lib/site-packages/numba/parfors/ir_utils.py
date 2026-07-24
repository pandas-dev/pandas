import numpy

from numba.core import types, ir
from numba.np import types as npy_types

from numba.core.ir_utils import (convert_size_to_var, mk_unique_var,
                                 get_np_ufunc_typ)


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
        typemap[attr_var.name] = get_np_ufunc_typ(numpy.empty, typingctx)
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
    if (
        isinstance(dtype, (npy_types.NPDatetime, npy_types.NPTimedelta))
        and dtype.unit != ''
    ):
        typename_const = ir.Const(dtype_str, loc)
        typ_var_assign = ir.Assign(typename_const, typ_var, loc)
    else:
        if dtype_str == 'bool':
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
        asfortranarray_attr_call = ir.Expr.getattr(g_np_var, "asfortranarray",
                                                   loc)
        afa_attr_var = ir.Var(scope, mk_unique_var("$asfortran_array_attr"),
                              loc)
        if typemap:
            typemap[afa_attr_var.name] = \
                get_np_ufunc_typ(numpy.asfortranarray, typingctx)
        afa_attr_assign = ir.Assign(asfortranarray_attr_call, afa_attr_var, loc)
        # call asfortranarray
        asfortranarray_call = ir.Expr.call(afa_attr_var, [empty_c_var], (), loc)
        if calltypes:
            calltypes[asfortranarray_call] = \
                typemap[afa_attr_var.name].get_call_type(
                    typingctx, [empty_c_typ], {})

        asfortranarray_assign = ir.Assign(asfortranarray_call, lhs, loc)

        out.extend([g_np_assign, attr_assign, typ_var_assign, empty_c_assign,
                    afa_attr_assign, asfortranarray_assign])
    else:
        alloc_assign = ir.Assign(alloc_call, lhs, loc)
        out.extend([g_np_assign, attr_assign, typ_var_assign, alloc_assign])

    return out
