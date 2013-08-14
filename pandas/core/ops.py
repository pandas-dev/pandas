"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
import operator
from functools import partial
import numpy as np
from pandas import compat, lib, tslib
import pandas.index as _index
from pandas.util.decorators import Appender
import pandas.core.common as com
import pandas.core.array as pa
import pandas.core.expressions as expressions
from pandas.core.common import(bind_method, is_list_like, notnull, isnull,
                               _values_from_object, _np_version_under1p6,
                               _np_version_under1p7,
                               _maybe_match_name)

# -----------------------------------------------------------------------------
# Functions that add arithmetic methods to objects, given arithmetic factory
# methods

def _create_methods(arith_method, radd_func, comp_method, bool_method,
                    use_numexpr, special=False, default_axis='columns'):
    # NOTE: Only frame cares about default_axis, specifically: special methods
    # have default axis None, whereas flex methods have default axis 'columns'
    # if we're not using numexpr, then don't pass a str_rep
    if use_numexpr:
        op = lambda x: x
    else:
        op = lambda x: None
    if special:
        def names(x):
            if x[-1] == "_":
                return "__%s_" % x
            else:
                return "__%s__" % x
    else:
        names = lambda x: x
    radd_func = radd_func or operator.add
    # Inframe, all special methods have default_axis=None, flex methods have default_axis set to the default (columns)
    new_methods = dict(
        add=arith_method(operator.add, names('add'), op('+'), default_axis=default_axis),
        radd=arith_method(radd_func, names('radd'), op('+'), default_axis=default_axis),
        sub=arith_method(operator.sub, names('sub'), op('-'), default_axis=default_axis),
        mul=arith_method(operator.mul, names('mul'), op('*'), default_axis=default_axis),
        truediv=arith_method(operator.truediv, names('truediv'), op('/'),
                             truediv=True, fill_zeros=np.inf, default_axis=default_axis),
        floordiv=arith_method(operator.floordiv, names('floordiv'), op('//'),
                              default_axis=default_axis, fill_zeros=np.inf),
        # Causes a floating point exception in the tests when numexpr
        # enabled, so for now no speedup
        mod=arith_method(operator.mod, names('mod'), default_axis=default_axis,
                         fill_zeros=np.nan),
        pow=arith_method(operator.pow, names('pow'), op('**'), default_axis=default_axis),
        # not entirely sure why this is necessary, but previously was included
        # so it's here to maintain compatibility
        rmul=arith_method(operator.mul, names('rmul'), default_axis=default_axis),
        rsub=arith_method(lambda x, y: y - x, names('rsub'), default_axis=default_axis),
        rtruediv=arith_method(lambda x, y: operator.truediv(y, x), names('rtruediv'),
                              truediv=True, fill_zeros=np.inf, default_axis=default_axis),
        rfloordiv=arith_method(lambda x, y: operator.floordiv(y, x), names('rfloordiv'),
                               default_axis=default_axis, fill_zeros=np.inf),
        rpow=arith_method(lambda x, y: y ** x, names('rpow'), default_axis=default_axis),
        rmod=arith_method(lambda x, y: y % x, names('rmod'), default_axis=default_axis),
    )
    if not compat.PY3:
        new_methods["div"] = arith_method(operator.div, names('div'), op('/'),
                                          truediv=False, fill_zeros=np.inf, default_axis=default_axis)
        new_methods["rdiv"] = arith_method(lambda x, y: operator.div(y, x), names('rdiv'),
                                           truediv=False, fill_zeros=np.inf, default_axis=default_axis)
    else:
        new_methods["div"] = arith_method(operator.truediv, names('div'), op('/'),
                                          truediv=True, fill_zeros=np.inf, default_axis=default_axis)
        new_methods["rdiv"] = arith_method(lambda x, y: operator.truediv(y, x), names('rdiv'),
                                           truediv=False, fill_zeros=np.inf, default_axis=default_axis)
        # Comp methods never had a default axis set
    if comp_method:
        new_methods.update(dict(
            eq=comp_method(operator.eq, names('eq'), op('==')),
            ne=comp_method(operator.ne, names('ne'), op('!='), masker=True),
            lt=comp_method(operator.lt, names('lt'), op('<')),
            gt=comp_method(operator.gt, names('gt'), op('>')),
            le=comp_method(operator.le, names('le'), op('<=')),
            ge=comp_method(operator.ge, names('ge'), op('>=')),
        ))
    if bool_method:
        new_methods.update(dict(
        and_=bool_method(operator.and_, names('and_ [&]'), op('&')),
        rand_=bool_method(lambda x, y: operator.and_(y, x), names('rand_[&]')),
        or_=bool_method(operator.or_, names('or_ [|]'), op('|')),
        ror_=bool_method(lambda x, y: operator.or_(y, x), names('ror_ [|]')),
        # For some reason ``^`` wasn't used in original.
        xor=bool_method(operator.xor, names('xor [^]')),
        rxor=bool_method(lambda x, y: operator.xor(y, x), names('rxor [^]'))
        ))

    new_methods = dict((names(k), v) for k, v in new_methods.items())
    return new_methods


def add_methods(cls, new_methods, force, select, exclude):
    if select and exclude:
        raise TypeError("May only pass either select or exclude")
    methods = new_methods
    if select:
        select = set(select)
        methods = {}
        for key, method in new_methods.items():
            if key in select:
                methods[key] = method
    if exclude:
        for k in exclude:
            new_methods.pop(k, None)

    for name, method in new_methods.items():
        if force or name not in cls.__dict__:
            bind_method(cls, name, method)

#----------------------------------------------------------------------
# Arithmetic
def add_special_arithmetic_methods(cls, arith_method=None, radd_func=None,
                                   comp_method=None, bool_method=None,
                                   use_numexpr=True, force=False, select=None,
                                   exclude=None):
    """
    Adds the full suite of special arithmetic methods (``__add__``, ``__sub__``, etc.) to the class.

    Parameters
    ----------
    arith_method : function (optional)
        factory for special arithmetic methods, with op string:
        f(op, name, str_rep, default_axis=None, fill_zeros=None, **eval_kwargs)
    radd_func :  function (optional)
        Possible replacement for ``operator.add`` for compatibility
    comp_method : function, optional,
        factory for rich comparison - signature: f(op, name, str_rep)
    use_numexpr : bool, default True
        whether to accelerate with numexpr, defaults to True
    force : bool, default False
        if False, checks whether function is defined **on ``cls.__dict__``** before defining
        if True, always defines functions on class base
    select : iterable of strings (optional)
        if passed, only sets functions with names in select
    exclude : iterable of strings (optional)
        if passed, will not set functions with names in exclude
    """
    radd_func = radd_func or operator.add
    # in frame, special methods have default_axis = None, comp methods use 'columns'
    new_methods = _create_methods(arith_method, radd_func, comp_method, bool_method, use_numexpr, default_axis=None,
                                  special=True)

    # inplace operators (I feel like these should get passed an `inplace=True`
    # or just be removed
    new_methods.update(dict(
        __iadd__=new_methods["__add__"],
        __isub__=new_methods["__sub__"],
        __imul__=new_methods["__mul__"],
        __itruediv__=new_methods["__truediv__"],
        __ipow__=new_methods["__pow__"]
    ))
    if not compat.PY3:
        new_methods["__idiv__"] = new_methods["__div__"]

    add_methods(cls, new_methods=new_methods, force=force, select=select, exclude=exclude)


def add_flex_arithmetic_methods(cls, flex_arith_method, radd_func=None,
                                flex_comp_method=None, flex_bool_method=None,
                                use_numexpr=True, force=False, select=None,
                                exclude=None):
    """
    Adds the full suite of flex arithmetic methods (``pow``, ``mul``, ``add``) to the class.

    Parameters
    ----------
    flex_arith_method : function (optional)
        factory for special arithmetic methods, with op string:
        f(op, name, str_rep, default_axis=None, fill_zeros=None, **eval_kwargs)
    radd_func :  function (optional)
        Possible replacement for ``lambda x, y: operator.add(y, x)`` for compatibility
    flex_comp_method : function, optional,
        factory for rich comparison - signature: f(op, name, str_rep)
    use_numexpr : bool, default True
        whether to accelerate with numexpr, defaults to True
    force : bool, default False
        if False, checks whether function is defined **on ``cls.__dict__``** before defining
        if True, always defines functions on class base
    select : iterable of strings (optional)
        if passed, only sets functions with names in select
    exclude : iterable of strings (optional)
        if passed, will not set functions with names in exclude
    """
    radd_func = radd_func or (lambda x, y: operator.add(y, x))
    # in frame, default axis is 'columns', doesn't matter for series and panel
    new_methods = _create_methods(
        flex_arith_method, radd_func, flex_comp_method, flex_bool_method,
        use_numexpr, default_axis='columns', special=False)
    new_methods.update(dict(
        multiply=new_methods['mul'],
        subtract=new_methods['sub'],
        divide=new_methods['div']
    ))

    add_methods(cls, new_methods=new_methods, force=force, select=select, exclude=exclude)

def cleanup_name(name):
    """cleanup special names
    >>> cleanup_name("__rsub__")
    sub
    >>> cleanup_name("rand_")
    and_
    """
    if name[:2] == "__":
        name = name[2:-2]
    if name[0] == "r":
        name = name[1:]
    # readd last _ for operator names.
    if name == "or":
        name = "or_"
    elif name == "and":
        name = "and_"
    return name

# ----------------------------------------------------------------------------
# Arithmetic factory methods

def mask_no_list(obj, mask):
    return obj[mask]

def mask_with_list(obj, mask):
    return np.array(list(obj[mask]))

def _standard_na_op(x, y, op, str_rep=None, fill_zeros=None,
                    convert_mask=False, masker=None, eval_kwargs={}):
    """standard internal operation handler, generally you want to wrap it
    in a partial so it can be called with just two arguments, e.g.
    na_op = partial(_standard_na_op, op=op, str_rep=str_rep, fill_zeros=fill_zeros)

    x : array-like
        object to be operated open
    y : object (e.g., scalar or ndarray-like)
        object to perform arithmetic with
    op : binary arithmetic function function
    str_rep : str (optional)
        string representation of operation to pass to expressions.evaluate.
    fill_zeros : object (optional)
        passes to com._fill_zeros (does nothing if None)
    convert_mask : bool (default False)
        if True, wraps mask with list and does not call upcast_putmask
        if False, does not wrap mask, but calls upcast_putmask
    masker : bool (optional)
        Value to fill in for masked values. Must explicitly pass if using
        convert_mask.
    eval_kwargs : dict (optional)
        Keyword arguments to pass to expressions.evaluate. Note that this must
        be a literal dict, you can't pass as ``**eval_kwargs``
    """
    # generally used for bool methods
    if convert_mask:
        do_mask = lambda obj, mask: np.array(list(obj[mask]))
        assert masker is not None, "With convert_mask, must also specify masker"
    else:
        do_mask = lambda obj, maks: obj[mask]

    try:
        result = expressions.evaluate(op, str_rep, x, y,
                                        raise_on_error=True, **eval_kwargs)
    except TypeError:
        xrav = x.ravel()
        result = np.empty(x.size, dtype=x.dtype)
        if isinstance(y, (np.ndarray, com.ABCSeries)):
            yrav = y.ravel()
            mask = notnull(xrav) & notnull(yrav)
            result[mask] = op(do_mask(xrav, mask), do_mask(yrav, mask))
        else:
            mask = notnull(xrav)
            result[mask] = op(do_mask(xrav, mask), y)
        if convert_mask:
            if not mask.all():
                result[-mask] = masker
        else:
            result, changed = com._maybe_upcast_putmask(result, -mask, np.nan)

        result = result.reshape(x.shape)

    # handles discrepancy between numpy and numexpr on division/mod by 0
    if fill_zeros is not None:
        result = com._fill_zeros(result, y, fill_zeros)
    return result


#----------------------------------------------------------------------
# Wrapper function for Series arithmetic methods
def _arith_method_SERIES(op, name, str_rep=None, fill_zeros=None, default_axis=None, **eval_kwargs):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    r_op = name.startswith("__r")
    na_op = partial(_standard_na_op, op=op, str_rep=str_rep,
                    fill_zeros=fill_zeros, convert_mask=False,
                    eval_kwargs=eval_kwargs)

    def wrapper(self, other, name=name):
        from pandas.core.frame import DataFrame
        dtype = None
        fill_value = tslib.iNaT
        wrap_results = lambda x: x

        lvalues, rvalues = self, other

        is_timedelta_lhs = com.is_timedelta64_dtype(self)
        is_datetime_lhs  = com.is_datetime64_dtype(self)
        is_integer_lhs   = lvalues.dtype.kind in ['i','u']

        if is_datetime_lhs or is_timedelta_lhs:
            if r_op:
                raise NotImplementedError('__r*__ operations for datetimes and timedeltas are not yet supported')

            coerce = 'compat' if _np_version_under1p7 else True

            # convert the argument to an ndarray
            def convert_to_array(values):
                if not is_list_like(values):
                    values = np.array([values])
                inferred_type = lib.infer_dtype(values)
                if inferred_type in set(['datetime64','datetime','date','time']):
                    # a datetlike
                    if not (isinstance(values, pa.Array) and com.is_datetime64_dtype(values)):
                        values = tslib.array_to_datetime(values)
                elif inferred_type in set(['timedelta']):
                    # have a timedelta, convert to to ns here
                    values = com._possibly_cast_to_timedelta(values, coerce=coerce)
                elif inferred_type in set(['timedelta64']):
                    # have a timedelta64, make sure dtype dtype is ns
                    values = com._possibly_cast_to_timedelta(values, coerce=coerce)
                elif inferred_type in set(['integer']):
                    # py3 compat where dtype is 'm' but is an integer
                    if values.dtype.kind == 'm':
                        values = values.astype('timedelta64[ns]')
                    elif name not in ['truediv','div','mul']:
                        raise TypeError("incompatible type for a datetime/timedelta operation [{0}]".format(name))
                elif isinstance(values[0],DateOffset):
                    # handle DateOffsets
                    os = pa.array([ getattr(v,'delta',None) for v in values ])
                    mask = isnull(os)
                    if mask.any():
                        raise TypeError("cannot use a non-absolute DateOffset in "
                                        "datetime/timedelta operations [{0}]".format(','.join([ com.pprint_thing(v) for v in values[mask] ])))
                    values = com._possibly_cast_to_timedelta(os, coerce=coerce)
                else:
                    raise TypeError("incompatible type [{0}] for a datetime/timedelta operation".format(pa.array(values).dtype))

                return values

            # convert lhs and rhs
            lvalues = convert_to_array(lvalues)
            rvalues = convert_to_array(rvalues)

            is_datetime_rhs  = com.is_datetime64_dtype(rvalues)
            is_timedelta_rhs = com.is_timedelta64_dtype(rvalues) or (not is_datetime_rhs and _np_version_under1p7)
            is_integer_rhs = rvalues.dtype.kind in ['i','u']
            mask = None

            # timedelta and integer mul/div
            if (is_timedelta_lhs and is_integer_rhs) or (is_integer_lhs and is_timedelta_rhs):

                if name not in ['truediv','div','mul']:
                    raise TypeError("can only operate on a timedelta and an integer for "
                                    "division, but the operator [%s] was passed" % name)
                dtype = 'timedelta64[ns]'
                mask = isnull(lvalues) | isnull(rvalues)
                lvalues = lvalues.astype(np.int64)
                rvalues = rvalues.astype(np.int64)

            # 2 datetimes
            elif is_datetime_lhs and is_datetime_rhs:
                if name != 'sub':
                    raise TypeError("can only operate on a datetimes for subtraction, "
                                    "but the operator [%s] was passed" % name)

                dtype = 'timedelta64[ns]'
                mask = isnull(lvalues) | isnull(rvalues)
                lvalues = lvalues.view('i8')
                rvalues = rvalues.view('i8')

            # 2 timedeltas
            elif is_timedelta_lhs and is_timedelta_rhs:
                mask = isnull(lvalues) | isnull(rvalues)

                # time delta division -> unit less
                if name in ['div','truediv']:
                    dtype = 'float64'
                    fill_value = np.nan
                    lvalues = lvalues.astype(np.int64).astype(np.float64)
                    rvalues = rvalues.astype(np.int64).astype(np.float64)

                # another timedelta
                elif name in ['add','sub']:
                    dtype = 'timedelta64[ns]'
                    lvalues = lvalues.astype(np.int64)
                    rvalues = rvalues.astype(np.int64)

                else:
                    raise TypeError("can only operate on a timedeltas for "
                                    "addition, subtraction, and division, but the operator [%s] was passed" % name)

            # datetime and timedelta
            elif is_timedelta_rhs and is_datetime_lhs:

                if name not in ['add','sub']:
                    raise TypeError("can only operate on a datetime with a rhs of a timedelta for "
                                    "addition and subtraction, but the operator [%s] was passed" % name)
                dtype = 'M8[ns]'
                lvalues = lvalues.view('i8')
                rvalues = rvalues.view('i8')

            elif is_timedelta_lhs and is_datetime_rhs:

                if name not in ['add']:
                    raise TypeError("can only operate on a timedelta and a datetime for "
                                    "addition, but the operator [%s] was passed" % name)
                dtype = 'M8[ns]'
                lvalues = lvalues.view('i8')
                rvalues = rvalues.view('i8')

            else:
                raise TypeError('cannot operate on a series with out a rhs '
                                'of a series/ndarray of type datetime64[ns] '
                                'or a timedelta')

            # if we need to mask the results
            if mask is not None:
                if mask.any():
                    def f(x):
                        x = pa.array(x,dtype=dtype)
                        np.putmask(x,mask,fill_value)
                        return x
                    wrap_results = f

        if isinstance(rvalues, com.ABCSeries):
            lvalues = getattr(lvalues, 'values', lvalues)
            rvalues = getattr(rvalues, 'values', rvalues)

            if self.index.equals(other.index):
                name = _maybe_match_name(self, other)
                return self._constructor(wrap_results(na_op(lvalues, rvalues)),
                                         index=self.index, dtype=dtype, name=name)

            join_idx, lidx, ridx = self.index.join(other.index, how='outer',
                                                   return_indexers=True)

            if lidx is not None:
                lvalues = com.take_1d(lvalues, lidx)

            if ridx is not None:
                rvalues = com.take_1d(rvalues, ridx)

            arr = na_op(lvalues, rvalues)

            name = _maybe_match_name(self, other)
            return self._constructor(wrap_results(arr), index=join_idx, name=name, dtype=dtype)
        elif isinstance(other, com.ABCDataFrame):
            return NotImplemented
        else:
            # scalars
            if hasattr(lvalues, 'values'):
                lvalues = lvalues.values
            return self._constructor(wrap_results(na_op(lvalues, rvalues)),
                                     index=self.index, name=self.name, dtype=dtype)
    wrapper.__name__ = name
    # handle r* methods, etc.
    name = cleanup_name(name)
    return wrapper


def _comp_method_SERIES(op, name, str_rep=None, masker=False):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def na_op(x, y):
        if x.dtype == np.object_:
            if isinstance(y, list):
                y = lib.list_to_object_array(y)

            if isinstance(y, (pa.Array, com.ABCSeries)):
                if y.dtype != np.object_:
                    result = lib.vec_compare(x, y.astype(np.object_), op)
                else:
                    result = lib.vec_compare(x, y, op)
            else:
                result = lib.scalar_compare(x, y, op)
        else:
            result = op(x, y)

        return result

    def wrapper(self, other):
        if isinstance(other, com.ABCSeries):
            name = _maybe_match_name(self, other)
            if len(self) != len(other):
                raise ValueError('Series lengths must match to compare')
            return self._constructor(na_op(self.values, other.values),
                                     index=self.index, name=name)
        elif isinstance(other, com.ABCDataFrame):  # pragma: no cover
            return NotImplemented
        elif isinstance(other, (pa.Array, com.ABCSeries)):
            if len(self) != len(other):
                raise ValueError('Lengths must match to compare')
            return self._constructor(na_op(self.values, np.asarray(other)),
                                     index=self.index, name=self.name)
        else:
            from pandas.core.series import Series

            mask = isnull(self)

            values = self.values
            other = _index.convert_scalar(values, other)

            if issubclass(values.dtype.type, np.datetime64):
                values = values.view('i8')

            # scalars
            res = na_op(values, other)
            if np.isscalar(res):
                raise TypeError('Could not compare %s type with Series'
                                % type(other))
            # always return a full value series here
            res = _values_from_object(res)

            res = Series(res, index=self.index, name=self.name, dtype='bool')

            # mask out the invalids
            if mask.any():
                res[mask.values] = masker

            return res

    wrapper.__name__ = name
    return wrapper


def _bool_method_SERIES(op, name, str_rep=None):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            if isinstance(y, list):
                y = lib.list_to_object_array(y)

            if isinstance(y, (pa.Array, com.ABCSeries)):
                if (x.dtype == np.bool_ and
                        y.dtype == np.bool_):  # pragma: no cover
                    result = op(x, y)  # when would this be hit?
                else:
                    x = com._ensure_object(x)
                    y = com._ensure_object(y)
                    result = lib.vec_binop(x, y, op)
            else:
                result = lib.scalar_binop(x, y, op)

        return result

    def wrapper(self, other):
        if isinstance(other, com.ABCSeries):
            name = _maybe_match_name(self, other)
            return self._constructor(na_op(self.values, other.values),
                                     index=self.index, name=name)
        elif isinstance(other, com.ABCDataFrame):
            return NotImplemented
        else:
            # scalars
            return self._constructor(na_op(self.values, other),
                                     index=self.index, name=self.name)

    wrapper.__name__ = name
    return wrapper


def _radd_compat_SERIES(left, right):
    radd = lambda x, y: y + x
    # GH #353, NumPy 1.5.1 workaround
    try:
        output = radd(left, right)
    except TypeError:
        cond = (_np_version_under1p6 and
                left.dtype == np.object_)
        if cond:  # pragma: no cover
            output = np.empty_like(left)
            output.flat[:] = [radd(x, right) for x in left.flat]
        else:
            raise

    return output


def _flex_method_SERIES(op, name, str_rep=None, default_axis=None, fill_zeros=None, **eval_kwargs):
    doc = """
    Binary operator %s with support to substitute a fill_value for missing data
    in one of the inputs

    Parameters
    ----------
    other: Series or scalar value
    fill_value : None or float value, default None (NaN)
        Fill missing (NaN) values with this value. If both Series are
        missing, the result will be missing
    level : int or name
        Broadcast across a level, matching Index values on the
        passed MultiIndex level

    Returns
    -------
    result : Series
    """ % name
    na_op = partial(_standard_na_op, op=op, str_rep=str_rep,
                    fill_zeros=fill_zeros, convert_mask=False,
                    eval_kwargs=eval_kwargs)

    @Appender(doc)
    def f(self, other, level=None, fill_value=None):
        if isinstance(other, com.ABCSeries):
            return self._binop(other, na_op, level=level, fill_value=fill_value)
        elif isinstance(other, (pa.Array, com.ABCSeries, list, tuple)):
            if len(other) != len(self):
                raise ValueError('Lengths must be equal')
            return self._binop(self._constructor(other, self.index), na_op,
                               level=level, fill_value=fill_value)
        else:
            return self._constructor(na_op(self.values, other), self.index,
                                     name=self.name)

    f.__name__ = name
    return f

series_flex_funcs = dict(flex_arith_method=_flex_method_SERIES,
                         radd_func=_radd_compat_SERIES,
                         flex_comp_method=_comp_method_SERIES)
series_special_funcs = dict(arith_method=_arith_method_SERIES,
                            radd_func=_radd_compat_SERIES,
                            comp_method=_comp_method_SERIES,
                            bool_method=_bool_method_SERIES)

_arith_doc_FRAME = """
Binary operator %s with support to substitute a fill_value for missing data in
one of the inputs

Parameters
----------
other : Series, DataFrame, or constant
axis : {0, 1, 'index', 'columns'}
    For Series input, axis to match Series index on
fill_value : None or float value, default None
    Fill missing (NaN) values with this value. If both DataFrame locations are
    missing, the result will be missing
level : int or name
    Broadcast across a level, matching Index values on the
    passed MultiIndex level

Notes
-----
Mismatched indices will be unioned together

Returns
-------
result : DataFrame
"""


def _arith_method_FRAME(op, name, str_rep=None, default_axis='columns', fill_zeros=None, **eval_kwargs):
    na_op = partial(_standard_na_op, op=op, str_rep=str_rep,
                    fill_zeros=fill_zeros, convert_mask=False,
                    eval_kwargs=eval_kwargs)

    @Appender(_arith_doc_FRAME % name)
    def f(self, other, axis=default_axis, level=None, fill_value=None):
        if isinstance(other, com.ABCDataFrame):    # Another DataFrame
            return self._combine_frame(other, na_op, fill_value, level)
        elif isinstance(other, com.ABCSeries):
            return self._combine_series(other, na_op, fill_value, axis, level)
        elif isinstance(other, (list, tuple)):
            if axis is not None and self._get_axis_name(axis) == 'index':
                casted = self._constructor_sliced(other, index=self.index)
            else:
                casted = self._constructor_sliced(other, index=self.columns)
            return self._combine_series(casted, na_op, fill_value, axis, level)
        elif isinstance(other, np.ndarray):
            if other.ndim == 1:
                if axis is not None and self._get_axis_name(axis) == 'index':
                    casted = self._constructor_sliced(other, index=self.index)
                else:
                    casted = self._constructor_sliced(other, index=self.columns)
                return self._combine_series(casted, na_op, fill_value,
                                            axis, level)
            elif other.ndim == 2:
                casted = self._constructor(other, index=self.index,
                                           columns=self.columns)
                return self._combine_frame(casted, na_op, fill_value, level)
            else:  # pragma: no cover
                raise ValueError("Bad argument shape")
        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


def _flex_comp_method_FRAME(op, name, str_rep=None, default_axis='columns',
                            masker=False):

    na_op = partial(_standard_na_op, op=op, str_rep=str_rep, convert_mask=True,
                    masker=masker)

    @Appender('Wrapper for flexible comparison methods %s' % name)
    def f(self, other, axis=default_axis, level=None):
        if isinstance(other, com.ABCDataFrame):    # Another DataFrame
            return self._flex_compare_frame(other, na_op, str_rep, level)

        elif isinstance(other, com.ABCSeries):
            return self._combine_series(other, na_op, None, axis, level)

        elif isinstance(other, (list, tuple)):
            if axis is not None and self._get_axis_name(axis) == 'index':
                casted = self._constructor_sliced(other, index=self.index)
            else:
                casted = self._constructor_sliced(other, index=self.columns)

            return self._combine_series(casted, na_op, None, axis, level)

        elif isinstance(other, np.ndarray):
            if other.ndim == 1:
                if axis is not None and self._get_axis_name(axis) == 'index':
                    casted = self._constructor_sliced(other, index=self.index)
                else:
                    casted = self._constructor_sliced(other, index=self.columns)

                return self._combine_series(casted, na_op, None, axis, level)

            elif other.ndim == 2:
                casted = self._constructor(other, index=self.index,
                                           columns=self.columns)

                return self._flex_compare_frame(casted, na_op, str_rep, level)

            else:  # pragma: no cover
                raise ValueError("Bad argument shape")

        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


def _comp_method_FRAME(func, name, str_rep, masker=False):
    @Appender('Wrapper for comparison method %s' % name)
    def f(self, other):
        if isinstance(other, com.ABCDataFrame):    # Another DataFrame
            return self._compare_frame(other, func, str_rep)
        elif isinstance(other, com.ABCSeries):
            return self._combine_series_infer(other, func)
        else:

            # straight boolean comparisions we want to allow all columns
            # (regardless of dtype to pass thru) See #4537 for discussion.
            return self._combine_const(other, func, raise_on_error=False).fillna(True).astype(bool)

    f.__name__ = name

    return f

frame_flex_funcs = dict(flex_arith_method=_arith_method_FRAME,
                        radd_func=_radd_compat_SERIES,
                        flex_comp_method=_flex_comp_method_FRAME)

frame_special_funcs = dict(arith_method=_arith_method_FRAME,
                           radd_func=_radd_compat_SERIES,
                           comp_method=_comp_method_FRAME,
                           bool_method=_arith_method_FRAME)

def _arith_method_PANEL(op, name, str_rep=None, fill_zeros=None, default_axis=None, **eval_kwargs):
    # work only for scalars
    def na_op(x, y):
        try:
            result = expressions.evaluate(op, str_rep, x, y, raise_on_error=True, **eval_kwargs)
        except TypeError:
            result = op(x, y)

        # handles discrepancy between numpy and numexpr on division/mod by 0
        result = com._fill_zeros(result,y,fill_zeros)
        return result

    def f(self, other):
        if not np.isscalar(other):
            raise ValueError('Simple arithmetic with %s can only be '
                             'done with scalar values' % self._constructor.__name__)

        return self._combine(other, na_op)
    f.__name__ = name
    return f


def _comp_method_PANEL(op, name, str_rep=None, masker=False):
    na_op = partial(_standard_na_op, op=op, str_rep=str_rep, convert_mask=True,
                    masker=masker)

    @Appender('Wrapper for comparison method %s' % name)
    def f(self, other):
        if isinstance(other, self._constructor):
            return self._compare_constructor(other, op)
        elif isinstance(other, (self._constructor_sliced, com.ABCDataFrame, com.ABCSeries)):
            raise Exception("input needs alignment for this object [%s]" %
                            self._constructor)
        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f

panel_special_funcs = dict(arith_method=_arith_method_PANEL,
                           comp_method=_comp_method_PANEL,
                           bool_method=_arith_method_PANEL)
