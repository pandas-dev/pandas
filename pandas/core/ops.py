"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
# necessary to enforce truediv in Python 2.X
from __future__ import division
import operator
import numpy as np
import pandas as pd
from pandas import compat, lib, tslib
import pandas.index as _index
from pandas.util.decorators import Appender
import pandas.core.common as com
import pandas.core.array as pa
import pandas.computation.expressions as expressions
from pandas.core.common import(bind_method, is_list_like, notnull, isnull,
                               _values_from_object, _maybe_match_name)

# -----------------------------------------------------------------------------
# Functions that add arithmetic methods to objects, given arithmetic factory
# methods


def _create_methods(arith_method, radd_func, comp_method, bool_method,
                    use_numexpr, special=False, default_axis='columns'):
    # creates actual methods based upon arithmetic, comp and bool method
    # constructors.

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
    # Inframe, all special methods have default_axis=None, flex methods have
    # default_axis set to the default (columns)
    new_methods = dict(
        add=arith_method(operator.add, names('add'), op('+'),
                         default_axis=default_axis),
        radd=arith_method(radd_func, names('radd'), op('+'),
                          default_axis=default_axis),
        sub=arith_method(operator.sub, names('sub'), op('-'),
                         default_axis=default_axis),
        mul=arith_method(operator.mul, names('mul'), op('*'),
                         default_axis=default_axis),
        truediv=arith_method(operator.truediv, names('truediv'), op('/'),
                             truediv=True, fill_zeros=np.inf,
                             default_axis=default_axis),
        floordiv=arith_method(operator.floordiv, names('floordiv'), op('//'),
                              default_axis=default_axis, fill_zeros=np.inf),
        # Causes a floating point exception in the tests when numexpr
        # enabled, so for now no speedup
        mod=arith_method(operator.mod, names('mod'), default_axis=default_axis,
                         fill_zeros=np.nan),
        pow=arith_method(operator.pow, names('pow'), op('**'),
                         default_axis=default_axis),
        # not entirely sure why this is necessary, but previously was included
        # so it's here to maintain compatibility
        rmul=arith_method(operator.mul, names('rmul'),
                          default_axis=default_axis),
        rsub=arith_method(lambda x, y: y - x, names('rsub'),
                          default_axis=default_axis),
        rtruediv=arith_method(lambda x, y: operator.truediv(y, x),
                              names('rtruediv'), truediv=True,
                              fill_zeros=np.inf, default_axis=default_axis),
        rfloordiv=arith_method(lambda x, y: operator.floordiv(y, x),
                               names('rfloordiv'), default_axis=default_axis,
                               fill_zeros=np.inf),
        rpow=arith_method(lambda x, y: y ** x, names('rpow'),
                          default_axis=default_axis),
        rmod=arith_method(lambda x, y: y % x, names('rmod'),
                          default_axis=default_axis),
    )
    new_methods['div'] = new_methods['truediv']
    new_methods['rdiv'] = new_methods['rtruediv']

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
            and_=bool_method(operator.and_, names('and_'), op('&')),
            or_=bool_method(operator.or_, names('or_'), op('|')),
            # For some reason ``^`` wasn't used in original.
            xor=bool_method(operator.xor, names('xor')),
            rand_=bool_method(lambda x, y: operator.and_(y, x),
                              names('rand_')),
            ror_=bool_method(lambda x, y: operator.or_(y, x), names('ror_')),
            rxor=bool_method(lambda x, y: operator.xor(y, x), names('rxor'))
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
    Adds the full suite of special arithmetic methods (``__add__``,
    ``__sub__``, etc.) to the class.

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
        if False, checks whether function is defined **on ``cls.__dict__``**
        before defining if True, always defines functions on class base
    select : iterable of strings (optional)
        if passed, only sets functions with names in select
    exclude : iterable of strings (optional)
        if passed, will not set functions with names in exclude
    """
    radd_func = radd_func or operator.add
    # in frame, special methods have default_axis = None, comp methods use
    # 'columns'
    new_methods = _create_methods(arith_method, radd_func, comp_method,
                                  bool_method, use_numexpr, default_axis=None,
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

    add_methods(cls, new_methods=new_methods, force=force, select=select,
                exclude=exclude)


def add_flex_arithmetic_methods(cls, flex_arith_method, radd_func=None,
                                flex_comp_method=None, flex_bool_method=None,
                                use_numexpr=True, force=False, select=None,
                                exclude=None):
    """
    Adds the full suite of flex arithmetic methods (``pow``, ``mul``, ``add``)
    to the class.

    Parameters
    ----------
    flex_arith_method : function (optional)
        factory for special arithmetic methods, with op string:
        f(op, name, str_rep, default_axis=None, fill_zeros=None, **eval_kwargs)
    radd_func :  function (optional)
        Possible replacement for ``lambda x, y: operator.add(y, x)`` for
        compatibility
    flex_comp_method : function, optional,
        factory for rich comparison - signature: f(op, name, str_rep)
    use_numexpr : bool, default True
        whether to accelerate with numexpr, defaults to True
    force : bool, default False
        if False, checks whether function is defined **on ``cls.__dict__``**
        before defining if True, always defines functions on class base
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
    # opt out of bool flex methods for now
    for k in ('ror_', 'rxor', 'rand_'):
        if k in new_methods:
            new_methods.pop(k)

    add_methods(cls, new_methods=new_methods, force=force, select=select,
                exclude=exclude)


class _TimeOp(object):

    """
    Wrapper around Series datetime/time/timedelta arithmetic operations.
    Generally, you should use classmethod ``maybe_convert_for_time_op`` as an
    entry point.
    """
    fill_value = tslib.iNaT
    wrap_results = staticmethod(lambda x: x)
    dtype = None

    def __init__(self, left, right, name):
        self.name = name

        lvalues = self._convert_to_array(left, name=name)
        rvalues = self._convert_to_array(right, name=name, other=lvalues)

        self.is_timedelta_lhs = com.is_timedelta64_dtype(left)
        self.is_datetime_lhs = com.is_datetime64_dtype(left)
        self.is_integer_lhs = left.dtype.kind in ['i', 'u']
        self.is_datetime_rhs = com.is_datetime64_dtype(rvalues)
        self.is_timedelta_rhs = (com.is_timedelta64_dtype(rvalues)
                                 or (not self.is_datetime_rhs
                                     and pd._np_version_under1p7))
        self.is_integer_rhs = rvalues.dtype.kind in ('i', 'u')

        self._validate()

        self._convert_for_datetime(lvalues, rvalues)

    def _validate(self):
        # timedelta and integer mul/div

        if (self.is_timedelta_lhs and self.is_integer_rhs) or\
           (self.is_integer_lhs and self.is_timedelta_rhs):

            if self.name not in ('__truediv__', '__div__', '__mul__'):
                raise TypeError("can only operate on a timedelta and an "
                                "integer for division, but the operator [%s]"
                                "was passed" % self.name)

        # 2 datetimes
        elif self.is_datetime_lhs and self.is_datetime_rhs:
            if self.name != '__sub__':
                raise TypeError("can only operate on a datetimes for"
                                " subtraction, but the operator [%s] was"
                                " passed" % self.name)

        # 2 timedeltas
        elif self.is_timedelta_lhs and self.is_timedelta_rhs:

            if self.name not in ('__div__', '__truediv__', '__add__',
                                 '__sub__'):
                raise TypeError("can only operate on a timedeltas for "
                                "addition, subtraction, and division, but the"
                                " operator [%s] was passed" % self.name)

        # datetime and timedelta
        elif self.is_datetime_lhs and self.is_timedelta_rhs:

            if self.name not in ('__add__', '__sub__'):
                raise TypeError("can only operate on a datetime with a rhs of"
                                " a timedelta for addition and subtraction, "
                                " but the operator [%s] was passed" %
                                self.name)

        elif self.is_timedelta_lhs and self.is_datetime_rhs:

            if self.name != '__add__':
                raise TypeError("can only operate on a timedelta and"
                                " a datetime for addition, but the operator"
                                " [%s] was passed" % self.name)
        else:
            raise TypeError('cannot operate on a series with out a rhs '
                            'of a series/ndarray of type datetime64[ns] '
                            'or a timedelta')

    def _convert_to_array(self, values, name=None, other=None):
        """converts values to ndarray"""
        from pandas.tseries.timedeltas import _possibly_cast_to_timedelta

        coerce = 'compat' if pd._np_version_under1p7 else True
        if not is_list_like(values):
            values = np.array([values])
        inferred_type = lib.infer_dtype(values)

        if inferred_type in ('datetime64', 'datetime', 'date', 'time'):
            # if we have a other of timedelta, but use pd.NaT here we
            # we are in the wrong path
            if (other is not None and other.dtype == 'timedelta64[ns]' and
                    all(isnull(v) for v in values)):
                values = np.empty(values.shape, dtype=other.dtype)
                values[:] = tslib.iNaT

            # a datetlike
            elif not (isinstance(values, (pa.Array, pd.Series)) and
                      com.is_datetime64_dtype(values)):
                values = tslib.array_to_datetime(values)
            elif isinstance(values, pd.DatetimeIndex):
                values = values.to_series()
        elif inferred_type in ('timedelta', 'timedelta64'):
            # have a timedelta, convert to to ns here
            values = _possibly_cast_to_timedelta(values, coerce=coerce)
        elif inferred_type == 'integer':
            # py3 compat where dtype is 'm' but is an integer
            if values.dtype.kind == 'm':
                values = values.astype('timedelta64[ns]')
            elif isinstance(values, pd.PeriodIndex):
                values = values.to_timestamp().to_series()
            elif name not in ('__truediv__', '__div__', '__mul__'):
                raise TypeError("incompatible type for a datetime/timedelta "
                                "operation [{0}]".format(name))
        elif isinstance(values[0], pd.DateOffset):
            # handle DateOffsets
            os = pa.array([getattr(v, 'delta', None) for v in values])
            mask = isnull(os)
            if mask.any():
                raise TypeError("cannot use a non-absolute DateOffset in "
                                "datetime/timedelta operations [{0}]".format(
                                    ', '.join([com.pprint_thing(v)
                                               for v in values[mask]])))
            values = _possibly_cast_to_timedelta(os, coerce=coerce)
        elif inferred_type == 'floating':

            # all nan, so ok, use the other dtype (e.g. timedelta or datetime)
            if isnull(values).all():
                values = np.empty(values.shape, dtype=other.dtype)
                values[:] = tslib.iNaT
            else:
                raise TypeError(
                    'incompatible type [{0}] for a datetime/timedelta '
                    'operation'.format(pa.array(values).dtype))
        else:
            raise TypeError("incompatible type [{0}] for a datetime/timedelta"
                            " operation".format(pa.array(values).dtype))

        return values

    def _convert_for_datetime(self, lvalues, rvalues):
        mask = None
        # datetimes require views
        if self.is_datetime_lhs or self.is_datetime_rhs:
            # datetime subtraction means timedelta
            if self.is_datetime_lhs and self.is_datetime_rhs:
                self.dtype = 'timedelta64[ns]'
            else:
                self.dtype = 'datetime64[ns]'
            mask = isnull(lvalues) | isnull(rvalues)
            lvalues = lvalues.view(np.int64)
            rvalues = rvalues.view(np.int64)

        # otherwise it's a timedelta
        else:
            self.dtype = 'timedelta64[ns]'
            mask = isnull(lvalues) | isnull(rvalues)
            lvalues = lvalues.astype(np.int64)
            rvalues = rvalues.astype(np.int64)

            # time delta division -> unit less
            # integer gets converted to timedelta in np < 1.6
            if (self.is_timedelta_lhs and self.is_timedelta_rhs) and\
               not self.is_integer_rhs and\
               not self.is_integer_lhs and\
               self.name in ('__div__', '__truediv__'):
                self.dtype = 'float64'
                self.fill_value = np.nan
                lvalues = lvalues.astype(np.float64)
                rvalues = rvalues.astype(np.float64)

        # if we need to mask the results
        if mask is not None:
            if mask.any():
                def f(x):
                    x = pa.array(x, dtype=self.dtype)
                    np.putmask(x, mask, self.fill_value)
                    return x
                self.wrap_results = f
        self.lvalues = lvalues
        self.rvalues = rvalues

    @classmethod
    def maybe_convert_for_time_op(cls, left, right, name):
        """
        if ``left`` and ``right`` are appropriate for datetime arithmetic with
        operation ``name``, processes them and returns a ``_TimeOp`` object
        that stores all the required values.  Otherwise, it will generate
        either a ``NotImplementedError`` or ``None``, indicating that the
        operation is unsupported for datetimes (e.g., an unsupported r_op) or
        that the data is not the right type for time ops.
        """
        # decide if we can do it
        is_timedelta_lhs = com.is_timedelta64_dtype(left)
        is_datetime_lhs = com.is_datetime64_dtype(left)
        if not (is_datetime_lhs or is_timedelta_lhs):
            return None
        # rops are allowed. No need for special checks, just strip off
        # r part.
        if name.startswith('__r'):
            name = "__" + name[3:]
        return cls(left, right, name)


def _arith_method_SERIES(op, name, str_rep=None, fill_zeros=None,
                         default_axis=None, **eval_kwargs):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def na_op(x, y):
        try:
            result = expressions.evaluate(op, str_rep, x, y,
                                          raise_on_error=True, **eval_kwargs)
        except TypeError:
            if isinstance(y, (pa.Array, pd.Series)):
                dtype = np.find_common_type([x.dtype, y.dtype], [])
                result = np.empty(x.size, dtype=dtype)
                mask = notnull(x) & notnull(y)
                result[mask] = op(x[mask], y[mask])
            else:
                result = pa.empty(len(x), dtype=x.dtype)
                mask = notnull(x)
                result[mask] = op(x[mask], y)

            result, changed = com._maybe_upcast_putmask(result, -mask, pa.NA)

        result = com._fill_zeros(result, y, fill_zeros)
        return result

    def wrapper(left, right, name=name):

        if isinstance(right, pd.DataFrame):
            return NotImplemented
        time_converted = _TimeOp.maybe_convert_for_time_op(left, right, name)

        if time_converted is None:
            lvalues, rvalues = left, right
            dtype = None
            wrap_results = lambda x: x
        elif time_converted == NotImplemented:
            return NotImplemented
        else:
            lvalues = time_converted.lvalues
            rvalues = time_converted.rvalues
            dtype = time_converted.dtype
            wrap_results = time_converted.wrap_results

        if isinstance(rvalues, pd.Series):
            join_idx, lidx, ridx = left.index.join(rvalues.index, how='outer',
                                                   return_indexers=True)
            rindex = rvalues.index
            name = _maybe_match_name(left, rvalues)
            lvalues = getattr(lvalues, 'values', lvalues)
            rvalues = getattr(rvalues, 'values', rvalues)
            if left.index.equals(rindex):
                index = left.index
            else:
                index = join_idx

                if lidx is not None:
                    lvalues = com.take_1d(lvalues, lidx)

                if ridx is not None:
                    rvalues = com.take_1d(rvalues, ridx)

            arr = na_op(lvalues, rvalues)

            return left._constructor(wrap_results(arr), index=index,
                                     name=name, dtype=dtype)
        else:
            # scalars
            if hasattr(lvalues, 'values'):
                lvalues = lvalues.values
            return left._constructor(wrap_results(na_op(lvalues, rvalues)),
                                     index=left.index, name=left.name,
                                     dtype=dtype)
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

            if isinstance(y, (pa.Array, pd.Series)):
                if y.dtype != np.object_:
                    result = lib.vec_compare(x, y.astype(np.object_), op)
                else:
                    result = lib.vec_compare(x, y, op)
            else:
                result = lib.scalar_compare(x, y, op)
        else:

            try:
                result = getattr(x, name)(y)
                if result is NotImplemented:
                    raise TypeError("invalid type comparison")
            except (AttributeError):
                result = op(x, y)

        return result

    def wrapper(self, other):
        if isinstance(other, pd.Series):
            name = _maybe_match_name(self, other)
            if len(self) != len(other):
                raise ValueError('Series lengths must match to compare')
            return self._constructor(na_op(self.values, other.values),
                                     index=self.index, name=name)
        elif isinstance(other, pd.DataFrame):  # pragma: no cover
            return NotImplemented
        elif isinstance(other, (pa.Array, pd.Series)):
            if len(self) != len(other):
                raise ValueError('Lengths must match to compare')
            return self._constructor(na_op(self.values, np.asarray(other)),
                                     index=self.index).__finalize__(self)
        else:

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

            res = pd.Series(res, index=self.index, name=self.name,
                            dtype='bool')

            # mask out the invalids
            if mask.any():
                res[mask] = masker

            return res
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

            if isinstance(y, (pa.Array, pd.Series)):
                if (x.dtype == np.bool_ and
                        y.dtype == np.bool_):  # pragma: no cover
                    result = op(x, y)  # when would this be hit?
                else:
                    x = com._ensure_object(x)
                    y = com._ensure_object(y)
                    result = lib.vec_binop(x, y, op)
            else:
                try:

                    # let null fall thru
                    if not isnull(y):
                        y = bool(y)
                    result = lib.scalar_binop(x, y, op)
                except:
                    raise TypeError("cannot compare a dtyped [{0}] array with "
                                    "a scalar of type [{1}]".format(
                                        x.dtype, type(y).__name__))

        return result

    def wrapper(self, other):
        if isinstance(other, pd.Series):
            name = _maybe_match_name(self, other)

            other = other.reindex_like(self).fillna(False).astype(bool)
            return self._constructor(na_op(self.values, other.values),
                                     index=self.index,
                                     name=name).fillna(False).astype(bool)
        elif isinstance(other, pd.DataFrame):
            return NotImplemented
        else:
            # scalars
            res = self._constructor(na_op(self.values, other),
                                    index=self.index).fillna(False)
            return res.astype(bool).__finalize__(self)
    return wrapper


def _radd_compat(left, right):
    radd = lambda x, y: y + x
    # GH #353, NumPy 1.5.1 workaround
    try:
        output = radd(left, right)
    except TypeError:
        cond = (pd._np_version_under1p6 and
                left.dtype == np.object_)
        if cond:  # pragma: no cover
            output = np.empty_like(left)
            output.flat[:] = [radd(x, right) for x in left.flat]
        else:
            raise

    return output


def _flex_method_SERIES(op, name, str_rep=None, default_axis=None,
                        fill_zeros=None, **eval_kwargs):
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

    @Appender(doc)
    def flex_wrapper(self, other, level=None, fill_value=None, axis=0):
        # validate axis
        self._get_axis_number(axis)
        if isinstance(other, pd.Series):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (pa.Array, pd.Series, list, tuple)):
            if len(other) != len(self):
                raise ValueError('Lengths must be equal')
            return self._binop(self._constructor(other, self.index), op,
                               level=level, fill_value=fill_value)
        else:
            return self._constructor(op(self.values, other),
                                     self.index).__finalize__(self)

    flex_wrapper.__name__ = name
    return flex_wrapper

series_flex_funcs = dict(flex_arith_method=_flex_method_SERIES,
                         radd_func=_radd_compat,
                         flex_comp_method=_comp_method_SERIES)

series_special_funcs = dict(arith_method=_arith_method_SERIES,
                            radd_func=_radd_compat,
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


def _arith_method_FRAME(op, name, str_rep=None, default_axis='columns',
                        fill_zeros=None, **eval_kwargs):
    def na_op(x, y):
        try:
            result = expressions.evaluate(
                op, str_rep, x, y, raise_on_error=True, **eval_kwargs)
        except TypeError:
            xrav = x.ravel()
            if isinstance(y, (np.ndarray, pd.Series)):
                dtype = np.find_common_type([x.dtype, y.dtype], [])
                result = np.empty(x.size, dtype=dtype)
                yrav = y.ravel()
                mask = notnull(xrav) & notnull(yrav)
                xrav = xrav[mask]
                yrav = yrav[mask]
                if np.prod(xrav.shape) and np.prod(yrav.shape):
                    result[mask] = op(xrav, yrav)
            else:
                result = np.empty(x.size, dtype=x.dtype)
                mask = notnull(xrav)
                xrav = xrav[mask]
                if np.prod(xrav.shape):
                    result[mask] = op(xrav, y)

            result, changed = com._maybe_upcast_putmask(result, -mask, np.nan)
            result = result.reshape(x.shape)

        result = com._fill_zeros(result, y, fill_zeros)

        return result

    @Appender(_arith_doc_FRAME % name)
    def f(self, other, axis=default_axis, level=None, fill_value=None):
        if isinstance(other, pd.DataFrame):    # Another DataFrame
            return self._combine_frame(other, na_op, fill_value, level)
        elif isinstance(other, pd.Series):
            return self._combine_series(other, na_op, fill_value, axis, level)
        elif isinstance(other, (list, tuple)):
            if axis is not None and self._get_axis_name(axis) == 'index':
                # TODO: Get all of these to use _constructor_sliced
                # casted = self._constructor_sliced(other, index=self.index)
                casted = pd.Series(other, index=self.index)
            else:
                # casted = self._constructor_sliced(other, index=self.columns)
                casted = pd.Series(other, index=self.columns)
            return self._combine_series(casted, na_op, fill_value, axis, level)
        elif isinstance(other, np.ndarray):
            if other.ndim == 1:
                if axis is not None and self._get_axis_name(axis) == 'index':
                    # casted = self._constructor_sliced(other,
                    #                                   index=self.index)
                    casted = pd.Series(other, index=self.index)
                else:
                    # casted = self._constructor_sliced(other,
                    #                                   index=self.columns)
                    casted = pd.Series(other, index=self.columns)
                return self._combine_series(casted, na_op, fill_value,
                                            axis, level)
            elif other.ndim == 2:
                # casted = self._constructor(other, index=self.index,
                #                            columns=self.columns)
                casted = pd.DataFrame(other, index=self.index,
                                      columns=self.columns)
                return self._combine_frame(casted, na_op, fill_value, level)
            else:
                raise ValueError("Incompatible argument shape: %s" %
                                 (other.shape, ))
        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


# Masker unused for now
def _flex_comp_method_FRAME(op, name, str_rep=None, default_axis='columns',
                            masker=False):

    def na_op(x, y):
        try:
            result = op(x, y)
        except TypeError:
            xrav = x.ravel()
            result = np.empty(x.size, dtype=x.dtype)
            if isinstance(y, (np.ndarray, pd.Series)):
                yrav = y.ravel()
                mask = notnull(xrav) & notnull(yrav)
                result[mask] = op(np.array(list(xrav[mask])),
                                  np.array(list(yrav[mask])))
            else:
                mask = notnull(xrav)
                result[mask] = op(np.array(list(xrav[mask])), y)

            if op == operator.ne:  # pragma: no cover
                np.putmask(result, -mask, True)
            else:
                np.putmask(result, -mask, False)
            result = result.reshape(x.shape)

        return result

    @Appender('Wrapper for flexible comparison methods %s' % name)
    def f(self, other, axis=default_axis, level=None):
        if isinstance(other, pd.DataFrame):    # Another DataFrame
            return self._flex_compare_frame(other, na_op, str_rep, level)

        elif isinstance(other, pd.Series):
            return self._combine_series(other, na_op, None, axis, level)

        elif isinstance(other, (list, tuple)):
            if axis is not None and self._get_axis_name(axis) == 'index':
                casted = pd.Series(other, index=self.index)
            else:
                casted = pd.Series(other, index=self.columns)

            return self._combine_series(casted, na_op, None, axis, level)

        elif isinstance(other, np.ndarray):
            if other.ndim == 1:
                if axis is not None and self._get_axis_name(axis) == 'index':
                    casted = pd.Series(other, index=self.index)
                else:
                    casted = pd.Series(other, index=self.columns)

                return self._combine_series(casted, na_op, None, axis, level)

            elif other.ndim == 2:
                casted = pd.DataFrame(other, index=self.index,
                                      columns=self.columns)

                return self._flex_compare_frame(casted, na_op, str_rep, level)

            else:
                raise ValueError("Incompatible argument shape: %s" %
                                 (other.shape, ))

        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


def _comp_method_FRAME(func, name, str_rep, masker=False):
    @Appender('Wrapper for comparison method %s' % name)
    def f(self, other):
        if isinstance(other, pd.DataFrame):    # Another DataFrame
            return self._compare_frame(other, func, str_rep)
        elif isinstance(other, pd.Series):
            return self._combine_series_infer(other, func)
        else:

            # straight boolean comparisions we want to allow all columns
            # (regardless of dtype to pass thru) See #4537 for discussion.
            res = self._combine_const(other, func, raise_on_error=False)
            return res.fillna(True).astype(bool)

    f.__name__ = name

    return f


frame_flex_funcs = dict(flex_arith_method=_arith_method_FRAME,
                        radd_func=_radd_compat,
                        flex_comp_method=_flex_comp_method_FRAME)


frame_special_funcs = dict(arith_method=_arith_method_FRAME,
                           radd_func=_radd_compat,
                           comp_method=_comp_method_FRAME,
                           bool_method=_arith_method_FRAME)


def _arith_method_PANEL(op, name, str_rep=None, fill_zeros=None,
                        default_axis=None, **eval_kwargs):
    # copied from Series na_op above, but without unnecessary branch for
    # non-scalar
    def na_op(x, y):
        try:
            result = expressions.evaluate(op, str_rep, x, y,
                                          raise_on_error=True, **eval_kwargs)
        except TypeError:

            # TODO: might need to find_common_type here?
            result = pa.empty(len(x), dtype=x.dtype)
            mask = notnull(x)
            result[mask] = op(x[mask], y)
            result, changed = com._maybe_upcast_putmask(result, -mask, pa.NA)

        result = com._fill_zeros(result, y, fill_zeros)
        return result

    # work only for scalars
    def f(self, other):
        if not np.isscalar(other):
            raise ValueError('Simple arithmetic with %s can only be '
                             'done with scalar values' %
                             self._constructor.__name__)

        return self._combine(other, op)
    f.__name__ = name
    return f


def _comp_method_PANEL(op, name, str_rep=None, masker=False):

    def na_op(x, y):
        try:
            result = expressions.evaluate(op, str_rep, x, y,
                                          raise_on_error=True)
        except TypeError:
            xrav = x.ravel()
            result = np.empty(x.size, dtype=bool)
            if isinstance(y, np.ndarray):
                yrav = y.ravel()
                mask = notnull(xrav) & notnull(yrav)
                result[mask] = op(np.array(list(xrav[mask])),
                                  np.array(list(yrav[mask])))
            else:
                mask = notnull(xrav)
                result[mask] = op(np.array(list(xrav[mask])), y)

            if op == operator.ne:  # pragma: no cover
                np.putmask(result, -mask, True)
            else:
                np.putmask(result, -mask, False)
            result = result.reshape(x.shape)

        return result

    @Appender('Wrapper for comparison method %s' % name)
    def f(self, other):
        if isinstance(other, self._constructor):
            return self._compare_constructor(other, na_op)
        elif isinstance(other, (self._constructor_sliced, pd.DataFrame,
                                pd.Series)):
            raise Exception("input needs alignment for this object [%s]" %
                            self._constructor)
        else:
            return self._combine_const(other, na_op)

    f.__name__ = name

    return f


panel_special_funcs = dict(arith_method=_arith_method_PANEL,
                           comp_method=_comp_method_PANEL,
                           bool_method=_arith_method_PANEL)
