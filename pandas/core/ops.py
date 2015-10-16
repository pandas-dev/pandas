"""
Arithmetic operations for PandasObjects

This is not a public API.
"""
# necessary to enforce truediv in Python 2.X
from __future__ import division
import operator
import warnings
import numpy as np
import pandas as pd
import datetime
from pandas import compat, lib, tslib
import pandas.index as _index
from pandas.util.decorators import Appender
import pandas.core.common as com
import pandas.computation.expressions as expressions
from pandas.lib import isscalar
from pandas.tslib import iNaT
from pandas.compat import bind_method
from pandas.core.common import(is_list_like, notnull, isnull,
                               _values_from_object, _maybe_match_name,
                               needs_i8_conversion, is_datetimelike_v_numeric,
                               is_integer_dtype, is_categorical_dtype, is_object_dtype,
                               is_timedelta64_dtype, is_datetime64_dtype, is_datetime64tz_dtype,
                               is_bool_dtype)
from pandas.io.common import PerformanceWarning

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
        mod=arith_method(operator.mod, names('mod'), None,
                         default_axis=default_axis, fill_zeros=np.nan),
        pow=arith_method(operator.pow, names('pow'), op('**'),
                         default_axis=default_axis),
        # not entirely sure why this is necessary, but previously was included
        # so it's here to maintain compatibility
        rmul=arith_method(operator.mul, names('rmul'), op('*'),
                          default_axis=default_axis, reversed=True),
        rsub=arith_method(lambda x, y: y - x, names('rsub'), op('-'),
                          default_axis=default_axis, reversed=True),
        rtruediv=arith_method(lambda x, y: operator.truediv(y, x),
                              names('rtruediv'), op('/'), truediv=True,
                              fill_zeros=np.inf, default_axis=default_axis,
                              reversed=True),
        rfloordiv=arith_method(lambda x, y: operator.floordiv(y, x),
                               names('rfloordiv'), op('//'),
                               default_axis=default_axis, fill_zeros=np.inf,
                               reversed=True),
        rpow=arith_method(lambda x, y: y ** x, names('rpow'), op('**'),
                          default_axis=default_axis, reversed=True),
        rmod=arith_method(lambda x, y: y % x, names('rmod'), op('%'),
                          default_axis=default_axis, fill_zeros=np.nan,
                          reversed=True),
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
            xor=bool_method(operator.xor, names('xor'), op('^')),
            rand_=bool_method(lambda x, y: operator.and_(y, x),
                              names('rand_'), op('&')),
            ror_=bool_method(lambda x, y: operator.or_(y, x), names('ror_'), op('|')),
            rxor=bool_method(lambda x, y: operator.xor(y, x), names('rxor'), op('^'))
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

    def _wrap_inplace_method(method):
        """
        return an inplace wrapper for this method
        """

        def f(self, other):
            result = method(self, other)

            # this makes sure that we are aligned like the input
            # we are updating inplace so we want to ignore is_copy
            self._update_inplace(result.reindex_like(self,copy=False)._data,
                                 verify_is_copy=False)

            return self
        return f

    new_methods.update(dict(
        __iadd__=_wrap_inplace_method(new_methods["__add__"]),
        __isub__=_wrap_inplace_method(new_methods["__sub__"]),
        __imul__=_wrap_inplace_method(new_methods["__mul__"]),
        __itruediv__=_wrap_inplace_method(new_methods["__truediv__"]),
        __ipow__=_wrap_inplace_method(new_methods["__pow__"]),
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
    flex_arith_method : function
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
    fill_value = iNaT
    wrap_results = staticmethod(lambda x: x)
    dtype = None

    def __init__(self, left, right, name, na_op):

        # need to make sure that we are aligning the data
        if isinstance(left, pd.Series) and isinstance(right, pd.Series):
            left, right = left.align(right,copy=False)

        lvalues = self._convert_to_array(left, name=name)
        rvalues = self._convert_to_array(right, name=name, other=lvalues)

        self.name = name
        self.na_op = na_op

        # left
        self.left = left
        self.is_offset_lhs = self._is_offset(left)
        self.is_timedelta_lhs = is_timedelta64_dtype(lvalues)
        self.is_datetime64_lhs = is_datetime64_dtype(lvalues)
        self.is_datetime64tz_lhs = is_datetime64tz_dtype(lvalues)
        self.is_datetime_lhs = self.is_datetime64_lhs or self.is_datetime64tz_lhs
        self.is_integer_lhs = left.dtype.kind in ['i', 'u']

        # right
        self.right = right
        self.is_offset_rhs = self._is_offset(right)
        self.is_datetime64_rhs = is_datetime64_dtype(rvalues)
        self.is_datetime64tz_rhs = is_datetime64tz_dtype(rvalues)
        self.is_datetime_rhs = self.is_datetime64_rhs or self.is_datetime64tz_rhs
        self.is_timedelta_rhs = is_timedelta64_dtype(rvalues)
        self.is_integer_rhs = rvalues.dtype.kind in ('i', 'u')

        self._validate(lvalues, rvalues, name)
        self.lvalues, self.rvalues = self._convert_for_datetime(lvalues, rvalues)

    def _validate(self, lvalues, rvalues, name):
        # timedelta and integer mul/div

        if (self.is_timedelta_lhs and self.is_integer_rhs) or (
            self.is_integer_lhs and self.is_timedelta_rhs):

            if name not in ('__div__', '__truediv__', '__mul__'):
                raise TypeError("can only operate on a timedelta and an "
                                "integer for division, but the operator [%s]"
                                "was passed" % name)

        # 2 datetimes
        elif self.is_datetime_lhs and self.is_datetime_rhs:

            if name not in ('__sub__','__rsub__'):
                raise TypeError("can only operate on a datetimes for"
                                " subtraction, but the operator [%s] was"
                                " passed" % name)

            # if tz's must be equal (same or None)
            if getattr(lvalues,'tz',None) != getattr(rvalues,'tz',None):
                raise ValueError("Incompatbile tz's on datetime subtraction ops")

        # 2 timedeltas
        elif ((self.is_timedelta_lhs and
               (self.is_timedelta_rhs or self.is_offset_rhs)) or
              (self.is_timedelta_rhs and
               (self.is_timedelta_lhs or self.is_offset_lhs))):

            if name not in ('__div__', '__rdiv__', '__truediv__', '__rtruediv__',
                            '__add__', '__radd__', '__sub__', '__rsub__'):
                raise TypeError("can only operate on a timedeltas for "
                                "addition, subtraction, and division, but the"
                                " operator [%s] was passed" % name)

        # datetime and timedelta/DateOffset
        elif (self.is_datetime_lhs and
              (self.is_timedelta_rhs or self.is_offset_rhs)):

            if name not in ('__add__', '__radd__', '__sub__'):
                raise TypeError("can only operate on a datetime with a rhs of"
                                " a timedelta/DateOffset for addition and subtraction,"
                                " but the operator [%s] was passed" %
                                name)

        elif ((self.is_timedelta_lhs or self.is_offset_lhs)
              and self.is_datetime_rhs):

            if name not in ('__add__', '__radd__'):
                raise TypeError("can only operate on a timedelta/DateOffset and"
                                " a datetime for addition, but the operator"
                                " [%s] was passed" % name)
        else:
            raise TypeError('cannot operate on a series with out a rhs '
                            'of a series/ndarray of type datetime64[ns] '
                            'or a timedelta')

    def _convert_to_array(self, values, name=None, other=None):
        """converts values to ndarray"""
        from pandas.tseries.timedeltas import to_timedelta

        ovalues = values
        if not is_list_like(values):
            values = np.array([values])

        inferred_type = lib.infer_dtype(values)

        if inferred_type in ('datetime64', 'datetime', 'date', 'time'):
            # if we have a other of timedelta, but use pd.NaT here we
            # we are in the wrong path
            if (other is not None and other.dtype == 'timedelta64[ns]' and
                    all(isnull(v) for v in values)):
                values = np.empty(values.shape, dtype=other.dtype)
                values[:] = iNaT

            # a datelike
            elif isinstance(values, pd.DatetimeIndex):
                values = values.to_series()
            # datetime with tz
            elif isinstance(ovalues, datetime.datetime) and hasattr(ovalues,'tz'):
                values = pd.DatetimeIndex(values)
            # datetime array with tz
            elif com.is_datetimetz(values):
                if isinstance(values, pd.Series):
                    values = values._values
            elif not (isinstance(values, (np.ndarray, pd.Series)) and
                      is_datetime64_dtype(values)):
                values = tslib.array_to_datetime(values)
        elif inferred_type in ('timedelta', 'timedelta64'):
            # have a timedelta, convert to to ns here
            values = to_timedelta(values, errors='coerce')
        elif inferred_type == 'integer':
            # py3 compat where dtype is 'm' but is an integer
            if values.dtype.kind == 'm':
                values = values.astype('timedelta64[ns]')
            elif isinstance(values, pd.PeriodIndex):
                values = values.to_timestamp().to_series()
            elif name not in ('__truediv__', '__div__', '__mul__'):
                raise TypeError("incompatible type for a datetime/timedelta "
                                "operation [{0}]".format(name))
        elif inferred_type == 'floating':
            # all nan, so ok, use the other dtype (e.g. timedelta or datetime)
            if isnull(values).all():
                values = np.empty(values.shape, dtype=other.dtype)
                values[:] = iNaT
            else:
                raise TypeError(
                    'incompatible type [{0}] for a datetime/timedelta '
                    'operation'.format(np.array(values).dtype))
        elif self._is_offset(values):
            return values
        else:
            raise TypeError("incompatible type [{0}] for a datetime/timedelta"
                            " operation".format(np.array(values).dtype))

        return values

    def _convert_for_datetime(self, lvalues, rvalues):
        from pandas.tseries.timedeltas import to_timedelta

        mask = isnull(lvalues) | isnull(rvalues)

        # datetimes require views
        if self.is_datetime_lhs or self.is_datetime_rhs:

            # datetime subtraction means timedelta
            if self.is_datetime_lhs and self.is_datetime_rhs:
                self.dtype = 'timedelta64[ns]'
            elif self.is_datetime64tz_lhs:
                self.dtype = lvalues.dtype
            elif self.is_datetime64tz_rhs:
                self.dtype = rvalues.dtype
            else:
                self.dtype = 'datetime64[ns]'

            # if adding single offset try vectorized path
            # in DatetimeIndex; otherwise elementwise apply
            def _offset(lvalues, rvalues):
                if len(lvalues) == 1:
                    rvalues = pd.DatetimeIndex(rvalues)
                    lvalues = lvalues[0]
                else:
                    warnings.warn("Adding/subtracting array of DateOffsets to Series not vectorized",
                                  PerformanceWarning)
                    rvalues = rvalues.astype('O')

                # pass thru on the na_op
                self.na_op = lambda x, y: getattr(x,self.name)(y)
                return lvalues, rvalues


            if self.is_offset_lhs:
                lvalues, rvalues = _offset(lvalues, rvalues)
            elif self.is_offset_rhs:
                rvalues, lvalues = _offset(rvalues, lvalues)
            else:

                # with tz, convert to UTC
                if self.is_datetime64tz_lhs:
                    lvalues = lvalues.tz_localize(None)
                if self.is_datetime64tz_rhs:
                    rvalues = rvalues.tz_localize(None)

                lvalues = lvalues.view(np.int64)
                rvalues = rvalues.view(np.int64)

        # otherwise it's a timedelta
        else:

            self.dtype = 'timedelta64[ns]'

            # convert Tick DateOffset to underlying delta
            if self.is_offset_lhs:
                lvalues = to_timedelta(lvalues)
            if self.is_offset_rhs:
                rvalues = to_timedelta(rvalues)

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
        if mask.any():
            def f(x):

                # datetime64[ns]/timedelta64[ns] masking
                try:
                    x = np.array(x, dtype=self.dtype)
                except TypeError:
                    x = np.array(x, dtype='datetime64[ns]')

                np.putmask(x, mask, self.fill_value)
                return x
            self.wrap_results = f

        return lvalues, rvalues


    def _is_offset(self, arr_or_obj):
        """ check if obj or all elements of list-like is DateOffset """
        if isinstance(arr_or_obj, pd.DateOffset):
            return True
        elif is_list_like(arr_or_obj):
            return all(isinstance(x, pd.DateOffset) for x in arr_or_obj)
        else:
            return False

    @classmethod
    def maybe_convert_for_time_op(cls, left, right, name, na_op):
        """
        if ``left`` and ``right`` are appropriate for datetime arithmetic with
        operation ``name``, processes them and returns a ``_TimeOp`` object
        that stores all the required values.  Otherwise, it will generate
        either a ``NotImplementedError`` or ``None``, indicating that the
        operation is unsupported for datetimes (e.g., an unsupported r_op) or
        that the data is not the right type for time ops.
        """
        # decide if we can do it
        is_timedelta_lhs = is_timedelta64_dtype(left)
        is_datetime_lhs = is_datetime64_dtype(left) or is_datetime64tz_dtype(left)

        if not (is_datetime_lhs or is_timedelta_lhs):
            return None

        return cls(left, right, name, na_op)


def _arith_method_SERIES(op, name, str_rep, fill_zeros=None,
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
            if isinstance(y, (np.ndarray, pd.Series, pd.Index)):
                dtype = np.find_common_type([x.dtype, y.dtype], [])
                result = np.empty(x.size, dtype=dtype)
                mask = notnull(x) & notnull(y)
                result[mask] = op(x[mask], _values_from_object(y[mask]))
            elif isinstance(x, np.ndarray):
                result = np.empty(len(x), dtype=x.dtype)
                mask = notnull(x)
                result[mask] = op(x[mask], y)
            else:
                raise TypeError("{typ} cannot perform the operation {op}".format(typ=type(x).__name__,op=str_rep))

            result, changed = com._maybe_upcast_putmask(result, ~mask, np.nan)

        result = com._fill_zeros(result, x, y, name, fill_zeros)
        return result

    def wrapper(left, right, name=name, na_op=na_op):

        if isinstance(right, pd.DataFrame):
            return NotImplemented

        time_converted = _TimeOp.maybe_convert_for_time_op(left, right, name, na_op)

        if time_converted is None:
            lvalues, rvalues = left, right
            dtype = None
            wrap_results = lambda x: x
        elif time_converted == NotImplemented:
            return NotImplemented
        else:
            left, right = time_converted.left, time_converted.right
            lvalues, rvalues = time_converted.lvalues, time_converted.rvalues
            dtype = time_converted.dtype
            wrap_results = time_converted.wrap_results
            na_op = time_converted.na_op

        if isinstance(rvalues, pd.Series):
            rindex = getattr(rvalues,'index',rvalues)
            name = _maybe_match_name(left, rvalues)
            lvalues = getattr(lvalues, 'values', lvalues)
            rvalues = getattr(rvalues, 'values', rvalues)
            if left.index.equals(rindex):
                index = left.index
            else:
                index, lidx, ridx = left.index.join(rindex, how='outer',
                                                       return_indexers=True)

                if lidx is not None:
                    lvalues = com.take_1d(lvalues, lidx)

                if ridx is not None:
                    rvalues = com.take_1d(rvalues, ridx)

            arr = na_op(lvalues, rvalues)

            return left._constructor(wrap_results(arr), index=index,
                                     name=name, dtype=dtype)
        else:
            # scalars
            if hasattr(lvalues, 'values') and not isinstance(lvalues, pd.DatetimeIndex):
                lvalues = lvalues.values

            return left._constructor(wrap_results(na_op(lvalues, rvalues)),
                                     index=left.index, name=left.name,
                                     dtype=dtype)
    return wrapper


def _comp_method_SERIES(op, name, str_rep, masker=False):
    """
    Wrapper function for Series arithmetic operations, to avoid
    code duplication.
    """
    def na_op(x, y):

        # dispatch to the categorical if we have a categorical
        # in either operand
        if is_categorical_dtype(x):
            return op(x,y)
        elif is_categorical_dtype(y) and not isscalar(y):
            return op(y,x)

        if is_object_dtype(x.dtype):
            if isinstance(y, list):
                y = lib.list_to_object_array(y)

            if isinstance(y, (np.ndarray, pd.Series)):
                if not is_object_dtype(y.dtype):
                    result = lib.vec_compare(x, y.astype(np.object_), op)
                else:
                    result = lib.vec_compare(x, y, op)
            else:
                result = lib.scalar_compare(x, y, op)
        else:

            # we want to compare like types
            # we only want to convert to integer like if
            # we are not NotImplemented, otherwise
            # we would allow datetime64 (but viewed as i8) against
            # integer comparisons
            if is_datetimelike_v_numeric(x, y):
                raise TypeError("invalid type comparison")

            # numpy does not like comparisons vs None
            if isscalar(y) and isnull(y):
                if name == '__ne__':
                    return np.ones(len(x), dtype=bool)
                else:
                    return np.zeros(len(x), dtype=bool)

            # we have a datetime/timedelta and may need to convert
            mask = None
            if needs_i8_conversion(x) or (not isscalar(y) and needs_i8_conversion(y)):

                if isscalar(y):
                    y = _index.convert_scalar(x,_values_from_object(y))
                else:
                    y = y.view('i8')

                mask = isnull(x)

                x = x.view('i8')

            try:
                result = getattr(x, name)(y)
                if result is NotImplemented:
                    raise TypeError("invalid type comparison")
            except AttributeError:
                result = op(x, y)

            if mask is not None and mask.any():
                result[mask] = masker

        return result

    def wrapper(self, other, axis=None):
        # Validate the axis parameter
        if axis is not None:
            self._get_axis_number(axis)

        if isinstance(other, pd.Series):
            name = _maybe_match_name(self, other)
            if len(self) != len(other):
                raise ValueError('Series lengths must match to compare')
            return self._constructor(na_op(self.values, other.values),
                                     index=self.index, name=name)
        elif isinstance(other, pd.DataFrame):  # pragma: no cover
            return NotImplemented
        elif isinstance(other, (np.ndarray, pd.Index)):
            if len(self) != len(other):
                raise ValueError('Lengths must match to compare')
            return self._constructor(na_op(self.values, np.asarray(other)),
                                     index=self.index).__finalize__(self)
        elif isinstance(other, pd.Categorical):
            if not is_categorical_dtype(self):
                msg = "Cannot compare a Categorical for op {op} with Series of dtype {typ}.\n"\
                      "If you want to compare values, use 'series <op> np.asarray(other)'."
                raise TypeError(msg.format(op=op,typ=self.dtype))


        if is_categorical_dtype(self):
            # cats are a special case as get_values() would return an ndarray, which would then
            # not take categories ordering into account
            # we can go directly to op, as the na_op would just test again and dispatch to it.
            res = op(self.values, other)
        else:
            values = self.get_values()
            if isinstance(other, (list, np.ndarray)):
                other = np.asarray(other)

            res = na_op(values, other)
            if isscalar(res):
                raise TypeError('Could not compare %s type with Series'
                                % type(other))

            # always return a full value series here
            res = _values_from_object(res)

        res = pd.Series(res, index=self.index, name=self.name,
                        dtype='bool')
        return res
    return wrapper


def _bool_method_SERIES(op, name, str_rep):
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

            if isinstance(y, (np.ndarray, pd.Series)):
                if (is_bool_dtype(x.dtype) and is_bool_dtype(y.dtype)):
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
        is_self_int_dtype = is_integer_dtype(self.dtype)

        fill_int = lambda x: x.fillna(0)
        fill_bool = lambda x: x.fillna(False).astype(bool)

        if isinstance(other, pd.Series):
            name = _maybe_match_name(self, other)
            other = other.reindex_like(self)
            is_other_int_dtype = is_integer_dtype(other.dtype)
            other = fill_int(other) if is_other_int_dtype else fill_bool(other)

            filler = fill_int if is_self_int_dtype and is_other_int_dtype else fill_bool
            return filler(self._constructor(na_op(self.values, other.values),
                                     index=self.index,
                                     name=name))

        elif isinstance(other, pd.DataFrame):
            return NotImplemented

        else:
            # scalars, list, tuple, np.array
            filler = fill_int if is_self_int_dtype and is_integer_dtype(np.asarray(other)) else fill_bool
            return filler(self._constructor(na_op(self.values, other),
                                    index=self.index)).__finalize__(self)

    return wrapper


def _radd_compat(left, right):
    radd = lambda x, y: y + x
    # GH #353, NumPy 1.5.1 workaround
    try:
        output = radd(left, right)
    except TypeError:
        raise

    return output

_op_descriptions = {'add': {'op': '+', 'desc': 'Addition', 'reversed': False, 'reverse': 'radd'},
                    'sub': {'op': '-', 'desc': 'Subtraction', 'reversed': False, 'reverse': 'rsub'},
                    'mul': {'op': '*', 'desc': 'Multiplication', 'reversed': False, 'reverse': 'rmul'},
                    'mod': {'op': '%', 'desc': 'Modulo', 'reversed': False, 'reverse': 'rmod'},
                    'pow': {'op': '**', 'desc': 'Exponential power', 'reversed': False, 'reverse': 'rpow'},
                    'truediv': {'op': '/', 'desc': 'Floating division', 'reversed': False, 'reverse': 'rtruediv'},
                    'floordiv': {'op': '//', 'desc': 'Integer division', 'reversed': False, 'reverse': 'rfloordiv'}}

_op_names = list(_op_descriptions.keys())
for k in _op_names:
    reverse_op = _op_descriptions[k]['reverse']
    _op_descriptions[reverse_op] = _op_descriptions[k].copy()
    _op_descriptions[reverse_op]['reversed'] = True
    _op_descriptions[reverse_op]['reverse'] = k

def _flex_method_SERIES(op, name, str_rep, default_axis=None,
                        fill_zeros=None, **eval_kwargs):
    op_name = name.replace('__', '')
    op_desc = _op_descriptions[op_name]
    if op_desc['reversed']:
        equiv = 'other ' + op_desc['op'] + ' series'
    else:
        equiv = 'series ' + op_desc['op'] + ' other'

    doc = """
    %s of series and other, element-wise (binary operator `%s`).

    Equivalent to ``%s``, but with support to substitute a fill_value for
    missing data in one of the inputs.

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

    See also
    --------
    Series.%s
    """ % (op_desc['desc'], op_name, equiv, op_desc['reverse'])

    @Appender(doc)
    def flex_wrapper(self, other, level=None, fill_value=None, axis=0):
        # validate axis
        self._get_axis_number(axis)
        if isinstance(other, pd.Series):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, pd.Series, list, tuple)):
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
            elif hasattr(x,'size'):
                result = np.empty(x.size, dtype=x.dtype)
                mask = notnull(xrav)
                xrav = xrav[mask]
                if np.prod(xrav.shape):
                    result[mask] = op(xrav, y)
            else:
                raise TypeError("cannot perform operation {op} between objects "
                                "of type {x} and {y}".format(op=name,x=type(x),y=type(y)))

            result, changed = com._maybe_upcast_putmask(result, ~mask, np.nan)
            result = result.reshape(x.shape)

        result = com._fill_zeros(result, x, y, name, fill_zeros)

        return result

    if name in _op_descriptions:
        op_name = name.replace('__', '')
        op_desc = _op_descriptions[op_name]
        if op_desc['reversed']:
            equiv = 'other ' + op_desc['op'] + ' dataframe'
        else:
            equiv = 'dataframe ' + op_desc['op'] + ' other'

        doc = """
        %s of dataframe and other, element-wise (binary operator `%s`).

        Equivalent to ``%s``, but with support to substitute a fill_value for
        missing data in one of the inputs.

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

        See also
        --------
        DataFrame.%s
        """ % (op_desc['desc'], op_name, equiv, op_desc['reverse'])
    else:
        doc = _arith_doc_FRAME % name

    @Appender(doc)
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
        elif isinstance(other, np.ndarray) and other.ndim:  # skips np scalar
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
                np.putmask(result, ~mask, True)
            else:
                np.putmask(result, ~mask, False)
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
            result = np.empty(len(x), dtype=x.dtype)
            mask = notnull(x)
            result[mask] = op(x[mask], y)
            result, changed = com._maybe_upcast_putmask(result, ~mask, np.nan)

        result = com._fill_zeros(result, x, y, name, fill_zeros)
        return result

    # work only for scalars
    def f(self, other):
        if not isscalar(other):
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
                np.putmask(result, ~mask, True)
            else:
                np.putmask(result, ~mask, False)
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
