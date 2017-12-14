import numpy as np
from pandas import compat
from pandas._libs import lib
from pandas.core.dtypes.common import (
    is_extension_type,
    is_sequence)

from pandas.io.formats.printing import pprint_thing


def frame_apply(obj, func, axis=0, broadcast=False,
                raw=False, reduce=None, args=(), **kwds):
    """ construct and return a row or column based frame apply object """

    axis = obj._get_axis_number(axis)
    if axis == 0:
        klass = FrameRowApply
    elif axis == 1:
        klass = FrameColumnApply

    return klass(obj, func, broadcast=broadcast,
                 raw=raw, reduce=reduce, args=args, kwds=kwds)


class FrameApply(object):

    def __init__(self, obj, func, broadcast, raw, reduce, args, kwds):
        self.obj = obj
        self.broadcast = broadcast
        self.raw = raw
        self.reduce = reduce
        self.args = args

        self.ignore_failures = kwds.pop('ignore_failures', False)
        self.kwds = kwds

        # curry if needed
        if kwds or args and not isinstance(func, np.ufunc):
            def f(x):
                return func(x, *args, **kwds)
        else:
            f = func

        self.f = f

    @property
    def columns(self):
        return self.obj.columns

    @property
    def index(self):
        return self.obj.index

    @property
    def values(self):
        return self.obj.values

    @property
    def agg_axis(self):
        return self.obj._get_agg_axis(self.axis)

    def get_result(self):
        """ compute the results """

        # all empty
        if len(self.columns) == 0 and len(self.index) == 0:
            return self.apply_empty_result()

        # string dispatch
        if isinstance(self.f, compat.string_types):
            if self.axis:
                self.kwds['axis'] = self.axis
            return getattr(self.obj, self.f)(*self.args, **self.kwds)

        # ufunc
        elif isinstance(self.f, np.ufunc):
            with np.errstate(all='ignore'):
                results = self.f(self.values)
            return self.obj._constructor(data=results, index=self.index,
                                         columns=self.columns, copy=False)

        # broadcasting
        if self.broadcast:
            return self.apply_broadcast()

        # one axis empty
        if not all(self.obj.shape):
            return self.apply_empty_result()

        # raw
        if self.raw and not self.obj._is_mixed_type:
            return self.apply_raw()

        return self.apply_standard()

    def apply_empty_result(self):
        from pandas import Series
        reduce = self.reduce

        if reduce is None:
            reduce = False

            EMPTY_SERIES = Series([])
            try:
                r = self.f(EMPTY_SERIES, *self.args, **self.kwds)
                reduce = not isinstance(r, Series)
            except Exception:
                pass

        if reduce:
            return Series(np.nan, index=self.agg_axis)
        else:
            return self.obj.copy()

    def apply_raw(self):
        try:
            result = lib.reduce(self.values, self.f, axis=self.axis)
        except Exception:
            result = np.apply_along_axis(self.f, self.axis, self.values)

        # TODO: mixed type case
        from pandas import DataFrame, Series
        if result.ndim == 2:
            return DataFrame(result, index=self.index, columns=self.columns)
        else:
            return Series(result, index=self.agg_axis)

    def apply_standard(self):
        from pandas import Series

        reduce = self.reduce
        if reduce is None:
            reduce = True

        # try to reduce first (by default)
        # this only matters if the reduction in values is of different dtype
        # e.g. if we want to apply to a SparseFrame, then can't directly reduce
        if reduce:
            values = self.values

            # we cannot reduce using non-numpy dtypes,
            # as demonstrated in gh-12244
            if not is_extension_type(values):

                # Create a dummy Series from an empty array
                index = self.obj._get_axis(self.axis)
                empty_arr = np.empty(len(index), dtype=values.dtype)

                dummy = Series(empty_arr, index=index, dtype=values.dtype)

                try:
                    labels = self.agg_axis
                    result = lib.reduce(values, self.f,
                                        axis=self.axis,
                                        dummy=dummy,
                                        labels=labels)
                    return Series(result, index=labels)
                except Exception:
                    pass

        # compute the result using the series generator
        results, res_index, res_columns = self._apply_series_generator()

        # wrap results
        return self.wrap_results(results, res_index, res_columns)

    def _apply_series_generator(self):
        series_gen = self.series_generator
        res_index = self.result_index
        res_columns = self.result_columns

        i = None
        keys = []
        results = {}
        if self.ignore_failures:
            successes = []
            for i, v in enumerate(series_gen):
                try:
                    results[i] = self.f(v)
                    keys.append(v.name)
                    successes.append(i)
                except Exception:
                    pass

            # so will work with MultiIndex
            if len(successes) < len(res_index):
                res_index = res_index.take(successes)

        else:
            try:
                for i, v in enumerate(series_gen):
                    results[i] = self.f(v)
                    keys.append(v.name)
            except Exception as e:
                if hasattr(e, 'args'):

                    # make sure i is defined
                    if i is not None:
                        k = res_index[i]
                        e.args = e.args + ('occurred at index %s' %
                                           pprint_thing(k), )
                raise

        return results, res_index, res_columns

    def wrap_results(self, results, res_index, res_columns):
        from pandas import Series

        if len(results) > 0 and is_sequence(results[0]):
            if not isinstance(results[0], Series):
                index = res_columns
            else:
                index = None

            result = self.obj._constructor(data=results, index=index)
            result.columns = res_index

            if self.axis == 1:
                result = result.T
            result = result._convert(
                datetime=True, timedelta=True, copy=False)

        else:

            result = Series(results)
            result.index = res_index

        return result

    def _apply_broadcast(self, target):
        result_values = np.empty_like(target.values)
        columns = target.columns
        for i, col in enumerate(columns):
            result_values[:, i] = self.f(target[col])

        result = self.obj._constructor(result_values, index=target.index,
                                       columns=target.columns)
        return result


class FrameRowApply(FrameApply):
    axis = 0

    def get_result(self):

        # dispatch to agg
        if isinstance(self.f, (list, dict)):
            return self.obj.aggregate(self.f, axis=self.axis,
                                      *self.args, **self.kwds)

        return super(FrameRowApply, self).get_result()

    def apply_broadcast(self):
        return self._apply_broadcast(self.obj)

    @property
    def series_generator(self):
        return (self.obj._ixs(i, axis=1)
                for i in range(len(self.columns)))

    @property
    def result_index(self):
        return self.columns

    @property
    def result_columns(self):
        return self.index


class FrameColumnApply(FrameApply):
    axis = 1

    def __init__(self, obj, func, broadcast, raw, reduce, args, kwds):
        super(FrameColumnApply, self).__init__(obj, func, broadcast,
                                               raw, reduce, args, kwds)

        # skip if we are mixed datelike and trying reduce across axes
        # GH6125
        if self.reduce:
            if self.obj._is_mixed_type and self.obj._is_datelike_mixed_type:
                self.reduce = False

    def apply_broadcast(self):
        return self._apply_broadcast(self.obj.T).T

    @property
    def series_generator(self):
        from pandas import Series
        dtype = object if self.obj._is_mixed_type else None
        return (Series._from_array(arr, index=self.columns, name=name,
                                   dtype=dtype)
                for i, (arr, name) in enumerate(zip(self.values,
                                                    self.index)))

    @property
    def result_index(self):
        return self.index

    @property
    def result_columns(self):
        return self.columns
