import abc
import inspect
from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional, Tuple, Type, Union

import numpy as np

from pandas._config import option_context

from pandas._typing import Axis
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    is_dict_like,
    is_extension_array_dtype,
    is_list_like,
    is_sequence,
)
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.construction import create_series_with_explicit_dtype

if TYPE_CHECKING:
    from pandas import DataFrame, Index, Series

ResType = Dict[int, Any]


def frame_apply(
    obj: "DataFrame",
    func,
    axis: Axis = 0,
    raw: bool = False,
    result_type: Optional[str] = None,
    ignore_failures: bool = False,
    args=None,
    kwds=None,
):
    """ construct and return a row or column based frame apply object """
    axis = obj._get_axis_number(axis)
    klass: Type[FrameApply]
    if axis == 0:
        klass = FrameRowApply
    elif axis == 1:
        klass = FrameColumnApply

    return klass(
        obj,
        func,
        raw=raw,
        result_type=result_type,
        ignore_failures=ignore_failures,
        args=args,
        kwds=kwds,
    )


class FrameApply(metaclass=abc.ABCMeta):

    # ---------------------------------------------------------------
    # Abstract Methods
    axis: int

    @property
    @abc.abstractmethod
    def result_index(self) -> "Index":
        pass

    @property
    @abc.abstractmethod
    def result_columns(self) -> "Index":
        pass

    @property
    @abc.abstractmethod
    def series_generator(self) -> Iterator["Series"]:
        pass

    @abc.abstractmethod
    def wrap_results_for_axis(
        self, results: ResType, res_index: "Index"
    ) -> Union["Series", "DataFrame"]:
        pass

    # ---------------------------------------------------------------

    def __init__(
        self,
        obj: "DataFrame",
        func,
        raw: bool,
        result_type: Optional[str],
        ignore_failures: bool,
        args,
        kwds,
    ):
        self.obj = obj
        self.raw = raw
        self.ignore_failures = ignore_failures
        self.args = args or ()
        self.kwds = kwds or {}

        if result_type not in [None, "reduce", "broadcast", "expand"]:
            raise ValueError(
                "invalid value for result_type, must be one "
                "of {None, 'reduce', 'broadcast', 'expand'}"
            )

        self.result_type = result_type

        # curry if needed
        if (kwds or args) and not isinstance(func, (np.ufunc, str)):

            def f(x):
                return func(x, *args, **kwds)

        else:
            f = func

        self.f = f

    @property
    def res_columns(self) -> "Index":
        return self.result_columns

    @property
    def columns(self) -> "Index":
        return self.obj.columns

    @property
    def index(self) -> "Index":
        return self.obj.index

    @cache_readonly
    def values(self):
        return self.obj.values

    @cache_readonly
    def dtypes(self) -> "Series":
        return self.obj.dtypes

    @property
    def agg_axis(self) -> "Index":
        return self.obj._get_agg_axis(self.axis)

    def get_result(self):
        """ compute the results """
        # dispatch to agg
        if is_list_like(self.f) or is_dict_like(self.f):
            return self.obj.aggregate(self.f, axis=self.axis, *self.args, **self.kwds)

        # all empty
        if len(self.columns) == 0 and len(self.index) == 0:
            return self.apply_empty_result()

        # string dispatch
        if isinstance(self.f, str):
            # Support for `frame.transform('method')`
            # Some methods (shift, etc.) require the axis argument, others
            # don't, so inspect and insert if necessary.
            func = getattr(self.obj, self.f)
            sig = inspect.getfullargspec(func)
            if "axis" in sig.args:
                self.kwds["axis"] = self.axis
            return func(*self.args, **self.kwds)

        # ufunc
        elif isinstance(self.f, np.ufunc):
            with np.errstate(all="ignore"):
                results = self.obj._mgr.apply("apply", func=self.f)
            # _constructor will retain self.index and self.columns
            return self.obj._constructor(data=results)

        # broadcasting
        if self.result_type == "broadcast":
            return self.apply_broadcast(self.obj)

        # one axis empty
        elif not all(self.obj.shape):
            return self.apply_empty_result()

        # raw
        elif self.raw:
            return self.apply_raw()

        return self.apply_standard()

    def apply_empty_result(self):
        """
        we have an empty result; at least 1 axis is 0

        we will try to apply the function to an empty
        series in order to see if this is a reduction function
        """
        # we are not asked to reduce or infer reduction
        # so just return a copy of the existing object
        if self.result_type not in ["reduce", None]:
            return self.obj.copy()

        # we may need to infer
        should_reduce = self.result_type == "reduce"

        from pandas import Series

        if not should_reduce:
            try:
                r = self.f(Series([], dtype=np.float64))
            except Exception:
                pass
            else:
                should_reduce = not isinstance(r, Series)

        if should_reduce:
            if len(self.agg_axis):
                r = self.f(Series([], dtype=np.float64))
            else:
                r = np.nan

            return self.obj._constructor_sliced(r, index=self.agg_axis)
        else:
            return self.obj.copy()

    def apply_raw(self):
        """ apply to the values as a numpy array """

        def wrap_function(func):
            """
            Wrap user supplied function to work around numpy issue.

            see https://github.com/numpy/numpy/issues/8352
            """

            def wrapper(*args, **kwargs):
                result = func(*args, **kwargs)
                if isinstance(result, str):
                    result = np.array(result, dtype=object)
                return result

            return wrapper

        result = np.apply_along_axis(wrap_function(self.f), self.axis, self.values)

        # TODO: mixed type case
        if result.ndim == 2:
            return self.obj._constructor(result, index=self.index, columns=self.columns)
        else:
            return self.obj._constructor_sliced(result, index=self.agg_axis)

    def apply_broadcast(self, target: "DataFrame") -> "DataFrame":
        result_values = np.empty_like(target.values)

        # axis which we want to compare compliance
        result_compare = target.shape[0]

        for i, col in enumerate(target.columns):
            res = self.f(target[col])
            ares = np.asarray(res).ndim

            # must be a scalar or 1d
            if ares > 1:
                raise ValueError("too many dims to broadcast")
            elif ares == 1:

                # must match return dim
                if result_compare != len(res):
                    raise ValueError("cannot broadcast result")

            result_values[:, i] = res

        # we *always* preserve the original index / columns
        result = self.obj._constructor(
            result_values, index=target.index, columns=target.columns
        )
        return result

    def apply_standard(self):
        results, res_index = self.apply_series_generator()

        # wrap results
        return self.wrap_results(results, res_index)

    def apply_series_generator(self) -> Tuple[ResType, "Index"]:
        series_gen = self.series_generator
        res_index = self.result_index

        results = {}

        if self.ignore_failures:
            successes = []
            for i, v in enumerate(series_gen):
                try:
                    results[i] = self.f(v)
                except Exception:
                    pass
                else:
                    successes.append(i)

            # so will work with MultiIndex
            if len(successes) < len(res_index):
                res_index = res_index.take(successes)

        else:
            with option_context("mode.chained_assignment", None):
                for i, v in enumerate(series_gen):
                    # ignore SettingWithCopy here in case the user mutates
                    results[i] = self.f(v)
                    if isinstance(results[i], ABCSeries):
                        # If we have a view on v, we need to make a copy because
                        #  series_generator will swap out the underlying data
                        results[i] = results[i].copy(deep=False)

        return results, res_index

    def wrap_results(
        self, results: ResType, res_index: "Index"
    ) -> Union["Series", "DataFrame"]:
        from pandas import Series

        # see if we can infer the results
        if len(results) > 0 and 0 in results and is_sequence(results[0]):
            return self.wrap_results_for_axis(results, res_index)

        # dict of scalars

        # the default dtype of an empty Series will be `object`, but this
        # code can be hit by df.mean() where the result should have dtype
        # float64 even if it's an empty Series.
        constructor_sliced = self.obj._constructor_sliced
        if constructor_sliced is Series:
            result = create_series_with_explicit_dtype(
                results, dtype_if_empty=np.float64
            )
        else:
            result = constructor_sliced(results)
        result.index = res_index

        return result


class FrameRowApply(FrameApply):
    axis = 0

    def apply_broadcast(self, target: "DataFrame") -> "DataFrame":
        return super().apply_broadcast(target)

    @property
    def series_generator(self):
        return (self.obj._ixs(i, axis=1) for i in range(len(self.columns)))

    @property
    def result_index(self) -> "Index":
        return self.columns

    @property
    def result_columns(self) -> "Index":
        return self.index

    def wrap_results_for_axis(
        self, results: ResType, res_index: "Index"
    ) -> Union["Series", "DataFrame"]:
        """ return the results for the rows """

        if self.result_type == "reduce":
            # e.g. test_apply_dict GH#8735
            res = self.obj._constructor_sliced(results)
            res.index = res_index
            return res

        elif self.result_type is None and all(
            isinstance(x, dict) for x in results.values()
        ):
            # Our operation was a to_dict op e.g.
            #  test_apply_dict GH#8735, test_apply_reduce_rows_to_dict GH#25196
            return self.obj._constructor_sliced(results)

        try:
            result = self.obj._constructor(data=results)
        except ValueError as err:
            if "arrays must all be same length" in str(err):
                # e.g. result = [[2, 3], [1.5], ['foo', 'bar']]
                #  see test_agg_listlike_result GH#29587
                res = self.obj._constructor_sliced(results)
                res.index = res_index
                return res
            else:
                raise

        if not isinstance(results[0], ABCSeries):
            if len(result.index) == len(self.res_columns):
                result.index = self.res_columns

        if len(result.columns) == len(res_index):
            result.columns = res_index

        return result


class FrameColumnApply(FrameApply):
    axis = 1

    def apply_broadcast(self, target: "DataFrame") -> "DataFrame":
        result = super().apply_broadcast(target.T)
        return result.T

    @property
    def series_generator(self):
        values = self.values
        assert len(values) > 0

        # We create one Series object, and will swap out the data inside
        #  of it.  Kids: don't do this at home.
        ser = self.obj._ixs(0, axis=0)
        mgr = ser._mgr
        blk = mgr.blocks[0]

        if is_extension_array_dtype(blk.dtype):
            # values will be incorrect for this block
            # TODO(EA2D): special case would be unnecessary with 2D EAs
            obj = self.obj
            for i in range(len(obj)):
                yield obj._ixs(i, axis=0)

        else:
            for (arr, name) in zip(values, self.index):
                # GH#35462 re-pin mgr in case setitem changed it
                ser._mgr = mgr
                blk.values = arr
                ser.name = name
                yield ser

    @property
    def result_index(self) -> "Index":
        return self.index

    @property
    def result_columns(self) -> "Index":
        return self.columns

    def wrap_results_for_axis(
        self, results: ResType, res_index: "Index"
    ) -> Union["Series", "DataFrame"]:
        """ return the results for the columns """
        result: Union["Series", "DataFrame"]

        # we have requested to expand
        if self.result_type == "expand":
            result = self.infer_to_same_shape(results, res_index)

        # we have a non-series and don't want inference
        elif not isinstance(results[0], ABCSeries):
            result = self.obj._constructor_sliced(results)
            result.index = res_index

        # we may want to infer results
        else:
            result = self.infer_to_same_shape(results, res_index)

        return result

    def infer_to_same_shape(self, results: ResType, res_index: "Index") -> "DataFrame":
        """ infer the results to the same shape as the input object """
        result = self.obj._constructor(data=results)
        result = result.T

        # set the index
        result.index = res_index

        # infer dtypes
        result = result.infer_objects()

        return result
