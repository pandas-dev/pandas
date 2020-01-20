import abc
import inspect
from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional, Tuple, Type, Union

import numpy as np

from pandas._libs import reduction as libreduction
from pandas._typing import Axis
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    is_dict_like,
    is_extension_array_dtype,
    is_list_like,
    is_sequence,
)
from pandas.core.dtypes.generic import ABCMultiIndex, ABCSeries

from pandas.core.construction import create_series_with_explicit_dtype

if TYPE_CHECKING:
    from pandas import DataFrame, Series, Index

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
                results = self.obj._data.apply("apply", func=self.f)
            return self.obj._constructor(
                data=results, index=self.index, columns=self.columns, copy=False
            )

        # broadcasting
        if self.result_type == "broadcast":
            return self.apply_broadcast(self.obj)

        # one axis empty
        elif not all(self.obj.shape):
            return self.apply_empty_result()

        # raw
        elif self.raw and not self.obj._is_mixed_type:
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
        try:
            result = libreduction.compute_reduction(self.values, self.f, axis=self.axis)
        except ValueError as err:
            if "Function does not reduce" not in str(err):
                # catch only ValueError raised intentionally in libreduction
                raise
            # We expect np.apply_along_axis to give a two-dimensional result, or
            #  also raise.
            result = np.apply_along_axis(self.f, self.axis, self.values)

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

        # try to reduce first (by default)
        # this only matters if the reduction in values is of different dtype
        # e.g. if we want to apply to a SparseFrame, then can't directly reduce

        # we cannot reduce using non-numpy dtypes,
        # as demonstrated in gh-12244
        if (
            self.result_type in ["reduce", None]
            and not self.dtypes.apply(is_extension_array_dtype).any()
            # Disallow complex_internals since libreduction shortcut
            #  cannot handle MultiIndex
            and not isinstance(self.agg_axis, ABCMultiIndex)
        ):

            values = self.values
            index = self.obj._get_axis(self.axis)
            labels = self.agg_axis
            empty_arr = np.empty(len(index), dtype=values.dtype)

            # Preserve subclass for e.g. test_subclassed_apply
            dummy = self.obj._constructor_sliced(
                empty_arr, index=index, dtype=values.dtype
            )

            try:
                result = libreduction.compute_reduction(
                    values, self.f, axis=self.axis, dummy=dummy, labels=labels
                )
            except ValueError as err:
                if "Function does not reduce" not in str(err):
                    # catch only ValueError raised intentionally in libreduction
                    raise
            except TypeError:
                # e.g. test_apply_ignore_failures we just ignore
                if not self.ignore_failures:
                    raise
            except ZeroDivisionError:
                # reached via numexpr; fall back to python implementation
                pass
            else:
                return self.obj._constructor_sliced(result, index=labels)

        # compute the result using the series generator
        results, res_index = self.apply_series_generator()

        # wrap results
        return self.wrap_results(results, res_index)

    def apply_series_generator(self) -> Tuple[ResType, "Index"]:
        series_gen = self.series_generator
        res_index = self.result_index

        keys = []
        results = {}
        if self.ignore_failures:
            successes = []
            for i, v in enumerate(series_gen):
                try:
                    results[i] = self.f(v)
                except Exception:
                    pass
                else:
                    keys.append(v.name)
                    successes.append(i)

            # so will work with MultiIndex
            if len(successes) < len(res_index):
                res_index = res_index.take(successes)

        else:
            for i, v in enumerate(series_gen):
                results[i] = self.f(v)
                keys.append(v.name)

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
    ) -> "DataFrame":
        """ return the results for the rows """

        result = self.obj._constructor(data=results)

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
        constructor = self.obj._constructor_sliced
        return (
            constructor(arr, index=self.columns, name=name)
            for i, (arr, name) in enumerate(zip(self.values, self.index))
        )

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
            from pandas import Series

            result = Series(results)
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
