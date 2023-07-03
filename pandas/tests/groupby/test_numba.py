import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Index,
    Series,
)
import pandas._testing as tm

pytestmark = pytest.mark.single_cpu


@td.skip_if_no("numba")
@pytest.mark.filterwarnings("ignore")
# Filter warnings when parallel=True and the function can't be parallelized by Numba
class TestEngine:
    def test_cython_vs_numba_frame(
        self, sort, nogil, parallel, nopython, numba_supported_reductions
    ):
        func, kwargs = numba_supported_reductions
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        gb = df.groupby("a", sort=sort)
        result = getattr(gb, func)(
            engine="numba", engine_kwargs=engine_kwargs, **kwargs
        )
        expected = getattr(gb, func)(**kwargs)
        tm.assert_frame_equal(result, expected)

    def test_cython_vs_numba_getitem(
        self, sort, nogil, parallel, nopython, numba_supported_reductions
    ):
        func, kwargs = numba_supported_reductions
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        gb = df.groupby("a", sort=sort)["c"]
        result = getattr(gb, func)(
            engine="numba", engine_kwargs=engine_kwargs, **kwargs
        )
        expected = getattr(gb, func)(**kwargs)
        tm.assert_series_equal(result, expected)

    def test_cython_vs_numba_series(
        self, sort, nogil, parallel, nopython, numba_supported_reductions
    ):
        func, kwargs = numba_supported_reductions
        ser = Series(range(3), index=[1, 2, 1], name="foo")
        engine_kwargs = {"nogil": nogil, "parallel": parallel, "nopython": nopython}
        gb = ser.groupby(level=0, sort=sort)
        result = getattr(gb, func)(
            engine="numba", engine_kwargs=engine_kwargs, **kwargs
        )
        expected = getattr(gb, func)(**kwargs)
        tm.assert_series_equal(result, expected)

    def test_as_index_false_unsupported(self, numba_supported_reductions):
        func, kwargs = numba_supported_reductions
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", as_index=False)
        with pytest.raises(NotImplementedError, match="as_index=False"):
            getattr(gb, func)(engine="numba", **kwargs)

    def test_axis_1_unsupported(self, numba_supported_reductions):
        func, kwargs = numba_supported_reductions
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", axis=1)
        with pytest.raises(NotImplementedError, match="axis=1"):
            getattr(gb, func)(engine="numba", **kwargs)

    @pytest.mark.parametrize("group_keys", [True, False])
    @pytest.mark.parametrize("as_index", [True, False])
    def test_cython_vs_numba_apply_identity(self, group_keys, as_index):
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", group_keys=group_keys, as_index=as_index)
        result = gb.apply(lambda group: group, engine="numba")
        expected = gb.apply(lambda group: group, engine="python")
        # TODO: Fix column anmes by exposing them to numba
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("group_keys", [True, False])
    @pytest.mark.parametrize("as_index", [True, False])
    def test_cython_vs_numba_apply_transform(self, group_keys, as_index):
        def f(group):
            new_values = (group.values - group.values.mean()) / group.values.std()
            return DataFrame(new_values, index=group.index)

        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", group_keys=group_keys, as_index=as_index)
        result = gb.apply(f, engine="numba")
        expected = gb.apply(f, engine="python")

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("group_keys", [True, False])
    @pytest.mark.parametrize("as_index", [True, False])
    def test_cython_vs_numba_apply_reduce_scalar(self, group_keys, as_index):
        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", group_keys=group_keys, as_index=as_index)
        result = gb.apply(lambda group: group.values.sum(), engine="numba")
        expected = gb.apply(lambda group: group.values.sum(axis=None), engine="python")
        if isinstance(expected, DataFrame):
            tm.assert_frame_equal(result, expected)
        else:
            tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("group_keys", [True, False])
    @pytest.mark.parametrize("as_index", [True, False])
    def test_cython_vs_numba_apply_reduce_axis_0(self, group_keys, as_index):
        def f(group):
            vals = group.values.sum(axis=0)
            return Series(vals, index=Index(np.arange(group.values.shape[1])))

        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", group_keys=group_keys, as_index=as_index)
        result = gb.apply(f, engine="numba")
        expected = gb.apply(lambda group: group.sum(axis=0))

        # TODO: Fix column names, needs to be exposed via the numba
        # wrappers
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("group_keys", [True, False])
    @pytest.mark.parametrize("as_index", [True, False])
    def test_cython_vs_numba_apply_reduce_axis_1(self, group_keys, as_index):
        def f(group):
            vals = group.values.sum(axis=1)
            return Series(vals, index=group.index)

        df = DataFrame({"a": [3, 2, 3, 2], "b": range(4), "c": range(1, 5)})
        gb = df.groupby("a", group_keys=group_keys, as_index=as_index)
        result = gb.apply(f, engine="numba")
        expected = gb.apply(lambda group: group.sum(axis=1))
        tm.assert_series_equal(result, expected)
