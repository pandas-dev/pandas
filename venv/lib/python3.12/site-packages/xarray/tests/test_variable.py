from __future__ import annotations

import warnings
from abc import ABC
from copy import copy, deepcopy
from datetime import datetime, timedelta
from textwrap import dedent
from typing import Any, Generic

import numpy as np
import pandas as pd
import pytest
import pytz

from xarray import DataArray, Dataset, IndexVariable, Variable, set_options
from xarray.core import dtypes, duck_array_ops, indexing
from xarray.core.common import full_like, ones_like, zeros_like
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.indexing import (
    BasicIndexer,
    CopyOnWriteArray,
    DaskIndexingAdapter,
    LazilyIndexedArray,
    MemoryCachedArray,
    NumpyIndexingAdapter,
    OuterIndexer,
    PandasIndexingAdapter,
    VectorizedIndexer,
)
from xarray.core.types import T_DuckArray
from xarray.core.utils import NDArrayMixin
from xarray.core.variable import as_compatible_data, as_variable
from xarray.namedarray.pycompat import array_type
from xarray.tests import (
    IndexableArray,
    assert_allclose,
    assert_array_equal,
    assert_equal,
    assert_identical,
    assert_no_warnings,
    has_dask_ge_2024_11_0,
    has_pandas_3,
    raise_if_dask_computes,
    requires_bottleneck,
    requires_cupy,
    requires_dask,
    requires_pint,
    requires_sparse,
    source_ndarray,
)
from xarray.tests.test_namedarray import NamedArraySubclassobjects

dask_array_type = array_type("dask")

_PAD_XR_NP_ARGS = [
    [{"x": (2, 1)}, ((2, 1), (0, 0), (0, 0))],
    [{"x": 1}, ((1, 1), (0, 0), (0, 0))],
    [{"y": (0, 3)}, ((0, 0), (0, 3), (0, 0))],
    [{"x": (3, 1), "z": (2, 0)}, ((3, 1), (0, 0), (2, 0))],
    [{"x": (3, 1), "z": 2}, ((3, 1), (0, 0), (2, 2))],
]


@pytest.fixture
def var():
    return Variable(dims=list("xyz"), data=np.random.rand(3, 4, 5))


@pytest.mark.parametrize(
    "data",
    [
        np.array(["a", "bc", "def"], dtype=object),
        np.array(["2019-01-01", "2019-01-02", "2019-01-03"], dtype="datetime64[ns]"),
    ],
)
def test_as_compatible_data_writeable(data):
    # In pandas 3 the mode.copy_on_write option defaults to True, so the option
    # setting logic can be removed once our minimum version of pandas is
    # greater than or equal to 3.
    if not has_pandas_3:
        pd.set_option("mode.copy_on_write", True)
    # GH8843, ensure writeable arrays for data_vars even with
    # pandas copy-on-write mode
    assert as_compatible_data(data).flags.writeable
    if not has_pandas_3:
        pd.reset_option("mode.copy_on_write")


class VariableSubclassobjects(NamedArraySubclassobjects, ABC):
    @pytest.fixture
    def target(self, data):
        data = 0.5 * np.arange(10).reshape(2, 5)
        return Variable(["x", "y"], data)

    def test_getitem_dict(self):
        v = self.cls(["x"], np.random.randn(5))
        actual = v[{"x": 0}]
        expected = v[0]
        assert_identical(expected, actual)

    def test_getitem_1d(self):
        data = np.array([0, 1, 2])
        v = self.cls(["x"], data)

        v_new = v[dict(x=[0, 1])]
        assert v_new.dims == ("x",)
        assert_array_equal(v_new, data[[0, 1]])

        v_new = v[dict(x=slice(None))]
        assert v_new.dims == ("x",)
        assert_array_equal(v_new, data)

        v_new = v[dict(x=Variable("a", [0, 1]))]
        assert v_new.dims == ("a",)
        assert_array_equal(v_new, data[[0, 1]])

        v_new = v[dict(x=1)]
        assert v_new.dims == ()
        assert_array_equal(v_new, data[1])

        # tuple argument
        v_new = v[slice(None)]
        assert v_new.dims == ("x",)
        assert_array_equal(v_new, data)

    def test_getitem_1d_fancy(self):
        v = self.cls(["x"], [0, 1, 2])
        # 1d-variable should be indexable by multi-dimensional Variable
        ind = Variable(("a", "b"), [[0, 1], [0, 1]])
        v_new = v[ind]
        assert v_new.dims == ("a", "b")
        expected = np.array(v._data)[([0, 1], [0, 1]), ...]
        assert_array_equal(v_new, expected)

        # boolean indexing
        ind = Variable(("x",), [True, False, True])
        v_new = v[ind]
        assert_identical(v[[0, 2]], v_new)
        v_new = v[[True, False, True]]
        assert_identical(v[[0, 2]], v_new)

        with pytest.raises(IndexError, match=r"Boolean indexer should"):
            ind = Variable(("a",), [True, False, True])
            v[ind]

    def test_getitem_with_mask(self):
        v = self.cls(["x"], [0, 1, 2])
        assert_identical(v._getitem_with_mask(-1), Variable((), np.nan))
        assert_identical(
            v._getitem_with_mask([0, -1, 1]), self.cls(["x"], [0, np.nan, 1])
        )
        assert_identical(v._getitem_with_mask(slice(2)), self.cls(["x"], [0, 1]))
        assert_identical(
            v._getitem_with_mask([0, -1, 1], fill_value=-99),
            self.cls(["x"], [0, -99, 1]),
        )

    def test_getitem_with_mask_size_zero(self):
        v = self.cls(["x"], [])
        assert_identical(v._getitem_with_mask(-1), Variable((), np.nan))
        assert_identical(
            v._getitem_with_mask([-1, -1, -1]),
            self.cls(["x"], [np.nan, np.nan, np.nan]),
        )

    def test_getitem_with_mask_nd_indexer(self):
        v = self.cls(["x"], [0, 1, 2])
        indexer = Variable(("x", "y"), [[0, -1], [-1, 2]])
        assert_identical(v._getitem_with_mask(indexer, fill_value=-1), indexer)

    def _assertIndexedLikeNDArray(self, variable, expected_value0, expected_dtype=None):
        """Given a 1-dimensional variable, verify that the variable is indexed
        like a numpy.ndarray.
        """
        assert variable[0].shape == ()
        assert variable[0].ndim == 0
        assert variable[0].size == 1
        # test identity
        assert variable.equals(variable.copy())
        assert variable.identical(variable.copy())
        # check value is equal for both ndarray and Variable
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "In the future, 'NAT == x'")
            np.testing.assert_equal(variable.values[0], expected_value0)
            np.testing.assert_equal(variable[0].values, expected_value0)
        # check type or dtype is consistent for both ndarray and Variable
        if expected_dtype is None:
            # check output type instead of array dtype
            assert type(variable.values[0]) is type(expected_value0)
            assert type(variable[0].values) is type(expected_value0)
        elif expected_dtype is not False:
            assert variable.values[0].dtype == expected_dtype
            assert variable[0].values.dtype == expected_dtype

    def test_index_0d_int(self):
        for value, dtype in [(0, np.int_), (np.int32(0), np.int32)]:
            x = self.cls(["x"], [value])
            self._assertIndexedLikeNDArray(x, value, dtype)

    def test_index_0d_float(self):
        for value, dtype in [(0.5, float), (np.float32(0.5), np.float32)]:
            x = self.cls(["x"], [value])
            self._assertIndexedLikeNDArray(x, value, dtype)

    def test_index_0d_string(self):
        value = "foo"
        dtype = np.dtype("U3")
        x = self.cls(["x"], [value])
        self._assertIndexedLikeNDArray(x, value, dtype)

    def test_index_0d_datetime(self):
        d = datetime(2000, 1, 1)
        x = self.cls(["x"], [d])
        self._assertIndexedLikeNDArray(x, np.datetime64(d))

        x = self.cls(["x"], [np.datetime64(d)])
        self._assertIndexedLikeNDArray(x, np.datetime64(d), "datetime64[us]")

        expected_unit = "us" if has_pandas_3 else "ns"
        x = self.cls(["x"], pd.DatetimeIndex([d]))
        self._assertIndexedLikeNDArray(
            x, np.datetime64(d), f"datetime64[{expected_unit}]"
        )

    def test_index_0d_timedelta64(self):
        td = timedelta(hours=1)
        # todo: discussion needed
        x = self.cls(["x"], [np.timedelta64(td)])
        self._assertIndexedLikeNDArray(
            x, np.timedelta64(td), np.dtype("timedelta64[us]")
        )

        x = self.cls(["x"], pd.to_timedelta([td]).as_unit("ns"))
        self._assertIndexedLikeNDArray(x, np.timedelta64(td), "timedelta64[ns]")

    def test_index_0d_not_a_time(self):
        d = np.datetime64("NaT", "ns")
        x = self.cls(["x"], [d])
        self._assertIndexedLikeNDArray(x, d)

    def test_index_0d_object(self):
        class HashableItemWrapper:
            def __init__(self, item):
                self.item = item

            def __eq__(self, other):
                return self.item == other.item

            def __hash__(self):
                return hash(self.item)

            def __repr__(self):
                return f"{type(self).__name__}(item={self.item!r})"

        item = HashableItemWrapper((1, 2, 3))
        x = self.cls("x", [item])
        self._assertIndexedLikeNDArray(x, item, expected_dtype=False)

    def test_0d_object_array_with_list(self):
        listarray = np.empty((1,), dtype=object)
        listarray[0] = [1, 2, 3]
        x = self.cls("x", listarray)
        assert_array_equal(x.data, listarray)
        assert_array_equal(x[0].data, listarray.squeeze())
        assert_array_equal(x.squeeze().data, listarray.squeeze())

    def test_index_and_concat_datetime(self):
        # regression test for #125
        date_range = pd.date_range("2011-09-01", periods=10)
        for dates in [date_range, date_range.values, date_range.to_pydatetime()]:
            expected = self.cls("t", dates)
            for times in [
                [expected[i] for i in range(10)],
                [expected[i : (i + 1)] for i in range(10)],
                [expected[[i]] for i in range(10)],
            ]:
                actual = Variable.concat(times, "t")
                assert expected.dtype == actual.dtype
                assert_array_equal(expected, actual)

    def test_0d_time_data(self):
        # regression test for #105
        x = self.cls("time", pd.date_range("2000-01-01", periods=5))
        expected = np.datetime64("2000-01-01", "ns")
        assert x[0].values == expected

    dt64_data = pd.date_range("1970-01-01", periods=3, unit="ns")

    @pytest.mark.parametrize(
        "values, unit",
        [
            (dt64_data, "ns"),
            (dt64_data.values, "ns"),
            (dt64_data.values.astype("datetime64[m]"), "s"),
            (dt64_data.values.astype("datetime64[s]"), "s"),
            (dt64_data.values.astype("datetime64[ps]"), "ns"),
            (
                dt64_data.to_pydatetime(),
                "us" if has_pandas_3 else "ns",
            ),
        ],
    )
    def test_datetime64_conversion(self, values, unit):
        v = self.cls(["t"], values)
        assert v.dtype == np.dtype(f"datetime64[{unit}]")
        assert_array_equal(v.values, self.dt64_data.values)
        assert v.values.dtype == np.dtype(f"datetime64[{unit}]")

    td64_data = pd.timedelta_range(start=0, periods=3)

    @pytest.mark.parametrize(
        "values, unit",
        [
            (td64_data, "ns"),
            (td64_data.values, "ns"),
            (td64_data.values.astype("timedelta64[m]"), "s"),
            (td64_data.values.astype("timedelta64[s]"), "s"),
            (td64_data.values.astype("timedelta64[ps]"), "ns"),
            (td64_data.to_pytimedelta(), "us" if has_pandas_3 else "ns"),
        ],
    )
    def test_timedelta64_conversion(self, values, unit):
        v = self.cls(["t"], values)
        assert v.dtype == np.dtype(f"timedelta64[{unit}]")
        assert_array_equal(v.values, self.td64_data.values)
        assert v.values.dtype == np.dtype(f"timedelta64[{unit}]")

    def test_object_conversion(self):
        data = np.arange(5).astype(str).astype(object)
        actual = self.cls("x", data)
        assert actual.dtype == data.dtype

    def test_pandas_data(self):
        v = self.cls(["x"], pd.Series([0, 1, 2], index=[3, 2, 1]))
        assert_identical(v, v[[0, 1, 2]])
        v = self.cls(["x"], pd.Index([0, 1, 2]))
        assert v[0].values == v.values[0]

    def test_pandas_period_index(self):
        v = self.cls(["x"], pd.period_range(start="2000", periods=20, freq="D"))
        v = v.load()  # for dask-based Variable
        assert v[0] == pd.Period("2000", freq="D")
        assert "PeriodArray" in repr(v)

    @pytest.mark.parametrize("dtype", [float, int])
    def test_1d_math(self, dtype: np.typing.DTypeLike | None) -> None:
        x = np.arange(5, dtype=dtype)
        y = np.ones(5, dtype=dtype)

        # should we need `.to_base_variable()`?
        # probably a break that `+v` changes type?
        v = self.cls(["x"], x)
        base_v = v.to_base_variable()
        # unary ops
        assert_identical(base_v, +v)
        assert_identical(base_v, abs(v))
        assert_array_equal((-v).values, -x)
        # binary ops with numbers
        assert_identical(base_v, v + 0)
        assert_identical(base_v, 0 + v)
        assert_identical(base_v, v * 1)
        if dtype is int:
            assert_identical(base_v, v << 0)
            assert_array_equal(v << 3, x << 3)
            assert_array_equal(v >> 2, x >> 2)
        # binary ops with numpy arrays
        assert_array_equal((v * x).values, x**2)
        assert_array_equal((x * v).values, x**2)
        assert_array_equal(v - y, v - 1)
        assert_array_equal(y - v, 1 - v)
        if dtype is int:
            assert_array_equal(v << x, x << x)
            assert_array_equal(v >> x, x >> x)
        # verify attributes are dropped
        v2 = self.cls(["x"], x, {"units": "meters"})
        with set_options(keep_attrs=False):
            assert_identical(base_v, +v2)
        # binary ops with all variables
        assert_array_equal(v + v, 2 * v)
        w = self.cls(["x"], y, {"foo": "bar"})
        # With drop_conflicts, v (no attrs) + w (has attrs) should keep w's attrs
        # Note: IndexVariable ops return Variable, not IndexVariable
        expected = self.cls(["x"], x + y, {"foo": "bar"}).to_base_variable()
        assert_identical(v + w, expected)
        assert_array_equal((v * w).values, x * y)

        # something complicated
        assert_array_equal((v**2 * w - 1 + x).values, x**2 * y - 1 + x)
        # make sure dtype is preserved (for Index objects)
        assert dtype == (+v).dtype
        assert dtype == (+v).values.dtype
        assert dtype == (0 + v).dtype
        assert dtype == (0 + v).values.dtype
        # check types of returned data
        assert isinstance(+v, Variable)
        assert not isinstance(+v, IndexVariable)
        assert isinstance(0 + v, Variable)
        assert not isinstance(0 + v, IndexVariable)

    def test_1d_reduce(self):
        x = np.arange(5)
        v = self.cls(["x"], x)
        actual = v.sum()
        expected = Variable((), 10)
        assert_identical(expected, actual)
        assert type(actual) is Variable

    def test_array_interface(self):
        x = np.arange(5)
        v = self.cls(["x"], x)
        assert_array_equal(np.asarray(v), x)
        # test patched in methods
        assert_array_equal(v.astype(float), x.astype(float))
        # think this is a break, that argsort changes the type
        assert_identical(v.argsort(), v.to_base_variable())
        assert_identical(v.clip(2, 3), self.cls("x", x.clip(2, 3)).to_base_variable())
        # test ufuncs
        assert_identical(np.sin(v), self.cls(["x"], np.sin(x)).to_base_variable())
        assert isinstance(np.sin(v), Variable)
        assert not isinstance(np.sin(v), IndexVariable)

    def example_1d_objects(self):
        for data in [
            range(3),
            0.5 * np.arange(3),
            0.5 * np.arange(3, dtype=np.float32),
            pd.date_range("2000-01-01", periods=3),
            np.array(["a", "b", "c"], dtype=object),
        ]:
            yield (self.cls("x", data), data)

    def test___array__(self):
        for v, data in self.example_1d_objects():
            assert_array_equal(v.values, np.asarray(data))
            assert_array_equal(np.asarray(v), np.asarray(data))
            assert v[0].values == np.asarray(data)[0]
            assert np.asarray(v[0]) == np.asarray(data)[0]

    def test_equals_all_dtypes(self):
        for v, _ in self.example_1d_objects():
            v2 = v.copy()
            assert v.equals(v2)
            assert v.identical(v2)
            assert v.no_conflicts(v2)
            assert v[0].equals(v2[0])
            assert v[0].identical(v2[0])
            assert v[0].no_conflicts(v2[0])
            assert v[:2].equals(v2[:2])
            assert v[:2].identical(v2[:2])
            assert v[:2].no_conflicts(v2[:2])

    def test_eq_all_dtypes(self):
        # ensure that we don't choke on comparisons for which numpy returns
        # scalars
        expected = Variable("x", 3 * [False])
        for v, _ in self.example_1d_objects():
            actual = "z" == v
            assert_identical(expected, actual)
            actual = ~("z" != v)
            assert_identical(expected, actual)

    def test_encoding_preserved(self):
        expected = self.cls("x", range(3), {"foo": 1}, {"bar": 2})
        for actual in [
            expected.T,
            expected[...],
            expected.squeeze(),
            expected.isel(x=slice(None)),
            expected.set_dims({"x": 3}),
            expected.copy(deep=True),
            expected.copy(deep=False),
        ]:
            assert_identical(expected.to_base_variable(), actual.to_base_variable())
            assert expected.encoding == actual.encoding

    def test_drop_encoding(self) -> None:
        encoding1 = {"scale_factor": 1}
        # encoding set via cls constructor
        v1 = self.cls(["a"], [0, 1, 2], encoding=encoding1)
        assert v1.encoding == encoding1
        v2 = v1.drop_encoding()
        assert v1.encoding == encoding1
        assert v2.encoding == {}

        # encoding set via setter
        encoding3 = {"scale_factor": 10}
        v3 = self.cls(["a"], [0, 1, 2], encoding=encoding3)
        assert v3.encoding == encoding3
        v4 = v3.drop_encoding()
        assert v3.encoding == encoding3
        assert v4.encoding == {}

    def test_concat(self):
        x = np.arange(5)
        y = np.arange(5, 10)
        v = self.cls(["a"], x)
        w = self.cls(["a"], y)
        assert_identical(
            Variable(["b", "a"], np.array([x, y])), Variable.concat([v, w], "b")
        )
        assert_identical(
            Variable(["b", "a"], np.array([x, y])), Variable.concat((v, w), "b")
        )
        assert_identical(
            Variable(["b", "a"], np.array([x, y])), Variable.concat((v, w), "b")
        )
        with pytest.raises(ValueError, match=r"Variable has dimensions"):
            Variable.concat([v, Variable(["c"], y)], "b")
        # test indexers
        actual = Variable.concat(
            [v, w], positions=[np.arange(0, 10, 2), np.arange(1, 10, 2)], dim="a"
        )
        expected = Variable("a", np.array([x, y]).ravel(order="F"))
        assert_identical(expected, actual)
        # test concatenating along a dimension
        v = Variable(["time", "x"], np.random.random((10, 8)))
        assert_identical(v, Variable.concat([v[:5], v[5:]], "time"))
        assert_identical(v, Variable.concat([v[:5], v[5:6], v[6:]], "time"))
        assert_identical(v, Variable.concat([v[:1], v[1:]], "time"))
        # test dimension order
        assert_identical(v, Variable.concat([v[:, :5], v[:, 5:]], "x"))
        with pytest.raises(ValueError, match=r"all input arrays must have"):
            Variable.concat([v[:, 0], v[:, 1:]], "x")

    def test_concat_attrs(self):
        # always keep attrs from first variable
        v = self.cls("a", np.arange(5), {"foo": "bar"})
        w = self.cls("a", np.ones(5))
        expected = self.cls(
            "a", np.concatenate([np.arange(5), np.ones(5)])
        ).to_base_variable()
        expected.attrs["foo"] = "bar"
        assert_identical(expected, Variable.concat([v, w], "a"))

    def test_concat_fixed_len_str(self):
        # regression test for #217
        for kind in ["S", "U"]:
            x = self.cls("animal", np.array(["horse"], dtype=kind))
            y = self.cls("animal", np.array(["aardvark"], dtype=kind))
            actual = Variable.concat([x, y], "animal")
            expected = Variable("animal", np.array(["horse", "aardvark"], dtype=kind))
            assert_equal(expected, actual)

    def test_concat_number_strings(self):
        # regression test for #305
        a = self.cls("x", ["0", "1", "2"])
        b = self.cls("x", ["3", "4"])
        actual = Variable.concat([a, b], dim="x")
        expected = Variable("x", np.arange(5).astype(str))
        assert_identical(expected, actual)
        assert actual.dtype.kind == expected.dtype.kind

    def test_concat_mixed_dtypes(self):
        a = self.cls("x", [0, 1])
        b = self.cls("x", ["two"])
        actual = Variable.concat([a, b], dim="x")
        expected = Variable("x", np.array([0, 1, "two"], dtype=object))
        assert_identical(expected, actual)
        assert actual.dtype == object

    @pytest.mark.parametrize("deep", [True, False])
    @pytest.mark.parametrize("astype", [float, int, str])
    def test_copy(self, deep: bool, astype: type[object]) -> None:
        v = self.cls("x", (0.5 * np.arange(10)).astype(astype), {"foo": "bar"})
        w = v.copy(deep=deep)
        assert type(v) is type(w)
        assert_identical(v, w)
        assert v.dtype == w.dtype
        if self.cls is Variable:
            if deep:
                assert source_ndarray(v.values) is not source_ndarray(w.values)
            else:
                assert source_ndarray(v.values) is source_ndarray(w.values)
        assert_identical(v, copy(v))

    def test_copy_deep_recursive(self) -> None:
        # GH:issue:7111

        # direct recursion
        v = self.cls("x", [0, 1])
        v.attrs["other"] = v

        # TODO: cannot use assert_identical on recursive Vars yet...
        # lets just ensure that deep copy works without RecursionError
        v.copy(deep=True)

        # indirect recursion
        v2 = self.cls("y", [2, 3])
        v.attrs["other"] = v2
        v2.attrs["other"] = v

        # TODO: cannot use assert_identical on recursive Vars yet...
        # lets just ensure that deep copy works without RecursionError
        v.copy(deep=True)
        v2.copy(deep=True)

    def test_copy_index(self):
        midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2], [-1, -2]], names=("one", "two", "three")
        )
        v = self.cls("x", midx)
        for deep in [True, False]:
            w = v.copy(deep=deep)
            assert isinstance(w._data, PandasIndexingAdapter)
            assert isinstance(w.to_index(), pd.MultiIndex)
            assert_array_equal(v._data.array, w._data.array)

    def test_copy_with_data(self) -> None:
        orig = Variable(("x", "y"), [[1.5, 2.0], [3.1, 4.3]], {"foo": "bar"})
        new_data = np.array([[2.5, 5.0], [7.1, 43]])
        actual = orig.copy(data=new_data)
        expected = orig.copy()
        expected.data = new_data
        assert_identical(expected, actual)

    def test_copy_with_data_errors(self) -> None:
        orig = Variable(("x", "y"), [[1.5, 2.0], [3.1, 4.3]], {"foo": "bar"})
        new_data = [2.5, 5.0]
        with pytest.raises(ValueError, match=r"must match shape of object"):
            orig.copy(data=new_data)  # type: ignore[arg-type]

    def test_copy_index_with_data(self) -> None:
        orig = IndexVariable("x", np.arange(5))
        new_data = np.arange(5, 10)
        actual = orig.copy(data=new_data)
        expected = IndexVariable("x", np.arange(5, 10))
        assert_identical(expected, actual)

    def test_copy_index_with_data_errors(self) -> None:
        orig = IndexVariable("x", np.arange(5))
        new_data = np.arange(5, 20)
        with pytest.raises(ValueError, match=r"must match shape of object"):
            orig.copy(data=new_data)
        with pytest.raises(ValueError, match=r"Cannot assign to the .data"):
            orig.data = new_data
        with pytest.raises(ValueError, match=r"Cannot assign to the .values"):
            orig.values = new_data

    def test_replace(self):
        var = Variable(("x", "y"), [[1.5, 2.0], [3.1, 4.3]], {"foo": "bar"})
        result = var._replace()
        assert_identical(result, var)

        new_data = np.arange(4).reshape(2, 2)
        result = var._replace(data=new_data)
        assert_array_equal(result.data, new_data)

    def test_real_and_imag(self):
        v = self.cls("x", np.arange(3) - 1j * np.arange(3), {"foo": "bar"})
        expected_re = self.cls("x", np.arange(3), {"foo": "bar"})
        assert_identical(v.real, expected_re)

        expected_im = self.cls("x", -np.arange(3), {"foo": "bar"})
        assert_identical(v.imag, expected_im)

        expected_abs = self.cls("x", np.sqrt(2 * np.arange(3) ** 2)).to_base_variable()
        assert_allclose(abs(v), expected_abs)

    def test_aggregate_complex(self):
        # should skip NaNs
        v = self.cls("x", [1, 2j, np.nan])
        expected = Variable((), 0.5 + 1j)
        assert_allclose(v.mean(), expected)

    def test_pandas_categorical_dtype(self):
        data = pd.Categorical(np.arange(10, dtype="int64"))
        v = self.cls("x", data)
        print(v)  # should not error
        assert v.dtype == data.dtype

    def test_pandas_datetime64_with_tz(self):
        data = pd.date_range(
            start="2000-01-01",
            tz=pytz.timezone("America/New_York"),
            periods=10,
            freq="1h",
        )
        v = self.cls("x", data)
        print(v)  # should not error
        if v.dtype == np.dtype("O"):
            import dask.array as da

            assert isinstance(v.data, da.Array)
        else:
            assert v.dtype == data.dtype

    def test_multiindex(self):
        idx = pd.MultiIndex.from_product([list("abc"), [0, 1]])
        v = self.cls("x", idx)
        assert_identical(Variable((), ("a", 0)), v[0])
        assert_identical(v, v[:])

    def test_load(self):
        array = self.cls("x", np.arange(5))
        orig_data = array._data
        copied = array.copy(deep=True)
        if array.chunks is None:
            array.load()
            assert type(array._data) is type(orig_data)
            assert type(copied._data) is type(orig_data)
            assert_identical(array, copied)

    def test_getitem_advanced(self):
        v = self.cls(["x", "y"], [[0, 1, 2], [3, 4, 5]])
        v_data = v.compute().data

        # orthogonal indexing
        v_new = v[([0, 1], [1, 0])]
        assert v_new.dims == ("x", "y")
        assert_array_equal(v_new, v_data[[0, 1]][:, [1, 0]])

        v_new = v[[0, 1]]
        assert v_new.dims == ("x", "y")
        assert_array_equal(v_new, v_data[[0, 1]])

        # with mixed arguments
        ind = Variable(["a"], [0, 1])
        v_new = v[dict(x=[0, 1], y=ind)]
        assert v_new.dims == ("x", "a")
        assert_array_equal(v_new, v_data[[0, 1]][:, [0, 1]])

        # boolean indexing
        v_new = v[dict(x=[True, False], y=[False, True, False])]
        assert v_new.dims == ("x", "y")
        assert_array_equal(v_new, v_data[0][1])

        # with scalar variable
        ind = Variable((), 2)
        v_new = v[dict(y=ind)]
        expected = v[dict(y=2)]
        assert_array_equal(v_new, expected)

        # with boolean variable with wrong shape
        ind2: np.ndarray[Any, np.dtype[np.bool_]] = np.array([True, False])
        with pytest.raises(IndexError, match=r"Boolean array size 2 is "):
            v[Variable(("a", "b"), [[0, 1]]), ind2]

        # boolean indexing with different dimension
        ind = Variable(["a"], [True, False, False])
        with pytest.raises(IndexError, match=r"Boolean indexer should be"):
            v[dict(y=ind)]

    def test_getitem_uint_1d(self):
        # regression test for #1405
        v = self.cls(["x"], [0, 1, 2])
        v_data = v.compute().data

        v_new = v[np.array([0])]
        assert_array_equal(v_new, v_data[0])
        v_new = v[np.array([0], dtype="uint64")]
        assert_array_equal(v_new, v_data[0])

    def test_getitem_uint(self):
        # regression test for #1405
        v = self.cls(["x", "y"], [[0, 1, 2], [3, 4, 5]])
        v_data = v.compute().data

        v_new = v[np.array([0])]
        assert_array_equal(v_new, v_data[[0], :])
        v_new = v[np.array([0], dtype="uint64")]
        assert_array_equal(v_new, v_data[[0], :])

        v_new = v[np.uint64(0)]
        assert_array_equal(v_new, v_data[0, :])

    def test_getitem_0d_array(self):
        # make sure 0d-np.array can be used as an indexer
        v = self.cls(["x"], [0, 1, 2])
        v_data = v.compute().data

        v_new = v[np.array([0])[0]]
        assert_array_equal(v_new, v_data[0])

        v_new = v[np.array(0)]
        assert_array_equal(v_new, v_data[0])

        v_new = v[Variable((), np.array(0))]
        assert_array_equal(v_new, v_data[0])

    def test_getitem_fancy(self):
        v = self.cls(["x", "y"], [[0, 1, 2], [3, 4, 5]])
        v_data = v.compute().data

        ind = Variable(["a", "b"], [[0, 1, 1], [1, 1, 0]])
        v_new = v[ind]
        assert v_new.dims == ("a", "b", "y")
        assert_array_equal(v_new, v_data[[[0, 1, 1], [1, 1, 0]], :])

        # It would be ok if indexed with the multi-dimensional array including
        # the same name
        ind = Variable(["x", "b"], [[0, 1, 1], [1, 1, 0]])
        v_new = v[ind]
        assert v_new.dims == ("x", "b", "y")
        assert_array_equal(v_new, v_data[[[0, 1, 1], [1, 1, 0]], :])

        ind = Variable(["a", "b"], [[0, 1, 2], [2, 1, 0]])
        v_new = v[dict(y=ind)]
        assert v_new.dims == ("x", "a", "b")
        assert_array_equal(v_new, v_data[:, ([0, 1, 2], [2, 1, 0])])

        ind = Variable(["a", "b"], [[0, 0], [1, 1]])
        v_new = v[dict(x=[1, 0], y=ind)]
        assert v_new.dims == ("x", "a", "b")
        assert_array_equal(v_new, v_data[[1, 0]][:, ind])

        # along diagonal
        ind = Variable(["a"], [0, 1])
        v_new = v[ind, ind]
        assert v_new.dims == ("a",)
        assert_array_equal(v_new, v_data[[0, 1], [0, 1]])

        # with integer
        ind = Variable(["a", "b"], [[0, 0], [1, 1]])
        v_new = v[dict(x=0, y=ind)]
        assert v_new.dims == ("a", "b")
        assert_array_equal(v_new[0], v_data[0][[0, 0]])
        assert_array_equal(v_new[1], v_data[0][[1, 1]])

        # with slice
        ind = Variable(["a", "b"], [[0, 0], [1, 1]])
        v_new = v[dict(x=slice(None), y=ind)]
        assert v_new.dims == ("x", "a", "b")
        assert_array_equal(v_new, v_data[:, [[0, 0], [1, 1]]])

        ind = Variable(["a", "b"], [[0, 0], [1, 1]])
        v_new = v[dict(x=ind, y=slice(None))]
        assert v_new.dims == ("a", "b", "y")
        assert_array_equal(v_new, v_data[[[0, 0], [1, 1]], :])

        ind = Variable(["a", "b"], [[0, 0], [1, 1]])
        v_new = v[dict(x=ind, y=slice(None, 1))]
        assert v_new.dims == ("a", "b", "y")
        assert_array_equal(v_new, v_data[[[0, 0], [1, 1]], slice(None, 1)])

        # slice matches explicit dimension
        ind = Variable(["y"], [0, 1])
        v_new = v[ind, :2]
        assert v_new.dims == ("y",)
        assert_array_equal(v_new, v_data[[0, 1], [0, 1]])

        # with multiple slices
        v = self.cls(["x", "y", "z"], [[[1, 2, 3], [4, 5, 6]]])
        ind = Variable(["a", "b"], [[0]])
        v_new = v[ind, :, :]
        expected = Variable(["a", "b", "y", "z"], v.data[np.newaxis, ...])
        assert_identical(v_new, expected)

        v = Variable(["w", "x", "y", "z"], [[[[1, 2, 3], [4, 5, 6]]]])
        ind = Variable(["y"], [0])
        v_new = v[ind, :, 1:2, 2]
        expected = Variable(["y", "x"], [[6]])
        assert_identical(v_new, expected)

        # slice and vector mixed indexing resulting in the same dimension
        v = Variable(["x", "y", "z"], np.arange(60).reshape(3, 4, 5))
        ind = Variable(["x"], [0, 1, 2])
        v_new = v[:, ind]
        expected = Variable(("x", "z"), np.zeros((3, 5)))
        expected[0] = v.data[0, 0]
        expected[1] = v.data[1, 1]
        expected[2] = v.data[2, 2]
        assert_identical(v_new, expected)

        v_new = v[:, ind.data]
        assert v_new.shape == (3, 3, 5)

    def test_getitem_error(self):
        v = self.cls(["x", "y"], [[0, 1, 2], [3, 4, 5]])

        with pytest.raises(IndexError, match=r"labeled multi-"):
            v[[[0, 1], [1, 2]]]

        ind_x = Variable(["a"], [0, 1, 1])
        ind_y = Variable(["a"], [0, 1])
        with pytest.raises(IndexError, match=r"Dimensions of indexers "):
            v[ind_x, ind_y]

        ind = Variable(["a", "b"], [[True, False], [False, True]])
        with pytest.raises(IndexError, match=r"2-dimensional boolean"):
            v[dict(x=ind)]

        v = Variable(["x", "y", "z"], np.arange(60).reshape(3, 4, 5))
        ind = Variable(["x"], [0, 1])
        with pytest.raises(IndexError, match=r"Dimensions of indexers mismatch"):
            v[:, ind]

    @pytest.mark.parametrize(
        "mode",
        [
            "mean",
            "median",
            "reflect",
            "edge",
            "linear_ramp",
            "maximum",
            "minimum",
            "symmetric",
            "wrap",
        ],
    )
    @pytest.mark.parametrize("xr_arg, np_arg", _PAD_XR_NP_ARGS)
    @pytest.mark.filterwarnings(
        r"ignore:dask.array.pad.+? converts integers to floats."
    )
    def test_pad(self, mode, xr_arg, np_arg):
        data = np.arange(4 * 3 * 2).reshape(4, 3, 2)
        v = self.cls(["x", "y", "z"], data)

        actual = v.pad(mode=mode, **xr_arg)
        expected = np.pad(data, np_arg, mode=mode)

        assert_array_equal(actual, expected)
        assert isinstance(actual._data, type(v._data))

    @pytest.mark.parametrize("xr_arg, np_arg", _PAD_XR_NP_ARGS)
    def test_pad_constant_values(self, xr_arg, np_arg):
        data = np.arange(4 * 3 * 2).reshape(4, 3, 2)
        v = self.cls(["x", "y", "z"], data)

        actual = v.pad(**xr_arg)
        expected = np.pad(
            np.array(duck_array_ops.astype(v.data, float)),
            np_arg,
            mode="constant",
            constant_values=np.nan,
        )
        assert_array_equal(actual, expected)
        assert isinstance(actual._data, type(v._data))

        # for the boolean array, we pad False
        data = np.full_like(data, False, dtype=bool).reshape(4, 3, 2)
        v = self.cls(["x", "y", "z"], data)

        actual = v.pad(mode="constant", constant_values=False, **xr_arg)
        expected = np.pad(
            np.array(v.data), np_arg, mode="constant", constant_values=False
        )
        assert_array_equal(actual, expected)

    @pytest.mark.parametrize(
        ["keep_attrs", "attrs", "expected"],
        [
            pytest.param(None, {"a": 1, "b": 2}, {"a": 1, "b": 2}, id="default"),
            pytest.param(False, {"a": 1, "b": 2}, {}, id="False"),
            pytest.param(True, {"a": 1, "b": 2}, {"a": 1, "b": 2}, id="True"),
        ],
    )
    def test_pad_keep_attrs(self, keep_attrs, attrs, expected):
        data = np.arange(10, dtype=float)
        v = self.cls(["x"], data, attrs)

        keep_attrs_ = "default" if keep_attrs is None else keep_attrs

        with set_options(keep_attrs=keep_attrs_):
            actual = v.pad({"x": (1, 1)}, mode="constant", constant_values=np.nan)

            assert actual.attrs == expected

        actual = v.pad(
            {"x": (1, 1)},
            mode="constant",
            constant_values=np.nan,
            keep_attrs=keep_attrs,
        )
        assert actual.attrs == expected

    @pytest.mark.parametrize("d, w", (("x", 3), ("y", 5)))
    def test_rolling_window(self, d, w):
        # Just a working test. See test_nputils for the algorithm validation
        v = self.cls(["x", "y", "z"], np.arange(40 * 30 * 2).reshape(40, 30, 2))
        v_rolling = v.rolling_window(d, w, d + "_window")
        assert v_rolling.dims == ("x", "y", "z", d + "_window")
        assert v_rolling.shape == v.shape + (w,)

        v_rolling = v.rolling_window(d, w, d + "_window", center=True)
        assert v_rolling.dims == ("x", "y", "z", d + "_window")
        assert v_rolling.shape == v.shape + (w,)

        # dask and numpy result should be the same
        v_loaded = v.load().rolling_window(d, w, d + "_window", center=True)
        assert_array_equal(v_rolling, v_loaded)

        # numpy backend should not be over-written
        if isinstance(v._data, np.ndarray):
            with pytest.raises(ValueError):
                v_loaded[0] = 1.0

    def test_rolling_1d(self):
        x = self.cls("x", np.array([1, 2, 3, 4], dtype=float))

        kwargs = dict(dim="x", window=3, window_dim="xw")
        actual = x.rolling_window(**kwargs, center=True, fill_value=np.nan)
        expected = Variable(
            ("x", "xw"),
            np.array(
                [[np.nan, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, np.nan]], dtype=float
            ),
        )
        assert_equal(actual, expected)

        actual = x.rolling_window(**kwargs, center=False, fill_value=0.0)
        expected = self.cls(
            ("x", "xw"),
            np.array([[0, 0, 1], [0, 1, 2], [1, 2, 3], [2, 3, 4]], dtype=float),
        )
        assert_equal(actual, expected)

        x = self.cls(("y", "x"), np.stack([x, x * 1.1]))
        actual = x.rolling_window(**kwargs, center=False, fill_value=0.0)
        expected = self.cls(
            ("y", "x", "xw"), np.stack([expected.data, expected.data * 1.1], axis=0)
        )
        assert_equal(actual, expected)

    @pytest.mark.parametrize("center", [[True, True], [False, False]])
    @pytest.mark.parametrize("dims", [("x", "y"), ("y", "z"), ("z", "x")])
    def test_nd_rolling(self, center, dims):
        x = self.cls(
            ("x", "y", "z"),
            np.arange(7 * 6 * 8).reshape(7, 6, 8).astype(float),
        )
        window = [3, 3]
        actual = x.rolling_window(
            dim=dims,
            window=window,
            window_dim=[f"{k}w" for k in dims],
            center=center,
            fill_value=np.nan,
        )
        expected = x
        for dim, win, cent in zip(dims, window, center, strict=True):
            expected = expected.rolling_window(
                dim=dim,
                window=win,
                window_dim=f"{dim}w",
                center=cent,
                fill_value=np.nan,
            )
        assert_equal(actual, expected)

    @pytest.mark.parametrize(
        ("dim, window, window_dim, center"),
        [
            ("x", [3, 3], "x_w", True),
            ("x", 3, ("x_w", "x_w"), True),
            ("x", 3, "x_w", [True, True]),
        ],
    )
    def test_rolling_window_errors(self, dim, window, window_dim, center):
        x = self.cls(
            ("x", "y", "z"),
            np.arange(7 * 6 * 8).reshape(7, 6, 8).astype(float),
        )
        with pytest.raises(ValueError):
            x.rolling_window(
                dim=dim,
                window=window,
                window_dim=window_dim,
                center=center,
            )


class TestVariable(VariableSubclassobjects):
    def cls(self, *args, **kwargs) -> Variable:
        return Variable(*args, **kwargs)

    @pytest.fixture(autouse=True)
    def setup(self):
        self.d = np.random.random((10, 3)).astype(np.float64)

    def test_values(self):
        v = Variable(["time", "x"], self.d)
        assert_array_equal(v.values, self.d)
        assert source_ndarray(v.values) is self.d
        with pytest.raises(ValueError):
            # wrong size
            v.values = np.random.random(5)
        d2 = np.random.random((10, 3))
        v.values = d2
        assert source_ndarray(v.values) is d2

    def test_numpy_same_methods(self):
        v = Variable([], np.float32(0.0))
        assert v.item() == 0  # type: ignore[attr-defined]
        assert type(v.item()) is float  # type: ignore[attr-defined]

        v = IndexVariable("x", np.arange(5))
        assert 2 == v.searchsorted(2)  # type: ignore[attr-defined]

    @pytest.mark.parametrize(
        "values, unit",
        [
            (np.datetime64("2000-01-01"), "s"),
            (
                pd.Timestamp("2000-01-01T00").as_unit("s"),
                "s" if has_pandas_3 else "ns",
            ),
            (
                datetime(2000, 1, 1),
                "us" if has_pandas_3 else "ns",
            ),
            (np.datetime64("2000-01-01T00:00:00.1234567891"), "ns"),
        ],
    )
    def test_datetime64_conversion_scalar(self, values, unit):
        v = Variable([], values)
        assert v.dtype == np.dtype(f"datetime64[{unit}]")
        assert np.issubdtype(v.values, "datetime64")
        assert v.values.dtype == np.dtype(f"datetime64[{unit}]")

    @pytest.mark.parametrize(
        "values, unit",
        [
            (np.timedelta64(1, "m"), "s"),
            (np.timedelta64(1, "D"), "s"),
            (np.timedelta64(1001, "ps"), "ns"),
            (pd.Timedelta("1 day").as_unit("ns"), "ns"),
            (timedelta(days=1), "us" if has_pandas_3 else "ns"),
        ],
    )
    def test_timedelta64_conversion_scalar(self, values, unit):
        v = Variable([], values)
        assert v.dtype == np.dtype(f"timedelta64[{unit}]")
        assert np.issubdtype(v.values, "timedelta64")
        assert v.values.dtype == np.dtype(f"timedelta64[{unit}]")

    def test_0d_str(self):
        v = Variable([], "foo")
        assert v.dtype == np.dtype("U3")
        assert v.values == "foo"

        v = Variable([], np.bytes_("foo"))
        assert v.dtype == np.dtype("S3")
        assert v.values == "foo".encode("ascii")

    def test_0d_datetime(self):
        v = Variable([], pd.Timestamp("2000-01-01").as_unit("s"))
        expected_unit = "s" if has_pandas_3 else "ns"
        assert v.dtype == np.dtype(f"datetime64[{expected_unit}]")
        assert v.values == np.datetime64("2000-01-01", expected_unit)  # type: ignore[call-overload]

    @pytest.mark.parametrize(
        "values, unit",
        [(pd.to_timedelta("1s").as_unit("ns"), "ns"), (np.timedelta64(1, "s"), "s")],
    )
    def test_0d_timedelta(self, values, unit):
        # todo: check, if this test is OK
        v = Variable([], values)
        assert v.dtype == np.dtype(f"timedelta64[{unit}]")
        assert v.values == np.timedelta64(10**9, "ns")

    def test_equals_and_identical(self):
        d = np.random.rand(10, 3)
        d[0, 0] = np.nan
        v1 = Variable(("dim1", "dim2"), data=d, attrs={"att1": 3, "att2": [1, 2, 3]})
        v2 = Variable(("dim1", "dim2"), data=d, attrs={"att1": 3, "att2": [1, 2, 3]})
        assert v1.equals(v2)
        assert v1.identical(v2)

        v3 = Variable(("dim1", "dim3"), data=d)
        assert not v1.equals(v3)

        v4 = Variable(("dim1", "dim2"), data=d)
        assert v1.equals(v4)
        assert not v1.identical(v4)

        v5 = deepcopy(v1)
        v5.values[:] = np.random.rand(10, 3)
        assert not v1.equals(v5)

        assert not v1.equals(None)
        assert not v1.equals(d)

        assert not v1.identical(None)
        assert not v1.identical(d)

    def test_broadcast_equals(self):
        v1 = Variable((), np.nan)
        v2 = Variable(("x"), [np.nan, np.nan])
        assert v1.broadcast_equals(v2)
        assert not v1.equals(v2)
        assert not v1.identical(v2)

        v3 = Variable(("x"), [np.nan])
        assert v1.broadcast_equals(v3)
        assert not v1.equals(v3)
        assert not v1.identical(v3)

        assert not v1.broadcast_equals(None)

        v4 = Variable(("x"), [np.nan] * 3)
        assert not v2.broadcast_equals(v4)

    def test_no_conflicts(self):
        v1 = Variable(("x"), [1, 2, np.nan, np.nan])
        v2 = Variable(("x"), [np.nan, 2, 3, np.nan])
        assert v1.no_conflicts(v2)
        assert not v1.equals(v2)
        assert not v1.broadcast_equals(v2)
        assert not v1.identical(v2)

        assert not v1.no_conflicts(None)

        v3 = Variable(("y"), [np.nan, 2, 3, np.nan])
        assert not v3.no_conflicts(v1)

        d = np.array([1, 2, np.nan, np.nan])
        assert not v1.no_conflicts(d)
        assert not v2.no_conflicts(d)

        v4 = Variable(("w", "x"), [d])
        assert v1.no_conflicts(v4)

    def test_as_variable(self):
        data = np.arange(10)
        expected = Variable("x", data)
        expected_extra = Variable(
            "x", data, attrs={"myattr": "val"}, encoding={"scale_factor": 1}
        )

        assert_identical(expected, as_variable(expected))

        ds = Dataset({"x": expected})
        var = as_variable(ds["x"]).to_base_variable()
        assert_identical(expected, var)
        assert not isinstance(ds["x"], Variable)
        assert isinstance(as_variable(ds["x"]), Variable)

        xarray_tuple = (
            expected_extra.dims,
            expected_extra.values,
            expected_extra.attrs,
            expected_extra.encoding,
        )
        assert_identical(expected_extra, as_variable(xarray_tuple))

        with pytest.raises(TypeError, match=r"tuple of form"):
            as_variable(tuple(data))
        with pytest.raises(ValueError, match=r"tuple of form"):  # GH1016
            as_variable(("five", "six", "seven"))
        with pytest.raises(TypeError, match=r"without an explicit list of dimensions"):
            as_variable(data)

        with pytest.warns(FutureWarning, match="IndexVariable"):
            actual = as_variable(data, name="x")
        assert_identical(expected.to_index_variable(), actual)

        actual = as_variable(0)
        expected = Variable([], 0)
        assert_identical(expected, actual)

        data2: np.ndarray[tuple[int, int], np.dtype[np.signedinteger[Any]]] = np.arange(
            9
        ).reshape((3, 3))
        expected = Variable(("x", "y"), data2)
        with pytest.raises(ValueError, match=r"without explicit dimension names"):
            as_variable(data2, name="x")

        # name of nD variable matches dimension name
        actual = as_variable(expected, name="x")
        assert_identical(expected, actual)

        # test datetime, timedelta conversion
        dt = np.array([datetime(1999, 1, 1) + timedelta(days=x) for x in range(10)])
        with pytest.warns(FutureWarning, match="IndexVariable"):
            assert as_variable(dt, "time").dtype.kind == "M"
        td = np.array([timedelta(days=x) for x in range(10)])
        with pytest.warns(FutureWarning, match="IndexVariable"):
            assert as_variable(td, "time").dtype.kind == "m"

        with pytest.raises(TypeError):
            as_variable(("x", DataArray([])))

    def test_repr(self):
        v = Variable(["time", "x"], [[1, 2, 3], [4, 5, 6]], {"foo": "bar"})
        v = v.astype(np.uint64)
        expected = dedent(
            """
        <xarray.Variable (time: 2, x: 3)> Size: 48B
        array([[1, 2, 3],
               [4, 5, 6]], dtype=uint64)
        Attributes:
            foo:      bar
        """
        ).strip()
        assert expected == repr(v)

    def test_repr_lazy_data(self):
        v = Variable("x", LazilyIndexedArray(np.arange(2e5)))
        assert "200000 values with dtype" in repr(v)
        assert isinstance(v._data, LazilyIndexedArray)

    def test_detect_indexer_type(self):
        """Tests indexer type was correctly detected."""
        data = np.random.random((10, 11))
        v = Variable(["x", "y"], data)

        _, ind, _ = v._broadcast_indexes((0, 1))
        assert type(ind) is indexing.BasicIndexer

        _, ind, _ = v._broadcast_indexes((0, slice(0, 8, 2)))
        assert type(ind) is indexing.BasicIndexer

        _, ind, _ = v._broadcast_indexes((0, [0, 1]))
        assert type(ind) is indexing.OuterIndexer

        _, ind, _ = v._broadcast_indexes(([0, 1], 1))
        assert type(ind) is indexing.OuterIndexer

        _, ind, _ = v._broadcast_indexes(([0, 1], [1, 2]))
        assert type(ind) is indexing.OuterIndexer

        _, ind, _ = v._broadcast_indexes(([0, 1], slice(0, 8, 2)))
        assert type(ind) is indexing.OuterIndexer

        vind = Variable(("a",), [0, 1])
        _, ind, _ = v._broadcast_indexes((vind, slice(0, 8, 2)))
        assert type(ind) is indexing.OuterIndexer

        vind = Variable(("y",), [0, 1])
        _, ind, _ = v._broadcast_indexes((vind, 3))
        assert type(ind) is indexing.OuterIndexer

        vind = Variable(("a",), [0, 1])
        _, ind, _ = v._broadcast_indexes((vind, vind))
        assert type(ind) is indexing.VectorizedIndexer

        vind = Variable(("a", "b"), [[0, 2], [1, 3]])
        _, ind, _ = v._broadcast_indexes((vind, 3))
        assert type(ind) is indexing.VectorizedIndexer

    def test_indexer_type(self):
        # GH:issue:1688. Wrong indexer type induces NotImplementedError
        data = np.random.random((10, 11))
        v = Variable(["x", "y"], data)

        def assert_indexer_type(key, object_type):
            _dims, index_tuple, _new_order = v._broadcast_indexes(key)
            assert isinstance(index_tuple, object_type)

        # should return BasicIndexer
        assert_indexer_type((0, 1), BasicIndexer)
        assert_indexer_type((0, slice(None, None)), BasicIndexer)
        assert_indexer_type((Variable([], 3), slice(None, None)), BasicIndexer)
        assert_indexer_type((Variable([], 3), (Variable([], 6))), BasicIndexer)

        # should return OuterIndexer
        assert_indexer_type(([0, 1], 1), OuterIndexer)
        assert_indexer_type(([0, 1], [1, 2]), OuterIndexer)
        assert_indexer_type((Variable(("x"), [0, 1]), 1), OuterIndexer)
        assert_indexer_type((Variable(("x"), [0, 1]), slice(None, None)), OuterIndexer)
        assert_indexer_type(
            (Variable(("x"), [0, 1]), Variable(("y"), [0, 1])), OuterIndexer
        )

        # should return VectorizedIndexer
        assert_indexer_type((Variable(("y"), [0, 1]), [0, 1]), VectorizedIndexer)
        assert_indexer_type(
            (Variable(("z"), [0, 1]), Variable(("z"), [0, 1])), VectorizedIndexer
        )
        assert_indexer_type(
            (
                Variable(("a", "b"), [[0, 1], [1, 2]]),
                Variable(("a", "b"), [[0, 1], [1, 2]]),
            ),
            VectorizedIndexer,
        )

    def test_items(self):
        data = np.random.random((10, 11))
        v = Variable(["x", "y"], data)
        # test slicing
        assert_identical(v, v[:])
        assert_identical(v, v[...])
        assert_identical(Variable(["y"], data[0]), v[0])
        assert_identical(Variable(["x"], data[:, 0]), v[:, 0])
        assert_identical(Variable(["x", "y"], data[:3, :2]), v[:3, :2])
        # test array indexing
        x = Variable(["x"], np.arange(10))
        y = Variable(["y"], np.arange(11))
        assert_identical(v, v[x.values])
        assert_identical(v, v[x])
        assert_identical(v[:3], v[x < 3])
        assert_identical(v[:, 3:], v[:, y >= 3])
        assert_identical(v[:3, 3:], v[x < 3, y >= 3])
        assert_identical(v[:3, :2], v[x[:3], y[:2]])
        assert_identical(v[:3, :2], v[range(3), range(2)])
        # test iteration
        for n, item in enumerate(v):
            assert_identical(Variable(["y"], data[n]), item)
        with pytest.raises(TypeError, match=r"iteration over a 0-d"):
            iter(Variable([], 0))
        # test setting
        v.values[:] = 0
        assert np.all(v.values == 0)
        # test orthogonal setting
        v[range(10), range(11)] = 1
        assert_array_equal(v.values, np.ones((10, 11)))

    def test_getitem_basic(self):
        v = self.cls(["x", "y"], [[0, 1, 2], [3, 4, 5]])

        # int argument
        v_new = v[0]
        assert v_new.dims == ("y",)
        assert_array_equal(v_new, v._data[0])

        # slice argument
        v_new = v[:2]
        assert v_new.dims == ("x", "y")
        assert_array_equal(v_new, v._data[:2])

        # list arguments
        v_new = v[[0]]
        assert v_new.dims == ("x", "y")
        assert_array_equal(v_new, v._data[[0]])  # type: ignore[call-overload]

        v_new = v[[]]
        assert v_new.dims == ("x", "y")
        assert_array_equal(v_new, v._data[[]])  # type: ignore[call-overload]

        # dict arguments
        v_new = v[dict(x=0)]
        assert v_new.dims == ("y",)
        assert_array_equal(v_new, v._data[0])

        v_new = v[dict(x=0, y=slice(None))]
        assert v_new.dims == ("y",)
        assert_array_equal(v_new, v._data[0])

        v_new = v[dict(x=0, y=1)]
        assert v_new.dims == ()
        assert_array_equal(v_new, v._data[0, 1])

        v_new = v[dict(y=1)]
        assert v_new.dims == ("x",)
        assert_array_equal(v_new, v._data[:, 1])

        # tuple argument
        v_new = v[(slice(None), 1)]
        assert v_new.dims == ("x",)
        assert_array_equal(v_new, v._data[:, 1])

        # test that we obtain a modifiable view when taking a 0d slice
        v_new = v[0, 0]
        v_new[...] += 99
        assert_array_equal(v_new, v._data[0, 0])

    def test_getitem_with_mask_2d_input(self):
        v = Variable(("x", "y"), [[0, 1, 2], [3, 4, 5]])
        assert_identical(
            v._getitem_with_mask(([-1, 0], [1, -1])),
            Variable(("x", "y"), [[np.nan, np.nan], [1, np.nan]]),
        )
        assert_identical(v._getitem_with_mask((slice(2), [0, 1, 2])), v)

    def test_isel(self):
        v = Variable(["time", "x"], self.d)
        assert_identical(v.isel(time=slice(None)), v)
        assert_identical(v.isel(time=0), v[0])
        assert_identical(v.isel(time=slice(0, 3)), v[:3])
        assert_identical(v.isel(x=0), v[:, 0])
        assert_identical(v.isel(x=[0, 2]), v[:, [0, 2]])
        assert_identical(v.isel(time=[]), v[[]])
        with pytest.raises(
            ValueError,
            match=r"Dimensions {'not_a_dim'} do not exist. Expected one or more of "
            r"\('time', 'x'\)",
        ):
            v.isel(not_a_dim=0)
        with pytest.warns(
            UserWarning,
            match=r"Dimensions {'not_a_dim'} do not exist. Expected one or more of "
            r"\('time', 'x'\)",
        ):
            v.isel(not_a_dim=0, missing_dims="warn")
        assert_identical(v, v.isel(not_a_dim=0, missing_dims="ignore"))

    def test_index_0d_numpy_string(self):
        # regression test to verify our work around for indexing 0d strings
        v = Variable([], np.bytes_("asdf"))
        assert_identical(v[()], v)

        v = Variable([], np.str_("asdf"))
        assert_identical(v[()], v)

    def test_indexing_0d_unicode(self):
        # regression test for GH568
        actual = Variable(("x"), ["tmax"])[0][()]
        expected = Variable((), "tmax")
        assert_identical(actual, expected)

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0])
    def test_shift(self, fill_value):
        v = Variable("x", [1, 2, 3, 4, 5])

        assert_identical(v, v.shift(x=0))
        assert v is not v.shift(x=0)

        expected = Variable("x", [np.nan, np.nan, 1, 2, 3])
        assert_identical(expected, v.shift(x=2))

        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_exp = np.nan
        else:
            fill_value_exp = fill_value

        expected = Variable("x", [fill_value_exp, 1, 2, 3, 4])
        assert_identical(expected, v.shift(x=1, fill_value=fill_value))

        expected = Variable("x", [2, 3, 4, 5, fill_value_exp])
        assert_identical(expected, v.shift(x=-1, fill_value=fill_value))

        expected = Variable("x", [fill_value_exp] * 5)
        assert_identical(expected, v.shift(x=5, fill_value=fill_value))
        assert_identical(expected, v.shift(x=6, fill_value=fill_value))

        with pytest.raises(ValueError, match=r"dimension"):
            v.shift(z=0)

        v = Variable("x", [1, 2, 3, 4, 5], {"foo": "bar"})
        assert_identical(v, v.shift(x=0))

        expected = Variable("x", [fill_value_exp, 1, 2, 3, 4], {"foo": "bar"})
        assert_identical(expected, v.shift(x=1, fill_value=fill_value))

    def test_shift2d(self):
        v = Variable(("x", "y"), [[1, 2], [3, 4]])
        expected = Variable(("x", "y"), [[np.nan, np.nan], [np.nan, 1]])
        assert_identical(expected, v.shift(x=1, y=1))

    def test_roll(self):
        v = Variable("x", [1, 2, 3, 4, 5])

        assert_identical(v, v.roll(x=0))
        assert v is not v.roll(x=0)

        expected = Variable("x", [5, 1, 2, 3, 4])
        assert_identical(expected, v.roll(x=1))
        assert_identical(expected, v.roll(x=-4))
        assert_identical(expected, v.roll(x=6))

        expected = Variable("x", [4, 5, 1, 2, 3])
        assert_identical(expected, v.roll(x=2))
        assert_identical(expected, v.roll(x=-3))

        with pytest.raises(ValueError, match=r"dimension"):
            v.roll(z=0)

    def test_roll_consistency(self):
        v = Variable(("x", "y"), np.random.randn(5, 6))

        for axis, dim in [(0, "x"), (1, "y")]:
            for shift in [-3, 0, 1, 7, 11]:
                expected = np.roll(v.values, shift, axis=axis)
                actual = v.roll(**{dim: shift}).values
                assert_array_equal(expected, actual)

    def test_transpose(self):
        v = Variable(["time", "x"], self.d)
        v2 = Variable(["x", "time"], self.d.T)
        assert_identical(v, v2.transpose())
        assert_identical(v.transpose(), v.T)
        x = np.random.randn(2, 3, 4, 5)
        w = Variable(["a", "b", "c", "d"], x)
        w2 = Variable(["d", "b", "c", "a"], np.einsum("abcd->dbca", x))
        assert w2.shape == (5, 3, 4, 2)
        assert_identical(w2, w.transpose("d", "b", "c", "a"))
        assert_identical(w2, w.transpose("d", ..., "a"))
        assert_identical(w2, w.transpose("d", "b", "c", ...))
        assert_identical(w2, w.transpose(..., "b", "c", "a"))
        assert_identical(w, w2.transpose("a", "b", "c", "d"))
        w3 = Variable(["b", "c", "d", "a"], np.einsum("abcd->bcda", x))
        assert_identical(w, w3.transpose("a", "b", "c", "d"))

        # test missing dimension, raise error
        with pytest.raises(ValueError):
            v.transpose(..., "not_a_dim")

        # test missing dimension, ignore error
        actual = v.transpose(..., "not_a_dim", missing_dims="ignore")
        expected_ell = v.transpose(...)
        assert_identical(expected_ell, actual)

        # test missing dimension, raise warning
        with pytest.warns(UserWarning):
            v.transpose(..., "not_a_dim", missing_dims="warn")
            assert_identical(expected_ell, actual)

    def test_transpose_0d(self):
        for value in [
            3.5,
            ("a", 1),
            np.datetime64("2000-01-01"),
            np.timedelta64(1, "h"),
            None,
            object(),
        ]:
            variable = Variable([], value)
            actual = variable.transpose()
            assert_identical(actual, variable)

    def test_pandas_categorical_dtype(self):
        data = pd.Categorical(np.arange(10, dtype="int64"))
        v = self.cls("x", data)
        print(v)  # should not error
        assert isinstance(v.dtype, pd.CategoricalDtype)

    def test_squeeze(self):
        v = Variable(["x", "y"], [[1]])
        assert_identical(Variable([], 1), v.squeeze())
        assert_identical(Variable(["y"], [1]), v.squeeze("x"))
        assert_identical(Variable(["y"], [1]), v.squeeze(["x"]))
        assert_identical(Variable(["x"], [1]), v.squeeze("y"))
        assert_identical(Variable([], 1), v.squeeze(["x", "y"]))

        v = Variable(["x", "y"], [[1, 2]])
        assert_identical(Variable(["y"], [1, 2]), v.squeeze())
        assert_identical(Variable(["y"], [1, 2]), v.squeeze("x"))
        with pytest.raises(ValueError, match=r"cannot select a dimension"):
            v.squeeze("y")

    def test_get_axis_num(self) -> None:
        v = Variable(["x", "y", "z"], np.random.randn(2, 3, 4))
        assert v.get_axis_num("x") == 0
        assert v.get_axis_num(["x"]) == (0,)
        assert v.get_axis_num(["x", "y"]) == (0, 1)
        assert v.get_axis_num(["z", "y", "x"]) == (2, 1, 0)
        with pytest.raises(ValueError, match=r"not found in array dim"):
            v.get_axis_num("foobar")
        # Test the type annotations: mypy will complain if the inferred
        # type is wrong
        v.get_axis_num("x") + 0
        v.get_axis_num(["x"]) + ()
        v.get_axis_num(("x", "y")) + ()

    def test_set_dims(self):
        v = Variable(["x"], [0, 1])
        actual = v.set_dims(["x", "y"])
        expected = Variable(["x", "y"], [[0], [1]])
        assert_identical(actual, expected)

        actual = v.set_dims(["y", "x"])
        assert_identical(actual, expected.T)

        actual = v.set_dims({"x": 2, "y": 2})
        expected = Variable(["x", "y"], [[0, 0], [1, 1]])
        assert_identical(actual, expected)

        v = Variable(["foo"], [0, 1])
        actual = v.set_dims("foo")
        expected = v
        assert_identical(actual, expected)

        with pytest.raises(ValueError, match=r"must be a superset"):
            v.set_dims(["z"])

    def test_set_dims_object_dtype(self):
        v = Variable([], ("a", 1))
        actual = v.set_dims(("x",), (3,))
        exp_values = np.empty((3,), dtype=object)
        for i in range(3):
            exp_values[i] = ("a", 1)
        expected = Variable(["x"], exp_values)
        assert_identical(actual, expected)

    def test_set_dims_without_broadcast(self):
        class ArrayWithoutBroadcastTo(NDArrayMixin, indexing.ExplicitlyIndexed):
            def __init__(self, array):
                self.array = array

            # Broadcasting with __getitem__ is "easier" to implement
            # especially for dims of 1
            def __getitem__(self, key):
                return self.array[key]

            def __array_function__(self, *args, **kwargs):
                raise NotImplementedError(
                    "Not we don't want to use broadcast_to here "
                    "https://github.com/pydata/xarray/issues/9462"
                )

        arr = ArrayWithoutBroadcastTo(np.zeros((3, 4)))
        # We should be able to add a new axis without broadcasting
        assert arr[np.newaxis, :, :].shape == (1, 3, 4)
        with pytest.raises(NotImplementedError):
            np.broadcast_to(arr, (1, 3, 4))

        v = Variable(["x", "y"], arr)
        v_expanded = v.set_dims(["z", "x", "y"])
        assert v_expanded.dims == ("z", "x", "y")
        assert v_expanded.shape == (1, 3, 4)

        v_expanded = v.set_dims(["x", "z", "y"])
        assert v_expanded.dims == ("x", "z", "y")
        assert v_expanded.shape == (3, 1, 4)

        v_expanded = v.set_dims(["x", "y", "z"])
        assert v_expanded.dims == ("x", "y", "z")
        assert v_expanded.shape == (3, 4, 1)

        # Explicitly asking for a shape of 1 triggers a different
        # codepath in set_dims
        # https://github.com/pydata/xarray/issues/9462
        v_expanded = v.set_dims(["z", "x", "y"], shape=(1, 3, 4))
        assert v_expanded.dims == ("z", "x", "y")
        assert v_expanded.shape == (1, 3, 4)

        v_expanded = v.set_dims(["x", "z", "y"], shape=(3, 1, 4))
        assert v_expanded.dims == ("x", "z", "y")
        assert v_expanded.shape == (3, 1, 4)

        v_expanded = v.set_dims(["x", "y", "z"], shape=(3, 4, 1))
        assert v_expanded.dims == ("x", "y", "z")
        assert v_expanded.shape == (3, 4, 1)

        v_expanded = v.set_dims({"z": 1, "x": 3, "y": 4})
        assert v_expanded.dims == ("z", "x", "y")
        assert v_expanded.shape == (1, 3, 4)

        v_expanded = v.set_dims({"x": 3, "z": 1, "y": 4})
        assert v_expanded.dims == ("x", "z", "y")
        assert v_expanded.shape == (3, 1, 4)

        v_expanded = v.set_dims({"x": 3, "y": 4, "z": 1})
        assert v_expanded.dims == ("x", "y", "z")
        assert v_expanded.shape == (3, 4, 1)

        with pytest.raises(NotImplementedError):
            v.set_dims({"z": 2, "x": 3, "y": 4})

        with pytest.raises(NotImplementedError):
            v.set_dims(["z", "x", "y"], shape=(2, 3, 4))

    def test_stack(self):
        v = Variable(["x", "y"], [[0, 1], [2, 3]], {"foo": "bar"})
        actual = v.stack(z=("x", "y"))
        expected = Variable("z", [0, 1, 2, 3], v.attrs)
        assert_identical(actual, expected)

        actual = v.stack(z=("x",))
        expected = Variable(("y", "z"), v.data.T, v.attrs)
        assert_identical(actual, expected)

        actual = v.stack(z=())
        assert_identical(actual, v)

        actual = v.stack(X=("x",), Y=("y",)).transpose("X", "Y")
        expected = Variable(("X", "Y"), v.data, v.attrs)
        assert_identical(actual, expected)

    def test_stack_errors(self):
        v = Variable(["x", "y"], [[0, 1], [2, 3]], {"foo": "bar"})

        with pytest.raises(ValueError, match=r"invalid existing dim"):
            v.stack(z=("x1",))
        with pytest.raises(ValueError, match=r"cannot create a new dim"):
            v.stack(x=("x",))

    def test_unstack(self):
        v = Variable("z", [0, 1, 2, 3], {"foo": "bar"})
        actual = v.unstack(z={"x": 2, "y": 2})
        expected = Variable(("x", "y"), [[0, 1], [2, 3]], v.attrs)
        assert_identical(actual, expected)

        actual = v.unstack(z={"x": 4, "y": 1})
        expected = Variable(("x", "y"), [[0], [1], [2], [3]], v.attrs)
        assert_identical(actual, expected)

        actual = v.unstack(z={"x": 4})
        expected = Variable("x", [0, 1, 2, 3], v.attrs)
        assert_identical(actual, expected)

    def test_unstack_errors(self):
        v = Variable("z", [0, 1, 2, 3])
        with pytest.raises(ValueError, match=r"invalid existing dim"):
            v.unstack(foo={"x": 4})
        with pytest.raises(ValueError, match=r"cannot create a new dim"):
            v.stack(z=("z",))
        with pytest.raises(ValueError, match=r"the product of the new dim"):
            v.unstack(z={"x": 5})

    def test_unstack_2d(self):
        v = Variable(["x", "y"], [[0, 1], [2, 3]])
        actual = v.unstack(y={"z": 2})
        expected = Variable(["x", "z"], v.data)
        assert_identical(actual, expected)

        actual = v.unstack(x={"z": 2})
        expected = Variable(["y", "z"], v.data.T)
        assert_identical(actual, expected)

    def test_stack_unstack_consistency(self):
        v = Variable(["x", "y"], [[0, 1], [2, 3]])
        actual = v.stack(z=("x", "y")).unstack(z={"x": 2, "y": 2})
        assert_identical(actual, v)

    @pytest.mark.filterwarnings("error::RuntimeWarning")
    def test_unstack_without_missing(self):
        v = Variable(["z"], [0, 1, 2, 3])
        expected = Variable(["x", "y"], [[0, 1], [2, 3]])

        actual = v.unstack(z={"x": 2, "y": 2})

        assert_identical(actual, expected)

    def test_broadcasting_math(self):
        x = np.random.randn(2, 3)
        v = Variable(["a", "b"], x)
        # 1d to 2d broadcasting
        assert_identical(v * v, Variable(["a", "b"], np.einsum("ab,ab->ab", x, x)))
        assert_identical(v * v[0], Variable(["a", "b"], np.einsum("ab,b->ab", x, x[0])))
        assert_identical(v[0] * v, Variable(["b", "a"], np.einsum("b,ab->ba", x[0], x)))
        assert_identical(
            v[0] * v[:, 0], Variable(["b", "a"], np.einsum("b,a->ba", x[0], x[:, 0]))
        )
        # higher dim broadcasting
        y = np.random.randn(3, 4, 5)
        w = Variable(["b", "c", "d"], y)
        assert_identical(
            v * w, Variable(["a", "b", "c", "d"], np.einsum("ab,bcd->abcd", x, y))
        )
        assert_identical(
            w * v, Variable(["b", "c", "d", "a"], np.einsum("bcd,ab->bcda", y, x))
        )
        assert_identical(
            v * w[0], Variable(["a", "b", "c", "d"], np.einsum("ab,cd->abcd", x, y[0]))
        )

    @pytest.mark.filterwarnings("ignore:Duplicate dimension names")
    def test_broadcasting_failures(self):
        a = Variable(["x"], np.arange(10))
        b = Variable(["x"], np.arange(5))
        c = Variable(["x", "x"], np.arange(100).reshape(10, 10))
        with pytest.raises(ValueError, match=r"mismatched lengths"):
            a + b
        with pytest.raises(ValueError, match=r"duplicate dimensions"):
            a + c

    def test_inplace_math(self):
        x = np.arange(5)
        v = Variable(["x"], x)
        v2 = v
        v2 += 1
        assert v is v2
        # since we provided an ndarray for data, it is also modified in-place
        assert source_ndarray(v.values) is x
        assert_array_equal(v.values, np.arange(5) + 1)

        with pytest.raises(ValueError, match=r"dimensions cannot change"):
            v += Variable("y", np.arange(5))

    def test_inplace_math_error(self):
        x = np.arange(5)
        v = IndexVariable(["x"], x)
        with pytest.raises(
            TypeError, match=r"Values of an IndexVariable are immutable"
        ):
            v += 1

    def test_reduce(self):
        v = Variable(["x", "y"], self.d, {"ignored": "attributes"})
        # Reduce keeps attrs by default
        expected = Variable(["y"], self.d.std(axis=0), {"ignored": "attributes"})
        assert_identical(v.reduce(np.std, "x"), expected)
        assert_identical(v.reduce(np.std, axis=0), v.reduce(np.std, dim="x"))
        assert_identical(
            v.reduce(np.std, ["y", "x"]),
            Variable([], self.d.std(axis=(0, 1)), {"ignored": "attributes"}),
        )
        assert_identical(
            v.reduce(np.std), Variable([], self.d.std(), {"ignored": "attributes"})
        )
        # Chained reductions both keep attrs
        expected_chained = Variable(
            [], self.d.mean(axis=0).std(), {"ignored": "attributes"}
        )
        assert_identical(
            v.reduce(np.mean, "x").reduce(np.std, "y"),
            expected_chained,
        )
        assert_allclose(v.mean("x"), v.reduce(np.mean, "x"))

        with pytest.raises(ValueError, match=r"cannot supply both"):
            v.mean(dim="x", axis=0)

    @requires_bottleneck
    @pytest.mark.parametrize("compute_backend", ["bottleneck"], indirect=True)
    def test_reduce_use_bottleneck(self, monkeypatch, compute_backend):
        def raise_if_called(*args, **kwargs):
            raise RuntimeError("should not have been called")

        import bottleneck as bn

        monkeypatch.setattr(bn, "nanmin", raise_if_called)

        v = Variable("x", [0.0, np.nan, 1.0])
        with pytest.raises(RuntimeError, match="should not have been called"):
            with set_options(use_bottleneck=True):
                v.min()

        with set_options(use_bottleneck=False):
            v.min()

    @pytest.mark.parametrize("skipna", [True, False, None])
    @pytest.mark.parametrize("q", [0.25, [0.50], [0.25, 0.75]])
    @pytest.mark.parametrize(
        "axis, dim",
        zip([None, 0, [0], [0, 1]], [None, "x", ["x"], ["x", "y"]], strict=True),
    )
    def test_quantile(self, q, axis, dim, skipna):
        d = self.d.copy()
        d[0, 0] = np.nan

        v = Variable(["x", "y"], d)
        actual = v.quantile(q, dim=dim, skipna=skipna)
        _percentile_func = np.nanpercentile if skipna in (True, None) else np.percentile
        expected = _percentile_func(d, np.array(q) * 100, axis=axis)
        np.testing.assert_allclose(actual.values, expected)

    @requires_dask
    @pytest.mark.parametrize("q", [0.25, [0.50], [0.25, 0.75]])
    @pytest.mark.parametrize("axis, dim", [[1, "y"], [[1], ["y"]]])
    def test_quantile_dask(self, q, axis, dim):
        v = Variable(["x", "y"], self.d).chunk({"x": 2})
        actual = v.quantile(q, dim=dim)
        assert isinstance(actual.data, dask_array_type)
        expected = np.nanpercentile(self.d, np.array(q) * 100, axis=axis)
        np.testing.assert_allclose(actual.values, expected)

    @pytest.mark.parametrize("method", ["midpoint", "lower"])
    @pytest.mark.parametrize(
        "use_dask", [pytest.param(True, marks=requires_dask), False]
    )
    def test_quantile_method(self, method, use_dask) -> None:
        v = Variable(["x", "y"], self.d)
        if use_dask:
            v = v.chunk({"x": 2})

        q = np.array([0.25, 0.5, 0.75])
        actual = v.quantile(q, dim="y", method=method)

        expected = np.nanquantile(self.d, q, axis=1, method=method)

        if use_dask:
            assert isinstance(actual.data, dask_array_type)

        np.testing.assert_allclose(actual.values, expected)

    @pytest.mark.filterwarnings(
        "default:The `interpolation` argument to quantile was renamed to `method`:FutureWarning"
    )
    @pytest.mark.parametrize("method", ["midpoint", "lower"])
    def test_quantile_interpolation_deprecation(self, method) -> None:
        v = Variable(["x", "y"], self.d)
        q = np.array([0.25, 0.5, 0.75])

        with pytest.warns(
            FutureWarning,
            match="`interpolation` argument to quantile was renamed to `method`",
        ):
            actual = v.quantile(q, dim="y", interpolation=method)

        expected = v.quantile(q, dim="y", method=method)

        np.testing.assert_allclose(actual.values, expected.values)

        with warnings.catch_warnings(record=True):
            with pytest.raises(TypeError, match="interpolation and method keywords"):
                v.quantile(q, dim="y", interpolation=method, method=method)

    @requires_dask
    def test_quantile_chunked_dim_error(self):
        v = Variable(["x", "y"], self.d).chunk({"x": 2})

        if has_dask_ge_2024_11_0:
            # Dask rechunks
            np.testing.assert_allclose(
                v.compute().quantile(0.5, dim="x"), v.quantile(0.5, dim="x")
            )

        else:
            # this checks for ValueError in dask.array.apply_gufunc
            with pytest.raises(ValueError, match=r"consists of multiple chunks"):
                v.quantile(0.5, dim="x")

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize("q", [-0.1, 1.1, [2], [0.25, 2]])
    def test_quantile_out_of_bounds(self, q, compute_backend):
        v = Variable(["x", "y"], self.d)

        # escape special characters
        with pytest.raises(
            ValueError,
            match=r"(Q|q)uantiles must be in the range \[0, 1\]",
        ):
            v.quantile(q, dim="x")

    @requires_dask
    @requires_bottleneck
    def test_rank_dask(self):
        # Instead of a single test here, we could parameterize the other tests for both
        # arrays. But this is sufficient.
        v = Variable(
            ["x", "y"], [[30.0, 1.0, np.nan, 20.0, 4.0], [30.0, 1.0, np.nan, 20.0, 4.0]]
        ).chunk(x=1)
        expected = Variable(
            ["x", "y"], [[4.0, 1.0, np.nan, 3.0, 2.0], [4.0, 1.0, np.nan, 3.0, 2.0]]
        )
        assert_equal(v.rank("y").compute(), expected)

        with pytest.raises(
            ValueError, match=r" with dask='parallelized' consists of multiple chunks"
        ):
            v.rank("x")

    def test_rank_use_bottleneck(self):
        v = Variable(["x"], [3.0, 1.0, np.nan, 2.0, 4.0])
        with set_options(use_bottleneck=False):
            with pytest.raises(RuntimeError):
                v.rank("x")

    @requires_bottleneck
    def test_rank(self):
        import bottleneck as bn

        # floats
        v = Variable(["x", "y"], [[3, 4, np.nan, 1]])
        expect_0 = bn.nanrankdata(v.data, axis=0)
        expect_1 = bn.nanrankdata(v.data, axis=1)
        np.testing.assert_allclose(v.rank("x").values, expect_0)
        np.testing.assert_allclose(v.rank("y").values, expect_1)
        # int
        v = Variable(["x"], [3, 2, 1])
        expect = bn.rankdata(v.data, axis=0)
        np.testing.assert_allclose(v.rank("x").values, expect)
        # str
        v = Variable(["x"], ["c", "b", "a"])
        expect = bn.rankdata(v.data, axis=0)
        np.testing.assert_allclose(v.rank("x").values, expect)
        # pct
        v = Variable(["x"], [3.0, 1.0, np.nan, 2.0, 4.0])
        v_expect = Variable(["x"], [0.75, 0.25, np.nan, 0.5, 1.0])
        assert_equal(v.rank("x", pct=True), v_expect)
        # invalid dim
        with pytest.raises(ValueError):
            # apply_ufunc error message isn't great here — `ValueError: tuple.index(x): x not in tuple`
            v.rank("y")

    def test_big_endian_reduce(self):
        # regression test for GH489
        data = np.ones(5, dtype=">f4")
        v = Variable(["x"], data)
        expected = Variable([], 5)
        assert_identical(expected, v.sum())

    def test_reduce_funcs(self):
        v = Variable("x", np.array([1, np.nan, 2, 3]))
        assert_identical(v.mean(), Variable([], 2))
        assert_identical(v.mean(skipna=True), Variable([], 2))
        assert_identical(v.mean(skipna=False), Variable([], np.nan))
        assert_identical(np.mean(v), Variable([], 2))

        assert_identical(v.prod(), Variable([], 6))
        assert_identical(v.cumsum(axis=0), Variable("x", np.array([1, 1, 3, 6])))
        assert_identical(v.cumprod(axis=0), Variable("x", np.array([1, 1, 2, 6])))
        assert_identical(v.var(), Variable([], 2.0 / 3))
        assert_identical(v.median(), Variable([], 2))

        v = Variable("x", [True, False, False])
        assert_identical(v.any(), Variable([], True))
        assert_identical(v.all(dim="x"), Variable([], False))

        v = Variable("t", pd.date_range("2000-01-01", periods=3))
        assert v.argmax(skipna=True, dim="t") == 2

        assert_identical(v.max(), Variable([], pd.Timestamp("2000-01-03")))

    def test_reduce_keepdims(self):
        v = Variable(["x", "y"], self.d)

        with set_options(use_numbagg=False):
            assert_identical(
                v.mean(keepdims=True), Variable(v.dims, np.mean(self.d, keepdims=True))
            )
            assert_identical(
                v.mean(dim="x", keepdims=True),
                Variable(v.dims, np.mean(self.d, axis=0, keepdims=True)),
            )
            assert_identical(
                v.mean(dim="y", keepdims=True),
                Variable(v.dims, np.mean(self.d, axis=1, keepdims=True)),
            )
            assert_identical(
                v.mean(dim=["y", "x"], keepdims=True),
                Variable(v.dims, np.mean(self.d, axis=(1, 0), keepdims=True)),
            )

            v = Variable([], 1.0)
            assert_identical(
                v.mean(keepdims=True), Variable([], np.mean(v.data, keepdims=True))
            )

    @requires_dask
    def test_reduce_keepdims_dask(self):
        import dask.array

        v = Variable(["x", "y"], self.d).chunk()

        actual = v.mean(keepdims=True)
        assert isinstance(actual.data, dask.array.Array)

        expected = Variable(v.dims, np.mean(self.d, keepdims=True))
        assert_identical(actual, expected)

        actual = v.mean(dim="y", keepdims=True)
        assert isinstance(actual.data, dask.array.Array)

        expected = Variable(v.dims, np.mean(self.d, axis=1, keepdims=True))
        assert_identical(actual, expected)

    def test_reduce_keep_attrs(self):
        _attrs = {"units": "test", "long_name": "testing"}

        v = Variable(["x", "y"], self.d, _attrs)

        # Test default behavior (keeps attrs for reduction operations)
        vm = v.mean()
        assert len(vm.attrs) == len(_attrs)
        assert vm.attrs == _attrs

        # Test explicitly keeping attrs
        vm = v.mean(keep_attrs=True)
        assert len(vm.attrs) == len(_attrs)
        assert vm.attrs == _attrs

        # Test explicitly dropping attrs
        vm = v.mean(keep_attrs=False)
        assert len(vm.attrs) == 0
        assert vm.attrs == {}

    def test_binary_ops_keep_attrs(self):
        _attrs = {"units": "test", "long_name": "testing"}
        a = Variable(["x", "y"], np.random.randn(3, 3), _attrs)
        b = Variable(["x", "y"], np.random.randn(3, 3), _attrs)
        # Test kept attrs (now default)
        d = a - b  # just one operation
        assert d.attrs == _attrs
        # Test dropped attrs
        with set_options(keep_attrs=False):
            d = a - b
        assert d.attrs == {}

    def test_binary_ops_attrs_drop_conflicts(self):
        # Test that binary operations combine attrs with drop_conflicts behavior
        attrs_a = {"units": "meters", "long_name": "distance", "source": "sensor_a"}
        attrs_b = {"units": "feet", "resolution": "high", "source": "sensor_b"}
        a = Variable(["x"], [1, 2, 3], attrs_a)
        b = Variable(["x"], [4, 5, 6], attrs_b)

        # With keep_attrs=True (default), should combine attrs dropping conflicts
        result = a + b
        # "units" and "source" conflict, so they're dropped
        # "long_name" only in a, "resolution" only in b, so they're kept
        assert result.attrs == {"long_name": "distance", "resolution": "high"}

        # Test with identical values for some attrs
        attrs_c = {"units": "meters", "type": "data", "source": "sensor_c"}
        c = Variable(["x"], [7, 8, 9], attrs_c)
        result2 = a + c
        # "units" has same value, so kept; "source" conflicts, so dropped
        # "long_name" from a, "type" from c
        assert result2.attrs == {
            "units": "meters",
            "long_name": "distance",
            "type": "data",
        }

        # With keep_attrs=False, attrs should be empty
        with set_options(keep_attrs=False):
            result3 = a + b
            assert result3.attrs == {}

    def test_count(self):
        expected = Variable([], 3)
        actual = Variable(["x"], [1, 2, 3, np.nan]).count()
        assert_identical(expected, actual)

        v = Variable(["x"], np.array(["1", "2", "3", np.nan], dtype=object))
        actual = v.count()
        assert_identical(expected, actual)

        actual = Variable(["x"], [True, False, True]).count()
        assert_identical(expected, actual)
        assert actual.dtype == int

        expected = Variable(["x"], [2, 3])
        actual = Variable(["x", "y"], [[1, 0, np.nan], [1, 1, 1]]).count("y")
        assert_identical(expected, actual)

    def test_setitem(self):
        v = Variable(["x", "y"], [[0, 3, 2], [3, 4, 5]])
        v[0, 1] = 1
        assert v[0, 1] == 1

        v = Variable(["x", "y"], [[0, 3, 2], [3, 4, 5]])
        v[dict(x=[0, 1])] = 1
        assert_array_equal(v[[0, 1]], np.ones_like(v[[0, 1]]))

        # boolean indexing
        v = Variable(["x", "y"], [[0, 3, 2], [3, 4, 5]])
        v[dict(x=[True, False])] = 1

        assert_array_equal(v[0], np.ones_like(v[0]))
        v = Variable(["x", "y"], [[0, 3, 2], [3, 4, 5]])
        v[dict(x=[True, False], y=[False, True, False])] = 1
        assert v[0, 1] == 1

    def test_setitem_fancy(self):
        # assignment which should work as np.ndarray does
        def assert_assigned_2d(array, key_x, key_y, values):
            expected = array.copy()
            expected[key_x, key_y] = values
            v = Variable(["x", "y"], array)
            v[dict(x=key_x, y=key_y)] = values
            assert_array_equal(expected, v)

        # 1d vectorized indexing
        assert_assigned_2d(
            np.random.randn(4, 3),
            key_x=Variable(["a"], [0, 1]),
            key_y=Variable(["a"], [0, 1]),
            values=0,
        )
        assert_assigned_2d(
            np.random.randn(4, 3),
            key_x=Variable(["a"], [0, 1]),
            key_y=Variable(["a"], [0, 1]),
            values=Variable((), 0),
        )
        assert_assigned_2d(
            np.random.randn(4, 3),
            key_x=Variable(["a"], [0, 1]),
            key_y=Variable(["a"], [0, 1]),
            values=Variable(("a"), [3, 2]),
        )
        assert_assigned_2d(
            np.random.randn(4, 3),
            key_x=slice(None),
            key_y=Variable(["a"], [0, 1]),
            values=Variable(("a"), [3, 2]),
        )

        # 2d-vectorized indexing
        assert_assigned_2d(
            np.random.randn(4, 3),
            key_x=Variable(["a", "b"], [[0, 1]]),
            key_y=Variable(["a", "b"], [[1, 0]]),
            values=0,
        )
        assert_assigned_2d(
            np.random.randn(4, 3),
            key_x=Variable(["a", "b"], [[0, 1]]),
            key_y=Variable(["a", "b"], [[1, 0]]),
            values=[0],
        )
        assert_assigned_2d(
            np.random.randn(5, 4),
            key_x=Variable(["a", "b"], [[0, 1], [2, 3]]),
            key_y=Variable(["a", "b"], [[1, 0], [3, 3]]),
            values=[2, 3],
        )

        # vindex with slice
        v = Variable(["x", "y", "z"], np.ones((4, 3, 2)))
        ind = Variable(["a"], [0, 1])
        v[dict(x=ind, z=ind)] = 0
        expected = Variable(["x", "y", "z"], np.ones((4, 3, 2)))
        expected[0, :, 0] = 0
        expected[1, :, 1] = 0
        assert_identical(expected, v)

        # dimension broadcast
        v = Variable(["x", "y"], np.ones((3, 2)))
        ind = Variable(["a", "b"], [[0, 1]])
        v[ind, :] = 0
        expected = Variable(["x", "y"], [[0, 0], [0, 0], [1, 1]])
        assert_identical(expected, v)

        with pytest.raises(ValueError, match=r"shape mismatch"):
            v[ind, ind] = np.zeros((1, 2, 1))

        v = Variable(["x", "y"], [[0, 3, 2], [3, 4, 5]])
        ind = Variable(["a"], [0, 1])
        v[dict(x=ind)] = Variable(["a", "y"], np.ones((2, 3), dtype=int) * 10)
        assert_array_equal(v[0], np.ones_like(v[0]) * 10)
        assert_array_equal(v[1], np.ones_like(v[1]) * 10)
        assert v.dims == ("x", "y")  # dimension should not change

        # increment
        v = Variable(["x", "y"], np.arange(6).reshape(3, 2))
        ind = Variable(["a"], [0, 1])
        v[dict(x=ind)] += 1
        expected = Variable(["x", "y"], [[1, 2], [3, 4], [4, 5]])
        assert_identical(v, expected)

        ind = Variable(["a"], [0, 0])
        v[dict(x=ind)] += 1
        expected = Variable(["x", "y"], [[2, 3], [3, 4], [4, 5]])
        assert_identical(v, expected)

    def test_coarsen(self):
        v = self.cls(["x"], [0, 1, 2, 3, 4])
        actual = v.coarsen({"x": 2}, boundary="pad", func="mean")
        expected = self.cls(["x"], [0.5, 2.5, 4])
        assert_identical(actual, expected)

        actual = v.coarsen({"x": 2}, func="mean", boundary="pad", side="right")
        expected = self.cls(["x"], [0, 1.5, 3.5])
        assert_identical(actual, expected)

        actual = v.coarsen({"x": 2}, func=np.mean, side="right", boundary="trim")
        expected = self.cls(["x"], [1.5, 3.5])
        assert_identical(actual, expected)

        # working test
        v = self.cls(["x", "y", "z"], np.arange(40 * 30 * 2).reshape(40, 30, 2))
        for windows, func, side, boundary in [
            ({"x": 2}, np.mean, "left", "trim"),
            ({"x": 2}, np.median, {"x": "left"}, "pad"),
            ({"x": 2, "y": 3}, np.max, "left", {"x": "pad", "y": "trim"}),
        ]:
            v.coarsen(windows, func, boundary, side)

    def test_coarsen_2d(self):
        # 2d-mean should be the same with the successive 1d-mean
        v = self.cls(["x", "y"], np.arange(6 * 12).reshape(6, 12))
        actual = v.coarsen({"x": 3, "y": 4}, func="mean")
        expected = v.coarsen({"x": 3}, func="mean").coarsen({"y": 4}, func="mean")
        assert_equal(actual, expected)

        v = self.cls(["x", "y"], np.arange(7 * 12).reshape(7, 12))
        actual = v.coarsen({"x": 3, "y": 4}, func="mean", boundary="trim")
        expected = v.coarsen({"x": 3}, func="mean", boundary="trim").coarsen(
            {"y": 4}, func="mean", boundary="trim"
        )
        assert_equal(actual, expected)

        # if there is nan, the two should be different
        v = self.cls(["x", "y"], 1.0 * np.arange(6 * 12).reshape(6, 12))
        v[2, 4] = np.nan
        v[3, 5] = np.nan
        actual = v.coarsen({"x": 3, "y": 4}, func="mean", boundary="trim")
        expected = (
            v.coarsen({"x": 3}, func="sum", boundary="trim").coarsen(
                {"y": 4}, func="sum", boundary="trim"
            )
            / 12
        )
        assert not actual.equals(expected)
        # adjusting the nan count
        expected[0, 1] *= 12 / 11
        expected[1, 1] *= 12 / 11
        assert_allclose(actual, expected)

        v = self.cls(("x", "y"), np.arange(4 * 4, dtype=np.float32).reshape(4, 4))
        actual = v.coarsen(dict(x=2, y=2), func="count", boundary="exact")
        expected = self.cls(("x", "y"), 4 * np.ones((2, 2)))
        assert_equal(actual, expected)

        v[0, 0] = np.nan
        v[-1, -1] = np.nan
        expected[0, 0] = 3
        expected[-1, -1] = 3
        actual = v.coarsen(dict(x=2, y=2), func="count", boundary="exact")
        assert_equal(actual, expected)

        actual = v.coarsen(dict(x=2, y=2), func="sum", boundary="exact", skipna=False)
        expected = self.cls(("x", "y"), [[np.nan, 18], [42, np.nan]])
        assert_equal(actual, expected)

        actual = v.coarsen(dict(x=2, y=2), func="sum", boundary="exact", skipna=True)
        expected = self.cls(("x", "y"), [[10, 18], [42, 35]])
        assert_equal(actual, expected)

    # perhaps @pytest.mark.parametrize("operation", [f for f in duck_array_ops])
    def test_coarsen_keep_attrs(self, operation="mean"):
        _attrs = {"units": "test", "long_name": "testing"}

        test_func = getattr(duck_array_ops, operation, None)

        # Test dropped attrs
        with set_options(keep_attrs=False):
            new = Variable(["coord"], np.linspace(1, 10, 100), attrs=_attrs).coarsen(
                windows={"coord": 1}, func=test_func, boundary="exact", side="left"
            )
        assert new.attrs == {}

        # Test kept attrs
        with set_options(keep_attrs=True):
            new = Variable(["coord"], np.linspace(1, 10, 100), attrs=_attrs).coarsen(
                windows={"coord": 1},
                func=test_func,
                boundary="exact",
                side="left",
            )
        assert new.attrs == _attrs


@requires_dask
class TestVariableWithDask(VariableSubclassobjects):
    def cls(self, *args, **kwargs) -> Variable:
        return Variable(*args, **kwargs).chunk()

    def test_chunk(self):
        unblocked = Variable(["dim_0", "dim_1"], np.ones((3, 4)))
        assert unblocked.chunks is None

        blocked = unblocked.chunk()
        assert blocked.chunks == ((3,), (4,))
        first_dask_name = blocked.data.name

        blocked = unblocked.chunk(chunks=((2, 1), (2, 2)))  # type: ignore[arg-type]
        assert blocked.chunks == ((2, 1), (2, 2))
        assert blocked.data.name != first_dask_name

        blocked = unblocked.chunk(chunks=(3, 3))
        assert blocked.chunks == ((3,), (3, 1))
        assert blocked.data.name != first_dask_name

        # name doesn't change when rechunking by same amount
        # this fails if ReprObject doesn't have __dask_tokenize__ defined
        assert unblocked.chunk(2).data.name == unblocked.chunk(2).data.name

        assert blocked.load().chunks is None

        # Check that kwargs are passed
        import dask.array as da

        blocked = unblocked.chunk(name="testname_")
        assert isinstance(blocked.data, da.Array)
        assert "testname_" in blocked.data.name

        # test kwargs form of chunks
        blocked = unblocked.chunk(dim_0=3, dim_1=3)
        assert blocked.chunks == ((3,), (3, 1))
        assert blocked.data.name != first_dask_name

    @pytest.mark.skip
    def test_0d_object_array_with_list(self):
        super().test_0d_object_array_with_list()

    @pytest.mark.skip
    def test_array_interface(self):
        # dask array does not have `argsort`
        super().test_array_interface()

    @pytest.mark.skip
    def test_copy_index(self):
        super().test_copy_index()

    @pytest.mark.skip
    @pytest.mark.filterwarnings("ignore:elementwise comparison failed.*:FutureWarning")
    def test_eq_all_dtypes(self):
        super().test_eq_all_dtypes()

    def test_getitem_fancy(self):
        super().test_getitem_fancy()

    def test_getitem_1d_fancy(self):
        super().test_getitem_1d_fancy()

    def test_getitem_with_mask_nd_indexer(self):
        import dask.array as da

        v = Variable(["x"], da.arange(3, chunks=3))
        indexer = Variable(("x", "y"), [[0, -1], [-1, 2]])
        assert_identical(
            v._getitem_with_mask(indexer, fill_value=-1),
            self.cls(("x", "y"), [[0, -1], [-1, 2]]),
        )

    @pytest.mark.parametrize("dim", ["x", "y"])
    @pytest.mark.parametrize("window", [3, 8, 11])
    @pytest.mark.parametrize("center", [True, False])
    def test_dask_rolling(self, dim, window, center):
        import dask
        import dask.array as da

        dask.config.set(scheduler="single-threaded")

        x = Variable(("x", "y"), np.array(np.random.randn(100, 40), dtype=float))
        dx = Variable(("x", "y"), da.from_array(x, chunks=[(6, 30, 30, 20, 14), 8]))

        expected = x.rolling_window(
            dim, window, "window", center=center, fill_value=np.nan
        )
        with raise_if_dask_computes():
            actual = dx.rolling_window(
                dim, window, "window", center=center, fill_value=np.nan
            )
        assert isinstance(actual.data, da.Array)
        assert actual.shape == expected.shape
        assert_equal(actual, expected)

    @pytest.mark.xfail(reason="https://github.com/dask/dask/issues/11585")
    def test_multiindex(self):
        super().test_multiindex()

    @pytest.mark.parametrize(
        "mode",
        [
            "mean",
            pytest.param(
                "median",
                marks=pytest.mark.xfail(reason="median is not implemented by Dask"),
            ),
            pytest.param(
                "reflect", marks=pytest.mark.xfail(reason="dask.array.pad bug")
            ),
            "edge",
            "linear_ramp",
            "maximum",
            "minimum",
            "symmetric",
            "wrap",
        ],
    )
    @pytest.mark.parametrize("xr_arg, np_arg", _PAD_XR_NP_ARGS)
    @pytest.mark.filterwarnings(
        r"ignore:dask.array.pad.+? converts integers to floats."
    )
    def test_pad(self, mode, xr_arg, np_arg):
        super().test_pad(mode, xr_arg, np_arg)

    @pytest.mark.skip(reason="dask doesn't support extension arrays")
    def test_pandas_period_index(self):
        super().test_pandas_period_index()

    @pytest.mark.skip(reason="dask doesn't support extension arrays")
    def test_pandas_datetime64_with_tz(self):
        super().test_pandas_datetime64_with_tz()

    @pytest.mark.skip(reason="dask doesn't support extension arrays")
    def test_pandas_categorical_dtype(self):
        super().test_pandas_categorical_dtype()


@requires_sparse
class TestVariableWithSparse:
    # TODO inherit VariableSubclassobjects to cover more tests

    def test_as_sparse(self):
        data = np.arange(12).reshape(3, 4)
        var = Variable(("x", "y"), data)._as_sparse(fill_value=-1)
        actual = var._to_dense()
        assert_identical(var, actual)


class TestIndexVariable(VariableSubclassobjects):
    def cls(self, *args, **kwargs) -> IndexVariable:
        return IndexVariable(*args, **kwargs)

    def test_init(self):
        with pytest.raises(ValueError, match=r"must be 1-dimensional"):
            IndexVariable((), 0)

    def test_to_index(self):
        data = 0.5 * np.arange(10)
        v = IndexVariable(["time"], data, {"foo": "bar"})
        assert pd.Index(data, name="time").identical(v.to_index())

    def test_to_index_multiindex_level(self):
        midx = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two"))
        with pytest.warns(FutureWarning):
            ds = Dataset(coords={"x": midx})
        assert ds.one.variable.to_index().equals(midx.get_level_values("one"))

    def test_multiindex_default_level_names(self):
        midx = pd.MultiIndex.from_product([["a", "b"], [1, 2]])
        v = IndexVariable(["x"], midx, {"foo": "bar"})
        assert v.to_index().names == ("x_level_0", "x_level_1")

    def test_data(self):  # type: ignore[override]
        x = IndexVariable("x", np.arange(3.0))
        assert isinstance(x._data, PandasIndexingAdapter)
        assert isinstance(x.data, np.ndarray)
        assert float == x.dtype
        assert_array_equal(np.arange(3), x)
        assert float == x.values.dtype
        with pytest.raises(TypeError, match=r"cannot be modified"):
            x[:] = 0

    def test_name(self):
        coord = IndexVariable("x", [10.0])
        assert coord.name == "x"

        with pytest.raises(AttributeError):
            coord.name = "y"

    def test_level_names(self):
        midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=["level_1", "level_2"]
        )
        x = IndexVariable("x", midx)
        assert x.level_names == midx.names

        assert IndexVariable("y", [10.0]).level_names is None

    def test_get_level_variable(self):
        midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=["level_1", "level_2"]
        )
        x = IndexVariable("x", midx)
        level_1 = IndexVariable("x", midx.get_level_values("level_1"))
        assert_identical(x.get_level_variable("level_1"), level_1)

        with pytest.raises(ValueError, match=r"has no MultiIndex"):
            IndexVariable("y", [10.0]).get_level_variable("level")

    def test_concat_periods(self):
        periods = pd.period_range("2000-01-01", periods=10)
        coords = [IndexVariable("t", periods[:5]), IndexVariable("t", periods[5:])]
        expected = IndexVariable("t", periods)
        actual = IndexVariable.concat(coords, dim="t")
        assert_identical(actual, expected)
        assert isinstance(actual.to_index(), pd.PeriodIndex)

        positions = [list(range(5)), list(range(5, 10))]
        actual = IndexVariable.concat(coords, dim="t", positions=positions)
        assert_identical(actual, expected)
        assert isinstance(actual.to_index(), pd.PeriodIndex)

    def test_concat_multiindex(self):
        idx = pd.MultiIndex.from_product([[0, 1, 2], ["a", "b"]])
        coords = [IndexVariable("x", idx[:2]), IndexVariable("x", idx[2:])]
        expected = IndexVariable("x", idx)
        actual = IndexVariable.concat(coords, dim="x")
        assert_identical(actual, expected)
        assert isinstance(actual.to_index(), pd.MultiIndex)

    @pytest.mark.parametrize("dtype", [str, bytes])
    def test_concat_str_dtype(self, dtype):
        a = IndexVariable("x", np.array(["a"], dtype=dtype))
        b = IndexVariable("x", np.array(["b"], dtype=dtype))
        expected = IndexVariable("x", np.array(["a", "b"], dtype=dtype))

        actual = IndexVariable.concat([a, b])
        assert actual.identical(expected)
        assert np.issubdtype(actual.dtype, dtype)

    def test_datetime64(self):
        # GH:1932  Make sure indexing keeps precision
        t = np.array([1518418799999986560, 1518418799999996560], dtype="datetime64[ns]")
        v = IndexVariable("t", t)
        assert v[0].data == t[0]

    # These tests make use of multi-dimensional variables, which are not valid
    # IndexVariable objects:
    @pytest.mark.skip
    def test_getitem_error(self):
        super().test_getitem_error()

    @pytest.mark.skip
    def test_getitem_advanced(self):
        super().test_getitem_advanced()

    @pytest.mark.skip
    def test_getitem_fancy(self):
        super().test_getitem_fancy()

    @pytest.mark.skip
    def test_getitem_uint(self):
        super().test_getitem_fancy()

    @pytest.mark.skip
    @pytest.mark.parametrize(
        "mode",
        [
            "mean",
            "median",
            "reflect",
            "edge",
            "linear_ramp",
            "maximum",
            "minimum",
            "symmetric",
            "wrap",
        ],
    )
    @pytest.mark.parametrize("xr_arg, np_arg", _PAD_XR_NP_ARGS)
    def test_pad(self, mode, xr_arg, np_arg):
        super().test_pad(mode, xr_arg, np_arg)

    @pytest.mark.skip
    def test_pad_constant_values(self, xr_arg, np_arg):
        super().test_pad_constant_values(xr_arg, np_arg)

    @pytest.mark.skip
    def test_rolling_window(self):
        super().test_rolling_window()

    @pytest.mark.skip
    def test_rolling_1d(self):
        super().test_rolling_1d()

    @pytest.mark.skip
    def test_nd_rolling(self):
        super().test_nd_rolling()

    @pytest.mark.skip
    def test_rolling_window_errors(self):
        super().test_rolling_window_errors()

    @pytest.mark.skip
    def test_coarsen_2d(self):
        super().test_coarsen_2d()  # type: ignore[misc]

    def test_to_index_variable_copy(self) -> None:
        # to_index_variable should return a copy
        # https://github.com/pydata/xarray/issues/6931
        a = IndexVariable("x", ["a"])
        b = a.to_index_variable()
        assert a is not b
        b.dims = ("y",)
        assert a.dims == ("x",)


class TestAsCompatibleData(Generic[T_DuckArray]):
    def test_unchanged_types(self):
        types = (np.asarray, PandasIndexingAdapter, LazilyIndexedArray)
        for t in types:
            for data in [
                np.arange(3),
                pd.date_range("2000-01-01", periods=3),
                pd.date_range("2000-01-01", periods=3).values,
            ]:
                x = t(data)  # type: ignore[arg-type]
                assert source_ndarray(x) is source_ndarray(as_compatible_data(x))

    def test_converted_types(self):
        for input_array in [
            [[0, 1, 2]],
            pd.DataFrame([[0, 1, 2]]),
            np.float64(1.4),
            np.str_("abc"),
        ]:
            actual = as_compatible_data(input_array)
            assert_array_equal(np.asarray(input_array), actual)
            assert np.ndarray is type(actual)
            assert np.asarray(input_array).dtype == actual.dtype

    def test_masked_array(self):
        original: Any = np.ma.MaskedArray(np.arange(5))
        expected = np.arange(5)
        actual: Any = as_compatible_data(original)
        assert_array_equal(expected, actual)
        assert np.dtype(int) == actual.dtype

        original1: Any = np.ma.MaskedArray(np.arange(5), mask=4 * [False] + [True])
        expected1: Any = np.arange(5.0)
        expected1[-1] = np.nan
        actual = as_compatible_data(original1)
        assert_array_equal(expected1, actual)
        assert np.dtype(float) == actual.dtype

        original2: Any = np.ma.MaskedArray([1.0, 2.0], mask=[True, False])
        original2.flags.writeable = False
        expected2: Any = [np.nan, 2.0]
        actual = as_compatible_data(original2)
        assert_array_equal(expected2, actual)
        assert np.dtype(float) == actual.dtype

        # GH2377
        actual_var: Any = Variable(dims=tuple(), data=np.ma.masked)
        expected_var = Variable(dims=tuple(), data=np.nan)
        assert_array_equal(expected_var, actual_var)
        assert actual_var.dtype == expected_var.dtype

    def test_datetime(self):
        expected = np.datetime64("2000-01-01")
        actual: Any = as_compatible_data(expected)
        assert expected == actual
        assert np.ndarray is type(actual)
        assert np.dtype("datetime64[s]") == actual.dtype

        expected_dt: Any = np.array([np.datetime64("2000-01-01")])
        actual = as_compatible_data(expected_dt)
        assert np.asarray(expected_dt) == actual
        assert np.ndarray is type(actual)
        assert np.dtype("datetime64[s]") == actual.dtype

        expected_dt_ns: Any = np.array([np.datetime64("2000-01-01", "ns")])
        actual = as_compatible_data(expected_dt_ns)
        assert np.asarray(expected_dt_ns) == actual
        assert np.ndarray is type(actual)
        assert np.dtype("datetime64[ns]") == actual.dtype
        assert expected_dt_ns is source_ndarray(np.asarray(actual))

        expected = np.datetime64(
            "2000-01-01",
            "us" if has_pandas_3 else "ns",
        )
        actual = as_compatible_data(datetime(2000, 1, 1))
        assert np.asarray(expected) == actual
        assert np.ndarray is type(actual)
        assert expected.dtype == actual.dtype

    def test_tz_datetime(self) -> None:
        tz = pytz.timezone("America/New_York")
        times_ns = pd.date_range("2000", periods=1, tz=tz)

        times_s = times_ns.astype(pd.DatetimeTZDtype("s", tz))  # type: ignore[arg-type]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            actual: T_DuckArray = as_compatible_data(times_s)
        assert actual.array == times_s
        assert actual.array.dtype == pd.DatetimeTZDtype("s", tz)  # type: ignore[arg-type]

        series = pd.Series(times_s)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            actual2: T_DuckArray = as_compatible_data(series)

        np.testing.assert_array_equal(actual2, np.asarray(series.values))
        assert actual2.dtype == np.dtype("datetime64[s]")

    def test_full_like(self) -> None:
        # For more thorough tests, see test_variable.py
        orig = Variable(
            dims=("x", "y"), data=[[1.5, 2.0], [3.1, 4.3]], attrs={"foo": "bar"}
        )

        expect = orig.copy(deep=True)
        # see https://github.com/python/mypy/issues/3004 for why we need to ignore type
        expect.values = [[2.0, 2.0], [2.0, 2.0]]  # type: ignore[assignment,unused-ignore]
        assert_identical(expect, full_like(orig, 2))

        # override dtype
        expect.values = [[True, True], [True, True]]  # type: ignore[assignment,unused-ignore]
        assert expect.dtype == bool
        assert_identical(expect, full_like(orig, True, dtype=bool))

        # raise error on non-scalar fill_value
        with pytest.raises(ValueError, match=r"must be scalar"):
            full_like(orig, [1.0, 2.0])

        with pytest.raises(ValueError, match="'dtype' cannot be dict-like"):
            full_like(orig, True, dtype={"x": bool})

    @requires_dask
    def test_full_like_dask(self) -> None:
        orig = Variable(
            dims=("x", "y"), data=[[1.5, 2.0], [3.1, 4.3]], attrs={"foo": "bar"}
        ).chunk(dict(x=(1, 1), y=(2,)))

        def check(actual, expect_dtype, expect_values):
            assert actual.dtype == expect_dtype
            assert actual.shape == orig.shape
            assert actual.dims == orig.dims
            assert actual.attrs == orig.attrs
            assert actual.chunks == orig.chunks
            assert_array_equal(actual.values, expect_values)

        check(full_like(orig, 2), orig.dtype, np.full_like(orig.values, 2))
        # override dtype
        check(
            full_like(orig, True, dtype=bool),
            bool,
            np.full_like(orig.values, True, dtype=bool),
        )

        # Check that there's no array stored inside dask
        # (e.g. we didn't create a numpy array and then we chunked it!)
        dsk = full_like(orig, 1).data.dask
        for v in dsk.values():
            if isinstance(v, tuple):
                for vi in v:
                    assert not isinstance(vi, np.ndarray)
            else:
                assert not isinstance(v, np.ndarray)

    def test_zeros_like(self) -> None:
        orig = Variable(
            dims=("x", "y"), data=[[1.5, 2.0], [3.1, 4.3]], attrs={"foo": "bar"}
        )
        assert_identical(zeros_like(orig), full_like(orig, 0))
        assert_identical(zeros_like(orig, dtype=int), full_like(orig, 0, dtype=int))

    def test_ones_like(self) -> None:
        orig = Variable(
            dims=("x", "y"), data=[[1.5, 2.0], [3.1, 4.3]], attrs={"foo": "bar"}
        )
        assert_identical(ones_like(orig), full_like(orig, 1))
        assert_identical(ones_like(orig, dtype=int), full_like(orig, 1, dtype=int))

    def test_numpy_ndarray_subclass(self):
        class SubclassedArray(np.ndarray):
            def __new__(cls, array, foo):
                obj = np.asarray(array).view(cls)
                obj.foo = foo  # type: ignore[attr-defined]
                return obj

        data = SubclassedArray([1, 2, 3], foo="bar")
        actual: Any = as_compatible_data(data)
        assert isinstance(actual, SubclassedArray)
        assert actual.foo == "bar"  # type: ignore[attr-defined]
        assert_array_equal(data, actual)

    def test_numpy_matrix(self):
        with pytest.warns(PendingDeprecationWarning):
            data = np.matrix([[1, 2], [3, 4]])
        actual: Any = as_compatible_data(data)
        assert isinstance(actual, np.ndarray)
        assert_array_equal(data, actual)

    def test_unsupported_type(self):
        # Non indexable type
        class CustomArray(NDArrayMixin):
            def __init__(self, array):
                self.array = array

        class CustomIndexable(CustomArray, indexing.ExplicitlyIndexed):
            pass

        # Type with data stored in values attribute
        class CustomWithValuesAttr:
            def __init__(self, array):
                self.values = array

        array = CustomArray(np.arange(3))
        orig = Variable(dims=("x"), data=array, attrs={"foo": "bar"})
        assert isinstance(orig._data, np.ndarray)  # should not be CustomArray

        array = CustomIndexable(np.arange(3))
        orig = Variable(dims=("x"), data=array, attrs={"foo": "bar"})
        assert isinstance(orig._data, CustomIndexable)

        array2: Any = CustomWithValuesAttr(np.arange(3))
        orig = Variable(dims=(), data=array2)
        assert isinstance(orig._data.item(), CustomWithValuesAttr)  # type: ignore[union-attr]


def test_raise_no_warning_for_nan_in_binary_ops():
    with assert_no_warnings():
        _ = Variable("x", [1, 2, np.nan]) > 0


class TestBackendIndexing:
    """Make sure all the array wrappers can be indexed."""

    @pytest.fixture(autouse=True)
    def setUp(self):
        self.d = np.random.random((10, 3)).astype(np.float64)
        self.cat = PandasExtensionArray(pd.Categorical(["a", "b"] * 5))

    async def check_orthogonal_indexing(self, v, load_async):
        expected = self.d[[8, 3]][:, [2, 1]]

        if load_async:
            result = await v.isel(x=[8, 3], y=[2, 1]).load_async()
        else:
            result = v.isel(x=[8, 3], y=[2, 1])

        assert np.allclose(result, expected)

    async def check_vectorized_indexing(self, v, load_async):
        ind_x = Variable("z", [0, 2])
        ind_y = Variable("z", [2, 1])
        expected = self.d[ind_x, ind_y]

        if load_async:
            result = await v.isel(x=ind_x, y=ind_y).load_async()
        else:
            result = v.isel(x=ind_x, y=ind_y).load()

        assert np.allclose(result, expected)

    @pytest.mark.asyncio
    @pytest.mark.parametrize("load_async", [True, False])
    async def test_NumpyIndexingAdapter(self, load_async):
        v = Variable(dims=("x", "y"), data=NumpyIndexingAdapter(self.d))
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)
        # could not doubly wrapping
        with pytest.raises(TypeError, match=r"NumpyIndexingAdapter only wraps "):
            v = Variable(
                dims=("x", "y"), data=NumpyIndexingAdapter(NumpyIndexingAdapter(self.d))
            )

    def test_extension_array_duck_array(self):
        lazy = LazilyIndexedArray(self.cat)
        assert (lazy.get_duck_array().array == self.cat).all()

    def test_extension_array_duck_indexed(self):
        lazy = Variable(dims=("x"), data=LazilyIndexedArray(self.cat))
        assert (lazy[[0, 1, 5]] == ["a", "b", "b"]).all()

    @pytest.mark.asyncio
    @pytest.mark.parametrize("load_async", [True, False])
    async def test_LazilyIndexedArray(self, load_async):
        v = Variable(dims=("x", "y"), data=LazilyIndexedArray(self.d))
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)
        # doubly wrapping
        v = Variable(
            dims=("x", "y"),
            data=LazilyIndexedArray(LazilyIndexedArray(self.d)),
        )
        await self.check_orthogonal_indexing(v, load_async)
        # hierarchical wrapping
        v = Variable(
            dims=("x", "y"), data=LazilyIndexedArray(NumpyIndexingAdapter(self.d))
        )
        await self.check_orthogonal_indexing(v, load_async)

    @pytest.mark.asyncio
    @pytest.mark.parametrize("load_async", [True, False])
    async def test_CopyOnWriteArray(self, load_async):
        v = Variable(dims=("x", "y"), data=CopyOnWriteArray(self.d))
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)
        # doubly wrapping
        v = Variable(dims=("x", "y"), data=CopyOnWriteArray(LazilyIndexedArray(self.d)))  # type: ignore[arg-type]
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)

    @pytest.mark.asyncio
    @pytest.mark.parametrize("load_async", [True, False])
    async def test_MemoryCachedArray(self, load_async):
        v = Variable(dims=("x", "y"), data=MemoryCachedArray(self.d))
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)
        # doubly wrapping
        v = Variable(dims=("x", "y"), data=CopyOnWriteArray(MemoryCachedArray(self.d)))  # type: ignore[arg-type]
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)

    @requires_dask
    @pytest.mark.asyncio
    @pytest.mark.parametrize("load_async", [True, False])
    async def test_DaskIndexingAdapter(self, load_async):
        import dask.array as da

        dask_array = da.asarray(self.d)
        v = Variable(dims=("x", "y"), data=DaskIndexingAdapter(dask_array))
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)
        # doubly wrapping
        v = Variable(
            dims=("x", "y"),
            data=CopyOnWriteArray(DaskIndexingAdapter(dask_array)),  # type: ignore[arg-type]
        )
        await self.check_orthogonal_indexing(v, load_async)
        await self.check_vectorized_indexing(v, load_async)


def test_clip(var):
    # Copied from test_dataarray (would there be a way to combine the tests?)
    result = var.clip(min=0.5)
    assert result.min(...) >= 0.5

    result = var.clip(max=0.5)
    assert result.max(...) <= 0.5

    result = var.clip(min=0.25, max=0.75)
    assert result.min(...) >= 0.25
    assert result.max(...) <= 0.75

    result = var.clip(min=var.mean("x"), max=var.mean("z"))
    assert result.dims == var.dims
    assert_array_equal(
        result.data,
        np.clip(
            var.data,
            var.mean("x").data[np.newaxis, :, :],
            var.mean("z").data[:, :, np.newaxis],
        ),
    )


@pytest.mark.parametrize("Var", [Variable, IndexVariable])
class TestNumpyCoercion:
    def test_from_numpy(self, Var):
        v = Var("x", [1, 2, 3])

        assert_identical(v.as_numpy(), v)
        np.testing.assert_equal(v.to_numpy(), np.array([1, 2, 3]))

    @requires_dask
    def test_from_dask(self, Var):
        v = Var("x", [1, 2, 3])
        v_chunked = v.chunk(1)

        assert_identical(v_chunked.as_numpy(), v.compute())
        np.testing.assert_equal(v.to_numpy(), np.array([1, 2, 3]))

    @requires_pint
    def test_from_pint(self, Var):
        import pint

        arr = np.array([1, 2, 3])

        # IndexVariable strips the unit
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=pint.UnitStrippedWarning)
            v = Var("x", pint.Quantity(arr, units="m"))

        assert_identical(v.as_numpy(), Var("x", arr))
        np.testing.assert_equal(v.to_numpy(), arr)

    @requires_sparse
    def test_from_sparse(self, Var):
        if Var is IndexVariable:
            pytest.skip("Can't have 2D IndexVariables")

        import sparse

        arr = np.diagflat([1, 2, 3])
        coords = np.array([[0, 1, 2], [0, 1, 2]])
        sparr = sparse.COO(coords=coords, data=[1, 2, 3], shape=(3, 3))
        v = Variable(["x", "y"], sparr)

        assert_identical(v.as_numpy(), Variable(["x", "y"], arr))
        np.testing.assert_equal(v.to_numpy(), arr)

    @requires_cupy
    def test_from_cupy(self, Var):
        if Var is IndexVariable:
            pytest.skip("cupy in default indexes is not supported at the moment")
        import cupy as cp

        arr = np.array([1, 2, 3])
        v = Var("x", cp.array(arr))

        assert_identical(v.as_numpy(), Var("x", arr))
        np.testing.assert_equal(v.to_numpy(), arr)

    @requires_dask
    @requires_pint
    def test_from_pint_wrapping_dask(self, Var):
        import dask
        import pint

        arr = np.array([1, 2, 3])
        d = dask.array.from_array(np.array([1, 2, 3]))

        # IndexVariable strips the unit
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=pint.UnitStrippedWarning)
            v = Var("x", pint.Quantity(d, units="m"))

        result = v.as_numpy()
        assert_identical(result, Var("x", arr))
        np.testing.assert_equal(v.to_numpy(), arr)


@pytest.mark.parametrize(
    ("values", "unit"),
    [
        (np.datetime64("2000-01-01", "ns"), "ns"),
        (np.datetime64("2000-01-01", "s"), "s"),
        (np.array([np.datetime64("2000-01-01", "ns")]), "ns"),
        (np.array([np.datetime64("2000-01-01", "s")]), "s"),
        (pd.date_range("2000", periods=1, unit="ns"), "ns"),
        (
            datetime(2000, 1, 1),
            "us" if has_pandas_3 else "ns",
        ),
        (
            np.array([datetime(2000, 1, 1)]),
            "us" if has_pandas_3 else "ns",
        ),
        (
            pd.date_range(
                "2000", periods=1, tz=pytz.timezone("America/New_York"), unit="ns"
            ),
            "ns",
        ),
        (
            pd.Series(
                pd.date_range(
                    "2000", periods=1, tz=pytz.timezone("America/New_York"), unit="ns"
                )
            ),
            "ns",
        ),
    ],
    ids=lambda x: f"{x}",
)
def test_datetime_conversion(values, unit) -> None:
    # todo: check for redundancy (suggested per review)
    dims = ["time"] if isinstance(values, np.ndarray | pd.Index | pd.Series) else []
    var = Variable(dims, values)
    if var.dtype.kind == "M" and isinstance(var.dtype, np.dtype):
        assert var.dtype == np.dtype(f"datetime64[{unit}]")
    else:
        # The only case where a non-datetime64 dtype can occur currently is in
        # the case that the variable is backed by a timezone-aware
        # DatetimeIndex, and thus is hidden within the PandasIndexingAdapter class.
        assert isinstance(var._data, PandasIndexingAdapter)
        assert var._data.array.dtype == pd.DatetimeTZDtype(
            "ns", pytz.timezone("America/New_York")
        )


tz_ny = pytz.timezone("America/New_York")


@pytest.mark.parametrize(
    ["data", "dtype"],
    [
        pytest.param(pd.date_range("2000", periods=1), "datetime64[s]", id="index-sec"),
        pytest.param(
            pd.Series(pd.date_range("2000", periods=1)),
            "datetime64[s]",
            id="series-sec",
        ),
        pytest.param(
            pd.date_range("2000", periods=1, tz=tz_ny),
            pd.DatetimeTZDtype("s", tz_ny),  # type: ignore[arg-type]
            id="index-timezone",
        ),
        pytest.param(
            pd.Series(pd.date_range("2000", periods=1, tz=tz_ny)),
            pd.DatetimeTZDtype("s", tz_ny),  # type: ignore[arg-type]
            id="series-timezone",
        ),
    ],
)
def test_pandas_two_only_datetime_conversion_warnings(
    data: pd.DatetimeIndex | pd.Series, dtype: str | pd.DatetimeTZDtype
) -> None:
    # todo: check for redundancy (suggested per review)
    var = Variable(["time"], data.astype(dtype))  # type: ignore[arg-type]

    # we internally convert series to numpy representations to avoid too much nastiness with extension arrays
    # when calling data.array e.g., with NumpyExtensionArrays
    if isinstance(data, pd.Series):
        assert var.dtype == np.dtype("datetime64[s]")
    elif var.dtype.kind == "M":
        assert var.dtype == dtype
    else:
        # The only case where a non-datetime64 dtype can occur currently is in
        # the case that the variable is backed by a timezone-aware
        # DatetimeIndex, and thus is hidden within the PandasIndexingAdapter class.
        assert isinstance(var._data, PandasIndexingAdapter)
        assert var._data.array.dtype == pd.DatetimeTZDtype("s", tz_ny)


@pytest.mark.parametrize(
    ("values", "unit"),
    [
        (np.timedelta64(10, "ns"), "ns"),
        (np.timedelta64(10, "s"), "s"),
        (np.array([np.timedelta64(10, "ns")]), "ns"),
        (np.array([np.timedelta64(10, "s")]), "s"),
        (pd.timedelta_range("1", periods=1), "ns"),
        (timedelta(days=1), "us" if has_pandas_3 else "ns"),
        (np.array([timedelta(days=1)]), "us" if has_pandas_3 else "ns"),
        (pd.timedelta_range("1", periods=1).astype("timedelta64[s]"), "s"),
    ],
    ids=lambda x: f"{x}",
)
def test_timedelta_conversion(values, unit) -> None:
    # todo: check for redundancy
    dims = ["time"] if isinstance(values, np.ndarray | pd.Index) else []
    var = Variable(dims, values)
    assert var.dtype == np.dtype(f"timedelta64[{unit}]")


def test_explicitly_indexed_array_preserved() -> None:
    """Test that methods using ._data preserve ExplicitlyIndexed arrays.

    Regression test for methods that should use ._data instead of .data
    to avoid loading lazy arrays into memory.
    """
    arr = IndexableArray(np.array([1, 2, 3]))
    var = Variable(["x"], arr)
    result = var.drop_encoding()
    assert isinstance(result._data, indexing.ExplicitlyIndexed)
