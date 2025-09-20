from __future__ import annotations

import pickle
import re
import sys
import warnings
from collections.abc import Hashable
from copy import copy, deepcopy
from io import StringIO
from textwrap import dedent
from typing import Any, Literal, cast

import numpy as np
import pandas as pd
import pytest
from packaging.version import Version
from pandas.core.indexes.datetimes import DatetimeIndex

# remove once numpy 2.0 is the oldest supported version
try:
    from numpy.exceptions import RankWarning
except ImportError:
    from numpy import RankWarning  # type: ignore[no-redef,attr-defined,unused-ignore]

import contextlib

from pandas.errors import UndefinedVariableError

import xarray as xr
from xarray import (
    AlignmentError,
    DataArray,
    Dataset,
    IndexVariable,
    MergeError,
    Variable,
    align,
    backends,
    broadcast,
    open_dataset,
    set_options,
)
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.core import dtypes, indexing, utils
from xarray.core.common import duck_array_ops, full_like
from xarray.core.coordinates import Coordinates, DatasetCoordinates
from xarray.core.indexes import Index, PandasIndex
from xarray.core.types import ArrayLike
from xarray.core.utils import is_scalar
from xarray.groupers import SeasonResampler, TimeResampler
from xarray.namedarray.pycompat import array_type, integer_types
from xarray.testing import _assert_internal_invariants
from xarray.tests import (
    DuckArrayWrapper,
    InaccessibleArray,
    UnexpectedDataAccess,
    assert_allclose,
    assert_array_equal,
    assert_equal,
    assert_identical,
    assert_no_warnings,
    assert_writeable,
    create_test_data,
    has_cftime,
    has_dask,
    has_pyarrow,
    raise_if_dask_computes,
    requires_bottleneck,
    requires_cftime,
    requires_cupy,
    requires_dask,
    requires_numexpr,
    requires_pint,
    requires_scipy,
    requires_sparse,
    source_ndarray,
)
from xarray.tests.indexes import ScalarIndex, XYIndex

with contextlib.suppress(ImportError):
    import dask.array as da

# from numpy version 2.0 trapz is deprecated and renamed to trapezoid
# remove once numpy 2.0 is the oldest supported version
try:
    from numpy import trapezoid  # type: ignore[attr-defined,unused-ignore]
except ImportError:
    from numpy import (  # type: ignore[arg-type,no-redef,attr-defined,unused-ignore]
        trapz as trapezoid,
    )

sparse_array_type = array_type("sparse")

pytestmark = [
    pytest.mark.filterwarnings("error:Mean of empty slice"),
    pytest.mark.filterwarnings("error:All-NaN (slice|axis) encountered"),
]


def create_append_test_data(seed=None) -> tuple[Dataset, Dataset, Dataset]:
    rs = np.random.default_rng(seed)

    lat = [2, 1, 0]
    lon = [0, 1, 2]
    nt1 = 3
    nt2 = 2
    time1 = pd.date_range("2000-01-01", periods=nt1).as_unit("ns")
    time2 = pd.date_range("2000-02-01", periods=nt2).as_unit("ns")
    string_var = np.array(["a", "bc", "def"], dtype=object)
    string_var_to_append = np.array(["asdf", "asdfg"], dtype=object)
    string_var_fixed_length = np.array(["aa", "bb", "cc"], dtype="|S2")
    string_var_fixed_length_to_append = np.array(["dd", "ee"], dtype="|S2")
    unicode_var = np.array(["áó", "áó", "áó"])
    datetime_var = np.array(
        ["2019-01-01", "2019-01-02", "2019-01-03"], dtype="datetime64[ns]"
    )
    datetime_var_to_append = np.array(
        ["2019-01-04", "2019-01-05"], dtype="datetime64[ns]"
    )
    bool_var = np.array([True, False, True], dtype=bool)
    bool_var_to_append = np.array([False, True], dtype=bool)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Converting non-default")
        ds = xr.Dataset(
            data_vars={
                "da": xr.DataArray(
                    rs.random((3, 3, nt1)),
                    coords=[lat, lon, time1],
                    dims=["lat", "lon", "time"],
                ),
                "string_var": ("time", string_var),
                "string_var_fixed_length": ("time", string_var_fixed_length),
                "unicode_var": ("time", unicode_var),
                "datetime_var": ("time", datetime_var),
                "bool_var": ("time", bool_var),
            }
        )

        ds_to_append = xr.Dataset(
            data_vars={
                "da": xr.DataArray(
                    rs.random((3, 3, nt2)),
                    coords=[lat, lon, time2],
                    dims=["lat", "lon", "time"],
                ),
                "string_var": ("time", string_var_to_append),
                "string_var_fixed_length": ("time", string_var_fixed_length_to_append),
                "unicode_var": ("time", unicode_var[:nt2]),
                "datetime_var": ("time", datetime_var_to_append),
                "bool_var": ("time", bool_var_to_append),
            }
        )

        ds_with_new_var = xr.Dataset(
            data_vars={
                "new_var": xr.DataArray(
                    rs.random((3, 3, nt1 + nt2)),
                    coords=[lat, lon, time1.append(time2)],
                    dims=["lat", "lon", "time"],
                )
            }
        )

    assert_writeable(ds)
    assert_writeable(ds_to_append)
    assert_writeable(ds_with_new_var)
    return ds, ds_to_append, ds_with_new_var


def create_append_string_length_mismatch_test_data(dtype) -> tuple[Dataset, Dataset]:
    def make_datasets(data, data_to_append) -> tuple[Dataset, Dataset]:
        ds = xr.Dataset(
            {"temperature": (["time"], data)},
            coords={"time": [0, 1, 2]},
        )
        ds_to_append = xr.Dataset(
            {"temperature": (["time"], data_to_append)}, coords={"time": [0, 1, 2]}
        )
        assert_writeable(ds)
        assert_writeable(ds_to_append)
        return ds, ds_to_append

    u2_strings = ["ab", "cd", "ef"]
    u5_strings = ["abc", "def", "ghijk"]

    s2_strings = np.array(["aa", "bb", "cc"], dtype="|S2")
    s3_strings = np.array(["aaa", "bbb", "ccc"], dtype="|S3")

    if dtype == "U":
        return make_datasets(u2_strings, u5_strings)
    elif dtype == "S":
        return make_datasets(s2_strings, s3_strings)
    else:
        raise ValueError(f"unsupported dtype {dtype}.")


def create_test_multiindex() -> Dataset:
    mindex = pd.MultiIndex.from_product(
        [["a", "b"], [1, 2]], names=("level_1", "level_2")
    )
    return Dataset({}, Coordinates.from_pandas_multiindex(mindex, "x"))


def create_test_stacked_array() -> tuple[DataArray, DataArray]:
    x = DataArray(pd.Index(np.r_[:10], name="x"))
    y = DataArray(pd.Index(np.r_[:20], name="y"))
    a = x * y
    b = x * y * y
    return a, b


class InaccessibleVariableDataStore(backends.InMemoryDataStore):
    """
    Store that does not allow any data access.
    """

    def __init__(self):
        super().__init__()
        self._indexvars = set()

    def store(self, variables, *args, **kwargs) -> None:
        super().store(variables, *args, **kwargs)
        for k, v in variables.items():
            if isinstance(v, IndexVariable):
                self._indexvars.add(k)

    def get_variables(self):
        def lazy_inaccessible(k, v):
            if k in self._indexvars:
                return v
            data = indexing.LazilyIndexedArray(InaccessibleArray(v.values))
            return Variable(v.dims, data, v.attrs)

        return {k: lazy_inaccessible(k, v) for k, v in self._variables.items()}


class DuckBackendArrayWrapper(backends.common.BackendArray):
    """Mimic a BackendArray wrapper around DuckArrayWrapper"""

    def __init__(self, array):
        self.array = DuckArrayWrapper(array)
        self.shape = array.shape
        self.dtype = array.dtype

    def get_array(self):
        return self.array

    def __getitem__(self, key):
        return self.array[key.tuple]


class AccessibleAsDuckArrayDataStore(backends.InMemoryDataStore):
    """
    Store that returns a duck array, not convertible to numpy array,
    on read. Modeled after nVIDIA's kvikio.
    """

    def __init__(self):
        super().__init__()
        self._indexvars = set()

    def store(self, variables, *args, **kwargs) -> None:
        super().store(variables, *args, **kwargs)
        for k, v in variables.items():
            if isinstance(v, IndexVariable):
                self._indexvars.add(k)

    def get_variables(self) -> dict[Any, xr.Variable]:
        def lazy_accessible(k, v) -> xr.Variable:
            if k in self._indexvars:
                return v
            data = indexing.LazilyIndexedArray(DuckBackendArrayWrapper(v.values))
            return Variable(v.dims, data, v.attrs)

        return {k: lazy_accessible(k, v) for k, v in self._variables.items()}


class TestDataset:
    def test_repr(self) -> None:
        data = create_test_data(seed=123, use_extension_array=True)
        data.attrs["foo"] = "bar"
        # need to insert str dtype at runtime to handle different endianness
        var5 = (
            "\n                var5     (dim1) int64[pyarrow] 64B 5 9 7 2 6 2 8 1"
            if has_pyarrow
            else ""
        )
        expected = dedent(
            f"""\
            <xarray.Dataset> Size: 2kB
            Dimensions:  (dim2: 9, dim3: 10, time: 20, dim1: 8)
            Coordinates:
              * dim2     (dim2) float64 72B 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0
              * dim3     (dim3) {data["dim3"].dtype} 40B 'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j'
              * time     (time) datetime64[ns] 160B 2000-01-01 2000-01-02 ... 2000-01-20
                numbers  (dim3) int64 80B 0 1 2 0 0 1 1 2 2 3
            Dimensions without coordinates: dim1
            Data variables:
                var1     (dim1, dim2) float64 576B -0.9891 -0.3678 1.288 ... -0.2116 0.364
                var2     (dim1, dim2) float64 576B 0.953 1.52 1.704 ... 0.1347 -0.6423
                var3     (dim3, dim1) float64 640B 0.4107 0.9941 0.1665 ... 0.716 1.555
                var4     (dim1) category 3{6 if Version(pd.__version__) >= Version("3.0.0dev0") else 2}B b c b a c a c a{var5}
            Attributes:
                foo:      bar"""
        )
        actual = "\n".join(x.rstrip() for x in repr(data).split("\n"))

        assert expected == actual

        with set_options(display_width=100):
            max_len = max(map(len, repr(data).split("\n")))
            assert 90 < max_len < 100

        expected = dedent(
            """\
            <xarray.Dataset> Size: 0B
            Dimensions:  ()
            Data variables:
                *empty*"""
        )
        actual = "\n".join(x.rstrip() for x in repr(Dataset()).split("\n"))
        print(actual)
        assert expected == actual

        # verify that ... doesn't appear for scalar coordinates
        data = Dataset({"foo": ("x", np.ones(10))}).mean()
        expected = dedent(
            """\
            <xarray.Dataset> Size: 8B
            Dimensions:  ()
            Data variables:
                foo      float64 8B 1.0"""
        )
        actual = "\n".join(x.rstrip() for x in repr(data).split("\n"))
        print(actual)
        assert expected == actual

        # verify long attributes are truncated
        data = Dataset(attrs={"foo": "bar" * 1000})
        assert len(repr(data)) < 1000

    def test_repr_multiindex(self) -> None:
        data = create_test_multiindex()
        obj_size = np.dtype("O").itemsize
        expected = dedent(
            f"""\
            <xarray.Dataset> Size: {8 * obj_size + 32}B
            Dimensions:  (x: 4)
            Coordinates:
              * x        (x) object {4 * obj_size}B MultiIndex
              * level_1  (x) object {4 * obj_size}B 'a' 'a' 'b' 'b'
              * level_2  (x) int64 32B 1 2 1 2
            Data variables:
                *empty*"""
        )
        actual = "\n".join(x.rstrip() for x in repr(data).split("\n"))
        print(actual)
        assert expected == actual

        # verify that long level names are not truncated
        midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=("a_quite_long_level_name", "level_2")
        )
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        data = Dataset({}, midx_coords)
        expected = dedent(
            f"""\
            <xarray.Dataset> Size: {8 * obj_size + 32}B
            Dimensions:                  (x: 4)
            Coordinates:
              * x                        (x) object {4 * obj_size}B MultiIndex
              * a_quite_long_level_name  (x) object {4 * obj_size}B 'a' 'a' 'b' 'b'
              * level_2                  (x) int64 32B 1 2 1 2
            Data variables:
                *empty*"""
        )
        actual = "\n".join(x.rstrip() for x in repr(data).split("\n"))
        print(actual)
        assert expected == actual

    def test_repr_period_index(self) -> None:
        data = create_test_data(seed=456)
        data.coords["time"] = pd.period_range("2000-01-01", periods=20, freq="D")

        # check that creating the repr doesn't raise an error #GH645
        repr(data)

    def test_unicode_data(self) -> None:
        # regression test for GH834
        data = Dataset({"foø": ["ba®"]}, attrs={"å": "∑"})
        repr(data)  # should not raise

        byteorder = "<" if sys.byteorder == "little" else ">"
        expected = dedent(
            f"""\
            <xarray.Dataset> Size: 12B
            Dimensions:  (foø: 1)
            Coordinates:
              * foø      (foø) {byteorder}U3 12B {"ba®"!r}
            Data variables:
                *empty*
            Attributes:
                å:        ∑"""
        )
        actual = str(data)
        assert expected == actual

    def test_repr_nep18(self) -> None:
        class Array:
            def __init__(self):
                self.shape = (2,)
                self.ndim = 1
                self.dtype = np.dtype(np.float64)

            def __array_function__(self, *args, **kwargs):
                return NotImplemented

            def __array_ufunc__(self, *args, **kwargs):
                return NotImplemented

            def __repr__(self):
                return "Custom\nArray"

        dataset = Dataset({"foo": ("x", Array())})
        expected = dedent(
            """\
            <xarray.Dataset> Size: 16B
            Dimensions:  (x: 2)
            Dimensions without coordinates: x
            Data variables:
                foo      (x) float64 16B Custom Array"""
        )
        assert expected == repr(dataset)

    def test_info(self) -> None:
        ds = create_test_data(seed=123)
        ds = ds.drop_vars("dim3")  # string type prints differently in PY2 vs PY3
        ds.attrs["unicode_attr"] = "ba®"
        ds.attrs["string_attr"] = "bar"

        buf = StringIO()
        ds.info(buf=buf)

        expected = dedent(
            """\
        xarray.Dataset {
        dimensions:
        \tdim2 = 9 ;
        \ttime = 20 ;
        \tdim1 = 8 ;
        \tdim3 = 10 ;

        variables:
        \tfloat64 dim2(dim2) ;
        \tdatetime64[ns] time(time) ;
        \tfloat64 var1(dim1, dim2) ;
        \t\tvar1:foo = variable ;
        \tfloat64 var2(dim1, dim2) ;
        \t\tvar2:foo = variable ;
        \tfloat64 var3(dim3, dim1) ;
        \t\tvar3:foo = variable ;
        \tint64 numbers(dim3) ;

        // global attributes:
        \t:unicode_attr = ba® ;
        \t:string_attr = bar ;
        }"""
        )
        actual = buf.getvalue()
        assert expected == actual
        buf.close()

    def test_constructor(self) -> None:
        x1 = ("x", 2 * np.arange(100))
        x2 = ("x", np.arange(1000))
        z = (["x", "y"], np.arange(1000).reshape(100, 10))

        with pytest.raises(ValueError, match=r"conflicting sizes"):
            Dataset({"a": x1, "b": x2})
        with pytest.raises(TypeError, match=r"tuple of form"):
            Dataset({"x": (1, 2, 3, 4, 5, 6, 7)})
        with pytest.raises(ValueError, match=r"already exists as a scalar"):
            Dataset({"x": 0, "y": ("x", [1, 2, 3])})

        # nD coordinate variable "x" sharing name with dimension
        actual = Dataset({"a": x1, "x": z})
        assert "x" not in actual.xindexes
        _assert_internal_invariants(actual, check_default_indexes=True)

        # verify handling of DataArrays
        expected = Dataset({"x": x1, "z": z})
        actual = Dataset({"z": expected["z"]})
        assert_identical(expected, actual)

    def test_constructor_1d(self) -> None:
        expected = Dataset({"x": (["x"], 5.0 + np.arange(5))})
        actual = Dataset({"x": 5.0 + np.arange(5)})
        assert_identical(expected, actual)

        actual = Dataset({"x": [5, 6, 7, 8, 9]})
        assert_identical(expected, actual)

    def test_constructor_0d(self) -> None:
        expected = Dataset({"x": ([], 1)})
        for arg in [1, np.array(1), expected["x"]]:
            actual = Dataset({"x": arg})
            assert_identical(expected, actual)

        class Arbitrary:
            pass

        d = pd.Timestamp("2000-01-01T12")
        args = [
            True,
            None,
            3.4,
            np.nan,
            "hello",
            b"raw",
            np.datetime64("2000-01-01"),
            d,
            d.to_pydatetime(),
            Arbitrary(),
        ]
        for arg in args:
            print(arg)
            expected = Dataset({"x": ([], arg)})
            actual = Dataset({"x": arg})
            assert_identical(expected, actual)

    def test_constructor_auto_align(self) -> None:
        a = DataArray([1, 2], [("x", [0, 1])])
        b = DataArray([3, 4], [("x", [1, 2])])

        # verify align uses outer join
        expected = Dataset(
            {"a": ("x", [1, 2, np.nan]), "b": ("x", [np.nan, 3, 4])}, {"x": [0, 1, 2]}
        )
        actual = Dataset({"a": a, "b": b})
        assert_identical(expected, actual)

        # regression test for GH346
        assert isinstance(actual.variables["x"], IndexVariable)

        # variable with different dimensions
        c = ("y", [3, 4])
        expected2 = expected.merge({"c": c})
        actual = Dataset({"a": a, "b": b, "c": c})
        assert_identical(expected2, actual)

        # variable that is only aligned against the aligned variables
        d = ("x", [3, 2, 1])
        expected3 = expected.merge({"d": d})
        actual = Dataset({"a": a, "b": b, "d": d})
        assert_identical(expected3, actual)

        e = ("x", [0, 0])
        with pytest.raises(ValueError, match=r"conflicting sizes"):
            Dataset({"a": a, "b": b, "e": e})

    def test_constructor_pandas_sequence(self) -> None:
        ds = self.make_example_math_dataset()
        pandas_objs = {
            var_name: ds[var_name].to_pandas() for var_name in ["foo", "bar"]
        }
        ds_based_on_pandas = Dataset(pandas_objs, ds.coords, attrs=ds.attrs)
        del ds_based_on_pandas["x"]
        assert_equal(ds, ds_based_on_pandas)

        # reindex pandas obj, check align works
        rearranged_index = reversed(pandas_objs["foo"].index)
        pandas_objs["foo"] = pandas_objs["foo"].reindex(rearranged_index)
        ds_based_on_pandas = Dataset(pandas_objs, ds.coords, attrs=ds.attrs)
        del ds_based_on_pandas["x"]
        assert_equal(ds, ds_based_on_pandas)

    def test_constructor_pandas_single(self) -> None:
        das = [
            DataArray(np.random.rand(4), dims=["a"]),  # series
            DataArray(np.random.rand(4, 3), dims=["a", "b"]),  # df
        ]

        for a in das:
            pandas_obj = a.to_pandas()
            ds_based_on_pandas = Dataset(pandas_obj)  # type: ignore[arg-type]  # TODO: improve typing of __init__
            for dim in ds_based_on_pandas.data_vars:
                assert isinstance(dim, int)
                assert_array_equal(ds_based_on_pandas[dim], pandas_obj[dim])

    def test_constructor_compat(self) -> None:
        data = {"x": DataArray(0, coords={"y": 1}), "y": ("z", [1, 1, 1])}
        expected = Dataset({"x": 0}, {"y": ("z", [1, 1, 1])})
        actual = Dataset(data)
        assert_identical(expected, actual)

        data = {"y": ("z", [1, 1, 1]), "x": DataArray(0, coords={"y": 1})}
        actual = Dataset(data)
        assert_identical(expected, actual)

        original = Dataset(
            {"a": (("x", "y"), np.ones((2, 3)))},
            {"c": (("x", "y"), np.zeros((2, 3))), "x": [0, 1]},
        )
        expected = Dataset(
            {"a": ("x", np.ones(2)), "b": ("y", np.ones(3))},
            {"c": (("x", "y"), np.zeros((2, 3))), "x": [0, 1]},
        )

        actual = Dataset(
            {"a": original["a"][:, 0], "b": original["a"][0].drop_vars("x")}
        )
        assert_identical(expected, actual)

        data = {"x": DataArray(0, coords={"y": 3}), "y": ("z", [1, 1, 1])}
        with pytest.raises(MergeError):
            Dataset(data)

        data = {"x": DataArray(0, coords={"y": 1}), "y": [1, 1]}
        actual = Dataset(data)
        expected = Dataset({"x": 0}, {"y": [1, 1]})
        assert_identical(expected, actual)

    def test_constructor_with_coords(self) -> None:
        with pytest.raises(ValueError, match=r"found in both data_vars and"):
            Dataset({"a": ("x", [1])}, {"a": ("x", [1])})

        ds = Dataset({}, {"a": ("x", [1])})
        assert not ds.data_vars
        assert list(ds.coords.keys()) == ["a"]

        mindex = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=("level_1", "level_2")
        )
        with pytest.raises(ValueError, match=r"conflicting MultiIndex"):
            with pytest.warns(
                FutureWarning,
                match=".*`pandas.MultiIndex`.*no longer be implicitly promoted.*",
            ):
                Dataset({}, {"x": mindex, "y": mindex})
                Dataset({}, {"x": mindex, "level_1": range(4)})

    def test_constructor_no_default_index(self) -> None:
        # explicitly passing a Coordinates object skips the creation of default index
        ds = Dataset(coords=Coordinates({"x": [1, 2, 3]}, indexes={}))
        assert "x" in ds
        assert "x" not in ds.xindexes

    def test_constructor_multiindex(self) -> None:
        midx = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("one", "two"))
        coords = Coordinates.from_pandas_multiindex(midx, "x")

        ds = Dataset(coords=coords)
        assert_identical(ds, coords.to_dataset())

        with pytest.warns(
            FutureWarning,
            match=".*`pandas.MultiIndex`.*no longer be implicitly promoted.*",
        ):
            Dataset(data_vars={"x": midx})

        with pytest.warns(
            FutureWarning,
            match=".*`pandas.MultiIndex`.*no longer be implicitly promoted.*",
        ):
            Dataset(coords={"x": midx})

    def test_constructor_custom_index(self) -> None:
        class CustomIndex(Index): ...

        coords = Coordinates(
            coords={"x": ("x", [1, 2, 3])}, indexes={"x": CustomIndex()}
        )
        ds = Dataset(coords=coords)
        assert isinstance(ds.xindexes["x"], CustomIndex)

        # test coordinate variables copied
        assert ds.variables["x"] is not coords.variables["x"]

    @pytest.mark.filterwarnings("ignore:return type")
    def test_properties(self) -> None:
        ds = create_test_data()

        # dims / sizes
        # These exact types aren't public API, but this makes sure we don't
        # change them inadvertently:
        assert isinstance(ds.dims, utils.Frozen)
        # TODO change after deprecation cycle in GH #8500 is complete
        assert isinstance(ds.dims.mapping, dict)
        assert type(ds.dims.mapping) is dict
        with pytest.warns(
            FutureWarning,
            match=" To access a mapping from dimension names to lengths, please use `Dataset.sizes`",
        ):
            assert ds.dims == ds.sizes
        assert ds.sizes == {"dim1": 8, "dim2": 9, "dim3": 10, "time": 20}

        # dtypes
        assert isinstance(ds.dtypes, utils.Frozen)
        assert isinstance(ds.dtypes.mapping, dict)
        assert ds.dtypes == {
            "var1": np.dtype("float64"),
            "var2": np.dtype("float64"),
            "var3": np.dtype("float64"),
        }

        # data_vars
        assert list(ds) == list(ds.data_vars)
        assert list(ds.keys()) == list(ds.data_vars)
        assert "aasldfjalskdfj" not in ds.variables
        assert "dim1" in repr(ds.variables)
        assert len(ds) == 3
        assert bool(ds)

        assert list(ds.data_vars) == ["var1", "var2", "var3"]
        assert list(ds.data_vars.keys()) == ["var1", "var2", "var3"]
        assert "var1" in ds.data_vars
        assert "dim1" not in ds.data_vars
        assert "numbers" not in ds.data_vars
        assert len(ds.data_vars) == 3

        # xindexes
        assert set(ds.xindexes) == {"dim2", "dim3", "time"}
        assert len(ds.xindexes) == 3
        assert "dim2" in repr(ds.xindexes)
        assert all(isinstance(idx, Index) for idx in ds.xindexes.values())

        # indexes
        assert set(ds.indexes) == {"dim2", "dim3", "time"}
        assert len(ds.indexes) == 3
        assert "dim2" in repr(ds.indexes)
        assert all(isinstance(idx, pd.Index) for idx in ds.indexes.values())

        # coords
        assert list(ds.coords) == ["dim2", "dim3", "time", "numbers"]
        assert "dim2" in ds.coords
        assert "numbers" in ds.coords
        assert "var1" not in ds.coords
        assert "dim1" not in ds.coords
        assert len(ds.coords) == 4

        # nbytes
        assert (
            Dataset({"x": np.int64(1), "y": np.array([1, 2], dtype=np.float32)}).nbytes
            == 16
        )

    def test_warn_ds_dims_deprecation(self) -> None:
        # TODO remove after deprecation cycle in GH #8500 is complete
        ds = create_test_data()

        with pytest.warns(FutureWarning, match="return type"):
            ds.dims["dim1"]

        with pytest.warns(FutureWarning, match="return type"):
            ds.dims.keys()

        with pytest.warns(FutureWarning, match="return type"):
            ds.dims.values()

        with pytest.warns(FutureWarning, match="return type"):
            ds.dims.items()

        with assert_no_warnings():
            len(ds.dims)
            ds.dims.__iter__()
            _ = "dim1" in ds.dims

    def test_asarray(self) -> None:
        ds = Dataset({"x": 0})
        with pytest.raises(TypeError, match=r"cannot directly convert"):
            np.asarray(ds)

    def test_get_index(self) -> None:
        ds = Dataset({"foo": (("x", "y"), np.zeros((2, 3)))}, coords={"x": ["a", "b"]})
        assert ds.get_index("x").equals(pd.Index(["a", "b"]))
        assert ds.get_index("y").equals(pd.Index([0, 1, 2]))
        with pytest.raises(KeyError):
            ds.get_index("z")

    def test_attr_access(self) -> None:
        ds = Dataset(
            {"tmin": ("x", [42], {"units": "Celsius"})}, attrs={"title": "My test data"}
        )
        assert_identical(ds.tmin, ds["tmin"])
        assert_identical(ds.tmin.x, ds.x)

        assert ds.title == ds.attrs["title"]
        assert ds.tmin.units == ds["tmin"].attrs["units"]

        assert {"tmin", "title"} <= set(dir(ds))
        assert "units" in set(dir(ds.tmin))

        # should defer to variable of same name
        ds.attrs["tmin"] = -999
        assert ds.attrs["tmin"] == -999
        assert_identical(ds.tmin, ds["tmin"])

    def test_variable(self) -> None:
        a = Dataset()
        d = np.random.random((10, 3))
        a["foo"] = (("time", "x"), d)
        assert "foo" in a.variables
        assert "foo" in a
        a["bar"] = (("time", "x"), d)
        # order of creation is preserved
        assert list(a.variables) == ["foo", "bar"]
        assert_array_equal(a["foo"].values, d)
        # try to add variable with dim (10,3) with data that's (3,10)
        with pytest.raises(ValueError):
            a["qux"] = (("time", "x"), d.T)

    def test_modify_inplace(self) -> None:
        a = Dataset()
        vec = np.random.random((10,))
        attributes = {"foo": "bar"}
        a["x"] = ("x", vec, attributes)
        assert "x" in a.coords
        assert isinstance(a.coords["x"].to_index(), pd.Index)
        assert_identical(a.coords["x"].variable, a.variables["x"])
        b = Dataset()
        b["x"] = ("x", vec, attributes)
        assert_identical(a["x"], b["x"])
        assert a.sizes == b.sizes
        # this should work
        a["x"] = ("x", vec[:5])
        a["z"] = ("x", np.arange(5))
        with pytest.raises(ValueError):
            # now it shouldn't, since there is a conflicting length
            a["x"] = ("x", vec[:4])
        arr = np.random.random((10, 1))
        scal = np.array(0)
        with pytest.raises(ValueError):
            a["y"] = ("y", arr)
        with pytest.raises(ValueError):
            a["y"] = ("y", scal)
        assert "y" not in a.dims

    def test_coords_properties(self) -> None:
        # use int64 for repr consistency on windows
        data = Dataset(
            {
                "x": ("x", np.array([-1, -2], "int64")),
                "y": ("y", np.array([0, 1, 2], "int64")),
                "foo": (["x", "y"], np.random.randn(2, 3)),
            },
            {"a": ("x", np.array([4, 5], "int64")), "b": np.int64(-10)},
        )

        coords = data.coords
        assert isinstance(coords, DatasetCoordinates)

        # len
        assert len(coords) == 4

        # iter
        assert list(coords) == ["x", "y", "a", "b"]

        assert_identical(coords["x"].variable, data["x"].variable)
        assert_identical(coords["y"].variable, data["y"].variable)

        assert "x" in coords
        assert "a" in coords
        assert 0 not in coords
        assert "foo" not in coords

        with pytest.raises(KeyError):
            coords["foo"]
        with pytest.raises(KeyError):
            coords[0]

        # repr
        expected = dedent(
            """\
        Coordinates:
          * x        (x) int64 16B -1 -2
          * y        (y) int64 24B 0 1 2
            a        (x) int64 16B 4 5
            b        int64 8B -10"""
        )
        actual = repr(coords)
        assert expected == actual

        # dims
        assert coords.sizes == {"x": 2, "y": 3}

        # dtypes
        assert coords.dtypes == {
            "x": np.dtype("int64"),
            "y": np.dtype("int64"),
            "a": np.dtype("int64"),
            "b": np.dtype("int64"),
        }

    def test_coords_modify(self) -> None:
        data = Dataset(
            {
                "x": ("x", [-1, -2]),
                "y": ("y", [0, 1, 2]),
                "foo": (["x", "y"], np.random.randn(2, 3)),
            },
            {"a": ("x", [4, 5]), "b": -10},
        )

        actual = data.copy(deep=True)
        actual.coords["x"] = ("x", ["a", "b"])
        assert_array_equal(actual["x"], ["a", "b"])

        actual = data.copy(deep=True)
        actual.coords["z"] = ("z", ["a", "b"])
        assert_array_equal(actual["z"], ["a", "b"])

        actual = data.copy(deep=True)
        with pytest.raises(ValueError, match=r"conflicting dimension sizes"):
            actual.coords["x"] = ("x", [-1])
        assert_identical(actual, data)  # should not be modified

        actual = data.copy()
        del actual.coords["b"]
        expected = data.reset_coords("b", drop=True)
        assert_identical(expected, actual)

        with pytest.raises(KeyError):
            del data.coords["not_found"]

        with pytest.raises(KeyError):
            del data.coords["foo"]

        actual = data.copy(deep=True)
        actual.coords.update({"c": 11})
        expected = data.merge({"c": 11}).set_coords("c")
        assert_identical(expected, actual)

        # regression test for GH3746
        del actual.coords["x"]
        assert "x" not in actual.xindexes

    def test_update_index(self) -> None:
        actual = Dataset(coords={"x": [1, 2, 3]})
        actual["x"] = ["a", "b", "c"]
        assert actual.xindexes["x"].to_pandas_index().equals(pd.Index(["a", "b", "c"]))

    def test_coords_setitem_with_new_dimension(self) -> None:
        actual = Dataset()
        actual.coords["foo"] = ("x", [1, 2, 3])
        expected = Dataset(coords={"foo": ("x", [1, 2, 3])})
        assert_identical(expected, actual)

    def test_coords_setitem_multiindex(self) -> None:
        data = create_test_multiindex()
        with pytest.raises(ValueError, match=r"cannot drop or update.*corrupt.*index "):
            data.coords["level_1"] = range(4)

    def test_coords_set(self) -> None:
        one_coord = Dataset({"x": ("x", [0]), "yy": ("x", [1]), "zzz": ("x", [2])})
        two_coords = Dataset({"zzz": ("x", [2])}, {"x": ("x", [0]), "yy": ("x", [1])})
        all_coords = Dataset(
            coords={"x": ("x", [0]), "yy": ("x", [1]), "zzz": ("x", [2])}
        )

        actual = one_coord.set_coords("x")
        assert_identical(one_coord, actual)
        actual = one_coord.set_coords(["x"])
        assert_identical(one_coord, actual)

        actual = one_coord.set_coords("yy")
        assert_identical(two_coords, actual)

        actual = one_coord.set_coords(["yy", "zzz"])
        assert_identical(all_coords, actual)

        actual = one_coord.reset_coords()
        assert_identical(one_coord, actual)
        actual = two_coords.reset_coords()
        assert_identical(one_coord, actual)
        actual = all_coords.reset_coords()
        assert_identical(one_coord, actual)

        actual = all_coords.reset_coords(["yy", "zzz"])
        assert_identical(one_coord, actual)
        actual = all_coords.reset_coords("zzz")
        assert_identical(two_coords, actual)

        with pytest.raises(ValueError, match=r"cannot remove index"):
            one_coord.reset_coords("x")

        actual = all_coords.reset_coords("zzz", drop=True)
        expected = all_coords.drop_vars("zzz")
        assert_identical(expected, actual)
        expected = two_coords.drop_vars("zzz")
        assert_identical(expected, actual)

    def test_coords_to_dataset(self) -> None:
        orig = Dataset({"foo": ("y", [-1, 0, 1])}, {"x": 10, "y": [2, 3, 4]})
        expected = Dataset(coords={"x": 10, "y": [2, 3, 4]})
        actual = orig.coords.to_dataset()
        assert_identical(expected, actual)

    def test_coords_merge(self) -> None:
        orig_coords = Dataset(coords={"a": ("x", [1, 2]), "x": [0, 1]}).coords
        other_coords = Dataset(coords={"b": ("x", ["a", "b"]), "x": [0, 1]}).coords
        expected = Dataset(
            coords={"a": ("x", [1, 2]), "b": ("x", ["a", "b"]), "x": [0, 1]}
        )
        actual = orig_coords.merge(other_coords)
        assert_identical(expected, actual)
        actual = other_coords.merge(orig_coords)
        assert_identical(expected, actual)

        other_coords = Dataset(coords={"x": ("x", ["a"])}).coords
        with pytest.raises(MergeError):
            orig_coords.merge(other_coords)
        other_coords = Dataset(coords={"x": ("x", ["a", "b"])}).coords
        with pytest.raises(MergeError):
            orig_coords.merge(other_coords)
        other_coords = Dataset(coords={"x": ("x", ["a", "b", "c"])}).coords
        with pytest.raises(MergeError):
            orig_coords.merge(other_coords)

        other_coords = Dataset(coords={"a": ("x", [8, 9])}).coords
        expected = Dataset(coords={"x": range(2)})
        actual = orig_coords.merge(other_coords)
        assert_identical(expected, actual)
        actual = other_coords.merge(orig_coords)
        assert_identical(expected, actual)

        other_coords = Dataset(coords={"x": np.nan}).coords
        actual = orig_coords.merge(other_coords)
        assert_identical(orig_coords.to_dataset(), actual)
        actual = other_coords.merge(orig_coords)
        assert_identical(orig_coords.to_dataset(), actual)

    def test_coords_merge_mismatched_shape(self) -> None:
        orig_coords = Dataset(coords={"a": ("x", [1, 1])}).coords
        other_coords = Dataset(coords={"a": 1}).coords
        expected = orig_coords.to_dataset()
        actual = orig_coords.merge(other_coords)
        assert_identical(expected, actual)

        other_coords = Dataset(coords={"a": ("y", [1])}).coords
        expected = Dataset(coords={"a": (["x", "y"], [[1], [1]])})
        actual = orig_coords.merge(other_coords)
        assert_identical(expected, actual)

        actual = other_coords.merge(orig_coords)
        assert_identical(expected.transpose(), actual)

        orig_coords = Dataset(coords={"a": ("x", [np.nan])}).coords
        other_coords = Dataset(coords={"a": np.nan}).coords
        expected = orig_coords.to_dataset()
        actual = orig_coords.merge(other_coords)
        assert_identical(expected, actual)

    def test_data_vars_properties(self) -> None:
        ds = Dataset()
        ds["foo"] = (("x",), [1.0])
        ds["bar"] = 2.0

        # iter
        assert set(ds.data_vars) == {"foo", "bar"}
        assert "foo" in ds.data_vars
        assert "x" not in ds.data_vars
        assert_identical(ds["foo"], ds.data_vars["foo"])

        # repr
        expected = dedent(
            """\
        Data variables:
            foo      (x) float64 8B 1.0
            bar      float64 8B 2.0"""
        )
        actual = repr(ds.data_vars)
        assert expected == actual

        # dtypes
        assert ds.data_vars.dtypes == {
            "foo": np.dtype("float64"),
            "bar": np.dtype("float64"),
        }

        # len
        ds.coords["x"] = [1]
        assert len(ds.data_vars) == 2

        # https://github.com/pydata/xarray/issues/7588
        with pytest.raises(
            AssertionError, match="something is wrong with Dataset._coord_names"
        ):
            ds._coord_names = {"w", "x", "y", "z"}
            len(ds.data_vars)

    def test_equals_and_identical(self) -> None:
        data = create_test_data(seed=42)
        assert data.equals(data)
        assert data.identical(data)

        data2 = create_test_data(seed=42)
        data2.attrs["foobar"] = "baz"
        assert data.equals(data2)
        assert not data.identical(data2)

        del data2["time"]
        assert not data.equals(data2)

        data = create_test_data(seed=42).rename({"var1": None})
        assert data.equals(data)
        assert data.identical(data)

        data2 = data.reset_coords()
        assert not data2.equals(data)
        assert not data2.identical(data)

    def test_equals_failures(self) -> None:
        data = create_test_data()
        assert not data.equals("foo")  # type: ignore[arg-type]
        assert not data.identical(123)  # type: ignore[arg-type]
        assert not data.broadcast_equals({1: 2})  # type: ignore[arg-type]

    def test_broadcast_equals(self) -> None:
        data1 = Dataset(coords={"x": 0})
        data2 = Dataset(coords={"x": [0]})
        assert data1.broadcast_equals(data2)
        assert not data1.equals(data2)
        assert not data1.identical(data2)

    def test_attrs(self) -> None:
        data = create_test_data(seed=42)
        data.attrs = {"foobar": "baz"}
        assert data.attrs["foobar"], "baz"
        assert isinstance(data.attrs, dict)

    def test_chunks_does_not_load_data(self) -> None:
        # regression test for GH6538
        store = InaccessibleVariableDataStore()
        create_test_data().dump_to_store(store)
        ds = open_dataset(store)
        assert ds.chunks == {}

    @requires_dask
    @pytest.mark.parametrize(
        "use_cftime,calendar",
        [
            (False, "standard"),
            (pytest.param(True, marks=pytest.mark.skipif(not has_cftime)), "standard"),
            (pytest.param(True, marks=pytest.mark.skipif(not has_cftime)), "noleap"),
            (pytest.param(True, marks=pytest.mark.skipif(not has_cftime)), "360_day"),
        ],
    )
    def test_chunk_by_season_resampler(self, use_cftime: bool, calendar: str) -> None:
        import dask.array

        N = 365 + 365  # 2 years - 1 day
        time = xr.date_range(
            "2000-01-01", periods=N, freq="D", use_cftime=use_cftime, calendar=calendar
        )

        ds = Dataset(
            {
                "pr": ("time", dask.array.random.random((N), chunks=(20))),
                "pr2d": (("x", "time"), dask.array.random.random((10, N), chunks=(20))),
                "ones": ("time", np.ones((N,))),
            },
            coords={"time": time},
        )

        # Standard seasons
        rechunked = ds.chunk(
            {"x": 2, "time": SeasonResampler(["DJF", "MAM", "JJA", "SON"])}
        )
        assert rechunked.chunksizes["x"] == (2,) * 5
        assert len(rechunked.chunksizes["time"]) == 9
        assert rechunked.chunksizes["x"] == (2,) * 5
        assert sum(rechunked.chunksizes["time"]) == ds.sizes["time"]

        if calendar == "standard":
            assert rechunked.chunksizes["time"] == (60, 92, 92, 91, 90, 92, 92, 91, 30)
        elif calendar == "noleap":
            assert rechunked.chunksizes["time"] == (59, 92, 92, 91, 90, 92, 92, 91, 31)
        elif calendar == "360_day":
            assert rechunked.chunksizes["time"] == (60, 90, 90, 90, 90, 90, 90, 90, 40)
        else:
            raise AssertionError("unreachable")

        # Custom seasons
        rechunked = ds.chunk(
            {"x": 2, "time": SeasonResampler(["DJFM", "AM", "JJA", "SON"])}
        )
        assert len(rechunked.chunksizes["time"]) == 9
        assert sum(rechunked.chunksizes["time"]) == ds.sizes["time"]
        assert rechunked.chunksizes["x"] == (2,) * 5

        if calendar == "standard":
            assert rechunked.chunksizes["time"] == (91, 61, 92, 91, 121, 61, 92, 91, 30)
        elif calendar == "noleap":
            assert rechunked.chunksizes["time"] == (90, 61, 92, 91, 121, 61, 92, 91, 31)
        elif calendar == "360_day":
            assert rechunked.chunksizes["time"] == (90, 60, 90, 90, 120, 60, 90, 90, 40)
        else:
            raise AssertionError("unreachable")

        # Test that drop_incomplete doesn't affect chunking
        rechunked_drop_true = ds.chunk(
            time=SeasonResampler(["DJF", "MAM", "JJA", "SON"], drop_incomplete=True)
        )
        rechunked_drop_false = ds.chunk(
            time=SeasonResampler(["DJF", "MAM", "JJA", "SON"], drop_incomplete=False)
        )
        assert (
            rechunked_drop_true.chunksizes["time"]
            == rechunked_drop_false.chunksizes["time"]
        )

    @requires_dask
    def test_chunk_by_season_resampler_errors(self):
        """Test error handling for SeasonResampler chunking."""
        # Test error on missing season (should fail with incomplete seasons)
        ds = Dataset(
            {"x": ("time", np.arange(12))},
            coords={"time": pd.date_range("2000-01-01", periods=12, freq="MS")},
        )
        with pytest.raises(ValueError, match="does not cover all 12 months"):
            ds.chunk(time=SeasonResampler(["DJF", "MAM", "SON"]))

        ds = Dataset({"foo": ("x", [1, 2, 3])})
        # Test error on virtual variable
        with pytest.raises(ValueError, match="virtual variable"):
            ds.chunk(x=SeasonResampler(["DJF", "MAM", "JJA", "SON"]))

        # Test error on non-datetime variable
        ds["x"] = ("x", [1, 2, 3])
        with pytest.raises(ValueError, match="datetime variables"):
            ds.chunk(x=SeasonResampler(["DJF", "MAM", "JJA", "SON"]))

        # Test successful case with 1D datetime variable
        ds["x"] = ("x", xr.date_range("2001-01-01", periods=3, freq="D"))
        # This should work
        result = ds.chunk(x=SeasonResampler(["DJF", "MAM", "JJA", "SON"]))
        assert result.chunks is not None

        # Test error on missing season (should fail with incomplete seasons)
        with pytest.raises(ValueError):
            ds.chunk(x=SeasonResampler(["DJF", "MAM", "SON"]))

    @requires_dask
    def test_chunk(self) -> None:
        data = create_test_data()
        for v in data.variables.values():
            assert isinstance(v.data, np.ndarray)
        assert data.chunks == {}

        reblocked = data.chunk()
        for k, v in reblocked.variables.items():
            if k in reblocked.dims:
                assert isinstance(v.data, np.ndarray)
            else:
                assert isinstance(v.data, da.Array)

        expected_chunks: dict[Hashable, tuple[int, ...]] = {
            "dim1": (8,),
            "dim2": (9,),
            "dim3": (10,),
        }
        assert reblocked.chunks == expected_chunks

        # test kwargs form of chunks
        assert data.chunk(expected_chunks).chunks == expected_chunks

        def get_dask_names(ds):
            return {k: v.data.name for k, v in ds.items()}

        orig_dask_names = get_dask_names(reblocked)

        reblocked = data.chunk({"time": 5, "dim1": 5, "dim2": 5, "dim3": 5})
        # time is not a dim in any of the data_vars, so it
        # doesn't get chunked
        expected_chunks = {"dim1": (5, 3), "dim2": (5, 4), "dim3": (5, 5)}
        assert reblocked.chunks == expected_chunks

        # make sure dask names change when rechunking by different amounts
        # regression test for GH3350
        new_dask_names = get_dask_names(reblocked)
        for k, v in new_dask_names.items():
            assert v != orig_dask_names[k]

        reblocked = data.chunk(expected_chunks)
        assert reblocked.chunks == expected_chunks

        # reblock on already blocked data
        orig_dask_names = get_dask_names(reblocked)
        reblocked = reblocked.chunk(expected_chunks)
        new_dask_names = get_dask_names(reblocked)
        assert reblocked.chunks == expected_chunks
        assert_identical(reblocked, data)
        # rechunking with same chunk sizes should not change names
        for k, v in new_dask_names.items():
            assert v == orig_dask_names[k]

        with pytest.raises(
            ValueError,
            match=re.escape(
                "chunks keys ('foo',) not found in data dimensions ('dim2', 'dim3', 'time', 'dim1')"
            ),
        ):
            data.chunk({"foo": 10})

    @requires_dask
    @pytest.mark.parametrize(
        "calendar",
        (
            "standard",
            pytest.param(
                "gregorian",
                marks=pytest.mark.skipif(not has_cftime, reason="needs cftime"),
            ),
        ),
    )
    @pytest.mark.parametrize("freq", ["D", "W", "5ME", "YE"])
    @pytest.mark.parametrize("add_gap", [True, False])
    def test_chunk_by_frequency(self, freq: str, calendar: str, add_gap: bool) -> None:
        import dask.array

        N = 365 * 2
        ΔN = 28  # noqa: PLC2401
        time = xr.date_range(
            "2001-01-01", periods=N + ΔN, freq="D", calendar=calendar
        ).to_numpy(copy=True)
        if add_gap:
            # introduce an empty bin
            time[31 : 31 + ΔN] = np.datetime64("NaT")
            time = time[~np.isnat(time)]
        else:
            time = time[:N]

        ds = Dataset(
            {
                "pr": ("time", dask.array.random.random((N), chunks=(20))),
                "pr2d": (("x", "time"), dask.array.random.random((10, N), chunks=(20))),
                "ones": ("time", np.ones((N,))),
            },
            coords={"time": time},
        )
        rechunked = ds.chunk(x=2, time=TimeResampler(freq))
        expected = tuple(
            ds.ones.resample(time=freq).sum().dropna("time").astype(int).data.tolist()
        )
        assert rechunked.chunksizes["time"] == expected
        assert rechunked.chunksizes["x"] == (2,) * 5

        rechunked = ds.chunk({"x": 2, "time": TimeResampler(freq)})
        assert rechunked.chunksizes["time"] == expected
        assert rechunked.chunksizes["x"] == (2,) * 5

    def test_chunk_by_frequency_errors(self):
        ds = Dataset({"foo": ("x", [1, 2, 3])})
        with pytest.raises(ValueError, match="virtual variable"):
            ds.chunk(x=TimeResampler("YE"))
        ds["x"] = ("x", [1, 2, 3])
        with pytest.raises(ValueError, match="datetime variables"):
            ds.chunk(x=TimeResampler("YE"))
        ds["x"] = ("x", xr.date_range("2001-01-01", periods=3, freq="D"))
        with pytest.raises(ValueError, match="Invalid frequency"):
            ds.chunk(x=TimeResampler("foo"))

    @requires_dask
    def test_dask_is_lazy(self) -> None:
        store = InaccessibleVariableDataStore()
        create_test_data().dump_to_store(store)
        ds = open_dataset(store).chunk()

        with pytest.raises(UnexpectedDataAccess):
            ds.load()
        with pytest.raises(UnexpectedDataAccess):
            _ = ds["var1"].values

        # these should not raise UnexpectedDataAccess:
        _ = ds.var1.data
        ds.isel(time=10)
        ds.isel(time=slice(10), dim1=[0]).isel(dim1=0, dim2=-1)
        ds.transpose()
        ds.mean()
        ds.fillna(0)
        ds.rename({"dim1": "foobar"})
        ds.set_coords("var1")
        ds.drop_vars("var1")

    def test_isel(self) -> None:
        data = create_test_data()
        slicers: dict[Hashable, slice] = {
            "dim1": slice(None, None, 2),
            "dim2": slice(0, 2),
        }
        ret = data.isel(slicers)

        # Verify that only the specified dimension was altered
        assert list(data.dims) == list(ret.dims)
        for d in data.dims:
            if d in slicers:
                assert ret.sizes[d] == np.arange(data.sizes[d])[slicers[d]].size
            else:
                assert data.sizes[d] == ret.sizes[d]
        # Verify that the data is what we expect
        for v in data.variables:
            assert data[v].dims == ret[v].dims
            assert data[v].attrs == ret[v].attrs
            slice_list = [slice(None)] * data[v].values.ndim
            for d, s in slicers.items():
                if d in data[v].dims:
                    inds = np.nonzero(np.array(data[v].dims) == d)[0]
                    for ind in inds:
                        slice_list[ind] = s
            expected = data[v].values[tuple(slice_list)]
            actual = ret[v].values
            np.testing.assert_array_equal(expected, actual)

        with pytest.raises(ValueError):
            data.isel(not_a_dim=slice(0, 2))
        with pytest.raises(
            ValueError,
            match=r"Dimensions {'not_a_dim'} do not exist. Expected "
            r"one or more of "
            r"[\w\W]*'dim\d'[\w\W]*'dim\d'[\w\W]*'time'[\w\W]*'dim\d'[\w\W]*",
        ):
            data.isel(not_a_dim=slice(0, 2))
        with pytest.warns(
            UserWarning,
            match=r"Dimensions {'not_a_dim'} do not exist. "
            r"Expected one or more of "
            r"[\w\W]*'dim\d'[\w\W]*'dim\d'[\w\W]*'time'[\w\W]*'dim\d'[\w\W]*",
        ):
            data.isel(not_a_dim=slice(0, 2), missing_dims="warn")
        assert_identical(data, data.isel(not_a_dim=slice(0, 2), missing_dims="ignore"))

        ret = data.isel(dim1=0)
        assert {"time": 20, "dim2": 9, "dim3": 10} == ret.sizes
        assert set(data.data_vars) == set(ret.data_vars)
        assert set(data.coords) == set(ret.coords)
        assert set(data.xindexes) == set(ret.xindexes)

        ret = data.isel(time=slice(2), dim1=0, dim2=slice(5))
        assert {"time": 2, "dim2": 5, "dim3": 10} == ret.sizes
        assert set(data.data_vars) == set(ret.data_vars)
        assert set(data.coords) == set(ret.coords)
        assert set(data.xindexes) == set(ret.xindexes)

        ret = data.isel(time=0, dim1=0, dim2=slice(5))
        assert {"dim2": 5, "dim3": 10} == ret.sizes
        assert set(data.data_vars) == set(ret.data_vars)
        assert set(data.coords) == set(ret.coords)
        assert set(data.xindexes) == set(list(ret.xindexes) + ["time"])

    def test_isel_fancy(self) -> None:
        # isel with fancy indexing.
        data = create_test_data()

        pdim1 = [1, 2, 3]
        pdim2 = [4, 5, 1]
        pdim3 = [1, 2, 3]
        actual = data.isel(
            dim1=(("test_coord",), pdim1),
            dim2=(("test_coord",), pdim2),
            dim3=(("test_coord",), pdim3),
        )
        assert "test_coord" in actual.dims
        assert actual.coords["test_coord"].shape == (len(pdim1),)

        # Should work with DataArray
        actual = data.isel(
            dim1=DataArray(pdim1, dims="test_coord"),
            dim2=(("test_coord",), pdim2),
            dim3=(("test_coord",), pdim3),
        )
        assert "test_coord" in actual.dims
        assert actual.coords["test_coord"].shape == (len(pdim1),)
        expected = data.isel(
            dim1=(("test_coord",), pdim1),
            dim2=(("test_coord",), pdim2),
            dim3=(("test_coord",), pdim3),
        )
        assert_identical(actual, expected)

        # DataArray with coordinate
        idx1 = DataArray(pdim1, dims=["a"], coords={"a": np.random.randn(3)})
        idx2 = DataArray(pdim2, dims=["b"], coords={"b": np.random.randn(3)})
        idx3 = DataArray(pdim3, dims=["c"], coords={"c": np.random.randn(3)})
        # Should work with DataArray
        actual = data.isel(dim1=idx1, dim2=idx2, dim3=idx3)
        assert "a" in actual.dims
        assert "b" in actual.dims
        assert "c" in actual.dims
        assert "time" in actual.coords
        assert "dim2" in actual.coords
        assert "dim3" in actual.coords
        expected = data.isel(
            dim1=(("a",), pdim1), dim2=(("b",), pdim2), dim3=(("c",), pdim3)
        )
        expected = expected.assign_coords(a=idx1["a"], b=idx2["b"], c=idx3["c"])
        assert_identical(actual, expected)

        idx1 = DataArray(pdim1, dims=["a"], coords={"a": np.random.randn(3)})
        idx2 = DataArray(pdim2, dims=["a"])
        idx3 = DataArray(pdim3, dims=["a"])
        # Should work with DataArray
        actual = data.isel(dim1=idx1, dim2=idx2, dim3=idx3)
        assert "a" in actual.dims
        assert "time" in actual.coords
        assert "dim2" in actual.coords
        assert "dim3" in actual.coords
        expected = data.isel(
            dim1=(("a",), pdim1), dim2=(("a",), pdim2), dim3=(("a",), pdim3)
        )
        expected = expected.assign_coords(a=idx1["a"])
        assert_identical(actual, expected)

        actual = data.isel(dim1=(("points",), pdim1), dim2=(("points",), pdim2))
        assert "points" in actual.dims
        assert "dim3" in actual.dims
        assert "dim3" not in actual.data_vars
        np.testing.assert_array_equal(data["dim2"][pdim2], actual["dim2"])

        # test that the order of the indexers doesn't matter
        assert_identical(
            data.isel(dim1=(("points",), pdim1), dim2=(("points",), pdim2)),
            data.isel(dim2=(("points",), pdim2), dim1=(("points",), pdim1)),
        )
        # make sure we're raising errors in the right places
        with pytest.raises(IndexError, match=r"Dimensions of indexers mismatch"):
            data.isel(dim1=(("points",), [1, 2]), dim2=(("points",), [1, 2, 3]))
        with pytest.raises(TypeError, match=r"cannot use a Dataset"):
            data.isel(dim1=Dataset({"points": [1, 2]}))

        # test to be sure we keep around variables that were not indexed
        ds = Dataset({"x": [1, 2, 3, 4], "y": 0})
        actual = ds.isel(x=(("points",), [0, 1, 2]))
        assert_identical(ds["y"], actual["y"])

        # tests using index or DataArray as indexers
        stations = Dataset()
        stations["station"] = (("station",), ["A", "B", "C"])
        stations["dim1s"] = (("station",), [1, 2, 3])
        stations["dim2s"] = (("station",), [4, 5, 1])

        actual = data.isel(dim1=stations["dim1s"], dim2=stations["dim2s"])
        assert "station" in actual.coords
        assert "station" in actual.dims
        assert_identical(actual["station"].drop_vars(["dim2"]), stations["station"])

        with pytest.raises(ValueError, match=r"conflicting values/indexes on "):
            data.isel(
                dim1=DataArray(
                    [0, 1, 2], dims="station", coords={"station": [0, 1, 2]}
                ),
                dim2=DataArray(
                    [0, 1, 2], dims="station", coords={"station": [0, 1, 3]}
                ),
            )

        # multi-dimensional selection
        stations = Dataset()
        stations["a"] = (("a",), ["A", "B", "C"])
        stations["b"] = (("b",), [0, 1])
        stations["dim1s"] = (("a", "b"), [[1, 2], [2, 3], [3, 4]])
        stations["dim2s"] = (("a",), [4, 5, 1])
        actual = data.isel(dim1=stations["dim1s"], dim2=stations["dim2s"])
        assert "a" in actual.coords
        assert "a" in actual.dims
        assert "b" in actual.coords
        assert "b" in actual.dims
        assert "dim2" in actual.coords
        assert "a" in actual["dim2"].dims

        assert_identical(actual["a"].drop_vars(["dim2"]), stations["a"])
        assert_identical(actual["b"], stations["b"])
        expected_var1 = data["var1"].variable[
            stations["dim1s"].variable, stations["dim2s"].variable
        ]
        expected_var2 = data["var2"].variable[
            stations["dim1s"].variable, stations["dim2s"].variable
        ]
        expected_var3 = data["var3"].variable[slice(None), stations["dim1s"].variable]
        assert_equal(actual["a"].drop_vars("dim2"), stations["a"])
        assert_array_equal(actual["var1"], expected_var1)
        assert_array_equal(actual["var2"], expected_var2)
        assert_array_equal(actual["var3"], expected_var3)

        # test that drop works
        ds = xr.Dataset({"a": (("x",), [1, 2, 3])}, coords={"b": (("x",), [5, 6, 7])})

        actual = ds.isel({"x": 1}, drop=False)
        expected = xr.Dataset({"a": 2}, coords={"b": 6})
        assert_identical(actual, expected)

        actual = ds.isel({"x": 1}, drop=True)
        expected = xr.Dataset({"a": 2})
        assert_identical(actual, expected)

        actual = ds.isel({"x": DataArray(1)}, drop=False)
        expected = xr.Dataset({"a": 2}, coords={"b": 6})
        assert_identical(actual, expected)

        actual = ds.isel({"x": DataArray(1)}, drop=True)
        expected = xr.Dataset({"a": 2})
        assert_identical(actual, expected)

    def test_isel_dataarray(self) -> None:
        """Test for indexing by DataArray"""
        data = create_test_data()
        # indexing with DataArray with same-name coordinates.
        indexing_da = DataArray(
            np.arange(1, 4), dims=["dim1"], coords={"dim1": np.random.randn(3)}
        )
        actual = data.isel(dim1=indexing_da)
        assert_identical(indexing_da["dim1"], actual["dim1"])
        assert_identical(data["dim2"], actual["dim2"])

        # Conflict in the dimension coordinate
        indexing_da = DataArray(
            np.arange(1, 4), dims=["dim2"], coords={"dim2": np.random.randn(3)}
        )
        with pytest.raises(IndexError, match=r"dimension coordinate 'dim2'"):
            data.isel(dim2=indexing_da)
        # Also the case for DataArray
        with pytest.raises(IndexError, match=r"dimension coordinate 'dim2'"):
            data["var2"].isel(dim2=indexing_da)
        with pytest.raises(IndexError, match=r"dimension coordinate 'dim2'"):
            data["dim2"].isel(dim2=indexing_da)

        # same name coordinate which does not conflict
        indexing_da = DataArray(
            np.arange(1, 4), dims=["dim2"], coords={"dim2": data["dim2"].values[1:4]}
        )
        actual = data.isel(dim2=indexing_da)
        assert_identical(actual["dim2"], indexing_da["dim2"])

        # Silently drop conflicted (non-dimensional) coordinate of indexer
        indexing_da = DataArray(
            np.arange(1, 4),
            dims=["dim2"],
            coords={
                "dim2": data["dim2"].values[1:4],
                "numbers": ("dim2", np.arange(2, 5)),
            },
        )
        actual = data.isel(dim2=indexing_da)
        assert_identical(actual["numbers"], data["numbers"])

        # boolean data array with coordinate with the same name
        indexing_da = DataArray(
            np.arange(1, 10), dims=["dim2"], coords={"dim2": data["dim2"].values}
        )
        indexing_da = indexing_da < 3
        actual = data.isel(dim2=indexing_da)
        assert_identical(actual["dim2"], data["dim2"][:2])

        # boolean data array with non-dimensioncoordinate
        indexing_da = DataArray(
            np.arange(1, 10),
            dims=["dim2"],
            coords={
                "dim2": data["dim2"].values,
                "non_dim": (("dim2",), np.random.randn(9)),
                "non_dim2": 0,
            },
        )
        indexing_da = indexing_da < 3
        actual = data.isel(dim2=indexing_da)
        assert_identical(
            actual["dim2"].drop_vars("non_dim").drop_vars("non_dim2"), data["dim2"][:2]
        )
        assert_identical(actual["non_dim"], indexing_da["non_dim"][:2])
        assert_identical(actual["non_dim2"], indexing_da["non_dim2"])

        # non-dimension coordinate will be also attached
        indexing_da = DataArray(
            np.arange(1, 4),
            dims=["dim2"],
            coords={"non_dim": (("dim2",), np.random.randn(3))},
        )
        actual = data.isel(dim2=indexing_da)
        assert "non_dim" in actual
        assert "non_dim" in actual.coords

        # Index by a scalar DataArray
        indexing_da = DataArray(3, dims=[], coords={"station": 2})
        actual = data.isel(dim2=indexing_da)
        assert "station" in actual
        actual = data.isel(dim2=indexing_da["station"])
        assert "station" in actual

        # indexer generated from coordinates
        indexing_ds = Dataset({}, coords={"dim2": [0, 1, 2]})
        with pytest.raises(IndexError, match=r"dimension coordinate 'dim2'"):
            actual = data.isel(dim2=indexing_ds["dim2"])

    def test_isel_fancy_convert_index_variable(self) -> None:
        # select index variable "x" with a DataArray of dim "z"
        # -> drop index and convert index variable to base variable
        ds = xr.Dataset({"foo": ("x", [1, 2, 3])}, coords={"x": [0, 1, 2]})
        idxr = xr.DataArray([1], dims="z", name="x")
        actual = ds.isel(x=idxr)
        assert "x" not in actual.xindexes
        assert not isinstance(actual.x.variable, IndexVariable)

    def test_isel_multicoord_index(self) -> None:
        # regression test https://github.com/pydata/xarray/issues/10063
        # isel on a multi-coordinate index should return a unique index associated
        # to each coordinate
        coords = xr.Coordinates(coords={"x": [0, 1], "y": [1, 2]}, indexes={})
        ds = xr.Dataset(coords=coords).set_xindex(["x", "y"], XYIndex)

        ds2 = ds.isel(x=slice(None), y=slice(None))
        assert ds2.xindexes["x"] is ds2.xindexes["y"]

    def test_sel(self) -> None:
        data = create_test_data()
        int_slicers = {"dim1": slice(None, None, 2), "dim2": slice(2), "dim3": slice(3)}
        loc_slicers = {
            "dim1": slice(None, None, 2),
            "dim2": slice(0, 0.5),
            "dim3": slice("a", "c"),
        }
        assert_equal(data.isel(int_slicers), data.sel(loc_slicers))
        data["time"] = ("time", pd.date_range("2000-01-01", periods=20))
        assert_equal(data.isel(time=0), data.sel(time="2000-01-01"))
        assert_equal(
            data.isel(time=slice(10)), data.sel(time=slice("2000-01-01", "2000-01-10"))
        )
        assert_equal(data, data.sel(time=slice("1999", "2005")))
        times = pd.date_range("2000-01-01", periods=3)
        assert_equal(data.isel(time=slice(3)), data.sel(time=times))
        assert_equal(
            data.isel(time=slice(3)), data.sel(time=(data["time.dayofyear"] <= 3))
        )

        td = pd.to_timedelta(np.arange(3), unit="days")
        data = Dataset({"x": ("td", np.arange(3)), "td": td})
        assert_equal(data, data.sel(td=td))
        assert_equal(data, data.sel(td=slice("3 days")))
        assert_equal(data.isel(td=0), data.sel(td=pd.Timedelta("0 days")))
        assert_equal(data.isel(td=0), data.sel(td=pd.Timedelta("0h")))
        assert_equal(data.isel(td=slice(1, 3)), data.sel(td=slice("1 days", "2 days")))

    def test_sel_dataarray(self) -> None:
        data = create_test_data()

        ind = DataArray([0.0, 0.5, 1.0], dims=["dim2"])
        actual = data.sel(dim2=ind)
        assert_equal(actual, data.isel(dim2=[0, 1, 2]))

        # with different dimension
        ind = DataArray([0.0, 0.5, 1.0], dims=["new_dim"])
        actual = data.sel(dim2=ind)
        expected = data.isel(dim2=Variable("new_dim", [0, 1, 2]))
        assert "new_dim" in actual.dims
        assert_equal(actual, expected)

        # Multi-dimensional
        ind = DataArray([[0.0], [0.5], [1.0]], dims=["new_dim", "new_dim2"])
        actual = data.sel(dim2=ind)
        expected = data.isel(dim2=Variable(("new_dim", "new_dim2"), [[0], [1], [2]]))
        assert "new_dim" in actual.dims
        assert "new_dim2" in actual.dims
        assert_equal(actual, expected)

        # with coordinate
        ind = DataArray(
            [0.0, 0.5, 1.0], dims=["new_dim"], coords={"new_dim": ["a", "b", "c"]}
        )
        actual = data.sel(dim2=ind)
        expected = data.isel(dim2=[0, 1, 2]).rename({"dim2": "new_dim"})
        assert "new_dim" in actual.dims
        assert "new_dim" in actual.coords
        assert_equal(
            actual.drop_vars("new_dim").drop_vars("dim2"), expected.drop_vars("new_dim")
        )
        assert_equal(actual["new_dim"].drop_vars("dim2"), ind["new_dim"])

        # with conflicted coordinate (silently ignored)
        ind = DataArray(
            [0.0, 0.5, 1.0], dims=["dim2"], coords={"dim2": ["a", "b", "c"]}
        )
        actual = data.sel(dim2=ind)
        expected = data.isel(dim2=[0, 1, 2])
        assert_equal(actual, expected)

        # with conflicted coordinate (silently ignored)
        ind = DataArray(
            [0.0, 0.5, 1.0],
            dims=["new_dim"],
            coords={"new_dim": ["a", "b", "c"], "dim2": 3},
        )
        actual = data.sel(dim2=ind)
        assert_equal(
            actual["new_dim"].drop_vars("dim2"), ind["new_dim"].drop_vars("dim2")
        )
        expected = data.isel(dim2=[0, 1, 2])
        expected["dim2"] = (("new_dim"), expected["dim2"].values)
        assert_equal(actual["dim2"].drop_vars("new_dim"), expected["dim2"])
        assert actual["var1"].dims == ("dim1", "new_dim")

        # with non-dimensional coordinate
        ind = DataArray(
            [0.0, 0.5, 1.0],
            dims=["dim2"],
            coords={
                "dim2": ["a", "b", "c"],
                "numbers": ("dim2", [0, 1, 2]),
                "new_dim": ("dim2", [1.1, 1.2, 1.3]),
            },
        )
        actual = data.sel(dim2=ind)
        expected = data.isel(dim2=[0, 1, 2])
        assert_equal(actual.drop_vars("new_dim"), expected)
        assert np.allclose(actual["new_dim"].values, ind["new_dim"].values)

    def test_sel_dataarray_mindex(self) -> None:
        midx = pd.MultiIndex.from_product([list("abc"), [0, 1]], names=("one", "two"))
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        midx_coords["y"] = range(3)

        mds = xr.Dataset(
            {"var": (("x", "y"), np.random.rand(6, 3))}, coords=midx_coords
        )

        actual_isel = mds.isel(x=xr.DataArray(np.arange(3), dims="x"))
        actual_sel = mds.sel(x=DataArray(midx[:3], dims="x"))
        assert actual_isel["x"].dims == ("x",)
        assert actual_sel["x"].dims == ("x",)
        assert_identical(actual_isel, actual_sel)

        actual_isel = mds.isel(x=xr.DataArray(np.arange(3), dims="z"))
        actual_sel = mds.sel(x=Variable("z", midx[:3]))
        assert actual_isel["x"].dims == ("z",)
        assert actual_sel["x"].dims == ("z",)
        assert_identical(actual_isel, actual_sel)

        # with coordinate
        actual_isel = mds.isel(
            x=xr.DataArray(np.arange(3), dims="z", coords={"z": [0, 1, 2]})
        )
        actual_sel = mds.sel(
            x=xr.DataArray(midx[:3], dims="z", coords={"z": [0, 1, 2]})
        )
        assert actual_isel["x"].dims == ("z",)
        assert actual_sel["x"].dims == ("z",)
        assert_identical(actual_isel, actual_sel)

        # Vectorized indexing with level-variables raises an error
        with pytest.raises(ValueError, match=r"Vectorized selection is "):
            mds.sel(one=["a", "b"])

        with pytest.raises(
            ValueError,
            match=r"Vectorized selection is not available along coordinate 'x' with a multi-index",
        ):
            mds.sel(
                x=xr.DataArray(
                    [np.array(midx[:2]), np.array(midx[-2:])], dims=["a", "b"]
                )
            )

    def test_sel_categorical(self) -> None:
        ind = pd.Series(["foo", "bar"], dtype="category")
        df = pd.DataFrame({"ind": ind, "values": [1, 2]})
        ds = df.set_index("ind").to_xarray()
        actual = ds.sel(ind="bar")
        expected = ds.isel(ind=1)
        assert_identical(expected, actual)

    def test_sel_categorical_error(self) -> None:
        ind = pd.Series(["foo", "bar"], dtype="category")
        df = pd.DataFrame({"ind": ind, "values": [1, 2]})
        ds = df.set_index("ind").to_xarray()
        with pytest.raises(ValueError):
            ds.sel(ind="bar", method="nearest")
        with pytest.raises(ValueError):
            ds.sel(ind="bar", tolerance="nearest")  # type: ignore[arg-type]

    def test_categorical_index(self) -> None:
        cat = pd.CategoricalIndex(
            ["foo", "bar", "foo"],
            categories=["foo", "bar", "baz", "qux", "quux", "corge"],
        )
        ds = xr.Dataset(
            {"var": ("cat", np.arange(3))},
            coords={"cat": ("cat", cat), "c": ("cat", [0, 1, 1])},
        )
        # test slice
        actual1 = ds.sel(cat="foo")
        expected1 = ds.isel(cat=[0, 2])
        assert_identical(expected1, actual1)
        # make sure the conversion to the array works
        actual2 = ds.sel(cat="foo")["cat"].values
        assert (actual2 == np.array(["foo", "foo"])).all()

        ds = ds.set_index(index=["cat", "c"])
        actual3 = ds.unstack("index")
        assert actual3["var"].shape == (2, 2)

    def test_categorical_index_reindex(self) -> None:
        cat = pd.CategoricalIndex(
            ["foo", "bar", "baz"],
            categories=["foo", "bar", "baz", "qux", "quux", "corge"],
        )
        ds = xr.Dataset(
            {"var": ("cat", np.arange(3))},
            coords={"cat": ("cat", cat), "c": ("cat", [0, 1, 2])},
        )
        actual = ds.reindex(cat=["foo"])["cat"].values
        assert (actual == np.array(["foo"])).all()

    @pytest.mark.parametrize("fill_value", [np.nan, pd.NA])
    def test_extensionarray_negative_reindex(self, fill_value) -> None:
        cat = pd.Categorical(
            ["foo", "bar", "baz"],
            categories=["foo", "bar", "baz", "qux", "quux", "corge"],
        )
        ds = xr.Dataset(
            {"cat": ("index", cat)},
            coords={"index": ("index", np.arange(3))},
        )
        reindexed_cat = cast(
            pd.api.extensions.ExtensionArray,
            (
                ds.reindex(index=[-1, 1, 1], fill_value=fill_value)["cat"]
                .to_pandas()
                .values
            ),
        )
        assert reindexed_cat.equals(pd.array([pd.NA, "bar", "bar"], dtype=cat.dtype))  # type: ignore[attr-defined]

    def test_extension_array_reindex_same(self) -> None:
        series = pd.Series([1, 2, pd.NA, 3], dtype=pd.Int32Dtype())
        test = xr.Dataset({"test": series})
        res = test.reindex(dim_0=series.index)
        align(res, test, join="exact")

    def test_categorical_multiindex(self) -> None:
        i1 = pd.Series([0, 0])
        cat = pd.CategoricalDtype(categories=["foo", "baz", "bar"])
        i2 = pd.Series(["baz", "bar"], dtype=cat)

        df = pd.DataFrame({"i1": i1, "i2": i2, "values": [1, 2]}).set_index(
            ["i1", "i2"]
        )
        actual = df.to_xarray()
        assert actual["values"].shape == (1, 2)

    def test_sel_drop(self) -> None:
        data = Dataset({"foo": ("x", [1, 2, 3])}, {"x": [0, 1, 2]})
        expected = Dataset({"foo": 1})
        selected = data.sel(x=0, drop=True)
        assert_identical(expected, selected)

        expected = Dataset({"foo": 1}, {"x": 0})
        selected = data.sel(x=0, drop=False)
        assert_identical(expected, selected)

        data = Dataset({"foo": ("x", [1, 2, 3])})
        expected = Dataset({"foo": 1})
        selected = data.sel(x=0, drop=True)
        assert_identical(expected, selected)

    def test_sel_drop_mindex(self) -> None:
        midx = pd.MultiIndex.from_arrays([["a", "a"], [1, 2]], names=("foo", "bar"))
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        data = Dataset(coords=midx_coords)

        actual = data.sel(foo="a", drop=True)
        assert "foo" not in actual.coords

        actual = data.sel(foo="a", drop=False)
        assert_equal(actual.foo, DataArray("a", coords={"foo": "a"}))

    def test_isel_drop(self) -> None:
        data = Dataset({"foo": ("x", [1, 2, 3])}, {"x": [0, 1, 2]})
        expected = Dataset({"foo": 1})
        selected = data.isel(x=0, drop=True)
        assert_identical(expected, selected)

        expected = Dataset({"foo": 1}, {"x": 0})
        selected = data.isel(x=0, drop=False)
        assert_identical(expected, selected)

    def test_head(self) -> None:
        data = create_test_data()

        expected = data.isel(time=slice(5), dim2=slice(6))
        actual = data.head(time=5, dim2=6)
        assert_equal(expected, actual)

        expected = data.isel(time=slice(0))
        actual = data.head(time=0)
        assert_equal(expected, actual)

        expected = data.isel({dim: slice(6) for dim in data.dims})
        actual = data.head(6)
        assert_equal(expected, actual)

        expected = data.isel({dim: slice(5) for dim in data.dims})
        actual = data.head()
        assert_equal(expected, actual)

        with pytest.raises(TypeError, match=r"either dict-like or a single int"):
            data.head([3])  # type: ignore[arg-type]
        with pytest.raises(TypeError, match=r"expected integer type"):
            data.head(dim2=3.1)
        with pytest.raises(ValueError, match=r"expected positive int"):
            data.head(time=-3)

    def test_tail(self) -> None:
        data = create_test_data()

        expected = data.isel(time=slice(-5, None), dim2=slice(-6, None))
        actual = data.tail(time=5, dim2=6)
        assert_equal(expected, actual)

        expected = data.isel(dim1=slice(0))
        actual = data.tail(dim1=0)
        assert_equal(expected, actual)

        expected = data.isel({dim: slice(-6, None) for dim in data.dims})
        actual = data.tail(6)
        assert_equal(expected, actual)

        expected = data.isel({dim: slice(-5, None) for dim in data.dims})
        actual = data.tail()
        assert_equal(expected, actual)

        with pytest.raises(TypeError, match=r"either dict-like or a single int"):
            data.tail([3])  # type: ignore[arg-type]
        with pytest.raises(TypeError, match=r"expected integer type"):
            data.tail(dim2=3.1)
        with pytest.raises(ValueError, match=r"expected positive int"):
            data.tail(time=-3)

    def test_thin(self) -> None:
        data = create_test_data()

        expected = data.isel(time=slice(None, None, 5), dim2=slice(None, None, 6))
        actual = data.thin(time=5, dim2=6)
        assert_equal(expected, actual)

        expected = data.isel({dim: slice(None, None, 6) for dim in data.dims})
        actual = data.thin(6)
        assert_equal(expected, actual)

        with pytest.raises(TypeError, match=r"either dict-like or a single int"):
            data.thin([3])  # type: ignore[arg-type]
        with pytest.raises(TypeError, match=r"expected integer type"):
            data.thin(dim2=3.1)
        with pytest.raises(ValueError, match=r"cannot be zero"):
            data.thin(time=0)
        with pytest.raises(ValueError, match=r"expected positive int"):
            data.thin(time=-3)

    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_sel_fancy(self) -> None:
        data = create_test_data()

        # add in a range() index
        data["dim1"] = data.dim1

        pdim1 = [1, 2, 3]
        pdim2 = [4, 5, 1]
        pdim3 = [1, 2, 3]
        expected = data.isel(
            dim1=Variable(("test_coord",), pdim1),
            dim2=Variable(("test_coord",), pdim2),
            dim3=Variable(("test_coord"), pdim3),
        )
        actual = data.sel(
            dim1=Variable(("test_coord",), data.dim1[pdim1]),
            dim2=Variable(("test_coord",), data.dim2[pdim2]),
            dim3=Variable(("test_coord",), data.dim3[pdim3]),
        )
        assert_identical(expected, actual)

        # DataArray Indexer
        idx_t = DataArray(
            data["time"][[3, 2, 1]].values, dims=["a"], coords={"a": ["a", "b", "c"]}
        )
        idx_2 = DataArray(
            data["dim2"][[3, 2, 1]].values, dims=["a"], coords={"a": ["a", "b", "c"]}
        )
        idx_3 = DataArray(
            data["dim3"][[3, 2, 1]].values, dims=["a"], coords={"a": ["a", "b", "c"]}
        )
        actual = data.sel(time=idx_t, dim2=idx_2, dim3=idx_3)
        expected = data.isel(
            time=Variable(("a",), [3, 2, 1]),
            dim2=Variable(("a",), [3, 2, 1]),
            dim3=Variable(("a",), [3, 2, 1]),
        )
        expected = expected.assign_coords(a=idx_t["a"])
        assert_identical(expected, actual)

        idx_t = DataArray(
            data["time"][[3, 2, 1]].values, dims=["a"], coords={"a": ["a", "b", "c"]}
        )
        idx_2 = DataArray(
            data["dim2"][[2, 1, 3]].values, dims=["b"], coords={"b": [0, 1, 2]}
        )
        idx_3 = DataArray(
            data["dim3"][[1, 2, 1]].values, dims=["c"], coords={"c": [0.0, 1.1, 2.2]}
        )
        actual = data.sel(time=idx_t, dim2=idx_2, dim3=idx_3)
        expected = data.isel(
            time=Variable(("a",), [3, 2, 1]),
            dim2=Variable(("b",), [2, 1, 3]),
            dim3=Variable(("c",), [1, 2, 1]),
        )
        expected = expected.assign_coords(a=idx_t["a"], b=idx_2["b"], c=idx_3["c"])
        assert_identical(expected, actual)

        # test from sel_points
        data = Dataset({"foo": (("x", "y"), np.arange(9).reshape(3, 3))})
        data.coords.update({"x": [0, 1, 2], "y": [0, 1, 2]})

        expected = Dataset(
            {"foo": ("points", [0, 4, 8])},
            coords={
                "x": Variable(("points",), [0, 1, 2]),
                "y": Variable(("points",), [0, 1, 2]),
            },
        )
        actual = data.sel(
            x=Variable(("points",), [0, 1, 2]), y=Variable(("points",), [0, 1, 2])
        )
        assert_identical(expected, actual)

        expected.coords.update({"x": ("points", [0, 1, 2]), "y": ("points", [0, 1, 2])})
        actual = data.sel(
            x=Variable(("points",), [0.1, 1.1, 2.5]),
            y=Variable(("points",), [0, 1.2, 2.0]),
            method="pad",
        )
        assert_identical(expected, actual)

        idx_x = DataArray([0, 1, 2], dims=["a"], coords={"a": ["a", "b", "c"]})
        idx_y = DataArray([0, 2, 1], dims=["b"], coords={"b": [0, 3, 6]})
        expected_ary = data["foo"][[0, 1, 2], [0, 2, 1]]
        actual = data.sel(x=idx_x, y=idx_y)
        assert_array_equal(expected_ary, actual["foo"])
        assert_identical(actual["a"].drop_vars("x"), idx_x["a"])
        assert_identical(actual["b"].drop_vars("y"), idx_y["b"])

        with pytest.raises(KeyError):
            data.sel(x=[2.5], y=[2.0], method="pad", tolerance=1e-3)

    def test_sel_method(self) -> None:
        data = create_test_data()

        expected = data.sel(dim2=1)
        actual = data.sel(dim2=0.95, method="nearest")
        assert_identical(expected, actual)

        actual = data.sel(dim2=0.95, method="nearest", tolerance=1)
        assert_identical(expected, actual)

        with pytest.raises(KeyError):
            actual = data.sel(dim2=np.pi, method="nearest", tolerance=0)

        expected = data.sel(dim2=[1.5])
        actual = data.sel(dim2=[1.45], method="backfill")
        assert_identical(expected, actual)

        with pytest.raises(NotImplementedError, match=r"slice objects"):
            data.sel(dim2=slice(1, 3), method="ffill")

        with pytest.raises(TypeError, match=r"``method``"):
            # this should not pass silently
            data.sel(dim2=1, method=data)  # type: ignore[arg-type]

        # cannot pass method if there is no associated coordinate
        with pytest.raises(ValueError, match=r"cannot supply"):
            data.sel(dim1=0, method="nearest")

    def test_loc(self) -> None:
        data = create_test_data()
        expected = data.sel(dim3="a")
        actual = data.loc[dict(dim3="a")]
        assert_identical(expected, actual)
        with pytest.raises(TypeError, match=r"can only lookup dict"):
            data.loc["a"]  # type: ignore[index]

    def test_selection_multiindex(self) -> None:
        midx = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2], [-1, -2]], names=("one", "two", "three")
        )
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        mdata = Dataset(data_vars={"var": ("x", range(8))}, coords=midx_coords)

        def test_sel(
            lab_indexer, pos_indexer, replaced_idx=False, renamed_dim=None
        ) -> None:
            ds = mdata.sel(x=lab_indexer)
            expected_ds = mdata.isel(x=pos_indexer)
            if not replaced_idx:
                assert_identical(ds, expected_ds)
            else:
                if renamed_dim:
                    assert ds["var"].dims[0] == renamed_dim
                    ds = ds.rename({renamed_dim: "x"})
                assert_identical(ds["var"].variable, expected_ds["var"].variable)
                assert not ds["x"].equals(expected_ds["x"])

        test_sel(("a", 1, -1), 0)
        test_sel(("b", 2, -2), -1)
        test_sel(("a", 1), [0, 1], replaced_idx=True, renamed_dim="three")
        test_sel(("a",), range(4), replaced_idx=True)
        test_sel("a", range(4), replaced_idx=True)
        test_sel([("a", 1, -1), ("b", 2, -2)], [0, 7])
        test_sel(slice("a", "b"), range(8))
        test_sel(slice(("a", 1), ("b", 1)), range(6))
        test_sel({"one": "a", "two": 1, "three": -1}, 0)
        test_sel({"one": "a", "two": 1}, [0, 1], replaced_idx=True, renamed_dim="three")
        test_sel({"one": "a"}, range(4), replaced_idx=True)

        assert_identical(mdata.loc[{"x": {"one": "a"}}], mdata.sel(x={"one": "a"}))
        assert_identical(mdata.loc[{"x": "a"}], mdata.sel(x="a"))
        assert_identical(mdata.loc[{"x": ("a", 1)}], mdata.sel(x=("a", 1)))
        assert_identical(mdata.loc[{"x": ("a", 1, -1)}], mdata.sel(x=("a", 1, -1)))

        assert_identical(mdata.sel(x={"one": "a", "two": 1}), mdata.sel(one="a", two=1))

    def test_broadcast_like(self) -> None:
        original1 = DataArray(
            np.random.randn(5), [("x", range(5))], name="a"
        ).to_dataset()

        original2 = DataArray(np.random.randn(6), [("y", range(6))], name="b")

        expected1, expected2 = broadcast(original1, original2)

        assert_identical(
            original1.broadcast_like(original2), expected1.transpose("y", "x")
        )

        assert_identical(original2.broadcast_like(original1), expected2)

    def test_to_pandas(self) -> None:
        # 0D -> series
        actual = Dataset({"a": 1, "b": 2}).to_pandas()
        expected = pd.Series([1, 2], ["a", "b"])
        assert_array_equal(actual, expected)

        # 1D -> dataframe
        x = np.random.randn(10)
        y = np.random.randn(10)
        t = list("abcdefghij")
        ds = Dataset({"a": ("t", x), "b": ("t", y), "t": ("t", t)})
        actual_df = ds.to_pandas()
        expected_df = ds.to_dataframe()
        assert expected_df.equals(actual_df), (expected_df, actual_df)

        # 2D -> error
        x2d = np.random.randn(10, 10)
        y2d = np.random.randn(10, 10)
        with pytest.raises(ValueError, match=r"cannot convert Datasets"):
            Dataset({"a": (["t", "r"], x2d), "b": (["t", "r"], y2d)}).to_pandas()

    def test_reindex_like(self) -> None:
        data = create_test_data()
        data["letters"] = ("dim3", 10 * ["a"])

        expected = data.isel(dim1=slice(10), time=slice(13))
        actual = data.reindex_like(expected)
        assert_identical(actual, expected)

        expected = data.copy(deep=True)
        expected["dim3"] = ("dim3", list("cdefghijkl"))
        expected["var3"][:-2] = expected["var3"][2:].values
        expected["var3"][-2:] = np.nan
        expected["letters"] = expected["letters"].astype(object)
        expected["letters"][-2:] = np.nan
        expected["numbers"] = expected["numbers"].astype(float)
        expected["numbers"][:-2] = expected["numbers"][2:].values
        expected["numbers"][-2:] = np.nan
        actual = data.reindex_like(expected)
        assert_identical(actual, expected)

    def test_reindex(self) -> None:
        data = create_test_data()
        assert_identical(data, data.reindex())

        expected = data.assign_coords(dim1=data["dim1"])
        actual = data.reindex(dim1=data["dim1"])
        assert_identical(actual, expected)

        actual = data.reindex(dim1=data["dim1"].values)
        assert_identical(actual, expected)

        actual = data.reindex(dim1=data["dim1"].to_index())
        assert_identical(actual, expected)

        with pytest.raises(
            ValueError, match=r"cannot reindex or align along dimension"
        ):
            data.reindex(dim1=data["dim1"][:5])

        expected = data.isel(dim2=slice(5))
        actual = data.reindex(dim2=data["dim2"][:5])
        assert_identical(actual, expected)

        # test dict-like argument
        actual = data.reindex({"dim2": data["dim2"]})
        expected = data
        assert_identical(actual, expected)
        with pytest.raises(ValueError, match=r"cannot specify both"):
            data.reindex({"x": 0}, x=0)
        with pytest.raises(ValueError, match=r"dictionary"):
            data.reindex("foo")  # type: ignore[arg-type]

        # invalid dimension
        # TODO: (benbovy - explicit indexes): uncomment?
        # --> from reindex docstrings: "any mismatched dimension is simply ignored"
        # with pytest.raises(ValueError, match=r"indexer keys.*not correspond.*"):
        #     data.reindex(invalid=0)

        # out of order
        expected = data.sel(dim2=data["dim2"][:5:-1])
        actual = data.reindex(dim2=data["dim2"][:5:-1])
        assert_identical(actual, expected)

        # multiple fill values
        expected = data.reindex(dim2=[0.1, 2.1, 3.1, 4.1]).assign(
            var1=lambda ds: ds.var1.copy(data=[[-10, -10, -10, -10]] * len(ds.dim1)),
            var2=lambda ds: ds.var2.copy(data=[[-20, -20, -20, -20]] * len(ds.dim1)),
        )
        actual = data.reindex(
            dim2=[0.1, 2.1, 3.1, 4.1], fill_value={"var1": -10, "var2": -20}
        )
        assert_identical(actual, expected)
        # use the default value
        expected = data.reindex(dim2=[0.1, 2.1, 3.1, 4.1]).assign(
            var1=lambda ds: ds.var1.copy(data=[[-10, -10, -10, -10]] * len(ds.dim1)),
            var2=lambda ds: ds.var2.copy(
                data=[[np.nan, np.nan, np.nan, np.nan]] * len(ds.dim1)
            ),
        )
        actual = data.reindex(dim2=[0.1, 2.1, 3.1, 4.1], fill_value={"var1": -10})
        assert_identical(actual, expected)

        # regression test for #279
        expected = Dataset({"x": ("time", np.random.randn(5))}, {"time": range(5)})
        time2 = DataArray(np.arange(5), dims="time2")
        with pytest.raises(ValueError):
            actual = expected.reindex(time=time2)

        # another regression test
        ds = Dataset(
            {"foo": (["x", "y"], np.zeros((3, 4)))}, {"x": range(3), "y": range(4)}
        )
        expected = Dataset(
            {"foo": (["x", "y"], np.zeros((3, 2)))}, {"x": [0, 1, 3], "y": [0, 1]}
        )
        expected["foo"][-1] = np.nan
        actual = ds.reindex(x=[0, 1, 3], y=[0, 1])
        assert_identical(expected, actual)

    def test_reindex_attrs_encoding(self) -> None:
        ds = Dataset(
            {"data": ("x", [1, 2, 3])},
            {"x": ("x", [0, 1, 2], {"foo": "bar"}, {"bar": "baz"})},
        )
        actual = ds.reindex(x=[0, 1])
        expected = Dataset(
            {"data": ("x", [1, 2])},
            {"x": ("x", [0, 1], {"foo": "bar"}, {"bar": "baz"})},
        )
        assert_identical(actual, expected)
        assert actual.x.encoding == expected.x.encoding

    def test_reindex_warning(self) -> None:
        data = create_test_data()

        with pytest.raises(ValueError):
            # DataArray with different dimension raises Future warning
            ind = xr.DataArray([0.0, 1.0], dims=["new_dim"], name="ind")
            data.reindex(dim2=ind)

        # Should not warn
        ind = xr.DataArray([0.0, 1.0], dims=["dim2"], name="ind")
        with warnings.catch_warnings(record=True) as ws:
            data.reindex(dim2=ind)
            assert len(ws) == 0

    def test_reindex_variables_copied(self) -> None:
        data = create_test_data()
        reindexed_data = data.reindex(copy=False)
        for k in data.variables:
            assert reindexed_data.variables[k] is not data.variables[k]

    def test_reindex_method(self) -> None:
        ds = Dataset({"x": ("y", [10, 20]), "y": [0, 1]})
        y = [-0.5, 0.5, 1.5]
        actual = ds.reindex(y=y, method="backfill")
        expected = Dataset({"x": ("y", [10, 20, np.nan]), "y": y})
        assert_identical(expected, actual)

        actual = ds.reindex(y=y, method="backfill", tolerance=0.1)
        expected = Dataset({"x": ("y", 3 * [np.nan]), "y": y})
        assert_identical(expected, actual)

        actual = ds.reindex(y=y, method="backfill", tolerance=[0.1, 0.5, 0.1])
        expected = Dataset({"x": ("y", [np.nan, 20, np.nan]), "y": y})
        assert_identical(expected, actual)

        actual = ds.reindex(y=[0.1, 0.1, 1], tolerance=[0, 0.1, 0], method="nearest")
        expected = Dataset({"x": ("y", [np.nan, 10, 20]), "y": [0.1, 0.1, 1]})
        assert_identical(expected, actual)

        actual = ds.reindex(y=y, method="pad")
        expected = Dataset({"x": ("y", [np.nan, 10, 20]), "y": y})
        assert_identical(expected, actual)

        alt = Dataset({"y": y})
        actual = ds.reindex_like(alt, method="pad")
        assert_identical(expected, actual)

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"x": 2, "z": 1}])
    def test_reindex_fill_value(self, fill_value) -> None:
        ds = Dataset({"x": ("y", [10, 20]), "z": ("y", [-20, -10]), "y": [0, 1]})
        y = [0, 1, 2]
        actual = ds.reindex(y=y, fill_value=fill_value)
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_x = fill_value_z = np.nan
        elif isinstance(fill_value, dict):
            fill_value_x = fill_value["x"]
            fill_value_z = fill_value["z"]
        else:
            fill_value_x = fill_value_z = fill_value
        expected = Dataset(
            {
                "x": ("y", [10, 20, fill_value_x]),
                "z": ("y", [-20, -10, fill_value_z]),
                "y": y,
            }
        )
        assert_identical(expected, actual)

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"x": 2, "z": 1}])
    def test_reindex_like_fill_value(self, fill_value) -> None:
        ds = Dataset({"x": ("y", [10, 20]), "z": ("y", [-20, -10]), "y": [0, 1]})
        y = [0, 1, 2]
        alt = Dataset({"y": y})
        actual = ds.reindex_like(alt, fill_value=fill_value)
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_x = fill_value_z = np.nan
        elif isinstance(fill_value, dict):
            fill_value_x = fill_value["x"]
            fill_value_z = fill_value["z"]
        else:
            fill_value_x = fill_value_z = fill_value
        expected = Dataset(
            {
                "x": ("y", [10, 20, fill_value_x]),
                "z": ("y", [-20, -10, fill_value_z]),
                "y": y,
            }
        )
        assert_identical(expected, actual)

    @pytest.mark.parametrize("dtype", [str, bytes])
    def test_reindex_str_dtype(self, dtype) -> None:
        data = Dataset({"data": ("x", [1, 2]), "x": np.array(["a", "b"], dtype=dtype)})

        actual = data.reindex(x=data.x)
        expected = data

        assert_identical(expected, actual)
        assert actual.x.dtype == expected.x.dtype

    def test_reindex_with_multiindex_level(self) -> None:
        # test for https://github.com/pydata/xarray/issues/10347
        mindex = pd.MultiIndex.from_product(
            [[100, 200, 300], [1, 2, 3, 4]], names=["x", "y"]
        )
        y_idx = PandasIndex(mindex.levels[1], "y")

        ds1 = xr.Dataset(coords={"y": [1, 2, 3]})
        ds2 = xr.Dataset(coords=xr.Coordinates.from_xindex(y_idx))

        actual = ds1.reindex(y=ds2.y)
        assert_identical(actual, ds2)

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"foo": 2, "bar": 1}])
    def test_align_fill_value(self, fill_value) -> None:
        x = Dataset({"foo": DataArray([1, 2], dims=["x"], coords={"x": [1, 2]})})
        y = Dataset({"bar": DataArray([1, 2], dims=["x"], coords={"x": [1, 3]})})
        x2, y2 = align(x, y, join="outer", fill_value=fill_value)
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_foo = fill_value_bar = np.nan
        elif isinstance(fill_value, dict):
            fill_value_foo = fill_value["foo"]
            fill_value_bar = fill_value["bar"]
        else:
            fill_value_foo = fill_value_bar = fill_value

        expected_x2 = Dataset(
            {
                "foo": DataArray(
                    [1, 2, fill_value_foo], dims=["x"], coords={"x": [1, 2, 3]}
                )
            }
        )
        expected_y2 = Dataset(
            {
                "bar": DataArray(
                    [1, fill_value_bar, 2], dims=["x"], coords={"x": [1, 2, 3]}
                )
            }
        )
        assert_identical(expected_x2, x2)
        assert_identical(expected_y2, y2)

    def test_align(self) -> None:
        left = create_test_data()
        right = left.copy(deep=True)
        right["dim3"] = ("dim3", list("cdefghijkl"))
        right["var3"][:-2] = right["var3"][2:].values
        right["var3"][-2:] = np.random.randn(*right["var3"][-2:].shape)
        right["numbers"][:-2] = right["numbers"][2:].values
        right["numbers"][-2:] = -10

        intersection = list("cdefghij")
        union = list("abcdefghijkl")

        left2, right2 = align(left, right, join="inner")
        assert_array_equal(left2["dim3"], intersection)
        assert_identical(left2, right2)

        left2, right2 = align(left, right, join="outer")

        assert_array_equal(left2["dim3"], union)
        assert_equal(left2["dim3"].variable, right2["dim3"].variable)

        assert_identical(left2.sel(dim3=intersection), right2.sel(dim3=intersection))
        assert np.isnan(left2["var3"][-2:]).all()
        assert np.isnan(right2["var3"][:2]).all()

        left2, right2 = align(left, right, join="left")
        assert_equal(left2["dim3"].variable, right2["dim3"].variable)
        assert_equal(left2["dim3"].variable, left["dim3"].variable)

        assert_identical(left2.sel(dim3=intersection), right2.sel(dim3=intersection))
        assert np.isnan(right2["var3"][:2]).all()

        left2, right2 = align(left, right, join="right")
        assert_equal(left2["dim3"].variable, right2["dim3"].variable)
        assert_equal(left2["dim3"].variable, right["dim3"].variable)

        assert_identical(left2.sel(dim3=intersection), right2.sel(dim3=intersection))

        assert np.isnan(left2["var3"][-2:]).all()

        with pytest.raises(ValueError, match=r"invalid value for join"):
            align(left, right, join="foobar")  # type: ignore[call-overload]
        with pytest.raises(TypeError):
            align(left, right, foo="bar")  # type: ignore[call-overload]

    def test_align_exact(self) -> None:
        left = xr.Dataset(coords={"x": [0, 1]})
        right = xr.Dataset(coords={"x": [1, 2]})

        left1, left2 = xr.align(left, left, join="exact")
        assert_identical(left1, left)
        assert_identical(left2, left)

        with pytest.raises(ValueError, match=r"cannot align.*join.*exact.*not equal.*"):
            xr.align(left, right, join="exact")

    def test_align_override(self) -> None:
        left = xr.Dataset(coords={"x": [0, 1, 2]})
        right = xr.Dataset(coords={"x": [0.1, 1.1, 2.1], "y": [1, 2, 3]})
        expected_right = xr.Dataset(coords={"x": [0, 1, 2], "y": [1, 2, 3]})

        new_left, new_right = xr.align(left, right, join="override")
        assert_identical(left, new_left)
        assert_identical(new_right, expected_right)

        new_left, new_right = xr.align(left, right, exclude="x", join="override")
        assert_identical(left, new_left)
        assert_identical(right, new_right)

        new_left, new_right = xr.align(
            left.isel(x=0, drop=True), right, exclude="x", join="override"
        )
        assert_identical(left.isel(x=0, drop=True), new_left)
        assert_identical(right, new_right)

        with pytest.raises(
            ValueError, match=r"cannot align.*join.*override.*same size"
        ):
            xr.align(left.isel(x=0).expand_dims("x"), right, join="override")

    def test_align_exclude(self) -> None:
        x = Dataset(
            {
                "foo": DataArray(
                    [[1, 2], [3, 4]], dims=["x", "y"], coords={"x": [1, 2], "y": [3, 4]}
                )
            }
        )
        y = Dataset(
            {
                "bar": DataArray(
                    [[1, 2], [3, 4]], dims=["x", "y"], coords={"x": [1, 3], "y": [5, 6]}
                )
            }
        )
        x2, y2 = align(x, y, exclude=["y"], join="outer")

        expected_x2 = Dataset(
            {
                "foo": DataArray(
                    [[1, 2], [3, 4], [np.nan, np.nan]],
                    dims=["x", "y"],
                    coords={"x": [1, 2, 3], "y": [3, 4]},
                )
            }
        )
        expected_y2 = Dataset(
            {
                "bar": DataArray(
                    [[1, 2], [np.nan, np.nan], [3, 4]],
                    dims=["x", "y"],
                    coords={"x": [1, 2, 3], "y": [5, 6]},
                )
            }
        )
        assert_identical(expected_x2, x2)
        assert_identical(expected_y2, y2)

    def test_align_nocopy(self) -> None:
        x = Dataset({"foo": DataArray([1, 2, 3], coords=[("x", [1, 2, 3])])})
        y = Dataset({"foo": DataArray([1, 2], coords=[("x", [1, 2])])})
        expected_x2 = x
        expected_y2 = Dataset(
            {"foo": DataArray([1, 2, np.nan], coords=[("x", [1, 2, 3])])}
        )

        x2, y2 = align(x, y, copy=False, join="outer")
        assert_identical(expected_x2, x2)
        assert_identical(expected_y2, y2)
        assert source_ndarray(x["foo"].data) is source_ndarray(x2["foo"].data)

        x2, y2 = align(x, y, copy=True, join="outer")
        assert source_ndarray(x["foo"].data) is not source_ndarray(x2["foo"].data)
        assert_identical(expected_x2, x2)
        assert_identical(expected_y2, y2)

    def test_align_indexes(self) -> None:
        x = Dataset({"foo": DataArray([1, 2, 3], dims="x", coords=[("x", [1, 2, 3])])})
        (x2,) = align(x, indexes={"x": [2, 3, 1]})
        expected_x2 = Dataset(
            {"foo": DataArray([2, 3, 1], dims="x", coords={"x": [2, 3, 1]})}
        )

        assert_identical(expected_x2, x2)

    def test_align_multiple_indexes_common_dim(self) -> None:
        a = Dataset(coords={"x": [1, 2], "xb": ("x", [3, 4])}).set_xindex("xb")
        b = Dataset(coords={"x": [1], "xb": ("x", [3])}).set_xindex("xb")

        (a2, b2) = align(a, b, join="inner")
        assert_identical(a2, b, check_default_indexes=False)
        assert_identical(b2, b, check_default_indexes=False)

        c = Dataset(coords={"x": [1, 3], "xb": ("x", [2, 4])}).set_xindex("xb")

        with pytest.raises(AlignmentError, match=".*conflicting re-indexers"):
            align(a, c)

    def test_align_conflicting_indexes(self) -> None:
        class CustomIndex(PandasIndex): ...

        a = Dataset(coords={"xb": ("x", [3, 4])}).set_xindex("xb")
        b = Dataset(coords={"xb": ("x", [3])}).set_xindex("xb", CustomIndex)

        with pytest.raises(AlignmentError, match="cannot align.*conflicting indexes"):
            align(a, b)

    def test_align_non_unique(self) -> None:
        x = Dataset({"foo": ("x", [3, 4, 5]), "x": [0, 0, 1]})
        x1, x2 = align(x, x)
        assert_identical(x1, x)
        assert_identical(x2, x)

        y = Dataset({"bar": ("x", [6, 7]), "x": [0, 1]})
        with pytest.raises(ValueError, match=r"cannot reindex or align"):
            align(x, y)

    def test_align_str_dtype(self) -> None:
        a = Dataset({"foo": ("x", [0, 1])}, coords={"x": ["a", "b"]})
        b = Dataset({"foo": ("x", [1, 2])}, coords={"x": ["b", "c"]})

        expected_a = Dataset(
            {"foo": ("x", [0, 1, np.nan])}, coords={"x": ["a", "b", "c"]}
        )
        expected_b = Dataset(
            {"foo": ("x", [np.nan, 1, 2])}, coords={"x": ["a", "b", "c"]}
        )

        actual_a, actual_b = xr.align(a, b, join="outer")

        assert_identical(expected_a, actual_a)
        assert expected_a.x.dtype == actual_a.x.dtype

        assert_identical(expected_b, actual_b)
        assert expected_b.x.dtype == actual_b.x.dtype

    @pytest.mark.parametrize("join", ["left", "override"])
    def test_align_index_var_attrs(self, join) -> None:
        # regression test https://github.com/pydata/xarray/issues/6852
        # aligning two objects should have no side effect on their index variable
        # metadata.

        ds = Dataset(coords={"x": ("x", [1, 2, 3], {"units": "m"})})
        ds_noattr = Dataset(coords={"x": ("x", [1, 2, 3])})

        xr.align(ds_noattr, ds, join=join)

        assert ds.x.attrs == {"units": "m"}
        assert ds_noattr.x.attrs == {}

    def test_align_scalar_index(self) -> None:
        # ensure that indexes associated with scalar coordinates are not ignored
        # during alignment
        ds1 = Dataset(coords={"x": 0}).set_xindex("x", ScalarIndex)
        ds2 = Dataset(coords={"x": 0}).set_xindex("x", ScalarIndex)

        actual = xr.align(ds1, ds2, join="exact")
        assert_identical(actual[0], ds1, check_default_indexes=False)
        assert_identical(actual[1], ds2, check_default_indexes=False)

        ds3 = Dataset(coords={"x": 1}).set_xindex("x", ScalarIndex)

        with pytest.raises(AlignmentError, match="cannot align objects"):
            xr.align(ds1, ds3, join="exact")

    def test_align_multi_dim_index_exclude_dims(self) -> None:
        ds1 = (
            Dataset(coords={"x": [1, 2], "y": [3, 4]})
            .drop_indexes(["x", "y"])
            .set_xindex(["x", "y"], XYIndex)
        )
        ds2 = (
            Dataset(coords={"x": [1, 2], "y": [5, 6]})
            .drop_indexes(["x", "y"])
            .set_xindex(["x", "y"], XYIndex)
        )

        for join in ("outer", "exact"):
            actual = xr.align(ds1, ds2, join=join, exclude="y")
            assert_identical(actual[0], ds1, check_default_indexes=False)
            assert_identical(actual[1], ds2, check_default_indexes=False)

        with pytest.raises(
            AlignmentError, match="cannot align objects.*index.*not equal"
        ):
            xr.align(ds1, ds2, join="exact")

        with pytest.raises(AlignmentError, match="cannot exclude dimension"):
            xr.align(ds1, ds2, join="override", exclude="y")

    def test_align_index_equals_future_warning(self) -> None:
        # TODO: remove this test once the deprecation cycle is completed
        class DeprecatedEqualsSignatureIndex(PandasIndex):
            def equals(self, other: Index) -> bool:  # type: ignore[override]
                return super().equals(other, exclude=None)

        ds = (
            Dataset(coords={"x": [1, 2]})
            .drop_indexes("x")
            .set_xindex("x", DeprecatedEqualsSignatureIndex)
        )

        with pytest.warns(FutureWarning, match="signature.*deprecated"):
            xr.align(ds, ds.copy(), join="exact")

    def test_broadcast(self) -> None:
        ds = Dataset(
            {"foo": 0, "bar": ("x", [1]), "baz": ("y", [2, 3])}, {"c": ("x", [4])}
        )
        expected = Dataset(
            {
                "foo": (("x", "y"), [[0, 0]]),
                "bar": (("x", "y"), [[1, 1]]),
                "baz": (("x", "y"), [[2, 3]]),
            },
            {"c": ("x", [4])},
        )
        (actual,) = broadcast(ds)
        assert_identical(expected, actual)

        ds_x = Dataset({"foo": ("x", [1])})
        ds_y = Dataset({"bar": ("y", [2, 3])})
        expected_x = Dataset({"foo": (("x", "y"), [[1, 1]])})
        expected_y = Dataset({"bar": (("x", "y"), [[2, 3]])})
        actual_x, actual_y = broadcast(ds_x, ds_y)
        assert_identical(expected_x, actual_x)
        assert_identical(expected_y, actual_y)

        array_y = ds_y["bar"]
        expected_y2 = expected_y["bar"]
        actual_x2, actual_y2 = broadcast(ds_x, array_y)
        assert_identical(expected_x, actual_x2)
        assert_identical(expected_y2, actual_y2)

    def test_broadcast_nocopy(self) -> None:
        # Test that data is not copied if not needed
        x = Dataset({"foo": (("x", "y"), [[1, 1]])})
        y = Dataset({"bar": ("y", [2, 3])})

        (actual_x,) = broadcast(x)
        assert_identical(x, actual_x)
        assert source_ndarray(actual_x["foo"].data) is source_ndarray(x["foo"].data)

        actual_x, actual_y = broadcast(x, y)
        assert_identical(x, actual_x)
        assert source_ndarray(actual_x["foo"].data) is source_ndarray(x["foo"].data)

    def test_broadcast_exclude(self) -> None:
        x = Dataset(
            {
                "foo": DataArray(
                    [[1, 2], [3, 4]], dims=["x", "y"], coords={"x": [1, 2], "y": [3, 4]}
                ),
                "bar": DataArray(5),
            }
        )
        y = Dataset(
            {
                "foo": DataArray(
                    [[1, 2]], dims=["z", "y"], coords={"z": [1], "y": [5, 6]}
                )
            }
        )
        x2, y2 = broadcast(x, y, exclude=["y"])

        expected_x2 = Dataset(
            {
                "foo": DataArray(
                    [[[1, 2]], [[3, 4]]],
                    dims=["x", "z", "y"],
                    coords={"z": [1], "x": [1, 2], "y": [3, 4]},
                ),
                "bar": DataArray(
                    [[5], [5]], dims=["x", "z"], coords={"x": [1, 2], "z": [1]}
                ),
            }
        )
        expected_y2 = Dataset(
            {
                "foo": DataArray(
                    [[[1, 2]], [[1, 2]]],
                    dims=["x", "z", "y"],
                    coords={"z": [1], "x": [1, 2], "y": [5, 6]},
                )
            }
        )
        assert_identical(expected_x2, x2)
        assert_identical(expected_y2, y2)

    def test_broadcast_misaligned(self) -> None:
        x = Dataset({"foo": DataArray([1, 2, 3], coords=[("x", [-1, -2, -3])])})
        y = Dataset(
            {
                "bar": DataArray(
                    [[1, 2], [3, 4]],
                    dims=["y", "x"],
                    coords={"y": [1, 2], "x": [10, -3]},
                )
            }
        )
        x2, y2 = broadcast(x, y)
        expected_x2 = Dataset(
            {
                "foo": DataArray(
                    [[3, 3], [2, 2], [1, 1], [np.nan, np.nan]],
                    dims=["x", "y"],
                    coords={"y": [1, 2], "x": [-3, -2, -1, 10]},
                )
            }
        )
        expected_y2 = Dataset(
            {
                "bar": DataArray(
                    [[2, 4], [np.nan, np.nan], [np.nan, np.nan], [1, 3]],
                    dims=["x", "y"],
                    coords={"y": [1, 2], "x": [-3, -2, -1, 10]},
                )
            }
        )
        assert_identical(expected_x2, x2)
        assert_identical(expected_y2, y2)

    def test_broadcast_multi_index(self) -> None:
        # GH6430
        ds = Dataset(
            {"foo": (("x", "y", "z"), np.ones((3, 4, 2)))},
            {"x": ["a", "b", "c"], "y": [1, 2, 3, 4]},
        )
        stacked = ds.stack(space=["x", "y"])
        broadcasted, _ = broadcast(stacked, stacked.space)

        assert broadcasted.xindexes["x"] is broadcasted.xindexes["space"]
        assert broadcasted.xindexes["y"] is broadcasted.xindexes["space"]

    def test_variable_indexing(self) -> None:
        data = create_test_data()
        v = data["var1"]
        d1 = data["dim1"]
        d2 = data["dim2"]
        assert_equal(v, v[d1.values])
        assert_equal(v, v[d1])
        assert_equal(v[:3], v[d1 < 3])
        assert_equal(v[:, 3:], v[:, d2 >= 1.5])
        assert_equal(v[:3, 3:], v[d1 < 3, d2 >= 1.5])
        assert_equal(v[:3, :2], v[range(3), range(2)])
        assert_equal(v[:3, :2], v.loc[d1[:3], d2[:2]])

    def test_drop_variables(self) -> None:
        data = create_test_data()

        assert_identical(data, data.drop_vars([]))

        expected = Dataset({k: data[k] for k in data.variables if k != "time"})
        actual = data.drop_vars("time")
        assert_identical(expected, actual)
        actual = data.drop_vars(["time"])
        assert_identical(expected, actual)

        with pytest.raises(
            ValueError,
            match=re.escape(
                "These variables cannot be found in this dataset: ['not_found_here']"
            ),
        ):
            data.drop_vars("not_found_here")

        actual = data.drop_vars("not_found_here", errors="ignore")
        assert_identical(data, actual)

        actual = data.drop_vars(["not_found_here"], errors="ignore")
        assert_identical(data, actual)

        actual = data.drop_vars(["time", "not_found_here"], errors="ignore")
        assert_identical(expected, actual)

        # deprecated approach with `drop` works (straight copy paste from above)

        with pytest.warns(DeprecationWarning):
            actual = data.drop("not_found_here", errors="ignore")
        assert_identical(data, actual)

        with pytest.warns(DeprecationWarning):
            actual = data.drop(["not_found_here"], errors="ignore")
        assert_identical(data, actual)

        with pytest.warns(DeprecationWarning):
            actual = data.drop(["time", "not_found_here"], errors="ignore")
        assert_identical(expected, actual)

        with pytest.warns(DeprecationWarning):
            actual = data.drop({"time", "not_found_here"}, errors="ignore")
        assert_identical(expected, actual)

    def test_drop_multiindex_level(self) -> None:
        data = create_test_multiindex()
        expected = data.drop_vars(["x", "level_1", "level_2"])
        with pytest.warns(DeprecationWarning):
            actual = data.drop_vars("level_1")
        assert_identical(expected, actual)

    def test_drop_index_labels(self) -> None:
        data = Dataset({"A": (["x", "y"], np.random.randn(2, 3)), "x": ["a", "b"]})

        with pytest.warns(DeprecationWarning):
            actual = data.drop(["a"], dim="x")
        expected = data.isel(x=[1])
        assert_identical(expected, actual)

        with pytest.warns(DeprecationWarning):
            actual = data.drop(["a", "b"], dim="x")
        expected = data.isel(x=slice(0, 0))
        assert_identical(expected, actual)

        with pytest.raises(KeyError):
            # not contained in axis
            with pytest.warns(DeprecationWarning):
                data.drop(["c"], dim="x")

        with pytest.warns(DeprecationWarning):
            actual = data.drop(["c"], dim="x", errors="ignore")
        assert_identical(data, actual)

        with pytest.raises(ValueError):
            data.drop(["c"], dim="x", errors="wrong_value")  # type: ignore[arg-type]

        with pytest.warns(DeprecationWarning):
            actual = data.drop(["a", "b", "c"], "x", errors="ignore")
        expected = data.isel(x=slice(0, 0))
        assert_identical(expected, actual)

        # DataArrays as labels are a nasty corner case as they are not
        # Iterable[Hashable] - DataArray.__iter__ yields scalar DataArrays.
        actual = data.drop_sel(x=DataArray(["a", "b", "c"]), errors="ignore")
        expected = data.isel(x=slice(0, 0))
        assert_identical(expected, actual)
        with pytest.warns(DeprecationWarning):
            data.drop(DataArray(["a", "b", "c"]), dim="x", errors="ignore")
        assert_identical(expected, actual)

        actual = data.drop_sel(y=[1])
        expected = data.isel(y=[0, 2])
        assert_identical(expected, actual)

        with pytest.raises(KeyError, match=r"not found in axis"):
            data.drop_sel(x=0)

    def test_drop_labels_by_keyword(self) -> None:
        data = Dataset(
            {"A": (["x", "y"], np.random.randn(2, 6)), "x": ["a", "b"], "y": range(6)}
        )
        # Basic functionality.
        assert len(data.coords["x"]) == 2

        with pytest.warns(DeprecationWarning):
            ds1 = data.drop(["a"], dim="x")
        ds2 = data.drop_sel(x="a")
        ds3 = data.drop_sel(x=["a"])
        ds4 = data.drop_sel(x=["a", "b"])
        ds5 = data.drop_sel(x=["a", "b"], y=range(0, 6, 2))

        arr = DataArray(range(3), dims=["c"])
        with pytest.warns(DeprecationWarning):
            data.drop(arr.coords)
        with pytest.warns(DeprecationWarning):
            data.drop(arr.xindexes)

        assert_array_equal(ds1.coords["x"], ["b"])
        assert_array_equal(ds2.coords["x"], ["b"])
        assert_array_equal(ds3.coords["x"], ["b"])
        assert ds4.coords["x"].size == 0
        assert ds5.coords["x"].size == 0
        assert_array_equal(ds5.coords["y"], [1, 3, 5])

        # Error handling if user tries both approaches.
        with pytest.raises(ValueError):
            data.drop(labels=["a"], x="a")
        with pytest.raises(ValueError):
            data.drop(labels=["a"], dim="x", x="a")
        warnings.filterwarnings("ignore", r"\W*drop")
        with pytest.raises(ValueError):
            data.drop(dim="x", x="a")

    def test_drop_labels_by_position(self) -> None:
        data = Dataset(
            {"A": (["x", "y"], np.random.randn(2, 6)), "x": ["a", "b"], "y": range(6)}
        )
        # Basic functionality.
        assert len(data.coords["x"]) == 2

        actual = data.drop_isel(x=0)
        expected = data.drop_sel(x="a")
        assert_identical(expected, actual)

        actual = data.drop_isel(x=[0])
        expected = data.drop_sel(x=["a"])
        assert_identical(expected, actual)

        actual = data.drop_isel(x=[0, 1])
        expected = data.drop_sel(x=["a", "b"])
        assert_identical(expected, actual)
        assert actual.coords["x"].size == 0

        actual = data.drop_isel(x=[0, 1], y=range(0, 6, 2))
        expected = data.drop_sel(x=["a", "b"], y=range(0, 6, 2))
        assert_identical(expected, actual)
        assert actual.coords["x"].size == 0

        with pytest.raises(KeyError):
            data.drop_isel(z=1)

    def test_drop_indexes(self) -> None:
        ds = Dataset(
            coords={
                "x": ("x", [0, 1, 2]),
                "y": ("y", [3, 4, 5]),
                "foo": ("x", ["a", "a", "b"]),
            }
        )

        actual = ds.drop_indexes("x")
        assert "x" not in actual.xindexes
        assert type(actual.x.variable) is Variable

        actual = ds.drop_indexes(["x", "y"])
        assert "x" not in actual.xindexes
        assert "y" not in actual.xindexes
        assert type(actual.x.variable) is Variable
        assert type(actual.y.variable) is Variable

        with pytest.raises(
            ValueError,
            match=r"The coordinates \('not_a_coord',\) are not found in the dataset coordinates",
        ):
            ds.drop_indexes("not_a_coord")

        with pytest.raises(ValueError, match="those coordinates do not have an index"):
            ds.drop_indexes("foo")

        actual = ds.drop_indexes(["foo", "not_a_coord"], errors="ignore")
        assert_identical(actual, ds)

        # test index corrupted
        midx = pd.MultiIndex.from_tuples([(1, 2), (3, 4)], names=["a", "b"])
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        ds = Dataset(coords=midx_coords)

        with pytest.raises(ValueError, match=".*would corrupt the following index.*"):
            ds.drop_indexes("a")

    def test_drop_dims(self) -> None:
        data = xr.Dataset(
            {
                "A": (["x", "y"], np.random.randn(2, 3)),
                "B": ("x", np.random.randn(2)),
                "x": ["a", "b"],
                "z": np.pi,
            }
        )

        actual = data.drop_dims("x")
        expected = data.drop_vars(["A", "B", "x"])
        assert_identical(expected, actual)

        actual = data.drop_dims("y")
        expected = data.drop_vars("A")
        assert_identical(expected, actual)

        actual = data.drop_dims(["x", "y"])
        expected = data.drop_vars(["A", "B", "x"])
        assert_identical(expected, actual)

        with pytest.raises((ValueError, KeyError)):
            data.drop_dims("z")  # not a dimension

        with pytest.raises((ValueError, KeyError)):
            data.drop_dims(None)  # type:ignore[arg-type]

        actual = data.drop_dims("z", errors="ignore")
        assert_identical(data, actual)

        # should this be allowed?
        actual = data.drop_dims(None, errors="ignore")  # type:ignore[arg-type]
        assert_identical(data, actual)

        with pytest.raises(ValueError):
            actual = data.drop_dims("z", errors="wrong_value")  # type: ignore[arg-type]

        actual = data.drop_dims(["x", "y", "z"], errors="ignore")
        expected = data.drop_vars(["A", "B", "x"])
        assert_identical(expected, actual)

    def test_copy(self) -> None:
        data = create_test_data()
        data.attrs["Test"] = [1, 2, 3]

        for copied in [data.copy(deep=False), copy(data)]:
            assert_identical(data, copied)
            assert data.encoding == copied.encoding
            # Note: IndexVariable objects with string dtype are always
            # copied because of xarray.core.indexes.safe_cast_to_index.
            # Limiting the test to data variables.
            for k in data.data_vars:
                v0 = data.variables[k]
                v1 = copied.variables[k]
                assert source_ndarray(v0.data) is source_ndarray(v1.data)
            copied["foo"] = ("z", np.arange(5))
            assert "foo" not in data

            copied.attrs["foo"] = "bar"
            assert "foo" not in data.attrs
            assert data.attrs["Test"] is copied.attrs["Test"]

        for copied in [data.copy(deep=True), deepcopy(data)]:
            assert_identical(data, copied)
            for k, v0 in data.variables.items():
                v1 = copied.variables[k]
                assert v0 is not v1

            assert data.attrs["Test"] is not copied.attrs["Test"]

    def test_copy_with_data(self) -> None:
        orig = create_test_data()
        new_data = {k: np.random.randn(*v.shape) for k, v in orig.data_vars.items()}
        actual = orig.copy(data=new_data)

        expected = orig.copy()
        for k, v in new_data.items():
            expected[k].data = v
        assert_identical(expected, actual)

    @pytest.mark.xfail(raises=AssertionError)
    @pytest.mark.parametrize(
        "deep, expected_orig",
        [
            [
                True,
                xr.DataArray(
                    xr.IndexVariable("a", np.array([1, 2])),
                    coords={"a": [1, 2]},
                    dims=["a"],
                ),
            ],
            [
                False,
                xr.DataArray(
                    xr.IndexVariable("a", np.array([999, 2])),
                    coords={"a": [999, 2]},
                    dims=["a"],
                ),
            ],
        ],
    )
    def test_copy_coords(self, deep, expected_orig) -> None:
        """The test fails for the shallow copy, and apparently only on Windows
        for some reason. In windows coords seem to be immutable unless it's one
        dataset deep copied from another."""
        ds = xr.DataArray(
            np.ones([2, 2, 2]),
            coords={"a": [1, 2], "b": ["x", "y"], "c": [0, 1]},
            dims=["a", "b", "c"],
            name="value",
        ).to_dataset()
        ds_cp = ds.copy(deep=deep)
        new_a = np.array([999, 2])
        ds_cp.coords["a"] = ds_cp.a.copy(data=new_a)

        expected_cp = xr.DataArray(
            xr.IndexVariable("a", new_a),
            coords={"a": [999, 2]},
            dims=["a"],
        )
        assert_identical(ds_cp.coords["a"], expected_cp)

        assert_identical(ds.coords["a"], expected_orig)

    def test_copy_with_data_errors(self) -> None:
        orig = create_test_data()
        new_var1 = np.arange(orig["var1"].size).reshape(orig["var1"].shape)
        with pytest.raises(ValueError, match=r"Data must be dict-like"):
            orig.copy(data=new_var1)  # type: ignore[arg-type]
        with pytest.raises(ValueError, match=r"only contain variables in original"):
            orig.copy(data={"not_in_original": new_var1})
        with pytest.raises(ValueError, match=r"contain all variables in original"):
            orig.copy(data={"var1": new_var1})

    def test_drop_encoding(self) -> None:
        orig = create_test_data()
        vencoding = {"scale_factor": 10}
        orig.encoding = {"foo": "bar"}

        for k in orig.variables.keys():
            orig[k].encoding = vencoding

        actual = orig.drop_encoding()
        assert actual.encoding == {}
        for v in actual.variables.values():
            assert v.encoding == {}

        assert_equal(actual, orig)

    def test_rename(self) -> None:
        data = create_test_data()
        newnames = {
            "var1": "renamed_var1",
            "dim2": "renamed_dim2",
        }
        renamed = data.rename(newnames)

        variables = dict(data.variables)
        for nk, nv in newnames.items():
            variables[nv] = variables.pop(nk)

        for k, v in variables.items():
            dims = list(v.dims)
            for name, newname in newnames.items():
                if name in dims:
                    dims[dims.index(name)] = newname

            assert_equal(
                Variable(dims, v.values, v.attrs),
                renamed[k].variable.to_base_variable(),
            )
            assert v.encoding == renamed[k].encoding
            assert type(v) is type(renamed.variables[k])

        assert "var1" not in renamed
        assert "dim2" not in renamed

        with pytest.raises(ValueError, match=r"cannot rename 'not_a_var'"):
            data.rename({"not_a_var": "nada"})

        with pytest.raises(ValueError, match=r"'var1' conflicts"):
            data.rename({"var2": "var1"})

        # verify that we can rename a variable without accessing the data
        var1 = data["var1"]
        data["var1"] = (var1.dims, InaccessibleArray(var1.values))
        renamed = data.rename(newnames)
        with pytest.raises(UnexpectedDataAccess):
            _ = renamed["renamed_var1"].values

        # https://github.com/python/mypy/issues/10008
        renamed_kwargs = data.rename(**newnames)  # type: ignore[arg-type]
        assert_identical(renamed, renamed_kwargs)

    def test_rename_old_name(self) -> None:
        # regtest for GH1477
        data = create_test_data()

        with pytest.raises(ValueError, match=r"'samecol' conflicts"):
            data.rename({"var1": "samecol", "var2": "samecol"})

        # This shouldn't cause any problems.
        data.rename({"var1": "var2", "var2": "var1"})

    def test_rename_same_name(self) -> None:
        data = create_test_data()
        newnames = {"var1": "var1", "dim2": "dim2"}
        renamed = data.rename(newnames)
        assert_identical(renamed, data)

    def test_rename_dims(self) -> None:
        original = Dataset({"x": ("x", [0, 1, 2]), "y": ("x", [10, 11, 12]), "z": 42})
        expected = Dataset(
            {"x": ("x_new", [0, 1, 2]), "y": ("x_new", [10, 11, 12]), "z": 42}
        )
        # TODO: (benbovy - explicit indexes) update when set_index supports
        # setting index for non-dimension variables
        expected = expected.set_coords("x")
        actual = original.rename_dims({"x": "x_new"})
        assert_identical(expected, actual, check_default_indexes=False)
        actual_2 = original.rename_dims(x="x_new")
        assert_identical(expected, actual_2, check_default_indexes=False)

        # Test to raise ValueError
        dims_dict_bad = {"x_bad": "x_new"}
        with pytest.raises(ValueError):
            original.rename_dims(dims_dict_bad)

        with pytest.raises(ValueError):
            original.rename_dims({"x": "z"})

    def test_rename_vars(self) -> None:
        original = Dataset({"x": ("x", [0, 1, 2]), "y": ("x", [10, 11, 12]), "z": 42})
        expected = Dataset(
            {"x_new": ("x", [0, 1, 2]), "y": ("x", [10, 11, 12]), "z": 42}
        )
        # TODO: (benbovy - explicit indexes) update when set_index supports
        # setting index for non-dimension variables
        expected = expected.set_coords("x_new")
        actual = original.rename_vars({"x": "x_new"})
        assert_identical(expected, actual, check_default_indexes=False)
        actual_2 = original.rename_vars(x="x_new")
        assert_identical(expected, actual_2, check_default_indexes=False)

        # Test to raise ValueError
        names_dict_bad = {"x_bad": "x_new"}
        with pytest.raises(ValueError):
            original.rename_vars(names_dict_bad)

    def test_rename_dimension_coord(self) -> None:
        # rename a dimension corodinate to a non-dimension coordinate
        # should preserve index
        original = Dataset(coords={"x": ("x", [0, 1, 2])})

        actual = original.rename_vars({"x": "x_new"})
        assert "x_new" in actual.xindexes

        actual_2 = original.rename_dims({"x": "x_new"})
        assert "x" in actual_2.xindexes

    def test_rename_dimension_coord_warnings(self) -> None:
        # create a dimension coordinate by renaming a dimension or coordinate
        # should raise a warning (no index created)
        ds = Dataset(coords={"x": ("y", [0, 1])})

        with pytest.warns(
            UserWarning, match="rename 'x' to 'y' does not create an index.*"
        ):
            ds.rename(x="y")

        ds = Dataset(coords={"y": ("x", [0, 1])})

        with pytest.warns(
            UserWarning, match="rename 'x' to 'y' does not create an index.*"
        ):
            ds.rename(x="y")

        # No operation should not raise a warning
        ds = Dataset(
            data_vars={"data": (("x", "y"), np.ones((2, 3)))},
            coords={"x": range(2), "y": range(3), "a": ("x", [3, 4])},
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ds.rename(x="x")

    def test_rename_multiindex(self) -> None:
        midx = pd.MultiIndex.from_tuples([(1, 2), (3, 4)], names=["a", "b"])
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        original = Dataset({}, midx_coords)

        midx_renamed = midx.rename(["a", "c"])
        midx_coords_renamed = Coordinates.from_pandas_multiindex(midx_renamed, "x")
        expected = Dataset({}, midx_coords_renamed)

        actual = original.rename({"b": "c"})
        assert_identical(expected, actual)

        with pytest.raises(ValueError, match=r"'a' conflicts"):
            with pytest.warns(UserWarning, match="does not create an index anymore"):
                original.rename({"x": "a"})

        with pytest.raises(ValueError, match=r"'x' conflicts"):
            with pytest.warns(UserWarning, match="does not create an index anymore"):
                original.rename({"a": "x"})

        with pytest.raises(ValueError, match=r"'b' conflicts"):
            original.rename({"a": "b"})

    def test_rename_preserve_attrs_encoding(self) -> None:
        # test propagate attrs/encoding to new variable(s) created from Index object
        original = Dataset(coords={"x": ("x", [0, 1, 2])})
        expected = Dataset(coords={"y": ("y", [0, 1, 2])})
        for ds, dim in zip([original, expected], ["x", "y"], strict=True):
            ds[dim].attrs = {"foo": "bar"}
            ds[dim].encoding = {"foo": "bar"}

        actual = original.rename({"x": "y"})
        assert_identical(actual, expected)

    @requires_cftime
    def test_rename_does_not_change_CFTimeIndex_type(self) -> None:
        # make sure CFTimeIndex is not converted to DatetimeIndex #3522

        time = xr.date_range(
            start="2000", periods=6, freq="2MS", calendar="noleap", use_cftime=True
        )
        orig = Dataset(coords={"time": time})

        renamed = orig.rename(time="time_new")
        assert "time_new" in renamed.xindexes
        # TODO: benbovy - flexible indexes: update when CFTimeIndex
        # inherits from xarray.Index
        assert isinstance(renamed.xindexes["time_new"].to_pandas_index(), CFTimeIndex)
        assert renamed.xindexes["time_new"].to_pandas_index().name == "time_new"

        # check original has not changed
        assert "time" in orig.xindexes
        assert isinstance(orig.xindexes["time"].to_pandas_index(), CFTimeIndex)
        assert orig.xindexes["time"].to_pandas_index().name == "time"

        # note: rename_dims(time="time_new") drops "ds.indexes"
        renamed = orig.rename_dims()
        assert isinstance(renamed.xindexes["time"].to_pandas_index(), CFTimeIndex)

        renamed = orig.rename_vars()
        assert isinstance(renamed.xindexes["time"].to_pandas_index(), CFTimeIndex)

    def test_rename_does_not_change_DatetimeIndex_type(self) -> None:
        # make sure DatetimeIndex is conderved on rename

        time = pd.date_range(start="2000", periods=6, freq="2MS")
        orig = Dataset(coords={"time": time})

        renamed = orig.rename(time="time_new")
        assert "time_new" in renamed.xindexes
        # TODO: benbovy - flexible indexes: update when DatetimeIndex
        # inherits from xarray.Index?
        assert isinstance(renamed.xindexes["time_new"].to_pandas_index(), DatetimeIndex)
        assert renamed.xindexes["time_new"].to_pandas_index().name == "time_new"

        # check original has not changed
        assert "time" in orig.xindexes
        assert isinstance(orig.xindexes["time"].to_pandas_index(), DatetimeIndex)
        assert orig.xindexes["time"].to_pandas_index().name == "time"

        # note: rename_dims(time="time_new") drops "ds.indexes"
        renamed = orig.rename_dims()
        assert isinstance(renamed.xindexes["time"].to_pandas_index(), DatetimeIndex)

        renamed = orig.rename_vars()
        assert isinstance(renamed.xindexes["time"].to_pandas_index(), DatetimeIndex)

    def test_swap_dims(self) -> None:
        original = Dataset({"x": [1, 2, 3], "y": ("x", list("abc")), "z": 42})
        expected = Dataset({"z": 42}, {"x": ("y", [1, 2, 3]), "y": list("abc")})
        actual = original.swap_dims({"x": "y"})
        assert_identical(expected, actual)
        assert isinstance(actual.variables["y"], IndexVariable)
        assert isinstance(actual.variables["x"], Variable)
        assert actual.xindexes["y"].equals(expected.xindexes["y"])

        roundtripped = actual.swap_dims({"y": "x"})
        assert_identical(original.set_coords("y"), roundtripped)

        with pytest.raises(ValueError, match=r"cannot swap"):
            original.swap_dims({"y": "x"})
        with pytest.raises(ValueError, match=r"replacement dimension"):
            original.swap_dims({"x": "z"})

        expected = Dataset(
            {"y": ("u", list("abc")), "z": 42}, coords={"x": ("u", [1, 2, 3])}
        )
        actual = original.swap_dims({"x": "u"})
        assert_identical(expected, actual)

        # as kwargs
        expected = Dataset(
            {"y": ("u", list("abc")), "z": 42}, coords={"x": ("u", [1, 2, 3])}
        )
        actual = original.swap_dims(x="u")
        assert_identical(expected, actual)

        # handle multiindex case
        midx = pd.MultiIndex.from_arrays([list("aab"), list("yzz")], names=["y1", "y2"])

        original = Dataset({"x": [1, 2, 3], "y": ("x", midx), "z": 42})

        midx_coords = Coordinates.from_pandas_multiindex(midx, "y")
        midx_coords["x"] = ("y", [1, 2, 3])
        expected = Dataset({"z": 42}, midx_coords)

        actual = original.swap_dims({"x": "y"})
        assert_identical(expected, actual)
        assert isinstance(actual.variables["y"], IndexVariable)
        assert isinstance(actual.variables["x"], Variable)
        assert actual.xindexes["y"].equals(expected.xindexes["y"])

    def test_expand_dims_error(self) -> None:
        original = Dataset(
            {
                "x": ("a", np.random.randn(3)),
                "y": (["b", "a"], np.random.randn(4, 3)),
                "z": ("a", np.random.randn(3)),
            },
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
            attrs={"key": "entry"},
        )

        with pytest.raises(ValueError, match=r"already exists"):
            original.expand_dims(dim=["x"])

        # Make sure it raises true error also for non-dimensional coordinates
        # which has dimension.
        original = original.set_coords("z")
        with pytest.raises(ValueError, match=r"already exists"):
            original.expand_dims(dim=["z"])

        original = Dataset(
            {
                "x": ("a", np.random.randn(3)),
                "y": (["b", "a"], np.random.randn(4, 3)),
                "z": ("a", np.random.randn(3)),
            },
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
            attrs={"key": "entry"},
        )
        with pytest.raises(TypeError, match=r"value of new dimension"):
            original.expand_dims({"d": 3.2})
        with pytest.raises(ValueError, match=r"both keyword and positional"):
            original.expand_dims({"d": 4}, e=4)

    def test_expand_dims_int(self) -> None:
        original = Dataset(
            {"x": ("a", np.random.randn(3)), "y": (["b", "a"], np.random.randn(4, 3))},
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
            attrs={"key": "entry"},
        )

        actual = original.expand_dims(["z"], [1])
        expected = Dataset(
            {
                "x": original["x"].expand_dims("z", 1),
                "y": original["y"].expand_dims("z", 1),
            },
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
            attrs={"key": "entry"},
        )
        assert_identical(expected, actual)
        # make sure squeeze restores the original data set.
        roundtripped = actual.squeeze("z")
        assert_identical(original, roundtripped)

        # another test with a negative axis
        actual = original.expand_dims(["z"], [-1])
        expected = Dataset(
            {
                "x": original["x"].expand_dims("z", -1),
                "y": original["y"].expand_dims("z", -1),
            },
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
            attrs={"key": "entry"},
        )
        assert_identical(expected, actual)
        # make sure squeeze restores the original data set.
        roundtripped = actual.squeeze("z")
        assert_identical(original, roundtripped)

    def test_expand_dims_coords(self) -> None:
        original = Dataset({"x": ("a", np.array([1, 2, 3]))})
        expected = Dataset(
            {"x": (("b", "a"), np.array([[1, 2, 3], [1, 2, 3]]))}, coords={"b": [1, 2]}
        )
        actual = original.expand_dims(dict(b=[1, 2]))
        assert_identical(expected, actual)
        assert "b" not in original._coord_names

    def test_expand_dims_existing_scalar_coord(self) -> None:
        original = Dataset({"x": 1}, {"a": 2})
        expected = Dataset({"x": (("a",), [1])}, {"a": [2]})
        actual = original.expand_dims("a")
        assert_identical(expected, actual)

    def test_isel_expand_dims_roundtrip(self) -> None:
        original = Dataset({"x": (("a",), [1])}, {"a": [2]})
        actual = original.isel(a=0).expand_dims("a")
        assert_identical(actual, original)

    def test_expand_dims_mixed_int_and_coords(self) -> None:
        # Test expanding one dimension to have size > 1 that doesn't have
        # coordinates, and also expanding another dimension to have size > 1
        # that DOES have coordinates.
        original = Dataset(
            {"x": ("a", np.random.randn(3)), "y": (["b", "a"], np.random.randn(4, 3))},
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
        )

        actual = original.expand_dims({"d": 4, "e": ["l", "m", "n"]})

        expected = Dataset(
            {
                "x": xr.DataArray(
                    original["x"].values * np.ones([4, 3, 3]),
                    coords=dict(d=range(4), e=["l", "m", "n"], a=np.linspace(0, 1, 3)),
                    dims=["d", "e", "a"],
                ).drop_vars("d"),
                "y": xr.DataArray(
                    original["y"].values * np.ones([4, 3, 4, 3]),
                    coords=dict(
                        d=range(4),
                        e=["l", "m", "n"],
                        b=np.linspace(0, 1, 4),
                        a=np.linspace(0, 1, 3),
                    ),
                    dims=["d", "e", "b", "a"],
                ).drop_vars("d"),
            },
            coords={"c": np.linspace(0, 1, 5)},
        )
        assert_identical(actual, expected)

    def test_expand_dims_kwargs_python36plus(self) -> None:
        original = Dataset(
            {"x": ("a", np.random.randn(3)), "y": (["b", "a"], np.random.randn(4, 3))},
            coords={
                "a": np.linspace(0, 1, 3),
                "b": np.linspace(0, 1, 4),
                "c": np.linspace(0, 1, 5),
            },
            attrs={"key": "entry"},
        )
        other_way = original.expand_dims(e=["l", "m", "n"])
        other_way_expected = Dataset(
            {
                "x": xr.DataArray(
                    original["x"].values * np.ones([3, 3]),
                    coords=dict(e=["l", "m", "n"], a=np.linspace(0, 1, 3)),
                    dims=["e", "a"],
                ),
                "y": xr.DataArray(
                    original["y"].values * np.ones([3, 4, 3]),
                    coords=dict(
                        e=["l", "m", "n"],
                        b=np.linspace(0, 1, 4),
                        a=np.linspace(0, 1, 3),
                    ),
                    dims=["e", "b", "a"],
                ),
            },
            coords={"c": np.linspace(0, 1, 5)},
            attrs={"key": "entry"},
        )
        assert_identical(other_way_expected, other_way)

    @pytest.mark.parametrize("create_index_for_new_dim_flag", [True, False])
    def test_expand_dims_create_index_data_variable(
        self, create_index_for_new_dim_flag
    ):
        # data variables should not gain an index ever
        ds = Dataset({"x": 0})

        if create_index_for_new_dim_flag:
            with pytest.warns(UserWarning, match="No index created"):
                expanded = ds.expand_dims(
                    "x", create_index_for_new_dim=create_index_for_new_dim_flag
                )
        else:
            expanded = ds.expand_dims(
                "x", create_index_for_new_dim=create_index_for_new_dim_flag
            )

        # TODO Can't just create the expected dataset directly using constructor because of GH issue 8959
        expected = Dataset({"x": ("x", [0])}).drop_indexes("x").reset_coords("x")

        assert_identical(expanded, expected, check_default_indexes=False)
        assert expanded.indexes == {}

    def test_expand_dims_create_index_coordinate_variable(self):
        # coordinate variables should gain an index only if create_index_for_new_dim is True (the default)
        ds = Dataset(coords={"x": 0})
        expanded = ds.expand_dims("x")
        expected = Dataset({"x": ("x", [0])})
        assert_identical(expanded, expected)

        expanded_no_index = ds.expand_dims("x", create_index_for_new_dim=False)

        # TODO Can't just create the expected dataset directly using constructor because of GH issue 8959
        expected = Dataset(coords={"x": ("x", [0])}).drop_indexes("x")

        assert_identical(expanded_no_index, expected, check_default_indexes=False)
        assert expanded_no_index.indexes == {}

    def test_expand_dims_create_index_from_iterable(self):
        ds = Dataset(coords={"x": 0})
        expanded = ds.expand_dims(x=[0, 1])
        expected = Dataset({"x": ("x", [0, 1])})
        assert_identical(expanded, expected)

        expanded_no_index = ds.expand_dims(x=[0, 1], create_index_for_new_dim=False)

        # TODO Can't just create the expected dataset directly using constructor because of GH issue 8959
        expected = Dataset(coords={"x": ("x", [0, 1])}).drop_indexes("x")

        assert_identical(expanded, expected, check_default_indexes=False)
        assert expanded_no_index.indexes == {}

    def test_expand_dims_non_nanosecond_conversion(self) -> None:
        # Regression test for https://github.com/pydata/xarray/issues/7493#issuecomment-1953091000
        # todo: test still needed?
        ds = Dataset().expand_dims({"time": [np.datetime64("2018-01-01", "m")]})
        assert ds.time.dtype == np.dtype("datetime64[s]")

    def test_set_index(self) -> None:
        expected = create_test_multiindex()
        mindex = expected["x"].to_index()
        indexes = [mindex.get_level_values(str(n)) for n in mindex.names]
        coords = {idx.name: ("x", idx) for idx in indexes}
        ds = Dataset({}, coords=coords)

        obj = ds.set_index(x=mindex.names)
        assert_identical(obj, expected)

        # ensure pre-existing indexes involved are removed
        # (level_2 should be a coordinate with no index)
        ds = create_test_multiindex()
        coords = {"x": coords["level_1"], "level_2": coords["level_2"]}
        expected = Dataset({}, coords=coords)

        obj = ds.set_index(x="level_1")
        assert_identical(obj, expected)

        # ensure set_index with no existing index and a single data var given
        # doesn't return multi-index
        ds = Dataset(data_vars={"x_var": ("x", [0, 1, 2])})
        expected = Dataset(coords={"x": [0, 1, 2]})
        assert_identical(ds.set_index(x="x_var"), expected)

        with pytest.raises(ValueError, match=r"bar variable\(s\) do not exist"):
            ds.set_index(foo="bar")

        with pytest.raises(ValueError, match=r"dimension mismatch.*"):
            ds.set_index(y="x_var")

        ds = Dataset(coords={"x": 1})
        with pytest.raises(
            ValueError, match=r".*cannot set a PandasIndex.*scalar variable.*"
        ):
            ds.set_index(x="x")

    def test_set_index_deindexed_coords(self) -> None:
        # test de-indexed coordinates are converted to base variable
        # https://github.com/pydata/xarray/issues/6969
        one = ["a", "a", "b", "b"]
        two = [1, 2, 1, 2]
        three = ["c", "c", "d", "d"]
        four = [3, 4, 3, 4]

        midx_12 = pd.MultiIndex.from_arrays([one, two], names=["one", "two"])
        midx_34 = pd.MultiIndex.from_arrays([three, four], names=["three", "four"])

        coords = Coordinates.from_pandas_multiindex(midx_12, "x")
        coords["three"] = ("x", three)
        coords["four"] = ("x", four)
        ds = xr.Dataset(coords=coords)
        actual = ds.set_index(x=["three", "four"])

        coords_expected = Coordinates.from_pandas_multiindex(midx_34, "x")
        coords_expected["one"] = ("x", one)
        coords_expected["two"] = ("x", two)
        expected = xr.Dataset(coords=coords_expected)

        assert_identical(actual, expected)

    def test_reset_index(self) -> None:
        ds = create_test_multiindex()
        mindex = ds["x"].to_index()
        indexes = [mindex.get_level_values(str(n)) for n in mindex.names]
        coords = {idx.name: ("x", idx) for idx in indexes}
        expected = Dataset({}, coords=coords)

        obj = ds.reset_index("x")
        assert_identical(obj, expected, check_default_indexes=False)
        assert len(obj.xindexes) == 0

        ds = Dataset(coords={"y": ("x", [1, 2, 3])})
        with pytest.raises(ValueError, match=r".*not coordinates with an index"):
            ds.reset_index("y")

    def test_reset_index_keep_attrs(self) -> None:
        coord_1 = DataArray([1, 2], dims=["coord_1"], attrs={"attrs": True})
        ds = Dataset({}, {"coord_1": coord_1})
        obj = ds.reset_index("coord_1")
        assert ds.coord_1.attrs == obj.coord_1.attrs
        assert len(obj.xindexes) == 0

    def test_reset_index_drop_dims(self) -> None:
        ds = Dataset(coords={"x": [1, 2]})
        reset = ds.reset_index("x", drop=True)
        assert len(reset.dims) == 0

    @pytest.mark.parametrize(
        ["arg", "drop", "dropped", "converted", "renamed"],
        [
            ("foo", False, [], [], {"bar": "x"}),
            ("foo", True, ["foo"], [], {"bar": "x"}),
            ("x", False, ["x"], ["foo", "bar"], {}),
            ("x", True, ["x", "foo", "bar"], [], {}),
            (["foo", "bar"], False, ["x"], ["foo", "bar"], {}),
            (["foo", "bar"], True, ["x", "foo", "bar"], [], {}),
            (["x", "foo"], False, ["x"], ["foo", "bar"], {}),
            (["foo", "x"], True, ["x", "foo", "bar"], [], {}),
        ],
    )
    def test_reset_index_drop_convert(
        self,
        arg: str | list[str],
        drop: bool,
        dropped: list[str],
        converted: list[str],
        renamed: dict[str, str],
    ) -> None:
        # regressions https://github.com/pydata/xarray/issues/6946 and
        # https://github.com/pydata/xarray/issues/6989
        # check that multi-index dimension or level coordinates are dropped, converted
        # from IndexVariable to Variable or renamed to dimension as expected
        midx = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=("foo", "bar"))
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        ds = xr.Dataset(coords=midx_coords)
        reset = ds.reset_index(arg, drop=drop)

        for name in dropped:
            assert name not in reset.variables
        for name in converted:
            assert_identical(reset[name].variable, ds[name].variable.to_base_variable())
        for old_name, new_name in renamed.items():
            assert_identical(ds[old_name].variable, reset[new_name].variable)

    def test_reorder_levels(self) -> None:
        ds = create_test_multiindex()
        mindex = ds["x"].to_index()
        assert isinstance(mindex, pd.MultiIndex)
        midx = mindex.reorder_levels(["level_2", "level_1"])
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        expected = Dataset({}, coords=midx_coords)

        # check attrs propagated
        ds["level_1"].attrs["foo"] = "bar"
        expected["level_1"].attrs["foo"] = "bar"

        reindexed = ds.reorder_levels(x=["level_2", "level_1"])
        assert_identical(reindexed, expected)

        ds = Dataset({}, coords={"x": [1, 2]})
        with pytest.raises(ValueError, match=r"has no MultiIndex"):
            ds.reorder_levels(x=["level_1", "level_2"])

    def test_set_xindex(self) -> None:
        ds = Dataset(
            coords={"foo": ("x", ["a", "a", "b", "b"]), "bar": ("x", [0, 1, 2, 3])}
        )

        actual = ds.set_xindex("foo")
        expected = ds.set_index(x="foo").rename_vars(x="foo")
        assert_identical(actual, expected, check_default_indexes=False)

        actual_mindex = ds.set_xindex(["foo", "bar"])
        expected_mindex = ds.set_index(x=["foo", "bar"])
        assert_identical(actual_mindex, expected_mindex)

        class NotAnIndex: ...

        with pytest.raises(TypeError, match=".*not a subclass of xarray.Index"):
            ds.set_xindex("foo", NotAnIndex)  # type: ignore[arg-type]

        with pytest.raises(ValueError, match="those variables don't exist"):
            ds.set_xindex("not_a_coordinate", PandasIndex)

        ds["data_var"] = ("x", [1, 2, 3, 4])

        with pytest.raises(ValueError, match="those variables are data variables"):
            ds.set_xindex("data_var", PandasIndex)

        ds2 = Dataset(coords={"x": ("x", [0, 1, 2, 3])})

        with pytest.raises(ValueError, match="those coordinates already have an index"):
            ds2.set_xindex("x", PandasIndex)

    def test_set_xindex_options(self) -> None:
        ds = Dataset(coords={"foo": ("x", ["a", "a", "b", "b"])})

        class IndexWithOptions(Index):
            def __init__(self, opt):
                self.opt = opt

            @classmethod
            def from_variables(cls, variables, options):
                return cls(options["opt"])

        indexed = ds.set_xindex("foo", IndexWithOptions, opt=1)
        assert indexed.xindexes["foo"].opt == 1  # type: ignore[attr-defined]

    def test_stack(self) -> None:
        ds = Dataset(
            data_vars={"b": (("x", "y"), [[0, 1], [2, 3]])},
            coords={"x": ("x", [0, 1]), "y": ["a", "b"]},
        )

        midx_expected = pd.MultiIndex.from_product(
            [[0, 1], ["a", "b"]], names=["x", "y"]
        )
        midx_coords_expected = Coordinates.from_pandas_multiindex(midx_expected, "z")
        expected = Dataset(
            data_vars={"b": ("z", [0, 1, 2, 3])}, coords=midx_coords_expected
        )
        # check attrs propagated
        ds["x"].attrs["foo"] = "bar"
        expected["x"].attrs["foo"] = "bar"

        actual = ds.stack(z=["x", "y"])
        assert_identical(expected, actual)
        assert list(actual.xindexes) == ["z", "x", "y"]

        actual = ds.stack(z=[...])
        assert_identical(expected, actual)

        # non list dims with ellipsis
        actual = ds.stack(z=(...,))
        assert_identical(expected, actual)

        # ellipsis with given dim
        actual = ds.stack(z=[..., "y"])
        assert_identical(expected, actual)

        midx_expected = pd.MultiIndex.from_product(
            [["a", "b"], [0, 1]], names=["y", "x"]
        )
        midx_coords_expected = Coordinates.from_pandas_multiindex(midx_expected, "z")
        expected = Dataset(
            data_vars={"b": ("z", [0, 2, 1, 3])}, coords=midx_coords_expected
        )
        expected["x"].attrs["foo"] = "bar"

        actual = ds.stack(z=["y", "x"])
        assert_identical(expected, actual)
        assert list(actual.xindexes) == ["z", "y", "x"]

    @pytest.mark.parametrize(
        "create_index,expected_keys",
        [
            (True, ["z", "x", "y"]),
            (False, []),
            (None, ["z", "x", "y"]),
        ],
    )
    def test_stack_create_index(self, create_index, expected_keys) -> None:
        ds = Dataset(
            data_vars={"b": (("x", "y"), [[0, 1], [2, 3]])},
            coords={"x": ("x", [0, 1]), "y": ["a", "b"]},
        )

        actual = ds.stack(z=["x", "y"], create_index=create_index)
        assert list(actual.xindexes) == expected_keys

        # TODO: benbovy (flexible indexes) - test error multiple indexes found
        # along dimension + create_index=True

    def test_stack_multi_index(self) -> None:
        # multi-index on a dimension to stack is discarded too
        midx = pd.MultiIndex.from_product([["a", "b"], [0, 1]], names=("lvl1", "lvl2"))
        coords = Coordinates.from_pandas_multiindex(midx, "x")
        coords["y"] = [0, 1]
        ds = xr.Dataset(
            data_vars={"b": (("x", "y"), [[0, 1], [2, 3], [4, 5], [6, 7]])},
            coords=coords,
        )
        expected = Dataset(
            data_vars={"b": ("z", [0, 1, 2, 3, 4, 5, 6, 7])},
            coords={
                "x": ("z", np.repeat(midx.values, 2)),
                "lvl1": ("z", np.repeat(midx.get_level_values("lvl1"), 2)),
                "lvl2": ("z", np.repeat(midx.get_level_values("lvl2"), 2)),
                "y": ("z", [0, 1, 0, 1] * 2),
            },
        )
        actual = ds.stack(z=["x", "y"], create_index=False)
        assert_identical(expected, actual)
        assert len(actual.xindexes) == 0

        with pytest.raises(ValueError, match=r"cannot create.*wraps a multi-index"):
            ds.stack(z=["x", "y"], create_index=True)

    def test_stack_non_dim_coords(self) -> None:
        ds = Dataset(
            data_vars={"b": (("x", "y"), [[0, 1], [2, 3]])},
            coords={"x": ("x", [0, 1]), "y": ["a", "b"]},
        ).rename_vars(x="xx")

        exp_index = pd.MultiIndex.from_product([[0, 1], ["a", "b"]], names=["xx", "y"])
        exp_coords = Coordinates.from_pandas_multiindex(exp_index, "z")
        expected = Dataset(data_vars={"b": ("z", [0, 1, 2, 3])}, coords=exp_coords)

        actual = ds.stack(z=["x", "y"])
        assert_identical(expected, actual)
        assert list(actual.xindexes) == ["z", "xx", "y"]

    def test_unstack(self) -> None:
        index = pd.MultiIndex.from_product([[0, 1], ["a", "b"]], names=["x", "y"])
        coords = Coordinates.from_pandas_multiindex(index, "z")
        ds = Dataset(data_vars={"b": ("z", [0, 1, 2, 3])}, coords=coords)
        expected = Dataset(
            {"b": (("x", "y"), [[0, 1], [2, 3]]), "x": [0, 1], "y": ["a", "b"]}
        )

        # check attrs propagated
        ds["x"].attrs["foo"] = "bar"
        expected["x"].attrs["foo"] = "bar"

        for dim in ["z", ["z"], None]:
            actual = ds.unstack(dim)
            assert_identical(actual, expected)

    def test_unstack_errors(self) -> None:
        ds = Dataset({"x": [1, 2, 3]})
        with pytest.raises(
            ValueError,
            match=re.escape("Dimensions ('foo',) not found in data dimensions ('x',)"),
        ):
            ds.unstack("foo")
        with pytest.raises(ValueError, match=r".*do not have exactly one multi-index"):
            ds.unstack("x")

        ds = Dataset({"da": [1, 2]}, coords={"y": ("x", [1, 1]), "z": ("x", [0, 0])})
        ds = ds.set_index(x=("y", "z"))

        with pytest.raises(
            ValueError, match="Cannot unstack MultiIndex containing duplicates"
        ):
            ds.unstack("x")

    def test_unstack_fill_value(self) -> None:
        ds = xr.Dataset(
            {"var": (("x",), np.arange(6)), "other_var": (("x",), np.arange(3, 9))},
            coords={"x": [0, 1, 2] * 2, "y": (("x",), ["a"] * 3 + ["b"] * 3)},
        )
        # make ds incomplete
        ds = ds.isel(x=[0, 2, 3, 4]).set_index(index=["x", "y"])
        # test fill_value
        actual1 = ds.unstack("index", fill_value=-1)
        expected1 = ds.unstack("index").fillna(-1).astype(int)
        assert actual1["var"].dtype == int
        assert_equal(actual1, expected1)

        actual2 = ds["var"].unstack("index", fill_value=-1)
        expected2 = ds["var"].unstack("index").fillna(-1).astype(int)
        assert_equal(actual2, expected2)

        actual3 = ds.unstack("index", fill_value={"var": -1, "other_var": 1})
        expected3 = ds.unstack("index").fillna({"var": -1, "other_var": 1}).astype(int)
        assert_equal(actual3, expected3)

    @requires_sparse
    def test_unstack_sparse(self) -> None:
        ds = xr.Dataset(
            {"var": (("x",), np.arange(6))},
            coords={"x": [0, 1, 2] * 2, "y": (("x",), ["a"] * 3 + ["b"] * 3)},
        )
        # make ds incomplete
        ds = ds.isel(x=[0, 2, 3, 4]).set_index(index=["x", "y"])
        # test fill_value
        actual1 = ds.unstack("index", sparse=True)
        expected1 = ds.unstack("index")
        assert isinstance(actual1["var"].data, sparse_array_type)
        assert actual1["var"].variable._to_dense().equals(expected1["var"].variable)
        assert actual1["var"].data.density < 1.0

        actual2 = ds["var"].unstack("index", sparse=True)
        expected2 = ds["var"].unstack("index")
        assert isinstance(actual2.data, sparse_array_type)
        assert actual2.variable._to_dense().equals(expected2.variable)
        assert actual2.data.density < 1.0

        midx = pd.MultiIndex.from_arrays([np.arange(3), np.arange(3)], names=["a", "b"])
        coords = Coordinates.from_pandas_multiindex(midx, "z")
        coords["foo"] = np.arange(4)
        coords["bar"] = np.arange(5)
        ds_eye = Dataset(
            {"var": (("z", "foo", "bar"), np.ones((3, 4, 5)))}, coords=coords
        )
        actual3 = ds_eye.unstack(sparse=True, fill_value=0)
        assert isinstance(actual3["var"].data, sparse_array_type)
        expected3 = xr.Dataset(
            {
                "var": (
                    ("foo", "bar", "a", "b"),
                    np.broadcast_to(np.eye(3, 3), (4, 5, 3, 3)),
                )
            },
            coords={
                "foo": np.arange(4),
                "bar": np.arange(5),
                "a": np.arange(3),
                "b": np.arange(3),
            },
        )
        actual3["var"].data = actual3["var"].data.todense()
        assert_equal(expected3, actual3)

    def test_stack_unstack_fast(self) -> None:
        ds = Dataset(
            {
                "a": ("x", [0, 1]),
                "b": (("x", "y"), [[0, 1], [2, 3]]),
                "x": [0, 1],
                "y": ["a", "b"],
            }
        )
        actual = ds.stack(z=["x", "y"]).unstack("z")
        assert actual.broadcast_equals(ds)

        actual = ds[["b"]].stack(z=["x", "y"]).unstack("z")
        assert actual.identical(ds[["b"]])

    def test_stack_unstack_slow(self) -> None:
        ds = Dataset(
            data_vars={
                "a": ("x", [0, 1]),
                "b": (("x", "y"), [[0, 1], [2, 3]]),
            },
            coords={"x": [0, 1], "y": ["a", "b"]},
        )
        stacked = ds.stack(z=["x", "y"])
        actual = stacked.isel(z=slice(None, None, -1)).unstack("z")
        assert actual.broadcast_equals(ds)

        stacked = ds[["b"]].stack(z=["x", "y"])
        actual = stacked.isel(z=slice(None, None, -1)).unstack("z")
        assert actual.identical(ds[["b"]])

    def test_to_stacked_array_invalid_sample_dims(self) -> None:
        data = xr.Dataset(
            data_vars={"a": (("x", "y"), [[0, 1, 2], [3, 4, 5]]), "b": ("x", [6, 7])},
            coords={"y": ["u", "v", "w"]},
        )
        with pytest.raises(
            ValueError,
            match=r"Variables in the dataset must contain all ``sample_dims`` \(\['y'\]\) but 'b' misses \['y'\]",
        ):
            data.to_stacked_array("features", sample_dims=["y"])

    def test_to_stacked_array_name(self) -> None:
        name = "adf9d"

        # make a two dimensional dataset
        a, b = create_test_stacked_array()
        D = xr.Dataset({"a": a, "b": b})
        sample_dims = ["x"]

        y = D.to_stacked_array("features", sample_dims, name=name)
        assert y.name == name

    def test_to_stacked_array_dtype_dims(self) -> None:
        # make a two dimensional dataset
        a, b = create_test_stacked_array()
        D = xr.Dataset({"a": a, "b": b})
        sample_dims = ["x"]
        y = D.to_stacked_array("features", sample_dims)
        mindex = y.xindexes["features"].to_pandas_index()
        assert isinstance(mindex, pd.MultiIndex)
        assert mindex.levels[1].dtype == D.y.dtype
        assert y.dims == ("x", "features")

    def test_to_stacked_array_to_unstacked_dataset(self) -> None:
        # single dimension: regression test for GH4049
        arr = xr.DataArray(np.arange(3), coords=[("x", [0, 1, 2])])
        data = xr.Dataset({"a": arr, "b": arr})
        stacked = data.to_stacked_array("y", sample_dims=["x"])
        unstacked = stacked.to_unstacked_dataset("y")
        assert_identical(unstacked, data)

        # make a two dimensional dataset
        a, b = create_test_stacked_array()
        D = xr.Dataset({"a": a, "b": b})
        sample_dims = ["x"]
        y = D.to_stacked_array("features", sample_dims).transpose("x", "features")

        x = y.to_unstacked_dataset("features")
        assert_identical(D, x)

        # test on just one sample
        x0 = y[0].to_unstacked_dataset("features")
        d0 = D.isel(x=0)
        assert_identical(d0, x0)

    def test_to_stacked_array_to_unstacked_dataset_different_dimension(self) -> None:
        # test when variables have different dimensionality
        a, b = create_test_stacked_array()
        sample_dims = ["x"]
        D = xr.Dataset({"a": a, "b": b.isel(y=0)})

        y = D.to_stacked_array("features", sample_dims)
        x = y.to_unstacked_dataset("features")
        assert_identical(D, x)

    def test_to_stacked_array_preserves_dtype(self) -> None:
        # regression test for bug found in https://github.com/pydata/xarray/pull/8872#issuecomment-2081218616
        ds = xr.Dataset(
            data_vars={
                "a": (("x", "y"), [[0, 1, 2], [3, 4, 5]]),
                "b": ("x", [6, 7]),
            },
            coords={"y": ["u", "v", "w"]},
        )
        stacked = ds.to_stacked_array("z", sample_dims=["x"])

        # coordinate created from variables names should be of string dtype
        data = np.array(["a", "a", "a", "b"], dtype="<U1")
        expected_stacked_variable = DataArray(name="variable", data=data, dims="z")
        assert_identical(
            stacked.coords["variable"].drop_vars(["z", "variable", "y"]),
            expected_stacked_variable,
        )

    def test_to_stacked_array_transposed(self) -> None:
        # test that to_stacked_array uses updated dim order after transposition
        ds = xr.Dataset(
            data_vars=dict(
                v1=(["d1", "d2"], np.arange(6).reshape((2, 3))),
            ),
            coords=dict(
                d1=(["d1"], np.arange(2)),
                d2=(["d2"], np.arange(3)),
            ),
        )
        da = ds.to_stacked_array(
            new_dim="new_dim",
            sample_dims=[],
            variable_dim="variable",
        )
        dsT = ds.transpose()
        daT = dsT.to_stacked_array(
            new_dim="new_dim",
            sample_dims=[],
            variable_dim="variable",
        )
        v1 = np.arange(6)
        v1T = np.arange(6).reshape((2, 3)).T.flatten()
        np.testing.assert_equal(da.to_numpy(), v1)
        np.testing.assert_equal(daT.to_numpy(), v1T)

    def test_update(self) -> None:
        data = create_test_data(seed=0)
        expected = data.copy()
        var2 = Variable("dim1", np.arange(8))
        actual = data
        actual.update({"var2": var2})
        expected["var2"] = var2
        assert_identical(expected, actual)

        actual = data.copy()
        actual.update(data)
        assert_identical(expected, actual)

        other = Dataset(attrs={"new": "attr"})
        actual = data.copy()
        actual.update(other)
        assert_identical(expected, actual)

    def test_update_overwrite_coords(self) -> None:
        data = Dataset({"a": ("x", [1, 2])}, {"b": 3})
        data.update(Dataset(coords={"b": 4}))
        expected = Dataset({"a": ("x", [1, 2])}, {"b": 4})
        assert_identical(data, expected)

        data = Dataset({"a": ("x", [1, 2])}, {"b": 3})
        data.update(Dataset({"c": 5}, coords={"b": 4}))
        expected = Dataset({"a": ("x", [1, 2]), "c": 5}, {"b": 4})
        assert_identical(data, expected)

        data = Dataset({"a": ("x", [1, 2])}, {"b": 3})
        data.update({"c": DataArray(5, coords={"b": 4})})
        expected = Dataset({"a": ("x", [1, 2]), "c": 5}, {"b": 3})
        assert_identical(data, expected)

    def test_update_multiindex_level(self) -> None:
        data = create_test_multiindex()

        with pytest.raises(
            ValueError, match=r"cannot set or update variable.*corrupt.*index "
        ):
            data.update({"level_1": range(4)})

    def test_update_auto_align(self) -> None:
        ds = Dataset({"x": ("t", [3, 4])}, {"t": [0, 1]})

        expected1 = Dataset(
            {"x": ("t", [3, 4]), "y": ("t", [np.nan, 5])}, {"t": [0, 1]}
        )
        actual1 = ds.copy()
        other1 = {"y": ("t", [5]), "t": [1]}
        with pytest.raises(ValueError, match=r"conflicting sizes"):
            actual1.update(other1)
        actual1.update(Dataset(other1))
        assert_identical(expected1, actual1)

        actual2 = ds.copy()
        other2 = Dataset({"y": ("t", [5]), "t": [100]})
        actual2.update(other2)
        expected2 = Dataset(
            {"x": ("t", [3, 4]), "y": ("t", [np.nan] * 2)}, {"t": [0, 1]}
        )
        assert_identical(expected2, actual2)

    def test_getitem(self) -> None:
        data = create_test_data()
        assert isinstance(data["var1"], DataArray)
        assert_equal(data["var1"].variable, data.variables["var1"])
        with pytest.raises(KeyError):
            data["notfound"]
        with pytest.raises(KeyError):
            data[["var1", "notfound"]]
        with pytest.raises(
            KeyError,
            match=r"Hint: use a list to select multiple variables, for example `ds\[\['var1', 'var2'\]\]`",
        ):
            data["var1", "var2"]

        actual1 = data[["var1", "var2"]]
        expected1 = Dataset({"var1": data["var1"], "var2": data["var2"]})
        assert_equal(expected1, actual1)

        actual2 = data["numbers"]
        expected2 = DataArray(
            data["numbers"].variable,
            {"dim3": data["dim3"], "numbers": data["numbers"]},
            dims="dim3",
            name="numbers",
        )
        assert_identical(expected2, actual2)

        actual3 = data[dict(dim1=0)]
        expected3 = data.isel(dim1=0)
        assert_identical(expected3, actual3)

    def test_getitem_hashable(self) -> None:
        data = create_test_data()
        data[(3, 4)] = data["var1"] + 1
        expected = data["var1"] + 1
        expected.name = (3, 4)
        assert_identical(expected, data[(3, 4)])
        with pytest.raises(KeyError, match=r"('var1', 'var2')"):
            data[("var1", "var2")]

    def test_getitem_multiple_dtype(self) -> None:
        keys = ["foo", 1]
        dataset = Dataset({key: ("dim0", range(1)) for key in keys})
        assert_identical(dataset, dataset[keys])

    def test_getitem_extra_dim_index_coord(self) -> None:
        class AnyIndex(Index):
            def should_add_coord_to_array(self, name, var, dims):
                return True

        idx = AnyIndex()
        coords = Coordinates(
            coords={
                "x": ("x", [1, 2]),
                "x_bounds": (("x", "x_bnds"), [(0.5, 1.5), (1.5, 2.5)]),
            },
            indexes={"x": idx, "x_bounds": idx},
        )

        ds = Dataset({"foo": (("x"), [1.0, 2.0])}, coords=coords)
        actual = ds["foo"]

        assert_identical(actual.coords, coords, check_default_indexes=False)
        assert "x_bnds" not in actual.dims

    def test_virtual_variables_default_coords(self) -> None:
        dataset = Dataset({"foo": ("x", range(10))})
        expected1 = DataArray(range(10), dims="x", name="x")
        actual1 = dataset["x"]
        assert_identical(expected1, actual1)
        assert isinstance(actual1.variable, IndexVariable)

        actual2 = dataset[["x", "foo"]]
        expected2 = dataset.assign_coords(x=range(10))
        assert_identical(expected2, actual2)

    def test_virtual_variables_time(self) -> None:
        # access virtual variables
        data = create_test_data()
        index = data.variables["time"].to_index()
        assert isinstance(index, pd.DatetimeIndex)
        assert_array_equal(data["time.month"].values, index.month)
        assert_array_equal(data["time.season"].values, "DJF")
        # test virtual variable math
        assert_array_equal(data["time.dayofyear"] + 1, 2 + np.arange(20))
        assert_array_equal(np.sin(data["time.dayofyear"]), np.sin(1 + np.arange(20)))
        # ensure they become coordinates
        expected = Dataset({}, {"dayofyear": data["time.dayofyear"]})
        actual = data[["time.dayofyear"]]
        assert_equal(expected, actual)
        # non-coordinate variables
        ds = Dataset({"t": ("x", pd.date_range("2000-01-01", periods=3))})
        assert (ds["t.year"] == 2000).all()

    def test_virtual_variable_same_name(self) -> None:
        # regression test for GH367
        times = pd.date_range("2000-01-01", freq="h", periods=5)
        data = Dataset({"time": times})
        actual = data["time.time"]
        expected = DataArray(times.time, [("time", times)], name="time")
        assert_identical(actual, expected)

    def test_time_season(self) -> None:
        time = xr.date_range("2000-01-01", periods=12, freq="ME", use_cftime=False)
        ds = Dataset({"t": time})
        seas = ["DJF"] * 2 + ["MAM"] * 3 + ["JJA"] * 3 + ["SON"] * 3 + ["DJF"]
        assert_array_equal(seas, ds["t.season"])

    def test_slice_virtual_variable(self) -> None:
        data = create_test_data()
        assert_equal(
            data["time.dayofyear"][:10].variable, Variable(["time"], 1 + np.arange(10))
        )
        assert_equal(data["time.dayofyear"][0].variable, Variable([], 1))

    def test_setitem(self) -> None:
        # assign a variable
        var = Variable(["dim1"], np.random.randn(8))
        data1 = create_test_data()
        data1["A"] = var
        data2 = data1.copy()
        data2["A"] = var
        assert_identical(data1, data2)
        # assign a dataset array
        dv = 2 * data2["A"]
        data1["B"] = dv.variable
        data2["B"] = dv
        assert_identical(data1, data2)
        # can't assign an ND array without dimensions
        with pytest.raises(ValueError, match=r"without explicit dimension names"):
            data2["C"] = var.values.reshape(2, 4)
        # but can assign a 1D array
        data1["C"] = var.values
        data2["C"] = ("C", var.values)
        assert_identical(data1, data2)
        # can assign a scalar
        data1["scalar"] = 0
        data2["scalar"] = ([], 0)
        assert_identical(data1, data2)
        # can't use the same dimension name as a scalar var
        with pytest.raises(ValueError, match=r"already exists as a scalar"):
            data1["newvar"] = ("scalar", [3, 4, 5])
        # can't resize a used dimension
        with pytest.raises(ValueError, match=r"conflicting dimension sizes"):
            data1["dim1"] = data1["dim1"][:5]
        # override an existing value
        data1["A"] = 3 * data2["A"]
        assert_equal(data1["A"], 3 * data2["A"])
        # can't assign a dataset to a single key
        with pytest.raises(TypeError, match="Cannot assign a Dataset to a single key"):
            data1["D"] = xr.Dataset()

        # test assignment with positional and label-based indexing
        data3 = data1[["var1", "var2"]]
        data3["var3"] = data3.var1.isel(dim1=0)
        data4 = data3.copy()
        err_msg = (
            "can only set locations defined by dictionaries from Dataset.loc. Got: a"
        )
        with pytest.raises(TypeError, match=err_msg):
            data1.loc["a"] = 0
        err_msg = r"Variables \['A', 'B', 'scalar'\] in new values not available in original dataset:"
        with pytest.raises(ValueError, match=err_msg):
            data4[{"dim2": 1}] = data1[{"dim2": 2}]
        err_msg = "Variable 'var3': indexer {'dim2': 0} not available"
        with pytest.raises(ValueError, match=err_msg):
            data1[{"dim2": 0}] = 0.0
        err_msg = "Variable 'var1': indexer {'dim2': 10} not available"
        with pytest.raises(ValueError, match=err_msg):
            data4[{"dim2": 10}] = data3[{"dim2": 2}]
        err_msg = "Variable 'var1': dimension 'dim2' appears in new values"
        with pytest.raises(KeyError, match=err_msg):
            data4[{"dim2": 2}] = data3[{"dim2": [2]}]
        err_msg = (
            "Variable 'var2': dimension order differs between original and new data"
        )
        data3["var2"] = data3["var2"].T
        with pytest.raises(ValueError, match=err_msg):
            data4[{"dim2": [2, 3]}] = data3[{"dim2": [2, 3]}]
        data3["var2"] = data3["var2"].T
        err_msg = r"cannot align objects.*not equal along these coordinates.*"
        with pytest.raises(ValueError, match=err_msg):
            data4[{"dim2": [2, 3]}] = data3[{"dim2": [2, 3, 4]}]
        err_msg = "Dataset assignment only accepts DataArrays, Datasets, and scalars."
        with pytest.raises(TypeError, match=err_msg):
            data4[{"dim2": [2, 3]}] = data3["var1"][{"dim2": [3, 4]}].values
        data5 = data4.astype(str)
        data5["var4"] = data4["var1"]
        # convert to `np.str_('a')` once `numpy<2.0` has been dropped
        err_msg = "could not convert string to float: .*'a'.*"
        with pytest.raises(ValueError, match=err_msg):
            data5[{"dim2": 1}] = "a"

        data4[{"dim2": 0}] = 0.0
        data4[{"dim2": 1}] = data3[{"dim2": 2}]
        data4.loc[{"dim2": 1.5}] = 1.0
        data4.loc[{"dim2": 2.0}] = data3.loc[{"dim2": 2.5}]
        for v, dat3 in data3.items():
            dat4 = data4[v]
            assert_array_equal(dat4[{"dim2": 0}], 0.0)
            assert_array_equal(dat4[{"dim2": 1}], dat3[{"dim2": 2}])
            assert_array_equal(dat4.loc[{"dim2": 1.5}], 1.0)
            assert_array_equal(dat4.loc[{"dim2": 2.0}], dat3.loc[{"dim2": 2.5}])
            unchanged = [1.0, 2.5, 3.0, 3.5, 4.0]
            assert_identical(
                dat4.loc[{"dim2": unchanged}], dat3.loc[{"dim2": unchanged}]
            )

    def test_setitem_pandas(self) -> None:
        ds = self.make_example_math_dataset()
        ds["x"] = np.arange(3)
        ds_copy = ds.copy()
        ds_copy["bar"] = ds["bar"].to_pandas()
        assert_equal(ds, ds_copy)

    def test_setitem_auto_align(self) -> None:
        ds = Dataset()
        ds["x"] = ("y", range(3))
        ds["y"] = 1 + np.arange(3)
        expected = Dataset({"x": ("y", range(3)), "y": 1 + np.arange(3)})
        assert_identical(ds, expected)

        ds["y"] = DataArray(range(3), dims="y")
        expected = Dataset({"x": ("y", range(3))}, {"y": range(3)})
        assert_identical(ds, expected)

        ds["x"] = DataArray([1, 2], coords=[("y", [0, 1])])
        expected = Dataset({"x": ("y", [1, 2, np.nan])}, {"y": range(3)})
        assert_identical(ds, expected)

        ds["x"] = 42
        expected = Dataset({"x": 42, "y": range(3)})
        assert_identical(ds, expected)

        ds["x"] = DataArray([4, 5, 6, 7], coords=[("y", [0, 1, 2, 3])])
        expected = Dataset({"x": ("y", [4, 5, 6])}, {"y": range(3)})
        assert_identical(ds, expected)

    def test_setitem_dimension_override(self) -> None:
        # regression test for GH-3377
        ds = xr.Dataset({"x": [0, 1, 2]})
        ds["x"] = ds["x"][:2]
        expected = Dataset({"x": [0, 1]})
        assert_identical(ds, expected)

        ds = xr.Dataset({"x": [0, 1, 2]})
        ds["x"] = np.array([0, 1])
        assert_identical(ds, expected)

        ds = xr.Dataset({"x": [0, 1, 2]})
        ds.coords["x"] = [0, 1]
        assert_identical(ds, expected)

    def test_setitem_with_coords(self) -> None:
        # Regression test for GH:2068
        ds = create_test_data()

        other = DataArray(
            np.arange(10), dims="dim3", coords={"numbers": ("dim3", np.arange(10))}
        )
        expected = ds.copy()
        expected["var3"] = other.drop_vars("numbers")
        actual = ds.copy()
        actual["var3"] = other
        assert_identical(expected, actual)
        assert "numbers" in other.coords  # should not change other

        # with alignment
        other = ds["var3"].isel(dim3=slice(1, -1))
        other["numbers"] = ("dim3", np.arange(8))
        actual = ds.copy()
        actual["var3"] = other
        assert "numbers" in other.coords  # should not change other
        expected = ds.copy()
        expected["var3"] = ds["var3"].isel(dim3=slice(1, -1))
        assert_identical(expected, actual)

        # with non-duplicate coords
        other = ds["var3"].isel(dim3=slice(1, -1))
        other["numbers"] = ("dim3", np.arange(8))
        other["position"] = ("dim3", np.arange(8))
        actual = ds.copy()
        actual["var3"] = other
        assert "position" in actual
        assert "position" in other.coords

        # assigning a coordinate-only dataarray
        actual = ds.copy()
        other = actual["numbers"]
        other[0] = 10
        actual["numbers"] = other
        assert actual["numbers"][0] == 10

        # GH: 2099
        ds = Dataset(
            {"var": ("x", [1, 2, 3])},
            coords={"x": [0, 1, 2], "z1": ("x", [1, 2, 3]), "z2": ("x", [1, 2, 3])},
        )
        ds["var"] = ds["var"] * 2
        assert np.allclose(ds["var"], [2, 4, 6])

    def test_setitem_align_new_indexes(self) -> None:
        ds = Dataset({"foo": ("x", [1, 2, 3])}, {"x": [0, 1, 2]})
        ds["bar"] = DataArray([2, 3, 4], [("x", [1, 2, 3])])
        expected = Dataset(
            {"foo": ("x", [1, 2, 3]), "bar": ("x", [np.nan, 2, 3])}, {"x": [0, 1, 2]}
        )
        assert_identical(ds, expected)

    def test_setitem_vectorized(self) -> None:
        # Regression test for GH:7030
        # Positional indexing
        da = xr.DataArray(np.r_[:120].reshape(2, 3, 4, 5), dims=["a", "b", "c", "d"])
        ds = xr.Dataset({"da": da})
        b = xr.DataArray([[0, 0], [1, 0]], dims=["u", "v"])
        c = xr.DataArray([[0, 1], [2, 3]], dims=["u", "v"])
        w = xr.DataArray([-1, -2], dims=["u"])
        index = dict(b=b, c=c)
        ds[index] = xr.Dataset({"da": w})
        assert (ds[index]["da"] == w).all()

        # Indexing with coordinates
        da = xr.DataArray(np.r_[:120].reshape(2, 3, 4, 5), dims=["a", "b", "c", "d"])
        ds = xr.Dataset({"da": da})
        ds.coords["b"] = [2, 4, 6]
        b = xr.DataArray([[2, 2], [4, 2]], dims=["u", "v"])
        c = xr.DataArray([[0, 1], [2, 3]], dims=["u", "v"])
        w = xr.DataArray([-1, -2], dims=["u"])
        index = dict(b=b, c=c)
        ds.loc[index] = xr.Dataset({"da": w}, coords={"b": ds.coords["b"]})
        assert (ds.loc[index]["da"] == w).all()

    @pytest.mark.parametrize("dtype", [str, bytes])
    def test_setitem_str_dtype(self, dtype) -> None:
        ds = xr.Dataset(coords={"x": np.array(["x", "y"], dtype=dtype)})
        # test Dataset update
        ds["foo"] = xr.DataArray(np.array([0, 0]), dims=["x"])

        assert np.issubdtype(ds.x.dtype, dtype)

    def test_setitem_using_list(self) -> None:
        # assign a list of variables
        var1 = Variable(["dim1"], np.random.randn(8))
        var2 = Variable(["dim1"], np.random.randn(8))
        actual = create_test_data()
        expected = actual.copy()
        expected["A"] = var1
        expected["B"] = var2
        actual[["A", "B"]] = [var1, var2]
        assert_identical(actual, expected)
        # assign a list of dataset arrays
        dv = 2 * expected[["A", "B"]]
        actual[["C", "D"]] = [d.variable for d in dv.data_vars.values()]
        expected[["C", "D"]] = dv
        assert_identical(actual, expected)

    @pytest.mark.parametrize(
        "var_list, data, error_regex",
        [
            (
                ["A", "B"],
                [Variable(["dim1"], np.random.randn(8))],
                r"Different lengths",
            ),
            ([], [Variable(["dim1"], np.random.randn(8))], r"Empty list of variables"),
            (["A", "B"], xr.DataArray([1, 2]), r"assign single DataArray"),
        ],
    )
    def test_setitem_using_list_errors(self, var_list, data, error_regex) -> None:
        actual = create_test_data()
        with pytest.raises(ValueError, match=error_regex):
            actual[var_list] = data

    def test_assign(self) -> None:
        ds = Dataset()
        actual = ds.assign(x=[0, 1, 2], y=2)
        expected = Dataset({"x": [0, 1, 2], "y": 2})
        assert_identical(actual, expected)
        assert list(actual.variables) == ["x", "y"]
        assert_identical(ds, Dataset())

        actual = actual.assign(y=lambda ds: ds.x**2)
        expected = Dataset({"y": ("x", [0, 1, 4]), "x": [0, 1, 2]})
        assert_identical(actual, expected)

        actual = actual.assign_coords(z=2)
        expected = Dataset({"y": ("x", [0, 1, 4])}, {"z": 2, "x": [0, 1, 2]})
        assert_identical(actual, expected)

    def test_assign_coords(self) -> None:
        ds = Dataset()

        actual = ds.assign(x=[0, 1, 2], y=2)
        actual = actual.assign_coords(x=list("abc"))
        expected = Dataset({"x": list("abc"), "y": 2})
        assert_identical(actual, expected)

        actual = ds.assign(x=[0, 1, 2], y=[2, 3])
        actual = actual.assign_coords({"y": [2.0, 3.0]})
        expected = ds.assign(x=[0, 1, 2], y=[2.0, 3.0])
        assert_identical(actual, expected)

    def test_assign_attrs(self) -> None:
        expected = Dataset(attrs=dict(a=1, b=2))
        new = Dataset()
        actual = new.assign_attrs(a=1, b=2)
        assert_identical(actual, expected)
        assert new.attrs == {}

        expected.attrs["c"] = 3
        new_actual = actual.assign_attrs({"c": 3})
        assert_identical(new_actual, expected)
        assert actual.attrs == dict(a=1, b=2)

    def test_drop_attrs(self) -> None:
        # Simple example
        ds = Dataset().assign_attrs(a=1, b=2)
        original = ds.copy()
        expected = Dataset()
        result = ds.drop_attrs()
        assert_identical(result, expected)

        # Doesn't change original
        assert_identical(ds, original)

        # Example with variables and coords with attrs, and a multiindex. (arguably
        # should have used a canonical dataset with all the features we're should
        # support...)
        var = Variable("x", [1, 2, 3], attrs=dict(x=1, y=2))
        idx = IndexVariable("y", [1, 2, 3], attrs=dict(c=1, d=2))
        mx = xr.Coordinates.from_pandas_multiindex(
            pd.MultiIndex.from_tuples([(1, 2), (3, 4)], names=["d", "e"]), "z"
        )
        ds = Dataset(dict(var1=var), coords=dict(y=idx, z=mx)).assign_attrs(a=1, b=2)
        assert ds.attrs != {}
        assert ds["var1"].attrs != {}
        assert ds["y"].attrs != {}
        assert ds.coords["y"].attrs != {}

        original = ds.copy(deep=True)
        result = ds.drop_attrs()

        assert result.attrs == {}
        assert result["var1"].attrs == {}
        assert result["y"].attrs == {}
        assert list(result.data_vars) == list(ds.data_vars)
        assert list(result.coords) == list(ds.coords)

        # Doesn't change original
        assert_identical(ds, original)
        # Specifically test that the attrs on the coords are still there. (The index
        # can't currently contain `attrs`, so we can't test those.)
        assert ds.coords["y"].attrs != {}

        # Test for deep=False
        result_shallow = ds.drop_attrs(deep=False)
        assert result_shallow.attrs == {}
        assert result_shallow["var1"].attrs != {}
        assert result_shallow["y"].attrs != {}
        assert list(result.data_vars) == list(ds.data_vars)
        assert list(result.coords) == list(ds.coords)

    def test_assign_multiindex_level(self) -> None:
        data = create_test_multiindex()
        with pytest.raises(ValueError, match=r"cannot drop or update.*corrupt.*index "):
            data.assign(level_1=range(4))
            data.assign_coords(level_1=range(4))

    def test_assign_new_multiindex(self) -> None:
        midx = pd.MultiIndex.from_arrays([["a", "a", "b", "b"], [0, 1, 0, 1]])
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")

        ds = Dataset(coords={"x": [1, 2]})
        expected = Dataset(coords=midx_coords)

        with pytest.warns(
            FutureWarning,
            match=".*`pandas.MultiIndex`.*no longer be implicitly promoted.*",
        ):
            actual = ds.assign(x=midx)
        assert_identical(actual, expected)

    @pytest.mark.parametrize("orig_coords", [{}, {"x": range(4)}])
    def test_assign_coords_new_multiindex(self, orig_coords) -> None:
        ds = Dataset(coords=orig_coords)
        midx = pd.MultiIndex.from_arrays(
            [["a", "a", "b", "b"], [0, 1, 0, 1]], names=("one", "two")
        )
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")

        expected = Dataset(coords=midx_coords)

        with pytest.warns(
            FutureWarning,
            match=".*`pandas.MultiIndex`.*no longer be implicitly promoted.*",
        ):
            actual = ds.assign_coords({"x": midx})
        assert_identical(actual, expected)

        actual = ds.assign_coords(midx_coords)
        assert_identical(actual, expected)

    def test_assign_coords_existing_multiindex(self) -> None:
        data = create_test_multiindex()
        with pytest.warns(
            FutureWarning, match=r"updating coordinate.*MultiIndex.*inconsistent"
        ):
            updated = data.assign_coords(x=range(4))
        # https://github.com/pydata/xarray/issues/7097 (coord names updated)
        assert len(updated.coords) == 1

        with pytest.warns(
            FutureWarning, match=r"updating coordinate.*MultiIndex.*inconsistent"
        ):
            updated = data.assign(x=range(4))
        # https://github.com/pydata/xarray/issues/7097 (coord names updated)
        assert len(updated.coords) == 1

    def test_assign_all_multiindex_coords(self) -> None:
        data = create_test_multiindex()
        actual = data.assign(x=range(4), level_1=range(4), level_2=range(4))
        # no error but multi-index dropped in favor of single indexes for each level
        assert (
            actual.xindexes["x"]
            is not actual.xindexes["level_1"]
            is not actual.xindexes["level_2"]
        )

    def test_assign_coords_custom_index_side_effect(self) -> None:
        # test that assigning new coordinates do not reset other dimension coord indexes
        # to default (pandas) index (https://github.com/pydata/xarray/issues/7346)
        class CustomIndex(PandasIndex):
            pass

        ds = (
            Dataset(coords={"x": [1, 2, 3]})
            .drop_indexes("x")
            .set_xindex("x", CustomIndex)
        )
        actual = ds.assign_coords(y=[4, 5, 6])
        assert isinstance(actual.xindexes["x"], CustomIndex)

    def test_assign_coords_custom_index(self) -> None:
        class CustomIndex(Index):
            pass

        coords = Coordinates(
            coords={"x": ("x", [1, 2, 3])}, indexes={"x": CustomIndex()}
        )
        ds = Dataset()
        actual = ds.assign_coords(coords)
        assert isinstance(actual.xindexes["x"], CustomIndex)

    def test_assign_coords_no_default_index(self) -> None:
        coords = Coordinates({"y": [1, 2, 3]}, indexes={})
        ds = Dataset()
        actual = ds.assign_coords(coords)
        expected = coords.to_dataset()
        assert_identical(expected, actual, check_default_indexes=False)
        assert "y" not in actual.xindexes

    def test_merge_multiindex_level(self) -> None:
        data = create_test_multiindex()

        other = Dataset({"level_1": ("x", [0, 1])})
        with pytest.raises(ValueError, match=r".*conflicting dimension sizes.*"):
            data.merge(other)

        other = Dataset({"level_1": ("x", range(4))})
        with pytest.raises(
            ValueError, match=r"unable to determine.*coordinates or not.*"
        ):
            data.merge(other)

        # `other` Dataset coordinates are ignored (bug or feature?)
        other = Dataset(coords={"level_1": ("x", range(4))})
        assert_identical(data.merge(other), data)

    def test_setitem_original_non_unique_index(self) -> None:
        # regression test for GH943
        original = Dataset({"data": ("x", np.arange(5))}, coords={"x": [0, 1, 2, 0, 1]})
        expected = Dataset({"data": ("x", np.arange(5))}, {"x": range(5)})

        actual = original.copy()
        actual["x"] = list(range(5))
        assert_identical(actual, expected)

        actual = original.copy()
        actual["x"] = ("x", list(range(5)))
        assert_identical(actual, expected)

        actual = original.copy()
        actual.coords["x"] = list(range(5))
        assert_identical(actual, expected)

    def test_setitem_both_non_unique_index(self) -> None:
        # regression test for GH956
        names = ["joaquin", "manolo", "joaquin"]
        values = np.random.randint(0, 256, (3, 4, 4))
        array = DataArray(
            values, dims=["name", "row", "column"], coords=[names, range(4), range(4)]
        )
        expected = Dataset({"first": array, "second": array})
        actual = array.rename("first").to_dataset()
        actual["second"] = array
        assert_identical(expected, actual)

    def test_setitem_multiindex_level(self) -> None:
        data = create_test_multiindex()
        with pytest.raises(
            ValueError, match=r"cannot set or update variable.*corrupt.*index "
        ):
            data["level_1"] = range(4)

    def test_delitem(self) -> None:
        data = create_test_data()
        all_items = set(data.variables)
        assert set(data.variables) == all_items
        del data["var1"]
        assert set(data.variables) == all_items - {"var1"}
        del data["numbers"]
        assert set(data.variables) == all_items - {"var1", "numbers"}
        assert "numbers" not in data.coords

        expected = Dataset()
        actual = Dataset({"y": ("x", [1, 2])})
        del actual["y"]
        assert_identical(expected, actual)

    def test_delitem_multiindex_level(self) -> None:
        data = create_test_multiindex()
        with pytest.raises(
            ValueError, match=r"cannot remove coordinate.*corrupt.*index "
        ):
            del data["level_1"]

    def test_squeeze(self) -> None:
        data = Dataset({"foo": (["x", "y", "z"], [[[1], [2]]])})
        test_args: list[list] = [[], [["x"]], [["x", "z"]]]
        for args in test_args:

            def get_args(args, v):
                return [set(args[0]) & set(v.dims)] if args else []

            expected = Dataset(
                {k: v.squeeze(*get_args(args, v)) for k, v in data.variables.items()}
            )
            expected = expected.set_coords(data.coords)
            assert_identical(expected, data.squeeze(*args))
        # invalid squeeze
        with pytest.raises(ValueError, match=r"cannot select a dimension"):
            data.squeeze("y")

    def test_squeeze_drop(self) -> None:
        data = Dataset({"foo": ("x", [1])}, {"x": [0]})
        expected = Dataset({"foo": 1})
        selected = data.squeeze(drop=True)
        assert_identical(expected, selected)

        expected = Dataset({"foo": 1}, {"x": 0})
        selected = data.squeeze(drop=False)
        assert_identical(expected, selected)

        data = Dataset({"foo": (("x", "y"), [[1]])}, {"x": [0], "y": [0]})
        expected = Dataset({"foo": 1})
        selected = data.squeeze(drop=True)
        assert_identical(expected, selected)

        expected = Dataset({"foo": ("x", [1])}, {"x": [0]})
        selected = data.squeeze(dim="y", drop=True)
        assert_identical(expected, selected)

        data = Dataset({"foo": (("x",), [])}, {"x": []})
        selected = data.squeeze(drop=True)
        assert_identical(data, selected)

    def test_to_dataarray(self) -> None:
        ds = Dataset(
            {"a": 1, "b": ("x", [1, 2, 3])},
            coords={"c": 42},
            attrs={"Conventions": "None"},
        )
        data = [[1, 1, 1], [1, 2, 3]]
        coords = {"c": 42, "variable": ["a", "b"]}
        dims = ("variable", "x")
        expected = DataArray(data, coords, dims, attrs=ds.attrs)
        actual = ds.to_dataarray()
        assert_identical(expected, actual)

        actual = ds.to_dataarray("abc", name="foo")
        expected = expected.rename({"variable": "abc"}).rename("foo")
        assert_identical(expected, actual)

    def test_to_and_from_dataframe(self) -> None:
        x = np.random.randn(10)
        y = np.random.randn(10)
        t = list("abcdefghij")
        cat = pd.Categorical(["a", "b"] * 5)
        ds = Dataset({"a": ("t", x), "b": ("t", y), "t": ("t", t), "cat": ("t", cat)})
        expected = pd.DataFrame(
            np.array([x, y]).T, columns=["a", "b"], index=pd.Index(t, name="t")
        )
        expected["cat"] = cat
        actual = ds.to_dataframe()
        # use the .equals method to check all DataFrame metadata
        assert expected.equals(actual), (expected, actual)

        # verify coords are included
        actual = ds.set_coords("b").to_dataframe()
        assert expected.equals(actual), (expected, actual)

        # check roundtrip
        assert_identical(ds, Dataset.from_dataframe(actual))
        assert isinstance(ds["cat"].variable.data.dtype, pd.CategoricalDtype)
        # test a case with a MultiIndex
        w = np.random.randn(2, 3)
        cat = pd.Categorical(["a", "a", "c"])
        ds = Dataset({"w": (("x", "y"), w), "cat": ("y", cat)})
        ds["y"] = ("y", list("abc"))
        exp_index = pd.MultiIndex.from_arrays(
            [[0, 0, 0, 1, 1, 1], ["a", "b", "c", "a", "b", "c"]], names=["x", "y"]
        )
        expected = pd.DataFrame(
            {"w": w.reshape(-1), "cat": pd.Categorical(["a", "a", "c", "a", "a", "c"])},
            index=exp_index,
        )
        actual = ds.to_dataframe()
        assert expected.equals(actual)

        # check roundtrip
        # from_dataframe attempts to broadcast across because it doesn't know better, so cat must be converted
        ds["cat"] = (("x", "y"), np.stack((ds["cat"].to_numpy(), ds["cat"].to_numpy())))
        assert_identical(ds.assign_coords(x=[0, 1]), Dataset.from_dataframe(actual))

        # Check multiindex reordering
        new_order = ["x", "y"]
        # revert broadcasting fix above for 1d arrays
        ds["cat"] = ("y", cat)
        actual = ds.to_dataframe(dim_order=new_order)
        assert expected.equals(actual)

        new_order = ["y", "x"]
        exp_index = pd.MultiIndex.from_arrays(
            [["a", "a", "b", "b", "c", "c"], [0, 1, 0, 1, 0, 1]], names=["y", "x"]
        )
        expected = pd.DataFrame(
            {
                "w": w.transpose().reshape(-1),
                "cat": pd.Categorical(["a", "a", "a", "a", "c", "c"]),
            },
            index=exp_index,
        )
        actual = ds.to_dataframe(dim_order=new_order)
        assert expected.equals(actual)

        invalid_order = ["x"]
        with pytest.raises(
            ValueError, match="does not match the set of dimensions of this"
        ):
            ds.to_dataframe(dim_order=invalid_order)

        invalid_order = ["x", "z"]
        with pytest.raises(
            ValueError, match="does not match the set of dimensions of this"
        ):
            ds.to_dataframe(dim_order=invalid_order)

        # check pathological cases
        df = pd.DataFrame([1])
        actual_ds = Dataset.from_dataframe(df)
        expected_ds = Dataset({0: ("index", [1])}, {"index": [0]})
        assert_identical(expected_ds, actual_ds)

        df = pd.DataFrame()
        actual_ds = Dataset.from_dataframe(df)
        expected_ds = Dataset(coords={"index": []})
        assert_identical(expected_ds, actual_ds)

        # GH697
        df = pd.DataFrame({"A": []})
        actual_ds = Dataset.from_dataframe(df)
        expected_ds = Dataset({"A": DataArray([], dims=("index",))}, {"index": []})
        assert_identical(expected_ds, actual_ds)

        # regression test for GH278
        # use int64 to ensure consistent results for the pandas .equals method
        # on windows (which requires the same dtype)
        ds = Dataset({"x": pd.Index(["bar"]), "a": ("y", np.array([1], "int64"))}).isel(
            x=0
        )
        # use .loc to ensure consistent results on Python 3
        actual = ds.to_dataframe().loc[:, ["a", "x"]]
        expected = pd.DataFrame(
            [[1, "bar"]], index=pd.Index([0], name="y"), columns=["a", "x"]
        )
        assert expected.equals(actual), (expected, actual)

        ds = Dataset({"x": np.array([0], "int64"), "y": np.array([1], "int64")})
        actual = ds.to_dataframe()
        idx = pd.MultiIndex.from_arrays([[0], [1]], names=["x", "y"])
        expected = pd.DataFrame([[]], index=idx)
        assert expected.equals(actual), (expected, actual)

    def test_from_dataframe_categorical_dtype_index(self) -> None:
        cat = pd.CategoricalIndex(list("abcd"))
        df = pd.DataFrame({"f": [0, 1, 2, 3]}, index=cat)
        ds = df.to_xarray()
        restored = ds.to_dataframe()
        df.index.name = (
            "index"  # restored gets the name because it has the coord with the name
        )
        pd.testing.assert_frame_equal(df, restored)

    def test_from_dataframe_categorical_index(self) -> None:
        cat = pd.CategoricalDtype(
            categories=["foo", "bar", "baz", "qux", "quux", "corge"]
        )
        i1 = pd.Series(["foo", "bar", "foo"], dtype=cat)
        i2 = pd.Series(["bar", "bar", "baz"], dtype=cat)

        df = pd.DataFrame({"i1": i1, "i2": i2, "values": [1, 2, 3]})
        ds = df.set_index("i1").to_xarray()
        assert len(ds["i1"]) == 3

        ds = df.set_index(["i1", "i2"]).to_xarray()
        assert len(ds["i1"]) == 2
        assert len(ds["i2"]) == 2

    def test_from_dataframe_categorical_index_string_categories(self) -> None:
        cat = pd.CategoricalIndex(
            pd.Categorical.from_codes(
                np.array([1, 1, 0, 2], dtype=np.int64),  # type: ignore[arg-type]
                categories=pd.Index(["foo", "bar", "baz"], dtype="string"),
            )
        )
        ser = pd.Series(1, index=cat)
        ds = ser.to_xarray()
        assert ds.coords.dtypes["index"] == ser.index.dtype

    @requires_sparse
    def test_from_dataframe_sparse(self) -> None:
        import sparse

        df_base = pd.DataFrame(
            {"x": range(10), "y": list("abcdefghij"), "z": np.arange(0, 100, 10)}
        )

        ds_sparse = Dataset.from_dataframe(df_base.set_index("x"), sparse=True)
        ds_dense = Dataset.from_dataframe(df_base.set_index("x"), sparse=False)
        assert isinstance(ds_sparse["y"].data, sparse.COO)
        assert isinstance(ds_sparse["z"].data, sparse.COO)
        ds_sparse["y"].data = ds_sparse["y"].data.todense()
        ds_sparse["z"].data = ds_sparse["z"].data.todense()
        assert_identical(ds_dense, ds_sparse)

        ds_sparse = Dataset.from_dataframe(df_base.set_index(["x", "y"]), sparse=True)
        ds_dense = Dataset.from_dataframe(df_base.set_index(["x", "y"]), sparse=False)
        assert isinstance(ds_sparse["z"].data, sparse.COO)
        ds_sparse["z"].data = ds_sparse["z"].data.todense()
        assert_identical(ds_dense, ds_sparse)

    def test_to_and_from_empty_dataframe(self) -> None:
        # GH697
        expected = pd.DataFrame({"foo": []})
        ds = Dataset.from_dataframe(expected)
        assert len(ds["foo"]) == 0
        actual = ds.to_dataframe()
        assert len(actual) == 0
        assert expected.equals(actual)

    def test_from_dataframe_multiindex(self) -> None:
        index = pd.MultiIndex.from_product([["a", "b"], [1, 2, 3]], names=["x", "y"])
        df = pd.DataFrame({"z": np.arange(6)}, index=index)

        expected = Dataset(
            {"z": (("x", "y"), [[0, 1, 2], [3, 4, 5]])},
            coords={"x": ["a", "b"], "y": [1, 2, 3]},
        )
        actual = Dataset.from_dataframe(df)
        assert_identical(actual, expected)

        df2 = df.iloc[[3, 2, 1, 0, 4, 5], :]
        actual = Dataset.from_dataframe(df2)
        assert_identical(actual, expected)

        df3 = df.iloc[:4, :]
        expected3 = Dataset(
            {"z": (("x", "y"), [[0, 1, 2], [3, np.nan, np.nan]])},
            coords={"x": ["a", "b"], "y": [1, 2, 3]},
        )
        actual = Dataset.from_dataframe(df3)
        assert_identical(actual, expected3)

        df_nonunique = df.iloc[[0, 0], :]
        with pytest.raises(ValueError, match=r"non-unique MultiIndex"):
            Dataset.from_dataframe(df_nonunique)

    def test_from_dataframe_unsorted_levels(self) -> None:
        # regression test for GH-4186
        index = pd.MultiIndex(
            levels=[["b", "a"], ["foo"]], codes=[[0, 1], [0, 0]], names=["lev1", "lev2"]
        )
        df = pd.DataFrame({"c1": [0, 2], "c2": [1, 3]}, index=index)
        expected = Dataset(
            {
                "c1": (("lev1", "lev2"), [[0], [2]]),
                "c2": (("lev1", "lev2"), [[1], [3]]),
            },
            coords={"lev1": ["b", "a"], "lev2": ["foo"]},
        )
        actual = Dataset.from_dataframe(df)
        assert_identical(actual, expected)

    def test_from_dataframe_non_unique_columns(self) -> None:
        # regression test for GH449
        df = pd.DataFrame(np.zeros((2, 2)))
        df.columns = ["foo", "foo"]  # type: ignore[assignment,unused-ignore]
        with pytest.raises(ValueError, match=r"non-unique columns"):
            Dataset.from_dataframe(df)

    def test_convert_dataframe_with_many_types_and_multiindex(self) -> None:
        # regression test for GH737
        df = pd.DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="America/New_York"),
            }
        )
        df.index = pd.MultiIndex.from_product([["a"], range(3)], names=["one", "two"])
        roundtripped = Dataset.from_dataframe(df).to_dataframe()
        # we can't do perfectly, but we should be at least as faithful as
        # np.asarray
        expected = df.apply(np.asarray)
        assert roundtripped.equals(expected)

    @pytest.mark.parametrize("encoding", [True, False])
    @pytest.mark.parametrize("data", [True, "list", "array"])
    def test_to_and_from_dict(
        self, encoding: bool, data: bool | Literal["list", "array"]
    ) -> None:
        # <xarray.Dataset>
        # Dimensions:  (t: 10)
        # Coordinates:
        #   * t        (t) <U1 'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j'
        # Data variables:
        #     a        (t) float64 0.6916 -1.056 -1.163 0.9792 -0.7865 ...
        #     b        (t) float64 1.32 0.1954 1.91 1.39 0.519 -0.2772 ...
        x = np.random.randn(10)
        y = np.random.randn(10)
        t = list("abcdefghij")
        ds = Dataset({"a": ("t", x), "b": ("t", y), "t": ("t", t)})
        expected: dict[str, dict[str, Any]] = {
            "coords": {"t": {"dims": ("t",), "data": t, "attrs": {}}},
            "attrs": {},
            "dims": {"t": 10},
            "data_vars": {
                "a": {"dims": ("t",), "data": x.tolist(), "attrs": {}},
                "b": {"dims": ("t",), "data": y.tolist(), "attrs": {}},
            },
        }
        if encoding:
            ds.t.encoding.update({"foo": "bar"})
            expected["encoding"] = {}
            expected["coords"]["t"]["encoding"] = ds.t.encoding
            for vvs in ["a", "b"]:
                expected["data_vars"][vvs]["encoding"] = {}

        actual = ds.to_dict(data=data, encoding=encoding)

        # check that they are identical
        np.testing.assert_equal(expected, actual)

        # check roundtrip
        ds_rt = Dataset.from_dict(actual)
        assert_identical(ds, ds_rt)
        if encoding:
            assert set(ds_rt.variables) == set(ds.variables)
            for vv in ds.variables:
                np.testing.assert_equal(ds_rt[vv].encoding, ds[vv].encoding)

        # check the data=False option
        expected_no_data = expected.copy()
        del expected_no_data["coords"]["t"]["data"]
        del expected_no_data["data_vars"]["a"]["data"]
        del expected_no_data["data_vars"]["b"]["data"]
        endiantype = "<U1" if sys.byteorder == "little" else ">U1"
        expected_no_data["coords"]["t"].update({"dtype": endiantype, "shape": (10,)})
        expected_no_data["data_vars"]["a"].update({"dtype": "float64", "shape": (10,)})
        expected_no_data["data_vars"]["b"].update({"dtype": "float64", "shape": (10,)})
        actual_no_data = ds.to_dict(data=False, encoding=encoding)
        assert expected_no_data == actual_no_data

        # verify coords are included roundtrip
        expected_ds = ds.set_coords("b")
        actual2 = Dataset.from_dict(expected_ds.to_dict(data=data, encoding=encoding))

        assert_identical(expected_ds, actual2)
        if encoding:
            assert set(expected_ds.variables) == set(actual2.variables)
            for vv in ds.variables:
                np.testing.assert_equal(expected_ds[vv].encoding, actual2[vv].encoding)

        # test some incomplete dicts:
        # this one has no attrs field, the dims are strings, and x, y are
        # np.arrays

        d = {
            "coords": {"t": {"dims": "t", "data": t}},
            "dims": "t",
            "data_vars": {"a": {"dims": "t", "data": x}, "b": {"dims": "t", "data": y}},
        }
        assert_identical(ds, Dataset.from_dict(d))

        # this is kind of a flattened version with no coords, or data_vars
        d = {
            "a": {"dims": "t", "data": x},
            "t": {"data": t, "dims": "t"},
            "b": {"dims": "t", "data": y},
        }
        assert_identical(ds, Dataset.from_dict(d))

        # this one is missing some necessary information
        d = {
            "a": {"data": x},
            "t": {"data": t, "dims": "t"},
            "b": {"dims": "t", "data": y},
        }
        with pytest.raises(
            ValueError, match=r"cannot convert dict without the key 'dims'"
        ):
            Dataset.from_dict(d)

    def test_to_and_from_dict_with_time_dim(self) -> None:
        x = np.random.randn(10, 3)
        y = np.random.randn(10, 3)
        t = pd.date_range("20130101", periods=10)
        lat = [77.7, 83.2, 76]
        ds = Dataset(
            {
                "a": (["t", "lat"], x),
                "b": (["t", "lat"], y),
                "t": ("t", t),
                "lat": ("lat", lat),
            }
        )
        roundtripped = Dataset.from_dict(ds.to_dict())
        assert_identical(ds, roundtripped)

    @pytest.mark.parametrize("data", [True, "list", "array"])
    def test_to_and_from_dict_with_nan_nat(
        self, data: bool | Literal["list", "array"]
    ) -> None:
        x = np.random.randn(10, 3)
        y = np.random.randn(10, 3)
        y[2] = np.nan
        t = pd.Series(pd.date_range("20130101", periods=10))
        t[2] = np.nan

        lat = [77.7, 83.2, 76]
        ds = Dataset(
            {
                "a": (["t", "lat"], x),
                "b": (["t", "lat"], y),
                "t": ("t", t),
                "lat": ("lat", lat),
            }
        )
        roundtripped = Dataset.from_dict(ds.to_dict(data=data))
        assert_identical(ds, roundtripped)

    def test_to_dict_with_numpy_attrs(self) -> None:
        # this doesn't need to roundtrip
        x = np.random.randn(10)
        y = np.random.randn(10)
        t = list("abcdefghij")
        attrs = {
            "created": np.float64(1998),
            "coords": np.array([37, -110.1, 100]),
            "maintainer": "bar",
        }
        ds = Dataset({"a": ("t", x, attrs), "b": ("t", y, attrs), "t": ("t", t)})
        expected_attrs = {
            "created": attrs["created"].item(),  # type: ignore[attr-defined]
            "coords": attrs["coords"].tolist(),  # type: ignore[attr-defined]
            "maintainer": "bar",
        }
        actual = ds.to_dict()

        # check that they are identical
        assert expected_attrs == actual["data_vars"]["a"]["attrs"]

    def test_pickle(self) -> None:
        data = create_test_data()
        roundtripped = pickle.loads(pickle.dumps(data))
        assert_identical(data, roundtripped)
        # regression test for #167:
        assert data.sizes == roundtripped.sizes

    def test_lazy_load(self) -> None:
        store = InaccessibleVariableDataStore()
        create_test_data().dump_to_store(store)

        for decode_cf in [True, False]:
            ds = open_dataset(store, decode_cf=decode_cf)
            with pytest.raises(UnexpectedDataAccess):
                ds.load()
            with pytest.raises(UnexpectedDataAccess):
                _ = ds["var1"].values

            # these should not raise UnexpectedDataAccess:
            ds.isel(time=10)
            ds.isel(time=slice(10), dim1=[0]).isel(dim1=0, dim2=-1)

    def test_lazy_load_duck_array(self) -> None:
        store = AccessibleAsDuckArrayDataStore()
        create_test_data().dump_to_store(store)

        for decode_cf in [True, False]:
            ds = open_dataset(store, decode_cf=decode_cf)
            with pytest.raises(UnexpectedDataAccess):
                _ = ds["var1"].values

            # these should not raise UnexpectedDataAccess:
            _ = ds.var1.data
            ds.isel(time=10)
            ds.isel(time=slice(10), dim1=[0]).isel(dim1=0, dim2=-1)
            repr(ds)

            # preserve the duck array type and don't cast to array
            assert isinstance(ds["var1"].load().data, DuckArrayWrapper)
            assert isinstance(
                ds["var1"].isel(dim2=0, dim1=0).load().data, DuckArrayWrapper
            )

            ds.close()

    def test_dropna(self) -> None:
        x = np.random.randn(4, 4)
        x[::2, 0] = np.nan
        y = np.random.randn(4)
        y[-1] = np.nan
        ds = Dataset({"foo": (("a", "b"), x), "bar": (("b", y))})

        expected = ds.isel(a=slice(1, None, 2))
        actual = ds.dropna("a")
        assert_identical(actual, expected)

        expected = ds.isel(b=slice(1, 3))
        actual = ds.dropna("b")
        assert_identical(actual, expected)

        actual = ds.dropna("b", subset=["foo", "bar"])
        assert_identical(actual, expected)

        expected = ds.isel(b=slice(1, None))
        actual = ds.dropna("b", subset=["foo"])
        assert_identical(actual, expected)

        expected = ds.isel(b=slice(3))
        actual = ds.dropna("b", subset=["bar"])
        assert_identical(actual, expected)

        actual = ds.dropna("a", subset=[])
        assert_identical(actual, ds)

        actual = ds.dropna("a", subset=["bar"])
        assert_identical(actual, ds)

        actual = ds.dropna("a", how="all")
        assert_identical(actual, ds)

        actual = ds.dropna("b", how="all", subset=["bar"])
        expected = ds.isel(b=[0, 1, 2])
        assert_identical(actual, expected)

        actual = ds.dropna("b", thresh=1, subset=["bar"])
        assert_identical(actual, expected)

        actual = ds.dropna("b", thresh=2)
        assert_identical(actual, ds)

        actual = ds.dropna("b", thresh=4)
        expected = ds.isel(b=[1, 2, 3])
        assert_identical(actual, expected)

        actual = ds.dropna("a", thresh=3)
        expected = ds.isel(a=[1, 3])
        assert_identical(actual, ds)

        with pytest.raises(
            ValueError,
            match=r"'foo' not found in data dimensions \('a', 'b'\)",
        ):
            ds.dropna("foo")
        with pytest.raises(ValueError, match=r"invalid how"):
            ds.dropna("a", how="somehow")  # type: ignore[arg-type]
        with pytest.raises(TypeError, match=r"must specify how or thresh"):
            ds.dropna("a", how=None)  # type: ignore[arg-type]

    def test_fillna(self) -> None:
        ds = Dataset({"a": ("x", [np.nan, 1, np.nan, 3])}, {"x": [0, 1, 2, 3]})

        # fill with -1
        actual1 = ds.fillna(-1)
        expected = Dataset({"a": ("x", [-1, 1, -1, 3])}, {"x": [0, 1, 2, 3]})
        assert_identical(expected, actual1)

        actual2 = ds.fillna({"a": -1})
        assert_identical(expected, actual2)

        other = Dataset({"a": -1})
        actual3 = ds.fillna(other)
        assert_identical(expected, actual3)

        actual4 = ds.fillna({"a": other.a})
        assert_identical(expected, actual4)

        # fill with range(4)
        b = DataArray(range(4), coords=[("x", range(4))])
        actual5 = ds.fillna(b)
        expected = b.rename("a").to_dataset()
        assert_identical(expected, actual5)

        actual6 = ds.fillna(expected)
        assert_identical(expected, actual6)

        actual7 = ds.fillna(np.arange(4))
        assert_identical(expected, actual7)

        actual8 = ds.fillna(b[:3])
        assert_identical(expected, actual8)

        # okay to only include some data variables
        ds["b"] = np.nan
        actual9 = ds.fillna({"a": -1})
        expected = Dataset(
            {"a": ("x", [-1, 1, -1, 3]), "b": np.nan}, {"x": [0, 1, 2, 3]}
        )
        assert_identical(expected, actual9)

        # but new data variables is not okay
        with pytest.raises(ValueError, match=r"must be contained"):
            ds.fillna({"x": 0})

        # empty argument should be OK
        result1 = ds.fillna({})
        assert_identical(ds, result1)

        result2 = ds.fillna(Dataset(coords={"c": 42}))
        expected = ds.assign_coords(c=42)
        assert_identical(expected, result2)

        da = DataArray(range(5), name="a", attrs={"attr": "da"})
        actual10 = da.fillna(1)
        assert actual10.name == "a"
        assert actual10.attrs == da.attrs

        ds = Dataset({"a": da}, attrs={"attr": "ds"})
        actual11 = ds.fillna({"a": 1})
        assert actual11.attrs == ds.attrs
        assert actual11.a.name == "a"
        assert actual11.a.attrs == ds.a.attrs

    @pytest.mark.parametrize(
        "func", [lambda x: x.clip(0, 1), lambda x: np.float64(1.0) * x, np.abs, abs]
    )
    def test_propagate_attrs(self, func) -> None:
        da = DataArray(range(5), name="a", attrs={"attr": "da"})
        ds = Dataset({"a": da}, attrs={"attr": "ds"})

        # test defaults
        assert func(ds).attrs == ds.attrs
        with set_options(keep_attrs=False):
            assert func(ds).attrs != ds.attrs
            assert func(ds).a.attrs != ds.a.attrs

        with set_options(keep_attrs=False):
            assert func(ds).attrs != ds.attrs
            assert func(ds).a.attrs != ds.a.attrs

        with set_options(keep_attrs=True):
            assert func(ds).attrs == ds.attrs
            assert func(ds).a.attrs == ds.a.attrs

    def test_where(self) -> None:
        ds = Dataset({"a": ("x", range(5))})
        expected1 = Dataset({"a": ("x", [np.nan, np.nan, 2, 3, 4])})
        actual1 = ds.where(ds > 1)
        assert_identical(expected1, actual1)

        actual2 = ds.where(ds.a > 1)
        assert_identical(expected1, actual2)

        actual3 = ds.where(ds.a.values > 1)
        assert_identical(expected1, actual3)

        actual4 = ds.where(True)
        assert_identical(ds, actual4)

        expected5 = ds.copy(deep=True)
        expected5["a"].values = np.array([np.nan] * 5)
        actual5 = ds.where(False)
        assert_identical(expected5, actual5)

        # 2d
        ds = Dataset({"a": (("x", "y"), [[0, 1], [2, 3]])})
        expected6 = Dataset({"a": (("x", "y"), [[np.nan, 1], [2, 3]])})
        actual6 = ds.where(ds > 0)
        assert_identical(expected6, actual6)

        # attrs
        da = DataArray(range(5), name="a", attrs={"attr": "da"})
        actual7 = da.where(da.values > 1)
        assert actual7.name == "a"
        assert actual7.attrs == da.attrs

        ds = Dataset({"a": da}, attrs={"attr": "ds"})
        actual8 = ds.where(ds > 0)
        assert actual8.attrs == ds.attrs
        assert actual8.a.name == "a"
        assert actual8.a.attrs == ds.a.attrs

        # lambda
        ds = Dataset({"a": ("x", range(5))})
        expected9 = Dataset({"a": ("x", [np.nan, np.nan, 2, 3, 4])})
        actual9 = ds.where(lambda x: x > 1)
        assert_identical(expected9, actual9)

    def test_where_other(self) -> None:
        ds = Dataset({"a": ("x", range(5))}, {"x": range(5)})
        expected = Dataset({"a": ("x", [-1, -1, 2, 3, 4])}, {"x": range(5)})
        actual = ds.where(ds > 1, -1)
        assert_equal(expected, actual)
        assert actual.a.dtype == int

        actual = ds.where(lambda x: x > 1, -1)
        assert_equal(expected, actual)

        actual = ds.where(ds > 1, other=-1, drop=True)
        expected_nodrop = ds.where(ds > 1, -1)
        _, expected = xr.align(actual, expected_nodrop, join="left")
        assert_equal(actual, expected)
        assert actual.a.dtype == int

        with pytest.raises(ValueError, match=r"cannot align .* are not equal"):
            ds.where(ds > 1, ds.isel(x=slice(3)))

        with pytest.raises(ValueError, match=r"exact match required"):
            ds.where(ds > 1, ds.assign(b=2))

    def test_where_drop(self) -> None:
        # if drop=True

        # 1d
        # data array case
        array = DataArray(range(5), coords=[range(5)], dims=["x"])
        expected1 = DataArray(range(5)[2:], coords=[range(5)[2:]], dims=["x"])
        actual1 = array.where(array > 1, drop=True)
        assert_identical(expected1, actual1)

        # dataset case
        ds = Dataset({"a": array})
        expected2 = Dataset({"a": expected1})

        actual2 = ds.where(ds > 1, drop=True)
        assert_identical(expected2, actual2)

        actual3 = ds.where(ds.a > 1, drop=True)
        assert_identical(expected2, actual3)

        with pytest.raises(TypeError, match=r"must be a"):
            ds.where(np.arange(5) > 1, drop=True)

        # 1d with odd coordinates
        array = DataArray(
            np.array([2, 7, 1, 8, 3]), coords=[np.array([3, 1, 4, 5, 9])], dims=["x"]
        )
        expected4 = DataArray(
            np.array([7, 8, 3]), coords=[np.array([1, 5, 9])], dims=["x"]
        )
        actual4 = array.where(array > 2, drop=True)
        assert_identical(expected4, actual4)

        # 1d multiple variables
        ds = Dataset({"a": (("x"), [0, 1, 2, 3]), "b": (("x"), [4, 5, 6, 7])})
        expected5 = Dataset(
            {"a": (("x"), [np.nan, 1, 2, 3]), "b": (("x"), [4, 5, 6, np.nan])}
        )
        actual5 = ds.where((ds > 0) & (ds < 7), drop=True)
        assert_identical(expected5, actual5)

        # 2d
        ds = Dataset({"a": (("x", "y"), [[0, 1], [2, 3]])})
        expected6 = Dataset({"a": (("x", "y"), [[np.nan, 1], [2, 3]])})
        actual6 = ds.where(ds > 0, drop=True)
        assert_identical(expected6, actual6)

        # 2d with odd coordinates
        ds = Dataset(
            {"a": (("x", "y"), [[0, 1], [2, 3]])},
            coords={
                "x": [4, 3],
                "y": [1, 2],
                "z": (["x", "y"], [[np.e, np.pi], [np.pi * np.e, np.pi * 3]]),
            },
        )
        expected7 = Dataset(
            {"a": (("x", "y"), [[3]])},
            coords={"x": [3], "y": [2], "z": (["x", "y"], [[np.pi * 3]])},
        )
        actual7 = ds.where(ds > 2, drop=True)
        assert_identical(expected7, actual7)

        # 2d multiple variables
        ds = Dataset(
            {"a": (("x", "y"), [[0, 1], [2, 3]]), "b": (("x", "y"), [[4, 5], [6, 7]])}
        )
        expected8 = Dataset(
            {
                "a": (("x", "y"), [[np.nan, 1], [2, 3]]),
                "b": (("x", "y"), [[4, 5], [6, 7]]),
            }
        )
        actual8 = ds.where(ds > 0, drop=True)
        assert_identical(expected8, actual8)

        # mixed dimensions: PR#6690, Issue#6227
        ds = xr.Dataset(
            {
                "a": ("x", [1, 2, 3]),
                "b": ("y", [2, 3, 4]),
                "c": (("x", "y"), np.arange(9).reshape((3, 3))),
            }
        )
        expected9 = xr.Dataset(
            {
                "a": ("x", [np.nan, 3]),
                "b": ("y", [np.nan, 3, 4]),
                "c": (("x", "y"), np.arange(3.0, 9.0).reshape((2, 3))),
            }
        )
        actual9 = ds.where(ds > 2, drop=True)
        assert actual9.sizes["x"] == 2
        assert_identical(expected9, actual9)

    def test_where_drop_empty(self) -> None:
        # regression test for GH1341
        array = DataArray(np.random.rand(100, 10), dims=["nCells", "nVertLevels"])
        mask = DataArray(np.zeros((100,), dtype="bool"), dims="nCells")
        actual = array.where(mask, drop=True)
        expected = DataArray(np.zeros((0, 10)), dims=["nCells", "nVertLevels"])
        assert_identical(expected, actual)

    def test_where_drop_no_indexes(self) -> None:
        ds = Dataset({"foo": ("x", [0.0, 1.0])})
        expected = Dataset({"foo": ("x", [1.0])})
        actual = ds.where(ds == 1, drop=True)
        assert_identical(expected, actual)

    def test_reduce(self) -> None:
        data = create_test_data()

        assert len(data.mean().coords) == 0

        actual = data.max()
        expected = Dataset({k: v.max() for k, v in data.data_vars.items()})
        assert_equal(expected, actual)

        assert_equal(data.min(dim=["dim1"]), data.min(dim="dim1"))

        for reduct, expected_dims in [
            ("dim2", ["dim3", "time", "dim1"]),
            (["dim2", "time"], ["dim3", "dim1"]),
            (("dim2", "time"), ["dim3", "dim1"]),
            ((), ["dim2", "dim3", "time", "dim1"]),
        ]:
            actual_dims = list(data.min(dim=reduct).dims)
            assert actual_dims == expected_dims

        assert_equal(data.mean(dim=[]), data)

        with pytest.raises(ValueError):
            data.mean(axis=0)

    def test_reduce_coords(self) -> None:
        # regression test for GH1470
        data = xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"b": 4})
        expected = xr.Dataset({"a": 2}, coords={"b": 4})
        actual = data.mean("x")
        assert_identical(actual, expected)

        # should be consistent
        actual = data["a"].mean("x").to_dataset()
        assert_identical(actual, expected)

    def test_mean_uint_dtype(self) -> None:
        data = xr.Dataset(
            {
                "a": (("x", "y"), np.arange(6).reshape(3, 2).astype("uint")),
                "b": (("x",), np.array([0.1, 0.2, np.nan])),
            }
        )
        actual = data.mean("x", skipna=True)
        expected = xr.Dataset(
            {"a": data["a"].mean("x"), "b": data["b"].mean("x", skipna=True)}
        )
        assert_identical(actual, expected)

    def test_reduce_bad_dim(self) -> None:
        data = create_test_data()
        with pytest.raises(
            ValueError,
            match=re.escape("Dimension(s) 'bad_dim' do not exist"),
        ):
            data.mean(dim="bad_dim")

    def test_reduce_cumsum(self) -> None:
        data = xr.Dataset(
            {"a": 1, "b": ("x", [1, 2]), "c": (("x", "y"), [[np.nan, 3], [0, 4]])}
        )
        assert_identical(data.fillna(0), data.cumsum("y"))

        expected = xr.Dataset(
            {"a": 1, "b": ("x", [1, 3]), "c": (("x", "y"), [[0, 3], [0, 7]])}
        )
        assert_identical(expected, data.cumsum())

    @pytest.mark.parametrize(
        "reduct, expected",
        [
            ("dim1", ["dim2", "dim3", "time", "dim1"]),
            ("dim2", ["dim3", "time", "dim1", "dim2"]),
            ("dim3", ["dim2", "time", "dim1", "dim3"]),
            ("time", ["dim2", "dim3", "dim1"]),
        ],
    )
    @pytest.mark.parametrize("func", ["cumsum", "cumprod"])
    def test_reduce_cumsum_test_dims(self, reduct, expected, func) -> None:
        data = create_test_data()
        with pytest.raises(
            ValueError,
            match=re.escape("Dimension(s) 'bad_dim' do not exist"),
        ):
            getattr(data, func)(dim="bad_dim")

        # ensure dimensions are correct
        actual = getattr(data, func)(dim=reduct).dims
        assert list(actual) == expected

    def test_reduce_non_numeric(self) -> None:
        data1 = create_test_data(seed=44, use_extension_array=True)
        data2 = create_test_data(seed=44)
        add_vars = {"var6": ["dim1", "dim2"], "var7": ["dim1"]}
        for v, dims in sorted(add_vars.items()):
            size = tuple(data1.sizes[d] for d in dims)
            data = np.random.randint(0, 100, size=size).astype(np.str_)
            data1[v] = (dims, data, {"foo": "variable"})
        # var4 and var5 are extension arrays and should be dropped
        assert (
            "var4" not in data1.mean()
            and "var5" not in data1.mean()
            and "var6" not in data1.mean()
            and "var7" not in data1.mean()
        )
        assert_equal(data1.mean(), data2.mean())
        assert_equal(data1.mean(dim="dim1"), data2.mean(dim="dim1"))
        assert "var6" not in data1.mean(dim="dim2") and "var7" in data1.mean(dim="dim2")

    @pytest.mark.filterwarnings(
        "ignore:Once the behaviour of DataArray:DeprecationWarning"
    )
    def test_reduce_strings(self) -> None:
        expected = Dataset({"x": "a"})
        ds = Dataset({"x": ("y", ["a", "b"])})
        ds.coords["y"] = [-10, 10]
        actual = ds.min()
        assert_identical(expected, actual)

        expected = Dataset({"x": "b"})
        actual = ds.max()
        assert_identical(expected, actual)

        expected = Dataset({"x": 0})
        actual = ds.argmin()
        assert_identical(expected, actual)

        expected = Dataset({"x": 1})
        actual = ds.argmax()
        assert_identical(expected, actual)

        expected = Dataset({"x": -10})
        actual = ds.idxmin()
        assert_identical(expected, actual)

        expected = Dataset({"x": 10})
        actual = ds.idxmax()
        assert_identical(expected, actual)

        expected = Dataset({"x": b"a"})
        ds = Dataset({"x": ("y", np.array(["a", "b"], "S1"))})
        actual = ds.min()
        assert_identical(expected, actual)

        expected = Dataset({"x": "a"})
        ds = Dataset({"x": ("y", np.array(["a", "b"], "U1"))})
        actual = ds.min()
        assert_identical(expected, actual)

    def test_reduce_dtypes(self) -> None:
        # regression test for GH342
        expected = Dataset({"x": 1})
        actual = Dataset({"x": True}).sum()
        assert_identical(expected, actual)

        # regression test for GH505
        expected = Dataset({"x": 3})
        actual = Dataset({"x": ("y", np.array([1, 2], "uint16"))}).sum()
        assert_identical(expected, actual)

        expected = Dataset({"x": 1 + 1j})
        actual = Dataset({"x": ("y", [1, 1j])}).sum()
        assert_identical(expected, actual)

    def test_reduce_keep_attrs(self) -> None:
        data = create_test_data()
        _attrs = {"attr1": "value1", "attr2": 2929}

        attrs = dict(_attrs)
        data.attrs = attrs

        # Test dropped attrs
        ds = data.mean()
        assert ds.attrs == {}
        for v in ds.data_vars.values():
            assert v.attrs == {}

        # Test kept attrs
        ds = data.mean(keep_attrs=True)
        assert ds.attrs == attrs
        for k, v in ds.data_vars.items():
            assert v.attrs == data[k].attrs

    @pytest.mark.filterwarnings(
        "ignore:Once the behaviour of DataArray:DeprecationWarning"
    )
    def test_reduce_argmin(self) -> None:
        # regression test for #205
        ds = Dataset({"a": ("x", [0, 1])})
        expected = Dataset({"a": ([], 0)})
        actual = ds.argmin()
        assert_identical(expected, actual)

        actual = ds.argmin("x")
        assert_identical(expected, actual)

    def test_reduce_scalars(self) -> None:
        ds = Dataset({"x": ("a", [2, 2]), "y": 2, "z": ("b", [2])})
        expected = Dataset({"x": 0, "y": 0, "z": 0})
        actual = ds.var()
        assert_identical(expected, actual)

        expected = Dataset({"x": 0, "y": 0, "z": ("b", [0])})
        actual = ds.var("a")
        assert_identical(expected, actual)

    def test_reduce_only_one_axis(self) -> None:
        def mean_only_one_axis(x, axis):
            if not isinstance(axis, integer_types):
                raise TypeError("non-integer axis")
            return x.mean(axis)

        ds = Dataset({"a": (["x", "y"], [[0, 1, 2, 3, 4]])})
        expected = Dataset({"a": ("x", [2])})
        actual = ds.reduce(mean_only_one_axis, "y")
        assert_identical(expected, actual)

        with pytest.raises(
            TypeError, match=r"missing 1 required positional argument: 'axis'"
        ):
            ds.reduce(mean_only_one_axis)

    def test_reduce_no_axis(self) -> None:
        def total_sum(x):
            return np.sum(x.flatten())

        ds = Dataset({"a": (["x", "y"], [[0, 1, 2, 3, 4]])})
        expected = Dataset({"a": ((), 10)})
        actual = ds.reduce(total_sum)
        assert_identical(expected, actual)

        with pytest.raises(TypeError, match=r"unexpected keyword argument 'axis'"):
            ds.reduce(total_sum, dim="x")

    def test_reduce_keepdims(self) -> None:
        ds = Dataset(
            {"a": (["x", "y"], [[0, 1, 2, 3, 4]])},
            coords={
                "y": [0, 1, 2, 3, 4],
                "x": [0],
                "lat": (["x", "y"], [[0, 1, 2, 3, 4]]),
                "c": -999.0,
            },
        )

        # Shape should match behaviour of numpy reductions with keepdims=True
        # Coordinates involved in the reduction should be removed
        actual = ds.mean(keepdims=True)
        expected = Dataset(
            {"a": (["x", "y"], np.mean(ds.a, keepdims=True).data)}, coords={"c": ds.c}
        )
        assert_identical(expected, actual)

        actual = ds.mean("x", keepdims=True)
        expected = Dataset(
            {"a": (["x", "y"], np.mean(ds.a, axis=0, keepdims=True).data)},
            coords={"y": ds.y, "c": ds.c},
        )
        assert_identical(expected, actual)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize("skipna", [True, False, None])
    @pytest.mark.parametrize("q", [0.25, [0.50], [0.25, 0.75]])
    def test_quantile(self, q, skipna, compute_backend) -> None:
        ds = create_test_data(seed=123)
        ds.var1.data[0, 0] = np.nan

        for dim in [None, "dim1", ["dim1"]]:
            ds_quantile = ds.quantile(q, dim=dim, skipna=skipna)
            if is_scalar(q):
                assert "quantile" not in ds_quantile.dims
            else:
                assert "quantile" in ds_quantile.dims

            for var, dar in ds.data_vars.items():
                assert var in ds_quantile
                assert_identical(
                    ds_quantile[var], dar.quantile(q, dim=dim, skipna=skipna)
                )
        dim = ["dim1", "dim2"]
        ds_quantile = ds.quantile(q, dim=dim, skipna=skipna)
        assert "dim3" in ds_quantile.dims
        assert all(d not in ds_quantile.dims for d in dim)

    @pytest.mark.parametrize("compute_backend", ["numbagg", None], indirect=True)
    @pytest.mark.parametrize("skipna", [True, False])
    def test_quantile_skipna(self, skipna, compute_backend) -> None:
        q = 0.1
        dim = "time"
        ds = Dataset({"a": ([dim], np.arange(0, 11))})
        ds = ds.where(ds >= 1)

        result = ds.quantile(q=q, dim=dim, skipna=skipna)

        value = 1.9 if skipna else np.nan
        expected = Dataset({"a": value}, coords={"quantile": q})
        assert_identical(result, expected)

    @pytest.mark.parametrize("method", ["midpoint", "lower"])
    def test_quantile_method(self, method) -> None:
        ds = create_test_data(seed=123)
        q = [0.25, 0.5, 0.75]

        result = ds.quantile(q, method=method)

        assert_identical(result.var1, ds.var1.quantile(q, method=method))
        assert_identical(result.var2, ds.var2.quantile(q, method=method))
        assert_identical(result.var3, ds.var3.quantile(q, method=method))

    @pytest.mark.filterwarnings(
        "default:The `interpolation` argument to quantile was renamed to `method`:FutureWarning"
    )
    @pytest.mark.parametrize("method", ["midpoint", "lower"])
    def test_quantile_interpolation_deprecated(self, method) -> None:
        ds = create_test_data(seed=123)
        q = [0.25, 0.5, 0.75]

        with pytest.warns(
            FutureWarning,
            match="`interpolation` argument to quantile was renamed to `method`",
        ):
            ds.quantile(q, interpolation=method)

        with warnings.catch_warnings(record=True):
            with pytest.raises(TypeError, match="interpolation and method keywords"):
                ds.quantile(q, method=method, interpolation=method)

    @requires_bottleneck
    def test_rank(self) -> None:
        ds = create_test_data(seed=1234)
        # only ds.var3 depends on dim3
        z = ds.rank("dim3")
        assert ["var3"] == list(z.data_vars)
        # same as dataarray version
        x = z.var3
        y = ds.var3.rank("dim3")
        assert_equal(x, y)
        # coordinates stick
        assert list(z.coords) == list(ds.coords)
        assert list(x.coords) == list(y.coords)
        # invalid dim
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Dimension 'invalid_dim' not found in data dimensions ('dim3', 'dim1')"
            ),
        ):
            x.rank("invalid_dim")

    def test_rank_use_bottleneck(self) -> None:
        ds = Dataset({"a": ("x", [0, np.nan, 2]), "b": ("y", [4, 6, 3, 4])})
        with xr.set_options(use_bottleneck=False):
            with pytest.raises(RuntimeError):
                ds.rank("x")

    def test_count(self) -> None:
        ds = Dataset({"x": ("a", [np.nan, 1]), "y": 0, "z": np.nan})
        expected = Dataset({"x": 1, "y": 1, "z": 0})
        actual = ds.count()
        assert_identical(expected, actual)

    def test_map(self) -> None:
        data = create_test_data()
        data.attrs["foo"] = "bar"

        assert_identical(data.map(np.mean), data.mean())

        expected = data.mean(keep_attrs=True)
        actual = data.map(lambda x: x.mean(keep_attrs=True), keep_attrs=True)
        assert_identical(expected, actual)

        assert_identical(data.map(lambda x: x, keep_attrs=True), data.drop_vars("time"))

        def scale(x, multiple=1):
            return multiple * x

        actual = data.map(scale, multiple=2)
        assert_equal(actual["var1"], 2 * data["var1"])
        assert_identical(actual["numbers"], data["numbers"])

        actual = data.map(np.asarray)
        expected = data.drop_vars("time")  # time is not used on a data var
        assert_equal(expected, actual)

    def test_apply_pending_deprecated_map(self) -> None:
        data = create_test_data()
        data.attrs["foo"] = "bar"

        with pytest.warns(PendingDeprecationWarning):
            assert_identical(data.apply(np.mean), data.mean())

    def make_example_math_dataset(self):
        variables = {
            "bar": ("x", np.arange(100, 400, 100)),
            "foo": (("x", "y"), 1.0 * np.arange(12).reshape(3, 4)),
        }
        coords = {"abc": ("x", ["a", "b", "c"]), "y": 10 * np.arange(4)}
        ds = Dataset(variables, coords)
        ds["foo"][0, 0] = np.nan
        return ds

    def test_dataset_number_math(self) -> None:
        ds = self.make_example_math_dataset()

        assert_identical(ds, +ds)
        assert_identical(ds, ds + 0)
        assert_identical(ds, 0 + ds)
        assert_identical(ds, ds + np.array(0))
        assert_identical(ds, np.array(0) + ds)

        actual = ds.copy(deep=True)
        actual += 0
        assert_identical(ds, actual)

    # casting nan warns
    @pytest.mark.filterwarnings("ignore:invalid value encountered in cast")
    def test_unary_ops(self) -> None:
        ds = self.make_example_math_dataset()

        assert_identical(ds.map(abs), abs(ds))
        assert_identical(ds.map(lambda x: x + 4), ds + 4)

        for func in [
            lambda x: x.isnull(),
            lambda x: x.round(),
            lambda x: x.astype(int),
        ]:
            assert_identical(ds.map(func), func(ds))

        assert_identical(ds.isnull(), ~ds.notnull())

        # don't actually patch these methods in
        with pytest.raises(AttributeError):
            _ = ds.item
        with pytest.raises(AttributeError):
            _ = ds.searchsorted

    def test_dataset_array_math(self) -> None:
        ds = self.make_example_math_dataset()

        expected = ds.map(lambda x: x - ds["foo"])
        assert_identical(expected, ds - ds["foo"])
        assert_identical(expected, -ds["foo"] + ds)
        assert_identical(expected, ds - ds["foo"].variable)
        assert_identical(expected, -ds["foo"].variable + ds)
        actual = ds.copy(deep=True)
        actual -= ds["foo"]
        assert_identical(expected, actual)

        expected = ds.map(lambda x: x + ds["bar"])
        assert_identical(expected, ds + ds["bar"])
        actual = ds.copy(deep=True)
        actual += ds["bar"]
        assert_identical(expected, actual)

        expected = Dataset({"bar": ds["bar"] + np.arange(3)})
        assert_identical(expected, ds[["bar"]] + np.arange(3))
        assert_identical(expected, np.arange(3) + ds[["bar"]])

    def test_dataset_dataset_math(self) -> None:
        ds = self.make_example_math_dataset()

        assert_identical(ds, ds + 0 * ds)
        assert_identical(ds, ds + {"foo": 0, "bar": 0})

        expected = ds.map(lambda x: 2 * x)
        assert_identical(expected, 2 * ds)
        assert_identical(expected, ds + ds)
        assert_identical(expected, ds + ds.data_vars)
        assert_identical(expected, ds + dict(ds.data_vars))

        actual = ds.copy(deep=True)
        expected_id = id(actual)
        actual += ds
        assert_identical(expected, actual)
        assert expected_id == id(actual)

        assert_identical(ds == ds, ds.notnull())

        subsampled = ds.isel(y=slice(2))
        expected = 2 * subsampled
        assert_identical(expected, subsampled + ds)
        assert_identical(expected, ds + subsampled)

    def test_dataset_math_auto_align(self) -> None:
        ds = self.make_example_math_dataset()
        subset = ds.isel(y=[1, 3])
        expected = 2 * subset
        actual = ds + subset
        assert_identical(expected, actual)

        actual = ds.isel(y=slice(1)) + ds.isel(y=slice(1, None))
        expected = 2 * ds.drop_sel(y=ds.y)
        assert_equal(actual, expected)

        actual = ds + ds[["bar"]]
        expected = (2 * ds[["bar"]]).merge(ds.coords, compat="override")
        assert_identical(expected, actual)

        assert_identical(ds + Dataset(), ds.coords.to_dataset())
        assert_identical(Dataset() + Dataset(), Dataset())

        ds2 = Dataset(coords={"bar": 42})
        assert_identical(ds + ds2, ds.coords.merge(ds2))

        # maybe unary arithmetic with empty datasets should raise instead?
        assert_identical(Dataset() + 1, Dataset())

        actual = ds.copy(deep=True)
        other = ds.isel(y=slice(2))
        actual += other
        expected = ds + other.reindex_like(ds)
        assert_identical(expected, actual)

    def test_dataset_math_errors(self) -> None:
        ds = self.make_example_math_dataset()

        with pytest.raises(TypeError):
            ds["foo"] += ds
        with pytest.raises(TypeError):
            ds["foo"].variable += ds
        with pytest.raises(ValueError, match=r"must have the same"):
            ds += ds[["bar"]]

        # verify we can rollback in-place operations if something goes wrong
        # nb. inplace datetime64 math actually will work with an integer array
        # but not floats thanks to numpy's inconsistent handling
        other = DataArray(np.datetime64("2000-01-01"), coords={"c": 2})
        actual = ds.copy(deep=True)
        with pytest.raises(TypeError):
            actual += other
        assert_identical(actual, ds)

    def test_dataset_transpose(self) -> None:
        ds = Dataset(
            {
                "a": (("x", "y"), np.random.randn(3, 4)),
                "b": (("y", "x"), np.random.randn(4, 3)),
            },
            coords={
                "x": range(3),
                "y": range(4),
                "xy": (("x", "y"), np.random.randn(3, 4)),
            },
        )

        actual = ds.transpose()
        expected = Dataset(
            {"a": (("y", "x"), ds.a.values.T), "b": (("x", "y"), ds.b.values.T)},
            coords={
                "x": ds.x.values,
                "y": ds.y.values,
                "xy": (("y", "x"), ds.xy.values.T),
            },
        )
        assert_identical(expected, actual)

        actual = ds.transpose(...)
        expected = ds
        assert_identical(expected, actual)

        actual = ds.transpose("x", "y")
        expected = ds.map(lambda x: x.transpose("x", "y", transpose_coords=True))
        assert_identical(expected, actual)

        ds = create_test_data()
        actual = ds.transpose()
        for k in ds.variables:
            assert actual[k].dims[::-1] == ds[k].dims

        new_order = ("dim2", "dim3", "dim1", "time")
        actual = ds.transpose(*new_order)
        for k in ds.variables:
            expected_dims = tuple(d for d in new_order if d in ds[k].dims)
            assert actual[k].dims == expected_dims

        # same as above but with ellipsis
        new_order = ("dim2", "dim3", "dim1", "time")
        actual = ds.transpose("dim2", "dim3", ...)
        for k in ds.variables:
            expected_dims = tuple(d for d in new_order if d in ds[k].dims)
            assert actual[k].dims == expected_dims

        # test missing dimension, raise error
        with pytest.raises(ValueError):
            ds.transpose(..., "not_a_dim")

        # test missing dimension, ignore error
        actual = ds.transpose(..., "not_a_dim", missing_dims="ignore")
        expected_ell = ds.transpose(...)
        assert_identical(expected_ell, actual)

        # test missing dimension, raise warning
        with pytest.warns(UserWarning):
            actual = ds.transpose(..., "not_a_dim", missing_dims="warn")
            assert_identical(expected_ell, actual)

        assert "T" not in dir(ds)

    def test_dataset_ellipsis_transpose_different_ordered_vars(self) -> None:
        # https://github.com/pydata/xarray/issues/1081#issuecomment-544350457
        ds = Dataset(
            dict(
                a=(("w", "x", "y", "z"), np.ones((2, 3, 4, 5))),
                b=(("x", "w", "y", "z"), np.zeros((3, 2, 4, 5))),
            )
        )
        result = ds.transpose(..., "z", "y")
        assert list(result["a"].dims) == list("wxzy")
        assert list(result["b"].dims) == list("xwzy")

    def test_dataset_retains_period_index_on_transpose(self) -> None:
        ds = create_test_data()
        ds["time"] = pd.period_range("2000-01-01", periods=20)

        transposed = ds.transpose()

        assert isinstance(transposed.time.to_index(), pd.PeriodIndex)

    def test_dataset_diff_n1_simple(self) -> None:
        ds = Dataset({"foo": ("x", [5, 5, 6, 6])})
        actual = ds.diff("x")
        expected = Dataset({"foo": ("x", [0, 1, 0])})
        assert_equal(expected, actual)

    def test_dataset_diff_n1_label(self) -> None:
        ds = Dataset({"foo": ("x", [5, 5, 6, 6])}, {"x": [0, 1, 2, 3]})
        actual = ds.diff("x", label="lower")
        expected = Dataset({"foo": ("x", [0, 1, 0])}, {"x": [0, 1, 2]})
        assert_equal(expected, actual)

        actual = ds.diff("x", label="upper")
        expected = Dataset({"foo": ("x", [0, 1, 0])}, {"x": [1, 2, 3]})
        assert_equal(expected, actual)

    def test_dataset_diff_n1(self) -> None:
        ds = create_test_data(seed=1)
        actual = ds.diff("dim2")
        expected_dict = {}
        expected_dict["var1"] = DataArray(
            np.diff(ds["var1"].values, axis=1),
            {"dim2": ds["dim2"].values[1:]},
            ["dim1", "dim2"],
        )
        expected_dict["var2"] = DataArray(
            np.diff(ds["var2"].values, axis=1),
            {"dim2": ds["dim2"].values[1:]},
            ["dim1", "dim2"],
        )
        expected_dict["var3"] = ds["var3"]
        expected = Dataset(expected_dict, coords={"time": ds["time"].values})
        expected.coords["numbers"] = ("dim3", ds["numbers"].values)
        assert_equal(expected, actual)

    def test_dataset_diff_n2(self) -> None:
        ds = create_test_data(seed=1)
        actual = ds.diff("dim2", n=2)
        expected_dict = {}
        expected_dict["var1"] = DataArray(
            np.diff(ds["var1"].values, axis=1, n=2),
            {"dim2": ds["dim2"].values[2:]},
            ["dim1", "dim2"],
        )
        expected_dict["var2"] = DataArray(
            np.diff(ds["var2"].values, axis=1, n=2),
            {"dim2": ds["dim2"].values[2:]},
            ["dim1", "dim2"],
        )
        expected_dict["var3"] = ds["var3"]
        expected = Dataset(expected_dict, coords={"time": ds["time"].values})
        expected.coords["numbers"] = ("dim3", ds["numbers"].values)
        assert_equal(expected, actual)

    def test_dataset_diff_exception_n_neg(self) -> None:
        ds = create_test_data(seed=1)
        with pytest.raises(ValueError, match=r"must be non-negative"):
            ds.diff("dim2", n=-1)

    def test_dataset_diff_exception_label_str(self) -> None:
        ds = create_test_data(seed=1)
        with pytest.raises(ValueError, match=r"'label' argument has to"):
            ds.diff("dim2", label="raise_me")  # type: ignore[arg-type]

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"foo": -10}])
    def test_shift(self, fill_value) -> None:
        coords = {"bar": ("x", list("abc")), "x": [-4, 3, 2]}
        attrs = {"meta": "data"}
        ds = Dataset({"foo": ("x", [1, 2, 3])}, coords, attrs)
        actual = ds.shift(x=1, fill_value=fill_value)
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value = np.nan
        elif isinstance(fill_value, dict):
            fill_value = fill_value.get("foo", np.nan)
        expected = Dataset({"foo": ("x", [fill_value, 1, 2])}, coords, attrs)
        assert_identical(expected, actual)

        with pytest.raises(ValueError, match=r"dimensions"):
            ds.shift(foo=123)

    def test_roll_coords(self) -> None:
        coords = {"bar": ("x", list("abc")), "x": [-4, 3, 2]}
        attrs = {"meta": "data"}
        ds = Dataset({"foo": ("x", [1, 2, 3])}, coords, attrs)
        actual = ds.roll(x=1, roll_coords=True)

        ex_coords = {"bar": ("x", list("cab")), "x": [2, -4, 3]}
        expected = Dataset({"foo": ("x", [3, 1, 2])}, ex_coords, attrs)
        assert_identical(expected, actual)

        with pytest.raises(ValueError, match=r"dimensions"):
            ds.roll(foo=123, roll_coords=True)

    def test_roll_no_coords(self) -> None:
        coords = {"bar": ("x", list("abc")), "x": [-4, 3, 2]}
        attrs = {"meta": "data"}
        ds = Dataset({"foo": ("x", [1, 2, 3])}, coords, attrs)
        actual = ds.roll(x=1)

        expected = Dataset({"foo": ("x", [3, 1, 2])}, coords, attrs)
        assert_identical(expected, actual)

        with pytest.raises(ValueError, match=r"dimensions"):
            ds.roll(abc=321)

    def test_roll_multidim(self) -> None:
        # regression test for 2445
        arr = xr.DataArray(
            [[1, 2, 3], [4, 5, 6]],
            coords={"x": range(3), "y": range(2)},
            dims=("y", "x"),
        )
        actual = arr.roll(x=1, roll_coords=True)
        expected = xr.DataArray(
            [[3, 1, 2], [6, 4, 5]], coords=[("y", [0, 1]), ("x", [2, 0, 1])]
        )
        assert_identical(expected, actual)

    def test_real_and_imag(self) -> None:
        attrs = {"foo": "bar"}
        ds = Dataset({"x": ((), 1 + 2j, attrs)}, attrs=attrs)

        expected_re = Dataset({"x": ((), 1, attrs)}, attrs=attrs)
        assert_identical(ds.real, expected_re)

        expected_im = Dataset({"x": ((), 2, attrs)}, attrs=attrs)
        assert_identical(ds.imag, expected_im)

    def test_setattr_raises(self) -> None:
        ds = Dataset({}, coords={"scalar": 1}, attrs={"foo": "bar"})
        with pytest.raises(AttributeError, match=r"cannot set attr"):
            ds.scalar = 2
        with pytest.raises(AttributeError, match=r"cannot set attr"):
            ds.foo = 2
        with pytest.raises(AttributeError, match=r"cannot set attr"):
            ds.other = 2

    def test_filter_by_attrs(self) -> None:
        precip = dict(standard_name="convective_precipitation_flux")
        temp0 = dict(standard_name="air_potential_temperature", height="0 m")
        temp10 = dict(standard_name="air_potential_temperature", height="10 m")
        ds = Dataset(
            {
                "temperature_0": (["t"], [0], temp0),
                "temperature_10": (["t"], [0], temp10),
                "precipitation": (["t"], [0], precip),
            },
            coords={"time": (["t"], [0], dict(axis="T", long_name="time_in_seconds"))},
        )

        # Test return empty Dataset.
        ds.filter_by_attrs(standard_name="invalid_standard_name")
        new_ds = ds.filter_by_attrs(standard_name="invalid_standard_name")
        assert not bool(new_ds.data_vars)

        # Test return one DataArray.
        new_ds = ds.filter_by_attrs(standard_name="convective_precipitation_flux")
        assert new_ds["precipitation"].standard_name == "convective_precipitation_flux"

        assert_equal(new_ds["precipitation"], ds["precipitation"])

        # Test filter coordinates
        new_ds = ds.filter_by_attrs(long_name="time_in_seconds")
        assert new_ds["time"].long_name == "time_in_seconds"
        assert not bool(new_ds.data_vars)

        # Test return more than one DataArray.
        new_ds = ds.filter_by_attrs(standard_name="air_potential_temperature")
        assert len(new_ds.data_vars) == 2
        for var in new_ds.data_vars:
            assert new_ds[var].standard_name == "air_potential_temperature"

        # Test callable.
        new_ds = ds.filter_by_attrs(height=lambda v: v is not None)
        assert len(new_ds.data_vars) == 2
        for var in new_ds.data_vars:
            assert new_ds[var].standard_name == "air_potential_temperature"

        new_ds = ds.filter_by_attrs(height="10 m")
        assert len(new_ds.data_vars) == 1
        for var in new_ds.data_vars:
            assert new_ds[var].height == "10 m"

        # Test return empty Dataset due to conflicting filters
        new_ds = ds.filter_by_attrs(
            standard_name="convective_precipitation_flux", height="0 m"
        )
        assert not bool(new_ds.data_vars)

        # Test return one DataArray with two filter conditions
        new_ds = ds.filter_by_attrs(
            standard_name="air_potential_temperature", height="0 m"
        )
        for var in new_ds.data_vars:
            assert new_ds[var].standard_name == "air_potential_temperature"
            assert new_ds[var].height == "0 m"
            assert new_ds[var].height != "10 m"

        # Test return empty Dataset due to conflicting callables
        new_ds = ds.filter_by_attrs(
            standard_name=lambda v: False, height=lambda v: True
        )
        assert not bool(new_ds.data_vars)

    def test_binary_op_propagate_indexes(self) -> None:
        ds = Dataset(
            {"d1": DataArray([1, 2, 3], dims=["x"], coords={"x": [10, 20, 30]})}
        )
        expected = ds.xindexes["x"]
        actual = (ds * 2).xindexes["x"]
        assert expected is actual

    def test_binary_op_join_setting(self) -> None:
        # arithmetic_join applies to data array coordinates
        missing_2 = xr.Dataset({"x": [0, 1]})
        missing_0 = xr.Dataset({"x": [1, 2]})
        with xr.set_options(arithmetic_join="outer"):
            actual = missing_2 + missing_0
        expected = xr.Dataset({"x": [0, 1, 2]})
        assert_equal(actual, expected)

        # arithmetic join also applies to data_vars
        ds1 = xr.Dataset({"foo": 1, "bar": 2})
        ds2 = xr.Dataset({"bar": 2, "baz": 3})
        expected = xr.Dataset({"bar": 4})  # default is inner joining
        actual = ds1 + ds2
        assert_equal(actual, expected)

        with xr.set_options(arithmetic_join="outer"):
            expected = xr.Dataset({"foo": np.nan, "bar": 4, "baz": np.nan})
            actual = ds1 + ds2
            assert_equal(actual, expected)

        with xr.set_options(arithmetic_join="left"):
            expected = xr.Dataset({"foo": np.nan, "bar": 4})
            actual = ds1 + ds2
            assert_equal(actual, expected)

        with xr.set_options(arithmetic_join="right"):
            expected = xr.Dataset({"bar": 4, "baz": np.nan})
            actual = ds1 + ds2
            assert_equal(actual, expected)

    @pytest.mark.parametrize(
        ["keep_attrs", "expected"],
        (
            pytest.param(False, {}, id="False"),
            pytest.param(True, {"foo": "a", "bar": "b"}, id="True"),
        ),
    )
    def test_binary_ops_keep_attrs(self, keep_attrs, expected) -> None:
        ds1 = xr.Dataset({"a": 1}, attrs={"foo": "a", "bar": "b"})
        ds2 = xr.Dataset({"a": 1}, attrs={"foo": "a", "baz": "c"})
        with xr.set_options(keep_attrs=keep_attrs):
            ds_result = ds1 + ds2

        assert ds_result.attrs == expected

    def test_full_like(self) -> None:
        # For more thorough tests, see test_variable.py
        # Note: testing data_vars with mismatched dtypes
        ds = Dataset(
            {
                "d1": DataArray([1, 2, 3], dims=["x"], coords={"x": [10, 20, 30]}),
                "d2": DataArray([1.1, 2.2, 3.3], dims=["y"]),
            },
            attrs={"foo": "bar"},
        )
        actual = full_like(ds, 2)

        expected = ds.copy(deep=True)
        # https://github.com/python/mypy/issues/3004
        expected["d1"].values = [2, 2, 2]  # type: ignore[assignment,unused-ignore]
        expected["d2"].values = [2.0, 2.0, 2.0]  # type: ignore[assignment,unused-ignore]
        assert expected["d1"].dtype == int
        assert expected["d2"].dtype == float
        assert_identical(expected, actual)

        # override dtype
        actual = full_like(ds, fill_value=True, dtype=bool)
        expected = ds.copy(deep=True)
        expected["d1"].values = [True, True, True]  # type: ignore[assignment,unused-ignore]
        expected["d2"].values = [True, True, True]  # type: ignore[assignment,unused-ignore]
        assert expected["d1"].dtype == bool
        assert expected["d2"].dtype == bool
        assert_identical(expected, actual)

        # with multiple fill values
        actual = full_like(ds, {"d1": 1, "d2": 2.3})
        expected = ds.assign(d1=("x", [1, 1, 1]), d2=("y", [2.3, 2.3, 2.3]))
        assert expected["d1"].dtype == int
        assert expected["d2"].dtype == float
        assert_identical(expected, actual)

        # override multiple dtypes
        actual = full_like(ds, fill_value={"d1": 1, "d2": 2.3}, dtype={"d1": bool})
        expected = ds.assign(d1=("x", [True, True, True]), d2=("y", [2.3, 2.3, 2.3]))
        assert expected["d1"].dtype == bool
        assert expected["d2"].dtype == float
        assert_identical(expected, actual)

    def test_combine_first(self) -> None:
        dsx0 = DataArray([0, 0], [("x", ["a", "b"])]).to_dataset(name="dsx0")
        dsx1 = DataArray([1, 1], [("x", ["b", "c"])]).to_dataset(name="dsx1")

        actual = dsx0.combine_first(dsx1)
        expected = Dataset(
            {"dsx0": ("x", [0, 0, np.nan]), "dsx1": ("x", [np.nan, 1, 1])},
            coords={"x": ["a", "b", "c"]},
        )
        assert_equal(actual, expected)
        assert_equal(actual, xr.merge([dsx0, dsx1], join="outer"))

        # works just like xr.merge([self, other])
        dsy2 = DataArray([2, 2, 2], [("x", ["b", "c", "d"])]).to_dataset(name="dsy2")
        actual = dsx0.combine_first(dsy2)
        expected = xr.merge([dsy2, dsx0], join="outer")
        assert_equal(actual, expected)

    def test_sortby(self) -> None:
        ds = Dataset(
            {
                "A": DataArray(
                    [[1, 2], [3, 4], [5, 6]], [("x", ["c", "b", "a"]), ("y", [1, 0])]
                ),
                "B": DataArray([[5, 6], [7, 8], [9, 10]], dims=["x", "y"]),
            }
        )

        sorted1d = Dataset(
            {
                "A": DataArray(
                    [[5, 6], [3, 4], [1, 2]], [("x", ["a", "b", "c"]), ("y", [1, 0])]
                ),
                "B": DataArray([[9, 10], [7, 8], [5, 6]], dims=["x", "y"]),
            }
        )

        sorted2d = Dataset(
            {
                "A": DataArray(
                    [[6, 5], [4, 3], [2, 1]], [("x", ["a", "b", "c"]), ("y", [0, 1])]
                ),
                "B": DataArray([[10, 9], [8, 7], [6, 5]], dims=["x", "y"]),
            }
        )

        expected = sorted1d
        dax = DataArray([100, 99, 98], [("x", ["c", "b", "a"])])
        actual = ds.sortby(dax)
        assert_equal(actual, expected)

        # test descending order sort
        actual = ds.sortby(dax, ascending=False)
        assert_equal(actual, ds)

        # test alignment (fills in nan for 'c')
        dax_short = DataArray([98, 97], [("x", ["b", "a"])])
        actual = ds.sortby(dax_short)
        assert_equal(actual, expected)

        # test 1-D lexsort
        # dax0 is sorted first to give indices of [1, 2, 0]
        # and then dax1 would be used to move index 2 ahead of 1
        dax0 = DataArray([100, 95, 95], [("x", ["c", "b", "a"])])
        dax1 = DataArray([0, 1, 0], [("x", ["c", "b", "a"])])
        actual = ds.sortby([dax0, dax1])  # lexsort underneath gives [2, 1, 0]
        assert_equal(actual, expected)

        expected = sorted2d
        # test multi-dim sort by 1D dataarray values
        day = DataArray([90, 80], [("y", [1, 0])])
        actual = ds.sortby([day, dax])
        assert_equal(actual, expected)

        # test exception-raising
        with pytest.raises(KeyError):
            actual = ds.sortby("z")

        with pytest.raises(ValueError) as excinfo:
            actual = ds.sortby(ds["A"])
        assert "DataArray is not 1-D" in str(excinfo.value)

        expected = sorted1d
        actual = ds.sortby("x")
        assert_equal(actual, expected)

        # test pandas.MultiIndex
        indices = (("b", 1), ("b", 0), ("a", 1), ("a", 0))
        midx = pd.MultiIndex.from_tuples(indices, names=["one", "two"])
        ds_midx = Dataset(
            {
                "A": DataArray(
                    [[1, 2], [3, 4], [5, 6], [7, 8]], [("x", midx), ("y", [1, 0])]
                ),
                "B": DataArray([[5, 6], [7, 8], [9, 10], [11, 12]], dims=["x", "y"]),
            }
        )
        actual = ds_midx.sortby("x")
        midx_reversed = pd.MultiIndex.from_tuples(
            tuple(reversed(indices)), names=["one", "two"]
        )
        expected = Dataset(
            {
                "A": DataArray(
                    [[7, 8], [5, 6], [3, 4], [1, 2]],
                    [("x", midx_reversed), ("y", [1, 0])],
                ),
                "B": DataArray([[11, 12], [9, 10], [7, 8], [5, 6]], dims=["x", "y"]),
            }
        )
        assert_equal(actual, expected)

        # multi-dim sort by coordinate objects
        expected = sorted2d
        actual = ds.sortby(["x", "y"])
        assert_equal(actual, expected)

        # test descending order sort
        actual = ds.sortby(["x", "y"], ascending=False)
        assert_equal(actual, ds)

    def test_attribute_access(self) -> None:
        ds = create_test_data(seed=1)
        for key in ["var1", "var2", "var3", "time", "dim1", "dim2", "dim3", "numbers"]:
            assert_equal(ds[key], getattr(ds, key))
            assert key in dir(ds)

        for key in ["dim3", "dim1", "numbers"]:
            assert_equal(ds["var3"][key], getattr(ds.var3, key))
            assert key in dir(ds["var3"])
        # attrs
        assert ds["var3"].attrs["foo"] == ds.var3.foo
        assert "foo" in dir(ds["var3"])

    def test_ipython_key_completion(self) -> None:
        ds = create_test_data(seed=1)
        actual = ds._ipython_key_completions_()
        expected = ["var1", "var2", "var3", "time", "dim1", "dim2", "dim3", "numbers"]
        for item in actual:
            ds[item]  # should not raise
        assert sorted(actual) == sorted(expected)

        # for dataarray
        actual = ds["var3"]._ipython_key_completions_()
        expected = ["dim3", "dim1", "numbers"]
        for item in actual:
            ds["var3"][item]  # should not raise
        assert sorted(actual) == sorted(expected)

        # MultiIndex
        ds_midx = ds.stack(dim12=["dim2", "dim3"])
        actual = ds_midx._ipython_key_completions_()
        expected = [
            "var1",
            "var2",
            "var3",
            "time",
            "dim1",
            "dim2",
            "dim3",
            "numbers",
            "dim12",
        ]
        for item in actual:
            ds_midx[item]  # should not raise
        assert sorted(actual) == sorted(expected)

        # coords
        actual = ds.coords._ipython_key_completions_()
        expected = ["time", "dim1", "dim2", "dim3", "numbers"]
        for item in actual:
            ds.coords[item]  # should not raise
        assert sorted(actual) == sorted(expected)

        actual = ds["var3"].coords._ipython_key_completions_()
        expected = ["dim1", "dim3", "numbers"]
        for item in actual:
            ds["var3"].coords[item]  # should not raise
        assert sorted(actual) == sorted(expected)

        coords = Coordinates(ds.coords)
        actual = coords._ipython_key_completions_()
        expected = ["time", "dim2", "dim3", "numbers"]
        for item in actual:
            coords[item]  # should not raise
        assert sorted(actual) == sorted(expected)

        # data_vars
        actual = ds.data_vars._ipython_key_completions_()
        expected = ["var1", "var2", "var3", "dim1"]
        for item in actual:
            ds.data_vars[item]  # should not raise
        assert sorted(actual) == sorted(expected)

    def test_polyfit_output(self) -> None:
        ds = create_test_data(seed=1)

        out = ds.polyfit("dim2", 2, full=False)
        assert "var1_polyfit_coefficients" in out

        out = ds.polyfit("dim1", 2, full=True)
        assert "var1_polyfit_coefficients" in out
        assert "dim1_matrix_rank" in out

        out = ds.polyfit("time", 2)
        assert len(out.data_vars) == 0

    def test_polyfit_weighted(self) -> None:
        ds = create_test_data(seed=1)
        ds = ds.broadcast_like(ds)  # test more than 2 dimensions (issue #9972)
        ds_copy = ds.copy(deep=True)

        expected = ds.polyfit("dim2", 2)
        actual = ds.polyfit("dim2", 2, w=np.ones(ds.sizes["dim2"]))
        xr.testing.assert_identical(expected, actual)

        # Make sure weighted polyfit does not change the original object (issue #5644)
        xr.testing.assert_identical(ds, ds_copy)

    def test_polyfit_coord(self) -> None:
        # Make sure polyfit works when given a non-dimension coordinate.
        ds = create_test_data(seed=1)

        out = ds.polyfit("numbers", 2, full=False)
        assert "var3_polyfit_coefficients" in out
        assert "dim1" in out.dims
        assert "dim2" not in out
        assert "dim3" not in out

    def test_polyfit_coord_output(self) -> None:
        da = xr.DataArray(
            [1, 3, 2], dims=["x"], coords=dict(x=["a", "b", "c"], y=("x", [0, 1, 2]))
        )
        out = da.polyfit("y", deg=1)["polyfit_coefficients"]
        assert out.sel(degree=0).item() == pytest.approx(1.5)
        assert out.sel(degree=1).item() == pytest.approx(0.5)

    def test_polyfit_warnings(self) -> None:
        ds = create_test_data(seed=1)

        with warnings.catch_warnings(record=True) as ws:
            ds.var1.polyfit("dim2", 10, full=False)
            assert len(ws) == 1
            assert ws[0].category == RankWarning
            ds.var1.polyfit("dim2", 10, full=True)
            assert len(ws) == 1

    def test_polyfit_polyval(self) -> None:
        da = xr.DataArray(
            np.arange(1, 10).astype(np.float64), dims=["x"], coords=dict(x=np.arange(9))
        )

        out = da.polyfit("x", 3, full=False)
        da_fitval = xr.polyval(da.x, out.polyfit_coefficients)
        # polyval introduces very small errors (1e-16 here)
        xr.testing.assert_allclose(da_fitval, da)

        da = da.assign_coords(x=xr.date_range("2001-01-01", periods=9, freq="YS"))
        out = da.polyfit("x", 3, full=False)
        da_fitval = xr.polyval(da.x, out.polyfit_coefficients)
        xr.testing.assert_allclose(da_fitval, da, rtol=1e-3)

    @requires_cftime
    def test_polyfit_polyval_cftime(self) -> None:
        da = xr.DataArray(
            np.arange(1, 10).astype(np.float64),
            dims=["x"],
            coords=dict(
                x=xr.date_range("2001-01-01", periods=9, freq="YS", calendar="noleap")
            ),
        )
        out = da.polyfit("x", 3, full=False)
        da_fitval = xr.polyval(da.x, out.polyfit_coefficients)
        np.testing.assert_allclose(da_fitval, da)

    @staticmethod
    def _test_data_var_interior(
        original_data_var, padded_data_var, padded_dim_name, expected_pad_values
    ):
        np.testing.assert_equal(
            np.unique(padded_data_var.isel({padded_dim_name: [0, -1]})),
            expected_pad_values,
        )
        np.testing.assert_array_equal(
            padded_data_var.isel({padded_dim_name: slice(1, -1)}), original_data_var
        )

    @pytest.mark.parametrize("padded_dim_name", ["dim1", "dim2", "dim3", "time"])
    @pytest.mark.parametrize(
        ["constant_values"],
        [
            pytest.param(None, id="default"),
            pytest.param(42, id="scalar"),
            pytest.param((42, 43), id="tuple"),
            pytest.param({"dim1": 42, "dim2": 43}, id="per dim scalar"),
            pytest.param({"dim1": (42, 43), "dim2": (43, 44)}, id="per dim tuple"),
            pytest.param({"var1": 42, "var2": (42, 43)}, id="per var"),
            pytest.param({"var1": 42, "dim1": (42, 43)}, id="mixed"),
        ],
    )
    def test_pad(self, padded_dim_name, constant_values) -> None:
        ds = create_test_data(seed=1)
        padded = ds.pad({padded_dim_name: (1, 1)}, constant_values=constant_values)

        # test padded dim values and size
        for ds_dim_name, ds_dim in ds.sizes.items():
            if ds_dim_name == padded_dim_name:
                np.testing.assert_equal(padded.sizes[ds_dim_name], ds_dim + 2)
                if ds_dim_name in padded.coords:
                    assert padded[ds_dim_name][[0, -1]].isnull().all()
            else:
                np.testing.assert_equal(padded.sizes[ds_dim_name], ds_dim)

        # check if coord "numbers" with dimension dim3 is padded correctly
        if padded_dim_name == "dim3":
            assert padded["numbers"][[0, -1]].isnull().all()
            # twarning: passes but dtype changes from int to float
            np.testing.assert_array_equal(padded["numbers"][1:-1], ds["numbers"])

        # test if data_vars are paded with correct values
        for data_var_name, data_var in padded.data_vars.items():
            if padded_dim_name in data_var.dims:
                if utils.is_dict_like(constant_values):
                    if (
                        expected := constant_values.get(data_var_name, None)
                    ) is not None or (
                        expected := constant_values.get(padded_dim_name, None)
                    ) is not None:
                        self._test_data_var_interior(
                            ds[data_var_name], data_var, padded_dim_name, expected
                        )
                    else:
                        self._test_data_var_interior(
                            ds[data_var_name], data_var, padded_dim_name, 0
                        )
                elif constant_values:
                    self._test_data_var_interior(
                        ds[data_var_name], data_var, padded_dim_name, constant_values
                    )
                else:
                    self._test_data_var_interior(
                        ds[data_var_name], data_var, padded_dim_name, np.nan
                    )
            else:
                assert_array_equal(data_var, ds[data_var_name])

    @pytest.mark.parametrize(
        ["keep_attrs", "attrs", "expected"],
        [
            pytest.param(None, {"a": 1, "b": 2}, {"a": 1, "b": 2}, id="default"),
            pytest.param(False, {"a": 1, "b": 2}, {}, id="False"),
            pytest.param(True, {"a": 1, "b": 2}, {"a": 1, "b": 2}, id="True"),
        ],
    )
    def test_pad_keep_attrs(self, keep_attrs, attrs, expected) -> None:
        ds = xr.Dataset(
            {"a": ("x", [1, 2], attrs), "b": ("y", [1, 2], attrs)},
            coords={"c": ("x", [-1, 1], attrs), "d": ("y", [-1, 1], attrs)},
            attrs=attrs,
        )
        expected = xr.Dataset(
            {"a": ("x", [0, 1, 2, 0], expected), "b": ("y", [1, 2], attrs)},
            coords={
                "c": ("x", [np.nan, -1, 1, np.nan], expected),
                "d": ("y", [-1, 1], attrs),
            },
            attrs=expected,
        )

        keep_attrs_ = "default" if keep_attrs is None else keep_attrs

        with set_options(keep_attrs=keep_attrs_):
            actual = ds.pad({"x": (1, 1)}, mode="constant", constant_values=0)
            xr.testing.assert_identical(actual, expected)

        actual = ds.pad(
            {"x": (1, 1)}, mode="constant", constant_values=0, keep_attrs=keep_attrs
        )
        xr.testing.assert_identical(actual, expected)

    def test_astype_attrs(self) -> None:
        data = create_test_data(seed=123)
        data.attrs["foo"] = "bar"

        assert data.attrs == data.astype(float).attrs
        assert data.var1.attrs == data.astype(float).var1.attrs
        assert not data.astype(float, keep_attrs=False).attrs
        assert not data.astype(float, keep_attrs=False).var1.attrs

    @pytest.mark.parametrize("parser", ["pandas", "python"])
    @pytest.mark.parametrize(
        "engine", ["python", None, pytest.param("numexpr", marks=[requires_numexpr])]
    )
    @pytest.mark.parametrize(
        "backend", ["numpy", pytest.param("dask", marks=[requires_dask])]
    )
    def test_query(self, backend, engine, parser) -> None:
        """Test querying a dataset."""

        # setup test data
        np.random.seed(42)
        a = np.arange(0, 10, 1)
        b = np.random.randint(0, 100, size=10)
        c = np.linspace(0, 1, 20)
        d = np.random.choice(["foo", "bar", "baz"], size=30, replace=True).astype(
            object
        )
        e = np.arange(0, 10 * 20).reshape(10, 20)
        f = np.random.normal(0, 1, size=(10, 20, 30))
        if backend == "numpy":
            ds = Dataset(
                {
                    "a": ("x", a),
                    "b": ("x", b),
                    "c": ("y", c),
                    "d": ("z", d),
                    "e": (("x", "y"), e),
                    "f": (("x", "y", "z"), f),
                },
                coords={
                    "a2": ("x", a),
                    "b2": ("x", b),
                    "c2": ("y", c),
                    "d2": ("z", d),
                    "e2": (("x", "y"), e),
                    "f2": (("x", "y", "z"), f),
                },
            )
        elif backend == "dask":
            ds = Dataset(
                {
                    "a": ("x", da.from_array(a, chunks=3)),
                    "b": ("x", da.from_array(b, chunks=3)),
                    "c": ("y", da.from_array(c, chunks=7)),
                    "d": ("z", da.from_array(d, chunks=12)),
                    "e": (("x", "y"), da.from_array(e, chunks=(3, 7))),
                    "f": (("x", "y", "z"), da.from_array(f, chunks=(3, 7, 12))),
                },
                coords={
                    "a2": ("x", a),
                    "b2": ("x", b),
                    "c2": ("y", c),
                    "d2": ("z", d),
                    "e2": (("x", "y"), e),
                    "f2": (("x", "y", "z"), f),
                },
            )

        # query single dim, single variable
        with raise_if_dask_computes():
            actual = ds.query(x="a2 > 5", engine=engine, parser=parser)
        expect = ds.isel(x=(a > 5))
        assert_identical(expect, actual)

        # query single dim, single variable, via dict
        with raise_if_dask_computes():
            actual = ds.query(dict(x="a2 > 5"), engine=engine, parser=parser)
        expect = ds.isel(dict(x=(a > 5)))
        assert_identical(expect, actual)

        # query single dim, single variable
        with raise_if_dask_computes():
            actual = ds.query(x="b2 > 50", engine=engine, parser=parser)
        expect = ds.isel(x=(b > 50))
        assert_identical(expect, actual)

        # query single dim, single variable
        with raise_if_dask_computes():
            actual = ds.query(y="c2 < .5", engine=engine, parser=parser)
        expect = ds.isel(y=(c < 0.5))
        assert_identical(expect, actual)

        # query single dim, single string variable
        if parser == "pandas":
            # N.B., this query currently only works with the pandas parser
            # xref https://github.com/pandas-dev/pandas/issues/40436
            with raise_if_dask_computes():
                actual = ds.query(z='d2 == "bar"', engine=engine, parser=parser)
            expect = ds.isel(z=(d == "bar"))
            assert_identical(expect, actual)

        # query single dim, multiple variables
        with raise_if_dask_computes():
            actual = ds.query(x="(a2 > 5) & (b2 > 50)", engine=engine, parser=parser)
        expect = ds.isel(x=((a > 5) & (b > 50)))
        assert_identical(expect, actual)

        # query single dim, multiple variables with computation
        with raise_if_dask_computes():
            actual = ds.query(x="(a2 * b2) > 250", engine=engine, parser=parser)
        expect = ds.isel(x=(a * b) > 250)
        assert_identical(expect, actual)

        # check pandas query syntax is supported
        if parser == "pandas":
            with raise_if_dask_computes():
                actual = ds.query(
                    x="(a2 > 5) and (b2 > 50)", engine=engine, parser=parser
                )
            expect = ds.isel(x=((a > 5) & (b > 50)))
            assert_identical(expect, actual)

        # query multiple dims via kwargs
        with raise_if_dask_computes():
            actual = ds.query(x="a2 > 5", y="c2 < .5", engine=engine, parser=parser)
        expect = ds.isel(x=(a > 5), y=(c < 0.5))
        assert_identical(expect, actual)

        # query multiple dims via kwargs
        if parser == "pandas":
            with raise_if_dask_computes():
                actual = ds.query(
                    x="a2 > 5",
                    y="c2 < .5",
                    z="d2 == 'bar'",
                    engine=engine,
                    parser=parser,
                )
            expect = ds.isel(x=(a > 5), y=(c < 0.5), z=(d == "bar"))
            assert_identical(expect, actual)

        # query multiple dims via dict
        with raise_if_dask_computes():
            actual = ds.query(
                dict(x="a2 > 5", y="c2 < .5"), engine=engine, parser=parser
            )
        expect = ds.isel(dict(x=(a > 5), y=(c < 0.5)))
        assert_identical(expect, actual)

        # query multiple dims via dict
        if parser == "pandas":
            with raise_if_dask_computes():
                actual = ds.query(
                    dict(x="a2 > 5", y="c2 < .5", z="d2 == 'bar'"),
                    engine=engine,
                    parser=parser,
                )
            expect = ds.isel(dict(x=(a > 5), y=(c < 0.5), z=(d == "bar")))
            assert_identical(expect, actual)

        # test error handling
        with pytest.raises(ValueError):
            ds.query("a > 5")  # type: ignore[arg-type] # must be dict or kwargs
        with pytest.raises(ValueError):
            ds.query(x=(a > 5))
        with pytest.raises(IndexError):
            ds.query(y="a > 5")  # wrong length dimension
        with pytest.raises(IndexError):
            ds.query(x="c < .5")  # wrong length dimension
        with pytest.raises(IndexError):
            ds.query(x="e > 100")  # wrong number of dimensions
        with pytest.raises(UndefinedVariableError):
            ds.query(x="spam > 50")  # name not present


# pytest tests — new tests should go here, rather than in the class.


@pytest.mark.parametrize("parser", ["pandas", "python"])
def test_eval(ds, parser) -> None:
    """Currently much more minimal testing that `query` above, and much of the setup
    isn't used. But the risks are fairly low — `query` shares much of the code, and
    the method is currently experimental."""

    actual = ds.eval("z1 + 5", parser=parser)
    expect = ds["z1"] + 5
    assert_identical(expect, actual)

    # check pandas query syntax is supported
    if parser == "pandas":
        actual = ds.eval("(z1 > 5) and (z2 > 0)", parser=parser)
        expect = (ds["z1"] > 5) & (ds["z2"] > 0)
        assert_identical(expect, actual)


@pytest.mark.parametrize("test_elements", ([1, 2], np.array([1, 2]), DataArray([1, 2])))
def test_isin(test_elements, backend) -> None:
    expected = Dataset(
        data_vars={
            "var1": (("dim1",), [0, 1]),
            "var2": (("dim1",), [1, 1]),
            "var3": (("dim1",), [0, 1]),
        }
    ).astype("bool")

    if backend == "dask":
        expected = expected.chunk()

    result = Dataset(
        data_vars={
            "var1": (("dim1",), [0, 1]),
            "var2": (("dim1",), [1, 2]),
            "var3": (("dim1",), [0, 1]),
        }
    ).isin(test_elements)

    assert_equal(result, expected)


def test_isin_dataset() -> None:
    ds = Dataset({"x": [1, 2]})
    with pytest.raises(TypeError):
        ds.isin(ds)


@pytest.mark.parametrize(
    "unaligned_coords",
    (
        {"x": [2, 1, 0]},
        {"x": (["x"], np.asarray([2, 1, 0]))},
        {"x": (["x"], np.asarray([1, 2, 0]))},
        {"x": pd.Index([2, 1, 0])},
        {"x": Variable(dims="x", data=[0, 2, 1])},
        {"x": IndexVariable(dims="x", data=[0, 1, 2])},
        {"y": 42},
        {"y": ("x", [2, 1, 0])},
        {"y": ("x", np.asarray([2, 1, 0]))},
        {"y": (["x"], np.asarray([2, 1, 0]))},
    ),
)
@pytest.mark.parametrize("coords", ({"x": ("x", [0, 1, 2])}, {"x": [0, 1, 2]}))
def test_dataset_constructor_aligns_to_explicit_coords(
    unaligned_coords, coords
) -> None:
    a = xr.DataArray([1, 2, 3], dims=["x"], coords=unaligned_coords)

    expected = xr.Dataset(coords=coords)
    expected["a"] = a

    result = xr.Dataset({"a": a}, coords=coords)

    assert_equal(expected, result)


def test_error_message_on_set_supplied() -> None:
    with pytest.raises(TypeError, match="has invalid type <class 'set'>"):
        xr.Dataset(dict(date=[1, 2, 3], sec={4}))


@pytest.mark.parametrize("unaligned_coords", ({"y": ("b", np.asarray([2, 1, 0]))},))
def test_constructor_raises_with_invalid_coords(unaligned_coords) -> None:
    with pytest.raises(ValueError, match="not a subset of the DataArray dimensions"):
        xr.DataArray([1, 2, 3], dims=["x"], coords=unaligned_coords)


@pytest.mark.parametrize("ds", [3], indirect=True)
def test_dir_expected_attrs(ds) -> None:
    some_expected_attrs = {"pipe", "mean", "isnull", "var1", "dim2", "numbers"}
    result = dir(ds)
    assert set(result) >= some_expected_attrs


def test_dir_non_string(ds) -> None:
    # add a numbered key to ensure this doesn't break dir
    ds[5] = "foo"
    result = dir(ds)
    assert 5 not in result

    # GH2172
    sample_data = np.random.uniform(size=[2, 2000, 10000])
    x = xr.Dataset({"sample_data": (sample_data.shape, sample_data)})
    x2 = x["sample_data"]
    dir(x2)


def test_dir_unicode(ds) -> None:
    ds["unicode"] = "uni"
    result = dir(ds)
    assert "unicode" in result


def test_raise_no_warning_for_nan_in_binary_ops() -> None:
    with assert_no_warnings():
        _ = Dataset(data_vars={"x": ("y", [1, 2, np.nan])}) > 0


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("ds", (2,), indirect=True)
def test_raise_no_warning_assert_close(ds) -> None:
    assert_allclose(ds, ds)


@pytest.mark.parametrize("dask", [True, False])
@pytest.mark.parametrize("edge_order", [1, 2])
def test_differentiate(dask, edge_order) -> None:
    rs = np.random.default_rng(42)
    coord = [0.2, 0.35, 0.4, 0.6, 0.7, 0.75, 0.76, 0.8]

    da = xr.DataArray(
        rs.random((8, 6)),
        dims=["x", "y"],
        coords={"x": coord, "z": 3, "x2d": (("x", "y"), rs.random((8, 6)))},
    )
    if dask and has_dask:
        da = da.chunk({"x": 4})

    ds = xr.Dataset({"var": da})

    # along x
    actual = da.differentiate("x", edge_order)
    expected_x = xr.DataArray(
        np.gradient(da, da["x"], axis=0, edge_order=edge_order),
        dims=da.dims,
        coords=da.coords,
    )
    assert_equal(expected_x, actual)
    assert_equal(
        ds["var"].differentiate("x", edge_order=edge_order),
        ds.differentiate("x", edge_order=edge_order)["var"],
    )
    # coordinate should not change
    assert_equal(da["x"], actual["x"])

    # along y
    actual = da.differentiate("y", edge_order)
    expected_y = xr.DataArray(
        np.gradient(da, da["y"], axis=1, edge_order=edge_order),
        dims=da.dims,
        coords=da.coords,
    )
    assert_equal(expected_y, actual)
    assert_equal(actual, ds.differentiate("y", edge_order=edge_order)["var"])
    assert_equal(
        ds["var"].differentiate("y", edge_order=edge_order),
        ds.differentiate("y", edge_order=edge_order)["var"],
    )

    with pytest.raises(ValueError):
        da.differentiate("x2d")


@pytest.mark.parametrize("dask", [True, False])
def test_differentiate_datetime(dask) -> None:
    rs = np.random.default_rng(42)
    coord = np.array(
        [
            "2004-07-13",
            "2006-01-13",
            "2010-08-13",
            "2010-09-13",
            "2010-10-11",
            "2010-12-13",
            "2011-02-13",
            "2012-08-13",
        ],
        dtype="datetime64",
    )

    da = xr.DataArray(
        rs.random((8, 6)),
        dims=["x", "y"],
        coords={"x": coord, "z": 3, "x2d": (("x", "y"), rs.random((8, 6)))},
    )
    if dask and has_dask:
        da = da.chunk({"x": 4})

    # along x
    actual = da.differentiate("x", edge_order=1, datetime_unit="D")
    expected_x = xr.DataArray(
        np.gradient(
            da, da["x"].variable._to_numeric(datetime_unit="D"), axis=0, edge_order=1
        ),
        dims=da.dims,
        coords=da.coords,
    )
    assert_equal(expected_x, actual)

    actual2 = da.differentiate("x", edge_order=1, datetime_unit="h")
    assert np.allclose(actual, actual2 * 24)

    # for datetime variable
    actual = da["x"].differentiate("x", edge_order=1, datetime_unit="D")
    assert np.allclose(actual, 1.0)

    # with different date unit
    da = xr.DataArray(coord.astype("datetime64[ms]"), dims=["x"], coords={"x": coord})
    actual = da.differentiate("x", edge_order=1)
    assert np.allclose(actual, 1.0)


@requires_cftime
@pytest.mark.parametrize("dask", [True, False])
def test_differentiate_cftime(dask) -> None:
    rs = np.random.default_rng(42)
    coord = xr.date_range("2000", periods=8, freq="2ME", use_cftime=True)

    da = xr.DataArray(
        rs.random((8, 6)),
        coords={"time": coord, "z": 3, "t2d": (("time", "y"), rs.random((8, 6)))},
        dims=["time", "y"],
    )

    if dask and has_dask:
        da = da.chunk({"time": 4})

    actual = da.differentiate("time", edge_order=1, datetime_unit="D")
    expected_data = np.gradient(
        da, da["time"].variable._to_numeric(datetime_unit="D"), axis=0, edge_order=1
    )
    expected = xr.DataArray(expected_data, coords=da.coords, dims=da.dims)
    assert_equal(expected, actual)

    actual2 = da.differentiate("time", edge_order=1, datetime_unit="h")
    assert_allclose(actual, actual2 * 24)

    # Test the differentiation of datetimes themselves
    actual = da["time"].differentiate("time", edge_order=1, datetime_unit="D")
    assert_allclose(actual, xr.ones_like(da["time"]).astype(float))


@pytest.mark.parametrize("dask", [True, False])
def test_integrate(dask) -> None:
    rs = np.random.default_rng(42)
    coord = [0.2, 0.35, 0.4, 0.6, 0.7, 0.75, 0.76, 0.8]

    da = xr.DataArray(
        rs.random((8, 6)),
        dims=["x", "y"],
        coords={
            "x": coord,
            "x2": (("x",), rs.random(8)),
            "z": 3,
            "x2d": (("x", "y"), rs.random((8, 6))),
        },
    )
    if dask and has_dask:
        da = da.chunk({"x": 4})

    ds = xr.Dataset({"var": da})

    # along x
    actual = da.integrate("x")
    # coordinate that contains x should be dropped.
    expected_x = xr.DataArray(
        trapezoid(da.compute(), da["x"], axis=0),
        dims=["y"],
        coords={k: v for k, v in da.coords.items() if "x" not in v.dims},
    )
    assert_allclose(expected_x, actual.compute())
    assert_equal(ds["var"].integrate("x"), ds.integrate("x")["var"])

    # make sure result is also a dask array (if the source is dask array)
    assert isinstance(actual.data, type(da.data))

    # along y
    actual = da.integrate("y")
    expected_y = xr.DataArray(
        trapezoid(da, da["y"], axis=1),
        dims=["x"],
        coords={k: v for k, v in da.coords.items() if "y" not in v.dims},
    )
    assert_allclose(expected_y, actual.compute())
    assert_equal(actual, ds.integrate("y")["var"])
    assert_equal(ds["var"].integrate("y"), ds.integrate("y")["var"])

    # along x and y
    actual = da.integrate(("y", "x"))
    assert actual.ndim == 0

    with pytest.raises(ValueError):
        da.integrate("x2d")


@requires_scipy
@pytest.mark.parametrize("dask", [True, False])
def test_cumulative_integrate(dask) -> None:
    rs = np.random.default_rng(43)
    coord = [0.2, 0.35, 0.4, 0.6, 0.7, 0.75, 0.76, 0.8]

    da = xr.DataArray(
        rs.random((8, 6)),
        dims=["x", "y"],
        coords={
            "x": coord,
            "x2": (("x",), rs.random(8)),
            "z": 3,
            "x2d": (("x", "y"), rs.random((8, 6))),
        },
    )
    if dask and has_dask:
        da = da.chunk({"x": 4})

    ds = xr.Dataset({"var": da})

    # along x
    actual = da.cumulative_integrate("x")

    from scipy.integrate import cumulative_trapezoid

    expected_x = xr.DataArray(
        cumulative_trapezoid(da.compute(), da["x"], axis=0, initial=0.0),  # type: ignore[call-overload,unused-ignore]
        dims=["x", "y"],
        coords=da.coords,
    )
    assert_allclose(expected_x, actual.compute())
    assert_equal(
        ds["var"].cumulative_integrate("x"),
        ds.cumulative_integrate("x")["var"],
    )

    # make sure result is also a dask array (if the source is dask array)
    assert isinstance(actual.data, type(da.data))

    # along y
    actual = da.cumulative_integrate("y")
    expected_y = xr.DataArray(
        cumulative_trapezoid(da, da["y"], axis=1, initial=0.0),  # type: ignore[call-overload,unused-ignore]
        dims=["x", "y"],
        coords=da.coords,
    )
    assert_allclose(expected_y, actual.compute())
    assert_equal(actual, ds.cumulative_integrate("y")["var"])
    assert_equal(
        ds["var"].cumulative_integrate("y"),
        ds.cumulative_integrate("y")["var"],
    )

    # along x and y
    actual = da.cumulative_integrate(("y", "x"))
    assert actual.ndim == 2

    with pytest.raises(ValueError):
        da.cumulative_integrate("x2d")


@pytest.mark.parametrize("dask", [True, False])
@pytest.mark.parametrize("which_datetime", ["np", "cftime"])
def test_trapezoid_datetime(dask, which_datetime) -> None:
    rs = np.random.default_rng(42)
    coord: ArrayLike
    if which_datetime == "np":
        coord = np.array(
            [
                "2004-07-13",
                "2006-01-13",
                "2010-08-13",
                "2010-09-13",
                "2010-10-11",
                "2010-12-13",
                "2011-02-13",
                "2012-08-13",
            ],
            dtype="datetime64",
        )
    else:
        if not has_cftime:
            pytest.skip("Test requires cftime.")
        coord = xr.date_range("2000", periods=8, freq="2D", use_cftime=True)

    da = xr.DataArray(
        rs.random((8, 6)),
        coords={"time": coord, "z": 3, "t2d": (("time", "y"), rs.random((8, 6)))},
        dims=["time", "y"],
    )

    if dask and has_dask:
        da = da.chunk({"time": 4})

    actual = da.integrate("time", datetime_unit="D")
    expected_data = trapezoid(
        da.compute().data,
        duck_array_ops.datetime_to_numeric(da["time"].data, datetime_unit="D"),
        axis=0,
    )
    expected = xr.DataArray(
        expected_data,
        dims=["y"],
        coords={k: v for k, v in da.coords.items() if "time" not in v.dims},
    )
    assert_allclose(expected, actual.compute())

    # make sure result is also a dask array (if the source is dask array)
    assert isinstance(actual.data, type(da.data))

    actual2 = da.integrate("time", datetime_unit="h")
    assert_allclose(actual, actual2 / 24.0)


def test_no_dict() -> None:
    d = Dataset()
    with pytest.raises(AttributeError):
        _ = d.__dict__


def test_subclass_slots() -> None:
    """Test that Dataset subclasses must explicitly define ``__slots__``.

    .. note::
       As of 0.13.0, this is actually mitigated into a FutureWarning for any class
       defined outside of the xarray package.
    """
    with pytest.raises(AttributeError) as e:

        class MyDS(Dataset):
            pass

    assert str(e.value) == "MyDS must explicitly define __slots__"


def test_weakref() -> None:
    """Classes with __slots__ are incompatible with the weakref module unless they
    explicitly state __weakref__ among their slots
    """
    from weakref import ref

    ds = Dataset()
    r = ref(ds)
    assert r() is ds


def test_deepcopy_obj_array() -> None:
    x0 = Dataset(dict(foo=DataArray(np.array([object()]))))
    x1 = deepcopy(x0)
    assert x0["foo"].values[0] is not x1["foo"].values[0]


def test_deepcopy_recursive() -> None:
    # GH:issue:7111

    # direct recursion
    ds = xr.Dataset({"a": (["x"], [1, 2])})
    ds.attrs["other"] = ds

    # TODO: cannot use assert_identical on recursive Vars yet...
    # lets just ensure that deep copy works without RecursionError
    ds.copy(deep=True)

    # indirect recursion
    ds2 = xr.Dataset({"b": (["y"], [3, 4])})
    ds.attrs["other"] = ds2
    ds2.attrs["other"] = ds

    # TODO: cannot use assert_identical on recursive Vars yet...
    # lets just ensure that deep copy works without RecursionError
    ds.copy(deep=True)
    ds2.copy(deep=True)


def test_clip(ds) -> None:
    result = ds.clip(min=0.5)
    assert all((result.min(...) >= 0.5).values())

    result = ds.clip(max=0.5)
    assert all((result.max(...) <= 0.5).values())

    result = ds.clip(min=0.25, max=0.75)
    assert all((result.min(...) >= 0.25).values())
    assert all((result.max(...) <= 0.75).values())

    result = ds.clip(min=ds.mean("y"), max=ds.mean("y"))
    assert result.sizes == ds.sizes


class TestDropDuplicates:
    @pytest.mark.parametrize("keep", ["first", "last", False])
    def test_drop_duplicates_1d(self, keep) -> None:
        ds = xr.Dataset(
            {"a": ("time", [0, 5, 6, 7]), "b": ("time", [9, 3, 8, 2])},
            coords={"time": [0, 0, 1, 2]},
        )

        if keep == "first":
            a = [0, 6, 7]
            b = [9, 8, 2]
            time = [0, 1, 2]
        elif keep == "last":
            a = [5, 6, 7]
            b = [3, 8, 2]
            time = [0, 1, 2]
        else:
            a = [6, 7]
            b = [8, 2]
            time = [1, 2]

        expected = xr.Dataset(
            {"a": ("time", a), "b": ("time", b)}, coords={"time": time}
        )
        result = ds.drop_duplicates("time", keep=keep)
        assert_equal(expected, result)

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Dimensions ('space',) not found in data dimensions ('time',)"
            ),
        ):
            ds.drop_duplicates("space", keep=keep)


class TestNumpyCoercion:
    def test_from_numpy(self) -> None:
        ds = xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"lat": ("x", [4, 5, 6])})

        assert_identical(ds.as_numpy(), ds)

    @requires_dask
    def test_from_dask(self) -> None:
        ds = xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"lat": ("x", [4, 5, 6])})
        ds_chunked = ds.chunk(1)

        assert_identical(ds_chunked.as_numpy(), ds.compute())

    @requires_pint
    def test_from_pint(self) -> None:
        from pint import Quantity

        arr = np.array([1, 2, 3])
        ds = xr.Dataset(
            {"a": ("x", Quantity(arr, units="Pa"))},
            coords={"lat": ("x", Quantity(arr + 3, units="m"))},
        )

        expected = xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"lat": ("x", arr + 3)})
        assert_identical(ds.as_numpy(), expected)

    @requires_sparse
    def test_from_sparse(self) -> None:
        import sparse

        arr = np.diagflat([1, 2, 3])
        sparr = sparse.COO.from_numpy(arr)
        ds = xr.Dataset(
            {"a": (["x", "y"], sparr)}, coords={"elev": (("x", "y"), sparr + 3)}
        )

        expected = xr.Dataset(
            {"a": (["x", "y"], arr)}, coords={"elev": (("x", "y"), arr + 3)}
        )
        assert_identical(ds.as_numpy(), expected)

    @requires_cupy
    def test_from_cupy(self) -> None:
        import cupy as cp

        arr = np.array([1, 2, 3])
        ds = xr.Dataset(
            {"a": ("x", cp.array(arr))}, coords={"lat": ("x", cp.array(arr + 3))}
        )

        expected = xr.Dataset({"a": ("x", [1, 2, 3])}, coords={"lat": ("x", arr + 3)})
        assert_identical(ds.as_numpy(), expected)

    @requires_dask
    @requires_pint
    def test_from_pint_wrapping_dask(self) -> None:
        import dask
        from pint import Quantity

        arr = np.array([1, 2, 3])
        d = dask.array.from_array(arr)
        ds = xr.Dataset(
            {"a": ("x", Quantity(d, units="Pa"))},
            coords={"lat": ("x", Quantity(d, units="m") * 2)},
        )

        result = ds.as_numpy()
        expected = xr.Dataset({"a": ("x", arr)}, coords={"lat": ("x", arr * 2)})
        assert_identical(result, expected)


def test_string_keys_typing() -> None:
    """Tests that string keys to `variables` are permitted by mypy"""

    da = xr.DataArray(np.arange(10), dims=["x"])
    ds = xr.Dataset(dict(x=da))
    mapping = {"y": da}
    ds.assign(variables=mapping)


def test_transpose_error() -> None:
    # Transpose dataset with list as argument
    # Should raise error
    ds = xr.Dataset({"foo": (("x", "y"), [[21]]), "bar": (("x", "y"), [[12]])})

    with pytest.raises(
        TypeError,
        match=re.escape(
            "transpose requires dim to be passed as multiple arguments. Expected `'y', 'x'`. Received `['y', 'x']` instead"
        ),
    ):
        ds.transpose(["y", "x"])  # type: ignore[arg-type]
