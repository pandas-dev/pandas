from __future__ import annotations

import contextlib
import warnings

import numpy as np
import pandas as pd
import pytest

from xarray import (
    Dataset,
    SerializationWarning,
    Variable,
    coding,
    conventions,
    date_range,
    open_dataset,
)
from xarray.backends.common import WritableCFDataStore
from xarray.backends.memory import InMemoryDataStore
from xarray.coders import CFDatetimeCoder, CFTimedeltaCoder
from xarray.conventions import decode_cf
from xarray.testing import assert_identical
from xarray.tests import (
    assert_array_equal,
    requires_cftime,
    requires_dask,
    requires_netCDF4,
)
from xarray.tests.test_backends import CFEncodedBase


class TestBoolTypeArray:
    def test_booltype_array(self) -> None:
        x = np.array([1, 0, 1, 1, 0], dtype="i1")
        bx = coding.variables.BoolTypeArray(x)
        assert bx.dtype == bool
        assert_array_equal(bx, np.array([True, False, True, True, False], dtype=bool))

        x = np.array([[1, 0, 1], [0, 1, 0]], dtype="i1")
        bx = coding.variables.BoolTypeArray(x)
        assert_array_equal(bx.transpose((1, 0)), x.transpose((1, 0)))


class TestNativeEndiannessArray:
    def test(self) -> None:
        x = np.arange(5, dtype=">i8")
        expected = np.arange(5, dtype="int64")
        a = coding.variables.NativeEndiannessArray(x)
        assert a.dtype == expected.dtype
        assert a.dtype == expected[:].dtype
        assert_array_equal(a, expected)

        y = np.arange(6, dtype=">i8").reshape((2, 3))
        b = coding.variables.NativeEndiannessArray(y)
        expected2 = np.arange(6, dtype="int64").reshape((2, 3))
        assert_array_equal(b.transpose((1, 0)), expected2.transpose((1, 0)))


def test_decode_cf_with_conflicting_fill_missing_value() -> None:
    expected = Variable(["t"], [np.nan, np.nan, 2], {"units": "foobar"})
    var = Variable(
        ["t"], np.arange(3), {"units": "foobar", "missing_value": 0, "_FillValue": 1}
    )
    with pytest.warns(SerializationWarning, match="has multiple fill"):
        actual = conventions.decode_cf_variable("t", var)
        assert_identical(actual, expected)

    expected = Variable(["t"], np.arange(10), {"units": "foobar"})

    var = Variable(
        ["t"],
        np.arange(10),
        {"units": "foobar", "missing_value": np.nan, "_FillValue": np.nan},
    )

    # the following code issues two warnings, so we need to check for both
    with pytest.warns(SerializationWarning) as winfo:
        actual = conventions.decode_cf_variable("t", var)
    for aw in winfo:
        assert "non-conforming" in str(aw.message)

    assert_identical(actual, expected)

    var = Variable(
        ["t"],
        np.arange(10),
        {
            "units": "foobar",
            "missing_value": np.float32(np.nan),
            "_FillValue": np.float32(np.nan),
        },
    )

    # the following code issues two warnings, so we need to check for both
    with pytest.warns(SerializationWarning) as winfo:
        actual = conventions.decode_cf_variable("t", var)
    for aw in winfo:
        assert "non-conforming" in str(aw.message)
    assert_identical(actual, expected)


def test_decode_cf_variable_with_mismatched_coordinates() -> None:
    # tests for decoding mismatched coordinates attributes
    # see GH #1809
    zeros1 = np.zeros((1, 5, 3))
    orig = Dataset(
        {
            "XLONG": (["x", "y"], zeros1.squeeze(0), {}),
            "XLAT": (["x", "y"], zeros1.squeeze(0), {}),
            "foo": (["time", "x", "y"], zeros1, {"coordinates": "XTIME XLONG XLAT"}),
            "time": ("time", [0.0], {"units": "hours since 2017-01-01"}),
        }
    )
    decoded = conventions.decode_cf(orig, decode_coords=True)
    assert decoded["foo"].encoding["coordinates"] == "XTIME XLONG XLAT"
    assert list(decoded.coords.keys()) == ["XLONG", "XLAT", "time"]

    decoded = conventions.decode_cf(orig, decode_coords=False)
    assert "coordinates" not in decoded["foo"].encoding
    assert decoded["foo"].attrs.get("coordinates") == "XTIME XLONG XLAT"
    assert list(decoded.coords.keys()) == ["time"]


@requires_cftime
class TestEncodeCFVariable:
    def test_incompatible_attributes(self) -> None:
        invalid_vars = [
            Variable(
                ["t"], pd.date_range("2000-01-01", periods=3), {"units": "foobar"}
            ),
            Variable(["t"], pd.to_timedelta(["1 day"]), {"units": "foobar"}),  # type: ignore[arg-type, unused-ignore]
            Variable(["t"], [0, 1, 2], {"add_offset": 0}, {"add_offset": 2}),
            Variable(["t"], [0, 1, 2], {"_FillValue": 0}, {"_FillValue": 2}),
        ]
        for var in invalid_vars:
            with pytest.raises(ValueError):
                conventions.encode_cf_variable(var)

    def test_missing_fillvalue(self) -> None:
        v = Variable(["x"], np.array([np.nan, 1, 2, 3]))
        v.encoding = {"dtype": "int16"}
        # Expect both the SerializationWarning and the RuntimeWarning from numpy
        with pytest.warns(Warning) as record:
            conventions.encode_cf_variable(v)
        # Check we got the expected warnings
        warning_messages = [str(w.message) for w in record]
        assert any(
            "floating point data as an integer" in msg for msg in warning_messages
        )
        assert any(
            "invalid value encountered in cast" in msg for msg in warning_messages
        )

    def test_multidimensional_coordinates(self) -> None:
        # regression test for GH1763
        # Set up test case with coordinates that have overlapping (but not
        # identical) dimensions.
        zeros1 = np.zeros((1, 5, 3))
        zeros2 = np.zeros((1, 6, 3))
        zeros3 = np.zeros((1, 5, 4))
        orig = Dataset(
            {
                "lon1": (["x1", "y1"], zeros1.squeeze(0), {}),
                "lon2": (["x2", "y1"], zeros2.squeeze(0), {}),
                "lon3": (["x1", "y2"], zeros3.squeeze(0), {}),
                "lat1": (["x1", "y1"], zeros1.squeeze(0), {}),
                "lat2": (["x2", "y1"], zeros2.squeeze(0), {}),
                "lat3": (["x1", "y2"], zeros3.squeeze(0), {}),
                "foo1": (["time", "x1", "y1"], zeros1, {"coordinates": "lon1 lat1"}),
                "foo2": (["time", "x2", "y1"], zeros2, {"coordinates": "lon2 lat2"}),
                "foo3": (["time", "x1", "y2"], zeros3, {"coordinates": "lon3 lat3"}),
                "time": ("time", [0.0], {"units": "hours since 2017-01-01"}),
            }
        )
        orig = conventions.decode_cf(orig)
        # Encode the coordinates, as they would be in a netCDF output file.
        enc, attrs = conventions.encode_dataset_coordinates(orig)
        # Make sure we have the right coordinates for each variable.
        foo1_coords = enc["foo1"].attrs.get("coordinates", "")
        foo2_coords = enc["foo2"].attrs.get("coordinates", "")
        foo3_coords = enc["foo3"].attrs.get("coordinates", "")
        assert foo1_coords == "lon1 lat1"
        assert foo2_coords == "lon2 lat2"
        assert foo3_coords == "lon3 lat3"
        # Should not have any global coordinates.
        assert "coordinates" not in attrs

    def test_var_with_coord_attr(self) -> None:
        # regression test for GH6310
        # don't overwrite user-defined "coordinates" attributes
        orig = Dataset(
            {"values": ("time", np.zeros(2), {"coordinates": "time lon lat"})},
            coords={
                "time": ("time", np.zeros(2)),
                "lat": ("time", np.zeros(2)),
                "lon": ("time", np.zeros(2)),
            },
        )
        # Encode the coordinates, as they would be in a netCDF output file.
        enc, attrs = conventions.encode_dataset_coordinates(orig)
        # Make sure we have the right coordinates for each variable.
        values_coords = enc["values"].attrs.get("coordinates", "")
        assert values_coords == "time lon lat"
        # Should not have any global coordinates.
        assert "coordinates" not in attrs

    def test_do_not_overwrite_user_coordinates(self) -> None:
        # don't overwrite user-defined "coordinates" encoding
        orig = Dataset(
            coords={"x": [0, 1, 2], "y": ("x", [5, 6, 7]), "z": ("x", [8, 9, 10])},
            data_vars={"a": ("x", [1, 2, 3]), "b": ("x", [3, 5, 6])},
        )
        orig["a"].encoding["coordinates"] = "y"
        orig["b"].encoding["coordinates"] = "z"
        enc, _ = conventions.encode_dataset_coordinates(orig)
        assert enc["a"].attrs["coordinates"] == "y"
        assert enc["b"].attrs["coordinates"] == "z"
        orig["a"].attrs["coordinates"] = "foo"
        with pytest.raises(ValueError, match=r"'coordinates' found in both attrs"):
            conventions.encode_dataset_coordinates(orig)

    def test_deterministic_coords_encoding(self) -> None:
        # the coordinates attribute is sorted when set by xarray.conventions ...
        # ... on a variable's coordinates attribute
        ds = Dataset({"foo": 0}, coords={"baz": 0, "bar": 0})
        vars, attrs = conventions.encode_dataset_coordinates(ds)
        assert vars["foo"].attrs["coordinates"] == "bar baz"
        assert attrs.get("coordinates") is None
        # ... on the global coordinates attribute
        ds = ds.drop_vars("foo")
        vars, attrs = conventions.encode_dataset_coordinates(ds)
        assert attrs["coordinates"] == "bar baz"

    def test_emit_coordinates_attribute_in_attrs(self) -> None:
        orig = Dataset(
            {"a": 1, "b": 1},
            coords={"t": np.array("2004-11-01T00:00:00", dtype=np.datetime64)},
        )

        orig["a"].attrs["coordinates"] = None
        enc, _ = conventions.encode_dataset_coordinates(orig)

        # check coordinate attribute emitted for 'a'
        assert "coordinates" not in enc["a"].attrs
        assert "coordinates" not in enc["a"].encoding

        # check coordinate attribute not emitted for 'b'
        assert enc["b"].attrs.get("coordinates") == "t"
        assert "coordinates" not in enc["b"].encoding

    def test_emit_coordinates_attribute_in_encoding(self) -> None:
        orig = Dataset(
            {"a": 1, "b": 1},
            coords={"t": np.array("2004-11-01T00:00:00", dtype=np.datetime64)},
        )

        orig["a"].encoding["coordinates"] = None
        enc, _ = conventions.encode_dataset_coordinates(orig)

        # check coordinate attribute emitted for 'a'
        assert "coordinates" not in enc["a"].attrs
        assert "coordinates" not in enc["a"].encoding

        # check coordinate attribute not emitted for 'b'
        assert enc["b"].attrs.get("coordinates") == "t"
        assert "coordinates" not in enc["b"].encoding


@requires_cftime
class TestDecodeCF:
    def test_dataset(self) -> None:
        original = Dataset(
            {
                "t": ("t", [0, 1, 2], {"units": "days since 2000-01-01"}),
                "foo": ("t", [0, 0, 0], {"coordinates": "y", "units": "bar"}),
                "y": ("t", [5, 10, -999], {"_FillValue": -999}),
            }
        )
        expected = Dataset(
            {"foo": ("t", [0, 0, 0], {"units": "bar"})},
            {
                "t": pd.date_range("2000-01-01", periods=3, unit="ns"),
                "y": ("t", [5.0, 10.0, np.nan]),
            },
        )
        actual = conventions.decode_cf(original)
        assert_identical(expected, actual)

    def test_invalid_coordinates(self) -> None:
        # regression test for GH308, GH1809
        original = Dataset({"foo": ("t", [1, 2], {"coordinates": "invalid"})})
        decoded = Dataset({"foo": ("t", [1, 2], {}, {"coordinates": "invalid"})})
        actual = conventions.decode_cf(original)
        assert_identical(decoded, actual)
        actual = conventions.decode_cf(original, decode_coords=False)
        assert_identical(original, actual)

    def test_decode_coordinates(self) -> None:
        # regression test for GH610
        original = Dataset(
            {"foo": ("t", [1, 2], {"coordinates": "x"}), "x": ("t", [4, 5])}
        )
        actual = conventions.decode_cf(original)
        assert actual.foo.encoding["coordinates"] == "x"

    def test_decode_coordinates_with_key_values(self) -> None:
        # regression test for GH9761
        original = Dataset(
            {
                "temp": (
                    ("y", "x"),
                    np.random.rand(2, 2),
                    {
                        "long_name": "temperature",
                        "units": "K",
                        "coordinates": "lat lon",
                        "grid_mapping": "crs",
                    },
                ),
                "x": (
                    ("x"),
                    np.arange(2),
                    {"standard_name": "projection_x_coordinate", "units": "m"},
                ),
                "y": (
                    ("y"),
                    np.arange(2),
                    {"standard_name": "projection_y_coordinate", "units": "m"},
                ),
                "lat": (
                    ("y", "x"),
                    np.random.rand(2, 2),
                    {"standard_name": "latitude", "units": "degrees_north"},
                ),
                "lon": (
                    ("y", "x"),
                    np.random.rand(2, 2),
                    {"standard_name": "longitude", "units": "degrees_east"},
                ),
                "crs": (
                    (),
                    None,
                    {
                        "grid_mapping_name": "transverse_mercator",
                        "longitude_of_central_meridian": -2.0,
                    },
                ),
                "crs2": (
                    (),
                    None,
                    {
                        "grid_mapping_name": "longitude_latitude",
                        "longitude_of_central_meridian": -2.0,
                    },
                ),
            },
        )

        original.temp.attrs["grid_mapping"] = "crs: x y"
        _vars, _attrs, coords = conventions.decode_cf_variables(
            original.variables, {}, decode_coords="all"
        )
        assert coords == {"lat", "lon", "crs"}

        original.temp.attrs["grid_mapping"] = "crs: x y crs2: lat lon"
        _vars, _attrs, coords = conventions.decode_cf_variables(
            original.variables, {}, decode_coords="all"
        )
        assert coords == {"lat", "lon", "crs", "crs2"}

        # stray colon
        original.temp.attrs["grid_mapping"] = "crs: x y crs2 : lat lon"
        _vars, _attrs, coords = conventions.decode_cf_variables(
            original.variables, {}, decode_coords="all"
        )
        assert coords == {"lat", "lon", "crs", "crs2"}

        original.temp.attrs["grid_mapping"] = "crs x y crs2: lat lon"
        with pytest.raises(ValueError, match="misses ':'"):
            conventions.decode_cf_variables(original.variables, {}, decode_coords="all")

        del original.temp.attrs["grid_mapping"]
        original.temp.attrs["formula_terms"] = "A: lat D: lon E: crs2"
        _vars, _attrs, coords = conventions.decode_cf_variables(
            original.variables, {}, decode_coords="all"
        )
        assert coords == {"lat", "lon", "crs2"}

        original.temp.attrs["formula_terms"] = "A: lat lon D: crs E: crs2"
        with pytest.warns(UserWarning, match="has malformed content"):
            _vars, _attrs, coords = conventions.decode_cf_variables(
                original.variables, {}, decode_coords="all"
            )
            assert coords == {"lat", "lon", "crs", "crs2"}

    def test_0d_int32_encoding(self) -> None:
        original = Variable((), np.int32(0), encoding={"dtype": "int64"})
        expected = Variable((), np.int64(0))
        actual = coding.variables.NonStringCoder().encode(original)
        assert_identical(expected, actual)

    def test_decode_cf_with_multiple_missing_values(self) -> None:
        original = Variable(["t"], [0, 1, 2], {"missing_value": np.array([0, 1])})
        expected = Variable(["t"], [np.nan, np.nan, 2], {})
        with pytest.warns(SerializationWarning, match="has multiple fill"):
            actual = conventions.decode_cf_variable("t", original)
            assert_identical(expected, actual)

    def test_decode_cf_with_drop_variables(self) -> None:
        original = Dataset(
            {
                "t": ("t", [0, 1, 2], {"units": "days since 2000-01-01"}),
                "x": ("x", [9, 8, 7], {"units": "km"}),
                "foo": (
                    ("t", "x"),
                    [[0, 0, 0], [1, 1, 1], [2, 2, 2]],
                    {"units": "bar"},
                ),
                "y": ("t", [5, 10, -999], {"_FillValue": -999}),
            }
        )
        expected = Dataset(
            {
                "t": pd.date_range("2000-01-01", periods=3, unit="ns"),
                "foo": (
                    ("t", "x"),
                    [[0, 0, 0], [1, 1, 1], [2, 2, 2]],
                    {"units": "bar"},
                ),
                "y": ("t", [5, 10, np.nan]),
            }
        )
        actual = conventions.decode_cf(original, drop_variables=("x",))
        actual2 = conventions.decode_cf(original, drop_variables="x")
        assert_identical(expected, actual)
        assert_identical(expected, actual2)

    @pytest.mark.filterwarnings("ignore:Ambiguous reference date string")
    def test_invalid_time_units_raises_eagerly(self) -> None:
        ds = Dataset({"time": ("time", [0, 1], {"units": "foobar since 123"})})
        with pytest.raises(ValueError, match=r"unable to decode time"):
            decode_cf(ds)

    @pytest.mark.parametrize("decode_times", [True, False])
    def test_invalid_timedelta_units_do_not_decode(self, decode_times) -> None:
        # regression test for #8269
        ds = Dataset(
            {"time": ("time", [0, 1, 20], {"units": "days invalid", "_FillValue": 20})}
        )
        expected = Dataset(
            {"time": ("time", [0.0, 1.0, np.nan], {"units": "days invalid"})}
        )
        assert_identical(expected, decode_cf(ds, decode_times=decode_times))

    @requires_cftime
    @pytest.mark.parametrize("time_unit", ["s", "ms", "us", "ns"])
    def test_dataset_repr_with_netcdf4_datetimes(self, time_unit) -> None:
        # regression test for #347
        attrs = {"units": "days since 0001-01-01", "calendar": "noleap"}
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "unable to decode time")
            ds = decode_cf(Dataset({"time": ("time", [0, 1], attrs)}))
            assert "(time) object" in repr(ds)

        attrs = {"units": "days since 1900-01-01"}
        ds = decode_cf(
            Dataset({"time": ("time", [0, 1], attrs)}),
            decode_times=CFDatetimeCoder(time_unit=time_unit),
        )
        assert f"(time) datetime64[{time_unit}]" in repr(ds)

    @requires_cftime
    def test_decode_cf_datetime_transition_to_invalid(self) -> None:
        # manually create dataset with not-decoded date
        from datetime import datetime

        ds = Dataset(coords={"time": [0, 266 * 365]})
        units = "days since 2000-01-01 00:00:00"
        ds.time.attrs = dict(units=units)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "unable to decode time")
            ds_decoded = conventions.decode_cf(ds)

        expected = np.array([datetime(2000, 1, 1, 0, 0), datetime(2265, 10, 28, 0, 0)])

        assert_array_equal(ds_decoded.time.values, expected)

    @requires_dask
    def test_decode_cf_with_dask(self) -> None:
        import dask.array as da

        original = Dataset(
            {
                "t": ("t", [0, 1, 2], {"units": "days since 2000-01-01"}),
                "foo": ("t", [0, 0, 0], {"coordinates": "y", "units": "bar"}),
                "bar": ("string2", [b"a", b"b"]),
                "baz": (("x"), [b"abc"], {"_Encoding": "utf-8"}),
                "y": ("t", [5, 10, -999], {"_FillValue": -999}),
            }
        ).chunk()
        decoded = conventions.decode_cf(original)
        assert all(
            isinstance(var.data, da.Array)
            for name, var in decoded.variables.items()
            if name not in decoded.xindexes
        )
        assert_identical(decoded, conventions.decode_cf(original).compute())

    @requires_dask
    def test_decode_dask_times(self) -> None:
        original = Dataset.from_dict(
            {
                "coords": {},
                "dims": {"time": 5},
                "data_vars": {
                    "average_T1": {
                        "dims": ("time",),
                        "attrs": {"units": "days since 1958-01-01 00:00:00"},
                        "data": [87659.0, 88024.0, 88389.0, 88754.0, 89119.0],
                    }
                },
            }
        )
        assert_identical(
            conventions.decode_cf(original.chunk()),
            conventions.decode_cf(original).chunk(),
        )

    @pytest.mark.parametrize("time_unit", ["s", "ms", "us", "ns"])
    def test_decode_cf_time_kwargs(self, time_unit) -> None:
        ds = Dataset.from_dict(
            {
                "coords": {
                    "timedelta": {
                        "data": np.array([1, 2, 3], dtype="int64"),
                        "dims": "timedelta",
                        "attrs": {"units": "days"},
                    },
                    "time": {
                        "data": np.array([1, 2, 3], dtype="int64"),
                        "dims": "time",
                        "attrs": {"units": "days since 2000-01-01"},
                    },
                },
                "dims": {"time": 3, "timedelta": 3},
                "data_vars": {
                    "a": {"dims": ("time", "timedelta"), "data": np.ones((3, 3))},
                },
            }
        )

        dsc = conventions.decode_cf(
            ds,
            decode_times=CFDatetimeCoder(time_unit=time_unit),
            decode_timedelta=CFTimedeltaCoder(time_unit=time_unit),
        )
        assert dsc.timedelta.dtype == np.dtype(f"m8[{time_unit}]")
        assert dsc.time.dtype == np.dtype(f"M8[{time_unit}]")
        dsc = conventions.decode_cf(ds, decode_times=False)
        assert dsc.timedelta.dtype == np.dtype("int64")
        assert dsc.time.dtype == np.dtype("int64")
        dsc = conventions.decode_cf(
            ds,
            decode_times=CFDatetimeCoder(time_unit=time_unit),
            decode_timedelta=False,
        )
        assert dsc.timedelta.dtype == np.dtype("int64")
        assert dsc.time.dtype == np.dtype(f"M8[{time_unit}]")
        dsc = conventions.decode_cf(ds, decode_times=False, decode_timedelta=True)
        assert dsc.timedelta.dtype == np.dtype("m8[ns]")
        assert dsc.time.dtype == np.dtype("int64")


class CFEncodedInMemoryStore(WritableCFDataStore, InMemoryDataStore):
    def encode_variable(self, var, name=None):
        """encode one variable"""
        coder = coding.strings.EncodedStringCoder(allows_unicode=True)
        var = coder.encode(var, name=name)
        return var


@requires_netCDF4
class TestCFEncodedDataStore(CFEncodedBase):
    @contextlib.contextmanager
    def create_store(self):
        yield CFEncodedInMemoryStore()

    @contextlib.contextmanager
    def roundtrip(
        self, data, save_kwargs=None, open_kwargs=None, allow_cleanup_failure=False
    ):
        if save_kwargs is None:
            save_kwargs = {}
        if open_kwargs is None:
            open_kwargs = {}
        store = CFEncodedInMemoryStore()
        data.dump_to_store(store, **save_kwargs)
        yield open_dataset(store, **open_kwargs)

    @pytest.mark.skip("cannot roundtrip coordinates yet for CFEncodedInMemoryStore")
    def test_roundtrip_coordinates(self) -> None:
        pass

    def test_invalid_dataarray_names_raise(self) -> None:
        # only relevant for on-disk file formats
        pass

    def test_encoding_kwarg(self) -> None:
        # we haven't bothered to raise errors yet for unexpected encodings in
        # this test dummy
        pass

    def test_encoding_kwarg_fixed_width_string(self) -> None:
        # CFEncodedInMemoryStore doesn't support explicit string encodings.
        pass

    def test_encoding_unlimited_dims(self) -> None:
        # CFEncodedInMemoryStore doesn't support unlimited_dims.
        pass


class TestDecodeCFVariableWithArrayUnits:
    def test_decode_cf_variable_with_array_units(self) -> None:
        v = Variable(["t"], [1, 2, 3], {"units": np.array(["foobar"], dtype=object)})
        v_decoded = conventions.decode_cf_variable("test2", v)
        assert_identical(v, v_decoded)


def test_decode_cf_variable_timedelta64():
    variable = Variable(["time"], pd.timedelta_range("1D", periods=2))
    decoded = conventions.decode_cf_variable("time", variable)
    assert decoded.encoding == {}
    assert_identical(decoded, variable)


def test_decode_cf_variable_datetime64():
    variable = Variable(["time"], pd.date_range("2000", periods=2))
    decoded = conventions.decode_cf_variable("time", variable)
    assert decoded.encoding == {}
    assert_identical(decoded, variable)


@requires_cftime
def test_decode_cf_variable_cftime():
    variable = Variable(["time"], date_range("2000", periods=2, use_cftime=True))
    decoded = conventions.decode_cf_variable("time", variable)
    assert decoded.encoding == {}
    assert_identical(decoded, variable)


def test_scalar_units() -> None:
    # test that scalar units does not raise an exception
    var = Variable(["t"], [np.nan, np.nan, 2], {"units": np.nan})

    actual = conventions.decode_cf_variable("t", var)
    assert_identical(actual, var)


def test_decode_cf_error_includes_variable_name():
    ds = Dataset({"my_invalid_var": ([], 1e36, {"units": "days since 2000-01-01"})})
    with pytest.raises(
        ValueError,
        match=r"unable to decode(?s:.*)my_invalid_var",
    ):
        decode_cf(ds)


def test_encode_cf_variable_with_vlen_dtype() -> None:
    v = Variable(
        ["x"], np.array(["a", "b"], dtype=coding.strings.create_vlen_dtype(str))
    )
    encoded_v = conventions.encode_cf_variable(v)
    assert encoded_v.data.dtype.kind == "O"
    assert coding.strings.check_vlen_dtype(encoded_v.data.dtype) is str

    # empty array
    v = Variable(["x"], np.array([], dtype=coding.strings.create_vlen_dtype(str)))
    encoded_v = conventions.encode_cf_variable(v)
    assert encoded_v.data.dtype.kind == "O"
    assert coding.strings.check_vlen_dtype(encoded_v.data.dtype) is str


def test_decode_cf_variables_decode_timedelta_warning() -> None:
    v = Variable(["time"], [1, 2], attrs={"units": "seconds"})
    variables = {"a": v}

    with warnings.catch_warnings():
        warnings.filterwarnings("error", "decode_timedelta", FutureWarning)
        conventions.decode_cf_variables(variables, {}, decode_timedelta=True)

    with pytest.warns(FutureWarning, match="decode_timedelta"):
        conventions.decode_cf_variables(variables, {})
