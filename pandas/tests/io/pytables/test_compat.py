from pathlib import Path

import pytest

from pandas.compat.numpy import np_version_gt2

import pandas as pd
import pandas._testing as tm
from pandas.tests.io.generate_legacy_storage_files import create_dataframe_all_types
from pandas.util.version import Version

tables = pytest.importorskip("tables")


@pytest.fixture
def pytables_hdf5_file(temp_h5_path):
    """
    Use PyTables to create a simple HDF5 file.
    """
    table_schema = {
        "c0": tables.Time64Col(pos=0),
        "c1": tables.StringCol(5, pos=1),
        "c2": tables.Int64Col(pos=2),
    }

    t0 = 1_561_105_000.0

    testsamples = [
        {"c0": t0, "c1": "aaaaa", "c2": 1},
        {"c0": t0 + 1, "c1": "bbbbb", "c2": 2},
        {"c0": t0 + 2, "c1": "ccccc", "c2": 10**5},
        {"c0": t0 + 3, "c1": "ddddd", "c2": 4_294_967_295},
    ]

    objname = "pandas_test_timeseries"

    with tables.open_file(temp_h5_path, mode="w") as f:
        t = f.create_table("/", name=objname, description=table_schema)
        for sample in testsamples:
            for key, value in sample.items():
                t.row[key] = value
            t.row.append()

    return temp_h5_path, objname, pd.DataFrame(testsamples)


class TestReadPyTablesHDF5:
    """
    A group of tests which covers reading HDF5 files written by plain PyTables
    (not written by pandas).

    Was introduced for regression-testing issue 11188.
    """

    def test_read_complete(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        result = pd.read_hdf(path, key=objname)
        expected = df
        tm.assert_frame_equal(result, expected, check_index_type=True)

    def test_read_with_start(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        # This is a regression test for pandas-dev/pandas/issues/11188
        result = pd.read_hdf(path, key=objname, start=1)
        expected = df[1:].reset_index(drop=True)
        tm.assert_frame_equal(result, expected, check_index_type=True)

    def test_read_with_stop(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        # This is a regression test for pandas-dev/pandas/issues/11188
        result = pd.read_hdf(path, key=objname, stop=1)
        expected = df[:1].reset_index(drop=True)
        tm.assert_frame_equal(result, expected, check_index_type=True)

    def test_read_with_startstop(self, pytables_hdf5_file):
        path, objname, df = pytables_hdf5_file
        # This is a regression test for pandas-dev/pandas/issues/11188
        result = pd.read_hdf(path, key=objname, start=1, stop=2)
        expected = df[1:2].reset_index(drop=True)
        tm.assert_frame_equal(result, expected, check_index_type=True)


_legacy_files = list(Path(__file__).parent.parent.glob("data/legacy_hdf/*/*.h5"))


@pytest.mark.parametrize("legacy_file", _legacy_files, ids=lambda x: x.name)
def test_legacy_files(datapath, legacy_file, using_infer_string, request):
    legacy_version = Version(legacy_file.parent.name)
    legacy_file = datapath(legacy_file)

    if not np_version_gt2 and legacy_file.endswith("fixed.h5"):
        # Files created for versions 2.0-3.0 used a numpy version >= 2.0, and
        # unpickling the object dtype column fails with older numpy
        pytest.skip("Fixed format pickle objects don't deserialize with numpy < 2.0")

    result = pd.read_hdf(legacy_file)

    expected = create_dataframe_all_types()

    # the fixed format doesn't include categorical columns (not supported)
    if legacy_file.endswith("fixed.h5"):
        expected = expected.drop(
            # columns=["categorical", "categorical_object", "categorical_int"]
            columns=["categorical_int"]
        )

    # # object dtype columns with strings get read as `str`
    # if using_infer_string:
    #     expected["object"] = expected["object"].astype("str")
    #     expected["object_nan"] = expected["object_nan"].astype("str")
    #     if legacy_file.endswith("table.h5"):
    #         expected["categorical_object"] = expected["categorical_object"].astype(
    #             pd.CategoricalDtype(
    #                 expected["categorical_object"].cat.categories.astype("str")
    #             )
    #         )
    # else:
    #     expected["string"] = expected["string"].astype("object")
    #     if legacy_file.endswith("table.h5"):
    #         expected["object"] = expected["object"].fillna(np.nan)
    #         expected["categorical"] = expected["categorical"].astype(
    #             pd.CategoricalDtype(
    #                 expected["categorical"].cat.categories.astype(object)
    #             )
    #         )
    #     else:
    #         expected["string"] = expected["string"].fillna("nan")

    if legacy_version < Version("2.2.0") or (
        legacy_version < Version("3.0.0") and legacy_file.endswith("fixed.h5")
    ):
        # timedelta columns gets read as nanoseconds, resulting in buggy values
        # (this also happened for direct roundtrips with those versions)
        assert not result["timedelta_us"].equals(expected["timedelta_us"])
        assert not result["timedelta_ms"].equals(expected["timedelta_ms"])
        assert not result["timedelta_s"].equals(expected["timedelta_s"])
        result = result.drop(columns=["timedelta_us", "timedelta_ms", "timedelta_s"])
        expected = expected.drop(
            columns=["timedelta_us", "timedelta_ms", "timedelta_s"]
        )

    if legacy_version < Version("2.2.0"):
        # datetime columns gets read as nanoseconds, resulting in buggy values
        # (this also happened for direct roundtrips with those versions)
        assert not result["datetime_us"].equals(expected["datetime_us"])
        assert not result["datetime_ms"].equals(expected["datetime_ms"])
        assert not result["datetime_s"].equals(expected["datetime_s"])
        assert not result["datetimetz_us"].equals(expected["datetimetz_us"])
        result = result.drop(
            columns=["datetime_us", "datetime_ms", "datetime_s", "datetimetz_us"]
        )
        expected = expected.drop(
            columns=["datetime_us", "datetime_ms", "datetime_s", "datetimetz_us"]
        )

    tm.assert_frame_equal(result, expected)
