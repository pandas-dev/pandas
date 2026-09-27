import pytest

import pandas as pd
import pandas._testing as tm

pytestmark = pytest.mark.single_cpu


def test_retain_index_attributes(temp_hdfstore, unit):
    # GH 3499, losing frequency info on index recreation
    dti = pd.date_range("2000-1-1", periods=3, freq="h", unit=unit)
    df = pd.DataFrame({"A": pd.Series(range(3), index=dti)})

    temp_hdfstore.put("data", df, format="table", track_times=False)

    result = temp_hdfstore.get("data")
    tm.assert_frame_equal(df, result)

    for attr in ["freq", "tz", "name"]:
        for idx in ["index", "columns"]:
            assert getattr(getattr(df, idx), attr, None) == getattr(
                getattr(result, idx), attr, None
            )

    dti2 = pd.date_range("2002-1-1", periods=3, freq="D", unit=unit)
    # try to append a table with a different frequency
    msg = r"the \[freq\] attribute of the existing index"
    with tm.assert_produces_warning(pd.errors.AttributeConflictWarning, match=msg):
        df2 = pd.DataFrame({"A": pd.Series(range(3), index=dti2)})
        temp_hdfstore.append("data", df2)

    assert temp_hdfstore.get_storer("data").info["index"]["freq"] is None

    # this is ok
    dti3 = pd.DatetimeIndex(
        ["2001-01-01", "2001-01-02", "2002-01-01"], dtype=f"M8[{unit}]"
    )
    df2 = pd.DataFrame(
        {
            "A": pd.Series(
                range(3),
                index=dti3,
            )
        }
    )
    temp_hdfstore.append("df2", df2)
    dti4 = pd.date_range("2002-1-1", periods=3, freq="D", unit=unit)
    df3 = pd.DataFrame({"A": pd.Series(range(3), index=dti4)})
    temp_hdfstore.append("df2", df3)


def test_retain_index_attributes2(temp_h5_path):
    msg = r"the \[freq\] attribute of the existing index"
    with tm.assert_produces_warning(pd.errors.AttributeConflictWarning, match=msg):
        df = pd.DataFrame(
            {
                "A": pd.Series(
                    range(3), index=pd.date_range("2000-1-1", periods=3, freq="h")
                )
            }
        )
        df.to_hdf(temp_h5_path, key="data", mode="w", append=True)
        df2 = pd.DataFrame(
            {
                "A": pd.Series(
                    range(3), index=pd.date_range("2002-1-1", periods=3, freq="D")
                )
            }
        )

        df2.to_hdf(temp_h5_path, key="data", append=True)

        idx = pd.date_range("2000-1-1", periods=3, freq="h")
        idx.name = "foo"
        df = pd.DataFrame({"A": pd.Series(range(3), index=idx)})
        df.to_hdf(temp_h5_path, key="data", mode="w", append=True)

    assert pd.read_hdf(temp_h5_path, key="data").index.name == "foo"

    msg = r"the \[index_name\] attribute of the existing index"
    with tm.assert_produces_warning(pd.errors.AttributeConflictWarning, match=msg):
        idx2 = pd.date_range("2001-1-1", periods=3, freq="h")
        idx2.name = "bar"
        df2 = pd.DataFrame({"A": pd.Series(range(3), index=idx2)})
        df2.to_hdf(temp_h5_path, key="data", append=True)

    assert pd.read_hdf(temp_h5_path, "data").index.name is None


def test_retain_datetime_attribute(temp_h5_path):
    ser = pd.Series(
        ["2024-08-26 15:13:14", "2024-08-26 15:14:14"],
        dtype="datetime64[us, UTC]",
    )
    dataframe = pd.DataFrame(ser)
    dataframe.to_hdf(temp_h5_path, key="Annotations", mode="w")

    recovered_dataframe = pd.read_hdf(temp_h5_path, key="Annotations")
    tm.assert_frame_equal(dataframe, recovered_dataframe)
