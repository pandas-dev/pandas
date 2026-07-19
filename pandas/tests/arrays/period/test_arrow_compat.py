import pytest

from pandas.core.dtypes.dtypes import PeriodDtype

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import (
    PeriodArray,
)

pytestmark = pytest.mark.filterwarnings(
    "ignore:Passing a BlockManager to DataFrame:DeprecationWarning"
)


pa = pytest.importorskip("pyarrow")


def test_arrow_extension_type():
    from pandas.core.arrays.arrow.extension_types import ArrowPeriodType

    p1 = ArrowPeriodType("D")
    p2 = ArrowPeriodType("D")
    p3 = ArrowPeriodType("M")

    assert p1.freq == "D"
    assert p1 == p2
    assert p1 != p3
    assert hash(p1) == hash(p2)
    assert hash(p1) != hash(p3)


@pytest.mark.single_cpu
def test_arrow_extension_type_registered_on_import(tmp_path):
    # GH#41432 the extension types are registered when pandas is imported, so
    # pyarrow reads them back correctly even without a prior pandas round-trip in
    # the reading process. Uses a subprocess for a clean import state, since
    # registration is a process-wide side effect other tests already trigger.
    import subprocess
    import sys

    path = tmp_path / "period.parquet"
    pd.DataFrame({"a": pd.period_range("2012", freq="Y", periods=3)}).to_parquet(path)

    code = (
        "import pandas  # noqa: F401  # registers the extension types on import\n"
        "import pyarrow.parquet as pq\n"
        f"typ = pq.read_table({str(path)!r}).schema.field('a').type\n"
        "assert str(typ).startswith('extension<pandas.period'), typ\n"
    )
    subprocess.check_call([sys.executable, "-c", code])


@pytest.mark.parametrize(
    "data, freq",
    [
        (pd.date_range("2017", periods=3), "D"),
        (pd.date_range("2017", periods=3, freq="YE"), "Y-DEC"),
    ],
)
def test_arrow_array(data, freq):
    from pandas.core.arrays.arrow.extension_types import ArrowPeriodType

    periods = pd.PeriodIndex(data, freq=freq).array
    result = pa.array(periods)
    assert isinstance(result.type, ArrowPeriodType)
    assert result.type.freq == freq
    expected = pa.array(periods.asi8, type="int64")
    assert result.storage.equals(expected)

    # convert to its storage type
    result = pa.array(periods, type=pa.int64())
    assert result.equals(expected)

    # unsupported conversions
    msg = "Not supported to convert PeriodArray to 'double' type"
    with pytest.raises(TypeError, match=msg):
        pa.array(periods, type="float64")


def test_arrow_array_missing():
    from pandas.core.arrays.arrow.extension_types import ArrowPeriodType

    arr = PeriodArray([1, 2, 3], dtype="period[D]")
    arr[1] = pd.NaT

    result = pa.array(arr)
    assert isinstance(result.type, ArrowPeriodType)
    assert result.type.freq == "D"
    expected = pa.array([1, None, 3], type="int64")
    assert result.storage.equals(expected)


def test_arrow_table_roundtrip():
    from pandas.core.arrays.arrow.extension_types import ArrowPeriodType

    arr = PeriodArray([1, 2, 3], dtype="period[D]")
    arr[1] = pd.NaT
    df = pd.DataFrame({"a": arr})

    table = pa.table(df)
    assert isinstance(table.field("a").type, ArrowPeriodType)
    result = table.to_pandas()
    assert isinstance(result["a"].dtype, PeriodDtype)
    tm.assert_frame_equal(result, df)

    table2 = pa.concat_tables([table, table])
    result = table2.to_pandas()
    expected = pd.concat([df, df], ignore_index=True)
    tm.assert_frame_equal(result, expected)


def test_arrow_load_from_zero_chunks():
    # GH-41040

    from pandas.core.arrays.arrow.extension_types import ArrowPeriodType

    arr = PeriodArray([], dtype="period[D]")
    df = pd.DataFrame({"a": arr})

    table = pa.table(df)
    assert isinstance(table.field("a").type, ArrowPeriodType)
    table = pa.table(
        [pa.chunked_array([], type=table.column(0).type)], schema=table.schema
    )

    result = table.to_pandas()
    assert isinstance(result["a"].dtype, PeriodDtype)
    tm.assert_frame_equal(result, df)


def test_arrow_table_roundtrip_without_metadata():
    arr = PeriodArray([1, 2, 3], dtype="period[h]")
    arr[1] = pd.NaT
    df = pd.DataFrame({"a": arr})

    table = pa.table(df)
    # remove the metadata
    table = table.replace_schema_metadata()
    assert table.schema.metadata is None

    result = table.to_pandas()
    assert isinstance(result["a"].dtype, PeriodDtype)
    tm.assert_frame_equal(result, df)
