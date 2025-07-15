from __future__ import annotations

import numpy as np
import pytest

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame(
        {
            "x": [None, 0, 1, 2, 3, 4] * 2,
            "ts": [
                pd.Timestamp("2017-05-09 00:00:00.006000"),
                pd.Timestamp("2017-05-09 00:00:00.006000"),
                pd.Timestamp("2017-05-09 07:56:23.858694"),
                pd.Timestamp("2017-05-09 05:59:58.938999"),
                None,
                None,
            ]
            * 2,
            "td": [
                np.timedelta64(3, "D"),
                np.timedelta64(1, "D"),
                None,
                None,
                np.timedelta64(3, "D"),
                np.timedelta64(1, "D"),
            ]
            * 2,
            "y": "a",
        }
    )
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=2)


def _drop_mean(df, col=None):
    """TODO: In pandas 2.0, mean is implemented for datetimes, but Dask returns None."""
    if isinstance(df, pd.DataFrame):
        df.at["mean", col] = np.nan
        df.dropna(how="all", inplace=True)
    elif isinstance(df, pd.Series):
        df.drop(labels=["mean"], inplace=True, errors="ignore")
    else:
        raise NotImplementedError("Expected Series or DataFrame with mean")
    return df


def test_describe_series(df, pdf):
    assert_eq(df.x.describe(), pdf.x.describe())
    assert_eq(df.y.describe(), pdf.y.describe())
    assert_eq(df.ts.describe(), _drop_mean(pdf.ts.describe()))
    assert_eq(df.td.describe(), pdf.td.describe())


def test_describe_df(df, pdf):
    assert_eq(df.describe(), _drop_mean(pdf.describe(), "ts"))
