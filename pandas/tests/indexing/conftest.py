import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
    Series,
    date_range,
)
from pandas.core.api import (
    Float64Index,
    UInt64Index,
)


@pytest.fixture(name="series_ints")
def fixture_series_ints():
    return Series(np.random.rand(4), index=np.arange(0, 8, 2))


@pytest.fixture(name="frame_ints")
def fixture_frame_ints():
    return DataFrame(
        np.random.randn(4, 4), index=np.arange(0, 8, 2), columns=np.arange(0, 12, 3)
    )


@pytest.fixture(name="series_uints")
def fixture_series_uints():
    return Series(np.random.rand(4), index=UInt64Index(np.arange(0, 8, 2)))


@pytest.fixture(name="frame_uints")
def fixture_frame_uints():
    return DataFrame(
        np.random.randn(4, 4),
        index=UInt64Index(range(0, 8, 2)),
        columns=UInt64Index(range(0, 12, 3)),
    )


@pytest.fixture(name="series_labels")
def fixture_series_labels():
    return Series(np.random.randn(4), index=list("abcd"))


@pytest.fixture(name="frame_labels")
def fixture_frame_labels():
    return DataFrame(np.random.randn(4, 4), index=list("abcd"), columns=list("ABCD"))


@pytest.fixture(name="series_ts")
def fixture_series_ts():
    return Series(np.random.randn(4), index=date_range("20130101", periods=4))


@pytest.fixture(name="frame_ts")
def fixture_frame_ts():
    return DataFrame(np.random.randn(4, 4), index=date_range("20130101", periods=4))


@pytest.fixture(name="series_floats")
def fixture_series_floats():
    return Series(np.random.rand(4), index=Float64Index(range(0, 8, 2)))


@pytest.fixture(name="frame_floats")
def fixture_frame_floats():
    return DataFrame(
        np.random.randn(4, 4),
        index=Float64Index(range(0, 8, 2)),
        columns=Float64Index(range(0, 12, 3)),
    )


@pytest.fixture(name="series_mixed")
def fixture_series_mixed():
    return Series(np.random.randn(4), index=[2, 4, "null", 8])


@pytest.fixture(name="frame_mixed")
def fixture_frame_mixed():
    return DataFrame(np.random.randn(4, 4), index=[2, 4, "null", 8])


@pytest.fixture(name="frame_empty")
def fixture_frame_empty():
    return DataFrame()


@pytest.fixture(name="series_empty")
def fixture_series_empty():
    return Series(dtype=object)


@pytest.fixture(name="frame_multi")
def fixture_frame_multi():
    return DataFrame(
        np.random.randn(4, 4),
        index=MultiIndex.from_product([[1, 2], [3, 4]]),
        columns=MultiIndex.from_product([[5, 6], [7, 8]]),
    )


@pytest.fixture(name="series_multi")
def fixture_series_multi():
    return Series(np.random.rand(4), index=MultiIndex.from_product([[1, 2], [3, 4]]))
