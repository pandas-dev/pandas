import numpy as np
import pytest

from pandas import (
    DataFrame,
    NaT,
    date_range,
)
import pandas._testing as tm


@pytest.fixture
def datetime_frame() -> DataFrame:
    """
    Fixture for DataFrame of floats with DatetimeIndex

    Columns are ['A', 'B', 'C', 'D']

                       A         B         C         D
    2000-01-03 -1.122153  0.468535  0.122226  1.693711
    2000-01-04  0.189378  0.486100  0.007864 -1.216052
    2000-01-05  0.041401 -0.835752 -0.035279 -0.414357
    2000-01-06  0.430050  0.894352  0.090719  0.036939
    2000-01-07 -0.620982 -0.668211 -0.706153  1.466335
    2000-01-10 -0.752633  0.328434 -0.815325  0.699674
    2000-01-11 -2.236969  0.615737 -0.829076 -1.196106
    ...              ...       ...       ...       ...
    2000-02-03  1.642618 -0.579288  0.046005  1.385249
    2000-02-04 -0.544873 -1.160962 -0.284071 -1.418351
    2000-02-07 -2.656149 -0.601387  1.410148  0.444150
    2000-02-08 -1.201881 -1.289040  0.772992 -1.445300
    2000-02-09  1.377373  0.398619  1.008453 -0.928207
    2000-02-10  0.473194 -0.636677  0.984058  0.511519
    2000-02-11 -0.965556  0.408313 -1.312844 -0.381948

    [30 rows x 4 columns]
    """
    return DataFrame(tm.getTimeSeriesData())


@pytest.fixture
def float_string_frame():
    """
    Fixture for DataFrame of floats and strings with index of unique strings

    Columns are ['A', 'B', 'C', 'D', 'foo'].

                       A         B         C         D  foo
    w3orJvq07g -1.594062 -1.084273 -1.252457  0.356460  bar
    PeukuVdmz2  0.109855 -0.955086 -0.809485  0.409747  bar
    ahp2KvwiM8 -1.533729 -0.142519 -0.154666  1.302623  bar
    3WSJ7BUCGd  2.484964  0.213829  0.034778 -2.327831  bar
    khdAmufk0U -0.193480 -0.743518 -0.077987  0.153646  bar
    LE2DZiFlrE -0.193566 -1.343194 -0.107321  0.959978  bar
    HJXSJhVn7b  0.142590  1.257603 -0.659409 -0.223844  bar
    ...              ...       ...       ...       ...  ...
    9a1Vypttgw -1.316394  1.601354  0.173596  1.213196  bar
    h5d1gVFbEy  0.609475  1.106738 -0.155271  0.294630  bar
    mK9LsTQG92  1.303613  0.857040 -1.019153  0.369468  bar
    oOLksd9gKH  0.558219 -0.134491 -0.289869 -0.951033  bar
    9jgoOjKyHg  0.058270 -0.496110 -0.413212 -0.852659  bar
    jZLDHclHAO  0.096298  1.267510  0.549206 -0.005235  bar
    lR0nxDp1C2 -2.119350 -0.794384  0.544118  0.145849  bar

    [30 rows x 5 columns]
    """
    df = DataFrame(tm.getSeriesData())
    df["foo"] = "bar"
    return df


@pytest.fixture
def mixed_float_frame():
    """
    Fixture for DataFrame of different float types with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].

                       A         B         C         D
    GI7bbDaEZe -0.237908 -0.246225 -0.468506  0.752993
    KGp9mFepzA -1.140809 -0.644046 -1.225586  0.801588
    VeVYLAb1l2 -1.154013 -1.677615  0.690430 -0.003731
    kmPME4WKhO  0.979578  0.998274 -0.776367  0.897607
    CPyopdXTiz  0.048119 -0.257174  0.836426  0.111266
    0kJZQndAj0  0.274357 -0.281135 -0.344238  0.834541
    tqdwQsaHG8 -0.979716 -0.519897  0.582031  0.144710
    ...              ...       ...       ...       ...
    7FhZTWILQj -2.906357  1.261039 -0.780273 -0.537237
    4pUDPM4eGq -2.042512 -0.464382 -0.382080  1.132612
    B8dUgUzwTi -1.506637 -0.364435  1.087891  0.297653
    hErlVYjVv9  1.477453 -0.495515 -0.713867  1.438427
    1BKN3o7YLs  0.127535 -0.349812 -0.881836  0.489827
    9S4Ekn7zga  1.445518 -2.095149  0.031982  0.373204
    xN1dNn6OV6  1.425017 -0.983995 -0.363281 -0.224502

    [30 rows x 4 columns]
    """
    df = DataFrame(tm.getSeriesData())
    df.A = df.A.astype("float32")
    df.B = df.B.astype("float32")
    df.C = df.C.astype("float16")
    df.D = df.D.astype("float64")
    return df


@pytest.fixture
def mixed_int_frame():
    """
    Fixture for DataFrame of different int types with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].

                A  B    C    D
    mUrCZ67juP  0  1    2    2
    rw99ACYaKS  0  1    0    0
    7QsEcpaaVU  0  1    1    1
    xkrimI2pcE  0  1    0    0
    dz01SuzoS8  0  1  255  255
    ccQkqOHX75 -1  1    0    0
    DN0iXaoDLd  0  1    0    0
    ...        .. ..  ...  ...
    Dfb141wAaQ  1  1  254  254
    IPD8eQOVu5  0  1    0    0
    CcaKulsCmv  0  1    0    0
    rIBa8gu7E5  0  1    0    0
    RP6peZmh5o  0  1    1    1
    NMb9pipQWQ  0  1    0    0
    PqgbJEzjib  0  1    3    3

    [30 rows x 4 columns]
    """
    df = DataFrame({k: v.astype(int) for k, v in tm.getSeriesData().items()})
    df.A = df.A.astype("int32")
    df.B = np.ones(len(df.B), dtype="uint64")
    df.C = df.C.astype("uint8")
    df.D = df.C.astype("int64")
    return df


@pytest.fixture
def timezone_frame():
    """
    Fixture for DataFrame of date_range Series with different time zones

    Columns are ['A', 'B', 'C']; some entries are missing

               A                         B                         C
    0 2013-01-01 2013-01-01 00:00:00-05:00 2013-01-01 00:00:00+01:00
    1 2013-01-02                       NaT                       NaT
    2 2013-01-03 2013-01-03 00:00:00-05:00 2013-01-03 00:00:00+01:00
    """
    df = DataFrame(
        {
            "A": date_range("20130101", periods=3),
            "B": date_range("20130101", periods=3, tz="US/Eastern"),
            "C": date_range("20130101", periods=3, tz="CET"),
        }
    )
    df.iloc[1, 1] = NaT
    df.iloc[1, 2] = NaT
    return df
