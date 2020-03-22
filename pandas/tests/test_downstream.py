"""
Testing that we work in the downstream packages
"""
import importlib
import subprocess
import sys

import numpy as np  # noqa
import pytest

from pandas import DataFrame
import pandas._testing as tm
from pandas.wesm import dataframe as dataframe_protocol


def import_module(name):
    # we *only* want to skip if the module is truly not available
    # and NOT just an actual import error because of pandas changes

    try:
        return importlib.import_module(name)
    except ModuleNotFoundError:  # noqa
        pytest.skip(f"skipping as {name} not available")


@pytest.fixture
def df():
    return DataFrame({"A": [1, 2, 3]})


def test_dask(df):

    toolz = import_module("toolz")  # noqa
    dask = import_module("dask")  # noqa

    import dask.dataframe as dd

    ddf = dd.from_pandas(df, npartitions=3)
    assert ddf.A is not None
    assert ddf.compute() is not None


@pytest.mark.filterwarnings("ignore:Panel class is removed")
def test_xarray(df):

    xarray = import_module("xarray")  # noqa

    assert df.to_xarray() is not None


def test_oo_optimizable():
    # GH 21071
    subprocess.check_call([sys.executable, "-OO", "-c", "import pandas"])


@tm.network
# Cython import warning
@pytest.mark.filterwarnings("ignore:can't:ImportWarning")
@pytest.mark.filterwarnings(
    # patsy needs to update their imports
    "ignore:Using or importing the ABCs from 'collections:DeprecationWarning"
)
def test_statsmodels():

    statsmodels = import_module("statsmodels")  # noqa
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    df = sm.datasets.get_rdataset("Guerry", "HistData").data
    smf.ols("Lottery ~ Literacy + np.log(Pop1831)", data=df).fit()


# Cython import warning
@pytest.mark.filterwarnings("ignore:can't:ImportWarning")
def test_scikit_learn(df):

    sklearn = import_module("sklearn")  # noqa
    from sklearn import svm, datasets

    digits = datasets.load_digits()
    clf = svm.SVC(gamma=0.001, C=100.0)
    clf.fit(digits.data[:-1], digits.target[:-1])
    clf.predict(digits.data[-1:])


# Cython import warning and traitlets
@tm.network
@pytest.mark.filterwarnings("ignore")
def test_seaborn():

    seaborn = import_module("seaborn")
    tips = seaborn.load_dataset("tips")
    seaborn.stripplot(x="day", y="total_bill", data=tips)


def test_pandas_gbq(df):

    pandas_gbq = import_module("pandas_gbq")  # noqa


@pytest.mark.xfail(reason="0.7.0 pending")
@tm.network
def test_pandas_datareader():

    pandas_datareader = import_module("pandas_datareader")  # noqa
    pandas_datareader.DataReader("F", "quandl", "2017-01-01", "2017-02-01")


# importing from pandas, Cython import warning
@pytest.mark.filterwarnings("ignore:can't resolve:ImportWarning")
def test_geopandas():

    geopandas = import_module("geopandas")  # noqa
    fp = geopandas.datasets.get_path("naturalearth_lowres")
    assert geopandas.read_file(fp) is not None


# Cython import warning
@pytest.mark.filterwarnings("ignore:can't resolve:ImportWarning")
@pytest.mark.filterwarnings("ignore:RangeIndex.* is deprecated:DeprecationWarning")
def test_pyarrow(df):

    pyarrow = import_module("pyarrow")  # noqa
    table = pyarrow.Table.from_pandas(df)
    result = table.to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.xfail(reason="pandas-wheels-50", strict=False)
def test_missing_required_dependency():
    # GH 23868
    # To ensure proper isolation, we pass these flags
    # -S : disable site-packages
    # -s : disable user site-packages
    # -E : disable PYTHON* env vars, especially PYTHONPATH
    # And, that's apparently not enough, so we give up.
    # https://github.com/MacPython/pandas-wheels/pull/50
    call = ["python", "-sSE", "-c", "import pandas"]

    msg = (
        r"Command '\['python', '-sSE', '-c', 'import pandas'\]' "
        "returned non-zero exit status 1."
    )

    with pytest.raises(subprocess.CalledProcessError, match=msg) as exc:
        subprocess.check_output(call, stderr=subprocess.STDOUT)

    output = exc.value.stdout.decode()
    for name in ["numpy", "pytz", "dateutil"]:
        assert name in output


# -----------------------------------------------------------------------------
# DataFrame interchange protocol
# -----------------------------------------------------------------------------


class TestDataFrameProtocol:
    def test_interface_smoketest(self):
        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        result = df.__dataframe__
        assert isinstance(result, dataframe_protocol.DataFrame)
        assert isinstance(result["a"], dataframe_protocol.Column)
        assert isinstance(result.column_by_index(0), dataframe_protocol.Column)
        # assert isinstance(result['a'].dtype, dataframe_protocol.DataType)

        assert result.num_rows == 3
        assert result.num_columns == 2
        assert result.column_names == ["a", "b"]
        assert list(result.iter_column_names()) == ["a", "b"]

        expected = np.array([1, 2, 3], dtype=np.int64)
        res = result["a"].to_numpy()
        tm.assert_numpy_array_equal(res, expected)
        res = result.column_by_index(0).to_numpy()
        tm.assert_numpy_array_equal(res, expected)

        assert result["a"].name == "a"
        assert result.column_by_index(0).name == 'a'

    def test_pandas_dataframe_constructor(self):
        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        result = DataFrame(df)
        tm.assert_frame_equal(result, df)
        assert result is not df

        result = DataFrame(df.__dataframe__)
        tm.assert_frame_equal(result, df)
        assert result is not df
