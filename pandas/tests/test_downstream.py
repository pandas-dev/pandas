"""
Testing that we work in the downstream packages
"""
import pytest
import numpy as np  # noqa
from pandas import DataFrame
from pandas.util import testing as tm


@pytest.fixture
def df():
    return DataFrame({'A': [1, 2, 3]})


def test_dask(df):

    toolz = pytest.importorskip('toolz')  # noqa
    dask = pytest.importorskip('dask')  # noqa

    import dask.dataframe as dd

    ddf = dd.from_pandas(df, npartitions=3)
    assert ddf.A is not None
    assert ddf.compute() is not None


def test_xarray(df):

    xarray = pytest.importorskip('xarray')  # noqa

    assert df.to_xarray() is not None


def test_statsmodels():

    statsmodels = pytest.importorskip('statsmodels')  # noqa
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    df = sm.datasets.get_rdataset("Guerry", "HistData").data
    smf.ols('Lottery ~ Literacy + np.log(Pop1831)', data=df).fit()


def test_scikit_learn(df):

    sklearn = pytest.importorskip('sklearn')  # noqa
    from sklearn import svm, datasets

    digits = datasets.load_digits()
    clf = svm.SVC(gamma=0.001, C=100.)
    clf.fit(digits.data[:-1], digits.target[:-1])
    clf.predict(digits.data[-1:])


def test_seaborn():

    seaborn = pytest.importorskip('seaborn')
    tips = seaborn.load_dataset("tips")
    seaborn.stripplot(x="day", y="total_bill", data=tips)


def test_pandas_gbq(df):

    pandas_gbq = pytest.importorskip('pandas-gbq')  # noqa


@tm.network
def test_pandas_datareader():

    pandas_datareader = pytest.importorskip('pandas-datareader')  # noqa
    pandas_datareader.get_data_yahoo('AAPL')


def test_geopandas():

    geopandas = pytest.importorskip('geopandas')  # noqa
    fp = geopandas.datasets.get_path('naturalearth_lowres')
    assert geopandas.read_file(fp) is not None


def test_pyarrow(df):

    pyarrow = pytest.importorskip('pyarrow')  # noqa
    table = pyarrow.Table.from_pandas(df)
    result = table.to_pandas()
    tm.assert_frame_equal(result, df)
