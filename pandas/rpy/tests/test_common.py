"""
Testing that functions from rpy work as expected
"""

# flake8: noqa

import pandas as pd
import numpy as np
import unittest
import nose
import warnings
import pandas.util.testing as tm

try:
    import pandas.rpy.common as com
    from rpy2.robjects import r
    import rpy2.robjects as robj
except ImportError:
    raise nose.SkipTest('R not installed')


class TestCommon(unittest.TestCase):
    def test_convert_list(self):
        obj = r('list(a=1, b=2, c=3)')

        converted = com.convert_robj(obj)
        expected = {'a': [1], 'b': [2], 'c': [3]}

        tm.assert_dict_equal(converted, expected)

    def test_convert_nested_list(self):
        obj = r('list(a=list(foo=1, bar=2))')

        converted = com.convert_robj(obj)
        expected = {'a': {'foo': [1], 'bar': [2]}}

        tm.assert_dict_equal(converted, expected)

    def test_convert_frame(self):
        # built-in dataset
        df = r['faithful']

        converted = com.convert_robj(df)

        assert np.array_equal(converted.columns, ['eruptions', 'waiting'])
        assert np.array_equal(converted.index, np.arange(1, 273))

    def _test_matrix(self):
        r('mat <- matrix(rnorm(9), ncol=3)')
        r('colnames(mat) <- c("one", "two", "three")')
        r('rownames(mat) <- c("a", "b", "c")')

        return r['mat']

    def test_convert_matrix(self):
        mat = self._test_matrix()

        converted = com.convert_robj(mat)

        assert np.array_equal(converted.index, ['a', 'b', 'c'])
        assert np.array_equal(converted.columns, ['one', 'two', 'three'])

    def test_convert_r_dataframe(self):

        is_na = robj.baseenv.get("is.na")

        seriesd = tm.getSeriesData()
        frame = pd.DataFrame(seriesd, columns=['D', 'C', 'B', 'A'])

        # Null data
        frame["E"] = [np.nan for item in frame["A"]]
        # Some mixed type data
        frame["F"] = ["text" if item %
                      2 == 0 else np.nan for item in range(30)]

        r_dataframe = com.convert_to_r_dataframe(frame)

        assert np.array_equal(
            com.convert_robj(r_dataframe.rownames), frame.index)
        assert np.array_equal(
            com.convert_robj(r_dataframe.colnames), frame.columns)
        assert all(is_na(item) for item in r_dataframe.rx2("E"))

        for column in frame[["A", "B", "C", "D"]]:
            coldata = r_dataframe.rx2(column)
            original_data = frame[column]
            assert np.array_equal(com.convert_robj(coldata), original_data)

        for column in frame[["D", "E"]]:
            for original, converted in zip(frame[column],
                                           r_dataframe.rx2(column)):

                if pd.isnull(original):
                    assert is_na(converted)
                else:
                    assert original == converted

    def test_convert_r_matrix(self):

        is_na = robj.baseenv.get("is.na")

        seriesd = tm.getSeriesData()
        frame = pd.DataFrame(seriesd, columns=['D', 'C', 'B', 'A'])
        # Null data
        frame["E"] = [np.nan for item in frame["A"]]

        r_dataframe = com.convert_to_r_matrix(frame)

        assert np.array_equal(
            com.convert_robj(r_dataframe.rownames), frame.index)
        assert np.array_equal(
            com.convert_robj(r_dataframe.colnames), frame.columns)
        assert all(is_na(item) for item in r_dataframe.rx(True, "E"))

        for column in frame[["A", "B", "C", "D"]]:
            coldata = r_dataframe.rx(True, column)
            original_data = frame[column]
            assert np.array_equal(com.convert_robj(coldata),
                                  original_data)

        # Pandas bug 1282
        frame["F"] = ["text" if item %
                      2 == 0 else np.nan for item in range(30)]

        try:
            wrong_matrix = com.convert_to_r_matrix(frame)
        except TypeError:
            pass
        except Exception:
            raise

    def test_dist(self):
        for name in ('eurodist',):
            df = com.load_data(name)
            dist = r[name]
            labels = r['labels'](dist)
            assert np.array_equal(df.index, labels)
            assert np.array_equal(df.columns, labels)

    def test_timeseries(self):
        """
        Test that the series has an informative index.
        Unfortunately the code currently does not build a DateTimeIndex
        """
        for name in (
            'austres', 'co2', 'fdeaths', 'freeny.y', 'JohnsonJohnson',
            'ldeaths', 'mdeaths', 'nottem', 'presidents', 'sunspot.month', 'sunspots',
            'UKDriverDeaths', 'UKgas', 'USAccDeaths',
            'airmiles', 'discoveries', 'EuStockMarkets',
            'LakeHuron', 'lh', 'lynx', 'nhtemp', 'Nile',
                'Seatbelts', 'sunspot.year', 'treering', 'uspop'):
            series = com.load_data(name)
            ts = r[name]
            assert np.array_equal(series.index, r['time'](ts))

    def test_numeric(self):
        for name in ('euro', 'islands', 'precip'):
            series = com.load_data(name)
            numeric = r[name]
            names = numeric.names
            assert np.array_equal(series.index, names)

    def test_table(self):
        iris3 = pd.DataFrame({'X0': {0: '0', 1: '1', 2: '2', 3: '3', 4: '4'},
                              'X1': {0: 'Sepal L.',
                                     1: 'Sepal L.',
                                     2: 'Sepal L.',
                                     3: 'Sepal L.',
                                     4: 'Sepal L.'},
                              'X2': {0: 'Setosa',
                                     1: 'Setosa',
                                     2: 'Setosa',
                                     3: 'Setosa',
                                     4: 'Setosa'},
                              'value': {0: '5.1', 1: '4.9', 2: '4.7', 3: '4.6', 4: '5.0'}})
        hec = pd.DataFrame(
            {
                'Eye': {0: 'Brown', 1: 'Brown', 2: 'Brown', 3: 'Brown', 4: 'Blue'},
                'Hair': {0: 'Black', 1: 'Brown', 2: 'Red', 3: 'Blond', 4: 'Black'},
                'Sex': {0: 'Male', 1: 'Male', 2: 'Male', 3: 'Male', 4: 'Male'},
                'value': {0: '32.0', 1: '53.0', 2: '10.0', 3: '3.0', 4: '11.0'}})
        titanic = pd.DataFrame(
            {
                'Age': {0: 'Child', 1: 'Child', 2: 'Child', 3: 'Child', 4: 'Child'},
                'Class': {0: '1st', 1: '2nd', 2: '3rd', 3: 'Crew', 4: '1st'},
                'Sex': {0: 'Male', 1: 'Male', 2: 'Male', 3: 'Male', 4: 'Female'},
                'Survived': {0: 'No', 1: 'No', 2: 'No', 3: 'No', 4: 'No'},
                'value': {0: '0.0', 1: '0.0', 2: '35.0', 3: '0.0', 4: '0.0'}})
        for name, expected in zip(('HairEyeColor', 'Titanic', 'iris3'),
                                 (hec, titanic, iris3)):
            df = com.load_data(name)
            table = r[name]
            names = r['dimnames'](table)
            try:
                columns = list(r['names'](names))[::-1]
            except TypeError:
                columns = ['X{:d}'.format(i) for i in range(len(names))][::-1]
            columns.append('value')
            assert np.array_equal(df.columns, columns)
            result = df.head()
            cond = ((result.sort(axis=1) == expected.sort(axis=1))).values
            assert np.all(cond)

    def test_factor(self):
        for name in ('state.division', 'state.region'):
            vector = r[name]
            factors = list(r['factor'](vector))
            level = list(r['levels'](vector))
            factors = [level[index - 1] for index in factors]
            result = com.load_data(name)
            assert np.equal(result, factors)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
