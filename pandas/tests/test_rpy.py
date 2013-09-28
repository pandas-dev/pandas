"""
Testing that functions from rpy work as expected
"""

import pandas as pd
import numpy as np
import unittest
import nose
import pandas.util.testing as tm
import pandas.rpy.common as com
from rpy2.robjects import r
import rpy2.robjects as robj

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
        frame["F"] = ["text" if item % 2 == 0 else np.nan for item in range(30)]

        r_dataframe = com.convert_to_r_dataframe(frame)

        assert np.array_equal(com.convert_robj(r_dataframe.rownames), frame.index)
        assert np.array_equal(com.convert_robj(r_dataframe.colnames), frame.columns)
        assert all(is_na(item) for item in r_dataframe.rx2("E"))

        for column in frame[["A", "B", "C", "D"]]:
            coldata = r_dataframe.rx2(column)
            original_data = frame[column]
            print(com.convert_robj(coldata))
            print(original_data)
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

        assert np.array_equal(com.convert_robj(r_dataframe.rownames), frame.index)
        assert np.array_equal(com.convert_robj(r_dataframe.colnames), frame.columns)
        assert all(is_na(item) for item in r_dataframe.rx(True, "E"))

        for column in frame[["A", "B", "C", "D"]]:
            coldata = r_dataframe.rx(True, column)
            original_data = frame[column]
            assert np.array_equal(com.convert_robj(coldata),
                                  original_data)

        # Pandas bug 1282
        frame["F"] = ["text" if item % 2 == 0 else np.nan for item in range(30)]

        try:
            wrong_matrix = com.convert_to_r_matrix(frame)
        except TypeError:
            pass
        except Exception:
            raise
    
if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)
