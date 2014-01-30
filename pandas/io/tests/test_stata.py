# pylint: disable=E1101

from datetime import datetime
import os
import warnings
import nose
import sys
from distutils.version import LooseVersion

import numpy as np

import pandas as pd
from pandas.core.frame import DataFrame, Series
from pandas.io.parsers import read_csv
from pandas.io.stata import read_stata, StataReader
import pandas.util.testing as tm
from pandas.util.misc import is_little_endian
from pandas import compat

def skip_if_not_little_endian():
    if not is_little_endian():
        raise nose.SkipTest("known failure of test on non-little endian")

class TestStata(tm.TestCase):

    def setUp(self):
        # Unit test datasets for dta7 - dta9 (old stata formats 104, 105 and 107) can be downloaded from:
        # http://stata-press.com/data/glmext.html
        self.dirpath = tm.get_data_path()
        self.dta1 = os.path.join(self.dirpath, 'stata1.dta')
        self.dta2 = os.path.join(self.dirpath, 'stata2.dta')
        self.dta3 = os.path.join(self.dirpath, 'stata3.dta')
        self.csv3 = os.path.join(self.dirpath, 'stata3.csv')
        self.dta4 = os.path.join(self.dirpath, 'stata4.dta')
        self.dta7 = os.path.join(self.dirpath, 'cancer.dta')
        self.csv7 = os.path.join(self.dirpath, 'cancer.csv')
        self.dta8 = os.path.join(self.dirpath, 'tbl19-3.dta')
        self.csv8 = os.path.join(self.dirpath, 'tbl19-3.csv')
        self.dta9 = os.path.join(self.dirpath, 'lbw.dta')
        self.csv9 = os.path.join(self.dirpath, 'lbw.csv')
        self.dta_encoding = os.path.join(self.dirpath, 'stata1_encoding.dta')
        self.dta1_13 = os.path.join(self.dirpath, 'stata1_v13.dta')
        self.dta2_13 = os.path.join(self.dirpath, 'stata2_v13.dta')
        self.dta3_13 = os.path.join(self.dirpath, 'stata3_v13.dta')
        self.dta4_13 = os.path.join(self.dirpath, 'stata4_v13.dta')

    def read_dta(self, file):
        return read_stata(file, convert_dates=True)

    def read_csv(self, file):
        return read_csv(file, parse_dates=True)

    def test_read_dta1(self):
        reader = StataReader(self.dta1)
        parsed = reader.data()
        reader_13 = StataReader(self.dta1_13)
        parsed_13 = reader_13.data()
        # Pandas uses np.nan as missing value.
        # Thus, all columns will be of type float, regardless of their name.
        expected = DataFrame([(np.nan, np.nan, np.nan, np.nan, np.nan)],
                             columns=['float_miss', 'double_miss', 'byte_miss',
                                      'int_miss', 'long_miss'])

        # this is an oddity as really the nan should be float64, but
        # the casting doesn't fail so need to match stata here
        expected['float_miss'] = expected['float_miss'].astype(np.float32)

        tm.assert_frame_equal(parsed, expected)
        tm.assert_frame_equal(parsed_13, expected)

    def test_read_dta2(self):
        if LooseVersion(sys.version) < '2.7':
            raise nose.SkipTest('datetime interp under 2.6 is faulty')

        expected = DataFrame.from_records(
            [
                (
                    datetime(2006, 11, 19, 23, 13, 20),
                    1479596223000,
                    datetime(2010, 1, 20),
                    datetime(2010, 1, 8),
                    datetime(2010, 1, 1),
                    datetime(1974, 7, 1),
                    datetime(2010, 1, 1),
                    datetime(2010, 1, 1)
                ),
                (
                    datetime(1959, 12, 31, 20, 3, 20),
                    -1479590,
                    datetime(1953, 10, 2),
                    datetime(1948, 6, 10),
                    datetime(1955, 1, 1),
                    datetime(1955, 7, 1),
                    datetime(1955, 1, 1),
                    datetime(2, 1, 1)
                ),
                (
                    pd.NaT,
                    pd.NaT,
                    pd.NaT,
                    pd.NaT,
                    pd.NaT,
                    pd.NaT,
                    pd.NaT,
                    pd.NaT,
                )
            ],
            columns=['datetime_c', 'datetime_big_c', 'date', 'weekly_date',
                     'monthly_date', 'quarterly_date', 'half_yearly_date',
                     'yearly_date']
        )

        with warnings.catch_warnings(record=True) as w:
            parsed = self.read_dta(self.dta2)
            parsed_13 = self.read_dta(self.dta2_13)
            np.testing.assert_equal(
                len(w), 1)  # should get a warning for that format.

        # buggy test because of the NaT comparison on certain platforms
        #
        #tm.assert_frame_equal(parsed, expected)
        #tm.assert_frame_equal(parsed_13, expected)

    def test_read_dta3(self):
        parsed = self.read_dta(self.dta3)
        parsed_13 = self.read_dta(self.dta3_13)

        # match stata here
        expected = self.read_csv(self.csv3)
        expected = expected.astype(np.float32)
        expected['year'] = expected['year'].astype(np.int32)
        expected['quarter'] = expected['quarter'].astype(np.int16)

        tm.assert_frame_equal(parsed, expected)
        tm.assert_frame_equal(parsed_13, expected)

    def test_read_dta4(self):
        parsed = self.read_dta(self.dta4)
        parsed_13 = self.read_dta(self.dta4_13)
        expected = DataFrame.from_records(
            [
                ["one", "ten", "one", "one", "one"],
                ["two", "nine", "two", "two", "two"],
                ["three", "eight", "three", "three", "three"],
                ["four", "seven", 4, "four", "four"],
                ["five", "six", 5, np.nan, "five"],
                ["six", "five", 6, np.nan, "six"],
                ["seven", "four", 7, np.nan, "seven"],
                ["eight", "three", 8, np.nan, "eight"],
                ["nine", "two", 9, np.nan, "nine"],
                ["ten", "one", "ten", np.nan, "ten"]
            ],
            columns=['fully_labeled', 'fully_labeled2', 'incompletely_labeled',
                     'labeled_with_missings', 'float_labelled'])

        tm.assert_frame_equal(parsed, expected)
        tm.assert_frame_equal(parsed_13, expected)

    def test_read_write_dta5(self):
        skip_if_not_little_endian()

        original = DataFrame([(np.nan, np.nan, np.nan, np.nan, np.nan)],
                             columns=['float_miss', 'double_miss', 'byte_miss',
                                      'int_miss', 'long_miss'])
        original.index.name = 'index'

        with tm.ensure_clean() as path:
            original.to_stata(path, None, False)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  original)

    def test_write_dta6(self):
        skip_if_not_little_endian()

        original = self.read_csv(self.csv3)
        original.index.name = 'index'

        with tm.ensure_clean() as path:
            original.to_stata(path, None, False)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  original)

    @nose.tools.nottest
    def test_read_dta7(self):
        expected = read_csv(self.csv7, parse_dates=True, sep='\t')
        parsed = self.read_dta(self.dta7)
        tm.assert_frame_equal(parsed, expected)

    @nose.tools.nottest
    def test_read_dta8(self):
        expected = read_csv(self.csv8, parse_dates=True, sep='\t')
        parsed = self.read_dta(self.dta8)
        tm.assert_frame_equal(parsed, expected)

    @nose.tools.nottest
    def test_read_dta9(self):
        expected = read_csv(self.csv9, parse_dates=True, sep='\t')
        parsed = self.read_dta(self.dta9)
        tm.assert_frame_equal(parsed, expected)

    def test_read_write_dta10(self):
        skip_if_not_little_endian()

        original = DataFrame(data=[["string", "object", 1, 1.1,
                                    np.datetime64('2003-12-25')]],
                             columns=['string', 'object', 'integer', 'float',
                                      'datetime'])
        original["object"] = Series(original["object"], dtype=object)
        original.index.name = 'index'

        with tm.ensure_clean() as path:
            original.to_stata(path, {'datetime': 'tc'}, False)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  original)

    def test_stata_doc_examples(self):
        with tm.ensure_clean() as path:
            df = DataFrame(np.random.randn(10, 2), columns=list('AB'))
            df.to_stata(path)

    def test_encoding(self):

        # GH 4626, proper encoding handling
        raw = read_stata(self.dta_encoding)
        encoded = read_stata(self.dta_encoding, encoding="latin-1")
        result = encoded.kreis1849[0]

        if compat.PY3:
            expected = raw.kreis1849[0]
            self.assert_(result == expected)
            self.assert_(isinstance(result, compat.string_types))
        else:
            expected = raw.kreis1849.str.decode("latin-1")[0]
            self.assert_(result == expected)
            self.assert_(isinstance(result, unicode))

    def test_read_write_dta11(self):
        skip_if_not_little_endian()

        original = DataFrame([(1, 2, 3, 4)],
                             columns=['good', compat.u('b\u00E4d'), '8number', 'astringwithmorethan32characters______'])
        formatted = DataFrame([(1, 2, 3, 4)],
                              columns=['good', 'b_d', '_8number', 'astringwithmorethan32characters_'])
        formatted.index.name = 'index'

        with tm.ensure_clean() as path:
            with warnings.catch_warnings(record=True) as w:
                original.to_stata(path, None, False)
                np.testing.assert_equal(
                    len(w), 1)  # should get a warning for that format.

            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'), formatted)

    def test_read_write_dta12(self):
        skip_if_not_little_endian()

        original = DataFrame([(1, 2, 3, 4)],
                             columns=['astringwithmorethan32characters_1', 'astringwithmorethan32characters_2', '+', '-'])
        formatted = DataFrame([(1, 2, 3, 4)],
                              columns=['astringwithmorethan32characters_', '_0astringwithmorethan32character', '_', '_1_'])
        formatted.index.name = 'index'

        with tm.ensure_clean() as path:
            with warnings.catch_warnings(record=True) as w:
                original.to_stata(path, None, False)
                np.testing.assert_equal(
                    len(w), 1)  # should get a warning for that format.

            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'), formatted)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
