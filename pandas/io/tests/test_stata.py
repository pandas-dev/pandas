# pylint: disable=E1101

from datetime import datetime
import datetime as dt
import os
import warnings
import nose
import sys
from distutils.version import LooseVersion

import numpy as np

import pandas as pd
from pandas.core.frame import DataFrame, Series
from pandas.io.parsers import read_csv
from pandas.io.stata import (read_stata, StataReader, InvalidColumnName,
    PossiblePrecisionLoss)
import pandas.util.testing as tm
from pandas.util.misc import is_little_endian
from pandas import compat

class TestStata(tm.TestCase):

    def setUp(self):
        # Unit test datasets for dta7 - dta9 (old stata formats 104, 105 and 107) can be downloaded from:
        # http://stata-press.com/data/glmext.html
        self.dirpath = tm.get_data_path()
        self.dta1_114 = os.path.join(self.dirpath, 'stata1_114.dta')
        self.dta1_117 = os.path.join(self.dirpath, 'stata1_117.dta')

        self.dta2_113 = os.path.join(self.dirpath, 'stata2_113.dta')
        self.dta2_114 = os.path.join(self.dirpath, 'stata2_114.dta')
        self.dta2_115 = os.path.join(self.dirpath, 'stata2_115.dta')
        self.dta2_117 = os.path.join(self.dirpath, 'stata2_117.dta')

        self.dta3_113 = os.path.join(self.dirpath, 'stata3_113.dta')
        self.dta3_114 = os.path.join(self.dirpath, 'stata3_114.dta')
        self.dta3_115 = os.path.join(self.dirpath, 'stata3_115.dta')
        self.dta3_117 = os.path.join(self.dirpath, 'stata3_117.dta')
        self.csv3 = os.path.join(self.dirpath, 'stata3.csv')

        self.dta4_113 = os.path.join(self.dirpath, 'stata4_113.dta')
        self.dta4_114 = os.path.join(self.dirpath, 'stata4_114.dta')
        self.dta4_115 = os.path.join(self.dirpath, 'stata4_115.dta')
        self.dta4_117 = os.path.join(self.dirpath, 'stata4_117.dta')

        self.dta7 = os.path.join(self.dirpath, 'cancer.dta')
        self.csv7 = os.path.join(self.dirpath, 'cancer.csv')

        self.dta8 = os.path.join(self.dirpath, 'tbl19-3.dta')

        self.csv8 = os.path.join(self.dirpath, 'tbl19-3.csv')

        self.dta9 = os.path.join(self.dirpath, 'lbw.dta')
        self.csv9 = os.path.join(self.dirpath, 'lbw.csv')

        self.dta_encoding = os.path.join(self.dirpath, 'stata1_encoding.dta')

        self.csv14 = os.path.join(self.dirpath, 'stata5.csv')
        self.dta14_113 = os.path.join(self.dirpath, 'stata5_113.dta')
        self.dta14_114 = os.path.join(self.dirpath, 'stata5_114.dta')
        self.dta14_115 = os.path.join(self.dirpath, 'stata5_115.dta')

        self.csv15 = os.path.join(self.dirpath, 'stata6.csv')
        self.dta15_113 = os.path.join(self.dirpath, 'stata6_113.dta')
        self.dta15_114 = os.path.join(self.dirpath, 'stata6_114.dta')
        self.dta15_115 = os.path.join(self.dirpath, 'stata6_115.dta')

    def read_dta(self, file):
        return read_stata(file, convert_dates=True)

    def read_csv(self, file):
        return read_csv(file, parse_dates=True)

    def test_read_dta1(self):
        reader_114 = StataReader(self.dta1_114)
        parsed_114 = reader_114.data()
        reader_117 = StataReader(self.dta1_117)
        parsed_117 = reader_117.data()
        # Pandas uses np.nan as missing value.
        # Thus, all columns will be of type float, regardless of their name.
        expected = DataFrame([(np.nan, np.nan, np.nan, np.nan, np.nan)],
                             columns=['float_miss', 'double_miss', 'byte_miss',
                                      'int_miss', 'long_miss'])

        # this is an oddity as really the nan should be float64, but
        # the casting doesn't fail so need to match stata here
        expected['float_miss'] = expected['float_miss'].astype(np.float32)

        tm.assert_frame_equal(parsed_114, expected)
        tm.assert_frame_equal(parsed_117, expected)

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
        expected['yearly_date'] = expected['yearly_date'].astype('O')

        with warnings.catch_warnings(record=True) as w:
            parsed_114 = self.read_dta(self.dta2_114)
            parsed_115 = self.read_dta(self.dta2_115)
            parsed_117 = self.read_dta(self.dta2_117)
            # 113 is buggy due ot limits date format support in Stata
            # parsed_113 = self.read_dta(self.dta2_113)
            tm.assert_equal(
                len(w), 1)  # should get a warning for that format.

        # buggy test because of the NaT comparison on certain platforms
        # Format 113 test fails since it does not support tc and tC formats
        # tm.assert_frame_equal(parsed_113, expected)
        tm.assert_frame_equal(parsed_114, expected)
        tm.assert_frame_equal(parsed_115, expected)
        tm.assert_frame_equal(parsed_117, expected)

    def test_read_dta3(self):
        parsed_113 = self.read_dta(self.dta3_113)
        parsed_114 = self.read_dta(self.dta3_114)
        parsed_115 = self.read_dta(self.dta3_115)
        parsed_117 = self.read_dta(self.dta3_117)

        # match stata here
        expected = self.read_csv(self.csv3)
        expected = expected.astype(np.float32)
        expected['year'] = expected['year'].astype(np.int16)
        expected['quarter'] = expected['quarter'].astype(np.int8)

        tm.assert_frame_equal(parsed_113, expected)
        tm.assert_frame_equal(parsed_114, expected)
        tm.assert_frame_equal(parsed_115, expected)
        tm.assert_frame_equal(parsed_117, expected)

    def test_read_dta4(self):
        parsed_113 = self.read_dta(self.dta4_113)
        parsed_114 = self.read_dta(self.dta4_114)
        parsed_115 = self.read_dta(self.dta4_115)
        parsed_117 = self.read_dta(self.dta4_117)

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

        tm.assert_frame_equal(parsed_113, expected)
        tm.assert_frame_equal(parsed_114, expected)
        tm.assert_frame_equal(parsed_115, expected)
        tm.assert_frame_equal(parsed_117, expected)

    def test_read_write_dta5(self):
        original = DataFrame([(np.nan, np.nan, np.nan, np.nan, np.nan)],
                             columns=['float_miss', 'double_miss', 'byte_miss',
                                      'int_miss', 'long_miss'])
        original.index.name = 'index'

        with tm.ensure_clean() as path:
            original.to_stata(path, None)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  original)

    def test_write_dta6(self):
        original = self.read_csv(self.csv3)
        original.index.name = 'index'
        original.index = original.index.astype(np.int32)
        original['year'] = original['year'].astype(np.int32)
        original['quarter'] = original['quarter'].astype(np.int32)

        with tm.ensure_clean() as path:
            original.to_stata(path, None)
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
        original = DataFrame(data=[["string", "object", 1, 1.1,
                                    np.datetime64('2003-12-25')]],
                             columns=['string', 'object', 'integer', 'floating',
                                      'datetime'])
        original["object"] = Series(original["object"], dtype=object)
        original.index.name = 'index'
        original.index = original.index.astype(np.int32)
        original['integer'] = original['integer'].astype(np.int32)

        with tm.ensure_clean() as path:
            original.to_stata(path, {'datetime': 'tc'})
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
            self.assertEqual(result, expected)
            self.assertIsInstance(result, compat.string_types)
        else:
            expected = raw.kreis1849.str.decode("latin-1")[0]
            self.assertEqual(result, expected)
            self.assertIsInstance(result, unicode)

    def test_read_write_dta11(self):
        original = DataFrame([(1, 2, 3, 4)],
                             columns=['good', compat.u('b\u00E4d'), '8number', 'astringwithmorethan32characters______'])
        formatted = DataFrame([(1, 2, 3, 4)],
                              columns=['good', 'b_d', '_8number', 'astringwithmorethan32characters_'])
        formatted.index.name = 'index'
        formatted = formatted.astype(np.int32)

        with tm.ensure_clean() as path:
            with warnings.catch_warnings(record=True) as w:
                original.to_stata(path, None)
                # should get a warning for that format.
            tm.assert_equal(len(w), 1)

            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'), formatted)

    def test_read_write_dta12(self):
        original = DataFrame([(1, 2, 3, 4, 5, 6)],
                             columns=['astringwithmorethan32characters_1',
                                      'astringwithmorethan32characters_2',
                                      '+',
                                      '-',
                                      'short',
                                      'delete'])
        formatted = DataFrame([(1, 2, 3, 4, 5, 6)],
                              columns=['astringwithmorethan32characters_',
                                       '_0astringwithmorethan32character',
                                       '_',
                                       '_1_',
                                       '_short',
                                       '_delete'])
        formatted.index.name = 'index'
        formatted = formatted.astype(np.int32)

        with tm.ensure_clean() as path:
            with warnings.catch_warnings(record=True) as w:
                original.to_stata(path, None)
                tm.assert_equal(len(w), 1)  # should get a warning for that format.

            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'), formatted)

    def test_read_write_dta13(self):
        s1 = Series(2**9, dtype=np.int16)
        s2 = Series(2**17, dtype=np.int32)
        s3 = Series(2**33, dtype=np.int64)
        original = DataFrame({'int16': s1, 'int32': s2, 'int64': s3})
        original.index.name = 'index'

        formatted = original
        formatted['int64'] = formatted['int64'].astype(np.float64)

        with tm.ensure_clean() as path:
            original.to_stata(path)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  formatted)

    def test_read_write_reread_dta14(self):
        expected = self.read_csv(self.csv14)
        cols = ['byte_', 'int_', 'long_', 'float_', 'double_']
        for col in cols:
            expected[col] = expected[col].convert_objects(convert_numeric=True)
        expected['float_'] = expected['float_'].astype(np.float32)
        expected['date_td'] = pd.to_datetime(expected['date_td'], coerce=True)

        parsed_113 = self.read_dta(self.dta14_113)
        parsed_113.index.name = 'index'
        parsed_114 = self.read_dta(self.dta14_114)
        parsed_114.index.name = 'index'
        parsed_115 = self.read_dta(self.dta14_115)
        parsed_115.index.name = 'index'

        tm.assert_frame_equal(parsed_114, parsed_113)
        tm.assert_frame_equal(parsed_114, parsed_115)

        with tm.ensure_clean() as path:
            parsed_114.to_stata(path, {'date_td': 'td'})
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'), parsed_114)

    def test_read_write_reread_dta15(self):
        expected = self.read_csv(self.csv15)
        expected['byte_'] = expected['byte_'].astype(np.int8)
        expected['int_'] = expected['int_'].astype(np.int16)
        expected['long_'] = expected['long_'].astype(np.int32)
        expected['float_'] = expected['float_'].astype(np.float32)
        expected['double_'] = expected['double_'].astype(np.float64)
        expected['date_td'] = expected['date_td'].apply(datetime.strptime, args=('%Y-%m-%d',))

        parsed_113 = self.read_dta(self.dta15_113)
        parsed_114 = self.read_dta(self.dta15_114)
        parsed_115 = self.read_dta(self.dta15_115)

        tm.assert_frame_equal(expected, parsed_114)
        tm.assert_frame_equal(parsed_113, parsed_114)
        tm.assert_frame_equal(parsed_114, parsed_115)

    def test_timestamp_and_label(self):
        original = DataFrame([(1,)], columns=['var'])
        time_stamp = datetime(2000, 2, 29, 14, 21)
        data_label = 'This is a data file.'
        with tm.ensure_clean() as path:
            original.to_stata(path, time_stamp=time_stamp, data_label=data_label)
            reader = StataReader(path)
            parsed_time_stamp = dt.datetime.strptime(reader.time_stamp, ('%d %b %Y %H:%M'))
            assert parsed_time_stamp == time_stamp
            assert reader.data_label == data_label

    def test_numeric_column_names(self):
        original = DataFrame(np.reshape(np.arange(25.0), (5, 5)))
        original.index.name = 'index'
        with tm.ensure_clean() as path:
            # should get a warning for that format.
            with warnings.catch_warnings(record=True) as w:
                tm.assert_produces_warning(original.to_stata(path), InvalidColumnName)
            # should produce a single warning
            tm.assert_equal(len(w), 1)

            written_and_read_again = self.read_dta(path)
            written_and_read_again = written_and_read_again.set_index('index')
            columns = list(written_and_read_again.columns)
            convert_col_name = lambda x: int(x[1])
            written_and_read_again.columns = map(convert_col_name, columns)
            tm.assert_frame_equal(original, written_and_read_again)

    def test_nan_to_missing_value(self):
        s1 = Series(np.arange(4.0), dtype=np.float32)
        s2 = Series(np.arange(4.0), dtype=np.float64)
        s1[::2] = np.nan
        s2[1::2] = np.nan
        original = DataFrame({'s1': s1, 's2': s2})
        original.index.name = 'index'
        with tm.ensure_clean() as path:
            original.to_stata(path)
            written_and_read_again = self.read_dta(path)
            written_and_read_again = written_and_read_again.set_index('index')
            tm.assert_frame_equal(written_and_read_again, original)

    def test_no_index(self):
        columns = ['x', 'y']
        original = DataFrame(np.reshape(np.arange(10.0), (5, 2)),
                             columns=columns)
        original.index.name = 'index_not_written'
        with tm.ensure_clean() as path:
            original.to_stata(path, write_index=False)
            written_and_read_again = self.read_dta(path)
            tm.assertRaises(KeyError,
                            lambda: written_and_read_again['index_not_written'])

    def test_string_no_dates(self):
        s1 = Series(['a', 'A longer string'])
        s2 = Series([1.0, 2.0], dtype=np.float64)
        original = DataFrame({'s1': s1, 's2': s2})
        original.index.name = 'index'
        with tm.ensure_clean() as path:
            original.to_stata(path)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  original)

    def test_large_value_conversion(self):
        s0 = Series([1, 99], dtype=np.int8)
        s1 = Series([1, 127], dtype=np.int8)
        s2 = Series([1, 2 ** 15 - 1], dtype=np.int16)
        s3 = Series([1, 2 ** 63 - 1], dtype=np.int64)
        original = DataFrame({'s0': s0, 's1': s1, 's2': s2, 's3': s3})
        original.index.name = 'index'
        with tm.ensure_clean() as path:
            with warnings.catch_warnings(record=True) as w:
                tm.assert_produces_warning(original.to_stata(path),
                                           PossiblePrecisionLoss)
            # should produce a single warning
            tm.assert_equal(len(w), 1)

            written_and_read_again = self.read_dta(path)
            modified = original.copy()
            modified['s1'] = Series(modified['s1'], dtype=np.int16)
            modified['s2'] = Series(modified['s2'], dtype=np.int32)
            modified['s3'] = Series(modified['s3'], dtype=np.float64)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  modified)

    def test_dates_invalid_column(self):
        original = DataFrame([datetime(2006, 11, 19, 23, 13, 20)])
        original.index.name = 'index'
        with tm.ensure_clean() as path:
            with warnings.catch_warnings(record=True) as w:
                tm.assert_produces_warning(original.to_stata(path, {0: 'tc'}),
                                           InvalidColumnName)
            tm.assert_equal(len(w), 1)

            written_and_read_again = self.read_dta(path)
            modified = original.copy()
            modified.columns = ['_0']
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  modified)

    def test_date_export_formats(self):
        columns = ['tc', 'td', 'tw', 'tm', 'tq', 'th', 'ty']
        conversions = dict(((c, c) for c in columns))
        data = [datetime(2006, 11, 20, 23, 13, 20)] * len(columns)
        original = DataFrame([data], columns=columns)
        original.index.name = 'index'
        expected_values = [datetime(2006, 11, 20, 23, 13, 20),  # Time
                           datetime(2006, 11, 20),  # Day
                           datetime(2006, 11, 19),  # Week
                           datetime(2006, 11, 1),  # Month
                           datetime(2006, 10, 1),  # Quarter year
                           datetime(2006, 7, 1),  # Half year
                           datetime(2006, 1, 1)]  # Year

        expected = DataFrame([expected_values], columns=columns)
        expected.index.name = 'index'
        with tm.ensure_clean() as path:
            original.to_stata(path, conversions)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  expected)

    def test_write_missing_strings(self):
        original = DataFrame([["1"], [None]], columns=["foo"])
        expected = DataFrame([["1"], [""]], columns=["foo"])
        expected.index.name = 'index'
        with tm.ensure_clean() as path:
            original.to_stata(path)
            written_and_read_again = self.read_dta(path)
            tm.assert_frame_equal(written_and_read_again.set_index('index'),
                                  expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)

