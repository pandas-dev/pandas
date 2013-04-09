# pylint: disable=E1101

from datetime import datetime
import os
import unittest

import warnings
import nose

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.categorical import Categorical
from pandas.io.parsers import (read_stata, read_csv, StataReader, StataWriter)
import pandas.util.testing as tm


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


class StataTests(unittest.TestCase):

    def setUp(self):
        self.dirpath = curpath()
        self.dta1 = os.path.join(self.dirpath, 'stata1.dta')
        self.dta2 = os.path.join(self.dirpath, 'stata2.dta')
        self.dta3 = os.path.join(self.dirpath, 'stata3.dta')
        self.csv3 = os.path.join(self.dirpath, 'stata3.csv')
        self.dta4 = os.path.join(self.dirpath, 'stata4.dta')
        self.dta5 = os.path.join(self.dirpath, 'stata5.dta')
        self.dta6 = os.path.join(self.dirpath, 'stata6.dta')

    def read_dta(self, file):
        return read_stata(file, convert_dates=True)

    def read_csv(self, file):
        return read_csv(file, parse_dates=True)

    def test_read_dta1(self):
        reader = StataReader(self.dta1)
        parsed = reader.data()
        # Pandas uses np.nan as missing value. Thus, all columns will be of type float, regardless of their name.
        expected = DataFrame([(np.nan, np.nan, np.nan, np.nan, np.nan)],
                             columns=['float_miss', 'double_miss', 'byte_miss', 'int_miss', 'long_miss'])

        for i, col in enumerate(parsed.columns):
            np.testing.assert_almost_equal(
                parsed[col],
                expected[expected.columns[i]]
            )

    def test_read_dta2(self):
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
                    np.datetime64('NaT'),
                    np.datetime64('NaT'),
                    np.datetime64('NaT'),
                    np.datetime64('NaT'),
                    np.datetime64('NaT'),
                    np.datetime64('NaT'),
                    np.datetime64('NaT'),
                    np.datetime64('NaT')
                )
            ],
            columns=['datetime_c', 'datetime_big_c', 'date', 'weekly_date', 'monthly_date', 'quarterly_date', 'half_yearly_date', 'yearly_date']
        )

        with warnings.catch_warnings(record=True) as w:
            parsed = self.read_dta(self.dta2)
            np.testing.assert_equal(
                len(w), 1)  # should get a warning for that format.

        tm.assert_frame_equal(parsed, expected)

    def test_read_dta3(self):
        parsed = self.read_dta(self.dta3)
        expected = self.read_csv(self.csv3)
        for i, col in enumerate(parsed.columns):
            np.testing.assert_almost_equal(
                parsed[col],
                expected[expected.columns[i]],
                decimal=3
            )

    def test_read_dta4(self):
        parsed = self.read_dta(self.dta4)
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
            columns=['fully_labeled', 'fully_labeled2', 'incompletely_labeled', 'labeled_with_missings', 'float_labelled'])

        tm.assert_frame_equal(parsed, expected)

    def test_write_dta5(self):
        original = DataFrame([(np.nan, np.nan, np.nan, np.nan, np.nan)],
                             columns=['float_miss', 'double_miss', 'byte_miss', 'int_miss', 'long_miss'])

        writer = StataWriter(self.dta5, original, None, False)
        writer.write_file()

        written_and_read_again = self.read_dta(self.dta5)
        tm.assert_frame_equal(written_and_read_again, original)

    def test_write_dta6(self):
        original = self.read_csv(self.csv3)

        writer = StataWriter(self.dta6, original, None, False)
        writer.write_file()

        written_and_read_again = self.read_dta(self.dta6)
        tm.assert_frame_equal(written_and_read_again, original)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
