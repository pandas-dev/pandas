# -*- coding: utf-8 -*-

import os
import nose

import pandas.util.testing as tm

from pandas import read_csv, read_table
from pandas.core.common import AbstractMethodError

from .common import ParserTests
from .header import HeaderTests
from .comment import CommentTests
from .usecols import UsecolsTests
from .skiprows import SkipRowsTests
from .index_col import IndexColTests
from .na_values import NAvaluesTests
from .converters import ConverterTests
from .c_parser_only import CParserTests
from .parse_dates import ParseDatesTests
from .compression import CompressionTests
from .multithread import MultithreadTests
from .python_parser_only import PythonParserTests


class BaseParser(CommentTests, CompressionTests,
                 ConverterTests, HeaderTests,
                 IndexColTests, MultithreadTests,
                 NAvaluesTests, ParseDatesTests,
                 ParserTests, SkipRowsTests,
                 UsecolsTests):
    def read_csv(self, *args, **kwargs):
        raise NotImplementedError

    def read_table(self, *args, **kwargs):
        raise NotImplementedError

    def float_precision_choices(self):
        raise AbstractMethodError(self)

    def setUp(self):
        self.dirpath = tm.get_data_path()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')


class TestCParserHighMemory(BaseParser, CParserTests, tm.TestCase):
    engine = 'c'
    low_memory = False
    float_precision_choices = [None, 'high', 'round_trip']

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = self.low_memory
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = self.low_memory
        return read_table(*args, **kwds)


class TestCParserLowMemory(BaseParser, CParserTests, tm.TestCase):
    engine = 'c'
    low_memory = True
    float_precision_choices = [None, 'high', 'round_trip']

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = self.low_memory
        kwds['buffer_lines'] = 2
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = True
        kwds['buffer_lines'] = 2
        return read_table(*args, **kwds)


class TestPythonParser(BaseParser, PythonParserTests, tm.TestCase):
    """
    Class for Python parser testing. Unless specifically stated
    as a PythonParser-specific issue, the goal is to eventually move
    as many of these tests into ParserTests as soon as the C parser
    can accept further specific arguments when parsing.
    """

    engine = 'python'
    float_precision_choices = [None]

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        return read_table(*args, **kwds)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
