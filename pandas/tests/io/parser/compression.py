# -*- coding: utf-8 -*-

"""
Tests compressed data parsing functionality for all
of the parsers defined in parsers.py
"""

import bz2
import gzip

import pytest

import pandas.compat as compat
import pandas.util._test_decorators as td

import pandas as pd
import pandas.util.testing as tm

try:
    lzma = compat.import_lzma()
except ImportError:
    lzma = None


class CompressionTests(object):

    def test_zip(self):
        import zipfile

        with open(self.csv1, 'rb') as data_file:
            data = data_file.read()
            expected = self.read_csv(self.csv1)

        with tm.ensure_clean('test_file.zip') as path:
            with zipfile.ZipFile(path, mode='w') as tmp:
                tmp.writestr('test_file', data)

            result = self.read_csv(path, compression='zip')
            tm.assert_frame_equal(result, expected)

            result = self.read_csv(path, compression='infer')
            tm.assert_frame_equal(result, expected)

            if self.engine is not 'python':
                with open(path, 'rb') as f:
                    result = self.read_csv(f, compression='zip')
                    tm.assert_frame_equal(result, expected)

        with tm.ensure_clean('combined_zip.zip') as path:
            inner_file_names = ['test_file', 'second_file']
            with zipfile.ZipFile(path, mode='w') as tmp:
                for file_name in inner_file_names:
                    tmp.writestr(file_name, data)

            tm.assert_raises_regex(ValueError, 'Multiple files',
                                   self.read_csv, path, compression='zip')

            tm.assert_raises_regex(ValueError, 'Multiple files',
                                   self.read_csv, path,
                                   compression='infer')

        with tm.ensure_clean() as path:
            with zipfile.ZipFile(path, mode='w') as tmp:
                pass

            tm.assert_raises_regex(ValueError, 'Zero files',
                                   self.read_csv, path, compression='zip')

        with tm.ensure_clean() as path:
            with open(path, 'wb') as f:
                pytest.raises(zipfile.BadZipfile, self.read_csv,
                              f, compression='zip')

    @pytest.mark.parametrize('compress_type, compress_method, ext', [
        ('gzip', gzip.GzipFile, 'gz'),
        ('bz2', bz2.BZ2File, 'bz2'),
        pytest.param('xz', getattr(lzma, 'LZMAFile', None), 'xz',
                     marks=td.skip_if_no_lzma)
    ])
    def test_other_compression(self, compress_type, compress_method, ext):

        with open(self.csv1, 'rb') as data_file:
            data = data_file.read()
            expected = self.read_csv(self.csv1)

        with tm.ensure_clean() as path:
            with compress_method(path, mode='wb') as tmp:
                tmp.write(data)

            result = self.read_csv(path, compression=compress_type)
            tm.assert_frame_equal(result, expected)

            if compress_type == 'bz2':
                pytest.raises(ValueError, self.read_csv,
                              path, compression='bz3')

            with open(path, 'rb') as fin:
                result = self.read_csv(fin, compression=compress_type)
                tm.assert_frame_equal(result, expected)

        with tm.ensure_clean('test.{}'.format(ext)) as path:
            with compress_method(path, mode='wb') as tmp:
                tmp.write(data)
            result = self.read_csv(path, compression='infer')
            tm.assert_frame_equal(result, expected)

    def test_read_csv_infer_compression(self):
        # see gh-9770
        expected = self.read_csv(self.csv1, index_col=0, parse_dates=True)

        with open(self.csv1) as f:
            inputs = [self.csv1, self.csv1 + '.gz',
                      self.csv1 + '.bz2', f]

            for inp in inputs:
                df = self.read_csv(inp, index_col=0, parse_dates=True,
                                   compression='infer')

                tm.assert_frame_equal(expected, df)

    def test_read_csv_compressed_utf16_example(self, datapath):
        # GH18071
        path = datapath('io', 'parser', 'data', 'utf16_ex_small.zip')

        result = self.read_csv(path, encoding='utf-16',
                               compression='zip', sep='\t')
        expected = pd.DataFrame({
            u'Country': [u'Venezuela', u'Venezuela'],
            u'Twitter': [u'Hugo Chávez Frías', u'Henrique Capriles R.']
        })

        tm.assert_frame_equal(result, expected)

    def test_invalid_compression(self):
        msg = 'Unrecognized compression type: sfark'
        with tm.assert_raises_regex(ValueError, msg):
            self.read_csv('test_file.zip', compression='sfark')
