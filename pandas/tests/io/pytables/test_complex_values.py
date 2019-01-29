from warnings import catch_warnings

import numpy as np
import pytest

from pandas import DataFrame, Panel, Series, concat
from pandas.util.testing import (
    assert_frame_equal, assert_panel_equal, assert_series_equal)

from pandas.io.pytables import read_hdf

from .base import Base, ensure_clean_path, ensure_clean_store


class TestHDFComplexValues(Base):
    # GH10447

    def test_complex_fixed(self):
        df = DataFrame(np.random.rand(4, 5).astype(np.complex64),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

        df = DataFrame(np.random.rand(4, 5).astype(np.complex128),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_table(self):
        df = DataFrame(np.random.rand(4, 5).astype(np.complex64),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

        df = DataFrame(np.random.rand(4, 5).astype(np.complex128),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table', mode='w')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_mixed_fixed(self):
        complex64 = np.array([1.0 + 1.0j, 1.0 + 1.0j,
                              1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64)
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j],
                              dtype=np.complex128)
        df = DataFrame({'A': [1, 2, 3, 4],
                        'B': ['a', 'b', 'c', 'd'],
                        'C': complex64,
                        'D': complex128,
                        'E': [1.0, 2.0, 3.0, 4.0]},
                       index=list('abcd'))
        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_mixed_table(self):
        complex64 = np.array([1.0 + 1.0j, 1.0 + 1.0j,
                              1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64)
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j],
                              dtype=np.complex128)
        df = DataFrame({'A': [1, 2, 3, 4],
                        'B': ['a', 'b', 'c', 'd'],
                        'C': complex64,
                        'D': complex128,
                        'E': [1.0, 2.0, 3.0, 4.0]},
                       index=list('abcd'))

        with ensure_clean_store(self.path) as store:
            store.append('df', df, data_columns=['A', 'B'])
            result = store.select('df', where='A>2')
            assert_frame_equal(df.loc[df.A > 2], result)

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    @pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
    def test_complex_across_dimensions_fixed(self):
        with catch_warnings(record=True):
            complex128 = np.array(
                [1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
            s = Series(complex128, index=list('abcd'))
            df = DataFrame({'A': s, 'B': s})
            p = Panel({'One': df, 'Two': df})

            objs = [s, df, p]
            comps = [assert_series_equal, assert_frame_equal,
                     assert_panel_equal]
            for obj, comp in zip(objs, comps):
                with ensure_clean_path(self.path) as path:
                    obj.to_hdf(path, 'obj', format='fixed')
                    reread = read_hdf(path, 'obj')
                    comp(obj, reread)

    @pytest.mark.filterwarnings("ignore:\\nPanel:FutureWarning")
    def test_complex_across_dimensions(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
        s = Series(complex128, index=list('abcd'))
        df = DataFrame({'A': s, 'B': s})

        with catch_warnings(record=True):
            p = Panel({'One': df, 'Two': df})

            objs = [df, p]
            comps = [assert_frame_equal, assert_panel_equal]
            for obj, comp in zip(objs, comps):
                with ensure_clean_path(self.path) as path:
                    obj.to_hdf(path, 'obj', format='table')
                    reread = read_hdf(path, 'obj')
                    comp(obj, reread)

    def test_complex_indexing_error(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j],
                              dtype=np.complex128)
        df = DataFrame({'A': [1, 2, 3, 4],
                        'B': ['a', 'b', 'c', 'd'],
                        'C': complex128},
                       index=list('abcd'))
        with ensure_clean_store(self.path) as store:
            pytest.raises(TypeError, store.append,
                          'df', df, data_columns=['C'])

    def test_complex_series_error(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
        s = Series(complex128, index=list('abcd'))

        with ensure_clean_path(self.path) as path:
            pytest.raises(TypeError, s.to_hdf, path, 'obj', format='t')

        with ensure_clean_path(self.path) as path:
            s.to_hdf(path, 'obj', format='t', index=False)
            reread = read_hdf(path, 'obj')
            assert_series_equal(s, reread)

    def test_complex_append(self):
        df = DataFrame({'a': np.random.randn(100).astype(np.complex128),
                        'b': np.random.randn(100)})

        with ensure_clean_store(self.path) as store:
            store.append('df', df, data_columns=['b'])
            store.append('df', df)
            result = store.select('df')
            assert_frame_equal(concat([df, df], 0), result)
