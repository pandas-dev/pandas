import pytest
from warnings import catch_warnings

import numpy as np
import pandas as pd
from pandas import (Series, DataFrame, compat)
from .common import (ensure_clean_store, ensure_clean_path,
                     _check_roundtrip, _maybe_remove)

import pandas.util.testing as tm
from pandas.compat import is_platform_little_endian, range, u
from pandas.io.pytables import Term, read_hdf
from pandas.core.dtypes.common import is_categorical_dtype


@pytest.mark.skipif(not is_platform_little_endian(),
                    reason="reason platform is not little endian")
def test_encoding():
    with ensure_clean_store() as store:
        df = DataFrame(dict(A='foo', B='bar'), index=range(5))
        df.loc[2, 'A'] = np.nan
        df.loc[3, 'B'] = np.nan
        _maybe_remove(store, 'df')
        store.append('df', df, encoding='ascii')
        tm.assert_frame_equal(store['df'], df)

        expected = df.reindex(columns=['A'])
        result = store.select('df', Term('columns=A', encoding='ascii'))
        tm.assert_frame_equal(result, expected)


def test_latin_encoding():
    if compat.PY2:
        tm.assert_raises_regex(
            TypeError, r'\[unicode\] is not implemented as a table column')
        return

    values = [[b'E\xc9, 17', b'', b'a', b'b', b'c'],
              [b'E\xc9, 17', b'a', b'b', b'c'],
              [b'EE, 17', b'', b'a', b'b', b'c'],
              [b'E\xc9, 17', b'\xf8\xfc', b'a', b'b', b'c'],
              [b'', b'a', b'b', b'c'],
              [b'\xf8\xfc', b'a', b'b', b'c'],
              [b'A\xf8\xfc', b'', b'a', b'b', b'c'],
              [np.nan, b'', b'b', b'c'],
              [b'A\xf8\xfc', np.nan, b'', b'b', b'c']]

    def _try_decode(x, encoding='latin-1'):
        try:
            return x.decode(encoding)
        except AttributeError:
            return x
    # not sure how to remove latin-1 from code in python 2 and 3
    values = [[_try_decode(x) for x in y] for y in values]

    examples = []
    for dtype in ['category', object]:
        for val in values:
            examples.append(pd.Series(val, dtype=dtype))

    def roundtrip(s, key='data', encoding='latin-1', nan_rep=''):
        with ensure_clean_path() as store:
            s.to_hdf(store, key, format='table', encoding=encoding,
                     nan_rep=nan_rep)
            retr = read_hdf(store, key)
            s_nan = s.replace(nan_rep, np.nan)
            if is_categorical_dtype(s_nan):
                assert is_categorical_dtype(retr)
                tm.assert_series_equal(s_nan, retr, check_dtype=False,
                                       check_categorical=False)
            else:
                tm.assert_series_equal(s_nan, retr)

    for s in examples:
        roundtrip(s)


def test_unicode_longer_encoded():
    # GH 11234
    char = '\u0394'
    df = pd.DataFrame({'A': [char]})
    with ensure_clean_store() as store:
        store.put('df', df, format='table', encoding='utf-8')
        result = store.get('df')
        tm.assert_frame_equal(result, df)

    df = pd.DataFrame({'A': ['a', char], 'B': ['b', 'b']})
    with ensure_clean_store() as store:
        store.put('df', df, format='table', encoding='utf-8')
        result = store.get('df')
        tm.assert_frame_equal(result, df)


def test_unicode_index():
    unicode_values = [u('\u03c3'), u('\u03c3\u03c3')]

    # PerformanceWarning
    with catch_warnings(record=True):
        s = Series(np.random.randn(len(unicode_values)), unicode_values)
        _check_roundtrip(s, tm.assert_series_equal)
