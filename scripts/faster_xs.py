import numpy as np

import pandas.util.testing as tm

from pandas.core.internals import _interleaved_dtype

df = tm.makeDataFrame()

df['E'] = 'foo'
df['F'] = 'foo'
df['G'] = 2
df['H'] = df['A'] > 0

blocks = df._data.blocks
items = df.columns
