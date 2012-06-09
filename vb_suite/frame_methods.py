from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# lookup

setup = common_setup + """
df = DataFrame(np.random.randn(10000, 8), columns=list('abcdefgh'))
df['foo'] = 'bar'

row_labels = list(df.index[::10])[:900]
col_labels = list(df.columns) * 100
row_labels_all = list(df.index) * len(df.columns)
col_labels_all = list(df.columns) * len(df.index)
"""

frame_fancy_lookup = Benchmark('df.lookup(row_labels, col_labels)', setup,
                               start_date=datetime(2012, 1, 12))

frame_fancy_lookup_all = Benchmark('df.lookup(row_labels_all, col_labels_all)',
                                   setup,
                                   start_date=datetime(2012, 1, 12))

#----------------------------------------------------------------------
# fillna in place

setup = common_setup + """
df = DataFrame(randn(10000, 100))
df.values[::2] = np.nan
"""

frame_fillna_inplace = Benchmark('df.fillna(0, inplace=True)', setup,
                                 start_date=datetime(2012, 4, 4))


#----------------------------------------------------------------------
# reindex both axes

setup = common_setup + """
df = DataFrame(randn(1000, 1000))
idx = range(400, 700)
"""

frame_reindex_axis0 = Benchmark('df.reindex(idx)', setup)

frame_reindex_axis1 = Benchmark('df.reindex(columns=idx)', setup)

frame_reindex_both_axes = Benchmark('df.reindex(index=idx, columns=idx)',
                                    setup, start_date=datetime(2011, 1, 1))

frame_reindex_both_axes_ix = Benchmark('df.ix[idx, idx]', setup,
                                       start_date=datetime(2011, 1, 1))

#----------------------------------------------------------------------
# boolean indexing

setup = common_setup + """
df = DataFrame(randn(10000, 100))
bool_arr = np.zeros(10000, dtype=bool)
bool_arr[:1000] = True
"""

frame_boolean_row_select = Benchmark('df[bool_arr]', setup,
                                     start_date=datetime(2011, 1, 1))
