from pandas import *
from pandas.util.testing import rands
from pandas.compat import range

i, j, k = 7, 771, 5532

panel = Panel(np.random.randn(i, j, k),
              items=[rands(10) for _ in range(i)],
              major_axis=DatetimeIndex('1/1/2000', periods=j,
                                       offset=offsets.Minute()),
              minor_axis=[rands(10) for _ in range(k)])


store = HDFStore('test.h5')
store.put('test_panel', panel, table=True)

retrieved = store['test_panel']
