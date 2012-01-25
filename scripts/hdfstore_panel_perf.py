from pandas import *
from pandas.util.testing import rands

i, j, k = 1, 10000, 500

panel = Panel(np.random.randn(i, j, k),
              items=[rands(10) for _ in xrange(i)],
              major_axis=[rands(10) for _ in xrange(j)],
              minor_axis=[rands(10) for _ in xrange(k)])


store = HDFStore('test.h5')
store.put('test_panel', panel, table=True)

retrieved = store['test_panel']
