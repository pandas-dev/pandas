from pandas import *
from pandas.util.testing import rands

i, j, k = 7, 771, 5532

panel = Panel(np.random.randn(i, j, k),
              items=[rands(10) for _ in xrange(i)],
              major_axis=DateRange('1/1/2000', periods=j,
                                   offset=datetools.Minute()),
              minor_axis=[rands(10) for _ in xrange(k)])


store = HDFStore('test.h5')
store.put('test_panel', panel, table=True)

retrieved = store['test_panel']
