import numpy as np
from pandas import DataFrame, date_range, read_msgpack
import pandas.util.testing as tm

from ..pandas_vb_common import BaseIO


class MSGPack(BaseIO):

    goal_time = 0.2

    def setup(self):
        self.fname = '__test__.msg'
        N = 100000
        C = 5
        self.df = DataFrame(np.random.randn(N, C),
                            columns=['float{}'.format(i) for i in range(C)],
                            index=date_range('20000101', periods=N, freq='H'))
        self.df['object'] = tm.makeStringIndex(N)
        self.df.to_msgpack(self.fname)

    def time_read_msgpack(self):
        read_msgpack(self.fname)

    def time_write_msgpack(self):
        self.df.to_msgpack(self.fname)


from ..pandas_vb_common import setup  # noqa: F401
