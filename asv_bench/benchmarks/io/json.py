import numpy as np
import pandas.util.testing as tm
from pandas import DataFrame, date_range, timedelta_range, concat, read_json

from ..pandas_vb_common import BaseIO


class ReadJSON(BaseIO):

    fname = "__test__.json"
    params = (['split', 'index', 'records'], ['int', 'datetime'])
    param_names = ['orient', 'index']

    def setup(self, orient, index):
        N = 100000
        indexes = {'int': np.arange(N),
                   'datetime': date_range('20000101', periods=N, freq='H')}
        df = DataFrame(np.random.randn(N, 5),
                       columns=['float_{}'.format(i) for i in range(5)],
                       index=indexes[index])
        df.to_json(self.fname, orient=orient)

    def time_read_json(self, orient, index):
        read_json(self.fname, orient=orient)


class ReadJSONLines(BaseIO):

    fname = "__test_lines__.json"
    params = ['int', 'datetime']
    param_names = ['index']

    def setup(self, index):
        N = 100000
        indexes = {'int': np.arange(N),
                   'datetime': date_range('20000101', periods=N, freq='H')}
        df = DataFrame(np.random.randn(N, 5),
                       columns=['float_{}'.format(i) for i in range(5)],
                       index=indexes[index])
        df.to_json(self.fname, orient='records', lines=True)

    def time_read_json_lines(self, index):
        read_json(self.fname, orient='records', lines=True)

    def time_read_json_lines_concat(self, index):
        concat(read_json(self.fname, orient='records', lines=True,
                         chunksize=25000))

    def peakmem_read_json_lines(self, index):
        read_json(self.fname, orient='records', lines=True)

    def peakmem_read_json_lines_concat(self, index):
        concat(read_json(self.fname, orient='records', lines=True,
                         chunksize=25000))


class ToJSON(BaseIO):

    fname = "__test__.json"
    params = ['split', 'columns', 'index']
    param_names = ['orient']

    def setup(self, lines_orient):
        N = 10**5
        ncols = 5
        index = date_range('20000101', periods=N, freq='H')
        timedeltas = timedelta_range(start=1, periods=N, freq='s')
        datetimes = date_range(start=1, periods=N, freq='s')
        ints = np.random.randint(100000000, size=N)
        floats = np.random.randn(N)
        strings = tm.makeStringIndex(N)
        self.df = DataFrame(np.random.randn(N, ncols), index=np.arange(N))
        self.df_date_idx = DataFrame(np.random.randn(N, ncols), index=index)
        self.df_td_int_ts = DataFrame({'td_1': timedeltas,
                                       'td_2': timedeltas,
                                       'int_1': ints,
                                       'int_2': ints,
                                       'ts_1': datetimes,
                                       'ts_2': datetimes},
                                      index=index)
        self.df_int_floats = DataFrame({'int_1': ints,
                                        'int_2': ints,
                                        'int_3': ints,
                                        'float_1': floats,
                                        'float_2': floats,
                                        'float_3': floats},
                                       index=index)
        self.df_int_float_str = DataFrame({'int_1': ints,
                                           'int_2': ints,
                                           'float_1': floats,
                                           'float_2': floats,
                                           'str_1': strings,
                                           'str_2': strings},
                                          index=index)

    def time_floats_with_int_index(self, orient):
        self.df.to_json(self.fname, orient=orient)

    def time_floats_with_dt_index(self, orient):
        self.df_date_idx.to_json(self.fname, orient=orient)

    def time_delta_int_tstamp(self, orient):
        self.df_td_int_ts.to_json(self.fname, orient=orient)

    def time_float_int(self, orient):
        self.df_int_floats.to_json(self.fname, orient=orient)

    def time_float_int_str(self, orient):
        self.df_int_float_str.to_json(self.fname, orient=orient)

    def time_floats_with_int_idex_lines(self, orient):
        self.df.to_json(self.fname, orient='records', lines=True)

    def time_floats_with_dt_index_lines(self, orient):
        self.df_date_idx.to_json(self.fname, orient='records', lines=True)

    def time_delta_int_tstamp_lines(self, orient):
        self.df_td_int_ts.to_json(self.fname, orient='records', lines=True)

    def time_float_int_lines(self, orient):
        self.df_int_floats.to_json(self.fname, orient='records', lines=True)

    def time_float_int_str_lines(self, orient):
        self.df_int_float_str.to_json(self.fname, orient='records', lines=True)


from ..pandas_vb_common import setup  # noqa: F401
