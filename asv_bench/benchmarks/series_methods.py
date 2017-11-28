from .pandas_vb_common import *


class series_constructor_no_data_datetime_index(object):
    goal_time = 0.2

    def setup(self):
        self.dr = pd.date_range(
            start=datetime(2015,10,26),
            end=datetime(2016,1,1),
            freq='50s'
        )  # ~100k long

    def time_series_constructor_no_data_datetime_index(self):
        Series(data=None, index=self.dr)


class series_constructor_dict_data_datetime_index(object):
    goal_time = 0.2

    def setup(self):
        self.dr = pd.date_range(
            start=datetime(2015, 10, 26),
            end=datetime(2016, 1, 1),
            freq='50s'
        )  # ~100k long
        self.data = {d: v for d, v in zip(self.dr, range(len(self.dr)))}

    def time_series_constructor_no_data_datetime_index(self):
        Series(data=self.data, index=self.dr)


class series_isin_int64(object):
    goal_time = 0.2

    def setup(self):
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.s4 = Series(np.random.randint(1, 100, 10000000)).astype('int64')
        self.values = [1, 2]

    def time_series_isin_int64(self):
        self.s3.isin(self.values)

    def time_series_isin_int64_large(self):
        self.s4.isin(self.values)


class series_isin_object(object):
    goal_time = 0.2

    def setup(self):
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_isin_object(self):
        self.s4.isin(self.values)


class series_nlargest1(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_nlargest1(self):
        self.s1.nlargest(3, keep='last')
        self.s1.nlargest(3, keep='first')


class series_nlargest2(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_nlargest2(self):
        self.s2.nlargest(3, keep='last')
        self.s2.nlargest(3, keep='first')


class series_nsmallest2(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_nsmallest2(self):
        self.s2.nsmallest(3, keep='last')
        self.s2.nsmallest(3, keep='first')


class series_dropna_int64(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randint(1, 10, 1000000))

    def time_series_dropna_int64(self):
        self.s.dropna()


class series_dropna_datetime(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(pd.date_range('2000-01-01', freq='S', periods=1000000))
        self.s[np.random.randint(1, 1000000, 100)] = pd.NaT

    def time_series_dropna_datetime(self):
        self.s.dropna()


class series_map_dict(object):
    goal_time = 0.2

    def setup(self):
        map_size = 1000
        self.s = Series(np.random.randint(0, map_size, 10000))
        self.map_dict = {i: map_size - i for i in range(map_size)}

    def time_series_map_dict(self):
        self.s.map(self.map_dict)


class series_map_series(object):
    goal_time = 0.2

    def setup(self):
        map_size = 1000
        self.s = Series(np.random.randint(0, map_size, 10000))
        self.map_series = Series(map_size - np.arange(map_size))

    def time_series_map_series(self):
        self.s.map(self.map_series)


class series_clip(object):
    goal_time = 0.2

    def setup(self):
        self.s = pd.Series(np.random.randn(50))

    def time_series_dropna_datetime(self):
        self.s.clip(0, 1)
