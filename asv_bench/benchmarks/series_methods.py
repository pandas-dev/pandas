from pandas_vb_common import *


class series_isin_int64(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_isin_int64(self):
        self.s3.isin(self.values)


class series_isin_object(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
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
        self.s1.nlargest(3, take_last=True)
        self.s1.nlargest(3, take_last=False)


class series_nlargest2(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_nlargest2(self):
        self.s2.nlargest(3, take_last=True)
        self.s2.nlargest(3, take_last=False)


class series_nsmallest2(object):
    goal_time = 0.2

    def setup(self):
        self.s1 = Series(np.random.randn(10000))
        self.s2 = Series(np.random.randint(1, 10, 10000))
        self.s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
        self.values = [1, 2]
        self.s4 = self.s3.astype('object')

    def time_series_nsmallest2(self):
        self.s2.nsmallest(3, take_last=True)
        self.s2.nsmallest(3, take_last=False)