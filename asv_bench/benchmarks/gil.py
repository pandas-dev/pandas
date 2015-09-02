from .pandas_vb_common import *
from pandas.core import common as com
try:
    from pandas.util.testing import test_parallel
    have_real_test_parallel = True
except ImportError:
    have_real_test_parallel = False

    def test_parallel(num_threads=1):

        def wrapper(fname):
            return fname
        return wrapper


class nogil_groupby_count_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_count_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].count()


class nogil_groupby_last_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_last_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].last()


class nogil_groupby_max_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_max_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].max()


class nogil_groupby_mean_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_mean_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].mean()


class nogil_groupby_min_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_min_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].min()


class nogil_groupby_prod_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_prod_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].prod()


class nogil_groupby_sum_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_sum_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].sum()


class nogil_groupby_sum_4(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_sum_4(self):
        self.pg4()

    def f(self):
        self.df.groupby('key')['data'].sum()

    def g2(self):
        for i in range(2):
            self.f()

    def g4(self):
        for i in range(4):
            self.f()

    def g8(self):
        for i in range(8):
            self.f()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.f()

    @test_parallel(num_threads=4)
    def pg4(self):
        self.f()

    @test_parallel(num_threads=8)
    def pg8(self):
        self.f()


class nogil_groupby_sum_8(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_sum_8(self):
        self.pg8()

    def f(self):
        self.df.groupby('key')['data'].sum()

    def g2(self):
        for i in range(2):
            self.f()

    def g4(self):
        for i in range(4):
            self.f()

    def g8(self):
        for i in range(8):
            self.f()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.f()

    @test_parallel(num_threads=4)
    def pg4(self):
        self.f()

    @test_parallel(num_threads=8)
    def pg8(self):
        self.f()


class nogil_groupby_var_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError

    def time_nogil_groupby_var_2(self):
        self.pg2()

    @test_parallel(num_threads=2)
    def pg2(self):
        self.df.groupby('key')['data'].var()


class nogil_take1d_float64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError
        self.N = 10000000.0
        self.df = DataFrame({'int64': np.arange(self.N, dtype='int64'), 'float64': np.arange(self.N, dtype='float64'), })
        self.indexer = np.arange(100, (len(self.df) - 100))

    def time_nogil_take1d_float64(self):
        self.take_1d_pg2_int64()

    @test_parallel(num_threads=2)
    def take_1d_pg2_int64(self):
        com.take_1d(self.df.int64.values, self.indexer)

    @test_parallel(num_threads=2)
    def take_1d_pg2_float64(self):
        com.take_1d(self.df.float64.values, self.indexer)


class nogil_take1d_int64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        if (not have_real_test_parallel):
            raise NotImplementedError
        self.N = 10000000.0
        self.df = DataFrame({'int64': np.arange(self.N, dtype='int64'), 'float64': np.arange(self.N, dtype='float64'), })
        self.indexer = np.arange(100, (len(self.df) - 100))

    def time_nogil_take1d_int64(self):
        self.take_1d_pg2_float64()

    @test_parallel(num_threads=2)
    def take_1d_pg2_int64(self):
        com.take_1d(self.df.int64.values, self.indexer)

    @test_parallel(num_threads=2)
    def take_1d_pg2_float64(self):
        com.take_1d(self.df.float64.values, self.indexer)


class nogil_kth_smallest(object):
    number = 1
    repeat = 5

    def setup(self):
        if (not have_real_test_parallel):
            raise NotImplementedError
        np.random.seed(1234)
        self.N = 10000000
        self.k = 500000
        self.a = np.random.randn(self.N)
        self.b = self.a.copy()
        self.kwargs_list = [{'arr': self.a}, {'arr': self.b}]

    def time_nogil_kth_smallest(self):
        @test_parallel(num_threads=2, kwargs_list=self.kwargs_list)
        def run(arr):
            algos.kth_smallest(arr, self.k)
        run()
