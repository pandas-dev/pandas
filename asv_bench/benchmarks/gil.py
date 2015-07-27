from pandas_vb_common import *
from pandas.core import common as com
from pandas.util.testing import test_parallel


class nogil_groupby_count_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].count()

    def time_nogil_groupby_count_2(self):
        pg2()


class nogil_groupby_last_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].last()

    def time_nogil_groupby_last_2(self):
        pg2()


class nogil_groupby_max_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].max()

    def time_nogil_groupby_max_2(self):
        pg2()


class nogil_groupby_mean_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].mean()

    def time_nogil_groupby_mean_2(self):
        pg2()


class nogil_groupby_min_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].min()

    def time_nogil_groupby_min_2(self):
        pg2()


class nogil_groupby_prod_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].prod()

    def time_nogil_groupby_prod_2(self):
        pg2()


class nogil_groupby_sum_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].sum()

    def time_nogil_groupby_sum_2(self):
        pg2()


class nogil_groupby_sum_4(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        def f():
            self.df.groupby('key')['data'].sum()

        def g2():
            for i in range(2):
                f()

        def g4():
            for i in range(4):
                f()

        def g8():
            for i in range(8):
                f()

        @test_parallel(num_threads=2)
        def pg2():
            f()

        @test_parallel(num_threads=4)
        def pg4():
            f()

        @test_parallel(num_threads=8)
        def pg8():
            f()

    def time_nogil_groupby_sum_4(self):
        pg4()


class nogil_groupby_sum_8(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        def f():
            self.df.groupby('key')['data'].sum()

        def g2():
            for i in range(2):
                f()

        def g4():
            for i in range(4):
                f()

        def g8():
            for i in range(8):
                f()

        @test_parallel(num_threads=2)
        def pg2():
            f()

        @test_parallel(num_threads=4)
        def pg4():
            f()

        @test_parallel(num_threads=8)
        def pg8():
            f()

    def time_nogil_groupby_sum_8(self):
        pg8()


class nogil_groupby_var_2(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })

        @test_parallel(num_threads=2)
        def pg2():
            self.df.groupby('key')['data'].var()

    def time_nogil_groupby_var_2(self):
        pg2()


class nogil_take1d_float64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        self.N = 10000000.0
        self.df = DataFrame({'int64': np.arange(self.N, dtype='int64'), 'float64': np.arange(self.N, dtype='float64'), })
        self.indexer = np.arange(100, (len(self.df) - 100))

        @test_parallel(num_threads=2)
        def take_1d_pg2_int64():
            com.take_1d(self.df.int64.values, self.indexer)

        @test_parallel(num_threads=2)
        def take_1d_pg2_float64():
            com.take_1d(self.df.float64.values, self.indexer)

    def time_nogil_take1d_float64(self):
        take_1d_pg2_int64()


class nogil_take1d_int64(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.ngroups = 1000
        np.random.seed(1234)
        self.df = DataFrame({'key': np.random.randint(0, self.ngroups, size=self.N), 'data': np.random.randn(self.N), })
        self.N = 10000000.0
        self.df = DataFrame({'int64': np.arange(self.N, dtype='int64'), 'float64': np.arange(self.N, dtype='float64'), })
        self.indexer = np.arange(100, (len(self.df) - 100))

        @test_parallel(num_threads=2)
        def take_1d_pg2_int64():
            com.take_1d(self.df.int64.values, self.indexer)

        @test_parallel(num_threads=2)
        def take_1d_pg2_float64():
            com.take_1d(self.df.float64.values, self.indexer)

    def time_nogil_take1d_int64(self):
        take_1d_pg2_float64()