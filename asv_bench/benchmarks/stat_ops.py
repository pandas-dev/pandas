from .pandas_vb_common import *


class stat_ops_frame_mean_float_axis_0(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_mean_float_axis_0(self):
        self.df.mean()


class stat_ops_frame_mean_float_axis_1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_mean_float_axis_1(self):
        self.df.mean(1)


class stat_ops_frame_mean_int_axis_0(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_mean_int_axis_0(self):
        self.dfi.mean()


class stat_ops_frame_mean_int_axis_1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_mean_int_axis_1(self):
        self.dfi.mean(1)


class stat_ops_frame_sum_float_axis_0(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_sum_float_axis_0(self):
        self.df.sum()


class stat_ops_frame_sum_float_axis_1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_sum_float_axis_1(self):
        self.df.sum(1)


class stat_ops_frame_sum_int_axis_0(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_sum_int_axis_0(self):
        self.dfi.sum()


class stat_ops_frame_sum_int_axis_1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 4))
        self.dfi = DataFrame(np.random.randint(1000, size=self.df.shape))

    def time_stat_ops_frame_sum_int_axis_1(self):
        self.dfi.sum(1)


class stat_ops_level_frame_sum(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4), index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4), index=self.index.levels[1])

    def time_stat_ops_level_frame_sum(self):
        self.df.sum(level=1)


class stat_ops_level_frame_sum_multiple(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4), index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4), index=self.index.levels[1])

    def time_stat_ops_level_frame_sum_multiple(self):
        self.df.sum(level=[0, 1])


class stat_ops_level_series_sum(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4), index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4), index=self.index.levels[1])

    def time_stat_ops_level_series_sum(self):
        self.df[1].sum(level=1)


class stat_ops_level_series_sum_multiple(object):
    goal_time = 0.2

    def setup(self):
        self.index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)], labels=[np.arange(10).repeat(10000), np.tile(np.arange(100).repeat(100), 10), np.tile(np.tile(np.arange(100), 100), 10)])
        random.shuffle(self.index.values)
        self.df = DataFrame(np.random.randn(len(self.index), 4), index=self.index)
        self.df_level = DataFrame(np.random.randn(100, 4), index=self.index.levels[1])

    def time_stat_ops_level_series_sum_multiple(self):
        self.df[1].sum(level=[0, 1])


class stat_ops_series_std(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randn(100000), index=np.arange(100000))
        self.s[::2] = np.nan

    def time_stat_ops_series_std(self):
        self.s.std()


class stats_corr_spearman(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 30))

    def time_stats_corr_spearman(self):
        self.df.corr(method='spearman')


class stats_rank2d_axis0_average(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(5000, 50))

    def time_stats_rank2d_axis0_average(self):
        self.df.rank()


class stats_rank2d_axis1_average(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(5000, 50))

    def time_stats_rank2d_axis1_average(self):
        self.df.rank(1)


class stats_rank_average(object):
    goal_time = 0.2

    def setup(self):
        self.values = np.concatenate([np.arange(100000), np.random.randn(100000), np.arange(100000)])
        self.s = Series(self.values)

    def time_stats_rank_average(self):
        self.s.rank()


class stats_rank_average_int(object):
    goal_time = 0.2

    def setup(self):
        self.values = np.random.randint(0, 100000, size=200000)
        self.s = Series(self.values)

    def time_stats_rank_average_int(self):
        self.s.rank()


class stats_rank_pct_average(object):
    goal_time = 0.2

    def setup(self):
        self.values = np.concatenate([np.arange(100000), np.random.randn(100000), np.arange(100000)])
        self.s = Series(self.values)

    def time_stats_rank_pct_average(self):
        self.s.rank(pct=True)


class stats_rank_pct_average_old(object):
    goal_time = 0.2

    def setup(self):
        self.values = np.concatenate([np.arange(100000), np.random.randn(100000), np.arange(100000)])
        self.s = Series(self.values)

    def time_stats_rank_pct_average_old(self):
        (self.s.rank() / len(self.s))


class stats_rolling_mean(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randn(100000)
        self.win = 100

    def time_rolling_mean(self):
        rolling_mean(self.arr, self.win)

    def time_rolling_median(self):
        rolling_median(self.arr, self.win)

    def time_rolling_min(self):
        rolling_min(self.arr, self.win)

    def time_rolling_max(self):
        rolling_max(self.arr, self.win)

    def time_rolling_sum(self):
        rolling_sum(self.arr, self.win)

    def time_rolling_std(self):
        rolling_std(self.arr, self.win)

    def time_rolling_var(self):
        rolling_var(self.arr, self.win)

    def time_rolling_skew(self):
        rolling_skew(self.arr, self.win)

    def time_rolling_kurt(self):
        rolling_kurt(self.arr, self.win)
