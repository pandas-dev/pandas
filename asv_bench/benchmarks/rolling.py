from .pandas_vb_common import *
import pandas as pd
import numpy as np


class DataframeRolling(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.Ns = 10000
        self.df = pd.DataFrame({'a': np.random.random(self.N)})
        self.dfs = pd.DataFrame({'a': np.random.random(self.Ns)})
        self.wins = 10
        self.winl = 1000

    def time_rolling_quantile_0(self):
        (self.df.rolling(self.wins).quantile(0.0))

    def time_rolling_quantile_1(self):
        (self.df.rolling(self.wins).quantile(1.0))

    def time_rolling_quantile_median(self):
        (self.df.rolling(self.wins).quantile(0.5))

    def time_rolling_median(self):
        (self.df.rolling(self.wins).median())

    def time_rolling_mean(self):
        (self.df.rolling(self.wins).mean())

    def time_rolling_max(self):
        (self.df.rolling(self.wins).max())

    def time_rolling_min(self):
        (self.df.rolling(self.wins).min())

    def time_rolling_std(self):
        (self.df.rolling(self.wins).std())

    def time_rolling_count(self):
        (self.df.rolling(self.wins).count())

    def time_rolling_skew(self):
        (self.df.rolling(self.wins).skew())

    def time_rolling_kurt(self):
        (self.df.rolling(self.wins).kurt())

    def time_rolling_sum(self):
        (self.df.rolling(self.wins).sum())

    def time_rolling_corr(self):
        (self.dfs.rolling(self.wins).corr())

    def time_rolling_cov(self):
        (self.dfs.rolling(self.wins).cov())
        
    def time_rolling_quantile_0_l(self):
        (self.df.rolling(self.winl).quantile(0.0))

    def time_rolling_quantile_1_l(self):
        (self.df.rolling(self.winl).quantile(1.0))

    def time_rolling_quantile_median_l(self):
        (self.df.rolling(self.winl).quantile(0.5))

    def time_rolling_median_l(self):
        (self.df.rolling(self.winl).median())

    def time_rolling_mean_l(self):
        (self.df.rolling(self.winl).mean())

    def time_rolling_max_l(self):
        (self.df.rolling(self.winl).max())

    def time_rolling_min_l(self):
        (self.df.rolling(self.winl).min())

    def time_rolling_std_l(self):
        (self.df.rolling(self.wins).std())

    def time_rolling_count_l(self):
        (self.df.rolling(self.wins).count())

    def time_rolling_skew_l(self):
        (self.df.rolling(self.wins).skew())

    def time_rolling_kurt_l(self):
        (self.df.rolling(self.wins).kurt())

    def time_rolling_sum_l(self):
        (self.df.rolling(self.wins).sum())


class SeriesRolling(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.Ns = 10000
        self.df = pd.DataFrame({'a': np.random.random(self.N)})
        self.dfs = pd.DataFrame({'a': np.random.random(self.Ns)})
        self.sr = self.df.a
        self.srs = self.dfs.a
        self.wins = 10
        self.winl = 1000

    def time_rolling_quantile_0(self):
        (self.sr.rolling(self.wins).quantile(0.0))

    def time_rolling_quantile_1(self):
        (self.sr.rolling(self.wins).quantile(1.0))

    def time_rolling_quantile_median(self):
        (self.sr.rolling(self.wins).quantile(0.5))

    def time_rolling_median(self):
        (self.sr.rolling(self.wins).median())

    def time_rolling_mean(self):
        (self.sr.rolling(self.wins).mean())

    def time_rolling_max(self):
        (self.sr.rolling(self.wins).max())

    def time_rolling_min(self):
        (self.sr.rolling(self.wins).min())

    def time_rolling_std(self):
        (self.sr.rolling(self.wins).std())

    def time_rolling_count(self):
        (self.sr.rolling(self.wins).count())

    def time_rolling_skew(self):
        (self.sr.rolling(self.wins).skew())

    def time_rolling_kurt(self):
        (self.sr.rolling(self.wins).kurt())

    def time_rolling_sum(self):
        (self.sr.rolling(self.wins).sum())

    def time_rolling_corr(self):
        (self.srs.rolling(self.wins).corr())

    def time_rolling_cov(self):
        (self.srs.rolling(self.wins).cov())
        
    def time_rolling_quantile_0_l(self):
        (self.sr.rolling(self.winl).quantile(0.0))

    def time_rolling_quantile_1_l(self):
        (self.sr.rolling(self.winl).quantile(1.0))

    def time_rolling_quantile_median_l(self):
        (self.sr.rolling(self.winl).quantile(0.5))

    def time_rolling_median_l(self):
        (self.sr.rolling(self.winl).median())

    def time_rolling_mean_l(self):
        (self.sr.rolling(self.winl).mean())

    def time_rolling_max_l(self):
        (self.sr.rolling(self.winl).max())

    def time_rolling_min_l(self):
        (self.sr.rolling(self.winl).min())

    def time_rolling_std_l(self):
        (self.sr.rolling(self.wins).std())

    def time_rolling_count_l(self):
        (self.sr.rolling(self.wins).count())

    def time_rolling_skew_l(self):
        (self.sr.rolling(self.wins).skew())

    def time_rolling_kurt_l(self):
        (self.sr.rolling(self.wins).kurt())

    def time_rolling_sum_l(self):
        (self.sr.rolling(self.wins).sum())
