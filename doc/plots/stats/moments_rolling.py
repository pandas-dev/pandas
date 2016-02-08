from moment_plots import *

ts = test_series()
s = ts.cumsum()

s[20:50] = np.NaN
s[120:150] = np.NaN
plot_timeseries(s,
                m.rolling_count(s, 50),
                m.rolling_sum(s, 50, min_periods=10),
                m.rolling_mean(s, 50, min_periods=10),
                m.rolling_std(s, 50, min_periods=10),
                m.rolling_skew(s, 50, min_periods=10),
                m.rolling_kurt(s, 50, min_periods=10),
                size=(10, 12),
                titles=('time series',
                        'rolling_count',
                        'rolling_sum',
                        'rolling_mean',
                        'rolling_std',
                        'rolling_skew',
                        'rolling_kurt'))
plt.show()
plt.close('all')
