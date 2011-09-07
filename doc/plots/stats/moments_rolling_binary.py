from moment_plots import *

np.random.seed(1)

ts = test_series()
s = ts.cumsum()
ts2 = test_series()
s2 = ts2.cumsum()

s[20:50] = np.NaN
s[120:150] = np.NaN
fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

ax0, ax1, ax2 = axes

ax0.plot(s.index, s.values)
ax0.plot(s2.index, s2.values)
ax0.set_title('time series')

ax1.plot(s.index, m.rolling_corr(s, s2, 50, min_periods=1).values)
ax1.set_title('rolling_corr')

ax2.plot(s.index, m.rolling_cov(s, s2, 50, min_periods=1).values)
ax2.set_title('rolling_cov')

fig.autofmt_xdate()
fig.subplots_adjust(bottom=0.10, top=0.95)

plt.show()
plt.close('all')
