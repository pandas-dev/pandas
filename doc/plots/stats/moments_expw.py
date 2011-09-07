from moment_plots import *

np.random.seed(1)

ts = test_series(500) * 10

# ts[::100] = 20

s = ts.cumsum()

fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

ax0, ax1, ax2 = axes

ax0.plot(s.index, s.values)
ax0.set_title('time series')

ax1.plot(s.index, m.ewma(s, span=50, min_periods=1).values, color='b')
ax1.plot(s.index, m.rolling_mean(s, 50, min_periods=1).values, color='r')
ax1.set_title('rolling_mean vs. ewma')

line1 = ax2.plot(s.index, m.ewmstd(s, span=50, min_periods=1).values, color='b')
line2 = ax2.plot(s.index, m.rolling_std(s, 50, min_periods=1).values, color='r')
ax2.set_title('rolling_std vs. ewmstd')

fig.legend((line1, line2),
           ('Exp-weighted', 'Equal-weighted'),
           loc='upper right')
fig.autofmt_xdate()
fig.subplots_adjust(bottom=0.10, top=0.95)

plt.show()
plt.close('all')
