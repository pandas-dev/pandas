import matplotlib.pyplot as plt
import pandas.util.testing as t
import pandas.stats.moments as m

t.N = 500
ts = t.makeTimeSeries()
ts[::100] = 20

s = ts.cumsum()


plt.figure(figsize=(10, 5))
plt.plot(s.index, m.ewmvol(s, span=50, min_periods=1).values, color='b')
plt.plot(s.index, m.rolling_std(s, 50, min_periods=1).values, color='r')

plt.title('Exp-weighted std with shocks')
plt.legend(('Exp-weighted', 'Equal-weighted'))

f = plt.gcf()
f.autofmt_xdate()

plt.show()
plt.close('all')
