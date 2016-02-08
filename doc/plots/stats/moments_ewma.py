import matplotlib.pyplot as plt
import pandas.util.testing as t
import pandas.stats.moments as m

t.N = 200
s = t.makeTimeSeries().cumsum()

plt.figure(figsize=(10, 5))
plt.plot(s.index, s.values)
plt.plot(s.index, m.ewma(s, 20, min_periods=1).values)
f = plt.gcf()
f.autofmt_xdate()

plt.show()
plt.close('all')
