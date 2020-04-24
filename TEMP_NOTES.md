In [20]: s = pd.Series(np.random.randn(100000))

In [21]: n = 10

In [22]: q = np.linspace(1./n, 1 - 1./n, n - 1)

In [23]: r = s.rolling(window=100)

In [24]: %time a = pd.DataFrame({qq: r.quantile(qq) for qq in q})
CPU times: user 333 ms, sys: 13.9 ms, total: 347 ms
Wall time: 347 ms

In [25]: %time b = r.quantile(q)
CPU times: user 97.9 ms, sys: 7.96 ms, total: 106 ms
Wall time: 106 ms

In [26]: r = s.rolling(window=1000)

In [27]: %time a = pd.DataFrame({qq: r.quantile(qq) for qq in q})
CPU times: user 477 ms, sys: 9.32 ms, total: 486 ms
Wall time: 486 ms

In [28]: %time b = r.quantile(q)
CPU times: user 148 ms, sys: 0 ns, total: 148 ms
Wall time: 156 ms

In [29]: n = 100

In [30]: q = np.linspace(1./n, 1 - 1./n, n - 1)

In [31]: %time a = pd.DataFrame({qq: r.quantile(qq) for qq in q})
CPU times: user 5.33 s, sys: 41.7 ms, total: 5.38 s
Wall time: 5.38 s

In [32]: %time b = r.quantile(q)
CPU times: user 1.25 s, sys: 25.6 ms, total: 1.28 s
Wall time: 1.28 s
