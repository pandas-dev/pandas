np.random.seed(5)
x = np.cumsum(np.random.normal(loc=1, scale=5, size=50))
s = pd.Series(x)
s.plot()
