spacing = np.linspace(-9 * np.pi, 9 * np.pi, num=1000)
s = pd.Series(0.7 * np.random.rand(1000) + 0.3 * np.sin(spacing))
pd.plotting.autocorrelation_plot(s)
