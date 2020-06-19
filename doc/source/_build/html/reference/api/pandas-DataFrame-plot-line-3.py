# An example with subplots, so an array of axes is returned.

axes = df.plot.line(subplots=True)
type(axes)
# <class 'numpy.ndarray'>
