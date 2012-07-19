import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

class RPlot:
	def __init__(self, data, x=None, y=None):
		self.data = data
		self.ax = plt.gca()
		self.aes = {}
		if x is not None and y is None:
			self.aes['x'] = x
		elif x is not None and y is not None:
			self.aes['x'] = x
			self.aes['y'] = y

	def __add__(self, other):
		other.plot(self)

class GeomPoint:
	def __init__(self, x=None, y=None, marker='o', colour='grey', size=20, alpha=1.0):
		self.x = x
		self.y = y
		self.marker = marker
		self.colour = colour
		self.size = size
		self.alpha = alpha

	def plot(self, rplot):
		aes = rplot.aes
		if self.x is not None:
			aes['x'] = self.x
		if self.y is not None:
			aes['y'] = self.y
		rplot.ax.scatter(rplot.data[aes['x']], rplot.data[aes['y']], c=self.colour, marker=self.marker, s=self.size, alpha=self.alpha)

class GeomDensity2d:
	def __init__(self, x=None, y=None, weight=1.0, colour='grey', size=0.5, linetype=1.0, alpha=1.0):
		self.x = x
		self.y = y
		self.weight = weight
		self.colour = colour
		self.size = size
		self.linetype = linetype
		self.alpha = alpha

	def plot(self, rplot):
		aes = rplot.aes
		if self.x is not None:
			aes['x'] = self.x
		if self.y is not None:
			aes['y'] = self.y
		x = rplot.data[aes['x']]
		y = rplot.data[aes['y']]
		rvs = np.array([x, y])
		x_min = x.min()
		x_max = x.max()
		y_min = y.min()
		y_max = y.max()
		X, Y = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
		positions = np.vstack([X.ravel(), Y.ravel()])
		values = np.vstack([x, y])
		kernel = stats.gaussian_kde(values)
		Z = np.reshape(kernel(positions).T, X.shape)
		rplot.ax.contour(Z, alpha=self.alpha, extent=[x_min, x_max, y_min, y_max])
