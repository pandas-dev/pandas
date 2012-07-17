import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

class RPlot:
	def __init__(self, data):
		self.data = data
		self.ax = plt.gca()
		self.aes = {}

	def __add__(self, other):
		other.plot(self.ax, self.data)

class GeomPoint:
	def __init__(self, x, y, marker='o', colour='grey', size=20, alpha=1.0):
		self.x = x
		self.y = y
		self.marker = marker
		self.colour = colour
		self.size = size
		self.alpha = alpha

	def plot(self, ax, data):
		ax.scatter(data[self.x], data[self.y], c=self.colour, marker=self.marker, s=self.size, alpha=self.alpha)


class GeomDensity2d:
	def __init__(self, x, y, weight=1.0, colour='grey', size=0.5, linetype=1.0, alpha=1.0):
		self.x = x
		self.y = y
		self.weight = weight
		self.colour = colour
		self.size = size
		self.linetype = linetype
		self.alpha = alpha

	def plot(self, ax, data):
		x = data[self.x]
		y = data[self.y]
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
		ax.contour(Z, alpha=alpha, extent=[x_min, x_max, y_min, y_max])