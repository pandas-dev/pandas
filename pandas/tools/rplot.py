import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import random
import pdb

def random_colour(name):
	"""Random colour from a string or other hashable value.

	Parameters:
	-----------
	name: A string value to use as a seed to the random number generator.

	Returns:
	--------
	(r, g, b): Where r, g, b are random numbers in the range (0, 1).
	"""
	random.seed(name)
	return [random.random() for _ in range(3)]

def filter_column(frame, column, filter_column, filter_value):
	"""Select only those values from column that have a specified value in another column.

	Parameters:
	-----------
	frame: pandas data frame object.
	column: column name from which to select values.
	filter_column: column used to filter.
	filter_value: only those rows with this value in filter_column will be selected.

	Returns:
	--------
	numpy array with filtered values.
	"""
	n = len(frame)
	vcol = frame[column]
	fcol = frame[filter_column]
	result = []
	for v, f in zip(vcol, fcol):
		if f == filter_value:
			result.append(v)
	return np.array(result)

def parse_facets(facet):
	"""Parse facets formula of the form 'lhs ~ rhs'.

	Parameters:
	-----------
	facets: facets formula.

	Returns:
	--------
	A list with LHS and RHS column names.
	"""
	lhs, rhs = [col.strip() for col in facet.split('~')]
	return (lhs, rhs)

def scale_size(column, categorical, min_size=1.0, max_size=50.0):
	def scaler(data, index):
		if categorical:
			pass
		else:
			x = data[column].iget(index)
			a = min(data[column])
			b = max(data[column])
			return min_size + ((x - a) / (b - a)) * (max_size - min_size)
	return scaler

class Layer:
	"""
	Layer object representing a single plot layer.
	"""
	def __init__(self, geom, geom_params, data, aesthetics):
		self.geom = geom
		self.geom_params = geom_params
		self.data = data
		self.aesthetics = aesthetics

	def render(self, ax):
		if self.geom == 'point':
			for index in range(len(self.data)):
				row = self.data.irow(index)
				x = row[self.aesthetics['x']]
				y = row[self.aesthetics['y']]
				size_scaler = self.aesthetics['size'](self.data, index)
				ax.scatter(x, y, size_scaler(self.data, index))

class RPlot:
	def __init__(self, data, x=None, y=None):
		seomlf.data = data
		self.ax = plt.gca()
		self.aes = {
			'x' : None,
			'y' : None,
		}
		if x is not None and y is None:
			self.aes['x'] = x
		elif x is not None and y is not None:
			self.aes['x'] = x
			self.aes['y'] = y

	def __add__(self, other):
		return other.plot(self)

class GeomPoint:
	def __init__(self, x=None, y=None, shape='o', colour='grey', size=20, alpha=1.0):
		self.x = x
		self.y = y
		self.shape = shape
		self.colour = colour
		self.size = size
		self.alpha = alpha

	def plot(self, rplot):
		aes = rplot.aes
		if self.x is not None:
			aes['x'] = self.x
		if self.y is not None:
			aes['y'] = self.y
		if type(self.colour) is not type(""):
			colours = list(set(self.colour))
			for colour in colours:
				xcol = filter_column(rplot.data, aes['x'], self.colour.name, colour)
				ycol = filter_column(rplot.data, aes['y'], self.colour.name, colour)
				rplot.ax.scatter(xcol, ycol, c=random_colour(colour), marker=self.shape, s=self.size, alpha=self.alpha, label=colour)
		else:
			rplot.ax.scatter(rplot.data[aes['x']], rplot.data[aes['y']], c=self.colour, marker=self.shape, s=self.size, alpha=self.alpha)
		rplot.ax.set_xlabel(aes['x'])
		rplot.ax.set_ylabel(aes['y'])
		rplot.ax.legend()
		return rplot

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
		return rplot