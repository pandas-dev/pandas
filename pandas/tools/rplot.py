import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import random
import pdb
from copy import deepcopy

#
# TODO:
# * Make sure trellis display works when two or one grouping variable is specified
# * Enable labelling for legend display
# * Expand RPlot class
#

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

def scale_size(column, categorical, min_size=1.0, max_size=80.0):
	"""Creates a function that converts between a data attribute to point size.

	Parameters:
	-----------
	column: string, a column name
	categorical: boolean, true if the column contains categorical data
	min_size: float, minimum point size
	max_size: float, maximum point size

	Returns:
	--------
	a function of two arguments that takes a data set and a row number, returns float
	"""
	def scaler(data, index):
		if categorical:
			pass
		else:
			x = data[column].iget(index)
			a = min(data[column])
			b = max(data[column])
			return min_size + ((x - a) / (b - a)) * (max_size - min_size)
	return scaler

def scale_gradient(column, categorical, colour1=(0.0, 0.0, 0.0), colour2=(1.0, 0.7, 0.8)):
	"""Create a function that converts between a data attribute value to a 
	point in colour space between two specified colours.

	Parameters:
	-----------
	column: string, a column name
	categorical: boolean, true if the column contains categorical data
	colour1: a tuple with three float values specifying rgb components
	colour2: a tuple with three float values specifying rgb components

	Returns:
	--------
	a function of two arguments that takes a data set and a row number, returns a
	tuple with three float values with rgb component values.
	"""
	def scaler(data, index):
		if categorical:
			pass
		else:
			x = data[column].iget(index)
			a = min(data[column])
			b = max(data[column])
			r1, g1, b1 = colour1
			r2, g2, b2 = colour2
			x_scaled = (x - a) / (b - a)
			return (r1 + (r2 - r1) * x_scaled,
					g1 + (g2 - g1) * x_scaled,
					b1 + (b2 - b1) * x_scaled)
	return scaler

def scale_gradient2(column, categorical, colour1=(0.0, 0.0, 0.0), colour2=(1.0, 0.7, 0.8), colour3=(0.2, 1.0, 0.5)):
	"""Create a function that converts between a data attribute value to a 
	point in colour space between three specified colours.

	Parameters:
	-----------
	column: string, a column name
	categorical: boolean, true if the column contains categorical data
	colour1: a tuple with three float values specifying rgb components
	colour2: a tuple with three float values specifying rgb components
	colour3: a tuple with three float values specifying rgb components

	Returns:
	--------
	a function of two arguments that takes a data set and a row number, returns a
	tuple with three float values with rgb component values.
	"""
	def scaler(data, index):
		if categorical:
			pass
		else:
			x = data[column].iget(index)
			a = min(data[column])
			b = max(data[column])
			r1, g1, b1 = colour1
			r2, g2, b2 = colour2
			r3, g3, b3 = colour3
			x_scaled = (x - a) / (b - a)
			if x_scaled < 0.5:
				x_scaled *= 2.0
				return (r1 + (r2 - r1) * x_scaled,
						g1 + (g2 - g1) * x_scaled,
						b1 + (b2 - b1) * x_scaled)
			else:
				x_scaled = (x_scaled - 0.5) * 2.0
				return (r2 + (r3 - r2) * x_scaled,
						g2 + (g3 - g2) * x_scaled,
						b2 + (b3 - b2) * x_scaled)
	return scaler

def scale_shape(column):
	"""Create a function that converts between a categorical value and a scatter plot shape.

	Parameters:
	-----------
	column: string, a column name to use

	Returns:
	--------
	a function of two arguments
	"""
	shapes = ['o', 'D', 'h', 'H', '_', '8', 'p', '+', '.', 's', '*', 'd', '^', '<', '>', 'v', '|', 'x']
	def scaler(data, index):
		values = list(set(data[column]))
		x = data[column].iget(index)
		return shapes[values.index(x)]
	return scaler

def scale_constant(constant):
	"""Create a function that always returns a specified constant value.

	Parameters:
	-----------
	constant: a Python object to be returned

	Returns:
	--------
	a two argument function
	"""
	def scaler(data, index):
		return constant
	return scaler

class Layer:
	"""
	Layer object representing a single plot layer.
	"""
	def __init__(self, data, aesthetics):
		self.data = data
		self.aesthetics = aesthetics

class GeomPoint(Layer):
	def work(self, ax=None, fig=None):
		"""Render the layer on a matplotlib axis.

		Parameters:
		-----------
		ax: matplotlib axis object to draw on

		Returns:
		--------
		ax: matplotlib axis object
		"""
		for index in range(len(self.data)):
			row = self.data.irow(index)
			x = row[self.aesthetics['x']]
			y = row[self.aesthetics['y']]
			size_scaler = self.aesthetics['size']
			colour_scaler = self.aesthetics['colour']
			shape_scaler = self.aesthetics['shape']
			alpha = self.aesthetics['alpha']
			ax.scatter(x, y, 
				s=size_scaler(self.data, index),
				c=colour_scaler(self.data, index),
				marker=shape_scaler(self.data, index),
				alpha=alpha)
		ax.set_xlabel(self.aesthetics['x'])
		ax.set_ylabel(self.aesthetics['y'])
	return ax

def display_grouped(grouped_data, x, y, fig):
	"""A test routine to display grouped data.

	Parameters:
	-----------
	grouped_data: data frame grouped by df.groupby pandas routine
	fig: matplotlib figure

	Returns:
	--------
	Nothing
	"""
	shingle1 = set([])
	shingle2 = set([])
	# Fill shingles.
	for name, group in grouped_data:
		if type(name) is type(()):
			shingle1.add(name[0])
			shingle2.add(name[1])
		else:
			shingle1.add(name)
	rows = len(shingle1)
	cols = len(shingle2)
	subplot_nr = 1
	for name, group, in grouped_data:
		ax = fig.add_subplot(rows, cols, subplot_nr)
		ax.scatter(group[x], group[y])
		subplot_nr += 1

class TrellisGrid:
	def __init__(self, layer, by):
		"""Initialize TreelisGrid instance.

		Parameters:
		-----------
		layer: a Layer object instance
		by: column names to group by
		"""
		self.data = layer.data
		self.grouped = self.data.groupby(by)
		self.groups = self.grouped.groups.keys()
		self.by = by
		self.shingle1 = set([g[0] for g in self.groups])
		self.shingle2 = set([g[1] for g in self.groups])
		self.rows = len(self.shingle1)
		self.cols = len(self.shingle2)
		self.grid = [[None for _ in range(self.cols)] for _ in range(self.rows)]
		self.group_grid = [[None for _ in range(self.cols)] for _ in range(self.rows)]
		row = 0
		col = 0
		for group, data in self.grouped:
			new_layer = deepcopy(layer)
			new_layer.data = data
			self.grid[row][col] = new_layer
			self.group_grid[row][col] = group
			col += 1
			if col >= self.cols:
				col = 0
				row += 1

	def get_layer(self, row, col):
		"""Get a layer associated with the specified row and col.

		Parameters:
		-----------
		row: integer row index
		col: integer column index

		Returns:
		--------
		layer object
		"""
		return self.grid[row][col]

	def get_group(self, row, col):
		"""Get a group tuple for the specified row and col.

		Parameters:
		-----------
		row: integer row index
		col: integer column index 

		Returns:
		--------
		a tuple
		"""
		return self.group_grid[row][col]

	def work(self, ax=None, fig=None):
		"""Render the trellis plot on a figure.

		Parameters:
		-----------
		fig: matplotlib figure to draw on

		Returns:
		--------
		matplotlib figure
		"""
		index = 1
		axes = []
		for row in self.grid:
			for layer in row:
				ax = fig.add_subplot(self.rows, self.cols, index)
				layer.render(ax)
				axes.append(ax)
				index += 1
		min_x = min([ax.get_xlim()[0] for ax in axes])
		max_x = max([ax.get_xlim()[1] for ax in axes])
		min_y = min([ax.get_ylim()[0] for ax in axes])
		max_y = max([ax.get_ylim()[1] for ax in axes])
		[ax.set_xlim(min_x, max_x) for ax in axes]
		[ax.set_ylim(min_y, max_y) for ax in axes]
		for index, axis in enumerate(axes):
			if index % self.cols == 0:
				pass
			else:
				axis.get_yaxis().set_ticks([])
				axis.set_ylabel('')
			if index / self.cols == self.rows - 1:
				pass
			else:
				axis.get_xaxis().set_ticks([])
				axis.set_xlabel('')
			label1 = "%s = %s" % (self.by[0], self.group_grid[index / self.cols][index % self.cols][0])
			label2 = "%s = %s" % (self.by[1], self.group_grid[index / self.cols][index % self.cols][1])
			if self.cols > 1:
				axis.table(cellText=[[label1], [label2]], 
					loc='top', cellLoc='center', 
					cellColours=[['lightgrey'], ['lightgrey']])
			else:
				axis.table(cellText=[[label1]], loc='top', cellLoc='center', cellColours=[['lightgrey']])
		fig.subplots_adjust(wspace=0.05, hspace=0.2)
		return fig

class RPlot:
	def __init__(self):
		self.layers = []

	def __add__(self, other):
		self.layers.append(other)

	def show(self, fig):
		for layer in self.layers:
			pass

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