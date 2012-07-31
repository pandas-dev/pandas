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

def default_aes(x=None, y=None):
	"""Create the default aesthetics dictionary.
	"""
	return {
		'x' : x,
		'y' : y,
		'size' : scale_constant(40.0),
		'colour' : scale_constant('grey'),
		'shape' : scale_constant('o'),
		'alpha' : scale_constant(1.0),
	}

def aes(x=None, y=None, size=None, colour=None, shape=None, alpha=None):
	"""Create an empty aesthetics dictionary.
	"""
	return {
		'x' : x,
		'y' : y,
		'size' : size,
		'colour' : colour,
		'shape' : shape,
		'alpha' : alpha,
	}

class Layer:
	"""
	Layer object representing a single plot layer.
	"""
	def __init__(self, data=None, aes_=None):
		self.data = data
		if aes_ is None:
			self.aes = aes()
		else:
			self.aes = aes_

	def work(self, fig=None, ax=None):
		pass

class GeomPoint(Layer):
	def work(self, fig=None, ax=None):
		"""Render the layer on a matplotlib axis.
		You can specify either a figure or an axis to draw on.

		Parameters:
		-----------
		fig: matplotlib figure object
		ax: matplotlib axis object to draw on

		Returns:
		--------
		fig, ax: matplotlib figure and axis objects
		"""
		if ax is None:
			if fig is None:
				return fig, ax
			else:
				ax = fig.gca()
		for index in range(len(self.data)):
			row = self.data.irow(index)
			x = row[self.aes['x']]
			y = row[self.aes['y']]
			size_scaler = self.aes['size']
			colour_scaler = self.aes['colour']
			shape_scaler = self.aes['shape']
			alpha = self.aes['alpha']
			ax.scatter(x, y, 
				s=size_scaler(self.data, index),
				c=colour_scaler(self.data, index),
				marker=shape_scaler(self.data, index),
				alpha=alpha(self.data, index))
		ax.set_xlabel(self.aes['x'])
		ax.set_ylabel(self.aes['y'])
		return fig, ax

class GeomDensity2D(Layer):
	def work(self, fig=None, ax=None):
		"""Render the layer on a matplotlib axis.
		You can specify either a figure or an axis to draw on.

		Parameters:
		-----------
		fig: matplotlib figure object
		ax: matplotlib axis object to draw on

		Returns:
		--------
		fig, ax: matplotlib figure and axis objects
		"""
		if ax is None:
			if fig is None:
				return fig, ax
			else:
				ax = fig.gca()
		x = self.data[self.aes['x']]
		y = self.data[self.aes['y']]
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
		ax.contour(Z, extent=[x_min, x_max, y_min, y_max])
		return fig, ax

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

class TrellisGrid(Layer):
	def __init__(self, by):
		"""Initialize TreelisGrid instance.

		Parameters:
		-----------
		by: column names to group by
		"""
		if len(by) != 2:
			raise ValueError("You must give a list of length 2 to group by")
		elif by[0] == '.' and by[1] == '.':
			raise ValueError("At least one of grouping attributes must be not a dot")
		self.by = by

	def trellis(self, layers):
		"""Create a trellis structure for a list of layers.
		Each layer will be cloned with different data in to a two dimensional grid.

		Parameters:
		-----------
		layers: a list of Layer objects

		Returns:
		--------
		trellised_layers: Clones of each layer in the list arranged in a trellised latice
		"""
		trellised_layers = []
		for layer in layers:
			data = layer.data
			if self.by[0] == '.':
				grouped = data.groupby(self.by[1])
			elif self.by[1] == '.':
				grouped = data.groupby(self.by[0])
			else:
				grouped = data.groupby(self.by)
			groups = grouped.groups.keys()
			if self.by[0] == '.' or self.by[1] == '.':
				shingle1 = set([g for g in groups])
			else:
				shingle1 = set([g[0] for g in groups])
				shingle2 = set([g[1] for g in groups])
			if self.by[0] == '.':
				self.rows = 1
				self.cols = len(shingle1)
			elif self.by[1] == '.':
				self.rows = len(shingle1)
				self.cols = 1
			else:
				self.rows = len(shingle1)
				self.cols = len(shingle2)
			trellised = [[None for _ in range(self.cols)] for _ in range(self.rows)]
			self.group_grid = [[None for _ in range(self.cols)] for _ in range(self.rows)]
			row = 0
			col = 0
			for group, data in grouped:
				new_layer = deepcopy(layer)
				new_layer.data = data
				trellised[row][col] = new_layer
				self.group_grid[row][col] = group
				col += 1
				if col >= self.cols:
					col = 0
					row += 1
			trellised_layers.append(trellised)
		return trellised_layers

def merge_aes(layer1, layer2):
	"""Merges the aesthetics dictionaries for the two layers. 
	Look up sequence_layers function. Which layer is first and which
	one is second is important.

	Parameters:
	-----------
	layer1: Layer object
	layer2: Layer object
	"""
	for key in layer2.aes.keys():
		if layer2.aes[key] is None:
			layer2.aes[key] = layer1.aes[key]

def sequence_layers(layers):
	"""Go through the list of layers and fill in the missing bits of information.
	The basic rules are this:
	* If the current layer has data set to None, take the data from previous layer.
	* For each aesthetic mapping, if that mapping is set to None, take it from previous layer.

	Parameters:
	-----------
	layers: a list of Layer objects
	"""
	for layer1, layer2 in zip(layers[:-1], layers[1:]):
		if layer2.data is None:
			layer2.data = layer1.data
		merge_aes(layer1, layer2)
	return layers

def sequence_grids(layer_grids):
	"""Go through the list of layer girds and perform the same thing as sequence_layers.

	Parameters:
	-----------
	layer_grids: a list of two dimensional layer grids
	"""
	for grid1, grid2 in zip(layer_grids[:-1], layer_grids[1:]):
		for row1, row2 in zip(grid1, grid2):
			for layer1, layer2 in zip(row1, row2):
				if layer2.data is None:
					layer2.data = layer1.data
				merge_aes(layer1, layer2)
	return layer_grids

def work_grid(grid, fig):
	"""Take a two dimensional grid, add subplots to a figure for each cell and do layer work.

	Parameters:
	-----------
	grid: a two dimensional grid of layers
	fig: matplotlib figure to draw on

	Returns:
	--------
	axes: a two dimensional list of matplotlib axes
	"""
	nrows = len(grid)
	ncols = len(grid[0])
	axes = [[None for _ in range(ncols)] for _ in range(nrows)]
	for row in range(nrows):
		for col in range(ncols):
			axes[row][col] = fig.add_subplot(nrows, ncols, ncols * row + col + 1)
			grid[row][col].work(ax=axes[row][col])
	return axes

def adjust_subplots(fig, axes, trellis):
	# Flatten the axes grid
	axes = [ax for row in axes for ax in row]
	min_x = min([ax.get_xlim()[0] for ax in axes])
	max_x = max([ax.get_xlim()[1] for ax in axes])
	min_y = min([ax.get_ylim()[0] for ax in axes])
	max_y = max([ax.get_ylim()[1] for ax in axes])
	[ax.set_xlim(min_x, max_x) for ax in axes]
	[ax.set_ylim(min_y, max_y) for ax in axes]
	for index, axis in enumerate(axes):
		if index % trellis.cols == 0:
			pass
		else:
			axis.get_yaxis().set_ticks([])
			axis.set_ylabel('')
		if index / trellis.cols == trellis.rows - 1:
			pass
		else:
			axis.get_xaxis().set_ticks([])
			axis.set_xlabel('')
		if trellis.by[0] == '.':
			label1 = "%s = %s" % (trellis.by[1], trellis.group_grid[index / trellis.cols][index % trellis.cols])
			label2 = None
		elif trellis.by[1] == '.':
			label1 = "%s = %s" % (trellis.by[0], trellis.group_grid[index / trellis.cols][index % trellis.cols])
			label2 = None
		else:
			label1 = "%s = %s" % (trellis.by[0], trellis.group_grid[index / trellis.cols][index % trellis.cols][0])
			label2 = "%s = %s" % (trellis.by[1], trellis.group_grid[index / trellis.cols][index % trellis.cols][1])
		if label2 is not None:
			axis.table(cellText=[[label1], [label2]], 
				loc='top', cellLoc='center', 
				cellColours=[['lightgrey'], ['lightgrey']])
		else:
			axis.table(cellText=[[label1]], loc='top', cellLoc='center', cellColours=[['lightgrey']])
	fig.subplots_adjust(wspace=0.05, hspace=0.2)
	return ax, fig

class RPlot:
	"""
	The main plot object. Add layers to an instance of this object to create a plot.
	"""
	def __init__(self, data, x=None, y=None):
		self.layers = [Layer(data, default_aes(x=x, y=y))]
		trellised = False

	def __add__(self, other):
		if not isinstance(other, Layer):
			raise TypeError("The operand on the right side of + must be a Layer instance")
		self.layers.append(other)

	def show(self, fig=None):
		if fig is None:
			fig = plt.gcf()
		# Look for the last TrellisGrid instance in the layer list
		last_trellis = None
		for layer in self.layers:
			if isinstance(layer, TrellisGrid):
				last_trellis = layer
		if last_trellis is None:
			# We have a simple, non-trellised plot
			new_layers = sequence_layers(self.layers)
			for layer in new_layers:
				layer.work(fig=fig)
		else:
			# We have a trellised plot.
			# First let's remove all other TrellisGrid instances from the layer list, 
			# including this one.
			new_layers = []
			for layer in self.layers:
				if not isinstance(layer, TrellisGrid):
					new_layers.append(layer)
			new_layers = sequence_layers(new_layers)
			# Now replace the old layers by their trellised versions
			new_layers = last_trellis.trellis(new_layers)
			# Prepare the subplots and draw on them
			new_layers = sequence_grids(new_layers)
			axes_grids = [work_grid(grid, fig) for grid in new_layers]
			axes_grid = axes_grids[-1]
			adjust_subplots(fig, axes_grid, last_trellis)
		# And we're done
		return fig