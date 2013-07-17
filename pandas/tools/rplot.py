import random
from copy import deepcopy
from pandas.core.common import _values_from_object

import numpy as np
from pandas.compat import range, zip
#
# TODO:
# * Make sure legends work properly
#

class Scale:
    """
    Base class for mapping between graphical and data attributes.
    """
    pass

class ScaleGradient(Scale):
    """
    A mapping between a data attribute value and a
    point in colour space between two specified colours.
    """
    def __init__(self, column, colour1, colour2):
        """Initialize ScaleGradient instance.

        Parameters:
        -----------
        column: string, pandas DataFrame column name
        colour1: tuple, 3 element tuple with float values representing an RGB colour
        colour2: tuple, 3 element tuple with float values representing an RGB colour
        """
        self.column = column
        self.colour1 = colour1
        self.colour2 = colour2
        self.categorical = False

    def __call__(self, data, index):
        """Return a colour corresponding to data attribute value.

        Parameters:
        -----------
        data: pandas DataFrame
        index: pandas DataFrame row index

        Returns:
        --------
        A three element tuple representing an RGB somewhere between colour1 and colour2
        """
        x = data[self.column].iget(index)
        a = min(data[self.column])
        b = max(data[self.column])
        r1, g1, b1 = self.colour1
        r2, g2, b2 = self.colour2
        x_scaled = (x - a) / (b - a)
        return (r1 + (r2 - r1) * x_scaled,
                g1 + (g2 - g1) * x_scaled,
                b1 + (b2 - b1) * x_scaled)

class ScaleGradient2(Scale):
    """
    Create a mapping between a data attribute value and a
    point in colour space in a line of three specified colours.
    """
    def __init__(self, column, colour1, colour2, colour3):
        """Initialize ScaleGradient2 instance.

        Parameters:
        -----------
        column: string, pandas DataFrame column name
        colour1: tuple, 3 element tuple with float values representing an RGB colour
        colour2: tuple, 3 element tuple with float values representing an RGB colour
        colour3: tuple, 3 element tuple with float values representing an RGB colour
        """
        self.column = column
        self.colour1 = colour1
        self.colour2 = colour2
        self.colour3 = colour3
        self.categorical = False

    def __call__(self, data, index):
        """Return a colour corresponding to data attribute value.

        Parameters:
        -----------
        data: pandas DataFrame
        index: pandas DataFrame row index

        Returns:
        --------
        A three element tuple representing an RGB somewhere along the line
        of colour1, colour2 and colour3
        """
        x = data[self.column].iget(index)
        a = min(data[self.column])
        b = max(data[self.column])
        r1, g1, b1 = self.colour1
        r2, g2, b2 = self.colour2
        r3, g3, b3 = self.colour3
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

class ScaleSize(Scale):
    """
    Provide a mapping between a DataFrame column and matplotlib
    scatter plot shape size.
    """
    def __init__(self, column, min_size=5.0, max_size=100.0, transform=lambda x: x):
        """Initialize ScaleSize instance.

        Parameters:
        -----------
        column: string, a column name
        min_size: float, minimum point size
        max_size: float, maximum point size
        transform: a one argument function of form float -> float (e.g. lambda x: log(x))
        """
        self.column = column
        self.min_size = min_size
        self.max_size = max_size
        self.transform = transform
        self.categorical = False

    def __call__(self, data, index):
        """Return matplotlib scatter plot marker shape size.

        Parameters:
        -----------
        data: pandas DataFrame
        index: pandas DataFrame row index
        """
        x = data[self.column].iget(index)
        a = float(min(data[self.column]))
        b = float(max(data[self.column]))
        return self.transform(self.min_size + ((x - a) / (b - a)) *
            (self.max_size - self.min_size))

class ScaleShape(Scale):
    """
    Provides a mapping between matplotlib marker shapes
    and attribute values.
    """
    def __init__(self, column):
        """Initialize ScaleShape instance.

        Parameters:
        -----------
        column: string, pandas DataFrame column name
        """
        self.column = column
        self.shapes = ['o', '+', 's', '*', '^', '<', '>', 'v', '|', 'x']
        self.legends = set([])
        self.categorical = True

    def __call__(self, data, index):
        """Returns a matplotlib marker identifier.

        Parameters:
        -----------
        data: pandas DataFrame
        index: pandas DataFrame row index

        Returns:
        --------
        a matplotlib marker identifier
        """
        values = sorted(list(set(data[self.column])))
        if len(values) > len(self.shapes):
            raise ValueError("Too many different values of the categorical attribute for ScaleShape")
        x = data[self.column].iget(index)
        return self.shapes[values.index(x)]

class ScaleRandomColour(Scale):
    """
    Maps a random colour to a DataFrame attribute.
    """
    def __init__(self, column):
        """Initialize ScaleRandomColour instance.

        Parameters:
        -----------
        column: string, pandas DataFrame column name
        """
        self.column = column
        self.categorical = True

    def __call__(self, data, index):
        """Return a tuple of three floats, representing
        an RGB colour.

        Parameters:
        -----------
        data: pandas DataFrame
        index: pandas DataFrame row index
        """
        random.seed(data[self.column].iget(index))
        return [random.random() for _ in range(3)]

class ScaleConstant(Scale):
    """
    Constant returning scale. Usually used automatically.
    """
    def __init__(self, value):
        """Initialize ScaleConstant instance.

        Parameters:
        -----------
        value: any Python value to be returned when called
        """
        self.value = value
        self.categorical = False

    def __call__(self, data, index):
        """Return the constant value.

        Parameters:
        -----------
        data: pandas DataFrame
        index: pandas DataFrame row index

        Returns:
        --------
        A constant value specified during initialisation
        """
        return self.value

def default_aes(x=None, y=None):
    """Create the default aesthetics dictionary.

    Parameters:
    -----------
    x: string, DataFrame column name
    y: string, DataFrame column name

    Returns:
    --------
    a dictionary with aesthetics bindings
    """
    return {
        'x' : x,
        'y' : y,
        'size' : ScaleConstant(40.0),
        'colour' : ScaleConstant('grey'),
        'shape' : ScaleConstant('o'),
        'alpha' : ScaleConstant(1.0),
    }

def make_aes(x=None, y=None, size=None, colour=None, shape=None, alpha=None):
    """Create an empty aesthetics dictionary.

    Parameters:
    -----------
    x: string, DataFrame column name
    y: string, DataFrame column name
    size: function, binding for size attribute of Geoms
    colour: function, binding for colour attribute of Geoms
    shape: function, binding for shape attribute of Geoms
    alpha: function, binding for alpha attribute of Geoms

    Returns:
    --------
    a dictionary with aesthetics bindings
    """
    if not hasattr(size, '__call__') and size is not None:
        size = ScaleConstant(size)
    if not hasattr(colour, '__call__') and colour is not None:
        colour = ScaleConstant(colour)
    if not hasattr(shape, '__call__') and shape is not None:
        shape = ScaleConstant(shape)
    if not hasattr(alpha, '__call__') and alpha is not None:
        alpha = ScaleConstant(alpha)
    if any([isinstance(size, scale) for scale in [ScaleConstant, ScaleSize]]) or size is None:
        pass
    else:
        raise ValueError('size mapping should be done through ScaleConstant or ScaleSize')
    if any([isinstance(colour, scale) for scale in [ScaleConstant, ScaleGradient, ScaleGradient2, ScaleRandomColour]]) or colour is None:
        pass
    else:
        raise ValueError('colour mapping should be done through ScaleConstant, ScaleRandomColour, ScaleGradient or ScaleGradient2')
    if any([isinstance(shape, scale) for scale in [ScaleConstant, ScaleShape]]) or shape is None:
        pass
    else:
        raise ValueError('shape mapping should be done through ScaleConstant or ScaleShape')
    if any([isinstance(alpha, scale) for scale in [ScaleConstant]]) or alpha is None:
        pass
    else:
        raise ValueError('alpha mapping should be done through ScaleConstant')
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
    def __init__(self, data=None, **kwds):
        """Initialize layer object.

        Parameters:
        -----------
        data: pandas DataFrame instance
        aes: aesthetics dictionary with bindings
        """
        self.data = data
        self.aes = make_aes(**kwds)
        self.legend = {}

    def work(self, fig=None, ax=None):
        """Do the drawing (usually) work.

        Parameters:
        -----------
        fig: matplotlib figure
        ax: matplotlib axis object

        Returns:
        --------
        a tuple with the same figure and axis instances
        """
        return fig, ax

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
            size_value = size_scaler(self.data, index)
            colour_value = colour_scaler(self.data, index)
            marker_value = shape_scaler(self.data, index)
            alpha_value = alpha(self.data, index)
            patch = ax.scatter(x, y,
                    s=size_value,
                    c=colour_value,
                    marker=marker_value,
                    alpha=alpha_value)
            label = []
            if colour_scaler.categorical:
                label += [colour_scaler.column, row[colour_scaler.column]]
            if shape_scaler.categorical:
                label += [shape_scaler.column, row[shape_scaler.column]]
            self.legend[tuple(label)] = patch
        ax.set_xlabel(self.aes['x'])
        ax.set_ylabel(self.aes['y'])
        return fig, ax

class GeomPolyFit(Layer):
    """
    Draw a polynomial fit of specified degree.
    """
    def __init__(self, degree, lw=2.0, colour='grey'):
        """Initialize GeomPolyFit object.

        Parameters:
        -----------
        degree: an integer, polynomial degree
        lw: line width
        colour: matplotlib colour
        """
        self.degree = degree
        self.lw = lw
        self.colour = colour
        Layer.__init__(self)

    def work(self, fig=None, ax=None):
        """Draw the polynomial fit on matplotlib figure or axis

        Parameters:
        -----------
        fig: matplotlib figure
        ax: matplotlib axis

        Returns:
        --------
        a tuple with figure and axis objects
        """
        if ax is None:
            if fig is None:
                return fig, ax
            else:
                ax = fig.gca()
        from numpy.polynomial.polynomial import polyfit
        from numpy.polynomial.polynomial import polyval
        x = self.data[self.aes['x']]
        y = self.data[self.aes['y']]
        min_x = min(x)
        max_x = max(x)
        c = polyfit(x, y, self.degree)
        x_ = np.linspace(min_x, max_x, len(x))
        y_ = polyval(x_, c)
        ax.plot(x_, y_, lw=self.lw, c=self.colour)
        return fig, ax

class GeomScatter(Layer):
    """
    An efficient scatter plot, use this instead of GeomPoint for speed.
    """
    def __init__(self, marker='o', colour='lightblue', alpha=1.0):
        """Initialize GeomScatter instance.

        Parameters:
        -----------
        marker: matplotlib marker string
        colour: matplotlib colour
        alpha: matplotlib alpha
        """
        self.marker = marker
        self.colour = colour
        self.alpha = alpha
        Layer.__init__(self)

    def work(self, fig=None, ax=None):
        """Draw a scatter plot on matplotlib figure or axis

        Parameters:
        -----------
        fig: matplotlib figure
        ax: matplotlib axis

        Returns:
        --------
        a tuple with figure and axis objects
        """
        if ax is None:
            if fig is None:
                return fig, ax
            else:
                ax = fig.gca()
        x = self.data[self.aes['x']]
        y = self.data[self.aes['y']]
        ax.scatter(x, y, marker=self.marker, c=self.colour, alpha=self.alpha)
        return fig, ax

class GeomHistogram(Layer):
    """
    An efficient histogram, use this instead of GeomBar for speed.
    """
    def __init__(self, bins=10, colour='lightblue'):
        """Initialize GeomHistogram instance.

        Parameters:
        -----------
        bins: integer, number of histogram bins
        colour: matplotlib colour
        """
        self.bins = bins
        self.colour = colour
        Layer.__init__(self)

    def work(self, fig=None, ax=None):
        """Draw a histogram on matplotlib figure or axis

        Parameters:
        -----------
        fig: matplotlib figure
        ax: matplotlib axis

        Returns:
        --------
        a tuple with figure and axis objects
        """
        if ax is None:
            if fig is None:
                return fig, ax
            else:
                ax = fig.gca()
        x = self.data[self.aes['x']]
        ax.hist(_values_from_object(x), self.bins, facecolor=self.colour)
        ax.set_xlabel(self.aes['x'])
        return fig, ax

class GeomDensity(Layer):
    """
    A kernel density estimation plot.
    """
    def work(self, fig=None, ax=None):
        """Draw a one dimensional kernel density plot.
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
        from scipy.stats import gaussian_kde
        x = self.data[self.aes['x']]
        gkde = gaussian_kde(x)
        ind = np.linspace(x.min(), x.max(), 200)
        ax.plot(ind, gkde.evaluate(ind))
        return fig, ax

class GeomDensity2D(Layer):
    def work(self, fig=None, ax=None):
        """Draw a two dimensional kernel density plot.
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
        import scipy.stats as stats
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)
        ax.contour(Z, extent=[x_min, x_max, y_min, y_max])
        return fig, ax

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
            groups = list(grouped.groups.keys())
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

def dictionary_union(dict1, dict2):
    """Take two dictionaries, return dictionary union.

    Parameters:
    -----------
    dict1: Python dictionary
    dict2: Python dictionary

    Returns:
    --------
    A union of the dictionaries. It assumes that values
    with the same keys are identical.
    """
    keys1 = list(dict1.keys())
    keys2 = list(dict2.keys())
    result = {}
    for key1 in keys1:
        result[key1] = dict1[key1]
    for key2 in keys2:
        result[key2] = dict2[key2]
    return result

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

def adjust_subplots(fig, axes, trellis, layers):
    """Adjust the subtplots on matplotlib figure with the
    fact that we have a trellis plot in mind.

    Parameters:
    -----------
    fig: matplotlib figure
    axes: a two dimensional grid of matplotlib axes
    trellis: TrellisGrid object
    layers: last grid of layers in the plot
    """
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
            label1 = "%s = %s" % (trellis.by[1], trellis.group_grid[index // trellis.cols][index % trellis.cols])
            label2 = None
        elif trellis.by[1] == '.':
            label1 = "%s = %s" % (trellis.by[0], trellis.group_grid[index // trellis.cols][index % trellis.cols])
            label2 = None
        else:
            label1 = "%s = %s" % (trellis.by[0], trellis.group_grid[index // trellis.cols][index % trellis.cols][0])
            label2 = "%s = %s" % (trellis.by[1], trellis.group_grid[index // trellis.cols][index % trellis.cols][1])
        if label2 is not None:
            axis.table(cellText=[[label1], [label2]],
                loc='top', cellLoc='center',
                cellColours=[['lightgrey'], ['lightgrey']])
        else:
            axis.table(cellText=[[label1]], loc='top', cellLoc='center', cellColours=[['lightgrey']])
    # Flatten the layer grid
    layers = [layer for row in layers for layer in row]
    legend = {}
    for layer in layers:
        legend = dictionary_union(legend, layer.legend)
    patches = []
    labels = []
    if len(list(legend.keys())) == 0:
        key_function = lambda tup: tup
    elif len(list(legend.keys())[0]) == 2:
        key_function = lambda tup: (tup[1])
    else:
        key_function = lambda tup: (tup[1], tup[3])
    for key in sorted(list(legend.keys()), key=key_function):
        value = legend[key]
        patches.append(value)
        if len(key) == 2:
            col, val = key
            labels.append("%s" % str(val))
        elif len(key) == 4:
            col1, val1, col2, val2 = key
            labels.append("%s, %s" % (str(val1), str(val2)))
        else:
            raise ValueError("Maximum 2 categorical attributes to display a lengend of")
    if len(legend):
        fig.legend(patches, labels, loc='upper right')
    fig.subplots_adjust(wspace=0.05, hspace=0.2)

class RPlot:
    """
    The main plot object. Add layers to an instance of this object to create a plot.
    """
    def __init__(self, data, x=None, y=None):
        """Initialize RPlot instance.

        Parameters:
        -----------
        data: pandas DataFrame instance
        x: string, DataFrame column name
        y: string, DataFrame column name
        """
        self.layers = [Layer(data, **default_aes(x=x, y=y))]
        trellised = False

    def add(self, layer):
        """Add a layer to RPlot instance.

        Parameters:
        -----------
        layer: Layer instance
        """
        if not isinstance(layer, Layer):
            raise TypeError("The operand on the right side of + must be a Layer instance")
        self.layers.append(layer)

    def render(self, fig=None):
        """Render all the layers on a matplotlib figure.

        Parameters:
        -----------
        fig: matplotlib figure
        """
        import matplotlib.pyplot as plt
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
            legend = {}
            for layer in new_layers:
                legend = dictionary_union(legend, layer.legend)
            patches = []
            labels = []
            if len(list(legend.keys())) == 0:
                key_function = lambda tup: tup
            elif len(list(legend.keys())[0]) == 2:
                key_function = lambda tup: (tup[1])
            else:
                key_function = lambda tup: (tup[1], tup[3])
            for key in sorted(list(legend.keys()), key=key_function):
                value = legend[key]
                patches.append(value)
                if len(key) == 2:
                    col, val = key
                    labels.append("%s" % str(val))
                elif len(key) == 4:
                    col1, val1, col2, val2 = key
                    labels.append("%s, %s" % (str(val1), str(val2)))
                else:
                    raise ValueError("Maximum 2 categorical attributes to display a lengend of")
            if len(legend):
                fig.legend(patches, labels, loc='upper right')
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
            adjust_subplots(fig, axes_grid, last_trellis, new_layers[-1])
        # And we're done
        return fig
