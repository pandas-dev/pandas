.. currentmodule:: pandas
.. _visualization:

.. ipython:: python
   :suppress:

   import numpy as np
   from numpy.random import randn, rand, randint
   np.random.seed(123456)
   from pandas import DataFrame, Series, date_range, options
   import pandas.util.testing as tm
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')
   options.display.mpl_style = 'default'

************************
Plotting with matplotlib
************************

.. note::

    We intend to build more plotting integration with `matplotlib
    <http://matplotlib.sourceforge.net>`__ as time goes on.

We use the standard convention for referencing the matplotlib API:

.. ipython:: python

   import matplotlib.pyplot as plt

.. _visualization.basic:

Basic plotting: ``plot``
------------------------

See the :ref:`cookbook<cookbook.plotting>` for some advanced strategies

The ``plot`` method on Series and DataFrame is just a simple wrapper around
``plt.plot``:

.. ipython:: python

   ts = Series(randn(1000), index=date_range('1/1/2000', periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png
   ts.plot()

If the index consists of dates, it calls ``gcf().autofmt_xdate()`` to try to
format the x-axis nicely as per above. The method takes a number of arguments
for controlling the look of the plot:

.. ipython:: python

   @savefig series_plot_basic2.png
   plt.figure(); ts.plot(style='k--', label='Series'); plt.legend()

On DataFrame, ``plot`` is a convenience to plot all of the columns with labels:

.. ipython:: python

   df = DataFrame(randn(1000, 4), index=ts.index, columns=list('ABCD'))
   df = df.cumsum()

   @savefig frame_plot_basic.png
   plt.figure(); df.plot(); plt.legend(loc='best')

You may set the ``legend`` argument to ``False`` to hide the legend, which is
shown by default.

.. ipython:: python

   @savefig frame_plot_basic_noleg.png
   df.plot(legend=False)

Some other options are available, like plotting each Series on a different axis:

.. ipython:: python

   @savefig frame_plot_subplots.png
   df.plot(subplots=True, figsize=(6, 6)); plt.legend(loc='best')

You may pass ``logy`` to get a log-scale Y axis.

.. ipython:: python

   plt.figure();

   ts = Series(randn(1000), index=date_range('1/1/2000', periods=1000))
   ts = np.exp(ts.cumsum())

   @savefig series_plot_logy.png
   ts.plot(logy=True)

You can plot one column versus another using the `x` and `y` keywords in
`DataFrame.plot`:

.. ipython:: python

   plt.figure()

   df3 = DataFrame(randn(1000, 2), columns=['B', 'C']).cumsum()
   df3['A'] = Series(range(len(df)))

   @savefig df_plot_xy.png
   df3.plot(x='A', y='B')


Plotting on a Secondary Y-axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To plot data on a secondary y-axis, use the ``secondary_y`` keyword:

.. ipython:: python

   plt.figure()

   df.A.plot()

   @savefig series_plot_secondary_y.png
   df.B.plot(secondary_y=True, style='g')


Selective Plotting on Secondary Y-axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To plot some columns in a DataFrame, give the column names to the ``secondary_y``
keyword:

.. ipython:: python

   plt.figure()
   ax = df.plot(secondary_y=['A', 'B'])
   ax.set_ylabel('CD scale')
   @savefig frame_plot_secondary_y.png
   ax.right_ax.set_ylabel('AB scale')



Note that the columns plotted on the secondary y-axis is automatically marked
with "(right)" in the legend. To turn off the automatic marking, use the
``mark_right=False`` keyword:

.. ipython:: python

   plt.figure()

   @savefig frame_plot_secondary_y_no_right.png
   df.plot(secondary_y=['A', 'B'], mark_right=False)


Suppressing tick resolution adjustment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pandas includes automatically tick resolution adjustment for regular frequency
time-series data. For limited cases where pandas cannot infer the frequency
information (e.g., in an externally created ``twinx``), you can choose to
suppress this behavior for alignment purposes.

Here is the default behavior, notice how the x-axis tick labelling is performed:

.. ipython:: python

   plt.figure()

   @savefig ser_plot_suppress.png
   df.A.plot()


Using the ``x_compat`` parameter, you can suppress this behavior:

.. ipython:: python

   plt.figure()

   @savefig ser_plot_suppress_parm.png
   df.A.plot(x_compat=True)


If you have more than one plot that needs to be suppressed, the ``use`` method
in ``pandas.plot_params`` can be used in a `with statement`:

.. ipython:: python

   import pandas as pd

   plt.figure()

   @savefig ser_plot_suppress_context.png
   with pd.plot_params.use('x_compat', True):
       df.A.plot(color='r')
       df.B.plot(color='g')
       df.C.plot(color='b')


Targeting different subplots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can pass an ``ax`` argument to ``Series.plot`` to plot on a particular axis:

.. ipython:: python
   :suppress:

   ts = Series(randn(1000), index=date_range('1/1/2000', periods=1000))
   ts = ts.cumsum()

   df = DataFrame(randn(1000, 4), index=ts.index, columns=list('ABCD'))
   df = df.cumsum()

.. ipython:: python

   fig, axes = plt.subplots(nrows=2, ncols=2)
   df['A'].plot(ax=axes[0,0]); axes[0,0].set_title('A')
   df['B'].plot(ax=axes[0,1]); axes[0,1].set_title('B')
   df['C'].plot(ax=axes[1,0]); axes[1,0].set_title('C')

   @savefig series_plot_multi.png
   df['D'].plot(ax=axes[1,1]); axes[1,1].set_title('D')


.. _visualization.other:

Other plotting features
-----------------------

.. _visualization.barplot:

Bar plots
~~~~~~~~~

For labeled, non-time series data, you may wish to produce a bar plot:

.. ipython:: python

   plt.figure();

   @savefig bar_plot_ex.png
   df.ix[5].plot(kind='bar'); plt.axhline(0, color='k')

Calling a DataFrame's ``plot`` method with ``kind='bar'`` produces a multiple
bar plot:

.. ipython:: python
   :suppress:

   plt.figure();

.. ipython:: python

   df2 = DataFrame(rand(10, 4), columns=['a', 'b', 'c', 'd'])

   @savefig bar_plot_multi_ex.png
   df2.plot(kind='bar');

To produce a stacked bar plot, pass ``stacked=True``:

.. ipython:: python
   :suppress:

   plt.figure();

.. ipython:: python

   @savefig bar_plot_stacked_ex.png
   df2.plot(kind='bar', stacked=True);

To get horizontal bar plots, pass ``kind='barh'``:

.. ipython:: python
   :suppress:

   plt.figure();

.. ipython:: python

   @savefig barh_plot_stacked_ex.png
   df2.plot(kind='barh', stacked=True);

Histograms
~~~~~~~~~~
.. ipython:: python

   plt.figure();

   @savefig hist_plot_ex.png
   df['A'].diff().hist()


For a DataFrame, ``hist`` plots the histograms of the columns on multiple
subplots:

.. ipython:: python

   plt.figure()

   @savefig frame_hist_ex.png
   df.diff().hist(color='k', alpha=0.5, bins=50)


New since 0.10.0, the ``by`` keyword can be specified to plot grouped histograms:

.. ipython:: python
   :suppress:

   plt.figure();

.. ipython:: python

   data = Series(randn(1000))

   @savefig grouped_hist.png
   data.hist(by=randint(0, 4, 1000), figsize=(6, 4))


.. _visualization.box:

Box-Plotting
~~~~~~~~~~~~

DataFrame has a ``boxplot`` method which allows you to visualize the
distribution of values within each column.

For instance, here is a boxplot representing five trials of 10 observations of
a uniform random variable on [0,1).

.. ipython:: python

   df = DataFrame(rand(10,5))
   plt.figure();

   @savefig box_plot_ex.png
   bp = df.boxplot()

You can create a stratified boxplot using the ``by`` keyword argument to create
groupings.  For instance,

.. ipython:: python

   df = DataFrame(rand(10,2), columns=['Col1', 'Col2'] )
   df['X'] = Series(['A','A','A','A','A','B','B','B','B','B'])

   plt.figure();

   @savefig box_plot_ex2.png
   bp = df.boxplot(by='X')

You can also pass a subset of columns to plot, as well as group by multiple
columns:

.. ipython:: python

   df = DataFrame(rand(10,3), columns=['Col1', 'Col2', 'Col3'])
   df['X'] = Series(['A','A','A','A','A','B','B','B','B','B'])
   df['Y'] = Series(['A','B','A','B','A','B','A','B','A','B'])

   plt.figure();

   @savefig box_plot_ex3.png
   bp = df.boxplot(column=['Col1','Col2'], by=['X','Y'])

.. _visualization.scatter_matrix:

Scatter plot matrix
~~~~~~~~~~~~~~~~~~~

*New in 0.7.3.* You can create a scatter plot matrix using the
 ``scatter_matrix`` method in ``pandas.tools.plotting``:

.. ipython:: python

   from pandas.tools.plotting import scatter_matrix
   df = DataFrame(randn(1000, 4), columns=['a', 'b', 'c', 'd'])

   @savefig scatter_matrix_kde.png
   scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal='kde')

.. _visualization.kde:

*New in 0.8.0* You can create density plots using the Series/DataFrame.plot and
setting `kind='kde'`:

.. ipython:: python
   :suppress:

   plt.figure();

.. ipython:: python

   ser = Series(randn(1000))

   @savefig kde_plot.png
   ser.plot(kind='kde')

.. _visualization.andrews_curves:

Andrews Curves
~~~~~~~~~~~~~~

Andrews curves allow one to plot multivariate data as a large number
of curves that are created using the attributes of samples as coefficients
for Fourier series. By coloring these curves differently for each class
it is possible to visualize data clustering. Curves belonging to samples
of the same class will usually be closer together and form larger structures.

**Note**: The "Iris" dataset is available `here <https://raw.github.com/pydata/pandas/master/pandas/tests/data/iris.csv>`__.

.. ipython:: python

   from pandas import read_csv
   from pandas.tools.plotting import andrews_curves

   data = read_csv('data/iris.data')

   plt.figure()

   @savefig andrews_curves.png
   andrews_curves(data, 'Name')

.. _visualization.parallel_coordinates:

Parallel Coordinates
~~~~~~~~~~~~~~~~~~~~

Parallel coordinates is a plotting technique for plotting multivariate data.
It allows one to see clusters in data and to estimate other statistics visually.
Using parallel coordinates points are represented as connected line segments.
Each vertical line represents one attribute. One set of connected line segments
represents one data point. Points that tend to cluster will appear closer together.

.. ipython:: python

   from pandas import read_csv
   from pandas.tools.plotting import parallel_coordinates

   data = read_csv('data/iris.data')

   plt.figure()

   @savefig parallel_coordinates.png
   parallel_coordinates(data, 'Name')

Lag Plot
~~~~~~~~

Lag plots are used to check if a data set or time series is random. Random
data should not exhibit any structure in the lag plot. Non-random structure
implies that the underlying data are not random.

.. ipython:: python

   from pandas.tools.plotting import lag_plot

   plt.figure()

   data = Series(0.1 * rand(1000) +
      0.9 * np.sin(np.linspace(-99 * np.pi, 99 * np.pi, num=1000)))

   @savefig lag_plot.png
   lag_plot(data)

Autocorrelation Plot
~~~~~~~~~~~~~~~~~~~~

Autocorrelation plots are often used for checking randomness in time series.
This is done by computing autocorrelations for data values at varying time lags.
If time series is random, such autocorrelations should be near zero for any and
all time-lag separations. If time series is non-random then one or more of the
autocorrelations will be significantly non-zero. The horizontal lines displayed
in the plot correspond to 95% and 99% confidence bands. The dashed line is 99%
confidence band.

.. ipython:: python

   from pandas.tools.plotting import autocorrelation_plot

   plt.figure()

   data = Series(0.7 * rand(1000) +
      0.3 * np.sin(np.linspace(-9 * np.pi, 9 * np.pi, num=1000)))

   @savefig autocorrelation_plot.png
   autocorrelation_plot(data)

.. _visualization.bootstrap:

Bootstrap Plot
~~~~~~~~~~~~~~

Bootstrap plots are used to visually assess the uncertainty of a statistic, such
as mean, median, midrange, etc. A random subset of a specified size is selected
from a data set, the statistic in question is computed for this subset and the
process is repeated a specified number of times. Resulting plots and histograms
are what constitutes the bootstrap plot.

.. ipython:: python

   from pandas.tools.plotting import bootstrap_plot

   data = Series(rand(1000))

   @savefig bootstrap_plot.png
   bootstrap_plot(data, size=50, samples=500, color='grey')

.. _visualization.radviz:

RadViz
~~~~~~

RadViz is a way of visualizing multi-variate data. It is based on a simple
spring tension minimization algorithm. Basically you set up a bunch of points in
a plane. In our case they are equally spaced on a unit circle. Each point
represents a single attribute. You then pretend that each sample in the data set
is attached to each of these points by a spring, the stiffness of which is
proportional to the numerical value of that attribute (they are normalized to
unit interval). The point in the plane, where our sample settles to (where the
forces acting on our sample are at an equilibrium) is where a dot representing
our sample will be drawn. Depending on which class that sample belongs it will
be colored differently.

**Note**: The "Iris" dataset is available `here <https://raw.github.com/pydata/pandas/master/pandas/tests/data/iris.csv>`__.

.. ipython:: python

   from pandas import read_csv
   from pandas.tools.plotting import radviz

   data = read_csv('data/iris.data')

   plt.figure()

   @savefig radviz.png
   radviz(data, 'Name')

.. _visualization.colormaps:

Colormaps
~~~~~~~~~

A potential issue when plotting a large number of columns is that it can be difficult to distinguish some series due to repetition in the default colors. To remedy this, DataFrame plotting supports the use of the ``colormap=`` argument, which accepts either a Matplotlib `colormap <http://matplotlib.org/api/cm_api.html>`__ or a string that is a name of a colormap registered with Matplotlib. A visualization of the default matplotlib colormaps is available `here <http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps>`__.

As matplotlib does not directly support colormaps for line-based plots, the colors are selected based on an even spacing determined by the number of columns in the DataFrame. There is no consideration made for background color, so some colormaps will produce lines that are not easily visible.

To use the jet colormap, we can simply pass ``'jet'`` to ``colormap=``

.. ipython:: python

   df = DataFrame(randn(1000, 10), index=ts.index)
   df = df.cumsum()

   plt.figure()

   @savefig jet.png
   df.plot(colormap='jet')

or we can pass the colormap itself

.. ipython:: python

   from matplotlib import cm

   plt.figure()

   @savefig jet_cm.png
   df.plot(colormap=cm.jet)

Colormaps can also be used other plot types, like bar charts:

.. ipython:: python

   dd = DataFrame(randn(10, 10)).applymap(abs)
   dd = dd.cumsum()

   plt.figure()

   @savefig greens.png
   dd.plot(kind='bar', colormap='Greens')

Parallel coordinates charts:

.. ipython:: python

   plt.figure()

   @savefig parallel_gist_rainbow.png
   parallel_coordinates(data, 'Name', colormap='gist_rainbow')

Andrews curves charts:

.. ipython:: python

   plt.figure()

   @savefig andrews_curve_winter.png
   andrews_curves(data, 'Name', colormap='winter')
