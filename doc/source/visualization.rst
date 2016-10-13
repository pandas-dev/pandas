.. currentmodule:: pandas
.. _visualization:

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows = 15
   import matplotlib
   matplotlib.style.use('ggplot')
   import matplotlib.pyplot as plt
   plt.close('all')

*************
Visualization
*************

We use the standard convention for referencing the matplotlib API:

.. ipython:: python

   import matplotlib.pyplot as plt

The plots in this document are made using matplotlib's ``ggplot`` style (new in version 1.4):

.. code-block:: python

   import matplotlib
   matplotlib.style.use('ggplot')

We provide the basics in pandas to easily create decent looking plots.
See the :ref:`ecosystem <ecosystem.visualization>` section for visualization
libraries that go beyond the basics documented here.

.. note::

   All calls to ``np.random`` are seeded with 123456.

.. _visualization.basic:

Basic Plotting: ``plot``
------------------------

See the :ref:`cookbook<cookbook.plotting>` for some advanced strategies

The ``plot`` method on Series and DataFrame is just a simple wrapper around
:meth:`plt.plot() <matplotlib.axes.Axes.plot>`:

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png
   ts.plot()

If the index consists of dates, it calls :meth:`gcf().autofmt_xdate() <matplotlib.figure.Figure.autofmt_xdate>`
to try to format the x-axis nicely as per above.

On DataFrame, :meth:`~DataFrame.plot` is a convenience to plot all of the columns with labels:

.. ipython:: python
   :suppress:

   plt.close('all')
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
   df = df.cumsum()

   @savefig frame_plot_basic.png
   plt.figure(); df.plot();

You can plot one column versus another using the `x` and `y` keywords in
:meth:`~DataFrame.plot`:

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df3 = pd.DataFrame(np.random.randn(1000, 2), columns=['B', 'C']).cumsum()
   df3['A'] = pd.Series(list(range(len(df))))

   @savefig df_plot_xy.png
   df3.plot(x='A', y='B')

.. note::

   For more formatting and styling options, see :ref:`below <visualization.formatting>`.

.. ipython:: python
    :suppress:

    plt.close('all')

.. _visualization.other:

Other Plots
-----------

Plotting methods allow for a handful of plot styles other than the
default Line plot. These methods can be provided as the ``kind``
keyword argument to :meth:`~DataFrame.plot`.
These include:

* :ref:`'bar' <visualization.barplot>` or :ref:`'barh' <visualization.barplot>` for bar plots
* :ref:`'hist' <visualization.hist>` for histogram
* :ref:`'box' <visualization.box>` for boxplot
* :ref:`'kde' <visualization.kde>` or ``'density'`` for density plots
* :ref:`'area' <visualization.area_plot>` for area plots
* :ref:`'scatter' <visualization.scatter>` for scatter plots
* :ref:`'hexbin' <visualization.hexbin>` for hexagonal bin plots
* :ref:`'pie' <visualization.pie>` for pie plots

For example, a bar plot can be created the following way:

.. ipython:: python

   plt.figure();

   @savefig bar_plot_ex.png
   df.ix[5].plot(kind='bar'); plt.axhline(0, color='k')

.. versionadded:: 0.17.0

You can also create these other plots using the methods ``DataFrame.plot.<kind>`` instead of providing the ``kind`` keyword argument. This makes it easier to discover plot methods and the specific arguments they use:

.. ipython::
    :verbatim:

    In [14]: df = pd.DataFrame()

    In [15]: df.plot.<TAB>
    df.plot.area     df.plot.barh     df.plot.density  df.plot.hist     df.plot.line     df.plot.scatter
    df.plot.bar      df.plot.box      df.plot.hexbin   df.plot.kde      df.plot.pie

In addition to these ``kind`` s, there are  the :ref:`DataFrame.hist() <visualization.hist>`,
and :ref:`DataFrame.boxplot() <visualization.box>` methods, which use a separate interface.

Finally, there are several :ref:`plotting functions <visualization.tools>` in ``pandas.tools.plotting``
that take a :class:`Series` or :class:`DataFrame` as an argument. These
include

* :ref:`Scatter Matrix <visualization.scatter_matrix>`
* :ref:`Andrews Curves <visualization.andrews_curves>`
* :ref:`Parallel Coordinates <visualization.parallel_coordinates>`
* :ref:`Lag Plot <visualization.lag>`
* :ref:`Autocorrelation Plot <visualization.autocorrelation>`
* :ref:`Bootstrap Plot <visualization.bootstrap>`
* :ref:`RadViz <visualization.radviz>`

Plots may also be adorned with :ref:`errorbars <visualization.errorbars>`
or :ref:`tables <visualization.table>`.

.. _visualization.barplot:

Bar plots
~~~~~~~~~

For labeled, non-time series data, you may wish to produce a bar plot:

.. ipython:: python

   plt.figure();

   @savefig bar_plot_ex.png
   df.ix[5].plot.bar(); plt.axhline(0, color='k')

Calling a DataFrame's :meth:`plot.bar() <DataFrame.plot.bar>` method produces a multiple
bar plot:

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df2 = pd.DataFrame(np.random.rand(10, 4), columns=['a', 'b', 'c', 'd'])

   @savefig bar_plot_multi_ex.png
   df2.plot.bar();

To produce a stacked bar plot, pass ``stacked=True``:

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()

.. ipython:: python

   @savefig bar_plot_stacked_ex.png
   df2.plot.bar(stacked=True);

To get horizontal bar plots, use the ``barh`` method:

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()

.. ipython:: python

   @savefig barh_plot_stacked_ex.png
   df2.plot.barh(stacked=True);

.. _visualization.hist:

Histograms
~~~~~~~~~~

.. versionadded:: 0.15.0

Histogram can be drawn by using the :meth:`DataFrame.plot.hist` and :meth:`Series.plot.hist` methods.

.. ipython:: python

   df4 = pd.DataFrame({'a': np.random.randn(1000) + 1, 'b': np.random.randn(1000),
                       'c': np.random.randn(1000) - 1}, columns=['a', 'b', 'c'])

   plt.figure();

   @savefig hist_new.png
   df4.plot.hist(alpha=0.5)


.. ipython:: python
   :suppress:

   plt.close('all')

Histogram can be stacked by ``stacked=True``. Bin size can be changed by ``bins`` keyword.

.. ipython:: python

   plt.figure();

   @savefig hist_new_stacked.png
   df4.plot.hist(stacked=True, bins=20)

.. ipython:: python
   :suppress:

   plt.close('all')

You can pass other keywords supported by matplotlib ``hist``. For example, horizontal and cumulative histgram can be drawn by ``orientation='horizontal'`` and ``cumulative='True'``.

.. ipython:: python

   plt.figure();

   @savefig hist_new_kwargs.png
   df4['a'].plot.hist(orientation='horizontal', cumulative=True)

.. ipython:: python
   :suppress:

   plt.close('all')

See the :meth:`hist <matplotlib.axes.Axes.hist>` method and the
`matplotlib hist documentation <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.hist>`__ for more.


The existing interface ``DataFrame.hist`` to plot histogram still can be used.

.. ipython:: python

   plt.figure();

   @savefig hist_plot_ex.png
   df['A'].diff().hist()

.. ipython:: python
   :suppress:

   plt.close('all')

:meth:`DataFrame.hist` plots the histograms of the columns on multiple
subplots:

.. ipython:: python

   plt.figure()

   @savefig frame_hist_ex.png
   df.diff().hist(color='k', alpha=0.5, bins=50)


.. versionadded:: 0.10.0

The ``by`` keyword can be specified to plot grouped histograms:

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   data = pd.Series(np.random.randn(1000))

   @savefig grouped_hist.png
   data.hist(by=np.random.randint(0, 4, 1000), figsize=(6, 4))


.. _visualization.box:

Box Plots
~~~~~~~~~

.. versionadded:: 0.15.0

Boxplot can be drawn calling :meth:`Series.plot.box` and :meth:`DataFrame.plot.box`,
or :meth:`DataFrame.boxplot` to visualize the distribution of values within each column.

For instance, here is a boxplot representing five trials of 10 observations of
a uniform random variable on [0,1).

.. ipython:: python
   :suppress:

   plt.close('all')
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.rand(10, 5), columns=['A', 'B', 'C', 'D', 'E'])

   @savefig box_plot_new.png
   df.plot.box()

Boxplot can be colorized by passing ``color`` keyword. You can pass a ``dict``
whose keys are ``boxes``, ``whiskers``, ``medians`` and ``caps``.
If some keys are missing in the ``dict``, default colors are used
for the corresponding artists. Also, boxplot has ``sym`` keyword to specify fliers style.

When you pass other type of arguments via ``color`` keyword, it will be directly
passed to matplotlib for all the ``boxes``, ``whiskers``, ``medians`` and ``caps``
colorization.

The colors are applied to every boxes to be drawn. If you want
more complicated colorization, you can get each drawn artists by passing
:ref:`return_type <visualization.box.return>`.

.. ipython:: python

   color = dict(boxes='DarkGreen', whiskers='DarkOrange',
                medians='DarkBlue', caps='Gray')

   @savefig box_new_colorize.png
   df.plot.box(color=color, sym='r+')

.. ipython:: python
   :suppress:

   plt.close('all')

Also, you can pass other keywords supported by matplotlib ``boxplot``.
For example, horizontal and custom-positioned boxplot can be drawn by
``vert=False`` and ``positions`` keywords.

.. ipython:: python

   @savefig box_new_kwargs.png
   df.plot.box(vert=False, positions=[1, 4, 5, 6, 8])


See the :meth:`boxplot <matplotlib.axes.Axes.boxplot>` method and the
`matplotlib boxplot documentation <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.boxplot>`__ for more.


The existing interface ``DataFrame.boxplot`` to plot boxplot still can be used.

.. ipython:: python
   :suppress:

   plt.close('all')
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10,5))
   plt.figure();

   @savefig box_plot_ex.png
   bp = df.boxplot()

You can create a stratified boxplot using the ``by`` keyword argument to create
groupings.  For instance,

.. ipython:: python
   :suppress:

   plt.close('all')
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10,2), columns=['Col1', 'Col2'] )
   df['X'] = pd.Series(['A','A','A','A','A','B','B','B','B','B'])

   plt.figure();

   @savefig box_plot_ex2.png
   bp = df.boxplot(by='X')

You can also pass a subset of columns to plot, as well as group by multiple
columns:

.. ipython:: python
   :suppress:

   plt.close('all')
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10,3), columns=['Col1', 'Col2', 'Col3'])
   df['X'] = pd.Series(['A','A','A','A','A','B','B','B','B','B'])
   df['Y'] = pd.Series(['A','B','A','B','A','B','A','B','A','B'])

   plt.figure();

   @savefig box_plot_ex3.png
   bp = df.boxplot(column=['Col1','Col2'], by=['X','Y'])

.. ipython:: python
   :suppress:

    plt.close('all')

.. _visualization.box.return:

.. warning::

   The default changed from ``'dict'`` to ``'axes'`` in version 0.19.0.

In ``boxplot``, the return type can be controlled by the ``return_type``, keyword. The valid choices are ``{"axes", "dict", "both", None}``.
Faceting, created by ``DataFrame.boxplot`` with the ``by``
keyword, will affect the output type as well:

================ ======= ==========================
``return_type=`` Faceted Output type
---------------- ------- --------------------------

``None``         No      axes
``None``         Yes     2-D ndarray of axes
``'axes'``       No      axes
``'axes'``       Yes     Series of axes
``'dict'``       No      dict of artists
``'dict'``       Yes     Series of dicts of artists
``'both'``       No      namedtuple
``'both'``       Yes     Series of namedtuples
================ ======= ==========================

``Groupby.boxplot`` always returns a Series of ``return_type``.

.. ipython:: python
   :okwarning:

   np.random.seed(1234)
   df_box = pd.DataFrame(np.random.randn(50, 2))
   df_box['g'] = np.random.choice(['A', 'B'], size=50)
   df_box.loc[df_box['g'] == 'B', 1] += 3

   @savefig boxplot_groupby.png
   bp = df_box.boxplot(by='g')

.. ipython:: python
   :suppress:

   plt.close('all')

Compare to:

.. ipython:: python
   :okwarning:

   @savefig groupby_boxplot_vis.png
   bp = df_box.groupby('g').boxplot()

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.area_plot:

Area Plot
~~~~~~~~~

.. versionadded:: 0.14

You can create area plots with :meth:`Series.plot.area` and :meth:`DataFrame.plot.area`.
Area plots are stacked by default. To produce stacked area plot, each column must be either all positive or all negative values.

When input data contains `NaN`, it will be automatically filled by 0. If you want to drop or fill by different values, use :func:`dataframe.dropna` or :func:`dataframe.fillna` before calling `plot`.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.figure()

.. ipython:: python

   df = pd.DataFrame(np.random.rand(10, 4), columns=['a', 'b', 'c', 'd'])

   @savefig area_plot_stacked.png
   df.plot.area();

To produce an unstacked plot, pass ``stacked=False``. Alpha value is set to 0.5 unless otherwise specified:

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()

.. ipython:: python

   @savefig area_plot_unstacked.png
   df.plot.area(stacked=False);

.. _visualization.scatter:

Scatter Plot
~~~~~~~~~~~~

.. versionadded:: 0.13

Scatter plot can be drawn by using the :meth:`DataFrame.plot.scatter` method.
Scatter plot requires numeric columns for x and y axis.
These can be specified by ``x`` and ``y`` keywords each.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.close('all')
   plt.figure()

.. ipython:: python

   df = pd.DataFrame(np.random.rand(50, 4), columns=['a', 'b', 'c', 'd'])

   @savefig scatter_plot.png
   df.plot.scatter(x='a', y='b');

To plot multiple column groups in a single axes, repeat ``plot`` method specifying target ``ax``.
It is recommended to specify ``color`` and ``label`` keywords to distinguish each groups.

.. ipython:: python

   ax = df.plot.scatter(x='a', y='b', color='DarkBlue', label='Group 1');
   @savefig scatter_plot_repeated.png
   df.plot.scatter(x='c', y='d', color='DarkGreen', label='Group 2', ax=ax);

.. ipython:: python
   :suppress:

   plt.close('all')

The keyword ``c`` may be given as the name of a column to provide colors for
each point:

.. ipython:: python

   @savefig scatter_plot_colored.png
   df.plot.scatter(x='a', y='b', c='c', s=50);


.. ipython:: python
   :suppress:

   plt.close('all')

You can pass other keywords supported by matplotlib ``scatter``.
Below example shows a bubble chart using a dataframe column values as bubble size.

.. ipython:: python

   @savefig scatter_plot_bubble.png
   df.plot.scatter(x='a', y='b', s=df['c']*200);

.. ipython:: python
   :suppress:

   plt.close('all')

See the :meth:`scatter <matplotlib.axes.Axes.scatter>` method and the
`matplotlib scatter documentation <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.scatter>`__ for more.

.. _visualization.hexbin:

Hexagonal Bin Plot
~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.14

You can create hexagonal bin plots with :meth:`DataFrame.plot.hexbin`.
Hexbin plots can be a useful alternative to scatter plots if your data are
too dense to plot each point individually.

.. ipython:: python
   :suppress:

   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 2), columns=['a', 'b'])
   df['b'] = df['b'] + np.arange(1000)

   @savefig hexbin_plot.png
   df.plot.hexbin(x='a', y='b', gridsize=25)


A useful keyword argument is ``gridsize``; it controls the number of hexagons
in the x-direction, and defaults to 100. A larger ``gridsize`` means more, smaller
bins.

By default, a histogram of the counts around each ``(x, y)`` point is computed.
You can specify alternative aggregations by passing values to the ``C`` and
``reduce_C_function`` arguments. ``C`` specifies the value at each ``(x, y)`` point
and ``reduce_C_function`` is a function of one argument that reduces all the
values in a bin to a single number (e.g. ``mean``, ``max``, ``sum``, ``std``).  In this
example the positions are given by columns ``a`` and ``b``, while the value is
given by column ``z``. The bins are aggregated with numpy's ``max`` function.

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 2), columns=['a', 'b'])
   df['b'] = df['b'] = df['b'] + np.arange(1000)
   df['z'] = np.random.uniform(0, 3, 1000)

   @savefig hexbin_plot_agg.png
   df.plot.hexbin(x='a', y='b', C='z', reduce_C_function=np.max,
           gridsize=25)

.. ipython:: python
   :suppress:

   plt.close('all')

See the :meth:`hexbin <matplotlib.axes.Axes.hexbin>` method and the
`matplotlib hexbin documentation <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.hexbin>`__ for more.

.. _visualization.pie:

Pie plot
~~~~~~~~

.. versionadded:: 0.14

You can create a pie plot with :meth:`DataFrame.plot.pie` or :meth:`Series.plot.pie`.
If your data includes any ``NaN``, they will be automatically filled with 0.
A ``ValueError`` will be raised if there are any negative values in your data.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.figure()

.. ipython:: python

   series = pd.Series(3 * np.random.rand(4), index=['a', 'b', 'c', 'd'], name='series')

   @savefig series_pie_plot.png
   series.plot.pie(figsize=(6, 6))

.. ipython:: python
   :suppress:

   plt.close('all')

For pie plots it's best to use square figures, one's with an equal aspect ratio. You can create the
figure with equal width and height, or force the aspect ratio to be equal after plotting by
calling ``ax.set_aspect('equal')`` on the returned ``axes`` object.

Note that pie plot with :class:`DataFrame` requires that you either specify a target column by the ``y``
argument or ``subplots=True``. When ``y`` is specified, pie plot of selected column
will be drawn. If ``subplots=True`` is specified, pie plots for each column are drawn as subplots.
A legend will be drawn in each pie plots by default; specify ``legend=False`` to hide it.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.figure()

.. ipython:: python

   df = pd.DataFrame(3 * np.random.rand(4, 2), index=['a', 'b', 'c', 'd'], columns=['x', 'y'])

   @savefig df_pie_plot.png
   df.plot.pie(subplots=True, figsize=(8, 4))

.. ipython:: python
   :suppress:

   plt.close('all')

You can use the ``labels`` and ``colors`` keywords to specify the labels and colors of each wedge.

.. warning::

   Most pandas plots use the the ``label`` and ``color`` arguments (note the lack of "s" on those).
   To be consistent with :func:`matplotlib.pyplot.pie` you must use ``labels`` and ``colors``.

If you want to hide wedge labels, specify ``labels=None``.
If ``fontsize`` is specified, the value will be applied to wedge labels.
Also, other keywords supported by :func:`matplotlib.pyplot.pie` can be used.


.. ipython:: python
   :suppress:

   plt.figure()

.. ipython:: python

   @savefig series_pie_plot_options.png
   series.plot.pie(labels=['AA', 'BB', 'CC', 'DD'], colors=['r', 'g', 'b', 'c'],
                   autopct='%.2f', fontsize=20, figsize=(6, 6))

If you pass values whose sum total is less than 1.0, matplotlib draws a semicircle.

.. ipython:: python
   :suppress:

   plt.close('all')
   plt.figure()

.. ipython:: python

   series = pd.Series([0.1] * 4, index=['a', 'b', 'c', 'd'], name='series2')

   @savefig series_pie_plot_semi.png
   series.plot.pie(figsize=(6, 6))

See the `matplotlib pie documentation <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.pie>`__ for more.

.. ipython:: python
    :suppress:

    plt.close('all')

.. _visualization.missing_data:

Plotting with Missing Data
--------------------------

Pandas tries to be pragmatic about plotting DataFrames or Series
that contain missing data. Missing values are dropped, left out, or filled
depending on the plot type.

+----------------+--------------------------------------+
| Plot Type      | NaN Handling                         |
+================+======================================+
| Line           | Leave gaps at NaNs                   |
+----------------+--------------------------------------+
| Line (stacked) | Fill 0's                             |
+----------------+--------------------------------------+
| Bar            | Fill 0's                             |
+----------------+--------------------------------------+
| Scatter        | Drop NaNs                            |
+----------------+--------------------------------------+
| Histogram      | Drop NaNs (column-wise)              |
+----------------+--------------------------------------+
| Box            | Drop NaNs (column-wise)              |
+----------------+--------------------------------------+
| Area           | Fill 0's                             |
+----------------+--------------------------------------+
| KDE            | Drop NaNs (column-wise)              |
+----------------+--------------------------------------+
| Hexbin         | Drop NaNs                            |
+----------------+--------------------------------------+
| Pie            | Fill 0's                             |
+----------------+--------------------------------------+

If any of these defaults are not what you want, or if you want to be
explicit about how missing values are handled, consider using
:meth:`~pandas.DataFrame.fillna` or :meth:`~pandas.DataFrame.dropna`
before plotting.

.. _visualization.tools:

Plotting Tools
--------------

These functions can be imported from ``pandas.tools.plotting``
and take a :class:`Series` or :class:`DataFrame` as an argument.

.. _visualization.scatter_matrix:

Scatter Matrix Plot
~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.7.3

You can create a scatter plot matrix using the
``scatter_matrix`` method in ``pandas.tools.plotting``:

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.tools.plotting import scatter_matrix
   df = pd.DataFrame(np.random.randn(1000, 4), columns=['a', 'b', 'c', 'd'])

   @savefig scatter_matrix_kde.png
   scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal='kde')

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.kde:

Density Plot
~~~~~~~~~~~~

.. versionadded:: 0.8.0

You can create density plots using the :meth:`Series.plot.kde` and :meth:`DataFrame.plot.kde` methods.

.. ipython:: python
   :suppress:

   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   ser = pd.Series(np.random.randn(1000))

   @savefig kde_plot.png
   ser.plot.kde()

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.andrews_curves:

Andrews Curves
~~~~~~~~~~~~~~

Andrews curves allow one to plot multivariate data as a large number
of curves that are created using the attributes of samples as coefficients
for Fourier series. By coloring these curves differently for each class
it is possible to visualize data clustering. Curves belonging to samples
of the same class will usually be closer together and form larger structures.

**Note**: The "Iris" dataset is available `here <https://raw.github.com/pandas-dev/pandas/master/pandas/tests/data/iris.csv>`__.

.. ipython:: python

   from pandas.tools.plotting import andrews_curves

   data = pd.read_csv('data/iris.data')

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

   from pandas.tools.plotting import parallel_coordinates

   data = pd.read_csv('data/iris.data')

   plt.figure()

   @savefig parallel_coordinates.png
   parallel_coordinates(data, 'Name')

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.lag:

Lag Plot
~~~~~~~~

Lag plots are used to check if a data set or time series is random. Random
data should not exhibit any structure in the lag plot. Non-random structure
implies that the underlying data are not random.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.tools.plotting import lag_plot

   plt.figure()

   data = pd.Series(0.1 * np.random.rand(1000) +
       0.9 * np.sin(np.linspace(-99 * np.pi, 99 * np.pi, num=1000)))

   @savefig lag_plot.png
   lag_plot(data)

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.autocorrelation:

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
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.tools.plotting import autocorrelation_plot

   plt.figure()

   data = pd.Series(0.7 * np.random.rand(1000) +
      0.3 * np.sin(np.linspace(-9 * np.pi, 9 * np.pi, num=1000)))

   @savefig autocorrelation_plot.png
   autocorrelation_plot(data)

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.bootstrap:

Bootstrap Plot
~~~~~~~~~~~~~~

Bootstrap plots are used to visually assess the uncertainty of a statistic, such
as mean, median, midrange, etc. A random subset of a specified size is selected
from a data set, the statistic in question is computed for this subset and the
process is repeated a specified number of times. Resulting plots and histograms
are what constitutes the bootstrap plot.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.tools.plotting import bootstrap_plot

   data = pd.Series(np.random.rand(1000))

   @savefig bootstrap_plot.png
   bootstrap_plot(data, size=50, samples=500, color='grey')

.. ipython:: python
   :suppress:

    plt.close('all')

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

**Note**: The "Iris" dataset is available `here <https://raw.github.com/pandas-dev/pandas/master/pandas/tests/data/iris.csv>`__.

.. ipython:: python

   from pandas.tools.plotting import radviz

   data = pd.read_csv('data/iris.data')

   plt.figure()

   @savefig radviz.png
   radviz(data, 'Name')

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.formatting:

Plot Formatting
---------------

Most plotting methods have a set of keyword arguments that control the
layout and formatting of the returned plot:

.. ipython:: python

   @savefig series_plot_basic2.png
   plt.figure(); ts.plot(style='k--', label='Series');

.. ipython:: python
   :suppress:

   plt.close('all')

For each kind of plot (e.g. `line`, `bar`, `scatter`) any additional arguments
keywords are passed along to the corresponding matplotlib function
(:meth:`ax.plot() <matplotlib.axes.Axes.plot>`,
:meth:`ax.bar() <matplotlib.axes.Axes.bar>`,
:meth:`ax.scatter() <matplotlib.axes.Axes.scatter>`). These can be used
to control additional styling, beyond what pandas provides.

Controlling the Legend
~~~~~~~~~~~~~~~~~~~~~~

You may set the ``legend`` argument to ``False`` to hide the legend, which is
shown by default.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
   df = df.cumsum()

   @savefig frame_plot_basic_noleg.png
   df.plot(legend=False)

.. ipython:: python
   :suppress:

   plt.close('all')

Scales
~~~~~~

You may pass ``logy`` to get a log-scale Y axis.

.. ipython:: python
   :suppress:

   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
   ts = np.exp(ts.cumsum())

   @savefig series_plot_logy.png
   ts.plot(logy=True)

.. ipython:: python
   :suppress:

   plt.close('all')

See also the ``logx`` and ``loglog`` keyword arguments.

Plotting on a Secondary Y-axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To plot data on a secondary y-axis, use the ``secondary_y`` keyword:

.. ipython:: python
   :suppress:

   plt.figure()

.. ipython:: python

   df.A.plot()

   @savefig series_plot_secondary_y.png
   df.B.plot(secondary_y=True, style='g')

.. ipython:: python
   :suppress:

   plt.close('all')

To plot some columns in a DataFrame, give the column names to the ``secondary_y``
keyword:

.. ipython:: python

   plt.figure()
   ax = df.plot(secondary_y=['A', 'B'])
   ax.set_ylabel('CD scale')
   @savefig frame_plot_secondary_y.png
   ax.right_ax.set_ylabel('AB scale')

.. ipython:: python
   :suppress:

   plt.close('all')

Note that the columns plotted on the secondary y-axis is automatically marked
with "(right)" in the legend. To turn off the automatic marking, use the
``mark_right=False`` keyword:

.. ipython:: python

   plt.figure()

   @savefig frame_plot_secondary_y_no_right.png
   df.plot(secondary_y=['A', 'B'], mark_right=False)

.. ipython:: python
   :suppress:

   plt.close('all')

Suppressing Tick Resolution Adjustment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pandas includes automatic tick resolution adjustment for regular frequency
time-series data. For limited cases where pandas cannot infer the frequency
information (e.g., in an externally created ``twinx``), you can choose to
suppress this behavior for alignment purposes.

Here is the default behavior, notice how the x-axis tick labelling is performed:

.. ipython:: python

   plt.figure()

   @savefig ser_plot_suppress.png
   df.A.plot()

.. ipython:: python
   :suppress:

   plt.close('all')

Using the ``x_compat`` parameter, you can suppress this behavior:

.. ipython:: python

   plt.figure()

   @savefig ser_plot_suppress_parm.png
   df.A.plot(x_compat=True)

.. ipython:: python
   :suppress:

   plt.close('all')

If you have more than one plot that needs to be suppressed, the ``use`` method
in ``pandas.plot_params`` can be used in a `with statement`:

.. ipython:: python

   plt.figure()

   @savefig ser_plot_suppress_context.png
   with pd.plot_params.use('x_compat', True):
       df.A.plot(color='r')
       df.B.plot(color='g')
       df.C.plot(color='b')

.. ipython:: python
   :suppress:

   plt.close('all')

Subplots
~~~~~~~~

Each Series in a DataFrame can be plotted on a different axis
with the ``subplots`` keyword:

.. ipython:: python

   @savefig frame_plot_subplots.png
   df.plot(subplots=True, figsize=(6, 6));

.. ipython:: python
   :suppress:

   plt.close('all')

Using Layout and Targeting Multiple Axes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The layout of subplots can be specified by ``layout`` keyword. It can accept
``(rows, columns)``. The ``layout`` keyword can be used in
``hist`` and ``boxplot`` also. If input is invalid, ``ValueError`` will be raised.

The number of axes which can be contained by rows x columns specified by ``layout`` must be
larger than the number of required subplots. If layout can contain more axes than required,
blank axes are not drawn. Similar to a numpy array's ``reshape`` method, you
can use ``-1`` for one dimension to automatically calculate the number of rows
or columns needed, given the other.

.. ipython:: python

   @savefig frame_plot_subplots_layout.png
   df.plot(subplots=True, layout=(2, 3), figsize=(6, 6), sharex=False);

.. ipython:: python
   :suppress:

   plt.close('all')

The above example is identical to using

.. ipython:: python

   df.plot(subplots=True, layout=(2, -1), figsize=(6, 6), sharex=False);

.. ipython:: python
   :suppress:

   plt.close('all')

The required number of columns (3) is inferred from the number of series to plot
and the given number of rows (2).

Also, you can pass multiple axes created beforehand as list-like via ``ax`` keyword.
This allows to use more complicated layout.
The passed axes must be the same number as the subplots being drawn.

When multiple axes are passed via ``ax`` keyword, ``layout``, ``sharex`` and ``sharey`` keywords
don't affect to the output. You should explicitly pass ``sharex=False`` and ``sharey=False``,
otherwise you will see a warning.

.. ipython:: python

   fig, axes = plt.subplots(4, 4, figsize=(6, 6));
   plt.subplots_adjust(wspace=0.5, hspace=0.5);
   target1 = [axes[0][0], axes[1][1], axes[2][2], axes[3][3]]
   target2 = [axes[3][0], axes[2][1], axes[1][2], axes[0][3]]

   df.plot(subplots=True, ax=target1, legend=False, sharex=False, sharey=False);
   @savefig frame_plot_subplots_multi_ax.png
   (-df).plot(subplots=True, ax=target2, legend=False, sharex=False, sharey=False);

.. ipython:: python
   :suppress:

   plt.close('all')

Another option is passing an ``ax`` argument to :meth:`Series.plot` to plot on a particular axis:

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
   ts = ts.cumsum()

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list('ABCD'))
   df = df.cumsum()

.. ipython:: python
   :suppress:

   plt.close('all')

.. ipython:: python

   fig, axes = plt.subplots(nrows=2, ncols=2)
   df['A'].plot(ax=axes[0,0]); axes[0,0].set_title('A');
   df['B'].plot(ax=axes[0,1]); axes[0,1].set_title('B');
   df['C'].plot(ax=axes[1,0]); axes[1,0].set_title('C');

   @savefig series_plot_multi.png
   df['D'].plot(ax=axes[1,1]); axes[1,1].set_title('D');

.. ipython:: python
   :suppress:

    plt.close('all')

.. _visualization.errorbars:

Plotting With Error Bars
~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.14

Plotting with error bars is now supported in the :meth:`DataFrame.plot` and :meth:`Series.plot`

Horizontal and vertical errorbars can be supplied to the ``xerr`` and ``yerr`` keyword arguments to :meth:`~DataFrame.plot()`. The error values can be specified using a variety of formats.

- As a :class:`DataFrame` or ``dict`` of errors with column names matching the ``columns`` attribute of the plotting :class:`DataFrame` or matching the ``name`` attribute of the :class:`Series`
- As a ``str`` indicating which of the columns of plotting :class:`DataFrame` contain the error values
- As raw values (``list``, ``tuple``, or ``np.ndarray``). Must be the same length as the plotting :class:`DataFrame`/:class:`Series`

Asymmetrical error bars are also supported, however raw error values must be provided in this case. For a ``M`` length :class:`Series`, a ``Mx2`` array should be provided indicating lower and upper (or left and right) errors. For a ``MxN`` :class:`DataFrame`, asymmetrical errors should be in a ``Mx2xN`` array.

Here is an example of one way to easily plot group means with standard deviations from the raw data.

.. ipython:: python

   # Generate the data
   ix3 = pd.MultiIndex.from_arrays([['a', 'a', 'a', 'a', 'b', 'b', 'b', 'b'], ['foo', 'foo', 'bar', 'bar', 'foo', 'foo', 'bar', 'bar']], names=['letter', 'word'])
   df3 = pd.DataFrame({'data1': [3, 2, 4, 3, 2, 4, 3, 2], 'data2': [6, 5, 7, 5, 4, 5, 6, 5]}, index=ix3)

   # Group by index labels and take the means and standard deviations for each group
   gp3 = df3.groupby(level=('letter', 'word'))
   means = gp3.mean()
   errors = gp3.std()
   means
   errors

   # Plot
   fig, ax = plt.subplots()
   @savefig errorbar_example.png
   means.plot.bar(yerr=errors, ax=ax)

.. ipython:: python
   :suppress:

   plt.close('all')

.. _visualization.table:

Plotting Tables
~~~~~~~~~~~~~~~

.. versionadded:: 0.14

Plotting with matplotlib table is now supported in  :meth:`DataFrame.plot` and :meth:`Series.plot` with a ``table`` keyword. The ``table`` keyword can accept ``bool``, :class:`DataFrame` or :class:`Series`. The simple way to draw a table is to specify ``table=True``. Data will be transposed to meet matplotlib's default layout.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   fig, ax = plt.subplots(1, 1)
   df = pd.DataFrame(np.random.rand(5, 3), columns=['a', 'b', 'c'])
   ax.get_xaxis().set_visible(False)   # Hide Ticks

   @savefig line_plot_table_true.png
   df.plot(table=True, ax=ax)

.. ipython:: python
   :suppress:

   plt.close('all')

Also, you can pass different :class:`DataFrame` or :class:`Series` for ``table`` keyword. The data will be drawn as displayed in print method (not transposed automatically). If required, it should be transposed manually as below example.

.. ipython:: python

   fig, ax = plt.subplots(1, 1)
   ax.get_xaxis().set_visible(False)   # Hide Ticks
   @savefig line_plot_table_data.png
   df.plot(table=np.round(df.T, 2), ax=ax)

.. ipython:: python
   :suppress:

   plt.close('all')

Finally, there is a helper function ``pandas.tools.plotting.table`` to create a table from :class:`DataFrame` and :class:`Series`, and add it to an ``matplotlib.Axes``. This function can accept keywords which matplotlib table has.

.. ipython:: python

   from pandas.tools.plotting import table
   fig, ax = plt.subplots(1, 1)

   table(ax, np.round(df.describe(), 2),
         loc='upper right', colWidths=[0.2, 0.2, 0.2])

   @savefig line_plot_table_describe.png
   df.plot(ax=ax, ylim=(0, 2), legend=None)

.. ipython:: python
   :suppress:

   plt.close('all')

**Note**: You can get table instances on the axes using ``axes.tables`` property for further decorations. See the `matplotlib table documentation <http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.table>`__ for more.

.. _visualization.colormaps:

Colormaps
~~~~~~~~~

A potential issue when plotting a large number of columns is that it can be
difficult to distinguish some series due to repetition in the default colors. To
remedy this, DataFrame plotting supports the use of the ``colormap=`` argument,
which accepts either a Matplotlib `colormap <http://matplotlib.org/api/cm_api.html>`__
or a string that is a name of a colormap registered with Matplotlib. A
visualization of the default matplotlib colormaps is available `here
<http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps>`__.

As matplotlib does not directly support colormaps for line-based plots, the
colors are selected based on an even spacing determined by the number of columns
in the DataFrame. There is no consideration made for background color, so some
colormaps will produce lines that are not easily visible.

To use the cubehelix colormap, we can simply pass ``'cubehelix'`` to ``colormap=``

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 10), index=ts.index)
   df = df.cumsum()

   plt.figure()

   @savefig cubehelix.png
   df.plot(colormap='cubehelix')

.. ipython:: python
   :suppress:

   plt.close('all')

or we can pass the colormap itself

.. ipython:: python

   from matplotlib import cm

   plt.figure()

   @savefig cubehelix_cm.png
   df.plot(colormap=cm.cubehelix)

.. ipython:: python
   :suppress:

   plt.close('all')

Colormaps can also be used other plot types, like bar charts:

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   dd = pd.DataFrame(np.random.randn(10, 10)).applymap(abs)
   dd = dd.cumsum()

   plt.figure()

   @savefig greens.png
   dd.plot.bar(colormap='Greens')

.. ipython:: python
   :suppress:

   plt.close('all')

Parallel coordinates charts:

.. ipython:: python

   plt.figure()

   @savefig parallel_gist_rainbow.png
   parallel_coordinates(data, 'Name', colormap='gist_rainbow')

.. ipython:: python
   :suppress:

   plt.close('all')

Andrews curves charts:

.. ipython:: python

   plt.figure()

   @savefig andrews_curve_winter.png
   andrews_curves(data, 'Name', colormap='winter')

.. ipython:: python
   :suppress:

   plt.close('all')

Plotting directly with matplotlib
---------------------------------

In some situations it may still be preferable or necessary to prepare plots
directly with matplotlib, for instance when a certain type of plot or
customization is not (yet) supported by pandas. Series and DataFrame objects
behave like arrays and can therefore be passed directly to matplotlib functions
without explicit casts.

pandas also automatically registers formatters and locators that recognize date
indices, thereby extending date and time support to practically all plot types
available in matplotlib. Although this formatting does not provide the same
level of refinement you would get when plotting via pandas, it can be faster
when plotting a large number of points.

.. note::

    The speed up for large data sets only applies to pandas 0.14.0 and later.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   price = pd.Series(np.random.randn(150).cumsum(),
                     index=pd.date_range('2000-1-1', periods=150, freq='B'))
   ma = price.rolling(20).mean()
   mstd = price.rolling(20).std()

   plt.figure()

   plt.plot(price.index, price, 'k')
   plt.plot(ma.index, ma, 'b')
   @savefig bollinger.png
   plt.fill_between(mstd.index, ma-2*mstd, ma+2*mstd, color='b', alpha=0.2)

.. ipython:: python
   :suppress:

    plt.close('all')


.. _rplot:


Trellis plotting interface
--------------------------

.. warning::

    The ``rplot`` trellis plotting interface has been **removed**. Please use
    external packages like `seaborn <https://github.com/mwaskom/seaborn>`_ for
    similar but more refined functionality and refer to our 0.18.1 documentation
    `here <http://pandas.pydata.org/pandas-docs/version/0.18.1/visualization.html>`__
    for how to convert to using it.
