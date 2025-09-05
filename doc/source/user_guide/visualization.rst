.. _visualization:

{{ header }}

*******************
Chart visualization
*******************


.. note::

   The examples below assume that you're using `Jupyter <https://jupyter.org/>`_.

This section demonstrates visualization through charting. For information on
visualization of tabular data please see the section on `Table Visualization <style.ipynb>`_.

We use the standard convention for referencing the matplotlib API:

.. ipython:: python

   import matplotlib.pyplot as plt

   plt.close("all")

We provide the basics in pandas to easily create decent looking plots.
See `the ecosystem page <https://pandas.pydata.org/community/ecosystem.html>`_ for visualization
libraries that go beyond the basics documented here.

.. note::

   All calls to ``np.random`` are seeded with 123456.

.. _visualization.basic:

Basic plotting: ``plot``
------------------------

We will demonstrate the basics, see the :ref:`cookbook<cookbook.plotting>` for
some advanced strategies.

The ``plot`` method on Series and DataFrame is just a simple wrapper around
:meth:`plt.plot() <matplotlib.axes.Axes.plot>`:

.. ipython:: python

   np.random.seed(123456)

   ts = pd.Series(np.random.randn(1000), index=pd.date_range("1/1/2000", periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png
   ts.plot();

If the index consists of dates, it calls :meth:`gcf().autofmt_xdate() <matplotlib.figure.Figure.autofmt_xdate>`
to try to format the x-axis nicely as per above.

On DataFrame, :meth:`~DataFrame.plot` is a convenience to plot all of the columns with labels:

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list("ABCD"))
   df = df.cumsum()

   plt.figure();
   @savefig frame_plot_basic.png
   df.plot();

You can plot one column versus another using the ``x`` and ``y`` keywords in
:meth:`~DataFrame.plot`:

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df3 = pd.DataFrame(np.random.randn(1000, 2), columns=["B", "C"]).cumsum()
   df3["A"] = pd.Series(list(range(len(df))))

   @savefig df_plot_xy.png
   df3.plot(x="A", y="B");

.. note::

   For more formatting and styling options, see
   :ref:`formatting <visualization.formatting>` below.

.. ipython:: python
    :suppress:

    plt.close("all")

.. _visualization.other:

Other plots
-----------

Plotting methods allow for a handful of plot styles other than the
default line plot. These methods can be provided as the ``kind``
keyword argument to :meth:`~DataFrame.plot`, and include:

* :ref:`'bar' <visualization.barplot>` or :ref:`'barh' <visualization.barplot>` for bar plots
* :ref:`'hist' <visualization.hist>` for histogram
* :ref:`'box' <visualization.box>` for boxplot
* :ref:`'kde' <visualization.kde>` or :ref:`'density' <visualization.kde>` for density plots
* :ref:`'area' <visualization.area_plot>` for area plots
* :ref:`'scatter' <visualization.scatter>` for scatter plots
* :ref:`'hexbin' <visualization.hexbin>` for hexagonal bin plots
* :ref:`'pie' <visualization.pie>` for pie plots

For example, a bar plot can be created the following way:

.. ipython:: python

   plt.figure();

   @savefig bar_plot_ex.png
   df.iloc[5].plot(kind="bar");

You can also create these other plots using the methods ``DataFrame.plot.<kind>`` instead of providing the ``kind`` keyword argument. This makes it easier to discover plot methods and the specific arguments they use:

.. ipython::
    :verbatim:

    In [14]: df = pd.DataFrame()

    In [15]: df.plot.<TAB>  # noqa: E225, E999
    df.plot.area     df.plot.barh     df.plot.density  df.plot.hist     df.plot.line     df.plot.scatter
    df.plot.bar      df.plot.box      df.plot.hexbin   df.plot.kde      df.plot.pie

In addition to these ``kind`` s, there are the :ref:`DataFrame.hist() <visualization.hist>`,
and :ref:`DataFrame.boxplot() <visualization.box>` methods, which use a separate interface.

Finally, there are several :ref:`plotting functions <visualization.tools>` in ``pandas.plotting``
that take a :class:`Series` or :class:`DataFrame` as an argument. These
include:

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
   df.iloc[5].plot.bar();
   plt.axhline(0, color="k");

Calling a DataFrame's :meth:`plot.bar() <DataFrame.plot.bar>` method produces a multiple
bar plot:

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df2 = pd.DataFrame(np.random.rand(10, 4), columns=["a", "b", "c", "d"])

   @savefig bar_plot_multi_ex.png
   df2.plot.bar();

To produce a stacked bar plot, pass ``stacked=True``:

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()

.. ipython:: python

   @savefig bar_plot_stacked_ex.png
   df2.plot.bar(stacked=True);

To get horizontal bar plots, use the ``barh`` method:

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()

.. ipython:: python

   @savefig barh_plot_stacked_ex.png
   df2.plot.barh(stacked=True);

.. _visualization.hist:

Histograms
~~~~~~~~~~

Histograms can be drawn by using the :meth:`DataFrame.plot.hist` and :meth:`Series.plot.hist` methods.

.. ipython:: python

   df4 = pd.DataFrame(
       {
           "a": np.random.randn(1000) + 1,
           "b": np.random.randn(1000),
           "c": np.random.randn(1000) - 1,
       },
       columns=["a", "b", "c"],
   )

   plt.figure();

   @savefig hist_new.png
   df4.plot.hist(alpha=0.5);


.. ipython:: python
   :suppress:

   plt.close("all")

A histogram can be stacked using ``stacked=True``. Bin size can be changed
using the ``bins`` keyword.

.. ipython:: python

   plt.figure();

   @savefig hist_new_stacked.png
   df4.plot.hist(stacked=True, bins=20);

.. ipython:: python
   :suppress:

   plt.close("all")

You can pass other keywords supported by matplotlib ``hist``. For example,
horizontal and cumulative histograms can be drawn by
``orientation='horizontal'`` and ``cumulative=True``.

.. ipython:: python

   plt.figure();

   @savefig hist_new_kwargs.png
   df4["a"].plot.hist(orientation="horizontal", cumulative=True);

.. ipython:: python
   :suppress:

   plt.close("all")

See the :meth:`hist <matplotlib.axes.Axes.hist>` method and the
`matplotlib hist documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html>`__ for more.


The existing interface ``DataFrame.hist`` to plot histogram still can be used.

.. ipython:: python

   plt.figure();

   @savefig hist_plot_ex.png
   df["A"].diff().hist();

.. ipython:: python
   :suppress:

   plt.close("all")

:meth:`DataFrame.hist` plots the histograms of the columns on multiple
subplots:

.. ipython:: python

   plt.figure();

   @savefig frame_hist_ex.png
   df.diff().hist(color="k", alpha=0.5, bins=50);


The ``by`` keyword can be specified to plot grouped histograms:

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   data = pd.Series(np.random.randn(1000))

   @savefig grouped_hist.png
   data.hist(by=np.random.randint(0, 4, 1000), figsize=(6, 4));

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

In addition, the ``by`` keyword can also be specified in :meth:`DataFrame.plot.hist`.

.. versionchanged:: 1.4.0

.. ipython:: python

   data = pd.DataFrame(
       {
           "a": np.random.choice(["x", "y", "z"], 1000),
           "b": np.random.choice(["e", "f", "g"], 1000),
           "c": np.random.randn(1000),
           "d": np.random.randn(1000) - 1,
       },
   )

   @savefig grouped_hist_by.png
   data.plot.hist(by=["a", "b"], figsize=(10, 5));

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.box:

Box plots
~~~~~~~~~

Boxplot can be drawn calling :meth:`Series.plot.box` and :meth:`DataFrame.plot.box`,
or :meth:`DataFrame.boxplot` to visualize the distribution of values within each column.

For instance, here is a boxplot representing five trials of 10 observations of
a uniform random variable on [0,1).

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.rand(10, 5), columns=["A", "B", "C", "D", "E"])

   @savefig box_plot_new.png
   df.plot.box();

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

   color = {
       "boxes": "DarkGreen",
       "whiskers": "DarkOrange",
       "medians": "DarkBlue",
       "caps": "Gray",
   }

   @savefig box_new_colorize.png
   df.plot.box(color=color, sym="r+");

.. ipython:: python
   :suppress:

   plt.close("all")

Also, you can pass other keywords supported by matplotlib ``boxplot``.
For example, horizontal and custom-positioned boxplot can be drawn by
``vert=False`` and ``positions`` keywords.

.. ipython:: python

   @savefig box_new_kwargs.png
   df.plot.box(vert=False, positions=[1, 4, 5, 6, 8]);


See the :meth:`boxplot <matplotlib.axes.Axes.boxplot>` method and the
`matplotlib boxplot documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html>`__ for more.


The existing interface ``DataFrame.boxplot`` to plot boxplot still can be used.

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10, 5))
   plt.figure();

   @savefig box_plot_ex.png
   bp = df.boxplot()

You can create a stratified boxplot using the ``by`` keyword argument to create
groupings.  For instance,

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10, 2), columns=["Col1", "Col2"])
   df["X"] = pd.Series(["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"])

   plt.figure();

   @savefig box_plot_ex2.png
   bp = df.boxplot(by="X")

You can also pass a subset of columns to plot, as well as group by multiple
columns:

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10, 3), columns=["Col1", "Col2", "Col3"])
   df["X"] = pd.Series(["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"])
   df["Y"] = pd.Series(["A", "B", "A", "B", "A", "B", "A", "B", "A", "B"])

   plt.figure();

   @savefig box_plot_ex3.png
   bp = df.boxplot(column=["Col1", "Col2"], by=["X", "Y"])

.. ipython:: python
   :suppress:

    plt.close("all")

You could also create groupings with :meth:`DataFrame.plot.box`, for instance:

.. versionchanged:: 1.4.0

.. ipython:: python
   :suppress:

   plt.close("all")
   np.random.seed(123456)

.. ipython:: python
   :okwarning:

   df = pd.DataFrame(np.random.rand(10, 3), columns=["Col1", "Col2", "Col3"])
   df["X"] = pd.Series(["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"])

   plt.figure();

   @savefig box_plot_ex4.png
   bp = df.plot.box(column=["Col1", "Col2"], by="X")

.. ipython:: python
   :suppress:

    plt.close("all")

.. _visualization.box.return:

In ``boxplot``, the return type can be controlled by the ``return_type``, keyword. The valid choices are ``{"axes", "dict", "both", None}``.
Faceting, created by ``DataFrame.boxplot`` with the ``by``
keyword, will affect the output type as well:

================ ======= ==========================
``return_type``  Faceted Output type
================ ======= ==========================
``None``         No      axes
``None``         Yes     2-D ndarray of axes
``'axes'``       No      axes
``'axes'``       Yes     Series of axes
``'dict'``       No      dict of artists
``'dict'``       Yes     Series of dicts of artists
``'both'``       No      namedtuple
``'both'``       Yes     Series of namedtuples
================ ======= ==========================

``Groupby.boxplot`` always returns a ``Series`` of ``return_type``.

.. ipython:: python
   :okwarning:

   np.random.seed(1234)
   df_box = pd.DataFrame(np.random.randn(50, 2))
   df_box["g"] = np.random.choice(["A", "B"], size=50)
   df_box.loc[df_box["g"] == "B", 1] += 3

   @savefig boxplot_groupby.png
   bp = df_box.boxplot(by="g")

.. ipython:: python
   :suppress:

   plt.close("all")

The subplots above are split by the numeric columns first, then the value of
the ``g`` column. Below the subplots are first split by the value of ``g``,
then by the numeric columns.

.. ipython:: python
   :okwarning:

   @savefig groupby_boxplot_vis.png
   bp = df_box.groupby("g").boxplot()

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.area_plot:

Area plot
~~~~~~~~~

You can create area plots with :meth:`Series.plot.area` and :meth:`DataFrame.plot.area`.
Area plots are stacked by default. To produce stacked area plot, each column must be either all positive or all negative values.

When input data contains ``NaN``, it will be automatically filled by 0. If you want to drop or fill by different values, use :func:`dataframe.dropna` or :func:`dataframe.fillna` before calling ``plot``.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.figure()

.. ipython:: python

   df = pd.DataFrame(np.random.rand(10, 4), columns=["a", "b", "c", "d"])

   @savefig area_plot_stacked.png
   df.plot.area();

To produce an unstacked plot, pass ``stacked=False``. Alpha value is set to 0.5 unless otherwise specified:

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()

.. ipython:: python

   @savefig area_plot_unstacked.png
   df.plot.area(stacked=False);

.. _visualization.scatter:

Scatter plot
~~~~~~~~~~~~

Scatter plot can be drawn by using the :meth:`DataFrame.plot.scatter` method.
Scatter plot requires numeric columns for the x and y axes.
These can be specified by the ``x`` and ``y`` keywords.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.close("all")
   plt.figure()

.. ipython:: python

   df = pd.DataFrame(np.random.rand(50, 4), columns=["a", "b", "c", "d"])
   df["species"] = pd.Categorical(
       ["setosa"] * 20 + ["versicolor"] * 20 + ["virginica"] * 10
   )

   @savefig scatter_plot.png
   df.plot.scatter(x="a", y="b");

To plot multiple column groups in a single axes, repeat ``plot`` method specifying target ``ax``.
It is recommended to specify ``color`` and ``label`` keywords to distinguish each groups.

.. ipython:: python
   :okwarning:

   ax = df.plot.scatter(x="a", y="b", color="DarkBlue", label="Group 1")
   @savefig scatter_plot_repeated.png
   df.plot.scatter(x="c", y="d", color="DarkGreen", label="Group 2", ax=ax);

.. ipython:: python
   :suppress:

   plt.close("all")

The keyword ``c`` may be given as the name of a column to provide colors for
each point:

.. ipython:: python

   @savefig scatter_plot_colored.png
   df.plot.scatter(x="a", y="b", c="c", s=50);


.. ipython:: python
   :suppress:

   plt.close("all")

If a categorical column is passed to ``c``, then a discrete colorbar will be produced:

.. versionadded:: 1.3.0

.. ipython:: python

   @savefig scatter_plot_categorical.png
   df.plot.scatter(x="a", y="b", c="species", cmap="viridis", s=50);


.. ipython:: python
   :suppress:

   plt.close("all")

You can pass other keywords supported by matplotlib
:meth:`scatter <matplotlib.axes.Axes.scatter>`. The example  below shows a
bubble chart using a column of the ``DataFrame`` as the bubble size.

.. ipython:: python

   @savefig scatter_plot_bubble.png
   df.plot.scatter(x="a", y="b", s=df["c"] * 200);

.. ipython:: python
   :suppress:

   plt.close("all")

See the :meth:`scatter <matplotlib.axes.Axes.scatter>` method and the
`matplotlib scatter documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html>`__ for more.

.. _visualization.hexbin:

Hexagonal bin plot
~~~~~~~~~~~~~~~~~~

You can create hexagonal bin plots with :meth:`DataFrame.plot.hexbin`.
Hexbin plots can be a useful alternative to scatter plots if your data are
too dense to plot each point individually.

.. ipython:: python
   :suppress:

   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 2), columns=["a", "b"])
   df["b"] = df["b"] + np.arange(1000)

   @savefig hexbin_plot.png
   df.plot.hexbin(x="a", y="b", gridsize=25);


A useful keyword argument is ``gridsize``; it controls the number of hexagons
in the x-direction, and defaults to 100. A larger ``gridsize`` means more, smaller
bins.

By default, a histogram of the counts around each ``(x, y)`` point is computed.
You can specify alternative aggregations by passing values to the ``C`` and
``reduce_C_function`` arguments. ``C`` specifies the value at each ``(x, y)`` point
and ``reduce_C_function`` is a function of one argument that reduces all the
values in a bin to a single number (e.g. ``mean``, ``max``, ``sum``, ``std``).  In this
example the positions are given by columns ``a`` and ``b``, while the value is
given by column ``z``. The bins are aggregated with NumPy's ``max`` function.

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 2), columns=["a", "b"])
   df["b"] = df["b"] + np.arange(1000)
   df["z"] = np.random.uniform(0, 3, 1000)

   @savefig hexbin_plot_agg.png
   df.plot.hexbin(x="a", y="b", C="z", reduce_C_function=np.max, gridsize=25);

.. ipython:: python
   :suppress:

   plt.close("all")

See the :meth:`hexbin <matplotlib.axes.Axes.hexbin>` method and the
`matplotlib hexbin documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hexbin.html>`__ for more.

.. _visualization.pie:

Pie plot
~~~~~~~~

You can create a pie plot with :meth:`DataFrame.plot.pie` or :meth:`Series.plot.pie`.
If your data includes any ``NaN``, they will be automatically filled with 0.
A ``ValueError`` will be raised if there are any negative values in your data.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.figure()

.. ipython:: python
   :okwarning:

   series = pd.Series(3 * np.random.rand(4), index=["a", "b", "c", "d"], name="series")

   @savefig series_pie_plot.png
   series.plot.pie(figsize=(6, 6));

.. ipython:: python
   :suppress:

   plt.close("all")

For pie plots it's best to use square figures, i.e. a figure aspect ratio 1.
You can create the figure with equal width and height, or force the aspect ratio
to be equal after plotting by calling ``ax.set_aspect('equal')`` on the returned
``axes`` object.

Note that pie plot with :class:`DataFrame` requires that you either specify a
target column by the ``y`` argument or ``subplots=True``. When ``y`` is
specified, pie plot of selected column will be drawn. If ``subplots=True`` is
specified, pie plots for each column are drawn as subplots. A legend will be
drawn in each pie plots by default; specify ``legend=False`` to hide it.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   plt.figure()

.. ipython:: python

   df = pd.DataFrame(
       3 * np.random.rand(4, 2), index=["a", "b", "c", "d"], columns=["x", "y"]
   )

   @savefig df_pie_plot.png
   df.plot.pie(subplots=True, figsize=(8, 4));

.. ipython:: python
   :suppress:

   plt.close("all")

You can use the ``labels`` and ``colors`` keywords to specify the labels and colors of each wedge.

.. warning::

   Most pandas plots use the ``label`` and ``color`` arguments (note the lack of "s" on those).
   To be consistent with :func:`matplotlib.pyplot.pie` you must use ``labels`` and ``colors``.

If you want to hide wedge labels, specify ``labels=None``.
If ``fontsize`` is specified, the value will be applied to wedge labels.
Also, other keywords supported by :func:`matplotlib.pyplot.pie` can be used.


.. ipython:: python
   :suppress:

   plt.figure()

.. ipython:: python

   @savefig series_pie_plot_options.png
   series.plot.pie(
       labels=["AA", "BB", "CC", "DD"],
       colors=["r", "g", "b", "c"],
       autopct="%.2f",
       fontsize=20,
       figsize=(6, 6),
   );

If you pass values whose sum total is less than 1.0 they will be rescaled so that they sum to 1.

.. ipython:: python
   :suppress:

   plt.close("all")
   plt.figure()

.. ipython:: python
   :okwarning:

   series = pd.Series([0.1] * 4, index=["a", "b", "c", "d"], name="series2")

   @savefig series_pie_plot_semi.png
   series.plot.pie(figsize=(6, 6));

See the `matplotlib pie documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pie.html>`__ for more.

.. ipython:: python
    :suppress:

    plt.close("all")

.. _visualization.missing_data:

Plotting with missing data
--------------------------

pandas tries to be pragmatic about plotting ``DataFrames`` or ``Series``
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

Plotting tools
--------------

These functions can be imported from ``pandas.plotting``
and take a :class:`Series` or :class:`DataFrame` as an argument.

.. _visualization.scatter_matrix:

Scatter matrix plot
~~~~~~~~~~~~~~~~~~~

You can create a scatter plot matrix using the
``scatter_matrix`` method in ``pandas.plotting``:

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.plotting import scatter_matrix

   df = pd.DataFrame(np.random.randn(1000, 4), columns=["a", "b", "c", "d"])

   @savefig scatter_matrix_kde.png
   scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal="kde");

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.kde:

Density plot
~~~~~~~~~~~~

You can create density plots using the :meth:`Series.plot.kde` and :meth:`DataFrame.plot.kde` methods.

.. ipython:: python
   :suppress:

   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   ser = pd.Series(np.random.randn(1000))

   @savefig kde_plot.png
   ser.plot.kde();

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.andrews_curves:

Andrews curves
~~~~~~~~~~~~~~

Andrews curves allow one to plot multivariate data as a large number
of curves that are created using the attributes of samples as coefficients
for Fourier series, see the `Wikipedia entry <https://en.wikipedia.org/wiki/Andrews_plot>`__
for more information. By coloring these curves differently for each class
it is possible to visualize data clustering. Curves belonging to samples
of the same class will usually be closer together and form larger structures.

**Note**: The "Iris" dataset is available `here <https://raw.githubusercontent.com/pandas-dev/pandas/main/pandas/tests/io/data/csv/iris.csv>`__.

.. ipython:: python

   from pandas.plotting import andrews_curves

   data = pd.read_csv("data/iris.data")

   plt.figure();

   @savefig andrews_curves.png
   andrews_curves(data, "Name");

.. _visualization.parallel_coordinates:

Parallel coordinates
~~~~~~~~~~~~~~~~~~~~

Parallel coordinates is a plotting technique for plotting multivariate data,
see the `Wikipedia entry <https://en.wikipedia.org/wiki/Parallel_coordinates>`__
for an introduction.
Parallel coordinates allows one to see clusters in data and to estimate other statistics visually.
Using parallel coordinates points are represented as connected line segments.
Each vertical line represents one attribute. One set of connected line segments
represents one data point. Points that tend to cluster will appear closer together.

.. ipython:: python

   from pandas.plotting import parallel_coordinates

   data = pd.read_csv("data/iris.data")

   plt.figure();

   @savefig parallel_coordinates.png
   parallel_coordinates(data, "Name");

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.lag:

Lag plot
~~~~~~~~

Lag plots are used to check if a data set or time series is random. Random
data should not exhibit any structure in the lag plot. Non-random structure
implies that the underlying data are not random. The ``lag`` argument may
be passed, and when ``lag=1`` the plot is essentially ``data[:-1]`` vs.
``data[1:]``.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.plotting import lag_plot

   plt.figure();

   spacing = np.linspace(-99 * np.pi, 99 * np.pi, num=1000)
   data = pd.Series(0.1 * np.random.rand(1000) + 0.9 * np.sin(spacing))

   @savefig lag_plot.png
   lag_plot(data);

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.autocorrelation:

Autocorrelation plot
~~~~~~~~~~~~~~~~~~~~

Autocorrelation plots are often used for checking randomness in time series.
This is done by computing autocorrelations for data values at varying time lags.
If time series is random, such autocorrelations should be near zero for any and
all time-lag separations. If time series is non-random then one or more of the
autocorrelations will be significantly non-zero. The horizontal lines displayed
in the plot correspond to 95% and 99% confidence bands. The dashed line is 99%
confidence band. See the
`Wikipedia entry <https://en.wikipedia.org/wiki/Correlogram>`__ for more about
autocorrelation plots.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   from pandas.plotting import autocorrelation_plot

   plt.figure();

   spacing = np.linspace(-9 * np.pi, 9 * np.pi, num=1000)
   data = pd.Series(0.7 * np.random.rand(1000) + 0.3 * np.sin(spacing))

   @savefig autocorrelation_plot.png
   autocorrelation_plot(data);

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.bootstrap:

Bootstrap plot
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

   from pandas.plotting import bootstrap_plot

   data = pd.Series(np.random.rand(1000))

   @savefig bootstrap_plot.png
   bootstrap_plot(data, size=50, samples=500, color="grey");

.. ipython:: python
   :suppress:

    plt.close("all")

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
See the R package `Radviz <https://cran.r-project.org/web/packages/Radviz/index.html>`__
for more information.

**Note**: The "Iris" dataset is available `here <https://raw.githubusercontent.com/pandas-dev/pandas/main/pandas/tests/io/data/csv/iris.csv>`__.

.. ipython:: python

   from pandas.plotting import radviz

   data = pd.read_csv("data/iris.data")

   plt.figure();

   @savefig radviz.png
   radviz(data, "Name");

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.formatting:

Plot formatting
---------------

Setting the plot style
~~~~~~~~~~~~~~~~~~~~~~

From version 1.5 and up, matplotlib offers a range of pre-configured plotting styles. Setting the
style can be used to easily give plots the general look that you want.
Setting the style is as easy as calling ``matplotlib.style.use(my_plot_style)`` before
creating your plot. For example you could write ``matplotlib.style.use('ggplot')`` for ggplot-style
plots.

You can see the various available style names at ``matplotlib.style.available`` and it's very
easy to try them out.

General plot style arguments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most plotting methods have a set of keyword arguments that control the
layout and formatting of the returned plot:

.. ipython:: python

   plt.figure();
   @savefig series_plot_basic2.png
   ts.plot(style="k--", label="Series");

.. ipython:: python
   :suppress:

   plt.close("all")

For each kind of plot (e.g. ``line``, ``bar``, ``scatter``) any additional arguments
keywords are passed along to the corresponding matplotlib function
(:meth:`ax.plot() <matplotlib.axes.Axes.plot>`,
:meth:`ax.bar() <matplotlib.axes.Axes.bar>`,
:meth:`ax.scatter() <matplotlib.axes.Axes.scatter>`). These can be used
to control additional styling, beyond what pandas provides.

Controlling the legend
~~~~~~~~~~~~~~~~~~~~~~

You may set the ``legend`` argument to ``False`` to hide the legend, which is
shown by default.

.. ipython:: python
   :suppress:

   np.random.seed(123456)

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list("ABCD"))
   df = df.cumsum()

   @savefig frame_plot_basic_noleg.png
   df.plot(legend=False);

.. ipython:: python
   :suppress:

   plt.close("all")


Controlling the labels
~~~~~~~~~~~~~~~~~~~~~~

You may set the ``xlabel`` and ``ylabel`` arguments to give the plot custom labels
for x and y axis. By default, pandas will pick up index name as xlabel, while leaving
it empty for ylabel.

.. ipython:: python

   df.plot();

   @savefig plot_xlabel_ylabel.png
   df.plot(xlabel="new x", ylabel="new y");

.. ipython:: python
   :suppress:

   plt.close("all")


Scales
~~~~~~

You may pass ``logy`` to get a log-scale Y axis.

.. ipython:: python
   :suppress:

   plt.figure()
   np.random.seed(123456)

.. ipython:: python

   ts = pd.Series(np.random.randn(1000), index=pd.date_range("1/1/2000", periods=1000))
   ts = np.exp(ts.cumsum())

   @savefig series_plot_logy.png
   ts.plot(logy=True);

.. ipython:: python
   :suppress:

   plt.close("all")

See also the ``logx`` and ``loglog`` keyword arguments.

Plotting on a secondary y-axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To plot data on a secondary y-axis, use the ``secondary_y`` keyword:

.. ipython:: python
   :suppress:

   plt.figure()

.. ipython:: python

   df["A"].plot();

   @savefig series_plot_secondary_y.png
   df["B"].plot(secondary_y=True, style="g");

.. ipython:: python
   :suppress:

   plt.close("all")

To plot some columns in a ``DataFrame``, give the column names to the ``secondary_y``
keyword:

.. ipython:: python

   plt.figure();
   ax = df.plot(secondary_y=["A", "B"])
   ax.set_ylabel("CD scale");
   @savefig frame_plot_secondary_y.png
   ax.right_ax.set_ylabel("AB scale");

.. ipython:: python
   :suppress:

   plt.close("all")

Note that the columns plotted on the secondary y-axis is automatically marked
with "(right)" in the legend. To turn off the automatic marking, use the
``mark_right=False`` keyword:

.. ipython:: python

   plt.figure();

   @savefig frame_plot_secondary_y_no_right.png
   df.plot(secondary_y=["A", "B"], mark_right=False);

.. ipython:: python
   :suppress:

   plt.close("all")

.. _plotting.formatters:

Custom formatters for timeseries plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pandas provides custom formatters for timeseries plots. These change the
formatting of the axis labels for dates and times. By default,
the custom formatters are applied only to plots created by pandas with
:meth:`DataFrame.plot` or :meth:`Series.plot`. To have them apply to all
plots, including those made by matplotlib, set the option
``pd.options.plotting.matplotlib.register_converters = True`` or use
:meth:`pandas.plotting.register_matplotlib_converters`.

Suppressing tick resolution adjustment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pandas includes automatic tick resolution adjustment for regular frequency
time-series data. For limited cases where pandas cannot infer the frequency
information (e.g., in an externally created ``twinx``), you can choose to
suppress this behavior for alignment purposes.

Here is the default behavior, notice how the x-axis tick labeling is performed:

.. ipython:: python

   plt.figure();

   @savefig ser_plot_suppress.png
   df["A"].plot();

.. ipython:: python
   :suppress:

   plt.close("all")

Using the ``x_compat`` parameter, you can suppress this behavior:

.. ipython:: python

   plt.figure();

   @savefig ser_plot_suppress_parm.png
   df["A"].plot(x_compat=True);

.. ipython:: python
   :suppress:

   plt.close("all")

If you have more than one plot that needs to be suppressed, the ``use`` method
in ``pandas.plotting.plot_params`` can be used in a ``with`` statement:

.. ipython:: python

   plt.figure();

   @savefig ser_plot_suppress_context.png
   with pd.plotting.plot_params.use("x_compat", True):
       df["A"].plot(color="r")
       df["B"].plot(color="g")
       df["C"].plot(color="b")

.. ipython:: python
   :suppress:

   plt.close("all")

Automatic date tick adjustment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``TimedeltaIndex`` now uses the native matplotlib
tick locator methods, it is useful to call the automatic
date tick adjustment from matplotlib for figures whose ticklabels overlap.

See the :meth:`autofmt_xdate <matplotlib.figure.autofmt_xdate>` method and the
`matplotlib documentation <https://matplotlib.org/2.0.2/users/recipes.html#fixing-common-date-annoyances>`__ for more.

Subplots
~~~~~~~~

Each ``Series`` in a ``DataFrame`` can be plotted on a different axis
with the ``subplots`` keyword:

.. ipython:: python

   @savefig frame_plot_subplots.png
   df.plot(subplots=True, figsize=(6, 6));

.. ipython:: python
   :suppress:

   plt.close("all")

Using layout and targeting multiple axes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The layout of subplots can be specified by the ``layout`` keyword. It can accept
``(rows, columns)``. The ``layout`` keyword can be used in
``hist`` and ``boxplot`` also. If the input is invalid, a ``ValueError`` will be raised.

The number of axes which can be contained by rows x columns specified by ``layout`` must be
larger than the number of required subplots. If layout can contain more axes than required,
blank axes are not drawn. Similar to a NumPy array's ``reshape`` method, you
can use ``-1`` for one dimension to automatically calculate the number of rows
or columns needed, given the other.

.. ipython:: python

   @savefig frame_plot_subplots_layout.png
   df.plot(subplots=True, layout=(2, 3), figsize=(6, 6), sharex=False);

.. ipython:: python
   :suppress:

   plt.close("all")

The above example is identical to using:

.. ipython:: python

   df.plot(subplots=True, layout=(2, -1), figsize=(6, 6), sharex=False);

.. ipython:: python
   :suppress:

   plt.close("all")

The required number of columns (3) is inferred from the number of series to plot
and the given number of rows (2).

You can pass multiple axes created beforehand as list-like via ``ax`` keyword.
This allows more complicated layouts.
The passed axes must be the same number as the subplots being drawn.

When multiple axes are passed via the ``ax`` keyword, ``layout``, ``sharex`` and ``sharey`` keywords
don't affect to the output. You should explicitly pass ``sharex=False`` and ``sharey=False``,
otherwise you will see a warning.

.. ipython:: python

   fig, axes = plt.subplots(4, 4, figsize=(9, 9))
   plt.subplots_adjust(wspace=0.5, hspace=0.5)
   target1 = [axes[0][0], axes[1][1], axes[2][2], axes[3][3]]
   target2 = [axes[3][0], axes[2][1], axes[1][2], axes[0][3]]

   df.plot(subplots=True, ax=target1, legend=False, sharex=False, sharey=False);
   @savefig frame_plot_subplots_multi_ax.png
   (-df).plot(subplots=True, ax=target2, legend=False, sharex=False, sharey=False);

.. ipython:: python
   :suppress:

   plt.close("all")

Another option is passing an ``ax`` argument to :meth:`Series.plot` to plot on a particular axis:

.. ipython:: python

   np.random.seed(123456)
   ts = pd.Series(np.random.randn(1000), index=pd.date_range("1/1/2000", periods=1000))
   ts = ts.cumsum()

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index, columns=list("ABCD"))
   df = df.cumsum()

.. ipython:: python
   :suppress:

   plt.close("all")

.. ipython:: python

   fig, axes = plt.subplots(nrows=2, ncols=2)
   plt.subplots_adjust(wspace=0.2, hspace=0.5)
   df["A"].plot(ax=axes[0, 0]);
   axes[0, 0].set_title("A");
   df["B"].plot(ax=axes[0, 1]);
   axes[0, 1].set_title("B");
   df["C"].plot(ax=axes[1, 0]);
   axes[1, 0].set_title("C");
   df["D"].plot(ax=axes[1, 1]);
   @savefig series_plot_multi.png
   axes[1, 1].set_title("D");

.. ipython:: python
   :suppress:

    plt.close("all")

.. _visualization.errorbars:

Plotting with error bars
~~~~~~~~~~~~~~~~~~~~~~~~

Plotting with error bars is supported in :meth:`DataFrame.plot` and :meth:`Series.plot`.

Horizontal and vertical error bars can be supplied to the ``xerr`` and ``yerr`` keyword arguments to :meth:`~DataFrame.plot`. The error values can be specified using a variety of formats:

* As a :class:`DataFrame` or ``dict`` of errors with column names matching the ``columns`` attribute of the plotting :class:`DataFrame` or matching the ``name`` attribute of the :class:`Series`.
* As a ``str`` indicating which of the columns of plotting :class:`DataFrame` contain the error values.
* As raw values (``list``, ``tuple``, or ``np.ndarray``). Must be the same length as the plotting :class:`DataFrame`/:class:`Series`.

Here is an example of one way to easily plot group means with standard deviations from the raw data.

.. ipython:: python

   # Generate the data
   ix3 = pd.MultiIndex.from_arrays(
       [
           ["a", "a", "a", "a", "a", "b", "b", "b", "b", "b"],
           ["foo", "foo", "foo", "bar", "bar", "foo", "foo", "bar", "bar", "bar"],
       ],
       names=["letter", "word"],
   )

   df3 = pd.DataFrame(
       {
           "data1": [9, 3, 2, 4, 3, 2, 4, 6, 3, 2],
           "data2": [9, 6, 5, 7, 5, 4, 5, 6, 5, 1],
       },
       index=ix3,
   )

   # Group by index labels and take the means and standard deviations
   # for each group
   gp3 = df3.groupby(level=("letter", "word"))
   means = gp3.mean()
   errors = gp3.std()
   means
   errors

   # Plot
   fig, ax = plt.subplots()
   @savefig errorbar_example.png
   means.plot.bar(yerr=errors, ax=ax, capsize=4, rot=0);

.. ipython:: python
   :suppress:

   plt.close("all")

Asymmetrical error bars are also supported, however raw error values must be provided in this case. For a ``N`` length :class:`Series`, a ``2xN`` array should be provided indicating lower and upper (or left and right) errors. For a ``MxN`` :class:`DataFrame`, asymmetrical errors should be in a ``Mx2xN`` array.

Here is an example of one way to plot the min/max range using asymmetrical error bars.

.. ipython:: python

   mins = gp3.min()
   maxs = gp3.max()

   # errors should be positive, and defined in the order of lower, upper
   errors = [[means[c] - mins[c], maxs[c] - means[c]] for c in df3.columns]

   # Plot
   fig, ax = plt.subplots()
   @savefig errorbar_asymmetrical_example.png
   means.plot.bar(yerr=errors, ax=ax, capsize=4, rot=0);

.. ipython:: python
   :suppress:

   plt.close("all")

.. _visualization.table:

Plotting tables
~~~~~~~~~~~~~~~

Plotting with matplotlib table is now supported in  :meth:`DataFrame.plot` and :meth:`Series.plot` with a ``table`` keyword. The ``table`` keyword can accept ``bool``, :class:`DataFrame` or :class:`Series`. The simple way to draw a table is to specify ``table=True``. Data will be transposed to meet matplotlib's default layout.

.. ipython:: python

   np.random.seed(123456)
   fig, ax = plt.subplots(1, 1, figsize=(7, 6.5))
   df = pd.DataFrame(np.random.rand(5, 3), columns=["a", "b", "c"])
   ax.xaxis.tick_top()  # Display x-axis ticks on top.

   @savefig line_plot_table_true.png
   df.plot(table=True, ax=ax);

.. ipython:: python
   :suppress:

   plt.close("all")

Also, you can pass a different :class:`DataFrame` or :class:`Series` to the
``table`` keyword. The data will be drawn as displayed in print method
(not transposed automatically). If required, it should be transposed manually
as seen in the example below.

.. ipython:: python

   fig, ax = plt.subplots(1, 1, figsize=(7, 6.75))
   ax.xaxis.tick_top()  # Display x-axis ticks on top.

   @savefig line_plot_table_data.png
   df.plot(table=np.round(df.T, 2), ax=ax);

.. ipython:: python
   :suppress:

   plt.close("all")

There also exists a helper function ``pandas.plotting.table``, which creates a
table from :class:`DataFrame` or :class:`Series`, and adds it to an
``matplotlib.Axes`` instance. This function can accept keywords which the
matplotlib `table <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.table.html>`__ has.

.. ipython:: python

   from pandas.plotting import table

   fig, ax = plt.subplots(1, 1)

   table(ax, np.round(df.describe(), 2), loc="upper right", colWidths=[0.2, 0.2, 0.2]);

   @savefig line_plot_table_describe.png
   df.plot(ax=ax, ylim=(0, 2), legend=None);

.. ipython:: python
   :suppress:

   plt.close("all")

**Note**: You can get table instances on the axes using ``axes.tables`` property for further decorations. See the `matplotlib table documentation <https://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.table>`__ for more.

.. _visualization.colormaps:

Colormaps
~~~~~~~~~

A potential issue when plotting a large number of columns is that it can be
difficult to distinguish some series due to repetition in the default colors. To
remedy this, ``DataFrame`` plotting supports the use of the ``colormap`` argument,
which accepts either a Matplotlib `colormap <https://matplotlib.org/api/cm_api.html>`__
or a string that is a name of a colormap registered with Matplotlib. A
visualization of the default matplotlib colormaps is available `here
<https://matplotlib.org/stable/gallery/color/colormap_reference.html>`__.

As matplotlib does not directly support colormaps for line-based plots, the
colors are selected based on an even spacing determined by the number of columns
in the ``DataFrame``. There is no consideration made for background color, so some
colormaps will produce lines that are not easily visible.

To use the cubehelix colormap, we can pass ``colormap='cubehelix'``.

.. ipython:: python

   np.random.seed(123456)
   df = pd.DataFrame(np.random.randn(1000, 10), index=ts.index)
   df = df.cumsum()

   plt.figure();

   @savefig cubehelix.png
   df.plot(colormap="cubehelix");

.. ipython:: python
   :suppress:

   plt.close("all")

Alternatively, we can pass the colormap itself:

.. ipython:: python

   from matplotlib import cm

   plt.figure();

   @savefig cubehelix_cm.png
   df.plot(colormap=cm.cubehelix);

.. ipython:: python
   :suppress:

   plt.close("all")

Colormaps can also be used other plot types, like bar charts:

.. ipython:: python

   np.random.seed(123456)
   dd = pd.DataFrame(np.random.randn(10, 10)).map(abs)
   dd = dd.cumsum()

   plt.figure();

   @savefig greens.png
   dd.plot.bar(colormap="Greens");

.. ipython:: python
   :suppress:

   plt.close("all")

Parallel coordinates charts:

.. ipython:: python

   plt.figure();

   @savefig parallel_gist_rainbow.png
   parallel_coordinates(data, "Name", colormap="gist_rainbow");

.. ipython:: python
   :suppress:

   plt.close("all")

Andrews curves charts:

.. ipython:: python

   plt.figure();

   @savefig andrews_curve_winter.png
   andrews_curves(data, "Name", colormap="winter");

.. ipython:: python
   :suppress:

   plt.close("all")

Plotting directly with Matplotlib
---------------------------------

In some situations it may still be preferable or necessary to prepare plots
directly with matplotlib, for instance when a certain type of plot or
customization is not (yet) supported by pandas. ``Series`` and ``DataFrame``
objects behave like arrays and can therefore be passed directly to
matplotlib functions without explicit casts.

pandas also automatically registers formatters and locators that recognize date
indices, thereby extending date and time support to practically all plot types
available in matplotlib. Although this formatting does not provide the same
level of refinement you would get when plotting via pandas, it can be faster
when plotting a large number of points.

.. ipython:: python

   np.random.seed(123456)
   price = pd.Series(
       np.random.randn(150).cumsum(),
       index=pd.date_range("2000-1-1", periods=150, freq="B"),
   )
   ma = price.rolling(20).mean()
   mstd = price.rolling(20).std()

   plt.figure();

   plt.plot(price.index, price, "k");
   plt.plot(ma.index, ma, "b");
   @savefig bollinger.png
   plt.fill_between(mstd.index, ma - 2 * mstd, ma + 2 * mstd, color="b", alpha=0.2);

.. ipython:: python
   :suppress:

    plt.close("all")

Plotting backends
-----------------

pandas can be extended with third-party plotting backends. The
main idea is letting users select a plotting backend different than the provided
one based on Matplotlib.

This can be done by passing 'backend.module' as the argument ``backend`` in ``plot``
function. For example:

.. code-block:: python

    >>> Series([1, 2, 3]).plot(backend="backend.module")

Alternatively, you can also set this option globally, do you don't need to specify
the keyword in each ``plot`` call. For example:

.. code-block:: python

    >>> pd.set_option("plotting.backend", "backend.module")
    >>> pd.Series([1, 2, 3]).plot()

Or:

.. code-block:: python

    >>> pd.options.plotting.backend = "backend.module"
    >>> pd.Series([1, 2, 3]).plot()

This would be more or less equivalent to:

.. code-block:: python

    >>> import backend.module
    >>> backend.module.plot(pd.Series([1, 2, 3]))

The backend module can then use other visualization tools (Bokeh, Altair, hvplot,...)
to generate the plots. Some libraries implementing a backend for pandas are listed
on `the ecosystem page <https://pandas.pydata.org/community/ecosystem.html>`_.

Developers guide can be found at
https://pandas.pydata.org/docs/dev/development/extending.html#plotting-backends
