.. currentmodule:: pandas
.. _rplot:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   tips_data = read_csv('data/tips.csv')
   iris_data = read_csv('data/iris.data')
   from pandas import read_csv
   from pandas.tools.plotting import radviz
   import pandas.tools.rplot as rplot
   plt.close('all')

**************************
Trellis plotting interface
**************************

.. note::

   The tips data set can be downloaded `here
   <http://wesmckinney.com/files/tips.csv>`_. Once you download it execute

   .. code-block:: python

      from pandas import read_csv
      tips_data = read_csv('tips.csv')

   from the directory where you downloaded the file.

We import the rplot API:

.. ipython:: python

   import pandas.tools.rplot as rplot

--------
Examples
--------

RPlot is a flexible API for producing Trellis plots. These plots allow you to arrange data in a rectangular grid by values of certain attributes. 

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot1_tips.png
   plot.render(plt.gcf())

In the example above, data from the tips data set is arranged by the attributes 'sex' and 'smoker'. Since both of those attributes can take on one of two values, the resulting grid has two columns and two rows. A histogram is displayed for each cell of the grid.

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomDensity())

   @savefig rplot2_tips.png
   plot.render(plt.gcf())

Example above is the same as previous except the plot is set to kernel density estimation. This shows how easy it is to have different plots for the same Trellis structure.

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomScatter())
   plot.add(rplot.GeomPolyFit(degree=2))

   @savefig rplot3_tips.png
   plot.render(plt.gcf())

The plot above shows that it is possible to have two or more plots for the same data displayed on the same Trellis grid cell.

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomScatter())
   plot.add(rplot.GeomDensity2D())

   @savefig rplot4_tips.png
   plot.render(plt.gcf())

Above is a similar plot but with 2D kernel desnity estimation plot superimposed.

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', '.']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot5_tips.png
   plot.render(plt.gcf())

It is possible to only use one attribute for grouping data. The example above only uses 'sex' attribute. If the second grouping attribute is not specified, the plots will be arranged in a column.

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['.', 'smoker']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot6_tips.png
   plot.render(plt.gcf())

If the first grouping attribute is not specified the plots will be arranged in a row.

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='total_bill', y='tip')
   plot.add(rplot.TrellisGrid(['.', 'smoker']))
   plot.add(rplot.GeomHistogram())

   plot = rplot.RPlot(tips_data, x='tip', y='total_bill')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomPoint(size=80.0, colour=rplot.ScaleRandomColour('day'), shape=rplot.ScaleShape('size'), alpha=1.0))

   @savefig rplot7_tips.png
   plot.render(plt.gcf())

As shown above, scatter plots are also possible. Scatter plots allow you to map various data attributes to graphical properties of the plot. In the example above the colour and shape of the scatter plot graphical objects is mapped to 'day' and 'size' attributes respectively. You use scale objects to specify these mappings. The list of scale classes is given below with initialization arguments for quick reference.

------
Scales
------

::

   ScaleGradient(column, colour1, colour2)

This one allows you to map an attribute (specified by parameter column) value to the colour of a graphical object. The larger the value of the attribute the closer the colour will be to colour2, the smaller the value, the closer it will be to colour1.

::

   ScaleGradient2(column, colour1, colour2, colour3)

The same as ScaleGradient but interpolates linearly between three colours instead of two.

::

   ScaleSize(column, min_size, max_size, transform)

Map attribute value to size of the graphical object. Parameter min_size (default 5.0) is the minimum size of the graphical object, max_size (default 100.0) is the maximum size and transform is a one argument function that will be used to transform the attribute value (defaults to lambda x: x).

::

   ScaleShape(column)

Map the shape of the object to attribute value. The attribute has to be categorical.

::

   ScaleRandomColour(column)

Assign a random colour to a value of categorical attribute specified by column.
