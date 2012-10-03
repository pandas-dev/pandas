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

RPlot is a flexible API for producing Trellis plots. These plots allow you to arrange data in a rectangular grid by values of certain attributes. 

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot1_tips.png width=8in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomDensity())

   @savefig rplot2_tips.png width=8in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomScatter())
   plot.add(rplot.GeomPolyFit(degree=2))

   @savefig rplot3_tips.png width=8in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomScatter())
   plot.add(rplot.GeomDensity2D())

   @savefig rplot4_tips.png width=8in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', '.']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot5_tips.png width=8in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['.', 'smoker']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot6_tips.png width=8in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['.', 'smoker']))
   plot.add(rplot.GeomHistogram())

   plot = rplot.RPlot(tips_data, x='tip', y='totbill')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomPoint(size=80.0, colour=rplot.ScaleRandomColour('day'), shape=rplot.ScaleShape('size'), alpha=1.0))

   @savefig rplot7_tips.png width=8in
   plot.render(plt.gcf())
