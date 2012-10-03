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
   from pandas import read_csv
   from pandas.tools.plotting import radviz
   import pandas.tools.rplot as rplot
   plt.close('all')

**************************
Trellis plotting interface
**************************

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomHistogram())

   @savefig rplot1_tips.png width=6in
   plot.render(plt.gcf())

.. ipython:: python

   plt.figure()

   plot = rplot.RPlot(tips_data, x='totbill', y='tip')
   plot.add(rplot.TrellisGrid(['sex', 'smoker']))
   plot.add(rplot.GeomDensity())

   @savefig rplot2_tips.png width=6in
   plot.render(plt.gcf())