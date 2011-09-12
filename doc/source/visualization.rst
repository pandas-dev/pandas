.. currentmodule:: pandas
.. _visualization:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')

************************
Plotting with matplotlib
************************

.. note::

    We intend to build more plotting integration with `matplotlib
    <http://matplotlib.sourceforge.net>`__ as time goes on.

We use the standard convention for referencing the matplotlib API:

.. ipython:: python

   import matplotlib.pyplot as plt

Basic plotting: ``plot``
------------------------

The ``plot`` method on Series and DataFrame is just a simple wrapper around
``plt.plot``:

.. ipython:: python

   ts = Series(randn(1000), index=DateRange('1/1/2000', periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png width=4.5in
   ts.plot()

If the index consists of dates, it calls ``gca().autofmt_xdate()`` to try to
format the x-axis nicely as per above. THe method takes a number of arguments
for controlling the look of the plot:

.. ipython:: python

   @savefig series_plot_basic2.png width=4.5in
   plt.figure(); ts.plot(style='k--', label='Series'); plt.legend()

On DataFrame, ``plot`` is a convenience to plot all of the columns with labels:

.. ipython:: python

   df = DataFrame(randn(1000, 4), index=ts.index,
                  columns=['A', 'B', 'C', 'D'])
   df = df.cumsum()

   @savefig frame_plot_basic.png width=4.5in
   plt.figure(); df.plot(); plt.legend(loc='best')

Some other options are available, like plotting each Series on a different axis:

.. ipython:: python

   @savefig frame_plot_subplots.png width=4.5in
   df.plot(subplots=True, figsize=(8, 8)); plt.legend(loc='best')

Targeting different subplots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can pass an ``ax`` argument to ``Series.plot`` to plot on a particular axis:

.. ipython:: python

   fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 5))
   df['A'].plot(ax=axes[0][0]); axes[0][0].set_title('A')
   df['B'].plot(ax=axes[0][1]); axes[0][1].set_title('B')
   df['C'].plot(ax=axes[1][0]); axes[1][0].set_title('C')

   @savefig series_plot_multi.png width=4.5in
   df['D'].plot(ax=axes[1][1]); axes[1][1].set_title('D')

Other plotting features
-----------------------

Plotting non-time series data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For labeled, non-time series data, you may wish to produce a bar plot:

.. ipython:: python

   plt.figure();

   @savefig bar_plot_ex.png width=4.5in
   df.ix[5].plot(kind='bar'); plt.axhline(0, color='k')

Histogramming
~~~~~~~~~~~~~
.. ipython:: python

   plt.figure();

   @savefig hist_plot_ex.png width=4.5in
   df['A'].diff().hist()

For a DataFrame, ``hist`` plots the histograms of the columns on multiple
subplots:

.. ipython:: python

   plt.figure()

   @savefig frame_hist_ex.png width=4.5in
   df.diff().hist(color='k', alpha=0.5, bins=50)
