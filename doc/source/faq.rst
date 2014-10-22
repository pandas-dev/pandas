.. currentmodule:: pandas
.. _faq:

********************************
Frequently Asked Questions (FAQ)
********************************

.. ipython:: python
   :suppress:

   from datetime import datetime
   import numpy as np
   np.random.seed(123456)
   from pandas import *
   options.display.max_rows=15
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   from dateutil.relativedelta import relativedelta
   from pandas.tseries.api import *
   from pandas.tseries.offsets import *
   import matplotlib.pyplot as plt
   plt.close('all')
   options.display.mpl_style='default'
   from pandas.compat import lrange


.. _df-memory-usage:

DataFrame memory usage
~~~~~~~~~~~~~~~~~~~~~~
As of pandas version 0.15.0, the memory usage of a dataframe (including
the index) is shown when accessing the ``info`` method of a dataframe. A
configuration option, ``display.memory_usage`` (see :ref:`options`),
specifies if the dataframe's memory usage will be displayed when
invoking the ``df.info()`` method.

For example, the memory usage of the dataframe below is shown
when calling ``df.info()``:

.. ipython:: python

    dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
               'complex128', 'object', 'bool']
    n = 5000
    data = dict([ (t, np.random.randint(100, size=n).astype(t))
                    for t in dtypes])
    df = DataFrame(data)
    df['categorical'] = df['object'].astype('category')

    df.info()

The ``+`` symbol indicates that the true memory usage could be higher, because
pandas does not count the memory used by values in columns with
``dtype=object``.

By default the display option is set to ``True`` but can be explicitly
overridden by passing the ``memory_usage`` argument when invoking ``df.info()``.

The memory usage of each column can be found by calling the ``memory_usage``
method. This returns a Series with an index represented by column names
and memory usage of each column shown in bytes. For the dataframe above,
the memory usage of each column and the total memory usage of the
dataframe can be found with the memory_usage method:

.. ipython:: python

    df.memory_usage()

    # total memory usage of dataframe
    df.memory_usage().sum()

By default the memory usage of the dataframe's index is not shown in the
returned Series, the memory usage of the index can be shown by passing
the ``index=True`` argument:

.. ipython:: python

    df.memory_usage(index=True)

The memory usage displayed by the ``info`` method utilizes the
``memory_usage`` method to determine the memory usage of a dataframe
while also formatting the output in human-readable units (base-2
representation; i.e., 1KB = 1024 bytes).

See also :ref:`Categorical Memory Usage <categorical.memory>`.

.. _ref-monkey-patching:

Adding Features to your pandas Installation
-------------------------------------------

pandas is a powerful tool and already has a plethora of data manipulation
operations implemented, most of them are very fast as well.
It's very possible however that certain functionality that would make your
life easier is missing. In that case you have several options:

1) Open an issue on `Github <https://github.com/pydata/pandas/issues/>`__ , explain your need and the sort of functionality you would like to see implemented.
2) Fork the repo, Implement the functionality yourself and open a PR
   on Github.
3) Write a method that performs the operation you are interested in and
   Monkey-patch the pandas class as part of your IPython profile startup
   or PYTHONSTARTUP file.

   For example, here is an example of adding an ``just_foo_cols()``
   method to the dataframe class:

::

   import pandas as pd
   def just_foo_cols(self):
       """Get a list of column names containing the string 'foo'

       """
       return [x for x in self.columns if 'foo' in x]

   pd.DataFrame.just_foo_cols = just_foo_cols # monkey-patch the DataFrame class
   df = pd.DataFrame([list(range(4))], columns=["A","foo","foozball","bar"])
   df.just_foo_cols()
   del pd.DataFrame.just_foo_cols # you can also remove the new method


Monkey-patching is usually frowned upon because it makes your code
less portable and can cause subtle bugs in some circumstances.
Monkey-patching existing methods is usually a bad idea in that respect.
When used with proper care, however, it's a very useful tool to have.


.. _ref-scikits-migration:

Migrating from scikits.timeseries to pandas >= 0.8.0
----------------------------------------------------

Starting with pandas 0.8.0, users of scikits.timeseries should have all of the
features that they need to migrate their code to use pandas. Portions of the
scikits.timeseries codebase for implementing calendar logic and timespan
frequency conversions (but **not** resampling, that has all been implemented
from scratch from the ground up) have been ported to the pandas codebase.

The scikits.timeseries notions of ``Date`` and ``DateArray`` are responsible
for implementing calendar logic:

::

    In [16]: dt = ts.Date('Q', '1984Q3')

    # sic
    In [17]: dt
    Out[17]: <Q-DEC : 1984Q1>

    In [18]: dt.asfreq('D', 'start')
    Out[18]: <D : 01-Jan-1984>

    In [19]: dt.asfreq('D', 'end')
    Out[19]: <D : 31-Mar-1984>

    In [20]: dt + 3
    Out[20]: <Q-DEC : 1984Q4>

``Date`` and ``DateArray`` from scikits.timeseries have been reincarnated in
pandas ``Period`` and ``PeriodIndex``:

.. ipython:: python

   pnow('D')  # scikits.timeseries.now()
   Period(year=2007, month=3, day=15, freq='D')
   p = Period('1984Q3')
   p
   p.asfreq('D', 'start')
   p.asfreq('D', 'end')
   (p + 3).asfreq('T') + 6 * 60 + 30
   rng = period_range('1990', '2010', freq='A')
   rng
   rng.asfreq('B', 'end') - 3

.. csv-table::
    :header: "scikits.timeseries", "pandas", "Notes"
    :widths: 20, 20, 60

    Date, Period, "A span of time, from yearly through to secondly"
    DateArray, PeriodIndex, "An array of timespans"
    convert, resample, "Frequency conversion in scikits.timeseries"
    convert_to_annual, pivot_annual, "currently supports up to daily frequency, see :issue:`736`"


PeriodIndex / DateArray properties and functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scikits.timeseries ``DateArray`` had a number of information
properties. Here are the pandas equivalents:

.. csv-table::
    :header: "scikits.timeseries", "pandas", "Notes"
    :widths: 20, 60, 20

    get_steps, ``np.diff(idx.values)``,
    has_missing_dates, ``not idx.is_full``,
    is_full, ``idx.is_full``,
    is_valid, ``idx.is_monotonic and idx.is_unique``,
    is_chronological, ``is_monotonic``,
    ``arr.sort_chronologically()``, ``idx.order()``,

Frequency conversion
~~~~~~~~~~~~~~~~~~~~

Frequency conversion is implemented using the ``resample`` method on TimeSeries
and DataFrame objects (multiple time series). ``resample`` also works on panels
(3D). Here is some code that resamples daily data to monthly:

.. ipython:: python

   rng = period_range('Jan-2000', periods=50, freq='M')
   data = Series(np.random.randn(50), index=rng)
   data
   data.resample('A', how=np.mean)

Plotting
~~~~~~~~

Much of the plotting functionality of scikits.timeseries has been ported and
adopted to pandas's data structures. For example:

.. ipython:: python

   rng = period_range('1987Q2', periods=10, freq='Q-DEC')
   data = Series(np.random.randn(10), index=rng)

   @savefig skts_ts_plot.png
   plt.figure(); data.plot()

Converting to and from period format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the ``to_timestamp`` and ``to_period`` instance methods.

Treatment of missing data
~~~~~~~~~~~~~~~~~~~~~~~~~

Unlike scikits.timeseries, pandas data structures are not based on NumPy's
``MaskedArray`` object. Missing data is represented as ``NaN`` in numerical
arrays and either as ``None`` or ``NaN`` in non-numerical arrays. Implementing
a version of pandas's data structures that use MaskedArray is possible but
would require the involvement of a dedicated maintainer. Active pandas
developers are not interested in this.

Resampling with timestamps and periods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``resample`` has a ``kind`` argument which allows you to resample time series
with a DatetimeIndex to PeriodIndex:

.. ipython:: python

   rng = date_range('1/1/2000', periods=200, freq='D')
   data = Series(np.random.randn(200), index=rng)
   data[:10]
   data.index
   data.resample('M', kind='period')

Similarly, resampling from periods to timestamps is possible with an optional
interval (``'start'`` or ``'end'``) convention:

.. ipython:: python

   rng = period_range('Jan-2000', periods=50, freq='M')
   data = Series(np.random.randn(50), index=rng)
   resampled = data.resample('A', kind='timestamp', convention='end')
   resampled.index


Byte-Ordering Issues
--------------------
Occasionally you may have to deal with data that were created on a machine with
a different byte order than the one on which you are running Python. To deal
with this issue you should convert the underlying NumPy array to the native
system byte order *before* passing it to Series/DataFrame/Panel constructors
using something similar to the following:

.. ipython:: python

   x = np.array(list(range(10)), '>i4') # big endian
   newx = x.byteswap().newbyteorder() # force native byteorder
   s = Series(newx)

See `the NumPy documentation on byte order
<http://docs.scipy.org/doc/numpy/user/basics.byteswapping.html>`__ for more
details.


Visualizing Data in Qt applications
-----------------------------------

There is experimental support for visualizing DataFrames in PyQt4 and PySide
applications. At the moment you can display and edit the values of the cells
in the DataFrame. Qt will take care of displaying just the portion of the
DataFrame that is currently visible and the edits will be immediately saved to
the underlying DataFrame

To demonstrate this we will create a simple PySide application that will switch
between two editable DataFrames. For this will use the ``DataFrameModel`` class
that handles the access to the DataFrame, and the ``DataFrameWidget``, which is
just a thin layer around the ``QTableView``.

.. code-block:: python

	import numpy as np
	import pandas as pd
	from pandas.sandbox.qtpandas import DataFrameModel, DataFrameWidget
	from PySide import QtGui, QtCore

	# Or if you use PyQt4:
	# from PyQt4 import QtGui, QtCore

	class MainWidget(QtGui.QWidget):
	    def __init__(self, parent=None):
	        super(MainWidget, self).__init__(parent)

	        # Create two DataFrames
	        self.df1 = pd.DataFrame(np.arange(9).reshape(3, 3),
	                                columns=['foo', 'bar', 'baz'])
	        self.df2 = pd.DataFrame({
	                'int': [1, 2, 3],
	                'float': [1.5, 2.5, 3.5],
	                'string': ['a', 'b', 'c'],
	                'nan': [np.nan, np.nan, np.nan]
	            }, index=['AAA', 'BBB', 'CCC'],
	            columns=['int', 'float', 'string', 'nan'])

	        # Create the widget and set the first DataFrame
	        self.widget = DataFrameWidget(self.df1)

	        # Create the buttons for changing DataFrames
	        self.button_first = QtGui.QPushButton('First')
	        self.button_first.clicked.connect(self.on_first_click)
	        self.button_second = QtGui.QPushButton('Second')
	        self.button_second.clicked.connect(self.on_second_click)

	        # Set the layout
	        vbox = QtGui.QVBoxLayout()
	        vbox.addWidget(self.widget)
	        hbox = QtGui.QHBoxLayout()
	        hbox.addWidget(self.button_first)
	        hbox.addWidget(self.button_second)
	        vbox.addLayout(hbox)
	        self.setLayout(vbox)

	    def on_first_click(self):
	    	'''Sets the first DataFrame'''
	        self.widget.setDataFrame(self.df1)

	    def on_second_click(self):
	    	'''Sets the second DataFrame'''
	        self.widget.setDataFrame(self.df2)

	if __name__ == '__main__':
	    import sys

	    # Initialize the application
	    app = QtGui.QApplication(sys.argv)
	    mw = MainWidget()
	    mw.show()
	    app.exec_()
