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
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   from dateutil.relativedelta import relativedelta
   from pandas.tseries.api import *
   from pandas.tseries.offsets import *
   import matplotlib.pyplot as plt
   plt.close('all')
   options.display.mpl_style='default'


.. _ref-repr-control:

How do I control the way my DataFrame is displayed?
---------------------------------------------------

Pandas users rely on a variety of environments for using pandas: scripts, terminal,
IPython qtconsole/ notebook, (IDLE, spyder, etc').
Each environment has it's own capabilities and limitations: HTML support,
horizontal scrolling, auto-detection of width/height.
To appropriately address all these environments, the display behavior is controlled
by several options, which you're encouraged to tweak to suit your setup.

As of 0.12, these are the relevant options, all under the `display` namespace,
(e.g. display.width,  etc'):

- notebook_repr_html: if True, IPython frontends with HTML support will display
  dataframes as HTML tables when possible.
- expand_repr (default True):  when the frame width cannot fit within the screen,
  the output will be broken into multiple pages to accomedate. This applies to
  textual (as opposed to HTML) display only.
- max_columns: max dataframe columns to display. a wider frame will trigger
  a summary view, unless `expand_repr` is True and HTML output is disabled.
- max_rows: max dataframe rows display. a longer frame will trigger a summary view.
- width: width of display screen in characters, used to determine the width of lines
  when expand_repr is active,  Setting this to None will trigger auto-detection of terminal
  width, this only works for proper terminals, not IPython frontends such as ipnb.
  width is ignored in IPython notebook, since the browser provides horizontal scrolling.

IPython users can use the IPython startup file to import pandas and set these
options automatically when starting up.


.. _ref-monkey-patching:

Adding Features to your Pandas Installation
-------------------------------------------

Pandas is a powerful tool and already has a plethora of data manipulation
operations implemented, most of them are very fast as well.
It's very possible however that certain functionality that would make your
life easier is missing. In that case you have several options:

1) Open an issue on `Github <https://github.com/pydata/pandas/issues/>`_ , explain your need and the sort of functionality you would like to see implemented.
2) Fork the repo, Implement the functionality yourself and open a PR
   on Github.
3) Write a method that performs the operation you are interested in and
   Monkey-patch the pandas class as part of your IPython profile startup
   or PYTHONSTARTUP file.

   For example, here is an example of adding an ``just_foo_cols()``
   method to the dataframe class:

.. ipython:: python

   import pandas as pd
   def just_foo_cols(self):
       """Get a list of column names containing the string 'foo'

       """
       return [x for x in self.columns if 'foo' in x]

   pd.DataFrame.just_foo_cols = just_foo_cols # monkey-patch the DataFrame class
   df = pd.DataFrame([range(4)],columns= ["A","foo","foozball","bar"])
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
(3D). Here is some code that resamples daily data to montly with
scikits.timeseries:

.. ipython:: python

   import scikits.timeseries as ts
   data = ts.time_series(np.random.randn(50), start_date='Jan-2000', freq='M')
   data
   data.convert('A', func=np.mean)

Here is the equivalent pandas code:

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

   x = np.array(range(10), '>i4') # big endian
   newx = x.byteswap().newbyteorder() # force native byteorder
   s = Series(newx)

See `the NumPy documentation on byte order
<http://docs.scipy.org/doc/numpy/user/basics.byteswapping.html>`__ for more
details.
