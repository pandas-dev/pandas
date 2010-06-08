.. _overview:

.. currentmodule:: pandas

****************
Package overview
****************

:mod:`pandas` is a library providing, among other things, a set of
convenient and powerful data structures for working with labeled
statistical (financial, economic, econometric) data sets. We will
refer to this data as *time series* and *cross-sectional* (or
*longitudinal*) which are common terms in statistics and
econometrics. pandas has multiple target audiences:

 * Users of R or MATLAB who wish to switch to Python for interactive
   data analysis and implementation of statistical models

 * NumPy users who are looking for richer data structures for working
   with time series and cross-sectional data.

 * System developers who wish to have a robust and well-tested library
   for building production applications involving such data sets.

History
-------

pandas development began at AQR Capital Management (a quantitative
hedge fund) in April 2008. It was open-sourced at the end of 2009 and
continues to be actively used and maintained.

Data structures at a glance
---------------------------

.. csv-table::
    :header: "Dimensions", "Name", "Description"
    :widths: 10, 15, 50

    1, Series, "Most generic 1D structure"
    1, TimeSeries, "Series indexed by datetimes"
    2, DataFrame, "General 2D indexed tabular structure"
    2, DataMatrix, "Same API as DataFrame but faster (for most operations)"
    3, WidePanel, "General 3D panel data"
    3, LongPanel, "Stacked (2D) format panel data"

Why more than 1 data structure?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best way to think about the pandas data tructures is as flexible
containers for lower dimensional data. For example, DataFrame /
DataMatrix are containers for Series, and WidePanel is a container for
DataFrame / DataMatrix objects. We would like to be able to insert and
remove objects from these containers in a dictionary-like fashion.

Also, we would like sensible default behaviors for the common API
functions which take into account the typical orientation of time
series and cross-sectional data sets. When using ndarrays to store 2-
and 3-dimensional data, a burden is placed on the user to consider the
orientation of the data set when writing functions; axes are
considered more or less equivalent (except when C- or
Fortran-contiguousness matters for performance). In pandas, the axes
are intended to lend more semantic meaning to the data; i.e., for a
particular data set there is likely to be a "right" way to orient the
data. The goal, then, is to reduce the amount of thought required to
code up data transformations in downstream functions.

Lest we be too hand-wavy, here are some common use cases to
illustrate:

 (A) :ref:`DataFrame <dataframe>` containing multiple related time series

  * **columns**: "data type" associated with each time series
  * **index**:  dates shared by time series

 (B) :ref:`DataFrame <dataframe>` containing multiple cross-sections

  * **columns**: "data type" associated with each cross-section
  * **index**:  individual / entity labels common to cross-sections

 (C) :ref:`WidePanel <panel>` containing panel data

  * **items**: "data type" associated with each collection of time series
  * **major_axis**: dates shared by time series
  * **minor_axis**: individual / entity labels common to time series

Lastly, particularly if you don't buy the above explanation, having a
specialized vocabulary to refer to types of data sets often serves as
a benefit when discussing a dataset with other users (or reading their
code).

A quick note on mutation
~~~~~~~~~~~~~~~~~~~~~~~~

Most instance methods on the pandas data structures return a new
object, rather than updating the original object in-place. However,
when working with the contents (e.g. a column in a DataFrame),
mutations **will** be reflected in the original structure. In general,
though, we like to "favor immutability" where sensible.

What else is in the package?
----------------------------



Installation
------------

You have the option to install an official release or to build from
source. If you choose to install from source and are on Windows, you
will have to ensure that you have a compatible C compiler (gcc)
installed (see below).

Binary installers
~~~~~~~~~~~~~~~~~

Available from the Google Code website and PyPI.

Dependencies
~~~~~~~~~~~~

  * `NumPy <http://www.numpy.org>`__: 1.3.0 or higher
  * `dateutil <http://labix.org/python-dateutil>`__

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

  * `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
  * `PyTables <http://www.pytables.org>`__: necessary for HDF5-based storage
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `scikits.statsmodels <http://statsmodels.sourceforge.net/>`__
     * Needed for parts of :mod:`pandas.stats`

.. note::

   Without the optional dependencies, many useful features will not
   work. Hence, it is highly recommended that you install these.

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

The source code is hosted at http://pandas.googlecode.com, it can be
checked out using SVN and compiled / installed like so:

::

  svn co http://pandas.googlecode.com/svn/trunk/ pandas

  cd pandas

  python setup.py install

On Windows, you will need to download and install `gcc / MinGW
<http://www.mingw.org/wiki/HOWTO_Install_the_MinGW_GCC_Compiler_Suite>`__.
After adding it to your system path, you can install pandas by typing
instead:

::

  python setup.py build --compiler=mingw32
  python setup.py install
