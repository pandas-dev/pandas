.. _overview:

.. currentmodule:: pandas

****************
Package overview
****************

:mod:`pandas` is a library providing a set of convenient and powerful
data structures for working with labeled statistical (financial,
economic, econometric) data sets. We will refer to this data as *time
series* and *cross-sectional* (or *longitudinal*) which are common
terms in statistics and econometrics. It has multiple target audiences:

  * Non-developers who wish to be able to easily manipulate data sets
in an interactive research environment.

History
-------

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



A quick note on mutation
~~~~~~~~~~~~~~~~~~~~~~~~

You will find that very few methods are capable of mutating a pandas
data structure like DataFrame. In general, the result of method calls
will return a new object (protecting the underlying data in the
calling object). So we like to "favor immutability" where sensible.

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
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `scikits.statsmodels <http://statsmodels.sourceforge.net/>`__
     * Needed for many parts of :mod:`pandas.stats`

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
After adding it to your system path , you can install pandas by typing
instead:

::

  python setup.py install --compiler=mingw32
