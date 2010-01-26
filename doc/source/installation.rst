************
Installation
************

You have the option to install an official release or to build from
source. If you choose to install from source and are on Windows, you
will have to ensure that you have a compatible C compiler (gcc)
installed (see below).

Binary installers
-----------------

Available from the Google Code website and PyPI.

Dependencies
------------
  * `NumPy <http://www.numpy.org>`__: 1.3.0 or higher
  * `dateutil <http://labix.org/python-dateutil>`__

Optional dependencies
---------------------

  * `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `scikits.statsmodels <http://statsmodels.sourceforge.net/>`__
     * Needed for many parts of :mod:`pandas.stats`

.. note::

   Without the optional dependencies, many useful features will not
   work. Hence, it is highly recommended that you install these.

Installing from source
----------------------

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
