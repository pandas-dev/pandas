.. _install:

.. currentmodule:: pandas

************
Installation
************

You have the option to install an `official release
<http://pypi.python.org/pypi/pandas>`__ or to build the `development version
<http://github.com/pydata/pandas>`__. If you choose to install from source and
are running Windows, you will have to ensure that you have a compatible C
compiler (MinGW or Visual Studio) installed. `How-to install MinGW on Windows
<http://docs.cython.org/src/tutorial/appendix.html>`__

Python version support
~~~~~~~~~~~~~~~~~~~~~~

Officially Python 2.5 to 2.7 and Python 3.1+, although Python 3 support is less
well tested. Python 2.4 support is being phased out since the userbase has
shrunk significantly. Continuing Python 2.4 support will require either monetary
development support or someone contributing to the project to maintain
compatibility.


Binary installers
~~~~~~~~~~~~~~~~~

Available on `PyPI <http://pypi.python.org/pypi/pandas>`__

Dependencies
~~~~~~~~~~~~

  * `NumPy <http://www.numpy.org>`__: 1.4.0 or higher. Recommend 1.5.1 or
    higher
  * `python-dateutil <http://labix.org/python-dateutil>`__ 1.5

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

  * `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
  * `PyTables <http://www.pytables.org>`__: necessary for HDF5-based storage
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `scikits.statsmodels <http://statsmodels.sourceforge.net/>`__
     * Needed for parts of :mod:`pandas.stats`
  * `pytz <http://pytz.sourceforge.net/>`__
     * Needed for time zone support with ``DateRange``

.. note::

   Without the optional dependencies, many useful features will not
   work. Hence, it is highly recommended that you install these. A packaged
   distribution like the `Enthought Python Distribution
   <http://enthought.com/products/epd.php>`__ may be worth considering.

Installing from source
~~~~~~~~~~~~~~~~~~~~~~
.. note::

   Installing from the git repository requires a recent installation of `Cython
   <http://cython.org>`__ as the cythonized C sources are no longer checked
   into source control. Released source distributions will contain the built C
   files. I recommend installing the latest Cython via ``easy_install -U
   Cython``

The source code is hosted at http://github.com/pydata/pandas, it can be checked
out using git and compiled / installed like so:

::

  git clone git://github.com/pydata/pandas.git
  cd pandas
  python setup.py install

On Windows, I suggest installing the MinGW compiler suite following the
directions linked to above. Once configured property, run the following on the
command line:

::

  python setup.py build --compiler=mingw32
  python setup.py install

Note that you will not be able to import pandas if you open an interpreter in
the source directory unless you build the C extensions in place:

::

  python setup.py build_ext --inplace


Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

pandas is equipped with an exhaustive set of unit tests covering about 97% of
the codebase as of this writing. To run it on your machine to verify that
everything is working (and you have all of the dependencies, soft and hard,
installed), make sure you have `nose
<http://readthedocs.org/docs/nose/en/latest/>`__ and run:

::

    $ nosetests pandas
    ..........................................................................
    .......................S..................................................
    ..........................................................................
    ..........................................................................
    ..........................................................................
    ..........................................................................
    ..........................................................................
    ..........................................................................
    ..........................................................................
    ..........................................................................
    .................S........................................................
    ....
    ----------------------------------------------------------------------
    Ran 818 tests in 21.631s

    OK (SKIP=2)
