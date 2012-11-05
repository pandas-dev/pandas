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

.. _all-platforms:

All platforms
_____________

Stable installers available on `PyPI <http://pypi.python.org/pypi/pandas>`__

Preliminary builds and installers on the `Pandas download page <http://pandas.pydata.org/getpandas.html>`__ .

Overview
___________



.. csv-table::
    :header: "Platform", "Distribution", "Status", "Download / Repository Link", "Install method"
    :widths: 10, 10, 10, 20, 50


    Windows, all, stable, :ref:`all-platforms`,  ``pip install pandas``
    Mac, all, stable, :ref:`all-platforms`,  ``pip install pandas``
    Linux, Debian, stable, `official Debian repository <http://packages.debian.org/search?keywords=pandas&searchon=names&suite=all&section=all>`_ , ``sudo apt-get install python-pandas``
    Linux, Debian & Ubuntu, unstable (latest packages), `NeuroDebian <http://neuro.debian.net/index.html#how-to-use-this-repository>`_ , ``sudo apt-get install python-pandas``
    Linux, Ubuntu, stable, `official Ubuntu repository <http://packages.ubuntu.com/search?keywords=pandas&searchon=names&suite=all&section=all>`_ , ``sudo apt-get install python-pandas``
    Linux, Ubuntu, unstable (daily builds), `PythonXY PPA  <https://code.launchpad.net/~pythonxy/+archive/pythonxy-devel>`_; activate by: ``sudo add-apt-repository ppa:pythonxy/pythonxy-devel && sudo apt-get update``, ``sudo apt-get install python-pandas``
	Linux, OpenSuse & Fedora, stable, `OpenSuse Repository  <http://software.opensuse.org/package/python-pandas?search_term=pandas>`_ , ``zypper in  python-pandas``










Dependencies
~~~~~~~~~~~~

  * `NumPy <http://www.numpy.org>`__: 1.6.1 or higher
  * `python-dateutil <http://labix.org/python-dateutil>`__ 1.5

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

  * `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
  * `PyTables <http://www.pytables.org>`__: necessary for HDF5-based storage
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `statsmodels <http://statsmodels.sourceforge.net/>`__: 0.4.0 or higher
     * Needed for parts of :mod:`pandas.stats`
  * `pytz <http://pytz.sourceforge.net/>`__
     * Needed for time zone support with ``date_range``

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

The most recent version of MinGW (any installer dated after 2011-08-03)
has removed the '-mno-cygwin' option but Distutils has not yet been updated to
reflect that. Thus, you may run into an error like "unrecognized command line
option '-mno-cygwin'". Until the bug is fixed in Distutils, you may need to
install a slightly older version of MinGW (2011-08-02 installer).

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
