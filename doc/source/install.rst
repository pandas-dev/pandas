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

Officially Python 2.6 to 2.7 and Python 3.1+, although Python 3 support is less
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
  * `pytz <http://pytz.sourceforge.net/>`__
     * Needed for time zone support

.. _install.recommended_dependencies:

Recommended Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

  * `numexpr <http://code.google.com/p/numexpr/>`__: for accelerating certain numerical operations.
    ``numexpr`` uses multiple cores as well as smart chunking and caching to achieve large speedups.
  * `bottleneck <http://berkeleyanalytics.com/bottleneck>`__: for accelerating certain types of ``nan``
    evaluations. ``bottleneck`` uses specialized cython routines to achieve large speedups.

.. note::

   You are highly encouraged to install these libraries, as they provide large speedups, especially
   if working with large data sets.


.. _install.optional_dependencies:

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~

  * `Cython <http://www.cython.org>`__: Only necessary to build development
    version. Version 0.17.1 or higher.
  * `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
  * `PyTables <http://www.pytables.org>`__: necessary for HDF5-based storage
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `statsmodels <http://statsmodels.sourceforge.net/>`__
     * Needed for parts of :mod:`pandas.stats`
  * `openpyxl <http://packages.python.org/openpyxl/>`__, `xlrd/xlwt <http://www.python-excel.org/>`__
     * openpyxl version 1.6.1 or higher
     * Needed for Excel I/O
  * `boto <https://pypi.python.org/pypi/boto>`__: necessary for Amazon S3
    access.
  * One of the following combinations of libraries is needed to use the
    top-level :func:`~pandas.io.html.read_html` function:

    * `BeautifulSoup4`_ and `html5lib`_ (Any recent version of `html5lib`_ is
      okay.)
    * `BeautifulSoup4`_ and `lxml`_ 
    * `BeautifulSoup4`_ and `html5lib`_ and `lxml`_ 
    * Only `lxml`_, although see :ref:`HTML reading gotchas <html-gotchas>`
      for reasons as to why you should probably **not** take this approach.

    .. warning::

       * if you install `BeautifulSoup4`_ you must install either
         `lxml`_ or `html5lib`_ or both.
         :func:`~pandas.io.html.read_html` will **not** work with *only*
         `BeautifulSoup4`_ installed.
       * You are highly encouraged to read :ref:`HTML reading gotchas
         <html-gotchas>`. It explains issues surrounding the installation and
         usage of the above three libraries
       * You may need to install an older version of `BeautifulSoup4`_:
           - Versions 4.2.1, 4.1.3 and 4.0.2 have been confirmed for 64 and
             32-bit Ubuntu/Debian
       * Additionally, if you're using `Anaconda`_ you should definitely
         read :ref:`the gotchas about HTML parsing libraries <html-gotchas>`

    .. note::

       * if you're on a system with ``apt-get`` you can do

         .. code-block:: sh

            sudo apt-get build-dep python-lxml

         to get the necessary dependencies for installation of `lxml`_. This
         will prevent further headaches down the line.


.. _html5lib: https://github.com/html5lib/html5lib-python
.. _BeautifulSoup4: http://www.crummy.com/software/BeautifulSoup
.. _lxml: http://lxml.de
.. _Anaconda: https://store.continuum.io/cshop/anaconda

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

Make sure you have Cython installed when installing from the repository,
rather then a tarball or pypi.

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
