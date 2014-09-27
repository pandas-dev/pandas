.. _install:

.. currentmodule:: pandas

============
Installation
============

The easiest way for the majority of users to install pandas is to install it
as part of the `Anaconda <http://docs.continuum.io/anaconda/>`__ distribution, a
cross platform distribution for data analysis and scientific computing.
This is the recommended installation method for most users.

Instructions for installing from source,
`PyPI <http://pypi.python.org/pypi/pandas>`__, various Linux distributions, or a
`development version <http://github.com/pydata/pandas>`__ are also provided.

Python version support
----------------------

Officially Python 2.6, 2.7, 3.2, 3.3, and 3.4.

Installing pandas
-----------------

Trying out pandas, no installation required!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to start experimenting with pandas doesn't involve installing
pandas at all.

`Wakari <https://wakari.io>`__ is a free service that provides a hosted
`IPython Notebook <http://ipython.org/notebook.html>`__ service in the cloud.

Simply create an account, and have access to pandas from within your brower via
an `IPython Notebook <http://ipython.org/notebook.html>`__ in a few minutes.

Installing pandas with Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing pandas and the rest of the `NumPy <http://www.numpy.org/>`__ and
`SciPy <http://www.scipy.org/>`__ stack can be a little
difficult for inexperienced users.

The simplest way to install not only pandas, but Python and the most popular
packages that make up the `SciPy <http://www.scipy.org/>`__ stack
(`IPython <http://ipython.org/>`__, `NumPy <http://www.numpy.org/>`__,
`Matplotlib <http://matplotlib.org/>`__, ...) is with
`Anaconda <http://docs.continuum.io/anaconda/>`__, a cross-platform
(Linux, Mac OS X, Windows) Python distribution for data analytics and
scientific computing.

After running a simple installer, the user will have access to pandas and the
rest of the `SciPy <http://www.scipy.org/>`__ stack without needing to install
anything else, and without needing to wait for any software to be compiled.

Installation instructions for `Anaconda <http://docs.continuum.io/anaconda/>`__
`can be found here <http://docs.continuum.io/anaconda/install.html>`__.

A full list of the packages available as part of the
`Anaconda <http://docs.continuum.io/anaconda/>`__ distribution
`can be found here <http://docs.continuum.io/anaconda/pkg-docs.html>`__.

An additional advantage of installing with Anaconda is that you don't require
admin rights to install it, it will install in the user's home directory, and
this also makes it trivial to delete Anaconda at a later date (just delete
that folder).

Installing pandas with Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous section outlined how to get pandas installed as part of the
`Anaconda <http://docs.continuum.io/anaconda/>`__ distribution.
However this approach means you will install well over one hundred packages
and involves downloading the installer which is a few hundred megabytes in size.

If you want to have more control on which packages, or have a limited internet
bandwidth, then installing pandas with
`Miniconda <http://conda.pydata.org/miniconda.html>`__ may be a better solution.

`Conda <http://conda.pydata.org/docs/>`__ is the package manager that the
`Anaconda <http://docs.continuum.io/anaconda/>`__ distribution is built upon.
It is a package manager that is both cross-platform and language agnostic
(it can play a similar role to a pip and virtualenv combination).

`Miniconda <http://conda.pydata.org/miniconda.html>`__ allows you to create a
minimal self contained Python installation, and then use the
`Conda <http://conda.pydata.org/docs/>`__ command to install additional packages.

First you will need `Conda <http://conda.pydata.org/docs/>`__ to be installed and
downloading and running the `Miniconda
<http://conda.pydata.org/miniconda.html>`__
will do this for you. The installer
`can be found here <http://conda.pydata.org/miniconda.html>`__

The next step is to create a new conda environment (these are analogous to a
virtualenv but they also allow you to specify precisely which Python version
to install also). Run the following commands from a terminal window::

    conda create -n name_of_my_env python

This will create a minimal environment with only Python installed in it.
To put your self inside this environment run::

    source activate name_of_my_env

On Windows the command is::

    activate name_of_my_env

The final step required is to install pandas. This can be done with the
following command::

    conda install pandas

To install a specific pandas version::

    conda install pandas=0.13.1

To install other packages, IPython for example::

    conda install ipython

To install the full `Anaconda <http://docs.continuum.io/anaconda/>`__
distribution::

    conda install anaconda

If you require any packages that are available to pip but not conda, simply
install pip, and use pip to install these packages::

    conda install pip
    pip install django

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

pandas can be installed via pip from
`PyPI <http://pypi.python.org/pypi/pandas>`__.

::

    pip install pandas

This will likely require the installation of a number of dependencies,
including NumPy, will require a compiler to compile required bits of code,
and can take a few minutes to complete.

Installing using your Linux distribution's package manager.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



.. csv-table::
    :header: "Distribution", "Status", "Download / Repository Link", "Install method"
    :widths: 10, 10, 20, 50


    Debian, stable, `official Debian repository <http://packages.debian.org/search?keywords=pandas&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install python-pandas``
    Debian & Ubuntu, unstable (latest packages), `NeuroDebian <http://neuro.debian.net/index.html#how-to-use-this-repository>`__ , ``sudo apt-get install python-pandas``
    Ubuntu, stable, `official Ubuntu repository <http://packages.ubuntu.com/search?keywords=pandas&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install python-pandas``
    Ubuntu, unstable (daily builds), `PythonXY PPA  <https://code.launchpad.net/~pythonxy/+archive/pythonxy-devel>`__; activate by: ``sudo add-apt-repository ppa:pythonxy/pythonxy-devel && sudo apt-get update``, ``sudo apt-get install python-pandas``
	OpenSuse & Fedora, stable, `OpenSuse Repository  <http://software.opensuse.org/package/python-pandas?search_term=pandas>`__ , ``zypper in  python-pandas``










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

Dependencies
------------

* `NumPy <http://www.numpy.org>`__: 1.7.0 or higher
* `python-dateutil <http://labix.org/python-dateutil>`__ 1.5
* `pytz <http://pytz.sourceforge.net/>`__
   * Needed for time zone support

.. _install.recommended_dependencies:

Recommended Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

* `numexpr <http://code.google.com/p/numexpr/>`__: for accelerating certain numerical operations.
  ``numexpr`` uses multiple cores as well as smart chunking and caching to achieve large speedups.
  If installed, must be Version 2.1 or higher.

* `bottleneck <http://berkeleyanalytics.com/bottleneck>`__: for accelerating certain types of ``nan``
  evaluations. ``bottleneck`` uses specialized cython routines to achieve large speedups.

.. note::

   You are highly encouraged to install these libraries, as they provide large speedups, especially
   if working with large data sets.


.. _install.optional_dependencies:

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~

* `Cython <http://www.cython.org>`__: Only necessary to build development
  version. Version 0.19.1 or higher.
* `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
* `PyTables <http://www.pytables.org>`__: necessary for HDF5-based storage. Version 3.0.0 or higher required.
* `SQLAlchemy <http://www.sqlalchemy.org>`__: for SQL database support. Version 0.8.1 or higher recommended.
* `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
* `statsmodels <http://statsmodels.sourceforge.net/>`__
   * Needed for parts of :mod:`pandas.stats`
* `openpyxl <http://packages.python.org/openpyxl/>`__, `xlrd/xlwt <http://www.python-excel.org/>`__
   * openpyxl version 1.6.1 or higher, but lower than 2.0.0
   * Needed for Excel I/O
* `XlsxWriter <https://pypi.python.org/pypi/XlsxWriter>`__
   * Alternative Excel writer.
* `boto <https://pypi.python.org/pypi/boto>`__: necessary for Amazon S3
  access.
* One of `PyQt4
  <http://www.riverbankcomputing.com/software/pyqt/download>`__, `PySide
  <http://qt-project.org/wiki/Category:LanguageBindings::PySide>`__, `pygtk
  <http://www.pygtk.org/>`__, `xsel
  <http://www.vergenet.net/~conrad/software/xsel/>`__, or `xclip
  <http://sourceforge.net/projects/xclip/>`__: necessary to use
  :func:`~pandas.io.clipboard.read_clipboard`. Most package managers on Linux
  distributions will have xclip and/or xsel immediately available for
  installation.
* Google's `python-gflags <http://code.google.com/p/python-gflags/>`__
  and `google-api-python-client <http://github.com/google/google-api-python-client>`__
   * Needed for :mod:`~pandas.io.gbq`
* `setuptools <https://pypi.python.org/pypi/setuptools/>`__
   * Needed for :mod:`~pandas.io.gbq` (specifically, it utilizes `pkg_resources`)
* `httplib2 <http://pypi.python.org/pypi/httplib2>`__
   * Needed for :mod:`~pandas.io.gbq`
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
   distribution like `Enthought Canopy
   <http://enthought.com/products/canopy>`__ may be worth considering.

