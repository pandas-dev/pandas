.. _install:

{{ header }}

============
Installation
============

The easiest way to install pandas is to install it
as part of the `Anaconda <https://docs.continuum.io/anaconda/>`__ distribution, a
cross platform distribution for data analysis and scientific computing.
This is the recommended installation method for most users.

Instructions for installing from source,
`PyPI <https://pypi.org/project/pandas>`__, `ActivePython <https://www.activestate.com/products/python/>`__, various Linux distributions, or a
`development version <https://github.com/pandas-dev/pandas>`__ are also provided.

.. _install.version:

Python version support
----------------------

Officially Python 3.8, 3.9 and 3.10.

Installing pandas
-----------------

.. _install.anaconda:

Installing with Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~

Installing pandas and the rest of the `NumPy <https://numpy.org/>`__ and
`SciPy <https://scipy.org/>`__ stack can be a little
difficult for inexperienced users.

The simplest way to install not only pandas, but Python and the most popular
packages that make up the `SciPy <https://scipy.org/>`__ stack
(`IPython <https://ipython.org/>`__, `NumPy <https://numpy.org/>`__,
`Matplotlib <https://matplotlib.org/>`__, ...) is with
`Anaconda <https://docs.continuum.io/anaconda/>`__, a cross-platform
(Linux, macOS, Windows) Python distribution for data analytics and
scientific computing.

After running the installer, the user will have access to pandas and the
rest of the `SciPy <https://scipy.org/>`__ stack without needing to install
anything else, and without needing to wait for any software to be compiled.

Installation instructions for `Anaconda <https://docs.continuum.io/anaconda/>`__
`can be found here <https://docs.continuum.io/anaconda/install/>`__.

A full list of the packages available as part of the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution
`can be found here <https://docs.continuum.io/anaconda/packages/pkg-docs/>`__.

Another advantage to installing Anaconda is that you don't need
admin rights to install it. Anaconda can install in the user's home directory,
which makes it trivial to delete Anaconda if you decide (just delete
that folder).

.. _install.miniconda:

Installing with Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~

The previous section outlined how to get pandas installed as part of the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution.
However this approach means you will install well over one hundred packages
and involves downloading the installer which is a few hundred megabytes in size.

If you want to have more control on which packages, or have a limited internet
bandwidth, then installing pandas with
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ may be a better solution.

`Conda <https://conda.io/en/latest/>`__ is the package manager that the
`Anaconda <https://docs.continuum.io/anaconda/>`__ distribution is built upon.
It is a package manager that is both cross-platform and language agnostic
(it can play a similar role to a pip and virtualenv combination).

`Miniconda <https://conda.pydata.org/miniconda.html>`__ allows you to create a
minimal self contained Python installation, and then use the
`Conda <https://conda.io/en/latest/>`__ command to install additional packages.

First you will need `Conda <https://conda.io/en/latest/>`__ to be installed and
downloading and running the `Miniconda
<https://conda.pydata.org/miniconda.html>`__
will do this for you. The installer
`can be found here <https://conda.pydata.org/miniconda.html>`__

The next step is to create a new conda environment. A conda environment is like a
virtualenv that allows you to specify a specific version of Python and set of libraries.
Run the following commands from a terminal window::

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

    conda install pandas=0.20.3

To install other packages, IPython for example::

    conda install ipython

To install the full `Anaconda <https://docs.continuum.io/anaconda/>`__
distribution::

    conda install anaconda

If you need packages that are available to pip but not conda, then
install pip, and then use pip to install those packages::

    conda install pip
    pip install django

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

pandas can be installed via pip from
`PyPI <https://pypi.org/project/pandas>`__.

.. note::
    You must have ``pip>=19.3`` to install from PyPI.

::

    pip install pandas

Installing with ActivePython
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation instructions for
`ActivePython <https://www.activestate.com/products/python/>`__ can be found
`here <https://www.activestate.com/products/python/>`__. Versions
2.7, 3.5 and 3.6 include pandas.

Installing using your Linux distribution's package manager.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The commands in this table will install pandas for Python 3 from your distribution.

.. csv-table::
    :header: "Distribution", "Status", "Download / Repository Link", "Install method"
    :widths: 10, 10, 20, 50


    Debian, stable, `official Debian repository <https://packages.debian.org/search?keywords=pandas&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install python3-pandas``
    Debian & Ubuntu, unstable (latest packages), `NeuroDebian <https://neuro.debian.net/index.html#how-to-use-this-repository>`__ , ``sudo apt-get install python3-pandas``
    Ubuntu, stable, `official Ubuntu repository <https://packages.ubuntu.com/search?keywords=pandas&searchon=names&suite=all&section=all>`__ , ``sudo apt-get install python3-pandas``
    OpenSuse, stable, `OpenSuse Repository  <https://software.opensuse.org/package/python-pandas?search_term=pandas>`__ , ``zypper in python3-pandas``
    Fedora, stable, `official Fedora repository  <https://src.fedoraproject.org/rpms/python-pandas>`__ , ``dnf install python3-pandas``
    Centos/RHEL, stable, `EPEL repository <https://admin.fedoraproject.org/pkgdb/package/rpms/python-pandas/>`__ , ``yum install python3-pandas``

**However**, the packages in the linux package managers are often a few versions behind, so
to get the newest version of pandas, it's recommended to install using the ``pip`` or ``conda``
methods described above.

Handling ImportErrors
~~~~~~~~~~~~~~~~~~~~~~

If you encounter an ImportError, it usually means that Python couldn't find pandas in the list of available
libraries. Python internally has a list of directories it searches through, to find packages. You can
obtain these directories with::

            import sys
            sys.path

One way you could be encountering this error is if you have multiple Python installations on your system
and you don't have pandas installed in the Python installation you're currently using.
In Linux/Mac you can run ``which python`` on your terminal and it will tell you which Python installation you're
using. If it's something like "/usr/bin/python", you're using the Python from the system, which is not recommended.

It is highly recommended to use ``conda``, for quick installation and for package and dependency updates.
You can find simple installation instructions for pandas in this document: ``installation instructions </getting_started.html>``.

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

See the :ref:`contributing guide <contributing>` for complete instructions on building from the git source tree. Further, see :ref:`creating a development environment <contributing_environment>` if you wish to create a pandas development environment.

Running the test suite
----------------------

pandas is equipped with an exhaustive set of unit tests, covering about 97% of
the code base as of this writing. To run it on your machine to verify that
everything is working (and that you have all of the dependencies, soft and hard,
installed), make sure you have `pytest
<https://docs.pytest.org/en/latest/>`__ >= 6.0 and `Hypothesis
<https://hypothesis.readthedocs.io/en/latest/>`__ >= 6.13.0, then run:

::

    >>> pd.test()
    running: pytest --skip-slow --skip-network --skip-db /home/user/anaconda3/lib/python3.9/site-packages/pandas

    ============================= test session starts ==============================
    platform linux -- Python 3.9.7, pytest-6.2.5, py-1.11.0, pluggy-1.0.0
    rootdir: /home/user
    plugins: dash-1.19.0, anyio-3.5.0, hypothesis-6.29.3
    collected 154975 items / 4 skipped / 154971 selected
    ........................................................................ [  0%]
    ........................................................................ [ 99%]
    .......................................                                  [100%]

    ==================================== ERRORS ====================================

    =================================== FAILURES ===================================

    =============================== warnings summary ===============================

    =========================== short test summary info ============================

    = 1 failed, 146194 passed, 7402 skipped, 1367 xfailed, 5 xpassed, 197 warnings, 10 errors in 1090.16s (0:18:10) =

This is just an example of what information is shown. You might see a slightly different result as what is shown above.

.. _install.dependencies:

Dependencies
------------

================================================================ ==========================
Package                                                          Minimum supported version
================================================================ ==========================
`NumPy <https://numpy.org>`__                                    1.19.5
`python-dateutil <https://dateutil.readthedocs.io/en/stable/>`__ 2.8.1
`pytz <https://pypi.org/project/pytz/>`__                        2020.1
================================================================ ==========================

.. _install.recommended_dependencies:

Recommended dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

* `numexpr <https://github.com/pydata/numexpr>`__: for accelerating certain numerical operations.
  ``numexpr`` uses multiple cores as well as smart chunking and caching to achieve large speedups.
  If installed, must be Version 2.7.3 or higher.

* `bottleneck <https://github.com/pydata/bottleneck>`__: for accelerating certain types of ``nan``
  evaluations. ``bottleneck`` uses specialized cython routines to achieve large speedups. If installed,
  must be Version 1.3.2 or higher.

.. note::

   You are highly encouraged to install these libraries, as they provide speed improvements, especially
   when working with large data sets.


.. _install.optional_dependencies:

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

pandas has many optional dependencies that are only used for specific methods.
For example, :func:`pandas.read_hdf` requires the ``pytables`` package, while
:meth:`DataFrame.to_markdown` requires the ``tabulate`` package. If the
optional dependency is not installed, pandas will raise an ``ImportError`` when
the method requiring that dependency is called.

Visualization
^^^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
matplotlib                3.3.2              Plotting library
Jinja2                    3.0.0              Conditional formatting with DataFrame.style
tabulate                  0.8.9              Printing in Markdown-friendly format (see `tabulate`_)
========================= ================== =============================================================

Computation
^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
SciPy                     1.7.1              Miscellaneous statistical functions
numba                     0.53.1             Alternative execution engine for rolling operations
                                             (see :ref:`Enhancing Performance <enhancingperf.numba>`)
xarray                    0.19.0             pandas-like API for N-dimensional data
========================= ================== =============================================================

Excel files
^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
xlrd                      2.0.1              Reading Excel
xlwt                      1.3.0              Writing Excel
xlsxwriter                1.4.3              Writing Excel
openpyxl                  3.0.7              Reading / writing for xlsx files
pyxlsb                    1.0.8              Reading for xlsb files
========================= ================== =============================================================

HTML
^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
BeautifulSoup4            4.9.3              HTML parser for read_html
html5lib                  1.1                HTML parser for read_html
lxml                      4.6.3              HTML parser for read_html
========================= ================== =============================================================

One of the following combinations of libraries is needed to use the
top-level :func:`~pandas.read_html` function:

* `BeautifulSoup4`_ and `html5lib`_
* `BeautifulSoup4`_ and `lxml`_
* `BeautifulSoup4`_ and `html5lib`_ and `lxml`_
* Only `lxml`_, although see :ref:`HTML Table Parsing <io.html.gotchas>`
  for reasons as to why you should probably **not** take this approach.

.. warning::

    * if you install `BeautifulSoup4`_ you must install either
      `lxml`_ or `html5lib`_ or both.
      :func:`~pandas.read_html` will **not** work with *only*
      `BeautifulSoup4`_ installed.
    * You are highly encouraged to read :ref:`HTML Table Parsing gotchas <io.html.gotchas>`.
      It explains issues surrounding the installation and
      usage of the above three libraries.

.. _html5lib: https://github.com/html5lib/html5lib-python
.. _BeautifulSoup4: https://www.crummy.com/software/BeautifulSoup
.. _lxml: https://lxml.de
.. _tabulate: https://github.com/astanin/python-tabulate

XML
^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
lxml                      4.5.0              XML parser for read_xml and tree builder for to_xml
========================= ================== =============================================================

SQL databases
^^^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
SQLAlchemy                1.4.16             SQL support for databases other than sqlite
psycopg2                  2.8.6              PostgreSQL engine for sqlalchemy
pymysql                   1.0.2              MySQL engine for sqlalchemy
========================= ================== =============================================================

Other data sources
^^^^^^^^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
PyTables                  3.6.1              HDF5-based reading / writing
blosc                     1.21.0             Compression for HDF5
zlib                                         Compression for HDF5
fastparquet               0.4.0              Parquet reading / writing
pyarrow                   1.0.1              Parquet, ORC, and feather reading / writing
pyreadstat                1.1.2              SPSS files (.sav) reading
========================= ================== =============================================================

.. _install.warn_orc:

.. warning::

    * If you want to use :func:`~pandas.read_orc`, it is highly recommended to install pyarrow using conda.
      The following is a summary of the environment in which :func:`~pandas.read_orc` can work.

      ========================= ================== =============================================================
      System                    Conda              PyPI
      ========================= ================== =============================================================
      Linux                     Successful         Failed(pyarrow==3.0 Successful)
      macOS                     Successful         Failed
      Windows                   Failed             Failed
      ========================= ================== =============================================================

Access data in the cloud
^^^^^^^^^^^^^^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
fsspec                    2021.5.0           Handling files aside from simple local and HTTP
gcsfs                     2021.5.0           Google Cloud Storage access
pandas-gbq                0.15.0             Google Big Query access
s3fs                      2021.05.0          Amazon S3 access
========================= ================== =============================================================

Clipboard
^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
PyQt4/PyQt5                                  Clipboard I/O
qtpy                                         Clipboard I/O
xclip                                        Clipboard I/O on linux
xsel                                         Clipboard I/O on linux
========================= ================== =============================================================


Compression
^^^^^^^^^^^

========================= ================== =============================================================
Dependency                Minimum Version    Notes
========================= ================== =============================================================
brotli                    0.7.0              Brotli compression
python-snappy             0.6.0              Snappy compression
Zstandard                 0.15.2             Zstandard compression
========================= ================== =============================================================
