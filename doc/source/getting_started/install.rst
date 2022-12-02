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

Officially Python 3.8, 3.9, 3.10 and 3.11.

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

pandas can also be installed with sets of optional dependencies to enable certain functionality. For example,
to install pandas with the optional dependencies to read Excel files.

::

    pip install "pandas[excel]"


The full list of extras that can be installed can be found in the :ref:`dependency section.<install.optional_dependencies>`

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

.. _install.required_dependencies:

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

pandas requires the following dependencies.

================================================================ ==========================
Package                                                          Minimum supported version
================================================================ ==========================
`NumPy <https://numpy.org>`__                                    1.20.3
`python-dateutil <https://dateutil.readthedocs.io/en/stable/>`__ 2.8.2
`pytz <https://pypi.org/project/pytz/>`__                        2020.1
================================================================ ==========================

.. _install.optional_dependencies:

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

pandas has many optional dependencies that are only used for specific methods.
For example, :func:`pandas.read_hdf` requires the ``pytables`` package, while
:meth:`DataFrame.to_markdown` requires the ``tabulate`` package. If the
optional dependency is not installed, pandas will raise an ``ImportError`` when
the method requiring that dependency is called.

If using pip, optional pandas dependencies can be installed or managed in a file (e.g. requirements.txt or pyproject.toml)
as optional extras (e.g.,``pandas[performance, aws]>=1.5.0``). All optional dependencies can be installed with ``pandas[all]``,
and specific sets of dependencies are listed in the sections below.

.. _install.recommended_dependencies:

Performance dependencies (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   You are highly encouraged to install these libraries, as they provide speed improvements, especially
   when working with large data sets.

Installable with ``pip install "pandas[performance]"``

===================================================== ================== ================== ===================================================================================================================================================================================
Dependency                                            Minimum Version    pip extra          Notes
===================================================== ================== ================== ===================================================================================================================================================================================
`numexpr <https://github.com/pydata/numexpr>`__       2.7.3              performance        Accelerates certain numerical operations by using uses multiple cores as well as smart chunking and caching to achieve large speedups
`bottleneck <https://github.com/pydata/bottleneck>`__ 1.3.2              performance        Accelerates certain types of ``nan`` by using specialized cython routines to achieve large speedup.
`numba <https://github.com/numba/numba>`__            0.53.1             performance        Alternative execution engine for operations that accept ``engine="numba"`` using a JIT compiler that translates Python functions to optimized machine code using the LLVM compiler.
===================================================== ================== ================== ===================================================================================================================================================================================

Timezones
^^^^^^^^^

Installable with ``pip install "pandas[timezone]"``

========================= ========================= =============== =============================================================
Dependency                Minimum Version           pip extra       Notes
========================= ========================= =============== =============================================================
tzdata                    2022.1(pypi)/             timezone        Allows the use of ``zoneinfo`` timezones with pandas.
                          2022a(for system tzdata)                  **Note**: You only need to install the pypi package if your
                                                                    system does not already provide the IANA tz database.
                                                                    However, the minimum tzdata version still applies, even if it
                                                                    is not enforced through an error.

                                                                    If you would like to keep your system tzdata version updated,
                                                                    it is recommended to use the ``tzdata`` package from
                                                                    conda-forge.
========================= ========================= =============== =============================================================

Visualization
^^^^^^^^^^^^^

Installable with ``pip install "pandas[plot, output_formatting]"``.

========================= ================== ================== =============================================================
Dependency                Minimum Version    pip extra          Notes
========================= ================== ================== =============================================================
matplotlib                3.6.1              plot               Plotting library
Jinja2                    3.0.0              output_formatting  Conditional formatting with DataFrame.style
tabulate                  0.8.9              output_formatting  Printing in Markdown-friendly format (see `tabulate`_)
========================= ================== ================== =============================================================

Computation
^^^^^^^^^^^

Installable with ``pip install "pandas[computation]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
SciPy                     1.7.1              computation     Miscellaneous statistical functions
xarray                    0.19.0             computation     pandas-like API for N-dimensional data
========================= ================== =============== =============================================================

Excel files
^^^^^^^^^^^

Installable with ``pip install "pandas[excel]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
xlrd                      2.0.1              excel           Reading Excel
xlsxwriter                1.4.3              excel           Writing Excel
openpyxl                  3.0.7              excel           Reading / writing for xlsx files
pyxlsb                    1.0.8              excel           Reading for xlsb files
========================= ================== =============== =============================================================

HTML
^^^^

Installable with ``pip install "pandas[html]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
BeautifulSoup4            4.9.3              html            HTML parser for read_html
html5lib                  1.1                html            HTML parser for read_html
lxml                      4.6.3              html            HTML parser for read_html
========================= ================== =============== =============================================================

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

Installable with ``pip install "pandas[xml]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
lxml                      4.6.3              xml             XML parser for read_xml and tree builder for to_xml
========================= ================== =============== =============================================================

SQL databases
^^^^^^^^^^^^^

Installable with ``pip install "pandas[postgresql, mysql, sql-other]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
SQLAlchemy                1.4.16             postgresql,     SQL support for databases other than sqlite
                                             mysql,
                                             sql-other
psycopg2                  2.8.6              postgresql      PostgreSQL engine for sqlalchemy
pymysql                   1.0.2              mysql           MySQL engine for sqlalchemy
========================= ================== =============== =============================================================

Other data sources
^^^^^^^^^^^^^^^^^^

Installable with ``pip install "pandas[hdf5, parquet, feather, spss, excel]"``

========================= ================== ================ =============================================================
Dependency                Minimum Version    pip extra        Notes
========================= ================== ================ =============================================================
PyTables                  3.6.1              hdf5             HDF5-based reading / writing
blosc                     1.21.0             hdf5             Compression for HDF5; only available on ``conda``
zlib                                         hdf5             Compression for HDF5
fastparquet               0.6.3              -                Parquet reading / writing (pyarrow is default)
pyarrow                   6.0.0              parquet, feather Parquet, ORC, and feather reading / writing
pyreadstat                1.1.2              spss             SPSS files (.sav) reading
odfpy                     1.4.1              excel            Open document format (.odf, .ods, .odt) reading / writing
========================= ================== ================ =============================================================

.. _install.warn_orc:

.. warning::

    * If you want to use :func:`~pandas.read_orc`, it is highly recommended to install pyarrow using conda.
      The following is a summary of the environment in which :func:`~pandas.read_orc` can work.

      ========================= ================== =============================================================
      System                    Conda              PyPI
      ========================= ================== =============================================================
      Linux                     Successful         Failed
      macOS                     Successful         Failed
      Windows                   Failed             Failed
      ========================= ================== =============================================================

Access data in the cloud
^^^^^^^^^^^^^^^^^^^^^^^^

Installable with ``pip install "pandas[fss, aws, gcp]"``

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
fsspec                    2021.7.0           fss, gcp, aws   Handling files aside from simple local and HTTP (required
                                                             dependency of s3fs, gcsfs).
gcsfs                     2021.7.0           gcp             Google Cloud Storage access
pandas-gbq                0.15.0             gcp             Google Big Query access
s3fs                      2021.08.0          aws             Amazon S3 access
========================= ================== =============== =============================================================

Clipboard
^^^^^^^^^

Installable with ``pip install "pandas[clipboard]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
PyQt4/PyQt5               5.15.1             clipboard       Clipboard I/O
qtpy                      2.2.0              clipboard       Clipboard I/O
========================= ================== =============== =============================================================

.. note::

   Depending on operating system, system-level packages may need to installed.
   For clipboard to operate on Linux one of the CLI tools ``xclip`` or ``xsel`` must be installed on your system.


Compression
^^^^^^^^^^^

Installable with ``pip install "pandas[compression]"``

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
brotli                    0.7.0              compression     Brotli compression
python-snappy             0.6.0              compression     Snappy compression
Zstandard                 0.15.2             compression     Zstandard compression
========================= ================== =============== =============================================================
