.. _install:

{{ header }}

============
Installation
============

The pandas development team officially distributes pandas for installation
through the following methods:

* Available on `conda-forge <https://anaconda.org/conda-forge/pandas>`__ for installation with the conda package manager.
* Available on `PyPI <https://pypi.org/project/pandas/>`__ for installation with pip.
* Available on `Github <https://github.com/pandas-dev/pandas>`__ for installation from source.

.. note::
    pandas may be installable from other sources besides the ones listed above,
    but they are **not** managed by the pandas development team.

.. _install.version:

Python version support
----------------------

See :ref:`Python support policy <policies.python_support>`.

Installing pandas
-----------------

.. _install.conda:

Installing with Conda
~~~~~~~~~~~~~~~~~~~~~

For users working with the `Conda <https://conda.io/en/latest/>`__ package manager,
pandas can be installed from the ``conda-forge`` channel.

.. code-block:: shell

    conda install -c conda-forge pandas

To install the Conda package manager on your system, the
`Miniforge distribution <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`__
is recommended.

Additionally, it is recommended to install and run pandas from a virtual environment.

.. code-block:: shell

    conda create -c conda-forge -n name_of_my_env python pandas
    # On Linux or MacOS
    source activate name_of_my_env
    # On Windows
    activate name_of_my_env

.. tip::
    For users that are new to Python, the easiest way to install Python, pandas, and the
    packages that make up the `PyData <https://pydata.org/>`__ stack such as
    `SciPy <https://scipy.org/>`__, `NumPy <https://numpy.org/>`__ and
    `Matplotlib <https://matplotlib.org/>`__
    is with `Anaconda <https://docs.anaconda.com/anaconda/install/>`__, a cross-platform
    (Linux, macOS, Windows) Python distribution for data analytics and
    scientific computing.

    However, pandas from Anaconda is **not** officially managed by the pandas development team.

.. _install.pip:

Installing with pip
~~~~~~~~~~~~~~~~~~~

For users working with the `pip <https://pip.pypa.io/en/stable/>`__ package manager,
pandas can be installed from `PyPI <https://pypi.org/project/pandas/>`__.

.. code-block:: shell

    pip install pandas

pandas can also be installed with sets of optional dependencies to enable certain functionality. For example,
to install pandas with the optional dependencies to read Excel files.

.. code-block:: shell

    pip install "pandas[excel]"

The full list of extras that can be installed can be found in the :ref:`dependency section.<install.optional_dependencies>`

Additionally, it is recommended to install and run pandas from a virtual environment, for example,
using the Python standard library's `venv <https://docs.python.org/3/library/venv.html>`__

.. _install.source:

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

See the :ref:`contributing guide <contributing>` for complete instructions on building from the git source tree.
Further, see :ref:`creating a development environment <contributing_environment>` if you wish to create
a pandas development environment.

.. _install.dev:

Installing the development version of pandas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing the development version is the quickest way to:

* Try a new feature that will be shipped in the next release (that is, a feature from a pull-request that was recently merged to the main branch).
* Check whether a bug you encountered has been fixed since the last release.

The development version is usually uploaded daily to the scientific-python-nightly-wheels
index from the PyPI registry of anaconda.org. You can install it by running.

.. code-block:: shell

    pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple pandas

.. note::
    You might be required to uninstall an existing version of pandas to install the development version.

    .. code-block:: shell

        pip uninstall pandas -y

Running the test suite
----------------------

If pandas has been installed :ref:`from source <install.source>`, running ``pytest pandas`` will run all of pandas unit tests.

The unit tests can also be run from the pandas module itself with the :func:`test` function. The packages required to run the tests
can be installed with ``pip install "pandas[test]"``.

.. note::

    Test failures are not necessarily indicative of a broken pandas installation.

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
`NumPy <https://numpy.org>`__                                    1.23.5
`python-dateutil <https://dateutil.readthedocs.io/en/stable/>`__ 2.8.2
`tzdata <https://pypi.org/project/tzdata/>`__                    2022.7
================================================================ ==========================

.. _install.optional_dependencies:

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

pandas has many optional dependencies that are only used for specific methods.
For example, :func:`pandas.read_hdf` requires the ``pytables`` package, while
:meth:`DataFrame.to_markdown` requires the ``tabulate`` package. If the
optional dependency is not installed, pandas will raise an ``ImportError`` when
the method requiring that dependency is called.

With pip, optional pandas dependencies can be installed or managed in a file (e.g. requirements.txt or pyproject.toml)
as optional extras (e.g. ``pandas[performance, aws]``). All optional dependencies can be installed with ``pandas[all]``,
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
`numexpr <https://github.com/pydata/numexpr>`__       2.8.4              performance        Accelerates certain numerical operations by using multiple cores as well as smart chunking and caching to achieve large speedups
`bottleneck <https://github.com/pydata/bottleneck>`__ 1.3.6              performance        Accelerates certain types of ``nan`` by using specialized cython routines to achieve large speedup.
`numba <https://github.com/numba/numba>`__            0.56.4             performance        Alternative execution engine for operations that accept ``engine="numba"`` using a JIT compiler that translates Python functions to optimized machine code using the LLVM compiler.
===================================================== ================== ================== ===================================================================================================================================================================================

Visualization
^^^^^^^^^^^^^

Installable with ``pip install "pandas[plot, output-formatting]"``.

========================= ================== ================== =============================================================
Dependency                Minimum Version    pip extra          Notes
========================= ================== ================== =============================================================
matplotlib                3.6.3              plot               Plotting library
Jinja2                    3.1.2              output-formatting  Conditional formatting with DataFrame.style
tabulate                  0.9.0              output-formatting  Printing in Markdown-friendly format (see `tabulate`_)
========================= ================== ================== =============================================================

Computation
^^^^^^^^^^^

Installable with ``pip install "pandas[computation]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
SciPy                     1.10.0             computation     Miscellaneous statistical functions
xarray                    2022.12.0          computation     pandas-like API for N-dimensional data
========================= ================== =============== =============================================================

.. _install.excel_dependencies:

Excel files
^^^^^^^^^^^

Installable with ``pip install "pandas[excel]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
xlrd                      2.0.1              excel           Reading for xls files
xlsxwriter                3.0.5              excel           Writing for xlsx files
openpyxl                  3.1.0              excel           Reading / writing for Excel 2010 xlsx/xlsm/xltx/xltm files
pyxlsb                    1.0.10             excel           Reading for xlsb files
python-calamine           0.1.7              excel           Reading for xls/xlsx/xlsm/xlsb/xla/xlam/ods files
odfpy                     1.4.1              excel           Reading / writing for OpenDocument 1.2 files
========================= ================== =============== =============================================================

HTML
^^^^

Installable with ``pip install "pandas[html]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
BeautifulSoup4            4.11.2             html            HTML parser for read_html
html5lib                  1.1                html            HTML parser for read_html
lxml                      4.9.2              html            HTML parser for read_html
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
lxml                      4.9.2              xml             XML parser for read_xml and tree builder for to_xml
========================= ================== =============== =============================================================

SQL databases
^^^^^^^^^^^^^

Traditional drivers are installable with ``pip install "pandas[postgresql, mysql, sql-other]"``

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
SQLAlchemy                2.0.0              postgresql,     SQL support for databases other than sqlite
                                             mysql,
                                             sql-other
psycopg2                  2.9.6              postgresql      PostgreSQL engine for sqlalchemy
pymysql                   1.0.2              mysql           MySQL engine for sqlalchemy
adbc-driver-postgresql    0.10.0             postgresql      ADBC Driver for PostgreSQL
adbc-driver-sqlite        0.8.0              sql-other       ADBC Driver for SQLite
========================= ================== =============== =============================================================

Other data sources
^^^^^^^^^^^^^^^^^^

Installable with ``pip install "pandas[hdf5, parquet, feather, spss, excel]"``

========================= ================== ================ =============================================================
Dependency                Minimum Version    pip extra        Notes
========================= ================== ================ =============================================================
PyTables                  3.8.0              hdf5             HDF5-based reading / writing
blosc                     1.21.3             hdf5             Compression for HDF5; only available on ``conda``
zlib                                         hdf5             Compression for HDF5
fastparquet               2023.10.0          -                Parquet reading / writing (pyarrow is default)
pyarrow                   10.0.1             parquet, feather Parquet, ORC, and feather reading / writing
pyreadstat                1.2.0              spss             SPSS files (.sav) reading
odfpy                     1.4.1              excel            Open document format (.odf, .ods, .odt) reading / writing
========================= ================== ================ =============================================================

.. _install.warn_orc:

.. warning::

    * If you want to use :func:`~pandas.read_orc`, it is highly recommended to install pyarrow using conda.
      :func:`~pandas.read_orc` may fail if pyarrow was installed from pypi, and :func:`~pandas.read_orc` is
      not compatible with Windows OS.

Access data in the cloud
^^^^^^^^^^^^^^^^^^^^^^^^

Installable with ``pip install "pandas[fss, aws, gcp]"``

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
fsspec                    2022.11.0          fss, gcp, aws   Handling files aside from simple local and HTTP (required
                                                             dependency of s3fs, gcsfs).
gcsfs                     2022.11.0          gcp             Google Cloud Storage access
s3fs                      2022.11.0          aws             Amazon S3 access
========================= ================== =============== =============================================================

Clipboard
^^^^^^^^^

Installable with ``pip install "pandas[clipboard]"``.

========================= ================== =============== =============================================================
Dependency                Minimum Version    pip extra       Notes
========================= ================== =============== =============================================================
PyQt4/PyQt5               5.15.9             clipboard       Clipboard I/O
qtpy                      2.3.0              clipboard       Clipboard I/O
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
Zstandard                 0.19.0             compression     Zstandard compression
========================= ================== =============== =============================================================

Timezone
^^^^^^^^

Installable with ``pip install "pandas[timezone]"``

========================= ================== =================== =============================================================
Dependency                Minimum Version    pip extra           Notes
========================= ================== =================== =============================================================
pytz                      2023.4             timezone            Alternative timezone library to ``zoneinfo``.
========================= ================== =================== =============================================================
