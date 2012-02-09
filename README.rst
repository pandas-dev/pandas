=============================================
pandas: powerful Python data analysis toolkit
=============================================

What is it
==========

**pandas** is a Python package providing fast, flexible, and expressive data
structures designed to make working with "relational" or "labeled" data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, **real world** data analysis in Python. Additionally, it has
the broader goal of becoming **the most powerful and flexible open source data
analysis / manipulation tool available in any language**. It is already well on
its way toward this goal.

Main Features
=============

Here are just a few of the things that pandas does well:

  - Easy handling of **missing data** (represented as NaN) in floating point as
    well as non-floating point data
  - Size mutability: columns can be **inserted and deleted** from DataFrame and
    higher dimensional objects
  - Automatic and explicit **data alignment**: objects can be explicitly
    aligned to a set of labels, or the user can simply ignore the labels and
    let `Series`, `DataFrame`, etc. automatically align the data for you in
    computations
  - Powerful, flexible **group by** functionality to perform
    split-apply-combine operations on data sets, for both aggregating and
    transforming data
  - Make it **easy to convert** ragged, differently-indexed data in other
    Python and NumPy data structures into DataFrame objects
  - Intelligent label-based **slicing**, **fancy indexing**, and **subsetting**
    of large data sets
  - Intuitive **merging** and **joining** data sets
  - Flexible **reshaping** and pivoting of data sets
  - **Hierarchical** labeling of axes (possible to have multiple labels per
    tick)
  - Robust IO tools for loading data from **flat files** (CSV and delimited),
    Excel files, databases, and saving / loading data from the ultrafast **HDF5
    format**
  - **Time series**-specific functionality: date range generation and frequency
    conversion, moving window statistics, moving window linear regressions,
    date shifting and lagging, etc.

Where to get it
===============

The source code is currently hosted on GitHub at: http://github.com/pydata/pandas

Binary installers for the latest released version are available at the Python
package index::

    http://pypi.python.org/pypi/pandas/

And via ``easy_install`` or ``pip``::

    easy_install pandas
    pip install pandas

Dependencies
============

  * `NumPy <http://www.numpy.org>`__: 1.6.1 or higher. Older versions will work
    but may not pass all of the unit tests. Bare minimum is NumPy 1.4.0.
  * `python-dateutil <http://labix.org/python-dateutil>`__ 1.5

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

  * `Cython <http://www.cython.org>`__: Only necessary to build development
    version
  * `SciPy <http://www.scipy.org>`__: miscellaneous statistical functions
  * `PyTables <http://www.pytables.org>`__: necessary for HDF5-based storage
  * `matplotlib <http://matplotlib.sourceforge.net/>`__: for plotting
  * `scikits.statsmodels <http://statsmodels.sourceforge.net/>`__
     * Needed for parts of :mod:`pandas.stats`
  * `pytz <http://pytz.sourceforge.net/>`__
     * Needed for time zone support with ``DateRange``

Installation from sources
=========================

In the ``pandas`` directory (same one where you found this file), execute::

    python setup.py install

On Windows, you will need to install MinGW and execute::

    python setup.py build --compiler=mingw32
    python setup.py install

See http://pandas.pydata.org/ for more information.

License
=======

BSD

Documentation
=============

The official documentation is hosted on PyData.org: http://pandas.pydata.org/

The Sphinx documentation should provide a good starting point for learning how
to use the library. Expect the docs to continue to expand as time goes on.

Background
==========

Work on ``pandas`` started at AQR (a quantitative hedge fund) in 2008 and
has been under active development since then.

Discussion and Development
==========================

Since ``pandas`` development is related to a number of other scientific
Python projects, questions are welcome on the scipy-user mailing
list. Specialized discussions or design issues should take place on
the pystatsmodels mailing list / Google group, where
``scikits.statsmodels`` and other libraries will also be discussed:

http://groups.google.com/group/pystatsmodels

  .. _NumPy: http://numpy.scipy.org/
