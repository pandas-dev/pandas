<div align="center">
  <img src="https://github.com/pydata/pandas/blob/master/doc/logo/pandas_logo.png"><br>
</div>

-----------------

# pandas: powerful Python data analysis toolkit

| | |
------------ | -------------
Latest Release | [![Latest Version](https://img.shields.io/pypi/v/pandas.svg)](https://pypi.python.org/pypi/pandas/)
Latest Release | [![Latest Version](https://anaconda.org/pandas/pandas/badges/version.svg)](https://anaconda.org/pandas/pandas/)
Package Status | [![Package Status](https://img.shields.io/pypi/status/pandas.svg)](https://pypi.python.org/pypi/pandas/)
License | [![License](https://img.shields.io/pypi/l/pandas.svg)](https://pypi.python.org/pypi/pandas/)
Build Status | [![Build Status](https://travis-ci.org/pydata/pandas.svg?branch=master)](https://travis-ci.org/pydata/pandas)
Build Status | [![Build Status](https://ci.appveyor.com/api/projects/status/iblk29s98quexwxi/branch/master?svg=true)](https://ci.appveyor.com/project/jreback/pandas-465)
Coverage | [![Coverage](https://codecov.io/github/pydata/pandas/coverage.svg?branch=master)](https://codecov.io/github/pydata/pandas?branch=master)
Conda | [![Conda](http://pubbadges.s3-website-us-east-1.amazonaws.com/pkgs-downloads-pandas.png)](http://conda.anaconda.org/pandas)
PyPI | [![PyPI](https://img.shields.io/pypi/dm/pandas.svg)](https://pypi.python.org/pypi/pandas/)
Gitter | [![https://gitter.im/pydata/pandas](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/pydata/pandas)

## What is it

**pandas** is a Python package providing fast, flexible, and expressive data
structures designed to make working with "relational" or "labeled" data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, **real world** data analysis in Python. Additionally, it has
the broader goal of becoming **the most powerful and flexible open source data
analysis / manipulation tool available in any language**. It is already well on
its way toward this goal.

## Main Features
Here are just a few of the things that pandas does well:

  - Easy handling of [**missing data**][missing-data] (represented as
    `NaN`) in floating point as well as non-floating point data
  - Size mutability: columns can be [**inserted and
    deleted**][insertion-deletion] from DataFrame and higher dimensional
    objects
  - Automatic and explicit [**data alignment**][alignment]: objects can
    be explicitly aligned to a set of labels, or the user can simply
    ignore the labels and let `Series`, `DataFrame`, etc. automatically
    align the data for you in computations
  - Powerful, flexible [**group by**][groupby] functionality to perform
    split-apply-combine operations on data sets, for both aggregating
    and transforming data
  - Make it [**easy to convert**][conversion] ragged,
    differently-indexed data in other Python and NumPy data structures
    into DataFrame objects
  - Intelligent label-based [**slicing**][slicing], [**fancy
    indexing**][fancy-indexing], and [**subsetting**][subsetting] of
    large data sets
  - Intuitive [**merging**][merging] and [**joining**][joining] data
    sets
  - Flexible [**reshaping**][reshape] and [**pivoting**][pivot-table] of
    data sets
  - [**Hierarchical**][mi] labeling of axes (possible to have multiple
    labels per tick)
  - Robust IO tools for loading data from [**flat files**][flat-files]
    (CSV and delimited), [**Excel files**][excel], [**databases**][db],
    and saving/loading data from the ultrafast [**HDF5 format**][hdfstore]
  - [**Time series**][timeseries]-specific functionality: date range
    generation and frequency conversion, moving window statistics,
    moving window linear regressions, date shifting and lagging, etc.


   [missing-data]: http://pandas.pydata.org/pandas-docs/stable/missing_data.html#working-with-missing-data
   [insertion-deletion]: http://pandas.pydata.org/pandas-docs/stable/dsintro.html#column-selection-addition-deletion
   [alignment]: http://pandas.pydata.org/pandas-docs/stable/dsintro.html?highlight=alignment#intro-to-data-structures
   [groupby]: http://pandas.pydata.org/pandas-docs/stable/groupby.html#group-by-split-apply-combine
   [conversion]: http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe
   [slicing]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#slicing-ranges
   [fancy-indexing]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#advanced-indexing-with-ix
   [subsetting]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#boolean-indexing
   [merging]: http://pandas.pydata.org/pandas-docs/stable/merging.html#database-style-dataframe-joining-merging
   [joining]: http://pandas.pydata.org/pandas-docs/stable/merging.html#joining-on-index
   [reshape]: http://pandas.pydata.org/pandas-docs/stable/reshaping.html#reshaping-and-pivot-tables
   [pivot-table]: http://pandas.pydata.org/pandas-docs/stable/reshaping.html#pivot-tables-and-cross-tabulations
   [mi]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#hierarchical-indexing-multiindex
   [flat-files]: http://pandas.pydata.org/pandas-docs/stable/io.html#csv-text-files
   [excel]: http://pandas.pydata.org/pandas-docs/stable/io.html#excel-files
   [db]: http://pandas.pydata.org/pandas-docs/stable/io.html#sql-queries
   [hdfstore]: http://pandas.pydata.org/pandas-docs/stable/io.html#hdf5-pytables
   [timeseries]: http://pandas.pydata.org/pandas-docs/stable/timeseries.html#time-series-date-functionality

## Where to get it
The source code is currently hosted on GitHub at:
http://github.com/pydata/pandas

Binary installers for the latest released version are available at the [Python
package index](http://pypi.python.org/pypi/pandas/) and on conda.

```sh
# conda
conda install pandas
```

```sh
# or PyPI
pip install pandas
```

## Dependencies
- [NumPy](http://www.numpy.org): 1.7.0 or higher
- [python-dateutil](http://labix.org/python-dateutil): 1.5 or higher
- [pytz](http://pytz.sourceforge.net)
    - Needed for time zone support with ``pandas.date_range``

See the [full installation instructions](http://pandas.pydata.org/pandas-docs/stable/install.html#dependencies)
for recommended and optional dependencies.

## Installation from sources
To install pandas from source you need Cython in addition to the normal
dependencies above. Cython can be installed from pypi:

```sh
pip install cython
```

In the `pandas` directory (same one where you found this file after
cloning the git repo), execute:

```sh
python setup.py install
```

or for installing in [development mode](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs):

```sh
python setup.py develop
```

Alternatively, you can use `pip` if you want all the dependencies pulled
in automatically (the `-e` option is for installing it in [development
mode](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs)):

```sh
pip install -e .
```

On Windows, you will need to install MinGW and execute:

```sh
python setup.py build --compiler=mingw32
python setup.py install
```

See http://pandas.pydata.org/ for more information.

## License
BSD

## Documentation
The official documentation is hosted on PyData.org: http://pandas.pydata.org/

The Sphinx documentation should provide a good starting point for learning how
to use the library. Expect the docs to continue to expand as time goes on.

## Background
Work on ``pandas`` started at AQR (a quantitative hedge fund) in 2008 and
has been under active development since then.

## Discussion and Development
Since pandas development is related to a number of other scientific
Python projects, questions are welcome on the scipy-user mailing
list. Specialized discussions or design issues should take place on
the PyData mailing list / Google group:

https://groups.google.com/forum/#!forum/pydata
