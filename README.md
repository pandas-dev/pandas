# pandas: powerful Python data analysis toolkit

![Travis-CI Build Status](https://travis-ci.org/pydata/pandas.svg)

[![Scatter-CI Status page](http://scatterci.github.io/scatterci48.jpg)](http://scatterci.github.io/pydata/pandas)

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

Binary installers for the latest released version are available at the Python
package index

    http://pypi.python.org/pypi/pandas/

And via `easy_install`:

```sh
easy_install pandas
```

or  `pip`:

```sh
pip install pandas
```

## Dependencies
- [NumPy](http://www.numpy.org): 1.6.1 or higher
- [python-dateutil](http://labix.org/python-dateutil): 1.5 or higher
- [pytz](http://pytz.sourceforge.net)
    - Needed for time zone support with ``pandas.date_range``

### Highly Recommended Dependencies
- [numexpr](http://code.google.com/p/numexpr/)
   - Needed to accelerate some expression evaluation operations
   - Required by PyTables
- [bottleneck](http://berkeleyanalytics.com/bottleneck)
   - Needed to accelerate certain numerical operations

### Optional dependencies
- [Cython](http://www.cython.org): Only necessary to build development version. Version 0.17.1 or higher.
- [SciPy](http://www.scipy.org): miscellaneous statistical functions
- [PyTables](http://www.pytables.org): necessary for HDF5-based storage
- [SQLAlchemy](http://www.sqlalchemy.org): for SQL database support. Version 0.8.1 or higher recommended.
- [matplotlib](http://matplotlib.sourceforge.net/): for plotting
- [statsmodels](http://statsmodels.sourceforge.net/)
   - Needed for parts of `pandas.stats`
- For Excel I/O:
  - [xlrd/xlwt](http://www.python-excel.org/)
     - Excel reading (xlrd) and writing (xlwt)
  - [openpyxl](http://packages.python.org/openpyxl/)
     - openpyxl version 1.6.1 or higher, but lower than 2.0.0, for
       writing .xlsx files
     - xlrd >= 0.9.0
  - [XlsxWriter](https://pypi.python.org/pypi/XlsxWriter)
     - Alternative Excel writer.
- [Google bq Command Line Tool](https://developers.google.com/bigquery/bq-command-line-tool/)
  - Needed for `pandas.io.gbq`
- [boto](https://pypi.python.org/pypi/boto): necessary for Amazon S3 access.
- One of the following combinations of libraries is needed to use the
  top-level [`pandas.read_html`][read-html-docs] function:
  - [BeautifulSoup4][BeautifulSoup4] and [html5lib][html5lib] (Any
    recent version of [html5lib][html5lib] is okay.)
  - [BeautifulSoup4][BeautifulSoup4] and [lxml][lxml]
  - [BeautifulSoup4][BeautifulSoup4] and [html5lib][html5lib] and [lxml][lxml]
  - Only [lxml][lxml], although see [HTML reading gotchas][html-gotchas]
    for reasons as to why you should probably **not** take this approach.

#### Notes about HTML parsing libraries
- If you install [BeautifulSoup4][BeautifulSoup4] you must install
  either [lxml][lxml] or [html5lib][html5lib] or both.
  `pandas.read_html` will **not** work with *only* `BeautifulSoup4`
  installed.
- You are strongly encouraged to read [HTML reading
  gotchas][html-gotchas]. It explains issues surrounding the
  installation and usage of the above three libraries.
- You may need to install an older version of
  [BeautifulSoup4][BeautifulSoup4]:
    - Versions 4.2.1, 4.1.3 and 4.0.2 have been confirmed for 64 and
      32-bit Ubuntu/Debian
- Additionally, if you're using [Anaconda][Anaconda] you should
  definitely read [the gotchas about HTML parsing][html-gotchas]
  libraries
- If you're on a system with `apt-get` you can do

  ```sh
  sudo apt-get build-dep python-lxml
  ```

  to get the necessary dependencies for installation of [lxml][lxml].
  This will prevent further headaches down the line.

   [html5lib]: https://github.com/html5lib/html5lib-python "html5lib"
   [BeautifulSoup4]: http://www.crummy.com/software/BeautifulSoup "BeautifulSoup4"
   [lxml]: http://lxml.de
   [Anaconda]: https://store.continuum.io/cshop/anaconda
   [NumPy]: http://numpy.scipy.org/
   [html-gotchas]: http://pandas.pydata.org/pandas-docs/stable/gotchas.html#html-table-parsing
   [read-html-docs]: http://pandas.pydata.org/pandas-docs/stable/generated/pandas.io.html.read_html.html#pandas.io.html.read_html

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

or for installing in [development mode](http://www.pip-installer.org/en/latest/usage.html):

```sh
python setup.py develop
```

Alternatively, you can use `pip` if you want all the dependencies pulled
in automatically (the `-e` option is for installing it in [development
mode](http://www.pip-installer.org/en/latest/usage.html)):

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
the pystatsmodels mailing list / Google group, where
``scikits.statsmodels`` and other libraries will also be discussed:

http://groups.google.com/group/pystatsmodels
