Installation from sources
=========================

In the pandas directory (same one where you found this file), execute:

python setup.py install

On Windows, you will need to install MinGW and execute

python setup.py install --compiler=mingw32

See

http://pandas.sourceforge.net/

For more information.

=============
Release Notes
=============

What it is
==========

pandas is a library for pan-el da-ta analysis, i.e. multidimensional
time series and cross-sectional data sets commonly found in
statistics, econometrics, or finance. It provides convenient and
easy-to-understand NumPy-based data structures for generic labeled
data, with focus on automatically aligning data based on its label(s)
and handling missing observations. One major goal of the library is to
simplify the implementation of statistical models on unreliable data.

Main Features
=============

* Data structures: for 1, 2, and 3 dimensional labeled data
  sets. Some of their main features include:

    * Automatically aligning data
    * Handling missing observations in calculations
    * Convenient slicing and reshaping ("reindexing") functions
    * Provide 'group by' aggregation or transformation functionality
    * Tools for merging / joining together data sets
    * Simple matplotlib integration for plotting

* Date tools: objects for expressing date offsets or generating date
  ranges; some functionality similar to scikits.timeseries

* Statistical models: convenient ordinary least squares and panel OLS
  implementations for in-sample or rolling time series /
  cross-sectional regressions. These will hopefully be the starting
  point for implementing other models

pandas is not necessarily intended as a standalone library but rather
as something which can be used in tandem with other NumPy-based
packages like scikits.statsmodels. Where possible wheel-reinvention
has largely been avoided. Also, its time series manipulation
capability is not as extensive as scikits.timeseries; pandas does have
its own time series object which fits into the unified data model.

Some other useful tools for time series data (moving average, standard
deviation, etc.) are available in the codebase but do not yet have a
convenient interface. These will be highlighted in a future release.

Where to get it
===============

The source code is currently hosted on googlecode at:

http://pandas.googlecode.com

Binary releases can be downloaded there, or alternately via the Python
package index or easy_install

PyPi: http://pypi.python.org/pypi/pandas/

License
=======

BSD

Documentation
=============

The official documentation is hosted on SourceForge.

http://pandas.sourceforge.net/

The sphinx documentation is still in an incomplete state, but it
should provide a good starting point for learning how to use the
library. Expect the docs to continue to expand as time goes on.

Background
==========

Work on pandas started at AQR (a quantitative hedge fund) in 2008 and
has been under active development since then.

Discussion and Development
==========================

Since pandas development is related to a number of other scientific
Python projects, questions are welcome on the scipy-user mailing
list. Specialized discussions or design issues should take place on
the pystatsmodels mailing list / google group, where
scikits.statsmodels and other libraries will also be discussed:

http://groups.google.com/group/pystatsmodels
