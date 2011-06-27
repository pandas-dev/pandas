
************************
pandas 0.4 Release Notes
************************

==========
What is it
==========

**pandas** is a library of powerful labeled data structures, statistical tools,
and general code for working with time series and cross-sectional data. It was
designed with the practical needs of statistical modeling and large,
inhomogeneous data sets in mind. It is particularly well suited for, among other
things, financial data analysis applications.

===============
Where to get it
===============

Source code: http://github.com/wesm/pandas
Binary installers on PyPI: http://pypi.python.org/pypi/pandas
Documentation: http://pandas.sourceforge.net

=============
Release notes
=============

**Release date:** NOT YET RELEASED

**New features / modules**

* `pandas.core.sparse` module: "Sparse" (mostly-NA, or some other fill value)
  versions of `Series`, `DataFrame`, and `WidePanel`. For low-density data, this
  will result in significant performance boosts, and smaller memory
  footprint. Added `to_sparse` methods to `Series`, `DataFrame`, and
  `WidePanel`. See online documentation for more on these
* `Series.describe`, `DataFrame.describe`: produces an R-like table of summary
  statistics about each data column
* `DataFrame.quantile`, `Series.quantile`
* Fancy indexing operator on Series / DataFrame, e.g. via .ix operator. Both
  getting and setting of values is supported; however, setting values will only
  currently work on homogeneously-typed DataFrame objects
  * series.ix[[d1, d2, d3]]
  * frame.ix[5:10, ['C', 'B', 'A']], frame.ix[5:10, 'A':'C']
  * frame.ix[date1:date2]
* `Series` arithmetic methods with optional fill_value for missing data,
  e.g. a.add(b, fill_value=0). If a location is missing for both it will still
  be missing in the result though.
* Boolean indexing with `DataFrame` objects: data[data > 0.1] = 0.1
* `pytz` / tzinfo support in `DateRange`
  * `tz_localize`, `tz_normalize`, and `tz_validate` methods added
* Added `ExcelFile` class to `pandas.io.parsers` for parsing multiple sheets out
  of a single Excel 2003 document
* `GroupBy` aggregations can now optionally *broadcast*, e.g. produce an object
  of the same size with the aggregated value propagated
* Added `select` function in all data structures: reindex axis based on
  arbitrary criterion (function returning boolean value),
  e.g. frame.select(lambda x: 'foo' in x, axis=1)
* `DataFrame.consolidate` method, API function relating to redesigned internals
* `HDFStore` class in `pandas.io.pytables` has been largely rewritten using
  patches from Jeff Reback from others. It now supports mixed-type `DataFrame`
  and `Series` data and can store `WidePanel` objects. It also has the option to
  query `DataFrame` and `WidePanel` data. Loading data from legacy `HDFStore`
  files is supported explicitly in the code

**Improvements to existing features**

* The 2-dimensional `DataFrame` and `DataMatrix` classes have been extensively
  redesigned internally into a single class `DataFrame`, preserving where
  possible their optimal performance characteristics. This should reduce
  confusion from users about which class to use.
  * Note that under ther hood there is a new essentially "lazy evaluation"
    scheme within respect to adding columns to DataFrame. During some
    operations, like-typed blocks will be "consolidated" but not before.
* Column ordering for mixed type data is now completely consistent in
  `DataFrame`. In prior releases, there was inconsistent column ordering in
  `DataMatrix`
* Improved console / string formatting of DataMatrix with negative numbers
* Added `skiprows` and `na_values` arguments to `pandas.io.parsers` functions
  for more flexible IO
* Can slice `DataFrame` and get a view of the data (when homogeneously typed),
  e.g. frame.xs(idx, copy=False) or frame.ix[idx]
* Many speed optimizations throughout `Series` and `DataFrame`
* Eager evaluation of groups when calling ``groupby`` functions, so if there is
  an exception with the grouping function it will raised immediately versus
  sometime later on when the groups are needed
* `datetools.WeekOfMonth` offset can be parameterized with `n` different than 1
  or -1.
* `parseCSV` / `read_csv` functions and others in `pandas.io.parsers` now can
  take a list of custom NA values, and also a list of rows to skip
* Statistical methods on DataFrame like `mean`, `std`, `var`, `skew` will now
  ignore non-numerical data. Before a not very useful error message was
  generated. A flag `numeric_only` has been added to `DataFrame.sum` and
  `DataFrame.count` to enable this behavior in those methods if so desired
  (disabled by default)
* `DataFrame.pivot` generalized to enable pivoting multiple columns into a
  `WidePanel`

**API Changes**

* The `DataMatrix` variable now refers to `DataFrame`, will be removed within
  two releases
* A `copy` argument has been added to the `DataFrame` constructor to avoid
  unnecessary copying of data. Data is no longer copied by default when passed
  into the constructor
* Handling of boolean dtype in `DataFrame` has been improved to support storage
  of boolean data with NA / NaN values. Before it was being converted to float64
  so this should not (in theory) cause API breakage
* Boolean indexing using Series must now have the same indices (labels)
* Backwards compatibility support for begin/end/nPeriods keyword arguments in
  DateRange class has been removed
* More intuitive / shorter filling aliases `ffill` (for `pad`) and `bfill` (for
  `backfill`) have been added to the functions that use them: `reindex`,
  `asfreq`, `fillna`.
* `pandas.core.mixins` code moved to `pandas.core.generic`
* `buffer` keyword arguments (e.g. `DataFrame.toString`) renamed to `buf` to
  avoid using Python built-in name
* `DataFrame.rows()` removed (use `DataFrame.index`)
* Added deprecation warning to `DataFrame.cols()`, to be removed in next release
* `DataFrame` deprecations and de-camelCasing: `merge`, `asMatrix`,
  `toDataMatrix`, `_firstTimeWithValue`, `_lastTimeWithValue`, `toRecords`,
  `fromRecords`
* `pandas.io.parsers` method deprecations
  * `parseCSV` is now `read_csv` and keyword arguments have been de-camelCased
  * `parseText` is now `read_table`
  * `parseExcel` is replaced by the `ExcelFile` class and its `parse` method
* `fillMethod` arguments (deprecated in prior release) removed, should be
  replaced with `method`
* `Series.fill`, `DataFrame.fill`, and `WidePanel.fill` removed, use `fillna`
  instead

**Bug fixes**

* Column ordering in `pandas.io.parsers.parseCSV` will match CSV in the presence
  of mixed-type data
* Fixed handling of Excel 2003 dates in `pandas.io.parsers`
* `DateRange` caching was happening with high resolution `DateOffset` objects,
  e.g. `DateOffset(seconds=1)`. This has been fixed
* Fixed __truediv__ issue in `DataFrame`
* Fixed `DataFrame.toCSV` bug preventing IO round trips in some cases
* Fixed bug in `Series.plot` causing matplotlib to barf in exceptional cases
* Disabled `Index` objects from being hashable, like ndarrays
* Added `__ne__` implementation to `Index` so that operations like ts[ts != idx]
  will work

Thanks
------
- Joon Ro
- Michael Pennington
- Chris Uga
- Chris Withers
- Jeff Reback
- William Ferreira
- Daniel Fortunov
- Martin Felder

************************
pandas 0.3 Release Notes
************************

=============
Release Notes
=============

This major release of pandas represents approximately 1 year of continuous
development work and brings with it many new features, bug fixes, speed
enhancements, and general quality-of-life improvements. The most significant
change from the 0.2 release has been the completion of a rigorous unit test
suite covering all of the core functionality.

==========
What is it
==========

**pandas** is a library of labeled data structures, statistical models, and
general code for working with time series and cross-sectional data. It was
designed with the practical needs of statistical modeling and large,
inhomogeneous data sets in mind.

===============
Where to get it
===============

Source code: http://github.com/wesm/pandas
Binary installers on PyPI: http://pypi.python.org/pypi/pandas
Documentation: http://pandas.sourceforge.net

pandas 0.3.0 release notes
==========================

**Release date:** February 20, 2011

**New features / modules**

* DataFrame / DataMatrix classes
 * `corrwith` function to compute column- or row-wise correlations between two
   objects
 * Can boolean-index DataFrame objects, e.g. df[df > 2] = 2, px[px > last_px] = 0
 * Added comparison magic methods (__lt__, __gt__, etc.)
 * Flexible explicit arithmetic methods (add, mul, sub, div, etc.)
 * Added `reindex_like` method

* WidePanel
 * Added `reindex_like` method

* `pandas.io`: IO utilities
  * `pandas.io.sql` module
   * Convenience functions for accessing SQL-like databases
  * `pandas.io.pytables` module
   * Added (still experimental) HDFStore class for storing pandas data
     structures using HDF5 / PyTables

* `pandas.core.datetools`
 * Added WeekOfMonth date offset

* `pandas.rpy` (experimental) module created, provide some interfacing /
  conversion between rpy2 and pandas

**Improvements**

* Unit test coverage: 100% line coverage of core data structures
* Speed enhancement to rolling_{median, max, min}
* Column ordering between DataFrame and DataMatrix is now consistent: before
  DataFrame would not respect column order
* Improved {Series, DataFrame}.plot methods to be more flexible (can pass
  matplotlib Axis arguments, plot DataFrame columns in multiple subplots, etc.)

**API Changes**

* Exponentially-weighted moment functions in `pandas.stats.moments`
  have a more consistent API and accept a min_periods argument like
  their regular moving counterparts.
* **fillMethod** argument in Series, DataFrame changed to **method**,
    `FutureWarning` added.
* **fill** method in Series, DataFrame/DataMatrix, WidePanel renamed to
    **fillna**, `FutureWarning` added to **fill**
* Renamed **DataFrame.getXS** to **xs**, `FutureWarning` added
* Removed **cap** and **floor** functions from DataFrame, renamed to
  **clip_upper** and **clip_lower** for consistency with NumPy

**Bug fixes**

* Fixed bug in IndexableSkiplist Cython code that was breaking
  rolling_max function
* Numerous numpy.int64-related indexing fixes
* Several NumPy 1.4.0 NaN-handling fixes
* Bug fixes to pandas.io.parsers.parseCSV
* Fixed `DateRange` caching issue with unusual date offsets
* Fixed bug in `DateRange.union`
* Fixed corner case in `IndexableSkiplist` implementation
