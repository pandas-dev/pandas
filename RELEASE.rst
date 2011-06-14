
************************
pandas 0.4 Release Notes
************************

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

=============
Release notes
=============

**Release date:** NOT YET RELEASED

**New features / modules**

* `Series.describe`, `DataFrame.describe`: produces an R-like table of summary
  statistics about each data column
* `DataFrame.quantile`, `Series.quantile`
* Fancy indexing operator on Series / DataFrame, e.g.:
  * frame.ix[5:10, ['C', 'B', 'A']]
  * frame.ix[date1:date2]
* Boolean indexing with DataFrame objects: df[df > 1] = 1
* `pytz` / tzinfo support in `DateRange`
  * `tz_localize`, `tz_normalize`, and `tz_validate` methods added

**Improvements to existing features**

* The 2-dimensional `DataFrame` and `DataMatrix` classes have been extensively
  refactored internally into a single class `DataFrame`, preserving where
  possible their optimal performance characteristics. This should reduce
  confusion from users about which class to use
* Column ordering for mixed type data is now completely consistent in
  `DataFrame`. In prior releases, there was inconsistent column ordering in
  `DataMatrix`

**API Changes**

* The `DataMatrix` variable now refers to `DataFrame`, will be removed in next
  release
* Handling of boolean dtype in `DataFrame` has been improved to support storage
  of boolean data with NA / NaN values. Before it was being converted to float64
  so this should not (in theory) cause API breakage
* Backwards compatibility support for begin/end/nPeriods keyword arguments in
  DateRange class has been removed
* `pandas.core.mixins` code moved to `pandas.core.generic`
* `buffer` keyword arguments (e.g. `DataFrame.toString`) renamed to `buf` to
  avoid using Python built-in name
* `DataFrame.rows()` removed (use `DataFrame.index`)
* Added deprecation warning to `DataFrame.cols()`, to be removed in next release
* `DataFrame` deprecations: `merge`, `asMatrix`, `toDataMatrix`,
  `_firstTimeWithValue`, `_lastTimeWithValue`
* `fillMethod` arguments (deprecated in prior release) removed, should be
  replaced with `method`
* `Series.fill`, `DataFrame.fill`, and `WidePanel.fill` removed, use `fillna`
  instead

**Bug fixes**

* Column ordering in pandas.io.parsers.

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
