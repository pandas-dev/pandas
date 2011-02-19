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
