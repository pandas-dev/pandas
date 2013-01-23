=============
Release Notes
=============

This is the list of changes to pandas between each release. For full details,
see the commit logs at http://github.com/pydata/pandas

What is it
----------

pandas is a Python package providing fast, flexible, and expressive data
structures designed to make working with “relational” or “labeled” data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, real world data analysis in Python. Additionally, it has the
broader goal of becoming the most powerful and flexible open source data
analysis / manipulation tool available in any language.

Where to get it
---------------

* Source code: http://github.com/pydata/pandas
* Binary installers on PyPI: http://pypi.python.org/pypi/pandas
* Documentation: http://pandas.pydata.org

pandas 0.10.1
=============

**Release date:** 2013-01-22

**New features**

  - Add data inferface to World Bank WDI pandas.io.wb (#2592)

**API Changes**

  - Restored inplace=True behavior returning self (same object) with
    deprecation warning until 0.11 (GH1893_)
  - ``HDFStore``
    - refactored HFDStore to deal with non-table stores as objects, will allow future enhancements
    - removed keyword ``compression`` from ``put`` (replaced by keyword
      ``complib`` to be consistent across library)
    - warn `PerformanceWarning` if you are attempting to store types that will be pickled by PyTables

**Improvements to existing features**

  - ``HDFStore``

    - enables storing of multi-index dataframes (closes GH1277_)
    - support data column indexing and selection, via ``data_columns`` keyword in append
    - support write chunking to reduce memory footprint, via ``chunksize``
      keyword to append
    - support automagic indexing via ``index`` keywork to append
    - support ``expectedrows`` keyword in append to inform ``PyTables`` about
      the expected tablesize
    - support ``start`` and ``stop`` keywords in select to limit the row
      selection space
    - added ``get_store`` context manager to automatically import with pandas
    - added column filtering via ``columns`` keyword in select
    - added methods append_to_multiple/select_as_multiple/select_as_coordinates
      to do multiple-table append/selection
    - added support for datetime64 in columns
    - added method ``unique`` to select the unique values in an indexable or data column
    - added method ``copy`` to copy an existing store (and possibly upgrade)
    - show the shape of the data on disk for non-table stores when printing the store
    - added ability to read PyTables flavor tables (allows compatiblity to other HDF5 systems)
  - Add ``logx`` option to DataFrame/Series.plot (GH2327_, #2565)
  - Support reading gzipped data from file-like object
  - ``pivot_table`` aggfunc can be anything used in GroupBy.aggregate (GH2643_)
  - Implement DataFrame merges in case where set cardinalities might overflow
    64-bit integer (GH2690_)
  - Raise exception in C file parser if integer dtype specified and have NA
    values. (GH2631_)
  - Attempt to parse ISO8601 format dates when parse_dates=True in read_csv for
    major performance boost in such cases (GH2698_)
  - Add methods ``neg`` and ``inv`` to Series
  - Implement ``kind`` option in ``ExcelFile`` to indicate whether it's an XLS
    or XLSX file (GH2613_)

**Bug fixes**

  - Fix read_csv/read_table multithreading issues (GH2608_)
  - ``HDFStore``

    - correctly handle ``nan`` elements in string columns; serialize via the
      ``nan_rep`` keyword to append
    - raise correctly on non-implemented column types (unicode/date)
    - handle correctly ``Term`` passed types (e.g. ``index<1000``, when index
      is ``Int64``), (closes GH512_)
    - handle Timestamp correctly in data_columns (closes GH2637_)
    - contains correctly matches on non-natural names
    - correctly store ``float32`` dtypes in tables (if not other float types in
      the same table)
  - Fix DataFrame.info bug with UTF8-encoded columns. (GH2576_)
  - Fix DatetimeIndex handling of FixedOffset tz (GH2604_)
  - More robust detection of being in IPython session for wide DataFrame
    console formatting (GH2585_)
  - Fix platform issues with ``file:///`` in unit test (#2564)
  - Fix bug and possible segfault when grouping by hierarchical level that
    contains NA values (GH2616_)
  - Ensure that MultiIndex tuples can be constructed with NAs (seen in #2616)
  - Fix int64 overflow issue when unstacking MultiIndex with many levels (#2616)
  - Exclude non-numeric data from DataFrame.quantile by default (GH2625_)
  - Fix a Cython C int64 boxing issue causing read_csv to return incorrect
    results (GH2599_)
  - Fix groupby summing performance issue on boolean data (GH2692_)
  - Don't bork Series containing datetime64 values with to_datetime (GH2699_)
  - Fix DataFrame.from_records corner case when passed columns, index column,
    but empty record list (GH2633_)
  - Fix C parser-tokenizer bug with trailing fields. (GH2668_)
  - Don't exclude non-numeric data from GroupBy.max/min (GH2700_)
  - Don't lose time zone when calling DatetimeIndex.drop (GH2621_)
  - Fix setitem on a Series with a boolean key and a non-scalar as value (GH2686_)
  - Box datetime64 values in Series.apply/map (GH2627_, GH2689_)
  - Upconvert datetime + datetime64 values when concatenating frames (GH2624_)
  - Raise a more helpful error message in merge operations when one DataFrame
    has duplicate columns (GH2649_)
  - Fix partial date parsing issue occuring only when code is run at EOM  (GH2618_)
  - Prevent MemoryError when using counting sort in sortlevel with
    high-cardinality MultiIndex objects (GH2684_)
  - Fix Period resampling bug when all values fall into a single bin (GH2070_)
  - Fix buggy interaction with usecols argument in read_csv when there is an
    implicit first index column (GH2654_)

.. _GH512: https://github.com/pydata/pandas/issues/512
.. _GH1277: https://github.com/pydata/pandas/issues/1277
.. _GH2070: https://github.com/pydata/pandas/issues/2070
.. _GH2327: https://github.com/pydata/pandas/issues/2327
.. _GH2585: https://github.com/pydata/pandas/issues/2585
.. _GH2599: https://github.com/pydata/pandas/issues/2599
.. _GH2604: https://github.com/pydata/pandas/issues/2604
.. _GH2576: https://github.com/pydata/pandas/issues/2576
.. _GH2608: https://github.com/pydata/pandas/issues/2608
.. _GH2613: https://github.com/pydata/pandas/issues/2613
.. _GH2616: https://github.com/pydata/pandas/issues/2616
.. _GH2621: https://github.com/pydata/pandas/issues/2621
.. _GH2624: https://github.com/pydata/pandas/issues/2624
.. _GH2625: https://github.com/pydata/pandas/issues/2625
.. _GH2627: https://github.com/pydata/pandas/issues/2627
.. _GH2631: https://github.com/pydata/pandas/issues/2631
.. _GH2633: https://github.com/pydata/pandas/issues/2633
.. _GH2637: https://github.com/pydata/pandas/issues/2637
.. _GH2643: https://github.com/pydata/pandas/issues/2643
.. _GH2649: https://github.com/pydata/pandas/issues/2649
.. _GH2654: https://github.com/pydata/pandas/issues/2654
.. _GH2668: https://github.com/pydata/pandas/issues/2668
.. _GH2684: https://github.com/pydata/pandas/issues/2684
.. _GH2689: https://github.com/pydata/pandas/issues/2689
.. _GH2690: https://github.com/pydata/pandas/issues/2690
.. _GH2692: https://github.com/pydata/pandas/issues/2692
.. _GH2698: https://github.com/pydata/pandas/issues/2698
.. _GH2699: https://github.com/pydata/pandas/issues/2699
.. _GH2700: https://github.com/pydata/pandas/issues/2700
.. _GH2694: https://github.com/pydata/pandas/issues/2694
.. _GH2686: https://github.com/pydata/pandas/issues/2686
.. _GH2618: https://github.com/pydata/pandas/issues/2618

pandas 0.10.0
=============

**Release date:** 2012-12-17

**New features**

  - Brand new high-performance delimited file parsing engine written in C and
    Cython. 50% or better performance in many standard use cases with a
    fraction as much memory usage. (GH407_, GH821_)
  - Many new file parser (read_csv, read_table) features:

    - Support for on-the-fly gzip or bz2 decompression (`compression` option)
    - Ability to get back numpy.recarray instead of DataFrame
      (`as_recarray=True`)
    - `dtype` option: explicit column dtypes
    - `usecols` option: specify list of columns to be read from a file. Good
      for reading very wide files with many irrelevant columns (GH1216_ GH926_, GH2465_)
    - Enhanced unicode decoding support via `encoding` option
    - `skipinitialspace` dialect option
    - Can specify strings to be recognized as True (`true_values`) or False
      (`false_values`)
    - High-performance `delim_whitespace` option for whitespace-delimited
      files; a preferred alternative to the '\s+' regular expression delimiter
    - Option to skip "bad" lines (wrong number of fields) that would otherwise
      have caused an error in the past (`error_bad_lines` and `warn_bad_lines`
      options)
    - Substantially improved performance in the parsing of integers with
      thousands markers and lines with comments
    - Easy of European (and other) decimal formats (`decimal` option) (GH584_, GH2466_)
    - Custom line terminators (e.g. lineterminator='~') (GH2457_)
    - Handling of no trailing commas in CSV files (GH2333_)
    - Ability to handle fractional seconds in date_converters (GH2209_)
    - read_csv allow scalar arg to na_values (GH1944_)
    - Explicit column dtype specification in read_* functions (GH1858_)
    - Easier CSV dialect specification (GH1743_)
    - Improve parser performance when handling special characters (GH1204_)

  - Google Analytics API integration with easy oauth2 workflow (GH2283_)
  - Add error handling to Series.str.encode/decode (GH2276_)
  - Add ``where`` and ``mask`` to Series (GH2337_)
  - Grouped histogram via `by` keyword in Series/DataFrame.hist (GH2186_)
  - Support optional ``min_periods`` keyword in ``corr`` and ``cov``
    for both Series and DataFrame (GH2002_)
  - Add ``duplicated`` and ``drop_duplicates`` functions to Series (GH1923_)
  - Add docs for ``HDFStore table`` format
  - 'density' property in `SparseSeries` (GH2384_)
  - Add ``ffill`` and ``bfill`` convenience functions for forward- and
    backfilling time series data (GH2284_)
  - New option configuration system and functions `set_option`, `get_option`,
    `describe_option`, and `reset_option`. Deprecate `set_printoptions` and
    `reset_printoptions` (GH2393_).
    You can also access options as attributes via ``pandas.options.X``
  - Wide DataFrames can be viewed more easily in the console with new
    `expand_frame_repr` and `line_width` configuration options. This is on by
    default now (GH2436_)
  - Scikits.timeseries-like moving window functions via ``rolling_window`` (GH1270_)

**Experimental Features**

  - Add support for Panel4D, a named 4 Dimensional stucture
  - Add support for ndpanel factory functions, to create custom,
    domain-specific N-Dimensional containers

**API Changes**

  - The default binning/labeling behavior for ``resample`` has been changed to
    `closed='left', label='left'` for daily and lower frequencies. This had
    been a large source of confusion for users. See "what's new" page for more
    on this. (GH2410_)
  - Methods with ``inplace`` option now return None instead of the calling
    (modified) object (GH1893_)
  - The special case DataFrame - TimeSeries doing column-by-column broadcasting
    has been deprecated. Users should explicitly do e.g. df.sub(ts, axis=0)
    instead. This is a legacy hack and can lead to subtle bugs.
  - inf/-inf are no longer considered as NA by isnull/notnull. To be clear, this
    is legacy cruft from early pandas. This behavior can be globally re-enabled
    using the new option ``mode.use_inf_as_null`` (GH2050_, GH1919_)
  - ``pandas.merge`` will now default to ``sort=False``. For many use cases
    sorting the join keys is not necessary, and doing it by default is wasteful
  - Specify ``header=0`` explicitly to replace existing column names in file in
    read_* functions.
  - Default column names for header-less parsed files (yielded by read_csv,
    etc.) are now the integers 0, 1, .... A new argument `prefix` has been
    added; to get the v0.9.x behavior specify ``prefix='X'`` (GH2034_). This API
    change was made to make the default column names more consistent with the
    DataFrame constructor's default column names when none are specified.
  - DataFrame selection using a boolean frame now preserves input shape
  - If function passed to Series.apply yields a Series, result will be a
    DataFrame (GH2316_)
  - Values like YES/NO/yes/no will not be considered as boolean by default any
    longer in the file parsers. This can be customized using the new
    ``true_values`` and ``false_values`` options (GH2360_)
  - `obj.fillna()` is no longer valid; make `method='pad'` no longer the
    default option, to be more explicit about what kind of filling to
    perform. Add `ffill/bfill` convenience functions per above (GH2284_)
  - `HDFStore.keys()` now returns an absolute path-name for each key
  - `to_string()` now always returns a unicode string. (GH2224_)
  - File parsers will not handle NA sentinel values arising from passed
    converter functions

**Improvements to existing features**

  - Add ``nrows`` option to DataFrame.from_records for iterators (GH1794_)
  - Unstack/reshape algorithm rewrite to avoid high memory use in cases where
    the number of observed key-tuples is much smaller than the total possible
    number that could occur (GH2278_). Also improves performance in most cases.
  - Support duplicate columns in DataFrame.from_records (GH2179_)
  - Add ``normalize`` option to Series/DataFrame.asfreq (GH2137_)
  - SparseSeries and SparseDataFrame construction from empty and scalar
    values now no longer create dense ndarrays unnecessarily (GH2322_)
  - ``HDFStore`` now supports hierarchial keys (GH2397_)
  - Support multiple query selection formats for ``HDFStore tables`` (GH1996_)
  - Support ``del store['df']`` syntax to delete HDFStores
  - Add multi-dtype support for ``HDFStore tables``
  - ``min_itemsize`` parameter can be specified in ``HDFStore table`` creation
  - Indexing support in ``HDFStore tables`` (GH698_)
  - Add `line_terminator` option to DataFrame.to_csv (GH2383_)
  - added implementation of str(x)/unicode(x)/bytes(x) to major pandas data
    structures, which should do the right thing on both py2.x and py3.x. (GH2224_)
  - Reduce groupby.apply overhead substantially by low-level manipulation of
    internal NumPy arrays in DataFrames (GH535_)
  - Implement ``value_vars`` in ``melt`` and add ``melt`` to pandas namespace
    (GH2412_)
  - Added boolean comparison operators to Panel
  - Enable ``Series.str.strip/lstrip/rstrip`` methods to take an argument (GH2411_)
  - The DataFrame ctor now respects column ordering when given
    an OrderedDict (GH2455_)
  - Assigning DatetimeIndex to Series changes the class to TimeSeries (GH2139_)
  - Improve performance of .value_counts method on non-integer data (GH2480_)
  - ``get_level_values`` method for MultiIndex return Index instead of ndarray (GH2449_)
  - ``convert_to_r_dataframe`` conversion for datetime values (GH2351_)
  - Allow ``DataFrame.to_csv`` to represent inf and nan differently (GH2026_)
  - Add ``min_i`` argument to ``nancorr`` to specify minimum required observations (GH2002_)
  - Add ``inplace`` option to ``sortlevel`` / ``sort`` functions on DataFrame (GH1873_)
  - Enable DataFrame to accept scalar constructor values like Series (GH1856_)
  - DataFrame.from_records now takes optional ``size`` parameter (GH1794_)
  - include iris dataset (GH1709_)
  - No datetime64 DataFrame column conversion of datetime.datetime with tzinfo (GH1581_)
  - Micro-optimizations in DataFrame for tracking state of internal consolidation (GH217_)
  - Format parameter in DataFrame.to_csv (GH1525_)
  - Partial string slicing for ``DatetimeIndex`` for daily and higher frequencies (GH2306_)
  - Implement ``col_space`` parameter in ``to_html`` and ``to_string`` in DataFrame (GH1000_)
  - Override ``Series.tolist`` and box datetime64 types (GH2447_)
  - Optimize ``unstack`` memory usage by compressing indices (GH2278_)
  - Fix HTML repr in IPython qtconsole if opening window is small (GH2275_)
  - Escape more special characters in console output (GH2492_)
  - df.select now invokes bool on the result of crit(x) (GH2487_)

**Bug fixes**

  - Fix major performance regression in DataFrame.iteritems (GH2273_)
  - Fixes bug when negative period passed to Series/DataFrame.diff (GH2266_)
  - Escape tabs in console output to avoid alignment issues (GH2038_)
  - Properly box datetime64 values when retrieving cross-section from
    mixed-dtype DataFrame (GH2272_)
  - Fix concatenation bug leading to GH2057_, GH2257_
  - Fix regression in Index console formatting (GH2319_)
  - Box Period data when assigning PeriodIndex to frame column (GH2243_, GH2281_)
  - Raise exception on calling reset_index on Series with inplace=True (GH2277_)
  - Enable setting multiple columns in DataFrame with hierarchical columns
    (GH2295_)
  - Respect dtype=object in DataFrame constructor (GH2291_)
  - Fix DatetimeIndex.join bug with tz-aware indexes and how='outer' (GH2317_)
  - pop(...) and del works with DataFrame with duplicate columns (GH2349_)
  - Treat empty strings as NA in date parsing (rather than let dateutil do
    something weird) (GH2263_)
  - Prevent uint64 -> int64 overflows (GH2355_)
  - Enable joins between MultiIndex and regular Index (GH2024_)
  - Fix time zone metadata issue when unioning non-overlapping DatetimeIndex
    objects (GH2367_)
  - Raise/handle int64 overflows in parsers (GH2247_)
  - Deleting of consecutive rows in ``HDFStore tables``` is much faster than before
  - Appending on a HDFStore would fail if the table was not first created via ``put``
  - Use `col_space` argument as minimum column width in DataFrame.to_html (GH2328_)
  - Fix tz-aware DatetimeIndex.to_period (GH2232_)
  - Fix DataFrame row indexing case with MultiIndex (GH2314_)
  - Fix to_excel exporting issues with Timestamp objects in index (GH2294_)
  - Fixes assigning scalars and array to hierarchical column chunk (GH1803_)
  - Fixed a UnicdeDecodeError with series tidy_repr (GH2225_)
  - Fixed issued with duplicate keys in an index (GH2347_, GH2380_)
  - Fixed issues re: Hash randomization, default on starting w/ py3.3 (GH2331_)
  - Fixed issue with missing attributes after loading a pickled dataframe (GH2431_)
  - Fix Timestamp formatting with tzoffset time zone in dateutil 2.1 (GH2443_)
  - Fix GroupBy.apply issue when using BinGrouper to do ts binning (GH2300_)
  - Fix issues resulting from datetime.datetime columns being converted to
    datetime64 when calling DataFrame.apply. (GH2374_)
  - Raise exception when calling to_panel on non uniquely-indexed frame (GH2441_)
  - Improved detection of console encoding on IPython zmq frontends (GH2458_)
  - Preserve time zone when .append-ing two time series (GH2260_)
  - Box timestamps when calling reset_index on time-zone-aware index rather
    than creating a tz-less datetime64 column (GH2262_)
  - Enable searching non-string columns in DataFrame.filter(like=...) (GH2467_)
  - Fixed issue with losing nanosecond precision upon conversion to DatetimeIndex(GH2252_)
  - Handle timezones in Datetime.normalize (GH2338_)
  - Fix test case where dtype specification with endianness causes
    failures on big endian machines (GH2318_)
  - Fix plotting bug where upsampling causes data to appear shifted in time (GH2448_)
  - Fix ``read_csv`` failure for UTF-16 with BOM and skiprows(GH2298_)
  - read_csv with names arg not implicitly setting header=None(GH2459_)
  - Unrecognized compression mode causes segfault in read_csv(GH2474_)
  - In read_csv, header=0 and passed names should discard first row(GH2269_)
  - Correctly route to stdout/stderr in read_table (GH2071_)
  - Fix exception when Timestamp.to_datetime is called on a Timestamp with tzoffset (GH2471_)
  - Fixed unintentional conversion of datetime64 to long in groupby.first() (GH2133_)
  - Union of empty DataFrames now return empty with concatenated index (GH2307_)
  - DataFrame.sort_index raises more helpful exception if sorting by column
    with duplicates (GH2488_)
  - DataFrame.to_string formatters can be list, too (GH2520_)
  - DataFrame.combine_first will always result in the union of the index and
    columns, even if one DataFrame is length-zero (GH2525_)
  - Fix several DataFrame.icol/irow with duplicate indices issues (GH2228_, GH2259_)
  - Use Series names for column names when using concat with axis=1 (GH2489_)
  - Raise Exception if start, end, periods all passed to date_range (GH2538_)
  - Fix Panel resampling issue (GH2537_)

.. _GH407: https://github.com/pydata/pandas/issues/407
.. _GH821: https://github.com/pydata/pandas/issues/821
.. _GH1216: https://github.com/pydata/pandas/issues/1216
.. _GH926: https://github.com/pydata/pandas/issues/926
.. _GH2465: https://github.com/pydata/pandas/issues/2465
.. _GH584: https://github.com/pydata/pandas/issues/584
.. _GH2466: https://github.com/pydata/pandas/issues/2466
.. _GH2457: https://github.com/pydata/pandas/issues/2457
.. _GH2333: https://github.com/pydata/pandas/issues/2333
.. _GH2209: https://github.com/pydata/pandas/issues/2209
.. _GH1944: https://github.com/pydata/pandas/issues/1944
.. _GH1858: https://github.com/pydata/pandas/issues/1858
.. _GH1743: https://github.com/pydata/pandas/issues/1743
.. _GH1204: https://github.com/pydata/pandas/issues/1204
.. _GH2283: https://github.com/pydata/pandas/issues/2283
.. _GH2276: https://github.com/pydata/pandas/issues/2276
.. _GH2337: https://github.com/pydata/pandas/issues/2337
.. _GH2186: https://github.com/pydata/pandas/issues/2186
.. _GH2002: https://github.com/pydata/pandas/issues/2002
.. _GH1923: https://github.com/pydata/pandas/issues/1923
.. _GH2384: https://github.com/pydata/pandas/issues/2384
.. _GH2284: https://github.com/pydata/pandas/issues/2284
.. _GH2393: https://github.com/pydata/pandas/issues/2393
.. _GH2436: https://github.com/pydata/pandas/issues/2436
.. _GH1270: https://github.com/pydata/pandas/issues/1270
.. _GH2410: https://github.com/pydata/pandas/issues/2410
.. _GH1893: https://github.com/pydata/pandas/issues/1893
.. _GH2050: https://github.com/pydata/pandas/issues/2050
.. _GH1919: https://github.com/pydata/pandas/issues/1919
.. _GH2034: https://github.com/pydata/pandas/issues/2034
.. _GH2316: https://github.com/pydata/pandas/issues/2316
.. _GH2360: https://github.com/pydata/pandas/issues/2360
.. _GH2224: https://github.com/pydata/pandas/issues/2224
.. _GH1794: https://github.com/pydata/pandas/issues/1794
.. _GH2278: https://github.com/pydata/pandas/issues/2278
.. _GH2179: https://github.com/pydata/pandas/issues/2179
.. _GH2137: https://github.com/pydata/pandas/issues/2137
.. _GH2322: https://github.com/pydata/pandas/issues/2322
.. _GH2397: https://github.com/pydata/pandas/issues/2397
.. _GH1996: https://github.com/pydata/pandas/issues/1996
.. _GH698: https://github.com/pydata/pandas/issues/698
.. _GH2383: https://github.com/pydata/pandas/issues/2383
.. _GH535: https://github.com/pydata/pandas/issues/535
.. _GH2412: https://github.com/pydata/pandas/issues/2412
.. _GH2411: https://github.com/pydata/pandas/issues/2411
.. _GH2455: https://github.com/pydata/pandas/issues/2455
.. _GH2139: https://github.com/pydata/pandas/issues/2139
.. _GH2480: https://github.com/pydata/pandas/issues/2480
.. _GH2449: https://github.com/pydata/pandas/issues/2449
.. _GH2351: https://github.com/pydata/pandas/issues/2351
.. _GH2026: https://github.com/pydata/pandas/issues/2026
.. _GH1873: https://github.com/pydata/pandas/issues/1873
.. _GH1856: https://github.com/pydata/pandas/issues/1856
.. _GH1709: https://github.com/pydata/pandas/issues/1709
.. _GH1581: https://github.com/pydata/pandas/issues/1581
.. _GH217: https://github.com/pydata/pandas/issues/217
.. _GH1525: https://github.com/pydata/pandas/issues/1525
.. _GH2306: https://github.com/pydata/pandas/issues/2306
.. _GH1000: https://github.com/pydata/pandas/issues/1000
.. _GH2447: https://github.com/pydata/pandas/issues/2447
.. _GH2275: https://github.com/pydata/pandas/issues/2275
.. _GH2492: https://github.com/pydata/pandas/issues/2492
.. _GH2487: https://github.com/pydata/pandas/issues/2487
.. _GH2273: https://github.com/pydata/pandas/issues/2273
.. _GH2266: https://github.com/pydata/pandas/issues/2266
.. _GH2038: https://github.com/pydata/pandas/issues/2038
.. _GH2272: https://github.com/pydata/pandas/issues/2272
.. _GH2057: https://github.com/pydata/pandas/issues/2057
.. _GH2257: https://github.com/pydata/pandas/issues/2257
.. _GH2319: https://github.com/pydata/pandas/issues/2319
.. _GH2243: https://github.com/pydata/pandas/issues/2243
.. _GH2281: https://github.com/pydata/pandas/issues/2281
.. _GH2277: https://github.com/pydata/pandas/issues/2277
.. _GH2295: https://github.com/pydata/pandas/issues/2295
.. _GH2291: https://github.com/pydata/pandas/issues/2291
.. _GH2317: https://github.com/pydata/pandas/issues/2317
.. _GH2349: https://github.com/pydata/pandas/issues/2349
.. _GH2263: https://github.com/pydata/pandas/issues/2263
.. _GH2355: https://github.com/pydata/pandas/issues/2355
.. _GH2024: https://github.com/pydata/pandas/issues/2024
.. _GH2367: https://github.com/pydata/pandas/issues/2367
.. _GH2247: https://github.com/pydata/pandas/issues/2247
.. _GH2328: https://github.com/pydata/pandas/issues/2328
.. _GH2232: https://github.com/pydata/pandas/issues/2232
.. _GH2314: https://github.com/pydata/pandas/issues/2314
.. _GH2294: https://github.com/pydata/pandas/issues/2294
.. _GH1803: https://github.com/pydata/pandas/issues/1803
.. _GH2225: https://github.com/pydata/pandas/issues/2225
.. _GH2347: https://github.com/pydata/pandas/issues/2347
.. _GH2380: https://github.com/pydata/pandas/issues/2380
.. _GH2331: https://github.com/pydata/pandas/issues/2331
.. _GH2431: https://github.com/pydata/pandas/issues/2431
.. _GH2443: https://github.com/pydata/pandas/issues/2443
.. _GH2300: https://github.com/pydata/pandas/issues/2300
.. _GH2374: https://github.com/pydata/pandas/issues/2374
.. _GH2441: https://github.com/pydata/pandas/issues/2441
.. _GH2458: https://github.com/pydata/pandas/issues/2458
.. _GH2260: https://github.com/pydata/pandas/issues/2260
.. _GH2262: https://github.com/pydata/pandas/issues/2262
.. _GH2467: https://github.com/pydata/pandas/issues/2467
.. _GH2252: https://github.com/pydata/pandas/issues/2252
.. _GH2338: https://github.com/pydata/pandas/issues/2338
.. _GH2318: https://github.com/pydata/pandas/issues/2318
.. _GH2448: https://github.com/pydata/pandas/issues/2448
.. _GH2298: https://github.com/pydata/pandas/issues/2298
.. _GH2459: https://github.com/pydata/pandas/issues/2459
.. _GH2474: https://github.com/pydata/pandas/issues/2474
.. _GH2269: https://github.com/pydata/pandas/issues/2269
.. _GH2071: https://github.com/pydata/pandas/issues/2071
.. _GH2471: https://github.com/pydata/pandas/issues/2471
.. _GH2133: https://github.com/pydata/pandas/issues/2133
.. _GH2307: https://github.com/pydata/pandas/issues/2307
.. _GH2488: https://github.com/pydata/pandas/issues/2488
.. _GH2520: https://github.com/pydata/pandas/issues/2520
.. _GH2525: https://github.com/pydata/pandas/issues/2525
.. _GH2228: https://github.com/pydata/pandas/issues/2228
.. _GH2259: https://github.com/pydata/pandas/issues/2259
.. _GH2489: https://github.com/pydata/pandas/issues/2489
.. _GH2538: https://github.com/pydata/pandas/issues/2538
.. _GH2537: https://github.com/pydata/pandas/issues/2537


pandas 0.9.1
============

**Release date:** 2012-11-14

**New features**

  - Can specify multiple sort orders in DataFrame/Series.sort/sort_index (GH928_)
  - New `top` and `bottom` options for handling NAs in rank (GH1508_, GH2159_)
  - Add `where` and `mask` functions to DataFrame (GH2109_, GH2151_)
  - Add `at_time` and `between_time` functions to DataFrame (GH2149_)
  - Add flexible `pow` and `rpow` methods to DataFrame (GH2190_)

**API Changes**

  - Upsampling period index "spans" intervals. Example: annual periods
    upsampled to monthly will span all months in each year
  - Period.end_time will yield timestamp at last nanosecond in the interval
    (GH2124_, GH2125_, GH1764_)
  - File parsers no longer coerce to float or bool for columns that have custom
    converters specified (GH2184_)

**Improvements to existing features**

  - Time rule inference for week-of-month (e.g. WOM-2FRI) rules (GH2140_)
  - Improve performance of datetime + business day offset with large number of
    offset periods
  - Improve HTML display of DataFrame objects with hierarchical columns
  - Enable referencing of Excel columns by their column names (GH1936_)
  - DataFrame.dot can accept ndarrays (GH2042_)
  - Support negative periods in Panel.shift (GH2164_)
  - Make .drop(...) work with non-unique indexes (GH2101_)
  - Improve performance of Series/DataFrame.diff (re: GH2087_)
  - Support unary ~ (__invert__) in DataFrame (GH2110_)
  - Turn off pandas-style tick locators and formatters (GH2205_)
  - DataFrame[DataFrame] uses DataFrame.where to compute masked frame (GH2230_)

**Bug fixes**

  - Fix some duplicate-column DataFrame constructor issues (GH2079_)
  - Fix bar plot color cycle issues (GH2082_)
  - Fix off-center grid for stacked bar plots (GH2157_)
  - Fix plotting bug if inferred frequency is offset with N > 1 (GH2126_)
  - Implement comparisons on date offsets with fixed delta (GH2078_)
  - Handle inf/-inf correctly in read_* parser functions (GH2041_)
  - Fix matplotlib unicode interaction bug
  - Make WLS r-squared match statsmodels 0.5.0 fixed value
  - Fix zero-trimming DataFrame formatting bug
  - Correctly compute/box datetime64 min/max values from Series.min/max (GH2083_)
  - Fix unstacking edge case with unrepresented groups (GH2100_)
  - Fix Series.str failures when using pipe pattern '|' (GH2119_)
  - Fix pretty-printing of dict entries in Series, DataFrame (GH2144_)
  - Cast other datetime64 values to nanoseconds in DataFrame ctor (GH2095_)
  - Alias Timestamp.astimezone to tz_convert, so will yield Timestamp (GH2060_)
  - Fix timedelta64 formatting from Series (GH2165_, GH2146_)
  - Handle None values gracefully in dict passed to Panel constructor (GH2075_)
  - Box datetime64 values as Timestamp objects in Series/DataFrame.iget (GH2148_)
  - Fix Timestamp indexing bug in DatetimeIndex.insert (GH2155_)
  - Use index name(s) (if any) in DataFrame.to_records (GH2161_)
  - Don't lose index names in Panel.to_frame/DataFrame.to_panel (GH2163_)
  - Work around length-0 boolean indexing NumPy bug (GH2096_)
  - Fix partial integer indexing bug in DataFrame.xs (GH2107_)
  - Fix variety of cut/qcut string-bin formatting bugs (GH1978_, GH1979_)
  - Raise Exception when xs view not possible of MultiIndex'd DataFrame (GH2117_)
  - Fix groupby(...).first() issue with datetime64 (GH2133_)
  - Better floating point error robustness in some rolling_* functions
    (GH2114_, GH2527_)
  - Fix ewma NA handling in the middle of Series (GH2128_)
  - Fix numerical precision issues in diff with integer data (GH2087_)
  - Fix bug in MultiIndex.__getitem__ with NA values (GH2008_)
  - Fix DataFrame.from_records dict-arg bug when passing columns (GH2179_)
  - Fix Series and DataFrame.diff for integer dtypes (GH2087_, GH2174_)
  - Fix bug when taking intersection of DatetimeIndex with empty index (GH2129_)
  - Pass through timezone information when calling DataFrame.align (GH2127_)
  - Properly sort when joining on datetime64 values (GH2196_)
  - Fix indexing bug in which False/True were being coerced to 0/1 (GH2199_)
  - Many unicode formatting fixes (GH2201_)
  - Fix improper MultiIndex conversion issue when assigning
    e.g. DataFrame.index (GH2200_)
  - Fix conversion of mixed-type DataFrame to ndarray with dup columns (GH2236_)
  - Fix duplicate columns issue (GH2218_, GH2219_)
  - Fix SparseSeries.__pow__ issue with NA input (GH2220_)
  - Fix icol with integer sequence failure (GH2228_)
  - Fixed resampling tz-aware time series issue (GH2245_)
  - SparseDataFrame.icol was not returning SparseSeries (GH2227_, GH2229_)
  - Enable ExcelWriter to handle PeriodIndex (GH2240_)
  - Fix issue constructing DataFrame from empty Series with name (GH2234_)
  - Use console-width detection in interactive sessions only (GH1610_)
  - Fix parallel_coordinates legend bug with mpl 1.2.0 (GH2237_)
  - Make tz_localize work in corner case of empty Series (GH2248_)

.. _GH928: https://github.com/pydata/pandas/issues/928
.. _GH1508: https://github.com/pydata/pandas/issues/1508
.. _GH2159: https://github.com/pydata/pandas/issues/2159
.. _GH2109: https://github.com/pydata/pandas/issues/2109
.. _GH2151: https://github.com/pydata/pandas/issues/2151
.. _GH2149: https://github.com/pydata/pandas/issues/2149
.. _GH2190: https://github.com/pydata/pandas/issues/2190
.. _GH2124: https://github.com/pydata/pandas/issues/2124
.. _GH2125: https://github.com/pydata/pandas/issues/2125
.. _GH1764: https://github.com/pydata/pandas/issues/1764
.. _GH2184: https://github.com/pydata/pandas/issues/2184
.. _GH2140: https://github.com/pydata/pandas/issues/2140
.. _GH1936: https://github.com/pydata/pandas/issues/1936
.. _GH2042: https://github.com/pydata/pandas/issues/2042
.. _GH2164: https://github.com/pydata/pandas/issues/2164
.. _GH2101: https://github.com/pydata/pandas/issues/2101
.. _GH2087: https://github.com/pydata/pandas/issues/2087
.. _GH2110: https://github.com/pydata/pandas/issues/2110
.. _GH2205: https://github.com/pydata/pandas/issues/2205
.. _GH2230: https://github.com/pydata/pandas/issues/2230
.. _GH2079: https://github.com/pydata/pandas/issues/2079
.. _GH2082: https://github.com/pydata/pandas/issues/2082
.. _GH2157: https://github.com/pydata/pandas/issues/2157
.. _GH2126: https://github.com/pydata/pandas/issues/2126
.. _GH2078: https://github.com/pydata/pandas/issues/2078
.. _GH2041: https://github.com/pydata/pandas/issues/2041
.. _GH2083: https://github.com/pydata/pandas/issues/2083
.. _GH2100: https://github.com/pydata/pandas/issues/2100
.. _GH2119: https://github.com/pydata/pandas/issues/2119
.. _GH2144: https://github.com/pydata/pandas/issues/2144
.. _GH2095: https://github.com/pydata/pandas/issues/2095
.. _GH2060: https://github.com/pydata/pandas/issues/2060
.. _GH2165: https://github.com/pydata/pandas/issues/2165
.. _GH2146: https://github.com/pydata/pandas/issues/2146
.. _GH2075: https://github.com/pydata/pandas/issues/2075
.. _GH2148: https://github.com/pydata/pandas/issues/2148
.. _GH2155: https://github.com/pydata/pandas/issues/2155
.. _GH2161: https://github.com/pydata/pandas/issues/2161
.. _GH2163: https://github.com/pydata/pandas/issues/2163
.. _GH2096: https://github.com/pydata/pandas/issues/2096
.. _GH2107: https://github.com/pydata/pandas/issues/2107
.. _GH1978: https://github.com/pydata/pandas/issues/1978
.. _GH1979: https://github.com/pydata/pandas/issues/1979
.. _GH2117: https://github.com/pydata/pandas/issues/2117
.. _GH2133: https://github.com/pydata/pandas/issues/2133
.. _GH2114: https://github.com/pydata/pandas/issues/2114
.. _GH2527: https://github.com/pydata/pandas/issues/2114
.. _GH2128: https://github.com/pydata/pandas/issues/2128
.. _GH2008: https://github.com/pydata/pandas/issues/2008
.. _GH2179: https://github.com/pydata/pandas/issues/2179
.. _GH2174: https://github.com/pydata/pandas/issues/2174
.. _GH2129: https://github.com/pydata/pandas/issues/2129
.. _GH2127: https://github.com/pydata/pandas/issues/2127
.. _GH2196: https://github.com/pydata/pandas/issues/2196
.. _GH2199: https://github.com/pydata/pandas/issues/2199
.. _GH2201: https://github.com/pydata/pandas/issues/2201
.. _GH2200: https://github.com/pydata/pandas/issues/2200
.. _GH2236: https://github.com/pydata/pandas/issues/2236
.. _GH2218: https://github.com/pydata/pandas/issues/2218
.. _GH2219: https://github.com/pydata/pandas/issues/2219
.. _GH2220: https://github.com/pydata/pandas/issues/2220
.. _GH2228: https://github.com/pydata/pandas/issues/2228
.. _GH2245: https://github.com/pydata/pandas/issues/2245
.. _GH2227: https://github.com/pydata/pandas/issues/2227
.. _GH2229: https://github.com/pydata/pandas/issues/2229
.. _GH2240: https://github.com/pydata/pandas/issues/2240
.. _GH2234: https://github.com/pydata/pandas/issues/2234
.. _GH1610: https://github.com/pydata/pandas/issues/1610
.. _GH2237: https://github.com/pydata/pandas/issues/2237
.. _GH2248: https://github.com/pydata/pandas/issues/2248


pandas 0.9.0
============

**Release date:** 10/7/2012

**New features**

  - Add ``str.encode`` and ``str.decode`` to Series (GH1706_)
  - Add `to_latex` method to DataFrame (GH1735_)
  - Add convenient expanding window equivalents of all rolling_* ops (GH1785_)
  - Add Options class to pandas.io.data for fetching options data from Yahoo!
    Finance (GH1748_, GH1739_)
  - Recognize and convert more boolean values in file parsing (Yes, No, TRUE,
    FALSE, variants thereof) (GH1691_, GH1295_)
  - Add Panel.update method, analogous to DataFrame.update (GH1999_, GH1988_)

**Improvements to existing features**

  - Proper handling of NA values in merge operations (GH1990_)
  - Add ``flags`` option for ``re.compile`` in some Series.str methods (GH1659_)
  - Parsing of UTC date strings in read_* functions (GH1693_)
  - Handle generator input to Series (GH1679_)
  - Add `na_action='ignore'` to Series.map to quietly propagate NAs (GH1661_)
  - Add args/kwds options to Series.apply (GH1829_)
  - Add inplace option to Series/DataFrame.reset_index (GH1797_)
  - Add ``level`` parameter to ``Series.reset_index``
  - Add quoting option for DataFrame.to_csv (GH1902_)
  - Indicate long column value truncation in DataFrame output with ... (GH1854_)
  - DataFrame.dot will not do data alignment, and also work with Series (GH1915_)
  - Add ``na`` option for missing data handling in some vectorized string
    methods (GH1689_)
  - If index_label=False in DataFrame.to_csv, do not print fields/commas in the
    text output. Results in easier importing into R (GH1583_)
  - Can pass tuple/list of axes to DataFrame.dropna to simplify repeated calls
    (dropping both columns and rows) (GH924_)
  - Improve DataFrame.to_html output for hierarchically-indexed rows (do not
    repeat levels) (GH1929_)
  - TimeSeries.between_time can now select times across midnight (GH1871_)
  - Enable `skip_footer` parameter in `ExcelFile.parse` (GH1843_)

**API Changes**

  - Change default header names in read_* functions to more Pythonic X0, X1,
    etc. instead of X.1, X.2. (GH2000_)
  - Deprecated ``day_of_year`` API removed from PeriodIndex, use ``dayofyear``
    (GH1723_)
  - Don't modify NumPy suppress printoption at import time
  - The internal HDF5 data arrangement for DataFrames has been
    transposed. Legacy files will still be readable by HDFStore (GH1834_, GH1824_)
  - Legacy cruft removed: pandas.stats.misc.quantileTS
  - Use ISO8601 format for Period repr: monthly, daily, and on down (GH1776_)
  - Empty DataFrame columns are now created as object dtype. This will prevent
    a class of TypeErrors that was occurring in code where the dtype of a
    column would depend on the presence of data or not (e.g. a SQL query having
    results) (GH1783_)
  - Setting parts of DataFrame/Panel using ix now aligns input Series/DataFrame
    (GH1630_)
  - `first` and `last` methods in `GroupBy` no longer drop non-numeric columns
    (GH1809_)
  - Resolved inconsistencies in specifying custom NA values in text parser.
    `na_values` of type dict no longer override default NAs unless
    `keep_default_na` is set to false explicitly (GH1657_)
  - Enable `skipfooter` parameter in text parsers as an alias for `skip_footer`

**Bug fixes**

  - Perform arithmetic column-by-column in mixed-type DataFrame to avoid type
    upcasting issues. Caused downstream DataFrame.diff bug (GH1896_)
  - Fix matplotlib auto-color assignment when no custom spectrum passed. Also
    respect passed color keyword argument (GH1711_)
  - Fix resampling logical error with closed='left' (GH1726_)
  - Fix critical DatetimeIndex.union bugs (GH1730_, GH1719_, GH1745_, GH1702_, GH1753_)
  - Fix critical DatetimeIndex.intersection bug with unanchored offsets (GH1708_)
  - Fix MM-YYYY time series indexing case (GH1672_)
  - Fix case where Categorical group key was not being passed into index in
    GroupBy result (GH1701_)
  - Handle Ellipsis in Series.__getitem__/__setitem__ (GH1721_)
  - Fix some bugs with handling datetime64 scalars of other units in NumPy 1.6
    and 1.7 (GH1717_)
  - Fix performance issue in MultiIndex.format (GH1746_)
  - Fixed GroupBy bugs interacting with DatetimeIndex asof / map methods (GH1677_)
  - Handle factors with NAs in pandas.rpy (GH1615_)
  - Fix statsmodels import in pandas.stats.var (GH1734_)
  - Fix DataFrame repr/info summary with non-unique columns (GH1700_)
  - Fix Series.iget_value for non-unique indexes (GH1694_)
  - Don't lose tzinfo when passing DatetimeIndex as DataFrame column (GH1682_)
  - Fix tz conversion with time zones that haven't had any DST transitions since
    first date in the array (GH1673_)
  - Fix field access with  UTC->local conversion on unsorted arrays (GH1756_)
  - Fix isnull handling of array-like (list) inputs (GH1755_)
  - Fix regression in handling of Series in Series constructor (GH1671_)
  - Fix comparison of Int64Index with DatetimeIndex (GH1681_)
  - Fix min_periods handling in new rolling_max/min at array start (GH1695_)
  - Fix errors with how='median' and generic NumPy resampling in some cases
    caused by SeriesBinGrouper (GH1648_, GH1688_)
  - When grouping by level, exclude unobserved levels (GH1697_)
  - Don't lose tzinfo in DatetimeIndex when shifting by different offset (GH1683_)
  - Hack to support storing data with a zero-length axis in HDFStore (GH1707_)
  - Fix DatetimeIndex tz-aware range generation issue (GH1674_)
  - Fix method='time' interpolation with intraday data (GH1698_)
  - Don't plot all-NA DataFrame columns as zeros (GH1696_)
  - Fix bug in scatter_plot with by option (GH1716_)
  - Fix performance problem in infer_freq with lots of non-unique stamps (GH1686_)
  - Fix handling of PeriodIndex as argument to create MultiIndex (GH1705_)
  - Fix re: unicode MultiIndex level names in Series/DataFrame repr (GH1736_)
  - Handle PeriodIndex in to_datetime instance method (GH1703_)
  - Support StaticTzInfo in DatetimeIndex infrastructure (GH1692_)
  - Allow MultiIndex setops with length-0 other type indexes (GH1727_)
  - Fix handling of DatetimeIndex in DataFrame.to_records (GH1720_)
  - Fix handling of general objects in isnull on which bool(...) fails (GH1749_)
  - Fix .ix indexing with MultiIndex ambiguity (GH1678_)
  - Fix .ix setting logic error with non-unique MultiIndex (GH1750_)
  - Basic indexing now works on MultiIndex with > 1000000 elements, regression
    from earlier version of pandas (GH1757_)
  - Handle non-float64 dtypes in fast DataFrame.corr/cov code paths (GH1761_)
  - Fix DatetimeIndex.isin to function properly (GH1763_)
  - Fix conversion of array of tz-aware datetime.datetime to DatetimeIndex with
    right time zone (GH1777_)
  - Fix DST issues with generating ancxhored date ranges (GH1778_)
  - Fix issue calling sort on result of Series.unique (GH1807_)
  - Fix numerical issue leading to square root of negative number in
    rolling_std (GH1840_)
  - Let Series.str.split accept no arguments (like str.split) (GH1859_)
  - Allow user to have dateutil 2.1 installed on a Python 2 system (GH1851_)
  - Catch ImportError less aggressively in pandas/__init__.py (GH1845_)
  - Fix pip source installation bug when installing from GitHub (GH1805_)
  - Fix error when window size > array size in rolling_apply (GH1850_)
  - Fix pip source installation issues via SSH from GitHub
  - Fix OLS.summary when column is a tuple (GH1837_)
  - Fix bug in __doc__ patching when -OO passed to interpreter
    (GH1792_ GH1741_ GH1774_)
  - Fix unicode console encoding issue in IPython notebook (GH1782_, GH1768_)
  - Fix unicode formatting issue with Series.name (GH1782_)
  - Fix bug in DataFrame.duplicated with datetime64 columns (GH1833_)
  - Fix bug in Panel internals resulting in error when doing fillna after
    truncate not changing size of panel (GH1823_)
  - Prevent segfault due to MultiIndex not being supported in HDFStore table
    format (GH1848_)
  - Fix UnboundLocalError in Panel.__setitem__ and add better error (GH1826_)
  - Fix to_csv issues with list of string entries. Isnull works on list of
    strings now too (GH1791_)
  - Fix Timestamp comparisons with datetime values outside the nanosecond range
    (1677-2262)
  - Revert to prior behavior of normalize_date with datetime.date objects
    (return datetime)
  - Fix broken interaction between np.nansum and Series.any/all
  - Fix bug with multiple column date parsers (GH1866_)
  - DatetimeIndex.union(Int64Index) was broken
  - Make plot x vs y interface consistent with integer indexing (GH1842_)
  - set_index inplace modified data even if unique check fails (GH1831_)
  - Only use Q-OCT/NOV/DEC in quarterly frequency inference (GH1789_)
  - Upcast to dtype=object when unstacking boolean DataFrame (GH1820_)
  - Fix float64/float32 merging bug (GH1849_)
  - Fixes to Period.start_time for non-daily frequencies (GH1857_)
  - Fix failure when converter used on index_col in read_csv (GH1835_)
  - Implement PeriodIndex.append so that pandas.concat works correctly (GH1815_)
  - Avoid Cython out-of-bounds access causing segfault sometimes in pad_2d,
    backfill_2d
  - Fix resampling error with intraday times and anchored target time (like
    AS-DEC) (GH1772_)
  - Fix .ix indexing bugs with mixed-integer indexes (GH1799_)
  - Respect passed color keyword argument in Series.plot (GH1890_)
  - Fix rolling_min/max when the window is larger than the size of the input
    array. Check other malformed inputs (GH1899_, GH1897_)
  - Rolling variance / standard deviation with only a single observation in
    window (GH1884_)
  - Fix unicode sheet name failure in to_excel (GH1828_)
  - Override DatetimeIndex.min/max to return Timestamp objects (GH1895_)
  - Fix column name formatting issue in length-truncated column (GH1906_)
  - Fix broken handling of copying Index metadata to new instances created by
    view(...) calls inside the NumPy infrastructure
  - Support datetime.date again in DateOffset.rollback/rollforward
  - Raise Exception if set passed to Series constructor (GH1913_)
  - Add TypeError when appending HDFStore table w/ wrong index type (GH1881_)
  - Don't raise exception on empty inputs in EW functions (e.g. ewma) (GH1900_)
  - Make asof work correctly with PeriodIndex (GH1883_)
  - Fix extlinks in doc build
  - Fill boolean DataFrame with NaN when calling shift (GH1814_)
  - Fix setuptools bug causing pip not to Cythonize .pyx files sometimes
  - Fix negative integer indexing regression in .ix from 0.7.x (GH1888_)
  - Fix error while retrieving timezone and utc offset from subclasses of
    datetime.tzinfo without .zone and ._utcoffset attributes (GH1922_)
  - Fix DataFrame formatting of small, non-zero FP numbers (GH1911_)
  - Various fixes by upcasting of date -> datetime (GH1395_)
  - Raise better exception when passing multiple functions with the same name,
    such as lambdas, to GroupBy.aggregate
  - Fix DataFrame.apply with axis=1 on a non-unique index (GH1878_)
  - Proper handling of Index subclasses in pandas.unique (GH1759_)
  - Set index names in DataFrame.from_records (GH1744_)
  - Fix time series indexing error with duplicates, under and over hash table
    size cutoff (GH1821_)
  - Handle list keys in addition to tuples in DataFrame.xs when
    partial-indexing a hierarchically-indexed DataFrame (GH1796_)
  - Support multiple column selection in DataFrame.__getitem__ with duplicate
    columns (GH1943_)
  - Fix time zone localization bug causing improper fields (e.g. hours) in time
    zones that have not had a UTC transition in a long time (GH1946_)
  - Fix errors when parsing and working with with fixed offset timezones
    (GH1922_, GH1928_)
  - Fix text parser bug when handling UTC datetime objects generated by
    dateutil (GH1693_)
  - Fix plotting bug when 'B' is the inferred frequency but index actually
    contains weekends (GH1668_, GH1669_)
  - Fix plot styling bugs (GH1666_, GH1665_, GH1658_)
  - Fix plotting bug with index/columns with unicode (GH1685_)
  - Fix DataFrame constructor bug when passed Series with datetime64 dtype
    in a dict (GH1680_)
  - Fixed regression in generating DatetimeIndex using timezone aware
    datetime.datetime (GH1676_)
  - Fix DataFrame bug when printing concatenated DataFrames with duplicated
    columns (GH1675_)
  - Fixed bug when plotting time series with multiple intraday frequencies
    (GH1732_)
  - Fix bug in DataFrame.duplicated to enable iterables other than list-types
    as input argument (GH1773_)
  - Fix resample bug when passed list of lambdas as `how` argument (GH1808_)
  - Repr fix for MultiIndex level with all NAs (GH1971_)
  - Fix PeriodIndex slicing bug when slice start/end are out-of-bounds (GH1977_)
  - Fix read_table bug when parsing unicode (GH1975_)
  - Fix BlockManager.iget bug when dealing with non-unique MultiIndex as columns
    (GH1970_)
  - Fix reset_index bug if both drop and level are specified (GH1957_)
  - Work around unsafe NumPy object->int casting with Cython function (GH1987_)
  - Fix datetime64 formatting bug in DataFrame.to_csv (GH1993_)
  - Default start date in pandas.io.data to 1/1/2000 as the docs say (GH2011_)


.. _GH1706: https://github.com/pydata/pandas/issues/1706
.. _GH1735: https://github.com/pydata/pandas/issues/1735
.. _GH1785: https://github.com/pydata/pandas/issues/1785
.. _GH1748: https://github.com/pydata/pandas/issues/1748
.. _GH1739: https://github.com/pydata/pandas/issues/1739
.. _GH1691: https://github.com/pydata/pandas/issues/1691
.. _GH1295: https://github.com/pydata/pandas/issues/1295
.. _GH1999: https://github.com/pydata/pandas/issues/1999
.. _GH1988: https://github.com/pydata/pandas/issues/1988
.. _GH1990: https://github.com/pydata/pandas/issues/1990
.. _GH1659: https://github.com/pydata/pandas/issues/1659
.. _GH1693: https://github.com/pydata/pandas/issues/1693
.. _GH1679: https://github.com/pydata/pandas/issues/1679
.. _GH1661: https://github.com/pydata/pandas/issues/1661
.. _GH1829: https://github.com/pydata/pandas/issues/1829
.. _GH1797: https://github.com/pydata/pandas/issues/1797
.. _GH1902: https://github.com/pydata/pandas/issues/1902
.. _GH1854: https://github.com/pydata/pandas/issues/1854
.. _GH1915: https://github.com/pydata/pandas/issues/1915
.. _GH1689: https://github.com/pydata/pandas/issues/1689
.. _GH1583: https://github.com/pydata/pandas/issues/1583
.. _GH924: https://github.com/pydata/pandas/issues/924
.. _GH1929: https://github.com/pydata/pandas/issues/1929
.. _GH1871: https://github.com/pydata/pandas/issues/1871
.. _GH1843: https://github.com/pydata/pandas/issues/1843
.. _GH2000: https://github.com/pydata/pandas/issues/2000
.. _GH1723: https://github.com/pydata/pandas/issues/1723
.. _GH1834: https://github.com/pydata/pandas/issues/1834
.. _GH1824: https://github.com/pydata/pandas/issues/1824
.. _GH1776: https://github.com/pydata/pandas/issues/1776
.. _GH1783: https://github.com/pydata/pandas/issues/1783
.. _GH1630: https://github.com/pydata/pandas/issues/1630
.. _GH1809: https://github.com/pydata/pandas/issues/1809
.. _GH1657: https://github.com/pydata/pandas/issues/1657
.. _GH1896: https://github.com/pydata/pandas/issues/1896
.. _GH1711: https://github.com/pydata/pandas/issues/1711
.. _GH1726: https://github.com/pydata/pandas/issues/1726
.. _GH1730: https://github.com/pydata/pandas/issues/1730
.. _GH1719: https://github.com/pydata/pandas/issues/1719
.. _GH1745: https://github.com/pydata/pandas/issues/1745
.. _GH1702: https://github.com/pydata/pandas/issues/1702
.. _GH1753: https://github.com/pydata/pandas/issues/1753
.. _GH1708: https://github.com/pydata/pandas/issues/1708
.. _GH1672: https://github.com/pydata/pandas/issues/1672
.. _GH1701: https://github.com/pydata/pandas/issues/1701
.. _GH1721: https://github.com/pydata/pandas/issues/1721
.. _GH1717: https://github.com/pydata/pandas/issues/1717
.. _GH1746: https://github.com/pydata/pandas/issues/1746
.. _GH1677: https://github.com/pydata/pandas/issues/1677
.. _GH1615: https://github.com/pydata/pandas/issues/1615
.. _GH1734: https://github.com/pydata/pandas/issues/1734
.. _GH1700: https://github.com/pydata/pandas/issues/1700
.. _GH1694: https://github.com/pydata/pandas/issues/1694
.. _GH1682: https://github.com/pydata/pandas/issues/1682
.. _GH1673: https://github.com/pydata/pandas/issues/1673
.. _GH1756: https://github.com/pydata/pandas/issues/1756
.. _GH1755: https://github.com/pydata/pandas/issues/1755
.. _GH1671: https://github.com/pydata/pandas/issues/1671
.. _GH1681: https://github.com/pydata/pandas/issues/1681
.. _GH1695: https://github.com/pydata/pandas/issues/1695
.. _GH1648: https://github.com/pydata/pandas/issues/1648
.. _GH1688: https://github.com/pydata/pandas/issues/1688
.. _GH1697: https://github.com/pydata/pandas/issues/1697
.. _GH1683: https://github.com/pydata/pandas/issues/1683
.. _GH1707: https://github.com/pydata/pandas/issues/1707
.. _GH1674: https://github.com/pydata/pandas/issues/1674
.. _GH1698: https://github.com/pydata/pandas/issues/1698
.. _GH1696: https://github.com/pydata/pandas/issues/1696
.. _GH1716: https://github.com/pydata/pandas/issues/1716
.. _GH1686: https://github.com/pydata/pandas/issues/1686
.. _GH1705: https://github.com/pydata/pandas/issues/1705
.. _GH1736: https://github.com/pydata/pandas/issues/1736
.. _GH1703: https://github.com/pydata/pandas/issues/1703
.. _GH1692: https://github.com/pydata/pandas/issues/1692
.. _GH1727: https://github.com/pydata/pandas/issues/1727
.. _GH1720: https://github.com/pydata/pandas/issues/1720
.. _GH1749: https://github.com/pydata/pandas/issues/1749
.. _GH1678: https://github.com/pydata/pandas/issues/1678
.. _GH1750: https://github.com/pydata/pandas/issues/1750
.. _GH1757: https://github.com/pydata/pandas/issues/1757
.. _GH1761: https://github.com/pydata/pandas/issues/1761
.. _GH1763: https://github.com/pydata/pandas/issues/1763
.. _GH1777: https://github.com/pydata/pandas/issues/1777
.. _GH1778: https://github.com/pydata/pandas/issues/1778
.. _GH1807: https://github.com/pydata/pandas/issues/1807
.. _GH1840: https://github.com/pydata/pandas/issues/1840
.. _GH1859: https://github.com/pydata/pandas/issues/1859
.. _GH1851: https://github.com/pydata/pandas/issues/1851
.. _GH1845: https://github.com/pydata/pandas/issues/1845
.. _GH1805: https://github.com/pydata/pandas/issues/1805
.. _GH1850: https://github.com/pydata/pandas/issues/1850
.. _GH1837: https://github.com/pydata/pandas/issues/1837
.. _GH1792: https://github.com/pydata/pandas/issues/1792
.. _GH1741: https://github.com/pydata/pandas/issues/1741
.. _GH1774: https://github.com/pydata/pandas/issues/1774
.. _GH1782: https://github.com/pydata/pandas/issues/1782
.. _GH1768: https://github.com/pydata/pandas/issues/1768
.. _GH1833: https://github.com/pydata/pandas/issues/1833
.. _GH1823: https://github.com/pydata/pandas/issues/1823
.. _GH1848: https://github.com/pydata/pandas/issues/1848
.. _GH1826: https://github.com/pydata/pandas/issues/1826
.. _GH1791: https://github.com/pydata/pandas/issues/1791
.. _GH1866: https://github.com/pydata/pandas/issues/1866
.. _GH1842: https://github.com/pydata/pandas/issues/1842
.. _GH1831: https://github.com/pydata/pandas/issues/1831
.. _GH1789: https://github.com/pydata/pandas/issues/1789
.. _GH1820: https://github.com/pydata/pandas/issues/1820
.. _GH1849: https://github.com/pydata/pandas/issues/1849
.. _GH1857: https://github.com/pydata/pandas/issues/1857
.. _GH1835: https://github.com/pydata/pandas/issues/1835
.. _GH1815: https://github.com/pydata/pandas/issues/1815
.. _GH1772: https://github.com/pydata/pandas/issues/1772
.. _GH1799: https://github.com/pydata/pandas/issues/1799
.. _GH1890: https://github.com/pydata/pandas/issues/1890
.. _GH1899: https://github.com/pydata/pandas/issues/1899
.. _GH1897: https://github.com/pydata/pandas/issues/1897
.. _GH1884: https://github.com/pydata/pandas/issues/1884
.. _GH1828: https://github.com/pydata/pandas/issues/1828
.. _GH1895: https://github.com/pydata/pandas/issues/1895
.. _GH1906: https://github.com/pydata/pandas/issues/1906
.. _GH1913: https://github.com/pydata/pandas/issues/1913
.. _GH1881: https://github.com/pydata/pandas/issues/1881
.. _GH1900: https://github.com/pydata/pandas/issues/1900
.. _GH1883: https://github.com/pydata/pandas/issues/1883
.. _GH1814: https://github.com/pydata/pandas/issues/1814
.. _GH1888: https://github.com/pydata/pandas/issues/1888
.. _GH1922: https://github.com/pydata/pandas/issues/1922
.. _GH1911: https://github.com/pydata/pandas/issues/1911
.. _GH1395: https://github.com/pydata/pandas/issues/1395
.. _GH1878: https://github.com/pydata/pandas/issues/1878
.. _GH1759: https://github.com/pydata/pandas/issues/1759
.. _GH1744: https://github.com/pydata/pandas/issues/1744
.. _GH1821: https://github.com/pydata/pandas/issues/1821
.. _GH1796: https://github.com/pydata/pandas/issues/1796
.. _GH1943: https://github.com/pydata/pandas/issues/1943
.. _GH1946: https://github.com/pydata/pandas/issues/1946
.. _GH1928: https://github.com/pydata/pandas/issues/1928
.. _GH1668: https://github.com/pydata/pandas/issues/1668
.. _GH1669: https://github.com/pydata/pandas/issues/1669
.. _GH1666: https://github.com/pydata/pandas/issues/1666
.. _GH1665: https://github.com/pydata/pandas/issues/1665
.. _GH1658: https://github.com/pydata/pandas/issues/1658
.. _GH1685: https://github.com/pydata/pandas/issues/1685
.. _GH1680: https://github.com/pydata/pandas/issues/1680
.. _GH1676: https://github.com/pydata/pandas/issues/1676
.. _GH1675: https://github.com/pydata/pandas/issues/1675
.. _GH1732: https://github.com/pydata/pandas/issues/1732
.. _GH1773: https://github.com/pydata/pandas/issues/1773
.. _GH1808: https://github.com/pydata/pandas/issues/1808
.. _GH1971: https://github.com/pydata/pandas/issues/1971
.. _GH1977: https://github.com/pydata/pandas/issues/1977
.. _GH1975: https://github.com/pydata/pandas/issues/1975
.. _GH1970: https://github.com/pydata/pandas/issues/1970
.. _GH1957: https://github.com/pydata/pandas/issues/1957
.. _GH1987: https://github.com/pydata/pandas/issues/1987
.. _GH1993: https://github.com/pydata/pandas/issues/1993
.. _GH2011: https://github.com/pydata/pandas/issues/2011


pandas 0.8.1
============

**Release date:** July 22, 2012

**New features**

  - Add vectorized, NA-friendly string methods to Series (GH1621_, GH620_)
  - Can pass dict of per-column line styles to DataFrame.plot (GH1559_)
  - Selective plotting to secondary y-axis on same subplot (GH1640_)
  - Add new ``bootstrap_plot`` plot function
  - Add new ``parallel_coordinates`` plot function (GH1488_)
  - Add ``radviz`` plot function (GH1566_)
  - Add ``multi_sparse`` option to ``set_printoptions`` to modify display of
    hierarchical indexes (GH1538_)
  - Add ``dropna`` method to Panel (GH171_)

**Improvements to existing features**

  - Use moving min/max algorithms from Bottleneck in rolling_min/rolling_max
    for > 100x speedup. (GH1504_, GH50_)
  - Add Cython group median method for >15x speedup (GH1358_)
  - Drastically improve ``to_datetime`` performance on ISO8601 datetime strings
    (with no time zones) (GH1571_)
  - Improve single-key groupby performance on large data sets, accelerate use of
    groupby with a Categorical variable
  - Add ability to append hierarchical index levels with ``set_index`` and to
    drop single levels with ``reset_index`` (GH1569_, GH1577_)
  - Always apply passed functions in ``resample``, even if upsampling (GH1596_)
  - Avoid unnecessary copies in DataFrame constructor with explicit dtype (GH1572_)
  - Cleaner DatetimeIndex string representation with 1 or 2 elements (GH1611_)
  - Improve performance of array-of-Period to PeriodIndex, convert such arrays
    to PeriodIndex inside Index (GH1215_)
  - More informative string representation for weekly Period objects (GH1503_)
  - Accelerate 3-axis multi data selection from homogeneous Panel (GH979_)
  - Add ``adjust`` option to ewma to disable adjustment factor (GH1584_)
  - Add new matplotlib converters for high frequency time series plotting (GH1599_)
  - Handling of tz-aware datetime.datetime objects in to_datetime; raise
    Exception unless utc=True given (GH1581_)

**Bug fixes**

  - Fix NA handling in DataFrame.to_panel (GH1582_)
  - Handle TypeError issues inside PyObject_RichCompareBool calls in khash
    (GH1318_)
  - Fix resampling bug to lower case daily frequency (GH1588_)
  - Fix kendall/spearman DataFrame.corr bug with no overlap (GH1595_)
  - Fix bug in DataFrame.set_index (GH1592_)
  - Don't ignore axes in boxplot if by specified (GH1565_)
  - Fix Panel .ix indexing with integers bug (GH1603_)
  - Fix Partial indexing bugs (years, months, ...) with PeriodIndex (GH1601_)
  - Fix MultiIndex console formatting issue (GH1606_)
  - Unordered index with duplicates doesn't yield scalar location for single
    entry (GH1586_)
  - Fix resampling of tz-aware time series with "anchored" freq (GH1591_)
  - Fix DataFrame.rank error on integer data (GH1589_)
  - Selection of multiple SparseDataFrame columns by list in __getitem__ (GH1585_)
  - Override Index.tolist for compatibility with MultiIndex (GH1576_)
  - Fix hierarchical summing bug with MultiIndex of length 1 (GH1568_)
  - Work around numpy.concatenate use/bug in Series.set_value (GH1561_)
  - Ensure Series/DataFrame are sorted before resampling (GH1580_)
  - Fix unhandled IndexError when indexing very large time series (GH1562_)
  - Fix DatetimeIndex intersection logic error with irregular indexes (GH1551_)
  - Fix unit test errors on Python 3 (GH1550_)
  - Fix .ix indexing bugs in duplicate DataFrame index (GH1201_)
  - Better handle errors with non-existing objects in HDFStore (GH1254_)
  - Don't copy int64 array data in DatetimeIndex when copy=False (GH1624_)
  - Fix resampling of conforming periods quarterly to annual (GH1622_)
  - Don't lose index name on resampling (GH1631_)
  - Support python-dateutil version 2.1 (GH1637_)
  - Fix broken scatter_matrix axis labeling, esp. with time series (GH1625_)
  - Fix cases where extra keywords weren't being passed on to matplotlib from
    Series.plot (GH1636_)
  - Fix BusinessMonthBegin logic for dates before 1st bday of month (GH1645_)
  - Ensure string alias converted (valid in DatetimeIndex.get_loc) in
    DataFrame.xs / __getitem__ (GH1644_)
  - Fix use of string alias timestamps with tz-aware time series (GH1647_)
  - Fix Series.max/min and Series.describe on len-0 series (GH1650_)
  - Handle None values in dict passed to concat (GH1649_)
  - Fix Series.interpolate with method='values' and DatetimeIndex (GH1646_)
  - Fix IndexError in left merges on a DataFrame with 0-length (GH1628_)
  - Fix DataFrame column width display with UTF-8 encoded characters (GH1620_)
  - Handle case in pandas.io.data.get_data_yahoo where Yahoo! returns duplicate
    dates for most recent business day
  - Avoid downsampling when plotting mixed frequencies on the same subplot (GH1619_)
  - Fix read_csv bug when reading a single line (GH1553_)
  - Fix bug in C code causing monthly periods prior to December 1969 to be off (GH1570_)

.. _GH1621: https://github.com/pydata/pandas/issues/1621
.. _GH620: https://github.com/pydata/pandas/issues/620
.. _GH1559: https://github.com/pydata/pandas/issues/1559
.. _GH1640: https://github.com/pydata/pandas/issues/1640
.. _GH1488: https://github.com/pydata/pandas/issues/1488
.. _GH1566: https://github.com/pydata/pandas/issues/1566
.. _GH1538: https://github.com/pydata/pandas/issues/1538
.. _GH171: https://github.com/pydata/pandas/issues/171
.. _GH1504: https://github.com/pydata/pandas/issues/1504
.. _GH50: https://github.com/pydata/pandas/issues/50
.. _GH1358: https://github.com/pydata/pandas/issues/1358
.. _GH1571: https://github.com/pydata/pandas/issues/1571
.. _GH1569: https://github.com/pydata/pandas/issues/1569
.. _GH1577: https://github.com/pydata/pandas/issues/1577
.. _GH1596: https://github.com/pydata/pandas/issues/1596
.. _GH1572: https://github.com/pydata/pandas/issues/1572
.. _GH1611: https://github.com/pydata/pandas/issues/1611
.. _GH1215: https://github.com/pydata/pandas/issues/1215
.. _GH1503: https://github.com/pydata/pandas/issues/1503
.. _GH979: https://github.com/pydata/pandas/issues/979
.. _GH1584: https://github.com/pydata/pandas/issues/1584
.. _GH1599: https://github.com/pydata/pandas/issues/1599
.. _GH1581: https://github.com/pydata/pandas/issues/1581
.. _GH1582: https://github.com/pydata/pandas/issues/1582
.. _GH1318: https://github.com/pydata/pandas/issues/1318
.. _GH1588: https://github.com/pydata/pandas/issues/1588
.. _GH1595: https://github.com/pydata/pandas/issues/1595
.. _GH1592: https://github.com/pydata/pandas/issues/1592
.. _GH1565: https://github.com/pydata/pandas/issues/1565
.. _GH1603: https://github.com/pydata/pandas/issues/1603
.. _GH1601: https://github.com/pydata/pandas/issues/1601
.. _GH1606: https://github.com/pydata/pandas/issues/1606
.. _GH1586: https://github.com/pydata/pandas/issues/1586
.. _GH1591: https://github.com/pydata/pandas/issues/1591
.. _GH1589: https://github.com/pydata/pandas/issues/1589
.. _GH1585: https://github.com/pydata/pandas/issues/1585
.. _GH1576: https://github.com/pydata/pandas/issues/1576
.. _GH1568: https://github.com/pydata/pandas/issues/1568
.. _GH1561: https://github.com/pydata/pandas/issues/1561
.. _GH1580: https://github.com/pydata/pandas/issues/1580
.. _GH1562: https://github.com/pydata/pandas/issues/1562
.. _GH1551: https://github.com/pydata/pandas/issues/1551
.. _GH1550: https://github.com/pydata/pandas/issues/1550
.. _GH1201: https://github.com/pydata/pandas/issues/1201
.. _GH1254: https://github.com/pydata/pandas/issues/1254
.. _GH1624: https://github.com/pydata/pandas/issues/1624
.. _GH1622: https://github.com/pydata/pandas/issues/1622
.. _GH1631: https://github.com/pydata/pandas/issues/1631
.. _GH1637: https://github.com/pydata/pandas/issues/1637
.. _GH1625: https://github.com/pydata/pandas/issues/1625
.. _GH1636: https://github.com/pydata/pandas/issues/1636
.. _GH1645: https://github.com/pydata/pandas/issues/1645
.. _GH1644: https://github.com/pydata/pandas/issues/1644
.. _GH1647: https://github.com/pydata/pandas/issues/1647
.. _GH1650: https://github.com/pydata/pandas/issues/1650
.. _GH1649: https://github.com/pydata/pandas/issues/1649
.. _GH1646: https://github.com/pydata/pandas/issues/1646
.. _GH1628: https://github.com/pydata/pandas/issues/1628
.. _GH1620: https://github.com/pydata/pandas/issues/1620
.. _GH1619: https://github.com/pydata/pandas/issues/1619
.. _GH1553: https://github.com/pydata/pandas/issues/1553
.. _GH1570: https://github.com/pydata/pandas/issues/1570


pandas 0.8.0
============

**Release date:** 6/29/2012

**New features**

  - New unified DatetimeIndex class for nanosecond-level timestamp data
  - New Timestamp datetime.datetime subclass with easy time zone conversions,
    and support for nanoseconds
  - New PeriodIndex class for timespans, calendar logic, and Period scalar object
  - High performance resampling of timestamp and period data. New `resample`
    method of all pandas data structures
  - New frequency names plus shortcut string aliases like '15h', '1h30min'
  - Time series string indexing shorthand (GH222_)
  - Add week, dayofyear array and other timestamp array-valued field accessor
    functions to DatetimeIndex
  - Add GroupBy.prod optimized aggregation function and 'prod' fast time series
    conversion method (GH1018_)
  - Implement robust frequency inference function and `inferred_freq` attribute
    on DatetimeIndex (GH391_)
  - New ``tz_convert`` and ``tz_localize`` methods in Series / DataFrame
  - Convert DatetimeIndexes to UTC if time zones are different in join/setops
    (GH864_)
  - Add limit argument for forward/backward filling to reindex, fillna,
    etc. (GH825_ and others)
  - Add support for indexes (dates or otherwise) with duplicates and common
    sense indexing/selection functionality
  - Series/DataFrame.update methods, in-place variant of combine_first (GH961_)
  - Add ``match`` function to API (GH502_)
  - Add Cython-optimized first, last, min, max, prod functions to GroupBy (GH994_,
    GH1043_)
  - Dates can be split across multiple columns (GH1227_, GH1186_)
  - Add experimental support for converting pandas DataFrame to R data.frame
    via rpy2 (GH350_, GH1212_)
  - Can pass list of (name, function) to GroupBy.aggregate to get aggregates in
    a particular order (GH610_)
  - Can pass dicts with lists of functions or dicts to GroupBy aggregate to do
    much more flexible multiple function aggregation (GH642_, GH610_)
  - New ordered_merge functions for merging DataFrames with ordered
    data. Also supports group-wise merging for panel data (GH813_)
  - Add keys() method to DataFrame
  - Add flexible replace method for replacing potentially values to Series and
    DataFrame (GH929_, GH1241_)
  - Add 'kde' plot kind for Series/DataFrame.plot (GH1059_)
  - More flexible multiple function aggregation with GroupBy
  - Add pct_change function to Series/DataFrame
  - Add option to interpolate by Index values in Series.interpolate (GH1206_)
  - Add ``max_colwidth`` option for DataFrame, defaulting to 50
  - Conversion of DataFrame through rpy2 to R data.frame (GH1282_, )
  - Add keys() method on DataFrame (GH1240_)
  - Add new ``match`` function to API (similar to R) (GH502_)
  - Add dayfirst option to parsers (GH854_)
  - Add ``method`` argument to ``align`` method for forward/backward fillin
    (GH216_)
  - Add Panel.transpose method for rearranging axes (GH695_)
  - Add new ``cut`` function (patterned after R) for discretizing data into
    equal range-length bins or arbitrary breaks of your choosing (GH415_)
  - Add new ``qcut`` for cutting with quantiles (GH1378_)
  - Add ``value_counts`` top level array method (GH1392_)
  - Added Andrews curves plot tupe (GH1325_)
  - Add lag plot (GH1440_)
  - Add autocorrelation_plot (GH1425_)
  - Add support for tox and Travis CI (GH1382_)
  - Add support for Categorical use in GroupBy (GH292_)
  - Add ``any`` and ``all`` methods to DataFrame (GH1416_)
  - Add ``secondary_y`` option to Series.plot
  - Add experimental ``lreshape`` function for reshaping wide to long

**Improvements to existing features**

  - Switch to klib/khash-based hash tables in Index classes for better
    performance in many cases and lower memory footprint
  - Shipping some functions from scipy.stats to reduce dependency,
    e.g. Series.describe and DataFrame.describe (GH1092_)
  - Can create MultiIndex by passing list of lists or list of arrays to Series,
    DataFrame constructor, etc. (GH831_)
  - Can pass arrays in addition to column names to DataFrame.set_index (GH402_)
  - Improve the speed of "square" reindexing of homogeneous DataFrame objects
    by significant margin (GH836_)
  - Handle more dtypes when passed MaskedArrays in DataFrame constructor (GH406_)
  - Improved performance of join operations on integer keys (GH682_)
  - Can pass multiple columns to GroupBy object, e.g. grouped[[col1, col2]] to
    only aggregate a subset of the value columns (GH383_)
  - Add histogram / kde plot options for scatter_matrix diagonals (GH1237_)
  - Add inplace option to Series/DataFrame.rename and sort_index,
    DataFrame.drop_duplicates (GH805_, GH207_)
  - More helpful error message when nothing passed to Series.reindex (GH1267_)
  - Can mix array and scalars as dict-value inputs to DataFrame ctor (GH1329_)
  - Use DataFrame columns' name for legend title in plots
  - Preserve frequency in DatetimeIndex when possible in boolean indexing
    operations
  - Promote datetime.date values in data alignment operations (GH867_)
  - Add ``order`` method to Index classes (GH1028_)
  - Avoid hash table creation in large monotonic hash table indexes (GH1160_)
  - Store time zones in HDFStore (GH1232_)
  - Enable storage of sparse data structures in HDFStore (#85)
  - Enable Series.asof to work with arrays of timestamp inputs
  - Cython implementation of DataFrame.corr speeds up by > 100x (GH1349_, GH1354_)
  - Exclude "nuisance" columns automatically in GroupBy.transform (GH1364_)
  - Support functions-as-strings in GroupBy.transform (GH1362_)
  - Use index name as xlabel/ylabel in plots (GH1415_)
  - Add ``convert_dtype`` option to Series.apply to be able to leave data as
    dtype=object (GH1414_)
  - Can specify all index level names in concat (GH1419_)
  - Add ``dialect`` keyword to parsers for quoting conventions (GH1363_)
  - Enable DataFrame[bool_DataFrame] += value (GH1366_)
  - Add ``retries`` argument to ``get_data_yahoo`` to try to prevent Yahoo! API
    404s (GH826_)
  - Improve performance of reshaping by using O(N) categorical sorting
  - Series names will be used for index of DataFrame if no index passed (GH1494_)
  - Header argument in DataFrame.to_csv can accept a list of column names to
    use instead of the object's columns (GH921_)
  - Add ``raise_conflict`` argument to DataFrame.update (GH1526_)
  - Support file-like objects in ExcelFile (GH1529_)

**API Changes**

  - Rename `pandas._tseries` to `pandas.lib`
  - Rename Factor to Categorical and add improvements. Numerous Categorical bug
    fixes
  - Frequency name overhaul, WEEKDAY/EOM and rules with @
    deprecated. get_legacy_offset_name backwards compatibility function added
  - Raise ValueError in DataFrame.__nonzero__, so "if df" no longer works
    (GH1073_)
  - Change BDay (business day) to not normalize dates by default (GH506_)
  - Remove deprecated DataMatrix name
  - Default merge suffixes for overlap now have underscores instead of periods
    to facilitate tab completion, etc. (GH1239_)
  - Deprecation of offset, time_rule timeRule parameters throughout codebase
  - Series.append and DataFrame.append no longer check for duplicate indexes
    by default, add verify_integrity parameter (GH1394_)
  - Refactor Factor class, old constructor moved to Factor.from_array
  - Modified internals of MultiIndex to use less memory (no longer represented
    as array of tuples) internally, speed up construction time and many methods
    which construct intermediate hierarchical indexes (GH1467_)

**Bug fixes**

  - Fix OverflowError from storing pre-1970 dates in HDFStore by switching to
    datetime64 (GH179_)
  - Fix logical error with February leap year end in YearEnd offset
  - Series([False, nan]) was getting casted to float64 (GH1074_)
  - Fix binary operations between boolean Series and object Series with
    booleans and NAs (GH1074_, GH1079_)
  - Couldn't assign whole array to column in mixed-type DataFrame via .ix
    (GH1142_)
  - Fix label slicing issues with float index values (GH1167_)
  - Fix segfault caused by empty groups passed to groupby (GH1048_)
  - Fix occasionally misbehaved reindexing in the presence of NaN labels (GH522_)
  - Fix imprecise logic causing weird Series results from .apply (GH1183_)
  - Unstack multiple levels in one shot, avoiding empty columns in some
    cases. Fix pivot table bug (GH1181_)
  - Fix formatting of MultiIndex on Series/DataFrame when index name coincides
    with label (GH1217_)
  - Handle Excel 2003 #N/A as NaN from xlrd (GH1213_, GH1225_)
  - Fix timestamp locale-related deserialization issues with HDFStore by moving
    to datetime64 representation (GH1081_, GH809_)
  - Fix DataFrame.duplicated/drop_duplicates NA value handling (GH557_)
  - Actually raise exceptions in fast reducer (GH1243_)
  - Fix various timezone-handling bugs from 0.7.3 (GH969_)
  - GroupBy on level=0 discarded index name (GH1313_)
  - Better error message with unmergeable DataFrames (GH1307_)
  - Series.__repr__ alignment fix with unicode index values (GH1279_)
  - Better error message if nothing passed to reindex (GH1267_)
  - More robust NA handling in DataFrame.drop_duplicates (GH557_)
  - Resolve locale-based and pre-epoch HDF5 timestamp deserialization issues
    (GH973_, GH1081_, GH179_)
  - Implement Series.repeat (GH1229_)
  - Fix indexing with namedtuple and other tuple subclasses (GH1026_)
  - Fix float64 slicing bug (GH1167_)
  - Parsing integers with commas (GH796_)
  - Fix groupby improper data type when group consists of one value (GH1065_)
  - Fix negative variance possibility in nanvar resulting from floating point
    error (GH1090_)
  - Consistently set name on groupby pieces (GH184_)
  - Treat dict return values as Series in GroupBy.apply (GH823_)
  - Respect column selection for DataFrame in in GroupBy.transform (GH1365_)
  - Fix MultiIndex partial indexing bug (GH1352_)
  - Enable assignment of rows in mixed-type DataFrame via .ix (GH1432_)
  - Reset index mapping when grouping Series in Cython (GH1423_)
  - Fix outer/inner DataFrame.join with non-unique indexes (GH1421_)
  - Fix MultiIndex groupby bugs with empty lower levels (GH1401_)
  - Calling fillna with a Series will have same behavior as with dict (GH1486_)
  - SparseSeries reduction bug (GH1375_)
  - Fix unicode serialization issue in HDFStore (GH1361_)
  - Pass keywords to pyplot.boxplot in DataFrame.boxplot (GH1493_)
  - Bug fixes in MonthBegin (GH1483_)
  - Preserve MultiIndex names in drop (GH1513_)
  - Fix Panel DataFrame slice-assignment bug (GH1533_)
  - Don't use locals() in read_* functions (GH1547_)

.. _GH222: https://github.com/pydata/pandas/issues/222
.. _GH1018: https://github.com/pydata/pandas/issues/1018
.. _GH391: https://github.com/pydata/pandas/issues/391
.. _GH864: https://github.com/pydata/pandas/issues/864
.. _GH825: https://github.com/pydata/pandas/issues/825
.. _GH961: https://github.com/pydata/pandas/issues/961
.. _GH502: https://github.com/pydata/pandas/issues/502
.. _GH994: https://github.com/pydata/pandas/issues/994
.. _GH1043: https://github.com/pydata/pandas/issues/1043
.. _GH1227: https://github.com/pydata/pandas/issues/1227
.. _GH1186: https://github.com/pydata/pandas/issues/1186
.. _GH350: https://github.com/pydata/pandas/issues/350
.. _GH1212: https://github.com/pydata/pandas/issues/1212
.. _GH610: https://github.com/pydata/pandas/issues/610
.. _GH642: https://github.com/pydata/pandas/issues/642
.. _GH813: https://github.com/pydata/pandas/issues/813
.. _GH929: https://github.com/pydata/pandas/issues/929
.. _GH1241: https://github.com/pydata/pandas/issues/1241
.. _GH1059: https://github.com/pydata/pandas/issues/1059
.. _GH1206: https://github.com/pydata/pandas/issues/1206
.. _GH1282: https://github.com/pydata/pandas/issues/1282
.. _GH1240: https://github.com/pydata/pandas/issues/1240
.. _GH854: https://github.com/pydata/pandas/issues/854
.. _GH216: https://github.com/pydata/pandas/issues/216
.. _GH695: https://github.com/pydata/pandas/issues/695
.. _GH415: https://github.com/pydata/pandas/issues/415
.. _GH1378: https://github.com/pydata/pandas/issues/1378
.. _GH1392: https://github.com/pydata/pandas/issues/1392
.. _GH1325: https://github.com/pydata/pandas/issues/1325
.. _GH1440: https://github.com/pydata/pandas/issues/1440
.. _GH1425: https://github.com/pydata/pandas/issues/1425
.. _GH1382: https://github.com/pydata/pandas/issues/1382
.. _GH292: https://github.com/pydata/pandas/issues/292
.. _GH1416: https://github.com/pydata/pandas/issues/1416
.. _GH1092: https://github.com/pydata/pandas/issues/1092
.. _GH831: https://github.com/pydata/pandas/issues/831
.. _GH402: https://github.com/pydata/pandas/issues/402
.. _GH836: https://github.com/pydata/pandas/issues/836
.. _GH406: https://github.com/pydata/pandas/issues/406
.. _GH682: https://github.com/pydata/pandas/issues/682
.. _GH383: https://github.com/pydata/pandas/issues/383
.. _GH1237: https://github.com/pydata/pandas/issues/1237
.. _GH805: https://github.com/pydata/pandas/issues/805
.. _GH207: https://github.com/pydata/pandas/issues/207
.. _GH1267: https://github.com/pydata/pandas/issues/1267
.. _GH1329: https://github.com/pydata/pandas/issues/1329
.. _GH867: https://github.com/pydata/pandas/issues/867
.. _GH1028: https://github.com/pydata/pandas/issues/1028
.. _GH1160: https://github.com/pydata/pandas/issues/1160
.. _GH1232: https://github.com/pydata/pandas/issues/1232
.. _GH1349: https://github.com/pydata/pandas/issues/1349
.. _GH1354: https://github.com/pydata/pandas/issues/1354
.. _GH1364: https://github.com/pydata/pandas/issues/1364
.. _GH1362: https://github.com/pydata/pandas/issues/1362
.. _GH1415: https://github.com/pydata/pandas/issues/1415
.. _GH1414: https://github.com/pydata/pandas/issues/1414
.. _GH1419: https://github.com/pydata/pandas/issues/1419
.. _GH1363: https://github.com/pydata/pandas/issues/1363
.. _GH1366: https://github.com/pydata/pandas/issues/1366
.. _GH826: https://github.com/pydata/pandas/issues/826
.. _GH1494: https://github.com/pydata/pandas/issues/1494
.. _GH921: https://github.com/pydata/pandas/issues/921
.. _GH1526: https://github.com/pydata/pandas/issues/1526
.. _GH1529: https://github.com/pydata/pandas/issues/1529
.. _GH1073: https://github.com/pydata/pandas/issues/1073
.. _GH506: https://github.com/pydata/pandas/issues/506
.. _GH1239: https://github.com/pydata/pandas/issues/1239
.. _GH1394: https://github.com/pydata/pandas/issues/1394
.. _GH1467: https://github.com/pydata/pandas/issues/1467
.. _GH179: https://github.com/pydata/pandas/issues/179
.. _GH1074: https://github.com/pydata/pandas/issues/1074
.. _GH1079: https://github.com/pydata/pandas/issues/1079
.. _GH1142: https://github.com/pydata/pandas/issues/1142
.. _GH1167: https://github.com/pydata/pandas/issues/1167
.. _GH1048: https://github.com/pydata/pandas/issues/1048
.. _GH522: https://github.com/pydata/pandas/issues/522
.. _GH1183: https://github.com/pydata/pandas/issues/1183
.. _GH1181: https://github.com/pydata/pandas/issues/1181
.. _GH1217: https://github.com/pydata/pandas/issues/1217
.. _GH1213: https://github.com/pydata/pandas/issues/1213
.. _GH1225: https://github.com/pydata/pandas/issues/1225
.. _GH1081: https://github.com/pydata/pandas/issues/1081
.. _GH809: https://github.com/pydata/pandas/issues/809
.. _GH557: https://github.com/pydata/pandas/issues/557
.. _GH1243: https://github.com/pydata/pandas/issues/1243
.. _GH969: https://github.com/pydata/pandas/issues/969
.. _GH1313: https://github.com/pydata/pandas/issues/1313
.. _GH1307: https://github.com/pydata/pandas/issues/1307
.. _GH1279: https://github.com/pydata/pandas/issues/1279
.. _GH973: https://github.com/pydata/pandas/issues/973
.. _GH1229: https://github.com/pydata/pandas/issues/1229
.. _GH1026: https://github.com/pydata/pandas/issues/1026
.. _GH796: https://github.com/pydata/pandas/issues/796
.. _GH1065: https://github.com/pydata/pandas/issues/1065
.. _GH1090: https://github.com/pydata/pandas/issues/1090
.. _GH184: https://github.com/pydata/pandas/issues/184
.. _GH823: https://github.com/pydata/pandas/issues/823
.. _GH1365: https://github.com/pydata/pandas/issues/1365
.. _GH1352: https://github.com/pydata/pandas/issues/1352
.. _GH1432: https://github.com/pydata/pandas/issues/1432
.. _GH1423: https://github.com/pydata/pandas/issues/1423
.. _GH1421: https://github.com/pydata/pandas/issues/1421
.. _GH1401: https://github.com/pydata/pandas/issues/1401
.. _GH1486: https://github.com/pydata/pandas/issues/1486
.. _GH1375: https://github.com/pydata/pandas/issues/1375
.. _GH1361: https://github.com/pydata/pandas/issues/1361
.. _GH1493: https://github.com/pydata/pandas/issues/1493
.. _GH1483: https://github.com/pydata/pandas/issues/1483
.. _GH1513: https://github.com/pydata/pandas/issues/1513
.. _GH1533: https://github.com/pydata/pandas/issues/1533
.. _GH1547: https://github.com/pydata/pandas/issues/1547


pandas 0.7.3
============

**Release date:** April 12, 2012

**New features / modules**

  - Support for non-unique indexes: indexing and selection, many-to-one and
    many-to-many joins (GH1306_)
  - Added fixed-width file reader, read_fwf (GH952_)
  - Add group_keys argument to groupby to not add group names to MultiIndex in
    result of apply (GH938_)
  - DataFrame can now accept non-integer label slicing (GH946_). Previously
    only DataFrame.ix was able to do so.
  - DataFrame.apply now retains name attributes on Series objects (GH983_)
  - Numeric DataFrame comparisons with non-numeric values now raises proper
    TypeError (GH943_). Previously raise "PandasError: DataFrame constructor
    not properly called!"
  - Add ``kurt`` methods to Series and DataFrame (GH964_)
  - Can pass dict of column -> list/set NA values for text parsers (GH754_)
  - Allows users specified NA values in text parsers (GH754_)
  - Parsers checks for openpyxl dependency and raises ImportError if not found
    (GH1007_)
  - New factory function to create HDFStore objects that can be used in a with
    statement so users do not have to explicitly call HDFStore.close (GH1005_)
  - pivot_table is now more flexible with same parameters as groupby (GH941_)
  - Added stacked bar plots (GH987_)
  - scatter_matrix method in pandas/tools/plotting.py (GH935_)
  - DataFrame.boxplot returns plot results for ex-post styling (GH985_)
  - Short version number accessible as pandas.version.short_version (GH930_)
  - Additional documentation in panel.to_frame (GH942_)
  - More informative Series.apply docstring regarding element-wise apply
    (GH977_)
  - Notes on rpy2 installation (GH1006_)
  - Add rotation and font size options to hist method (GH1012_)
  - Use exogenous / X variable index in result of OLS.y_predict. Add
    OLS.predict method (GH1027_, GH1008_)

**API Changes**

  - Calling apply on grouped Series, e.g. describe(), will no longer yield
    DataFrame by default. Will have to call unstack() to get prior behavior
  - NA handling in non-numeric comparisons has been tightened up (GH933_, GH953_)
  - No longer assign dummy names key_0, key_1, etc. to groupby index (GH1291_)

**Bug fixes**

  - Fix logic error when selecting part of a row in a DataFrame with a
    MultiIndex index (GH1013_)
  - Series comparison with Series of differing length causes crash (GH1016_).
  - Fix bug in indexing when selecting section of hierarchically-indexed row
    (GH1013_)
  - DataFrame.plot(logy=True) has no effect (GH1011_).
  - Broken arithmetic operations between SparsePanel-Panel (GH1015_)
  - Unicode repr issues in MultiIndex with non-ascii characters (GH1010_)
  - DataFrame.lookup() returns inconsistent results if exact match not present
    (GH1001_)
  - DataFrame arithmetic operations not treating None as NA (GH992_)
  - DataFrameGroupBy.apply returns incorrect result (GH991_)
  - Series.reshape returns incorrect result for multiple dimensions (GH989_)
  - Series.std and Series.var ignores ddof parameter (GH934_)
  - DataFrame.append loses index names (GH980_)
  - DataFrame.plot(kind='bar') ignores color argument (GH958_)
  - Inconsistent Index comparison results (GH948_)
  - Improper int dtype DataFrame construction from data with NaN (GH846_)
  - Removes default 'result' name in grouby results (GH995_)
  - DataFrame.from_records no longer mutate input columns (GH975_)
  - Use Index name when grouping by it (GH1313_)

.. _GH1306: https://github.com/pydata/pandas/issues/1306
.. _GH952: https://github.com/pydata/pandas/issues/952
.. _GH938: https://github.com/pydata/pandas/issues/938
.. _GH946: https://github.com/pydata/pandas/issues/946
.. _GH983: https://github.com/pydata/pandas/issues/983
.. _GH943: https://github.com/pydata/pandas/issues/943
.. _GH964: https://github.com/pydata/pandas/issues/964
.. _GH754: https://github.com/pydata/pandas/issues/754
.. _GH1007: https://github.com/pydata/pandas/issues/1007
.. _GH1005: https://github.com/pydata/pandas/issues/1005
.. _GH941: https://github.com/pydata/pandas/issues/941
.. _GH987: https://github.com/pydata/pandas/issues/987
.. _GH935: https://github.com/pydata/pandas/issues/935
.. _GH985: https://github.com/pydata/pandas/issues/985
.. _GH930: https://github.com/pydata/pandas/issues/930
.. _GH942: https://github.com/pydata/pandas/issues/942
.. _GH977: https://github.com/pydata/pandas/issues/977
.. _GH1006: https://github.com/pydata/pandas/issues/1006
.. _GH1012: https://github.com/pydata/pandas/issues/1012
.. _GH1027: https://github.com/pydata/pandas/issues/1027
.. _GH1008: https://github.com/pydata/pandas/issues/1008
.. _GH933: https://github.com/pydata/pandas/issues/933
.. _GH953: https://github.com/pydata/pandas/issues/953
.. _GH1291: https://github.com/pydata/pandas/issues/1291
.. _GH1013: https://github.com/pydata/pandas/issues/1013
.. _GH1016: https://github.com/pydata/pandas/issues/1016
.. _GH1011: https://github.com/pydata/pandas/issues/1011
.. _GH1015: https://github.com/pydata/pandas/issues/1015
.. _GH1010: https://github.com/pydata/pandas/issues/1010
.. _GH1001: https://github.com/pydata/pandas/issues/1001
.. _GH992: https://github.com/pydata/pandas/issues/992
.. _GH991: https://github.com/pydata/pandas/issues/991
.. _GH989: https://github.com/pydata/pandas/issues/989
.. _GH934: https://github.com/pydata/pandas/issues/934
.. _GH980: https://github.com/pydata/pandas/issues/980
.. _GH958: https://github.com/pydata/pandas/issues/958
.. _GH948: https://github.com/pydata/pandas/issues/948
.. _GH846: https://github.com/pydata/pandas/issues/846
.. _GH995: https://github.com/pydata/pandas/issues/995
.. _GH975: https://github.com/pydata/pandas/issues/975
.. _GH1313: https://github.com/pydata/pandas/issues/1313


pandas 0.7.2
============

**Release date:** March 16, 2012

**New features / modules**

  - Add additional tie-breaking methods in DataFrame.rank (GH874_)
  - Add ascending parameter to rank in Series, DataFrame (GH875_)
  - Add coerce_float option to DataFrame.from_records (GH893_)
  - Add sort_columns parameter to allow unsorted plots (GH918_)
  - IPython tab completion on GroupBy objects

**API Changes**

  - Series.sum returns 0 instead of NA when called on an empty
    series. Analogously for a DataFrame whose rows or columns are length 0
    (GH844_)

**Improvements to existing features**

  - Don't use groups dict in Grouper.size (GH860_)
  - Use khash for Series.value_counts, add raw function to algorithms.py (GH861_)
  - Enable column access via attributes on GroupBy (GH882_)
  - Enable setting existing columns (only) via attributes on DataFrame, Panel
    (GH883_)
  - Intercept __builtin__.sum in groupby (GH885_)
  - Can pass dict to DataFrame.fillna to use different values per column (GH661_)
  - Can select multiple hierarchical groups by passing list of values in .ix
    (GH134_)
  - Add level keyword to ``drop`` for dropping values from a level (GH159_)
  - Add ``coerce_float`` option on DataFrame.from_records (# 893)
  - Raise exception if passed date_parser fails in ``read_csv``
  - Add ``axis`` option to DataFrame.fillna (GH174_)
  - Fixes to Panel to make it easier to subclass (GH888_)

**Bug fixes**

  - Fix overflow-related bugs in groupby (GH850_, GH851_)
  - Fix unhelpful error message in parsers (GH856_)
  - Better err msg for failed boolean slicing of dataframe (GH859_)
  - Series.count cannot accept a string (level name) in the level argument (GH869_)
  - Group index platform int check (GH870_)
  - concat on axis=1 and ignore_index=True raises TypeError (GH871_)
  - Further unicode handling issues resolved (GH795_)
  - Fix failure in multiindex-based access in Panel (GH880_)
  - Fix DataFrame boolean slice assignment failure (GH881_)
  - Fix combineAdd NotImplementedError for SparseDataFrame (GH887_)
  - Fix DataFrame.to_html encoding and columns (GH890_, GH891_, GH909_)
  - Fix na-filling handling in mixed-type DataFrame (GH910_)
  - Fix to DataFrame.set_value with non-existant row/col (GH911_)
  - Fix malformed block in groupby when excluding nuisance columns (GH916_)
  - Fix inconsistant NA handling in dtype=object arrays (GH925_)
  - Fix missing center-of-mass computation in ewmcov (GH862_)
  - Don't raise exception when opening read-only HDF5 file (GH847_)
  - Fix possible out-of-bounds memory access in 0-length Series (GH917_)

.. _GH874: https://github.com/pydata/pandas/issues/874
.. _GH875: https://github.com/pydata/pandas/issues/875
.. _GH893: https://github.com/pydata/pandas/issues/893
.. _GH918: https://github.com/pydata/pandas/issues/918
.. _GH844: https://github.com/pydata/pandas/issues/844
.. _GH860: https://github.com/pydata/pandas/issues/860
.. _GH861: https://github.com/pydata/pandas/issues/861
.. _GH882: https://github.com/pydata/pandas/issues/882
.. _GH883: https://github.com/pydata/pandas/issues/883
.. _GH885: https://github.com/pydata/pandas/issues/885
.. _GH661: https://github.com/pydata/pandas/issues/661
.. _GH134: https://github.com/pydata/pandas/issues/134
.. _GH159: https://github.com/pydata/pandas/issues/159
.. _GH174: https://github.com/pydata/pandas/issues/174
.. _GH888: https://github.com/pydata/pandas/issues/888
.. _GH850: https://github.com/pydata/pandas/issues/850
.. _GH851: https://github.com/pydata/pandas/issues/851
.. _GH856: https://github.com/pydata/pandas/issues/856
.. _GH859: https://github.com/pydata/pandas/issues/859
.. _GH869: https://github.com/pydata/pandas/issues/869
.. _GH870: https://github.com/pydata/pandas/issues/870
.. _GH871: https://github.com/pydata/pandas/issues/871
.. _GH795: https://github.com/pydata/pandas/issues/795
.. _GH880: https://github.com/pydata/pandas/issues/880
.. _GH881: https://github.com/pydata/pandas/issues/881
.. _GH887: https://github.com/pydata/pandas/issues/887
.. _GH890: https://github.com/pydata/pandas/issues/890
.. _GH891: https://github.com/pydata/pandas/issues/891
.. _GH909: https://github.com/pydata/pandas/issues/909
.. _GH910: https://github.com/pydata/pandas/issues/910
.. _GH911: https://github.com/pydata/pandas/issues/911
.. _GH916: https://github.com/pydata/pandas/issues/916
.. _GH925: https://github.com/pydata/pandas/issues/925
.. _GH862: https://github.com/pydata/pandas/issues/862
.. _GH847: https://github.com/pydata/pandas/issues/847
.. _GH917: https://github.com/pydata/pandas/issues/917


pandas 0.7.1
============

**Release date:** February 29, 2012

**New features / modules**

  - Add ``to_clipboard`` function to pandas namespace for writing objects to
    the system clipboard (GH774_)
  - Add ``itertuples`` method to DataFrame for iterating through the rows of a
    dataframe as tuples (GH818_)
  - Add ability to pass fill_value and method to DataFrame and Series align
    method (GH806_, GH807_)
  - Add fill_value option to reindex, align methods (GH784_)
  - Enable concat to produce DataFrame from Series (GH787_)
  - Add ``between`` method to Series (GH802_)
  - Add HTML representation hook to DataFrame for the IPython HTML notebook
    (GH773_)
  - Support for reading Excel 2007 XML documents using openpyxl

**Improvements to existing features**

  - Improve performance and memory usage of fillna on DataFrame
  - Can concatenate a list of Series along axis=1 to obtain a DataFrame (GH787_)

**Bug fixes**

  - Fix memory leak when inserting large number of columns into a single
    DataFrame (GH790_)
  - Appending length-0 DataFrame with new columns would not result in those new
    columns being part of the resulting concatenated DataFrame (GH782_)
  - Fixed groupby corner case when passing dictionary grouper and as_index is
    False (GH819_)
  - Fixed bug whereby bool array sometimes had object dtype (GH820_)
  - Fix exception thrown on np.diff (GH816_)
  - Fix to_records where columns are non-strings (GH822_)
  - Fix Index.intersection where indices have incomparable types (GH811_)
  - Fix ExcelFile throwing an exception for two-line file (GH837_)
  - Add clearer error message in csv parser (GH835_)
  - Fix loss of fractional seconds in HDFStore (GH513_)
  - Fix DataFrame join where columns have datetimes (GH787_)
  - Work around numpy performance issue in take (GH817_)
  - Improve comparison operations for NA-friendliness (GH801_)
  - Fix indexing operation for floating point values (GH780_, GH798_)
  - Fix groupby case resulting in malformed dataframe (GH814_)
  - Fix behavior of reindex of Series dropping name (GH812_)
  - Improve on redudant groupby computation (GH775_)
  - Catch possible NA assignment to int/bool series with exception (GH839_)

.. _GH774: https://github.com/pydata/pandas/issues/774
.. _GH818: https://github.com/pydata/pandas/issues/818
.. _GH806: https://github.com/pydata/pandas/issues/806
.. _GH807: https://github.com/pydata/pandas/issues/807
.. _GH784: https://github.com/pydata/pandas/issues/784
.. _GH787: https://github.com/pydata/pandas/issues/787
.. _GH802: https://github.com/pydata/pandas/issues/802
.. _GH773: https://github.com/pydata/pandas/issues/773
.. _GH790: https://github.com/pydata/pandas/issues/790
.. _GH782: https://github.com/pydata/pandas/issues/782
.. _GH819: https://github.com/pydata/pandas/issues/819
.. _GH820: https://github.com/pydata/pandas/issues/820
.. _GH816: https://github.com/pydata/pandas/issues/816
.. _GH822: https://github.com/pydata/pandas/issues/822
.. _GH811: https://github.com/pydata/pandas/issues/811
.. _GH837: https://github.com/pydata/pandas/issues/837
.. _GH835: https://github.com/pydata/pandas/issues/835
.. _GH513: https://github.com/pydata/pandas/issues/513
.. _GH817: https://github.com/pydata/pandas/issues/817
.. _GH801: https://github.com/pydata/pandas/issues/801
.. _GH780: https://github.com/pydata/pandas/issues/780
.. _GH798: https://github.com/pydata/pandas/issues/798
.. _GH814: https://github.com/pydata/pandas/issues/814
.. _GH812: https://github.com/pydata/pandas/issues/812
.. _GH775: https://github.com/pydata/pandas/issues/775
.. _GH839: https://github.com/pydata/pandas/issues/839


pandas 0.7.0
============

**Release date:** 2/9/2012

**New features / modules**

  - New ``merge`` function for efficiently performing full gamut of database /
    relational-algebra operations. Refactored existing join methods to use the
    new infrastructure, resulting in substantial performance gains (GH220_,
    GH249_, GH267_)
  - New ``concat`` function for concatenating DataFrame or Panel objects along
    an axis. Can form union or intersection of the other axes. Improves
    performance of ``DataFrame.append`` (GH468_, GH479_, GH273_)
  - Handle differently-indexed output values in ``DataFrame.apply`` (GH498_)
  - Can pass list of dicts (e.g., a list of shallow JSON objects) to DataFrame
    constructor (GH526_)
  - Add ``reorder_levels`` method to Series and DataFrame (GH534_)
  - Add dict-like ``get`` function to DataFrame and Panel (GH521_)
  - ``DataFrame.iterrows`` method for efficiently iterating through the rows of
    a DataFrame
  - Added ``DataFrame.to_panel`` with code adapted from ``LongPanel.to_long``
  - ``reindex_axis`` method added to DataFrame
  - Add ``level`` option to binary arithmetic functions on ``DataFrame`` and
    ``Series``
  - Add ``level`` option to the ``reindex`` and ``align`` methods on Series and
    DataFrame for broadcasting values across a level (GH542_, GH552_, others)
  - Add attribute-based item access to ``Panel`` and add IPython completion (PR
    GH554_)
  - Add ``logy`` option to ``Series.plot`` for log-scaling on the Y axis
  - Add ``index``, ``header``, and ``justify`` options to
    ``DataFrame.to_string``. Add option to   (GH570_, GH571_)
  - Can pass multiple DataFrames to ``DataFrame.join`` to join on index (GH115_)
  - Can pass multiple Panels to ``Panel.join`` (GH115_)
  - Can pass multiple DataFrames to `DataFrame.append` to concatenate (stack)
    and multiple Series to ``Series.append`` too
  - Added ``justify`` argument to ``DataFrame.to_string`` to allow different
    alignment of column headers
  - Add ``sort`` option to GroupBy to allow disabling sorting of the group keys
    for potential speedups (GH595_)
  - Can pass MaskedArray to Series constructor (GH563_)
  - Add Panel item access via attributes and IPython completion (GH554_)
  - Implement ``DataFrame.lookup``, fancy-indexing analogue for retrieving
    values given a sequence of row and column labels (GH338_)
  - Add ``verbose`` option to ``read_csv`` and ``read_table`` to show number of
    NA values inserted in non-numeric columns (GH614_)
  - Can pass a list of dicts or Series to ``DataFrame.append`` to concatenate
    multiple rows (GH464_)
  - Add ``level`` argument to ``DataFrame.xs`` for selecting data from other
    MultiIndex levels. Can take one or more levels with potentially a tuple of
    keys for flexible retrieval of data (GH371_, GH629_)
  - New ``crosstab`` function for easily computing frequency tables (GH170_)
  - Can pass a list of functions to aggregate with groupby on a DataFrame,
    yielding an aggregated result with hierarchical columns (GH166_)
  - Add integer-indexing functions ``iget`` in Series and ``irow`` / ``iget``
    in DataFrame (GH628_)
  - Add new ``Series.unique`` function, significantly faster than
    ``numpy.unique`` (GH658_)
  - Add new ``cummin`` and ``cummax`` instance methods to ``Series`` and
    ``DataFrame`` (GH647_)
  - Add new ``value_range`` function to return min/max of a dataframe (GH288_)
  - Add ``drop`` parameter to ``reset_index`` method of ``DataFrame`` and added
    method to ``Series`` as well (GH699_)
  - Add ``isin`` method to Index objects, works just like ``Series.isin`` (GH
    GH657_)
  - Implement array interface on Panel so that ufuncs work (re: GH740_)
  - Add ``sort`` option to ``DataFrame.join`` (GH731_)
  - Improved handling of NAs (propagation) in binary operations with
    dtype=object arrays (GH737_)
  - Add ``abs`` method to Pandas objects
  - Added ``algorithms`` module to start collecting central algos

**API Changes**

  - Label-indexing with integer indexes now raises KeyError if a label is not
    found instead of falling back on location-based indexing (GH700_)
  - Label-based slicing via ``ix`` or ``[]`` on Series will now only work if
    exact matches for the labels are found or if the index is monotonic (for
    range selections)
  - Label-based slicing and sequences of labels can be passed to ``[]`` on a
    Series for both getting and setting (GH #86)
  - `[]` operator (``__getitem__`` and ``__setitem__``) will raise KeyError
    with integer indexes when an index is not contained in the index. The prior
    behavior would fall back on position-based indexing if a key was not found
    in the index which would lead to subtle bugs. This is now consistent with
    the behavior of ``.ix`` on DataFrame and friends (GH328_)
  - Rename ``DataFrame.delevel`` to ``DataFrame.reset_index`` and add
    deprecation warning
  - `Series.sort` (an in-place operation) called on a Series which is a view on
    a larger array (e.g. a column in a DataFrame) will generate an Exception to
    prevent accidentally modifying the data source (GH316_)
  - Refactor to remove deprecated ``LongPanel`` class (GH552_)
  - Deprecated ``Panel.to_long``, renamed to ``to_frame``
  - Deprecated ``colSpace`` argument in ``DataFrame.to_string``, renamed to
    ``col_space``
  - Rename ``precision`` to ``accuracy`` in engineering float formatter (GH
    GH395_)
  - The default delimiter for ``read_csv`` is comma rather than letting
    ``csv.Sniffer`` infer it
  - Rename ``col_or_columns`` argument in ``DataFrame.drop_duplicates`` (GH
    GH734_)

**Improvements to existing features**

  - Better error message in DataFrame constructor when passed column labels
    don't match data (GH497_)
  - Substantially improve performance of multi-GroupBy aggregation when a
    Python function is passed, reuse ndarray object in Cython (GH496_)
  - Can store objects indexed by tuples and floats in HDFStore (GH492_)
  - Don't print length by default in Series.to_string, add `length` option (GH
    GH489_)
  - Improve Cython code for multi-groupby to aggregate without having to sort
    the data (GH #93)
  - Improve MultiIndex reindexing speed by storing tuples in the MultiIndex,
    test for backwards unpickling compatibility
  - Improve column reindexing performance by using specialized Cython take
    function
  - Further performance tweaking of Series.__getitem__ for standard use cases
  - Avoid Index dict creation in some cases (i.e. when getting slices, etc.),
    regression from prior versions
  - Friendlier error message in setup.py if NumPy not installed
  - Use common set of NA-handling operations (sum, mean, etc.) in Panel class
    also (GH536_)
  - Default name assignment when calling ``reset_index`` on DataFrame with a
    regular (non-hierarchical) index (GH476_)
  - Use Cythonized groupers when possible in Series/DataFrame stat ops with
    ``level`` parameter passed (GH545_)
  - Ported skiplist data structure to C to speed up ``rolling_median`` by about
    5-10x in most typical use cases (GH374_)
  - Some performance enhancements in constructing a Panel from a dict of
    DataFrame objects
  - Made ``Index._get_duplicates`` a public method by removing the underscore
  - Prettier printing of floats, and column spacing fix (GH395_, GH571_)
  - Add ``bold_rows`` option to DataFrame.to_html (GH586_)
  - Improve the performance of ``DataFrame.sort_index`` by up to 5x or more
    when sorting by multiple columns
  - Substantially improve performance of DataFrame and Series constructors when
    passed a nested dict or dict, respectively (GH540_, GH621_)
  - Modified setup.py so that pip / setuptools will install dependencies (GH
    GH507_, various pull requests)
  - Unstack called on DataFrame with non-MultiIndex will return Series (GH
    GH477_)
  - Improve DataFrame.to_string and console formatting to be more consistent in
    the number of displayed digits (GH395_)
  - Use bottleneck if available for performing NaN-friendly statistical
    operations that it implemented (GH #91)
  - Monkey-patch context to traceback in ``DataFrame.apply`` to indicate which
    row/column the function application failed on (GH614_)
  - Improved ability of read_table and read_clipboard to parse
    console-formatted DataFrames (can read the row of index names, etc.)
  - Can pass list of group labels (without having to convert to an ndarray
    yourself) to ``groupby`` in some cases (GH659_)
  - Use ``kind`` argument to Series.order for selecting different sort kinds
    (GH668_)
  - Add option to Series.to_csv to omit the index (GH684_)
  - Add ``delimiter`` as an alternative to ``sep`` in ``read_csv`` and other
    parsing functions
  - Substantially improved performance of groupby on DataFrames with many
    columns by aggregating blocks of columns all at once (GH745_)
  - Can pass a file handle or StringIO to Series/DataFrame.to_csv (GH765_)
  - Can pass sequence of integers to DataFrame.irow(icol) and Series.iget, (GH
    GH654_)
  - Prototypes for some vectorized string functions
  - Add float64 hash table to solve the Series.unique problem with NAs (GH714_)
  - Memoize objects when reading from file to reduce memory footprint
  - Can get and set a column of a DataFrame with hierarchical columns
    containing "empty" ('') lower levels without passing the empty levels (PR
    GH768_)

**Bug fixes**

  - Raise exception in out-of-bounds indexing of Series instead of
    seg-faulting, regression from earlier releases (GH495_)
  - Fix error when joining DataFrames of different dtypes within the same
    typeclass (e.g. float32 and float64) (GH486_)
  - Fix bug in Series.min/Series.max on objects like datetime.datetime (GH
    GH487_)
  - Preserve index names in Index.union (GH501_)
  - Fix bug in Index joining causing subclass information (like DateRange type)
    to be lost in some cases (GH500_)
  - Accept empty list as input to DataFrame constructor, regression from 0.6.0
    (GH491_)
  - Can output DataFrame and Series with ndarray objects in a dtype=object
    array (GH490_)
  - Return empty string from Series.to_string when called on empty Series (GH
    GH488_)
  - Fix exception passing empty list to DataFrame.from_records
  - Fix Index.format bug (excluding name field) with datetimes with time info
  - Fix scalar value access in Series to always return NumPy scalars,
    regression from prior versions (GH510_)
  - Handle rows skipped at beginning of file in read_* functions (GH505_)
  - Handle improper dtype casting in ``set_value`` methods
  - Unary '-' / __neg__ operator on DataFrame was returning integer values
  - Unbox 0-dim ndarrays from certain operators like all, any in Series
  - Fix handling of missing columns (was combine_first-specific) in
    DataFrame.combine for general case (GH529_)
  - Fix type inference logic with boolean lists and arrays in DataFrame indexing
  - Use centered sum of squares in R-square computation if entity_effects=True
    in panel regression
  - Handle all NA case in Series.{corr, cov}, was raising exception (GH548_)
  - Aggregating by multiple levels with ``level`` argument to DataFrame, Series
    stat method, was broken (GH545_)
  - Fix Cython buf when converter passed to read_csv produced a numeric array
    (buffer dtype mismatch when passed to Cython type inference function) (GH
    GH546_)
  - Fix exception when setting scalar value using .ix on a DataFrame with a
    MultiIndex (GH551_)
  - Fix outer join between two DateRanges with different offsets that returned
    an invalid DateRange
  - Cleanup DataFrame.from_records failure where index argument is an integer
  - Fix Data.from_records failure when passed a dictionary
  - Fix NA handling in {Series, DataFrame}.rank with non-floating point dtypes
  - Fix bug related to integer type-checking in .ix-based indexing
  - Handle non-string index name passed to DataFrame.from_records
  - DataFrame.insert caused the columns name(s) field to be discarded (GH527_)
  - Fix erroneous in monotonic many-to-one left joins
  - Fix DataFrame.to_string to remove extra column white space (GH571_)
  - Format floats to default to same number of digits (GH395_)
  - Added decorator to copy docstring from one function to another (GH449_)
  - Fix error in monotonic many-to-one left joins
  - Fix __eq__ comparison between DateOffsets with different relativedelta
    keywords passed
  - Fix exception caused by parser converter returning strings (GH583_)
  - Fix MultiIndex formatting bug with integer names (GH601_)
  - Fix bug in handling of non-numeric aggregates in Series.groupby (GH612_)
  - Fix TypeError with tuple subclasses (e.g. namedtuple) in
    DataFrame.from_records (GH611_)
  - Catch misreported console size when running IPython within Emacs
  - Fix minor bug in pivot table margins, loss of index names and length-1
    'All' tuple in row labels
  - Add support for legacy WidePanel objects to be read from HDFStore
  - Fix out-of-bounds segfault in pad_object and backfill_object methods when
    either source or target array are empty
  - Could not create a new column in a DataFrame from a list of tuples
  - Fix bugs preventing SparseDataFrame and SparseSeries working with groupby
    (GH666_)
  - Use sort kind in Series.sort / argsort (GH668_)
  - Fix DataFrame operations on non-scalar, non-pandas objects (GH672_)
  - Don't convert DataFrame column to integer type when passing integer to
    __setitem__ (GH669_)
  - Fix downstream bug in pivot_table caused by integer level names in
    MultiIndex (GH678_)
  - Fix SparseSeries.combine_first when passed a dense Series (GH687_)
  - Fix performance regression in HDFStore loading when DataFrame or Panel
    stored in table format with datetimes
  - Raise Exception in DateRange when offset with n=0 is passed (GH683_)
  - Fix get/set inconsistency with .ix property and integer location but
    non-integer index (GH707_)
  - Use right dropna function for SparseSeries. Return dense Series for NA fill
    value (GH730_)
  - Fix Index.format bug causing incorrectly string-formatted Series with
    datetime indexes (# 726, 758)
  - Fix errors caused by object dtype arrays passed to ols (GH759_)
  - Fix error where column names lost when passing list of labels to
    DataFrame.__getitem__, (GH662_)
  - Fix error whereby top-level week iterator overwrote week instance
  - Fix circular reference causing memory leak in sparse array / series /
    frame, (GH663_)
  - Fix integer-slicing from integers-as-floats (GH670_)
  - Fix zero division errors in nanops from object dtype arrays in all NA case
    (GH676_)
  - Fix csv encoding when using unicode (GH705_, GH717_, GH738_)
  - Fix assumption that each object contains every unique block type in concat,
    (GH708_)
  - Fix sortedness check of multiindex in to_panel (GH719_, 720)
  - Fix that None was not treated as NA in PyObjectHashtable
  - Fix hashing dtype because of endianness confusion (GH747_, GH748_)
  - Fix SparseSeries.dropna to return dense Series in case of NA fill value (GH
    GH730_)
  - Use map_infer instead of np.vectorize. handle NA sentinels if converter
    yields numeric array, (GH753_)
  - Fixes and improvements to DataFrame.rank (GH742_)
  - Fix catching AttributeError instead of NameError for bottleneck
  - Try to cast non-MultiIndex to better dtype when calling reset_index (GH726_
    GH440_)
  - Fix #1.QNAN0' float bug on 2.6/win64
  - Allow subclasses of dicts in DataFrame constructor, with tests
  - Fix problem whereby set_index destroys column multiindex (GH764_)
  - Hack around bug in generating DateRange from naive DateOffset (GH770_)
  - Fix bug in DateRange.intersection causing incorrect results with some
    overlapping ranges (GH771_)

Thanks
------
- Craig Austin
- Chris Billington
- Marius Cobzarenco
- Mario Gamboa-Cavazos
- Hans-Martin Gaudecker
- Arthur Gerigk
- Yaroslav Halchenko
- Jeff Hammerbacher
- Matt Harrison
- Andreas Hilboll
- Luc Kesters
- Adam Klein
- Gregg Lind
- Solomon Negusse
- Wouter Overmeire
- Christian Prinoth
- Jeff Reback
- Sam Reckoner
- Craig Reeson
- Jan Schulz
- Skipper Seabold
- Ted Square
- Graham Taylor
- Aman Thakral
- Chris Uga
- Dieter Vandenbussche
- Texas P.
- Pinxing Ye
- ... and everyone I forgot

.. _GH220: https://github.com/pydata/pandas/issues/220
.. _GH249: https://github.com/pydata/pandas/issues/249
.. _GH267: https://github.com/pydata/pandas/issues/267
.. _GH468: https://github.com/pydata/pandas/issues/468
.. _GH479: https://github.com/pydata/pandas/issues/479
.. _GH273: https://github.com/pydata/pandas/issues/273
.. _GH498: https://github.com/pydata/pandas/issues/498
.. _GH526: https://github.com/pydata/pandas/issues/526
.. _GH534: https://github.com/pydata/pandas/issues/534
.. _GH521: https://github.com/pydata/pandas/issues/521
.. _GH542: https://github.com/pydata/pandas/issues/542
.. _GH552: https://github.com/pydata/pandas/issues/552
.. _GH554: https://github.com/pydata/pandas/issues/554
.. _GH570: https://github.com/pydata/pandas/issues/570
.. _GH571: https://github.com/pydata/pandas/issues/571
.. _GH115: https://github.com/pydata/pandas/issues/115
.. _GH595: https://github.com/pydata/pandas/issues/595
.. _GH563: https://github.com/pydata/pandas/issues/563
.. _GH338: https://github.com/pydata/pandas/issues/338
.. _GH614: https://github.com/pydata/pandas/issues/614
.. _GH464: https://github.com/pydata/pandas/issues/464
.. _GH371: https://github.com/pydata/pandas/issues/371
.. _GH629: https://github.com/pydata/pandas/issues/629
.. _GH170: https://github.com/pydata/pandas/issues/170
.. _GH166: https://github.com/pydata/pandas/issues/166
.. _GH628: https://github.com/pydata/pandas/issues/628
.. _GH658: https://github.com/pydata/pandas/issues/658
.. _GH647: https://github.com/pydata/pandas/issues/647
.. _GH288: https://github.com/pydata/pandas/issues/288
.. _GH699: https://github.com/pydata/pandas/issues/699
.. _GH657: https://github.com/pydata/pandas/issues/657
.. _GH740: https://github.com/pydata/pandas/issues/740
.. _GH731: https://github.com/pydata/pandas/issues/731
.. _GH737: https://github.com/pydata/pandas/issues/737
.. _GH700: https://github.com/pydata/pandas/issues/700
.. _GH328: https://github.com/pydata/pandas/issues/328
.. _GH316: https://github.com/pydata/pandas/issues/316
.. _GH395: https://github.com/pydata/pandas/issues/395
.. _GH734: https://github.com/pydata/pandas/issues/734
.. _GH497: https://github.com/pydata/pandas/issues/497
.. _GH496: https://github.com/pydata/pandas/issues/496
.. _GH492: https://github.com/pydata/pandas/issues/492
.. _GH489: https://github.com/pydata/pandas/issues/489
.. _GH536: https://github.com/pydata/pandas/issues/536
.. _GH476: https://github.com/pydata/pandas/issues/476
.. _GH545: https://github.com/pydata/pandas/issues/545
.. _GH374: https://github.com/pydata/pandas/issues/374
.. _GH586: https://github.com/pydata/pandas/issues/586
.. _GH540: https://github.com/pydata/pandas/issues/540
.. _GH621: https://github.com/pydata/pandas/issues/621
.. _GH507: https://github.com/pydata/pandas/issues/507
.. _GH477: https://github.com/pydata/pandas/issues/477
.. _GH659: https://github.com/pydata/pandas/issues/659
.. _GH668: https://github.com/pydata/pandas/issues/668
.. _GH684: https://github.com/pydata/pandas/issues/684
.. _GH745: https://github.com/pydata/pandas/issues/745
.. _GH765: https://github.com/pydata/pandas/issues/765
.. _GH654: https://github.com/pydata/pandas/issues/654
.. _GH714: https://github.com/pydata/pandas/issues/714
.. _GH768: https://github.com/pydata/pandas/issues/768
.. _GH495: https://github.com/pydata/pandas/issues/495
.. _GH486: https://github.com/pydata/pandas/issues/486
.. _GH487: https://github.com/pydata/pandas/issues/487
.. _GH501: https://github.com/pydata/pandas/issues/501
.. _GH500: https://github.com/pydata/pandas/issues/500
.. _GH491: https://github.com/pydata/pandas/issues/491
.. _GH490: https://github.com/pydata/pandas/issues/490
.. _GH488: https://github.com/pydata/pandas/issues/488
.. _GH510: https://github.com/pydata/pandas/issues/510
.. _GH505: https://github.com/pydata/pandas/issues/505
.. _GH529: https://github.com/pydata/pandas/issues/529
.. _GH548: https://github.com/pydata/pandas/issues/548
.. _GH546: https://github.com/pydata/pandas/issues/546
.. _GH551: https://github.com/pydata/pandas/issues/551
.. _GH527: https://github.com/pydata/pandas/issues/527
.. _GH449: https://github.com/pydata/pandas/issues/449
.. _GH583: https://github.com/pydata/pandas/issues/583
.. _GH601: https://github.com/pydata/pandas/issues/601
.. _GH612: https://github.com/pydata/pandas/issues/612
.. _GH611: https://github.com/pydata/pandas/issues/611
.. _GH666: https://github.com/pydata/pandas/issues/666
.. _GH672: https://github.com/pydata/pandas/issues/672
.. _GH669: https://github.com/pydata/pandas/issues/669
.. _GH678: https://github.com/pydata/pandas/issues/678
.. _GH687: https://github.com/pydata/pandas/issues/687
.. _GH683: https://github.com/pydata/pandas/issues/683
.. _GH707: https://github.com/pydata/pandas/issues/707
.. _GH730: https://github.com/pydata/pandas/issues/730
.. _GH759: https://github.com/pydata/pandas/issues/759
.. _GH662: https://github.com/pydata/pandas/issues/662
.. _GH663: https://github.com/pydata/pandas/issues/663
.. _GH670: https://github.com/pydata/pandas/issues/670
.. _GH676: https://github.com/pydata/pandas/issues/676
.. _GH705: https://github.com/pydata/pandas/issues/705
.. _GH717: https://github.com/pydata/pandas/issues/717
.. _GH738: https://github.com/pydata/pandas/issues/738
.. _GH708: https://github.com/pydata/pandas/issues/708
.. _GH719: https://github.com/pydata/pandas/issues/719
.. _GH747: https://github.com/pydata/pandas/issues/747
.. _GH748: https://github.com/pydata/pandas/issues/748
.. _GH753: https://github.com/pydata/pandas/issues/753
.. _GH742: https://github.com/pydata/pandas/issues/742
.. _GH726: https://github.com/pydata/pandas/issues/726
.. _GH440: https://github.com/pydata/pandas/issues/440
.. _GH764: https://github.com/pydata/pandas/issues/764
.. _GH770: https://github.com/pydata/pandas/issues/770
.. _GH771: https://github.com/pydata/pandas/issues/771


pandas 0.6.1
============

**Release date:** 12/13/2011

**API Changes**

  - Rename `names` argument in DataFrame.from_records to `columns`. Add
    deprecation warning
  - Boolean get/set operations on Series with boolean Series will reindex
    instead of requiring that the indexes be exactly equal (GH429_)

**New features / modules**

  - Can pass Series to DataFrame.append with ignore_index=True for appending a
    single row (GH430_)
  - Add Spearman and Kendall correlation options to Series.corr and
    DataFrame.corr (GH428_)
  - Add new `get_value` and `set_value` methods to Series, DataFrame, and Panel
    to very low-overhead access to scalar elements. df.get_value(row, column)
    is about 3x faster than df[column][row] by handling fewer cases (GH437_,
    GH438_). Add similar methods to sparse data structures for compatibility
  - Add Qt table widget to sandbox (GH435_)
  - DataFrame.align can accept Series arguments, add axis keyword (GH461_)
  - Implement new SparseList and SparseArray data structures. SparseSeries now
    derives from SparseArray (GH463_)
  - max_columns / max_rows options in set_printoptions (GH453_)
  - Implement Series.rank and DataFrame.rank, fast versions of
    scipy.stats.rankdata (GH428_)
  - Implement DataFrame.from_items alternate constructor (GH444_)
  - DataFrame.convert_objects method for inferring better dtypes for object
    columns (GH302_)
  - Add rolling_corr_pairwise function for computing Panel of correlation
    matrices (GH189_)
  - Add `margins` option to `pivot_table` for computing subgroup aggregates (GH
    GH114_)
  - Add `Series.from_csv` function (GH482_)

**Improvements to existing features**

  - Improve memory usage of `DataFrame.describe` (do not copy data
    unnecessarily) (GH425_)
  - Use same formatting function for outputting floating point Series to console
    as in DataFrame (GH420_)
  - DataFrame.delevel will try to infer better dtype for new columns (GH440_)
  - Exclude non-numeric types in DataFrame.{corr, cov}
  - Override Index.astype to enable dtype casting (GH412_)
  - Use same float formatting function for Series.__repr__ (GH420_)
  - Use available console width to output DataFrame columns (GH453_)
  - Accept ndarrays when setting items in Panel (GH452_)
  - Infer console width when printing __repr__ of DataFrame to console (PR
    GH453_)
  - Optimize scalar value lookups in the general case by 25% or more in Series
    and DataFrame
  - Can pass DataFrame/DataFrame and DataFrame/Series to
    rolling_corr/rolling_cov (GH462_)
  - Fix performance regression in cross-sectional count in DataFrame, affecting
    DataFrame.dropna speed
  - Column deletion in DataFrame copies no data (computes views on blocks) (GH
    GH158_)
  - MultiIndex.get_level_values can take the level name
  - More helpful error message when DataFrame.plot fails on one of the columns
    (GH478_)
  - Improve performance of DataFrame.{index, columns} attribute lookup

**Bug fixes**

  - Fix O(K^2) memory leak caused by inserting many columns without
    consolidating, had been present since 0.4.0 (GH467_)
  - `DataFrame.count` should return Series with zero instead of NA with length-0
    axis (GH423_)
  - Fix Yahoo! Finance API usage in pandas.io.data (GH419_, GH427_)
  - Fix upstream bug causing failure in Series.align with empty Series (GH434_)
  - Function passed to DataFrame.apply can return a list, as long as it's the
    right length. Regression from 0.4 (GH432_)
  - Don't "accidentally" upcast scalar values when indexing using .ix (GH431_)
  - Fix groupby exception raised with as_index=False and single column selected
    (GH421_)
  - Implement DateOffset.__ne__ causing downstream bug (GH456_)
  - Fix __doc__-related issue when converting py -> pyo with py2exe
  - Bug fix in left join Cython code with duplicate monotonic labels
  - Fix bug when unstacking multiple levels described in GH451_
  - Exclude NA values in dtype=object arrays, regression from 0.5.0 (GH469_)
  - Use Cython map_infer function in DataFrame.applymap to properly infer
    output type, handle tuple return values and other things that were breaking
    (GH465_)
  - Handle floating point index values in HDFStore (GH454_)
  - Fixed stale column reference bug (cached Series object) caused by type
    change / item deletion in DataFrame (GH473_)
  - Index.get_loc should always raise Exception when there are duplicates
  - Handle differently-indexed Series input to DataFrame constructor (GH475_)
  - Omit nuisance columns in multi-groupby with Python function
  - Buglet in handling of single grouping in general apply
  - Handle type inference properly when passing list of lists or tuples to
    DataFrame constructor (GH484_)
  - Preserve Index / MultiIndex names in GroupBy.apply concatenation step (GH
    GH481_)

Thanks
------
- Ralph Bean
- Luca Beltrame
- Marius Cobzarenco
- Andreas Hilboll
- Jev Kuznetsov
- Adam Lichtenstein
- Wouter Overmeire
- Fernando Perez
- Nathan Pinger
- Christian Prinoth
- Alex Reyfman
- Joon Ro
- Chang She
- Ted Square
- Chris Uga
- Dieter Vandenbussche

.. _GH429: https://github.com/pydata/pandas/issues/429
.. _GH430: https://github.com/pydata/pandas/issues/430
.. _GH428: https://github.com/pydata/pandas/issues/428
.. _GH437: https://github.com/pydata/pandas/issues/437
.. _GH438: https://github.com/pydata/pandas/issues/438
.. _GH435: https://github.com/pydata/pandas/issues/435
.. _GH461: https://github.com/pydata/pandas/issues/461
.. _GH463: https://github.com/pydata/pandas/issues/463
.. _GH453: https://github.com/pydata/pandas/issues/453
.. _GH444: https://github.com/pydata/pandas/issues/444
.. _GH302: https://github.com/pydata/pandas/issues/302
.. _GH189: https://github.com/pydata/pandas/issues/189
.. _GH114: https://github.com/pydata/pandas/issues/114
.. _GH482: https://github.com/pydata/pandas/issues/482
.. _GH425: https://github.com/pydata/pandas/issues/425
.. _GH420: https://github.com/pydata/pandas/issues/420
.. _GH440: https://github.com/pydata/pandas/issues/440
.. _GH412: https://github.com/pydata/pandas/issues/412
.. _GH452: https://github.com/pydata/pandas/issues/452
.. _GH462: https://github.com/pydata/pandas/issues/462
.. _GH158: https://github.com/pydata/pandas/issues/158
.. _GH478: https://github.com/pydata/pandas/issues/478
.. _GH467: https://github.com/pydata/pandas/issues/467
.. _GH423: https://github.com/pydata/pandas/issues/423
.. _GH419: https://github.com/pydata/pandas/issues/419
.. _GH427: https://github.com/pydata/pandas/issues/427
.. _GH434: https://github.com/pydata/pandas/issues/434
.. _GH432: https://github.com/pydata/pandas/issues/432
.. _GH431: https://github.com/pydata/pandas/issues/431
.. _GH421: https://github.com/pydata/pandas/issues/421
.. _GH456: https://github.com/pydata/pandas/issues/456
.. _GH451: https://github.com/pydata/pandas/issues/451
.. _GH469: https://github.com/pydata/pandas/issues/469
.. _GH465: https://github.com/pydata/pandas/issues/465
.. _GH454: https://github.com/pydata/pandas/issues/454
.. _GH473: https://github.com/pydata/pandas/issues/473
.. _GH475: https://github.com/pydata/pandas/issues/475
.. _GH484: https://github.com/pydata/pandas/issues/484
.. _GH481: https://github.com/pydata/pandas/issues/481


pandas 0.6.0
============

**Release date:** 11/25/2011

**API Changes**

  - Arithmetic methods like `sum` will attempt to sum dtype=object values by
    default instead of excluding them (GH382_)

**New features / modules**

  - Add `melt` function to `pandas.core.reshape`
  - Add `level` parameter to group by level in Series and DataFrame
    descriptive statistics (GH313_)
  - Add `head` and `tail` methods to Series, analogous to to DataFrame (PR
    GH296_)
  - Add `Series.isin` function which checks if each value is contained in a
    passed sequence (GH289_)
  - Add `float_format` option to `Series.to_string`
  - Add `skip_footer` (GH291_) and `converters` (GH343_) options to
    `read_csv` and `read_table`
  - Add proper, tested weighted least squares to standard and panel OLS (GH
    GH303_)
  - Add `drop_duplicates` and `duplicated` functions for removing duplicate
    DataFrame rows and checking for duplicate rows, respectively (GH319_)
  - Implement logical (boolean) operators &, |, ^ on DataFrame (GH347_)
  - Add `Series.mad`, mean absolute deviation, matching DataFrame
  - Add `QuarterEnd` DateOffset (GH321_)
  - Add matrix multiplication function `dot` to DataFrame (GH #65)
  - Add `orient` option to `Panel.from_dict` to ease creation of mixed-type
    Panels (GH359_, GH301_)
  - Add `DataFrame.from_dict` with similar `orient` option
  - Can now pass list of tuples or list of lists to `DataFrame.from_records`
    for fast conversion to DataFrame (GH357_)
  - Can pass multiple levels to groupby, e.g. `df.groupby(level=[0, 1])` (GH
    GH103_)
  - Can sort by multiple columns in `DataFrame.sort_index` (GH #92, GH362_)
  - Add fast `get_value` and `put_value` methods to DataFrame and
    micro-performance tweaks (GH360_)
  - Add `cov` instance methods to Series and DataFrame (GH194_, GH362_)
  - Add bar plot option to `DataFrame.plot` (GH348_)
  - Add `idxmin` and `idxmax` functions to Series and DataFrame for computing
    index labels achieving maximum and minimum values (GH286_)
  - Add `read_clipboard` function for parsing DataFrame from OS clipboard,
    should work across platforms (GH300_)
  - Add `nunique` function to Series for counting unique elements (GH297_)
  - DataFrame constructor will use Series name if no columns passed (GH373_)
  - Support regular expressions and longer delimiters in read_table/read_csv,
    but does not handle quoted strings yet (GH364_)
  - Add `DataFrame.to_html` for formatting DataFrame to HTML (GH387_)
  - MaskedArray can be passed to DataFrame constructor and masked values will be
    converted to NaN (GH396_)
  - Add `DataFrame.boxplot` function (GH368_, others)
  - Can pass extra args, kwds to DataFrame.apply (GH376_)

**Improvements to existing features**

  - Raise more helpful exception if date parsing fails in DateRange (GH298_)
  - Vastly improved performance of GroupBy on axes with a MultiIndex (GH299_)
  - Print level names in hierarchical index in Series repr (GH305_)
  - Return DataFrame when performing GroupBy on selected column and
    as_index=False (GH308_)
  - Can pass vector to `on` argument in `DataFrame.join` (GH312_)
  - Don't show Series name if it's None in the repr, also omit length for short
    Series (GH317_)
  - Show legend by default in `DataFrame.plot`, add `legend` boolean flag (GH
    GH324_)
  - Significantly improved performance of `Series.order`, which also makes
    np.unique called on a Series faster (GH327_)
  - Faster cythonized count by level in Series and DataFrame (GH341_)
  - Raise exception if dateutil 2.0 installed on Python 2.x runtime (GH346_)
  - Significant GroupBy performance enhancement with multiple keys with many
    "empty" combinations
  - New Cython vectorized function `map_infer` speeds up `Series.apply` and
    `Series.map` significantly when passed elementwise Python function,
    motivated by GH355_
  - Cythonized `cache_readonly`, resulting in substantial micro-performance
    enhancements throughout the codebase (GH361_)
  - Special Cython matrix iterator for applying arbitrary reduction operations
    with 3-5x better performance than `np.apply_along_axis` (GH309_)
  - Add `raw` option to `DataFrame.apply` for getting better performance when
    the passed function only requires an ndarray (GH309_)
  - Improve performance of `MultiIndex.from_tuples`
  - Can pass multiple levels to `stack` and `unstack` (GH370_)
  - Can pass multiple values columns to `pivot_table` (GH381_)
  - Can call `DataFrame.delevel` with standard Index with name set (GH393_)
  - Use Series name in GroupBy for result index (GH363_)
  - Refactor Series/DataFrame stat methods to use common set of NaN-friendly
    function
  - Handle NumPy scalar integers at C level in Cython conversion routines

**Bug fixes**

  - Fix bug in `DataFrame.to_csv` when writing a DataFrame with an index
    name (GH290_)
  - DataFrame should clear its Series caches on consolidation, was causing
    "stale" Series to be returned in some corner cases (GH304_)
  - DataFrame constructor failed if a column had a list of tuples (GH293_)
  - Ensure that `Series.apply` always returns a Series and implement
    `Series.round` (GH314_)
  - Support boolean columns in Cythonized groupby functions (GH315_)
  - `DataFrame.describe` should not fail if there are no numeric columns,
    instead return categorical describe (GH323_)
  - Fixed bug which could cause columns to be printed in wrong order in
    `DataFrame.to_string` if specific list of columns passed (GH325_)
  - Fix legend plotting failure if DataFrame columns are integers (GH326_)
  - Shift start date back by one month for Yahoo! Finance API in pandas.io.data
    (GH329_)
  - Fix `DataFrame.join` failure on unconsolidated inputs (GH331_)
  - DataFrame.min/max will no longer fail on mixed-type DataFrame (GH337_)
  - Fix `read_csv` / `read_table` failure when passing list to index_col that is
    not in ascending order (GH349_)
  - Fix failure passing Int64Index to Index.union when both are monotonic
  - Fix error when passing SparseSeries to (dense) DataFrame constructor
  - Added missing bang at top of setup.py (GH352_)
  - Change `is_monotonic` on MultiIndex so it properly compares the tuples
  - Fix MultiIndex outer join logic (GH351_)
  - Set index name attribute with single-key groupby (GH358_)
  - Bug fix in reflexive binary addition in Series and DataFrame for
    non-commutative operations (like string concatenation) (GH353_)
  - setupegg.py will invoke Cython (GH192_)
  - Fix block consolidation bug after inserting column into MultiIndex (GH366_)
  - Fix bug in join operations between Index and Int64Index (GH367_)
  - Handle min_periods=0 case in moving window functions (GH365_)
  - Fixed corner cases in DataFrame.apply/pivot with empty DataFrame (GH378_)
  - Fixed repr exception when Series name is a tuple
  - Always return DateRange from `asfreq` (GH390_)
  - Pass level names to `swaplavel` (GH379_)
  - Don't lose index names in `MultiIndex.droplevel` (GH394_)
  - Infer more proper return type in `DataFrame.apply` when no columns or rows
    depending on whether the passed function is a reduction (GH389_)
  - Always return NA/NaN from Series.min/max and DataFrame.min/max when all of a
    row/column/values are NA (GH384_)
  - Enable partial setting with .ix / advanced indexing (GH397_)
  - Handle mixed-type DataFrames correctly in unstack, do not lose type
    information (GH403_)
  - Fix integer name formatting bug in Index.format and in Series.__repr__
  - Handle label types other than string passed to groupby (GH405_)
  - Fix bug in .ix-based indexing with partial retrieval when a label is not
    contained in a level
  - Index name was not being pickled (GH408_)
  - Level name should be passed to result index in GroupBy.apply (GH416_)

Thanks
------

- Craig Austin
- Marius Cobzarenco
- Joel Cross
- Jeff Hammerbacher
- Adam Klein
- Thomas Kluyver
- Jev Kuznetsov
- Kieran O'Mahony
- Wouter Overmeire
- Nathan Pinger
- Christian Prinoth
- Skipper Seabold
- Chang She
- Ted Square
- Aman Thakral
- Chris Uga
- Dieter Vandenbussche
- carljv
- rsamson

.. _GH382: https://github.com/pydata/pandas/issues/382
.. _GH313: https://github.com/pydata/pandas/issues/313
.. _GH296: https://github.com/pydata/pandas/issues/296
.. _GH289: https://github.com/pydata/pandas/issues/289
.. _GH291: https://github.com/pydata/pandas/issues/291
.. _GH343: https://github.com/pydata/pandas/issues/343
.. _GH303: https://github.com/pydata/pandas/issues/303
.. _GH319: https://github.com/pydata/pandas/issues/319
.. _GH347: https://github.com/pydata/pandas/issues/347
.. _GH321: https://github.com/pydata/pandas/issues/321
.. _GH359: https://github.com/pydata/pandas/issues/359
.. _GH301: https://github.com/pydata/pandas/issues/301
.. _GH357: https://github.com/pydata/pandas/issues/357
.. _GH103: https://github.com/pydata/pandas/issues/103
.. _GH362: https://github.com/pydata/pandas/issues/362
.. _GH360: https://github.com/pydata/pandas/issues/360
.. _GH194: https://github.com/pydata/pandas/issues/194
.. _GH348: https://github.com/pydata/pandas/issues/348
.. _GH286: https://github.com/pydata/pandas/issues/286
.. _GH300: https://github.com/pydata/pandas/issues/300
.. _GH297: https://github.com/pydata/pandas/issues/297
.. _GH373: https://github.com/pydata/pandas/issues/373
.. _GH364: https://github.com/pydata/pandas/issues/364
.. _GH387: https://github.com/pydata/pandas/issues/387
.. _GH396: https://github.com/pydata/pandas/issues/396
.. _GH368: https://github.com/pydata/pandas/issues/368
.. _GH376: https://github.com/pydata/pandas/issues/376
.. _GH298: https://github.com/pydata/pandas/issues/298
.. _GH299: https://github.com/pydata/pandas/issues/299
.. _GH305: https://github.com/pydata/pandas/issues/305
.. _GH308: https://github.com/pydata/pandas/issues/308
.. _GH312: https://github.com/pydata/pandas/issues/312
.. _GH317: https://github.com/pydata/pandas/issues/317
.. _GH324: https://github.com/pydata/pandas/issues/324
.. _GH327: https://github.com/pydata/pandas/issues/327
.. _GH341: https://github.com/pydata/pandas/issues/341
.. _GH346: https://github.com/pydata/pandas/issues/346
.. _GH355: https://github.com/pydata/pandas/issues/355
.. _GH361: https://github.com/pydata/pandas/issues/361
.. _GH309: https://github.com/pydata/pandas/issues/309
.. _GH370: https://github.com/pydata/pandas/issues/370
.. _GH381: https://github.com/pydata/pandas/issues/381
.. _GH393: https://github.com/pydata/pandas/issues/393
.. _GH363: https://github.com/pydata/pandas/issues/363
.. _GH290: https://github.com/pydata/pandas/issues/290
.. _GH304: https://github.com/pydata/pandas/issues/304
.. _GH293: https://github.com/pydata/pandas/issues/293
.. _GH314: https://github.com/pydata/pandas/issues/314
.. _GH315: https://github.com/pydata/pandas/issues/315
.. _GH323: https://github.com/pydata/pandas/issues/323
.. _GH325: https://github.com/pydata/pandas/issues/325
.. _GH326: https://github.com/pydata/pandas/issues/326
.. _GH329: https://github.com/pydata/pandas/issues/329
.. _GH331: https://github.com/pydata/pandas/issues/331
.. _GH337: https://github.com/pydata/pandas/issues/337
.. _GH349: https://github.com/pydata/pandas/issues/349
.. _GH352: https://github.com/pydata/pandas/issues/352
.. _GH351: https://github.com/pydata/pandas/issues/351
.. _GH358: https://github.com/pydata/pandas/issues/358
.. _GH353: https://github.com/pydata/pandas/issues/353
.. _GH192: https://github.com/pydata/pandas/issues/192
.. _GH366: https://github.com/pydata/pandas/issues/366
.. _GH367: https://github.com/pydata/pandas/issues/367
.. _GH365: https://github.com/pydata/pandas/issues/365
.. _GH378: https://github.com/pydata/pandas/issues/378
.. _GH390: https://github.com/pydata/pandas/issues/390
.. _GH379: https://github.com/pydata/pandas/issues/379
.. _GH394: https://github.com/pydata/pandas/issues/394
.. _GH389: https://github.com/pydata/pandas/issues/389
.. _GH384: https://github.com/pydata/pandas/issues/384
.. _GH397: https://github.com/pydata/pandas/issues/397
.. _GH403: https://github.com/pydata/pandas/issues/403
.. _GH405: https://github.com/pydata/pandas/issues/405
.. _GH408: https://github.com/pydata/pandas/issues/408
.. _GH416: https://github.com/pydata/pandas/issues/416


pandas 0.5.0
============

**Release date:** 10/24/2011

This release of pandas includes a number of API changes (see below) and cleanup
of deprecated APIs from pre-0.4.0 releases. There are also bug fixes, new
features, numerous significant performance enhancements, and includes a new
IPython completer hook to enable tab completion of DataFrame columns accesses
as attributes (a new feature).

In addition to the changes listed here from 0.4.3 to 0.5.0, the minor releases
0.4.1, 0.4.2, and 0.4.3 brought some significant new functionality and
performance improvements that are worth taking a look at.

Thanks to all for bug reports, contributed patches and generally providing
feedback on the library.

**API Changes**

  - `read_table`, `read_csv`, and `ExcelFile.parse` default arguments for
    `index_col` is now None. To use one or more of the columns as the resulting
    DataFrame's index, these must be explicitly specified now
  - Parsing functions like `read_csv` no longer parse dates by default (GH
    GH225_)
  - Removed `weights` option in panel regression which was not doing anything
    principled (GH155_)
  - Changed `buffer` argument name in `Series.to_string` to `buf`
  - `Series.to_string` and `DataFrame.to_string` now return strings by default
    instead of printing to sys.stdout
  - Deprecated `nanRep` argument in various `to_string` and `to_csv` functions
    in favor of `na_rep`. Will be removed in 0.6 (GH275_)
  - Renamed `delimiter` to `sep` in `DataFrame.from_csv` for consistency
  - Changed order of `Series.clip` arguments to match those of `numpy.clip` and
    added (unimplemented) `out` argument so `numpy.clip` can be called on a
    Series (GH272_)
  - Series functions renamed (and thus deprecated) in 0.4 series have been
    removed:

    * `asOf`, use `asof`
    * `toDict`, use `to_dict`
    * `toString`, use `to_string`
    * `toCSV`, use `to_csv`
    * `merge`, use `map`
    * `applymap`, use `apply`
    * `combineFirst`, use `combine_first`
    * `_firstTimeWithValue` use `first_valid_index`
    * `_lastTimeWithValue` use `last_valid_index`

  - DataFrame functions renamed / deprecated in 0.4 series have been removed:

    * `asMatrix` method, use `as_matrix` or `values` attribute
    * `combineFirst`, use `combine_first`
    * `getXS`, use `xs`
    * `merge`, use `join`
    * `fromRecords`, use `from_records`
    * `fromcsv`, use `from_csv`
    * `toRecords`, use `to_records`
    * `toDict`, use `to_dict`
    * `toString`, use `to_string`
    * `toCSV`, use `to_csv`
    * `_firstTimeWithValue` use `first_valid_index`
    * `_lastTimeWithValue` use `last_valid_index`
    * `toDataMatrix` is no longer needed
    * `rows()` method, use `index` attribute
    * `cols()` method, use `columns` attribute
    * `dropEmptyRows()`, use `dropna(how='all')`
    * `dropIncompleteRows()`, use `dropna()`
    * `tapply(f)`, use `apply(f, axis=1)`
    * `tgroupby(keyfunc, aggfunc)`, use `groupby` with `axis=1`

  - Other outstanding deprecations have been removed:

    * `indexField` argument in `DataFrame.from_records`
    * `missingAtEnd` argument in `Series.order`. Use `na_last` instead
    * `Series.fromValue` classmethod, use regular `Series` constructor instead
    * Functions `parseCSV`, `parseText`, and `parseExcel` methods in
      `pandas.io.parsers` have been removed
    * `Index.asOfDate` function
    * `Panel.getMinorXS` (use `minor_xs`) and `Panel.getMajorXS` (use
      `major_xs`)
    * `Panel.toWide`, use `Panel.to_wide` instead

**New features / modules**

  - Added `DataFrame.align` method with standard join options
  - Added `parse_dates` option to `read_csv` and `read_table` methods to
    optionally try to parse dates in the index columns
  - Add `nrows`, `chunksize`, and `iterator` arguments to `read_csv` and
    `read_table`. The last two return a new `TextParser` class capable of
    lazily iterating through chunks of a flat file (GH242_)
  - Added ability to join on multiple columns in `DataFrame.join` (GH214_)
  - Added private `_get_duplicates` function to `Index` for identifying
    duplicate values more easily
  - Added column attribute access to DataFrame, e.g. df.A equivalent to df['A']
    if 'A' is a column in the DataFrame (GH213_)
  - Added IPython tab completion hook for DataFrame columns. (GH233_, GH230_)
  - Implement `Series.describe` for Series containing objects (GH241_)
  - Add inner join option to `DataFrame.join` when joining on key(s) (GH248_)
  - Can select set of DataFrame columns by passing a list to `__getitem__` (GH
    GH253_)
  - Can use & and | to intersection / union Index objects, respectively (GH
    GH261_)
  - Added `pivot_table` convenience function to pandas namespace (GH234_)
  - Implemented `Panel.rename_axis` function (GH243_)
  - DataFrame will show index level names in console output
  - Implemented `Panel.take`
  - Add `set_eng_float_format` function for setting alternate DataFrame
    floating point string formatting
  - Add convenience `set_index` function for creating a DataFrame index from
    its existing columns

**Improvements to existing features**

  - Major performance improvements in file parsing functions `read_csv` and
    `read_table`
  - Added Cython function for converting tuples to ndarray very fast. Speeds up
    many MultiIndex-related operations
  - File parsing functions like `read_csv` and `read_table` will explicitly
    check if a parsed index has duplicates and raise a more helpful exception
    rather than deferring the check until later
  - Refactored merging / joining code into a tidy class and disabled unnecessary
    computations in the float/object case, thus getting about 10% better
    performance (GH211_)
  - Improved speed of `DataFrame.xs` on mixed-type DataFrame objects by about
    5x, regression from 0.3.0 (GH215_)
  - With new `DataFrame.align` method, speeding up binary operations between
    differently-indexed DataFrame objects by 10-25%.
  - Significantly sped up conversion of nested dict into DataFrame (GH212_)
  - Can pass hierarchical index level name to `groupby` instead of the level
    number if desired (GH223_)
  - Add support for different delimiters in `DataFrame.to_csv` (GH244_)
  - Add more helpful error message when importing pandas post-installation from
    the source directory (GH250_)
  - Significantly speed up DataFrame `__repr__` and `count` on large mixed-type
    DataFrame objects
  - Better handling of pyx file dependencies in Cython module build (GH271_)

**Bug fixes**

  - `read_csv` / `read_table` fixes
    - Be less aggressive about converting float->int in cases of floating point
      representations of integers like 1.0, 2.0, etc.
    - "True"/"False" will not get correctly converted to boolean
    - Index name attribute will get set when specifying an index column
    - Passing column names should force `header=None` (GH257_)
    - Don't modify passed column names when `index_col` is not
      None (GH258_)
    - Can sniff CSV separator in zip file (since seek is not supported, was
      failing before)
  - Worked around matplotlib "bug" in which series[:, np.newaxis] fails. Should
    be reported upstream to matplotlib (GH224_)
  - DataFrame.iteritems was not returning Series with the name attribute
    set. Also neither was DataFrame._series
  - Can store datetime.date objects in HDFStore (GH231_)
  - Index and Series names are now stored in HDFStore
  - Fixed problem in which data would get upcasted to object dtype in
    GroupBy.apply operations (GH237_)
  - Fixed outer join bug with empty DataFrame (GH238_)
  - Can create empty Panel (GH239_)
  - Fix join on single key when passing list with 1 entry (GH246_)
  - Don't raise Exception on plotting DataFrame with an all-NA column (GH251_,
    GH254_)
  - Bug min/max errors when called on integer DataFrames (GH241_)
  - `DataFrame.iteritems` and `DataFrame._series` not assigning name attribute
  - Panel.__repr__ raised exception on length-0 major/minor axes
  - `DataFrame.join` on key with empty DataFrame produced incorrect columns
  - Implemented `MultiIndex.diff` (GH260_)
  - `Int64Index.take` and `MultiIndex.take` lost name field, fix downstream
    issue GH262_
  - Can pass list of tuples to `Series` (GH270_)
  - Can pass level name to `DataFrame.stack`
  - Support set operations between MultiIndex and Index
  - Fix many corner cases in MultiIndex set operations
    - Fix MultiIndex-handling bug with GroupBy.apply when returned groups are not
    indexed the same
  - Fix corner case bugs in DataFrame.apply
  - Setting DataFrame index did not cause Series cache to get cleared
  - Various int32 -> int64 platform-specific issues
  - Don't be too aggressive converting to integer when parsing file with
    MultiIndex (GH285_)
  - Fix bug when slicing Series with negative indices before beginning

Thanks
------

- Thomas Kluyver
- Daniel Fortunov
- Aman Thakral
- Luca Beltrame
- Wouter Overmeire

.. _GH225: https://github.com/pydata/pandas/issues/225
.. _GH155: https://github.com/pydata/pandas/issues/155
.. _GH275: https://github.com/pydata/pandas/issues/275
.. _GH272: https://github.com/pydata/pandas/issues/272
.. _GH242: https://github.com/pydata/pandas/issues/242
.. _GH214: https://github.com/pydata/pandas/issues/214
.. _GH213: https://github.com/pydata/pandas/issues/213
.. _GH233: https://github.com/pydata/pandas/issues/233
.. _GH230: https://github.com/pydata/pandas/issues/230
.. _GH241: https://github.com/pydata/pandas/issues/241
.. _GH248: https://github.com/pydata/pandas/issues/248
.. _GH253: https://github.com/pydata/pandas/issues/253
.. _GH261: https://github.com/pydata/pandas/issues/261
.. _GH234: https://github.com/pydata/pandas/issues/234
.. _GH243: https://github.com/pydata/pandas/issues/243
.. _GH211: https://github.com/pydata/pandas/issues/211
.. _GH215: https://github.com/pydata/pandas/issues/215
.. _GH212: https://github.com/pydata/pandas/issues/212
.. _GH223: https://github.com/pydata/pandas/issues/223
.. _GH244: https://github.com/pydata/pandas/issues/244
.. _GH250: https://github.com/pydata/pandas/issues/250
.. _GH271: https://github.com/pydata/pandas/issues/271
.. _GH257: https://github.com/pydata/pandas/issues/257
.. _GH258: https://github.com/pydata/pandas/issues/258
.. _GH224: https://github.com/pydata/pandas/issues/224
.. _GH231: https://github.com/pydata/pandas/issues/231
.. _GH237: https://github.com/pydata/pandas/issues/237
.. _GH238: https://github.com/pydata/pandas/issues/238
.. _GH239: https://github.com/pydata/pandas/issues/239
.. _GH246: https://github.com/pydata/pandas/issues/246
.. _GH251: https://github.com/pydata/pandas/issues/251
.. _GH254: https://github.com/pydata/pandas/issues/254
.. _GH260: https://github.com/pydata/pandas/issues/260
.. _GH262: https://github.com/pydata/pandas/issues/262
.. _GH270: https://github.com/pydata/pandas/issues/270
.. _GH285: https://github.com/pydata/pandas/issues/285


pandas 0.4.3
============

Release notes
-------------

**Release date:** 10/9/2011

This is largely a bugfix release from 0.4.2 but also includes a handful of new
and enhanced features. Also, pandas can now be installed and used on Python 3
(thanks Thomas Kluyver!).

**New features / modules**

  - Python 3 support using 2to3 (GH200_, Thomas Kluyver)
  - Add `name` attribute to `Series` and added relevant logic and tests. Name
    now prints as part of `Series.__repr__`
  - Add `name` attribute to standard Index so that stacking / unstacking does
    not discard names and so that indexed DataFrame objects can be reliably
    round-tripped to flat files, pickle, HDF5, etc.
  - Add `isnull` and `notnull` as instance methods on Series (GH209_, GH203_)

**Improvements to existing features**

  - Skip xlrd-related unit tests if not installed
  - `Index.append` and `MultiIndex.append` can accept a list of Index objects to
    concatenate together
  - Altered binary operations on differently-indexed SparseSeries objects to use
    the integer-based (dense) alignment logic which is faster with a larger
    number of blocks (GH205_)
  - Refactored `Series.__repr__` to be a bit more clean and consistent

**API Changes**

  - `Series.describe` and `DataFrame.describe` now bring the 25% and 75%
    quartiles instead of the 10% and 90% deciles. The other outputs have not
    changed
  - `Series.toString` will print deprecation warning, has been de-camelCased to
    `to_string`

**Bug fixes**

  - Fix broken interaction between `Index` and `Int64Index` when calling
    intersection. Implement `Int64Index.intersection`
  - `MultiIndex.sortlevel` discarded the level names (GH202_)
  - Fix bugs in groupby, join, and append due to improper concatenation of
    `MultiIndex` objects (GH201_)
  - Fix regression from 0.4.1, `isnull` and `notnull` ceased to work on other
    kinds of Python scalar objects like `datetime.datetime`
  - Raise more helpful exception when attempting to write empty DataFrame or
    LongPanel to `HDFStore` (GH204_)
  - Use stdlib csv module to properly escape strings with commas in
    `DataFrame.to_csv` (GH206_, Thomas Kluyver)
  - Fix Python ndarray access in Cython code for sparse blocked index integrity
    check
  - Fix bug writing Series to CSV in Python 3 (GH209_)
  - Miscellaneous Python 3 bugfixes

Thanks
------

  - Thomas Kluyver
  - rsamson

.. _GH200: https://github.com/pydata/pandas/issues/200
.. _GH209: https://github.com/pydata/pandas/issues/209
.. _GH203: https://github.com/pydata/pandas/issues/203
.. _GH205: https://github.com/pydata/pandas/issues/205
.. _GH202: https://github.com/pydata/pandas/issues/202
.. _GH201: https://github.com/pydata/pandas/issues/201
.. _GH204: https://github.com/pydata/pandas/issues/204
.. _GH206: https://github.com/pydata/pandas/issues/206


pandas 0.4.2
============

Release notes
-------------

**Release date:** 10/3/2011

This is a performance optimization release with several bug fixes. The new
Int64Index and new merging / joining Cython code and related Python
infrastructure are the main new additions

**New features / modules**

  - Added fast `Int64Index` type with specialized join, union,
    intersection. Will result in significant performance enhancements for
    int64-based time series (e.g. using NumPy's datetime64 one day) and also
    faster operations on DataFrame objects storing record array-like data.
  - Refactored `Index` classes to have a `join` method and associated data
    alignment routines throughout the codebase to be able to leverage optimized
    joining / merging routines.
  - Added `Series.align` method for aligning two series with choice of join
    method
  - Wrote faster Cython data alignment / merging routines resulting in
    substantial speed increases
  - Added `is_monotonic` property to `Index` classes with associated Cython
    code to evaluate the monotonicity of the `Index` values
  - Add method `get_level_values` to `MultiIndex`
  - Implemented shallow copy of `BlockManager` object in `DataFrame` internals

**Improvements to existing features**

  - Improved performance of `isnull` and `notnull`, a regression from v0.3.0
    (GH187_)
  - Wrote templating / code generation script to auto-generate Cython code for
    various functions which need to be available for the 4 major data types
    used in pandas (float64, bool, object, int64)
  - Refactored code related to `DataFrame.join` so that intermediate aligned
    copies of the data in each `DataFrame` argument do not need to be
    created. Substantial performance increases result (GH176_)
  - Substantially improved performance of generic `Index.intersection` and
    `Index.union`
  - Improved performance of `DateRange.union` with overlapping ranges and
    non-cacheable offsets (like Minute). Implemented analogous fast
    `DateRange.intersection` for overlapping ranges.
  - Implemented `BlockManager.take` resulting in significantly faster `take`
    performance on mixed-type `DataFrame` objects (GH104_)
  - Improved performance of `Series.sort_index`
  - Significant groupby performance enhancement: removed unnecessary integrity
    checks in DataFrame internals that were slowing down slicing operations to
    retrieve groups
  - Added informative Exception when passing dict to DataFrame groupby
    aggregation with axis != 0

**API Changes**

None

**Bug fixes**

  - Fixed minor unhandled exception in Cython code implementing fast groupby
    aggregation operations
  - Fixed bug in unstacking code manifesting with more than 3 hierarchical
    levels
  - Throw exception when step specified in label-based slice (GH185_)
  - Fix isnull to correctly work with np.float32. Fix upstream bug described in
    GH182_
  - Finish implementation of as_index=False in groupby for DataFrame
    aggregation (GH181_)
  - Raise SkipTest for pre-epoch HDFStore failure. Real fix will be sorted out
    via datetime64 dtype

Thanks
------

- Uri Laserson
- Scott Sinclair

.. _GH187: https://github.com/pydata/pandas/issues/187
.. _GH176: https://github.com/pydata/pandas/issues/176
.. _GH104: https://github.com/pydata/pandas/issues/104
.. _GH185: https://github.com/pydata/pandas/issues/185
.. _GH182: https://github.com/pydata/pandas/issues/182
.. _GH181: https://github.com/pydata/pandas/issues/181


pandas 0.4.1
============

Release notes
-------------

**Release date:** 9/25/2011

This is primarily a bug fix release but includes some new features and
improvements

**New features / modules**

  - Added new `DataFrame` methods `get_dtype_counts` and property `dtypes`
  - Setting of values using ``.ix`` indexing attribute in mixed-type DataFrame
    objects has been implemented (fixes GH135_)
  - `read_csv` can read multiple columns into a `MultiIndex`. DataFrame's
    `to_csv` method will properly write out a `MultiIndex` which can be read
    back (GH151_, thanks to Skipper Seabold)
  - Wrote fast time series merging / joining methods in Cython. Will be
    integrated later into DataFrame.join and related functions
  - Added `ignore_index` option to `DataFrame.append` for combining unindexed
    records stored in a DataFrame

**Improvements to existing features**

  - Some speed enhancements with internal Index type-checking function
  - `DataFrame.rename` has a new `copy` parameter which can rename a DataFrame
    in place
  - Enable unstacking by level name (GH142_)
  - Enable sortlevel to work by level name (GH141_)
  - `read_csv` can automatically "sniff" other kinds of delimiters using
    `csv.Sniffer` (GH146_)
  - Improved speed of unit test suite by about 40%
  - Exception will not be raised calling `HDFStore.remove` on non-existent node
    with where clause
  - Optimized `_ensure_index` function resulting in performance savings in
    type-checking Index objects

**API Changes**

None

**Bug fixes**

  - Fixed DataFrame constructor bug causing downstream problems (e.g. .copy()
    failing) when passing a Series as the values along with a column name and
    index
  - Fixed single-key groupby on DataFrame with as_index=False (GH160_)
  - `Series.shift` was failing on integer Series (GH154_)
  - `unstack` methods were producing incorrect output in the case of duplicate
    hierarchical labels. An exception will now be raised (GH147_)
  - Calling `count` with level argument caused reduceat failure or segfault in
    earlier NumPy (GH169_)
  - Fixed `DataFrame.corrwith` to automatically exclude non-numeric data (GH
    GH144_)
  - Unicode handling bug fixes in `DataFrame.to_string` (GH138_)
  - Excluding OLS degenerate unit test case that was causing platform specific
    failure (GH149_)
  - Skip blosc-dependent unit tests for PyTables < 2.2 (GH137_)
  - Calling `copy` on `DateRange` did not copy over attributes to the new object
    (GH168_)
  - Fix bug in `HDFStore` in which Panel data could be appended to a Table with
    different item order, thus resulting in an incorrect result read back

Thanks
------
- Yaroslav Halchenko
- Jeff Reback
- Skipper Seabold
- Dan Lovell
- Nick Pentreath

.. _GH135: https://github.com/pydata/pandas/issues/135
.. _GH151: https://github.com/pydata/pandas/issues/151
.. _GH142: https://github.com/pydata/pandas/issues/142
.. _GH141: https://github.com/pydata/pandas/issues/141
.. _GH146: https://github.com/pydata/pandas/issues/146
.. _GH160: https://github.com/pydata/pandas/issues/160
.. _GH154: https://github.com/pydata/pandas/issues/154
.. _GH147: https://github.com/pydata/pandas/issues/147
.. _GH169: https://github.com/pydata/pandas/issues/169
.. _GH144: https://github.com/pydata/pandas/issues/144
.. _GH138: https://github.com/pydata/pandas/issues/138
.. _GH149: https://github.com/pydata/pandas/issues/149
.. _GH137: https://github.com/pydata/pandas/issues/137
.. _GH168: https://github.com/pydata/pandas/issues/168


pandas 0.4.0
============

Release notes
-------------

**Release date:** 9/12/2011

**New features / modules**

  - `pandas.core.sparse` module: "Sparse" (mostly-NA, or some other fill value)
    versions of `Series`, `DataFrame`, and `Panel`. For low-density data, this
    will result in significant performance boosts, and smaller memory
    footprint. Added `to_sparse` methods to `Series`, `DataFrame`, and
    `Panel`. See online documentation for more on these
  - Fancy indexing operator on Series / DataFrame, e.g. via .ix operator. Both
    getting and setting of values is supported; however, setting values will only
    currently work on homogeneously-typed DataFrame objects. Things like:

    * series.ix[[d1, d2, d3]]
    * frame.ix[5:10, ['C', 'B', 'A']], frame.ix[5:10, 'A':'C']
    * frame.ix[date1:date2]

  - Significantly enhanced `groupby` functionality

    * Can groupby multiple keys, e.g. df.groupby(['key1', 'key2']). Iteration with
      multiple groupings products a flattened tuple
    * "Nuisance" columns (non-aggregatable) will automatically be excluded from
      DataFrame aggregation operations
    * Added automatic "dispatching to Series / DataFrame methods to more easily
      invoke methods on groups. e.g. s.groupby(crit).std() will work even though
      `std` is not implemented on the `GroupBy` class

  - Hierarchical / multi-level indexing

    * New the `MultiIndex` class. Integrated `MultiIndex` into `Series` and
      `DataFrame` fancy indexing, slicing, __getitem__ and __setitem,
      reindexing, etc. Added `level` keyword argument to `groupby` to enable
      grouping by a level of a `MultiIndex`

  - New data reshaping functions: `stack` and `unstack` on DataFrame and Series

    * Integrate with MultiIndex to enable sophisticated reshaping of data

  - `Index` objects (labels for axes) are now capable of holding tuples
  - `Series.describe`, `DataFrame.describe`: produces an R-like table of summary
    statistics about each data column
  - `DataFrame.quantile`, `Series.quantile` for computing sample quantiles of data
    across requested axis
  - Added general `DataFrame.dropna` method to replace `dropIncompleteRows` and
    `dropEmptyRows`, deprecated those.
  - `Series` arithmetic methods with optional fill_value for missing data,
    e.g. a.add(b, fill_value=0). If a location is missing for both it will still
    be missing in the result though.
  - fill_value option has been added to `DataFrame`.{add, mul, sub, div} methods
    similar to `Series`
  - Boolean indexing with `DataFrame` objects: data[data > 0.1] = 0.1 or
    data[data> other] = 1.
  - `pytz` / tzinfo support in `DateRange`

    * `tz_localize`, `tz_normalize`, and `tz_validate` methods added

  - Added `ExcelFile` class to `pandas.io.parsers` for parsing multiple sheets out
    of a single Excel 2003 document
  - `GroupBy` aggregations can now optionally *broadcast*, e.g. produce an object
    of the same size with the aggregated value propagated
  - Added `select` function in all data structures: reindex axis based on
    arbitrary criterion (function returning boolean value),
    e.g. frame.select(lambda x: 'foo' in x, axis=1)
  - `DataFrame.consolidate` method, API function relating to redesigned internals
  - `DataFrame.insert` method for inserting column at a specified location rather
    than the default __setitem__ behavior (which puts it at the end)
  - `HDFStore` class in `pandas.io.pytables` has been largely rewritten using
    patches from Jeff Reback from others. It now supports mixed-type `DataFrame`
    and `Series` data and can store `Panel` objects. It also has the option to
    query `DataFrame` and `Panel` data. Loading data from legacy `HDFStore`
    files is supported explicitly in the code
  - Added `set_printoptions` method to modify appearance of DataFrame tabular
    output
  - `rolling_quantile` functions; a moving version of `Series.quantile` /
    `DataFrame.quantile`
  - Generic `rolling_apply` moving window function
  - New `drop` method added to `Series`, `DataFrame`, etc. which can drop a set of
    labels from an axis, producing a new object
  - `reindex` methods now sport a `copy` option so that data is not forced to be
    copied then the resulting object is indexed the same
  - Added `sort_index` methods to Series and Panel. Renamed `DataFrame.sort`
    to `sort_index`. Leaving `DataFrame.sort` for now.
  - Added ``skipna`` option to statistical instance methods on all the data
    structures
  - `pandas.io.data` module providing a consistent interface for reading time
    series data from several different sources

**Improvements to existing features**

  * The 2-dimensional `DataFrame` and `DataMatrix` classes have been extensively
    redesigned internally into a single class `DataFrame`, preserving where
    possible their optimal performance characteristics. This should reduce
    confusion from users about which class to use.

    * Note that under the hood there is a new essentially "lazy evaluation"
      scheme within respect to adding columns to DataFrame. During some
      operations, like-typed blocks will be "consolidated" but not before.

  * `DataFrame` accessing columns repeatedly is now significantly faster than
    `DataMatrix` used to be in 0.3.0 due to an internal Series caching mechanism
    (which are all views on the underlying data)
  * Column ordering for mixed type data is now completely consistent in
    `DataFrame`. In prior releases, there was inconsistent column ordering in
    `DataMatrix`
  * Improved console / string formatting of DataMatrix with negative numbers
  * Improved tabular data parsing functions, `read_table` and `read_csv`:

    * Added `skiprows` and `na_values` arguments to `pandas.io.parsers` functions
      for more flexible IO
    * `parseCSV` / `read_csv` functions and others in `pandas.io.parsers` now can
      take a list of custom NA values, and also a list of rows to skip

  * Can slice `DataFrame` and get a view of the data (when homogeneously typed),
    e.g. frame.xs(idx, copy=False) or frame.ix[idx]
  * Many speed optimizations throughout `Series` and `DataFrame`
  * Eager evaluation of groups when calling ``groupby`` functions, so if there is
    an exception with the grouping function it will raised immediately versus
    sometime later on when the groups are needed
  * `datetools.WeekOfMonth` offset can be parameterized with `n` different than 1
    or -1.
  * Statistical methods on DataFrame like `mean`, `std`, `var`, `skew` will now
    ignore non-numerical data. Before a not very useful error message was
    generated. A flag `numeric_only` has been added to `DataFrame.sum` and
    `DataFrame.count` to enable this behavior in those methods if so desired
    (disabled by default)
  * `DataFrame.pivot` generalized to enable pivoting multiple columns into a
    `DataFrame` with hierarchical columns
  * `DataFrame` constructor can accept structured / record arrays
  * `Panel` constructor can accept a dict of DataFrame-like objects. Do not
    need to use `from_dict` anymore (`from_dict` is there to stay, though).

**API Changes**

  * The `DataMatrix` variable now refers to `DataFrame`, will be removed within
    two releases
  * `WidePanel` is now known as `Panel`. The `WidePanel` variable in the pandas
    namespace now refers to the renamed `Panel` class
  * `LongPanel` and `Panel` / `WidePanel` now no longer have a common
    subclass. `LongPanel` is now a subclass of `DataFrame` having a number of
    additional methods and a hierarchical index instead of the old
    `LongPanelIndex` object, which has been removed. Legacy `LongPanel` pickles
    may not load properly
  * Cython is now required to build `pandas` from a development branch. This was
    done to avoid continuing to check in cythonized C files into source
    control. Builds from released source distributions will not require Cython
  * Cython code has been moved up to a top level `pandas/src` directory. Cython
    extension modules have been renamed and promoted from the `lib` subpackage to
    the top level, i.e.

    * `pandas.lib.tseries` -> `pandas._tseries`
    * `pandas.lib.sparse` -> `pandas._sparse`

  * `DataFrame` pickling format has changed. Backwards compatibility for legacy
    pickles is provided, but it's recommended to consider PyTables-based
    `HDFStore` for storing data with a longer expected shelf life
  * A `copy` argument has been added to the `DataFrame` constructor to avoid
    unnecessary copying of data. Data is no longer copied by default when passed
    into the constructor
  * Handling of boolean dtype in `DataFrame` has been improved to support storage
    of boolean data with NA / NaN values. Before it was being converted to float64
    so this should not (in theory) cause API breakage
  * To optimize performance, Index objects now only check that their labels are
    unique when uniqueness matters (i.e. when someone goes to perform a
    lookup). This is a potentially dangerous tradeoff, but will lead to much
    better performance in many places (like groupby).
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
    `fromRecords`, `tgroupby`, `toString`
  * `pandas.io.parsers` method deprecations

    * `parseCSV` is now `read_csv` and keyword arguments have been de-camelCased
    * `parseText` is now `read_table`
    * `parseExcel` is replaced by the `ExcelFile` class and its `parse` method

  * `fillMethod` arguments (deprecated in prior release) removed, should be
    replaced with `method`
  * `Series.fill`, `DataFrame.fill`, and `Panel.fill` removed, use `fillna`
    instead
  * `groupby` functions now exclude NA / NaN values from the list of groups. This
    matches R behavior with NAs in factors e.g. with the `tapply` function
  * Removed `parseText`, `parseCSV` and `parseExcel` from pandas namespace
  * `Series.combineFunc` renamed to `Series.combine` and made a bit more general
    with a `fill_value` keyword argument defaulting to NaN
  * Removed `pandas.core.pytools` module. Code has been moved to
    `pandas.core.common`
  * Tacked on `groupName` attribute for groups in GroupBy renamed to `name`
  * Panel/LongPanel `dims` attribute renamed to `shape` to be more conformant
  * Slicing a `Series` returns a view now
  * More Series deprecations / renaming: `toCSV` to `to_csv`, `asOf` to `asof`,
    `merge` to `map`, `applymap` to `apply`, `toDict` to `to_dict`,
    `combineFirst` to `combine_first`. Will print `FutureWarning`.
  * `DataFrame.to_csv` does not write an "index" column label by default
    anymore since the output file can be read back without it. However, there
    is a new ``index_label`` argument. So you can do ``index_label='index'`` to
    emulate the old behavior
  * `datetools.Week` argument renamed from `dayOfWeek` to `weekday`
  * `timeRule` argument in `shift` has been deprecated in favor of using the
    `offset` argument for everything. So you can still pass a time rule string
    to `offset`
  * Added optional `encoding` argument to `read_csv`, `read_table`, `to_csv`,
    `from_csv` to handle unicode in python 2.x

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
  * Added `__ne__` implementation to `DataFrame`
  * Bug / unintuitive result when calling `fillna` on unordered labels
  * Bug calling `sum` on boolean DataFrame
  * Bug fix when creating a DataFrame from a dict with scalar values
  * Series.{sum, mean, std, ...} now return NA/NaN when the whole Series is NA
  * NumPy 1.4 through 1.6 compatibility fixes
  * Fixed bug in bias correction in `rolling_cov`, was affecting `rolling_corr`
    too
  * R-square value was incorrect in the presence of fixed and time effects in
    the `PanelOLS` classes
  * `HDFStore` can handle duplicates in table format, will take

Thanks
------
  - Joon Ro
  - Michael Pennington
  - Chris Uga
  - Chris Withers
  - Jeff Reback
  - Ted Square
  - Craig Austin
  - William Ferreira
  - Daniel Fortunov
  - Tony Roberts
  - Martin Felder
  - John Marino
  - Tim McNamara
  - Justin Berka
  - Dieter Vandenbussche
  - Shane Conway
  - Skipper Seabold
  - Chris Jordan-Squire

pandas 0.3.0
============

Release notes
-------------

**Release date:** February 20, 2011

**New features / modules**

  - `corrwith` function to compute column- or row-wise correlations between two
	DataFrame objects
  - Can boolean-index DataFrame objects, e.g. df[df > 2] = 2, px[px > last_px] = 0
  - Added comparison magic methods (__lt__, __gt__, etc.)
  - Flexible explicit arithmetic methods (add, mul, sub, div, etc.)
  - Added `reindex_like` method
  - Added `reindex_like` method to WidePanel
  - Convenience functions for accessing SQL-like databases in `pandas.io.sql`
	module
  - Added (still experimental) HDFStore class for storing pandas data
	structures using HDF5 / PyTables in `pandas.io.pytables` module
  - Added WeekOfMonth date offset
  - `pandas.rpy` (experimental) module created, provide some interfacing /
   conversion between rpy2 and pandas

**Improvements**

  - Unit test coverage: 100% line coverage of core data structures
  - Speed enhancement to rolling_{median, max, min}
  - Column ordering between DataFrame and DataMatrix is now consistent: before
	DataFrame would not respect column order
  - Improved {Series, DataFrame}.plot methods to be more flexible (can pass
	matplotlib Axis arguments, plot DataFrame columns in multiple subplots,
	etc.)

**API Changes**

  - Exponentially-weighted moment functions in `pandas.stats.moments` have a
	more consistent API and accept a min_periods argument like their regular
	moving counterparts.
  - **fillMethod** argument in Series, DataFrame changed to **method**,
	`FutureWarning` added.
  - **fill** method in Series, DataFrame/DataMatrix, WidePanel renamed to
	**fillna**, `FutureWarning` added to **fill**
  - Renamed **DataFrame.getXS** to **xs**, `FutureWarning` added
  - Removed **cap** and **floor** functions from DataFrame, renamed to
	**clip_upper** and **clip_lower** for consistency with NumPy

**Bug fixes**

  - Fixed bug in IndexableSkiplist Cython code that was breaking
	rolling_max function
  - Numerous numpy.int64-related indexing fixes
  - Several NumPy 1.4.0 NaN-handling fixes
  - Bug fixes to pandas.io.parsers.parseCSV
  - Fixed `DateRange` caching issue with unusual date offsets
  - Fixed bug in `DateRange.union`
  - Fixed corner case in `IndexableSkiplist` implementation
