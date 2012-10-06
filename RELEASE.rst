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

pandas 0.9.0
============

**Release date:** NOT YET RELEASED

**New features**

  - Add ``str.encode`` and ``str.decode`` to Series (#1706)
  - Add `to_latex` method to DataFrame (#1735)
  - Add convenient expanding window equivalents of all rolling_* ops (#1785)
  - Add Options class to pandas.io.data for fetching options data from Yahoo!
    Finance (#1748, #1739)
  - Recognize and convert more boolean values in file parsing (Yes, No, TRUE,
    FALSE, variants thereof) (#1691, #1295)
  - Add Panel.update method, analogous to DataFrame.update (#1999, #1988)

**Improvements to existing features**

  - Proper handling of NA values in merge operations (#1990)
  - Add ``flags`` option for ``re.compile`` in some Series.str methods (#1659)
  - Parsing of UTC date strings in read_* functions (#1693)
  - Handle generator input to Series (#1679)
  - Add `na_action='ignore'` to Series.map to quietly propagate NAs (#1661)
  - Add args/kwds options to Series.apply (#1829)
  - Add inplace option to Series/DataFrame.reset_index (#1797)
  - Add ``level`` parameter to ``Series.reset_index``
  - Add quoting option for DataFrame.to_csv (#1902)
  - Indicate long column value truncation in DataFrame output with ... (#1854)
  - DataFrame.dot will not do data alignment, and also work with Series (#1915)
  - Add ``na`` option for missing data handling in some vectorized string
    methods (#1689)
  - If index_label=False in DataFrame.to_csv, do not print fields/commas in the
    text output. Results in easier importing into R (#1583)
  - Can pass tuple/list of axes to DataFrame.dropna to simplify repeated calls
    (dropping both columns and rows) (#924)
  - Improve DataFrame.to_html output for hierarchically-indexed rows (do not
    repeat levels) (#1929)
  - TimeSeries.between_time can now select times across midnight (#1871)
  - Enable `skip_footer` parameter in `ExcelFile.parse` (#1843)

**API Changes**

  - Change default header names in read_* functions to more Pythonic X0, X1,
    etc. instead of X.1, X.2. (#2000)
  - Deprecated ``day_of_year`` API removed from PeriodIndex, use ``dayofyear``
    (#1723)
  - Don't modify NumPy suppress printoption at import time
  - The internal HDF5 data arrangement for DataFrames has been
    transposed. Legacy files will still be readable by HDFStore (#1834, #1824)
  - Legacy cruft removed: pandas.stats.misc.quantileTS
  - Use ISO8601 format for Period repr: monthly, daily, and on down (#1776)
  - Empty DataFrame columns are now created as object dtype. This will prevent
    a class of TypeErrors that was occurring in code where the dtype of a
    column would depend on the presence of data or not (e.g. a SQL query having
    results) (#1783)
  - Setting parts of DataFrame/Panel using ix now aligns input Series/DataFrame
    (#1630)
  - `first` and `last` methods in `GroupBy` no longer drop non-numeric columns
    (#1809)
  - Resolved inconsistencies in specifying custom NA values in text parser.
    `na_values` of type dict no longer override default NAs unless
    `keep_default_na` is set to false explicitly (#1657)
  - Enable `skipfooter` parameter in text parsers as an alias for `skip_footer`

**Bug fixes**

  - Perform arithmetic column-by-column in mixed-type DataFrame to avoid type
    upcasting issues. Caused downstream DataFrame.diff bug (#1896)
  - Fix matplotlib auto-color assignment when no custom spectrum passed. Also
    respect passed color keyword argument (#1711)
  - Fix resampling logical error with closed='left' (#1726)
  - Fix critical DatetimeIndex.union bugs (#1730, #1719, #1745, #1702, #1753)
  - Fix critical DatetimeIndex.intersection bug with unanchored offsets (#1708)
  - Fix MM-YYYY time series indexing case (#1672)
  - Fix case where Categorical group key was not being passed into index in
    GroupBy result (#1701)
  - Handle Ellipsis in Series.__getitem__/__setitem__ (#1721)
  - Fix some bugs with handling datetime64 scalars of other units in NumPy 1.6
    and 1.7 (#1717)
  - Fix performance issue in MultiIndex.format (#1746)
  - Fixed GroupBy bugs interacting with DatetimeIndex asof / map methods (#1677)
  - Handle factors with NAs in pandas.rpy (#1615)
  - Fix statsmodels import in pandas.stats.var (#1734)
  - Fix DataFrame repr/info summary with non-unique columns (#1700)
  - Fix Series.iget_value for non-unique indexes (#1694)
  - Don't lose tzinfo when passing DatetimeIndex as DataFrame column (#1682)
  - Fix tz conversion with time zones that haven't had any DST transitions since
    first date in the array (#1673)
  - Fix field access with  UTC->local conversion on unsorted arrays (#1756)
  - Fix isnull handling of array-like (list) inputs (#1755)
  - Fix regression in handling of Series in Series constructor (#1671)
  - Fix comparison of Int64Index with DatetimeIndex (#1681)
  - Fix min_periods handling in new rolling_max/min at array start (#1695)
  - Fix errors with how='median' and generic NumPy resampling in some cases
    caused by SeriesBinGrouper (#1648, #1688)
  - When grouping by level, exclude unobserved levels (#1697)
  - Don't lose tzinfo in DatetimeIndex when shifting by different offset (#1683)
  - Hack to support storing data with a zero-length axis in HDFStore (#1707)
  - Fix DatetimeIndex tz-aware range generation issue (#1674)
  - Fix method='time' interpolation with intraday data (#1698)
  - Don't plot all-NA DataFrame columns as zeros (#1696)
  - Fix bug in scatter_plot with by option (#1716)
  - Fix performance problem in infer_freq with lots of non-unique stamps (#1686)
  - Fix handling of PeriodIndex as argument to create MultiIndex (#1705)
  - Fix re: unicode MultiIndex level names in Series/DataFrame repr (#1736)
  - Handle PeriodIndex in to_datetime instance method (#1703)
  - Support StaticTzInfo in DatetimeIndex infrastructure (#1692)
  - Allow MultiIndex setops with length-0 other type indexes (#1727)
  - Fix handling of DatetimeIndex in DataFrame.to_records (#1720)
  - Fix handling of general objects in isnull on which bool(...) fails (#1749)
  - Fix .ix indexing with MultiIndex ambiguity (#1678)
  - Fix .ix setting logic error with non-unique MultiIndex (#1750)
  - Basic indexing now works on MultiIndex with > 1000000 elements, regression
    from earlier version of pandas (#1757)
  - Handle non-float64 dtypes in fast DataFrame.corr/cov code paths (#1761)
  - Fix DatetimeIndex.isin to function properly (#1763)
  - Fix conversion of array of tz-aware datetime.datetime to DatetimeIndex with
    right time zone (#1777)
  - Fix DST issues with generating ancxhored date ranges (#1778)
  - Fix issue calling sort on result of Series.unique (#1807)
  - Fix numerical issue leading to square root of negative number in
    rolling_std (#1840)
  - Let Series.str.split accept no arguments (like str.split) (#1859)
  - Allow user to have dateutil 2.1 installed on a Python 2 system (#1851)
  - Catch ImportError less aggressively in pandas/__init__.py (#1845)
  - Fix pip source installation bug when installing from GitHub (#1805)
  - Fix error when window size > array size in rolling_apply (#1850)
  - Fix pip source installation issues via SSH from GitHub
  - Fix OLS.summary when column is a tuple (#1837)
  - Fix bug in __doc__ patching when -OO passed to interpreter
    (#1792 #1741 #1774)
  - Fix unicode console encoding issue in IPython notebook (#1782, #1768)
  - Fix unicode formatting issue with Series.name (#1782)
  - Fix bug in DataFrame.duplicated with datetime64 columns (#1833)
  - Fix bug in Panel internals resulting in error when doing fillna after
    truncate not changing size of panel (#1823)
  - Prevent segfault due to MultiIndex not being supported in HDFStore table
    format (#1848)
  - Fix UnboundLocalError in Panel.__setitem__ and add better error (#1826)
  - Fix to_csv issues with list of string entries. Isnull works on list of
    strings now too (#1791)
  - Fix Timestamp comparisons with datetime values outside the nanosecond range
    (1677-2262)
  - Revert to prior behavior of normalize_date with datetime.date objects
    (return datetime)
  - Fix broken interaction between np.nansum and Series.any/all
  - Fix bug with multiple column date parsers (#1866)
  - DatetimeIndex.union(Int64Index) was broken
  - Make plot x vs y interface consistent with integer indexing (#1842)
  - set_index inplace modified data even if unique check fails (#1831)
  - Only use Q-OCT/NOV/DEC in quarterly frequency inference (#1789)
  - Upcast to dtype=object when unstacking boolean DataFrame (#1820)
  - Fix float64/float32 merging bug (#1849)
  - Fixes to Period.start_time for non-daily frequencies (#1857)
  - Fix failure when converter used on index_col in read_csv (#1835)
  - Implement PeriodIndex.append so that pandas.concat works correctly (#1815)
  - Avoid Cython out-of-bounds access causing segfault sometimes in pad_2d,
    backfill_2d
  - Fix resampling error with intraday times and anchored target time (like
    AS-DEC) (#1772)
  - Fix .ix indexing bugs with mixed-integer indexes (#1799)
  - Respect passed color keyword argument in Series.plot (#1890)
  - Fix rolling_min/max when the window is larger than the size of the input
    array. Check other malformed inputs (#1899, #1897)
  - Rolling variance / standard deviation with only a single observation in
    window (#1884)
  - Fix unicode sheet name failure in to_excel (#1828)
  - Override DatetimeIndex.min/max to return Timestamp objects (#1895)
  - Fix column name formatting issue in length-truncated column (#1906)
  - Fix broken handling of copying Index metadata to new instances created by
    view(...) calls inside the NumPy infrastructure
  - Support datetime.date again in DateOffset.rollback/rollforward
  - Raise Exception if set passed to Series constructor (#1913)
  - Add TypeError when appending HDFStore table w/ wrong index type (#1881)
  - Don't raise exception on empty inputs in EW functions (e.g. ewma) (#1900)
  - Make asof work correctly with PeriodIndex (#1883)
  - Fix extlinks in doc build
  - Fill boolean DataFrame with NaN when calling shift (#1814)
  - Fix setuptools bug causing pip not to Cythonize .pyx files sometimes
  - Fix negative integer indexing regression in .ix from 0.7.x (#1888)
  - Fix error while retrieving timezone and utc offset from subclasses of
    datetime.tzinfo without .zone and ._utcoffset attributes (#1922)
  - Fix DataFrame formatting of small, non-zero FP numbers (#1911)
  - Various fixes by upcasting of date -> datetime (#1395)
  - Raise better exception when passing multiple functions with the same name,
    such as lambdas, to GroupBy.aggregate
  - Fix DataFrame.apply with axis=1 on a non-unique index (#1878)
  - Proper handling of Index subclasses in pandas.unique (#1759)
  - Set index names in DataFrame.from_records (#1744)
  - Fix time series indexing error with duplicates, under and over hash table
    size cutoff (#1821)
  - Handle list keys in addition to tuples in DataFrame.xs when
    partial-indexing a hierarchically-indexed DataFrame (#1796)
  - Support multiple column selection in DataFrame.__getitem__ with duplicate
    columns (#1943)
  - Fix time zone localization bug causing improper fields (e.g. hours) in time
    zones that have not had a UTC transition in a long time (#1946)
  - Fix errors when parsing and working with with fixed offset timezones
    (#1922, #1928)
  - Fix text parser bug when handling UTC datetime objects generated by
    dateutil (#1693)
  - Fix plotting bug when 'B' is the inferred frequency but index actually
    contains weekends (#1668, #1669)
  - Fix plot styling bugs (#1666, #1665, #1658)
  - Fix plotting bug with index/columns with unicode (#1685)
  - Fix DataFrame constructor bug when passed Series with datetime64 dtype
    in a dict (#1680)
  - Fixed regression in generating DatetimeIndex using timezone aware
    datetime.datetime (#1676)
  - Fix DataFrame bug when printing concatenated DataFrames with duplicated
    columns (#1675)
  - Fixed bug when plotting time series with multiple intraday frequencies
    (#1732)
  - Fix bug in DataFrame.duplicated to enable iterables other than list-types
    as input argument (#1773)
  - Fix resample bug when passed list of lambdas as `how` argument (#1808)
  - Repr fix for MultiIndex level with all NAs (#1971)
  - Fix PeriodIndex slicing bug when slice start/end are out-of-bounds (#1977)
  - Fix read_table bug when parsing unicode (#1975)
  - Fix BlockManager.iget bug when dealing with non-unique MultiIndex as columns
    (#1970)
  - Fix reset_index bug if both drop and level are specified (#1957)
  - Work around unsafe NumPy object->int casting with Cython function (#1987)
  - Fix datetime64 formatting bug in DataFrame.to_csv (#1993)
  - Default start date in pandas.io.data to 1/1/2000 as the docs say (#2011)


pandas 0.8.1
============

**Release date:** July 22, 2012

**New features**

  - Add vectorized, NA-friendly string methods to Series (#1621, #620)
  - Can pass dict of per-column line styles to DataFrame.plot (#1559)
  - Selective plotting to secondary y-axis on same subplot (PR #1640)
  - Add new ``bootstrap_plot`` plot function
  - Add new ``parallel_coordinates`` plot function (#1488)
  - Add ``radviz`` plot function (#1566)
  - Add ``multi_sparse`` option to ``set_printoptions`` to modify display of
    hierarchical indexes (#1538)
  - Add ``dropna`` method to Panel (#171)

**Improvements to existing features**

  - Use moving min/max algorithms from Bottleneck in rolling_min/rolling_max
    for > 100x speedup. (#1504, #50)
  - Add Cython group median method for >15x speedup (#1358)
  - Drastically improve ``to_datetime`` performance on ISO8601 datetime strings
    (with no time zones) (#1571)
  - Improve single-key groupby performance on large data sets, accelerate use of
    groupby with a Categorical variable
  - Add ability to append hierarchical index levels with ``set_index`` and to
    drop single levels with ``reset_index`` (#1569, #1577)
  - Always apply passed functions in ``resample``, even if upsampling (#1596)
  - Avoid unnecessary copies in DataFrame constructor with explicit dtype (#1572)
  - Cleaner DatetimeIndex string representation with 1 or 2 elements (#1611)
  - Improve performance of array-of-Period to PeriodIndex, convert such arrays
    to PeriodIndex inside Index (#1215)
  - More informative string representation for weekly Period objects (#1503)
  - Accelerate 3-axis multi data selection from homogeneous Panel (#979)
  - Add ``adjust`` option to ewma to disable adjustment factor (#1584)
  - Add new matplotlib converters for high frequency time series plotting (#1599)
  - Handling of tz-aware datetime.datetime objects in to_datetime; raise
    Exception unless utc=True given (#1581)

**Bug fixes**

  - Fix NA handling in DataFrame.to_panel (#1582)
  - Handle TypeError issues inside PyObject_RichCompareBool calls in khash
    (#1318)
  - Fix resampling bug to lower case daily frequency (#1588)
  - Fix kendall/spearman DataFrame.corr bug with no overlap (#1595)
  - Fix bug in DataFrame.set_index (#1592)
  - Don't ignore axes in boxplot if by specified (#1565)
  - Fix Panel .ix indexing with integers bug (#1603)
  - Fix Partial indexing bugs (years, months, ...) with PeriodIndex (#1601)
  - Fix MultiIndex console formatting issue (#1606)
  - Unordered index with duplicates doesn't yield scalar location for single
    entry (#1586)
  - Fix resampling of tz-aware time series with "anchored" freq (#1591)
  - Fix DataFrame.rank error on integer data (#1589)
  - Selection of multiple SparseDataFrame columns by list in __getitem__ (#1585)
  - Override Index.tolist for compatibility with MultiIndex (#1576)
  - Fix hierarchical summing bug with MultiIndex of length 1 (#1568)
  - Work around numpy.concatenate use/bug in Series.set_value (#1561)
  - Ensure Series/DataFrame are sorted before resampling (#1580)
  - Fix unhandled IndexError when indexing very large time series (#1562)
  - Fix DatetimeIndex intersection logic error with irregular indexes (#1551)
  - Fix unit test errors on Python 3 (#1550)
  - Fix .ix indexing bugs in duplicate DataFrame index (#1201)
  - Better handle errors with non-existing objects in HDFStore (#1254)
  - Don't copy int64 array data in DatetimeIndex when copy=False (#1624)
  - Fix resampling of conforming periods quarterly to annual (#1622)
  - Don't lose index name on resampling (#1631)
  - Support python-dateutil version 2.1 (#1637)
  - Fix broken scatter_matrix axis labeling, esp. with time series (#1625)
  - Fix cases where extra keywords weren't being passed on to matplotlib from
    Series.plot (#1636)
  - Fix BusinessMonthBegin logic for dates before 1st bday of month (#1645)
  - Ensure string alias converted (valid in DatetimeIndex.get_loc) in
    DataFrame.xs / __getitem__ (#1644)
  - Fix use of string alias timestamps with tz-aware time series (#1647)
  - Fix Series.max/min and Series.describe on len-0 series (#1650)
  - Handle None values in dict passed to concat (#1649)
  - Fix Series.interpolate with method='values' and DatetimeIndex (#1646)
  - Fix IndexError in left merges on a DataFrame with 0-length (#1628)
  - Fix DataFrame column width display with UTF-8 encoded characters (#1620)
  - Handle case in pandas.io.data.get_data_yahoo where Yahoo! returns duplicate
    dates for most recent business day
  - Avoid downsampling when plotting mixed frequencies on the same subplot (#1619)
  - Fix read_csv bug when reading a single line (#1553)
  - Fix bug in C code causing monthly periods prior to December 1969 to be off (#1570)

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
  - Time series string indexing shorthand (#222)
  - Add week, dayofyear array and other timestamp array-valued field accessor
    functions to DatetimeIndex
  - Add GroupBy.prod optimized aggregation function and 'prod' fast time series
    conversion method (#1018)
  - Implement robust frequency inference function and `inferred_freq` attribute
    on DatetimeIndex (#391)
  - New ``tz_convert`` and ``tz_localize`` methods in Series / DataFrame
  - Convert DatetimeIndexes to UTC if time zones are different in join/setops
    (#864)
  - Add limit argument for forward/backward filling to reindex, fillna,
    etc. (#825 and others)
  - Add support for indexes (dates or otherwise) with duplicates and common
    sense indexing/selection functionality
  - Series/DataFrame.update methods, in-place variant of combine_first (#961)
  - Add ``match`` function to API (#502)
  - Add Cython-optimized first, last, min, max, prod functions to GroupBy (#994,
    #1043)
  - Dates can be split across multiple columns (#1227, #1186)
  - Add experimental support for converting pandas DataFrame to R data.frame
    via rpy2 (#350, #1212)
  - Can pass list of (name, function) to GroupBy.aggregate to get aggregates in
    a particular order (#610)
  - Can pass dicts with lists of functions or dicts to GroupBy aggregate to do
    much more flexible multiple function aggregation (#642, #610)
  - New ordered_merge functions for merging DataFrames with ordered
    data. Also supports group-wise merging for panel data (#813)
  - Add keys() method to DataFrame
  - Add flexible replace method for replacing potentially values to Series and
    DataFrame (#929, #1241)
  - Add 'kde' plot kind for Series/DataFrame.plot (#1059)
  - More flexible multiple function aggregation with GroupBy
  - Add pct_change function to Series/DataFrame
  - Add option to interpolate by Index values in Series.interpolate (#1206)
  - Add ``max_colwidth`` option for DataFrame, defaulting to 50
  - Conversion of DataFrame through rpy2 to R data.frame (#1282, )
  - Add keys() method on DataFrame (#1240)
  - Add new ``match`` function to API (similar to R) (#502)
  - Add dayfirst option to parsers (#854)
  - Add ``method`` argument to ``align`` method for forward/backward fillin
    (#216)
  - Add Panel.transpose method for rearranging axes (#695)
  - Add new ``cut`` function (patterned after R) for discretizing data into
    equal range-length bins or arbitrary breaks of your choosing (#415)
  - Add new ``qcut`` for cutting with quantiles (#1378)
  - Add ``value_counts`` top level array method (#1392)
  - Added Andrews curves plot tupe (#1325)
  - Add lag plot (#1440)
  - Add autocorrelation_plot (#1425)
  - Add support for tox and Travis CI (#1382)
  - Add support for Categorical use in GroupBy (#292)
  - Add ``any`` and ``all`` methods to DataFrame (#1416)
  - Add ``secondary_y`` option to Series.plot
  - Add experimental ``lreshape`` function for reshaping wide to long

**Improvements to existing features**

  - Switch to klib/khash-based hash tables in Index classes for better
    performance in many cases and lower memory footprint
  - Shipping some functions from scipy.stats to reduce dependency,
    e.g. Series.describe and DataFrame.describe (GH #1092)
  - Can create MultiIndex by passing list of lists or list of arrays to Series,
    DataFrame constructor, etc. (#831)
  - Can pass arrays in addition to column names to DataFrame.set_index (#402)
  - Improve the speed of "square" reindexing of homogeneous DataFrame objects
    by significant margin (#836)
  - Handle more dtypes when passed MaskedArrays in DataFrame constructor (#406)
  - Improved performance of join operations on integer keys (#682)
  - Can pass multiple columns to GroupBy object, e.g. grouped[[col1, col2]] to
    only aggregate a subset of the value columns (#383)
  - Add histogram / kde plot options for scatter_matrix diagonals (#1237)
  - Add inplace option to Series/DataFrame.rename and sort_index,
    DataFrame.drop_duplicates (#805, #207)
  - More helpful error message when nothing passed to Series.reindex (#1267)
  - Can mix array and scalars as dict-value inputs to DataFrame ctor (#1329)
  - Use DataFrame columns' name for legend title in plots
  - Preserve frequency in DatetimeIndex when possible in boolean indexing
    operations
  - Promote datetime.date values in data alignment operations (#867)
  - Add ``order`` method to Index classes (#1028)
  - Avoid hash table creation in large monotonic hash table indexes (#1160)
  - Store time zones in HDFStore (#1232)
  - Enable storage of sparse data structures in HDFStore (#85)
  - Enable Series.asof to work with arrays of timestamp inputs
  - Cython implementation of DataFrame.corr speeds up by > 100x (#1349, #1354)
  - Exclude "nuisance" columns automatically in GroupBy.transform (#1364)
  - Support functions-as-strings in GroupBy.transform (#1362)
  - Use index name as xlabel/ylabel in plots (#1415)
  - Add ``convert_dtype`` option to Series.apply to be able to leave data as
    dtype=object (#1414)
  - Can specify all index level names in concat (#1419)
  - Add ``dialect`` keyword to parsers for quoting conventions (#1363)
  - Enable DataFrame[bool_DataFrame] += value (#1366)
  - Add ``retries`` argument to ``get_data_yahoo`` to try to prevent Yahoo! API
    404s (#826)
  - Improve performance of reshaping by using O(N) categorical sorting
  - Series names will be used for index of DataFrame if no index passed (#1494)
  - Header argument in DataFrame.to_csv can accept a list of column names to
    use instead of the object's columns (#921)
  - Add ``raise_conflict`` argument to DataFrame.update (#1526)
  - Support file-like objects in ExcelFile (#1529)

**API Changes**

  - Rename `pandas._tseries` to `pandas.lib`
  - Rename Factor to Categorical and add improvements. Numerous Categorical bug
    fixes
  - Frequency name overhaul, WEEKDAY/EOM and rules with @
    deprecated. get_legacy_offset_name backwards compatibility function added
  - Raise ValueError in DataFrame.__nonzero__, so "if df" no longer works
    (#1073)
  - Change BDay (business day) to not normalize dates by default (#506)
  - Remove deprecated DataMatrix name
  - Default merge suffixes for overlap now have underscores instead of periods
    to facilitate tab completion, etc. (#1239)
  - Deprecation of offset, time_rule timeRule parameters throughout codebase
  - Series.append and DataFrame.append no longer check for duplicate indexes
    by default, add verify_integrity parameter (#1394)
  - Refactor Factor class, old constructor moved to Factor.from_array
  - Modified internals of MultiIndex to use less memory (no longer represented
    as array of tuples) internally, speed up construction time and many methods
    which construct intermediate hierarchical indexes (#1467)

**Bug fixes**

  - Fix OverflowError from storing pre-1970 dates in HDFStore by switching to
    datetime64 (GH #179)
  - Fix logical error with February leap year end in YearEnd offset
  - Series([False, nan]) was getting casted to float64 (GH #1074)
  - Fix binary operations between boolean Series and object Series with
    booleans and NAs (GH #1074, #1079)
  - Couldn't assign whole array to column in mixed-type DataFrame via .ix
    (#1142)
  - Fix label slicing issues with float index values (#1167)
  - Fix segfault caused by empty groups passed to groupby (#1048)
  - Fix occasionally misbehaved reindexing in the presence of NaN labels (#522)
  - Fix imprecise logic causing weird Series results from .apply (#1183)
  - Unstack multiple levels in one shot, avoiding empty columns in some
    cases. Fix pivot table bug (#1181)
  - Fix formatting of MultiIndex on Series/DataFrame when index name coincides
    with label (#1217)
  - Handle Excel 2003 #N/A as NaN from xlrd (#1213, #1225)
  - Fix timestamp locale-related deserialization issues with HDFStore by moving
    to datetime64 representation (#1081, #809)
  - Fix DataFrame.duplicated/drop_duplicates NA value handling (#557)
  - Actually raise exceptions in fast reducer (#1243)
  - Fix various timezone-handling bugs from 0.7.3 (#969)
  - GroupBy on level=0 discarded index name (#1313)
  - Better error message with unmergeable DataFrames (#1307)
  - Series.__repr__ alignment fix with unicode index values (#1279)
  - Better error message if nothing passed to reindex (#1267)
  - More robust NA handling in DataFrame.drop_duplicates (#557)
  - Resolve locale-based and pre-epoch HDF5 timestamp deserialization issues
    (#973, #1081, #179)
  - Implement Series.repeat (#1229)
  - Fix indexing with namedtuple and other tuple subclasses (#1026)
  - Fix float64 slicing bug (#1167)
  - Parsing integers with commas (#796)
  - Fix groupby improper data type when group consists of one value (#1065)
  - Fix negative variance possibility in nanvar resulting from floating point
    error (#1090)
  - Consistently set name on groupby pieces (#184)
  - Treat dict return values as Series in GroupBy.apply (#823)
  - Respect column selection for DataFrame in in GroupBy.transform (#1365)
  - Fix MultiIndex partial indexing bug (#1352)
  - Enable assignment of rows in mixed-type DataFrame via .ix (#1432)
  - Reset index mapping when grouping Series in Cython (#1423)
  - Fix outer/inner DataFrame.join with non-unique indexes (#1421)
  - Fix MultiIndex groupby bugs with empty lower levels (#1401)
  - Calling fillna with a Series will have same behavior as with dict (#1486)
  - SparseSeries reduction bug (#1375)
  - Fix unicode serialization issue in HDFStore (#1361)
  - Pass keywords to pyplot.boxplot in DataFrame.boxplot (#1493)
  - Bug fixes in MonthBegin (#1483)
  - Preserve MultiIndex names in drop (#1513)
  - Fix Panel DataFrame slice-assignment bug (#1533)
  - Don't use locals() in read_* functions (#1547)

pandas 0.7.3
============

**Release date:** April 12, 2012

**New features / modules**

  - Support for non-unique indexes: indexing and selection, many-to-one and
    many-to-many joins (#1306)
  - Added fixed-width file reader, read_fwf (PR #952)
  - Add group_keys argument to groupby to not add group names to MultiIndex in
    result of apply (GH #938)
  - DataFrame can now accept non-integer label slicing (GH #946). Previously
    only DataFrame.ix was able to do so.
  - DataFrame.apply now retains name attributes on Series objects (GH #983)
  - Numeric DataFrame comparisons with non-numeric values now raises proper
    TypeError (GH #943). Previously raise "PandasError: DataFrame constructor
    not properly called!"
  - Add ``kurt`` methods to Series and DataFrame (PR #964)
  - Can pass dict of column -> list/set NA values for text parsers (GH #754)
  - Allows users specified NA values in text parsers (GH #754)
  - Parsers checks for openpyxl dependency and raises ImportError if not found
    (PR #1007)
  - New factory function to create HDFStore objects that can be used in a with
    statement so users do not have to explicitly call HDFStore.close (PR #1005)
  - pivot_table is now more flexible with same parameters as groupby (GH #941)
  - Added stacked bar plots (GH #987)
  - scatter_matrix method in pandas/tools/plotting.py (PR #935)
  - DataFrame.boxplot returns plot results for ex-post styling (GH #985)
  - Short version number accessible as pandas.version.short_version (GH #930)
  - Additional documentation in panel.to_frame (GH #942)
  - More informative Series.apply docstring regarding element-wise apply
    (GH #977)
  - Notes on rpy2 installation (GH #1006)
  - Add rotation and font size options to hist method (#1012)
  - Use exogenous / X variable index in result of OLS.y_predict. Add
    OLS.predict method (PR #1027, #1008)

**API Changes**

  - Calling apply on grouped Series, e.g. describe(), will no longer yield
    DataFrame by default. Will have to call unstack() to get prior behavior
  - NA handling in non-numeric comparisons has been tightened up (#933, #953)
  - No longer assign dummy names key_0, key_1, etc. to groupby index (#1291)

**Bug fixes**

  - Fix logic error when selecting part of a row in a DataFrame with a
    MultiIndex index (GH #1013)
  - Series comparison with Series of differing length causes crash (GH #1016).
  - Fix bug in indexing when selecting section of hierarchically-indexed row
    (GH #1013)
  - DataFrame.plot(logy=True) has no effect (GH #1011).
  - Broken arithmetic operations between SparsePanel-Panel (GH #1015)
  - Unicode repr issues in MultiIndex with non-ascii characters (GH #1010)
  - DataFrame.lookup() returns inconsistent results if exact match not present
    (GH #1001)
  - DataFrame arithmetic operations not treating None as NA (GH #992)
  - DataFrameGroupBy.apply returns incorrect result (GH #991)
  - Series.reshape returns incorrect result for multiple dimensions (GH #989)
  - Series.std and Series.var ignores ddof parameter (GH #934)
  - DataFrame.append loses index names (GH #980)
  - DataFrame.plot(kind='bar') ignores color argument (GH #958)
  - Inconsistent Index comparison results (GH #948)
  - Improper int dtype DataFrame construction from data with NaN (GH #846)
  - Removes default 'result' name in grouby results (GH #995)
  - DataFrame.from_records no longer mutate input columns (PR #975)
  - Use Index name when grouping by it (#1313)

pandas 0.7.2
============

**Release date:** March 16, 2012

**New features / modules**

  - Add additional tie-breaking methods in DataFrame.rank (#874)
  - Add ascending parameter to rank in Series, DataFrame (#875)
  - Add coerce_float option to DataFrame.from_records (#893)
  - Add sort_columns parameter to allow unsorted plots (#918)
  - IPython tab completion on GroupBy objects

**API Changes**

  - Series.sum returns 0 instead of NA when called on an empty
    series. Analogously for a DataFrame whose rows or columns are length 0
    (#844)

**Improvements to existing features**

  - Don't use groups dict in Grouper.size (#860)
  - Use khash for Series.value_counts, add raw function to algorithms.py (#861)
  - Enable column access via attributes on GroupBy (#882)
  - Enable setting existing columns (only) via attributes on DataFrame, Panel
    (#883)
  - Intercept __builtin__.sum in groupby (#885)
  - Can pass dict to DataFrame.fillna to use different values per column (#661)
  - Can select multiple hierarchical groups by passing list of values in .ix
    (#134)
  - Add level keyword to ``drop`` for dropping values from a level (GH #159)
  - Add ``coerce_float`` option on DataFrame.from_records (# 893)
  - Raise exception if passed date_parser fails in ``read_csv``
  - Add ``axis`` option to DataFrame.fillna (#174)
  - Fixes to Panel to make it easier to subclass (PR #888)

**Bug fixes**

  - Fix overflow-related bugs in groupby (#850, #851)
  - Fix unhelpful error message in parsers (#856)
  - Better err msg for failed boolean slicing of dataframe (#859)
  - Series.count cannot accept a string (level name) in the level argument (#869)
  - Group index platform int check (#870)
  - concat on axis=1 and ignore_index=True raises TypeError (#871)
  - Further unicode handling issues resolved (#795)
  - Fix failure in multiindex-based access in Panel (#880)
  - Fix DataFrame boolean slice assignment failure (#881)
  - Fix combineAdd NotImplementedError for SparseDataFrame (#887)
  - Fix DataFrame.to_html encoding and columns (#890, #891, #909)
  - Fix na-filling handling in mixed-type DataFrame (#910)
  - Fix to DataFrame.set_value with non-existant row/col (#911)
  - Fix malformed block in groupby when excluding nuisance columns (#916)
  - Fix inconsistant NA handling in dtype=object arrays (#925)
  - Fix missing center-of-mass computation in ewmcov (#862)
  - Don't raise exception when opening read-only HDF5 file (#847)
  - Fix possible out-of-bounds memory access in 0-length Series (#917)

pandas 0.7.1
============

**Release date:** February 29, 2012

**New features / modules**

  - Add ``to_clipboard`` function to pandas namespace for writing objects to
    the system clipboard (#774)
  - Add ``itertuples`` method to DataFrame for iterating through the rows of a
    dataframe as tuples (#818)
  - Add ability to pass fill_value and method to DataFrame and Series align
    method (#806, #807)
  - Add fill_value option to reindex, align methods (#784)
  - Enable concat to produce DataFrame from Series (#787)
  - Add ``between`` method to Series (#802)
  - Add HTML representation hook to DataFrame for the IPython HTML notebook
    (#773)
  - Support for reading Excel 2007 XML documents using openpyxl

**Improvements to existing features**

  - Improve performance and memory usage of fillna on DataFrame
  - Can concatenate a list of Series along axis=1 to obtain a DataFrame (#787)

**Bug fixes**

  - Fix memory leak when inserting large number of columns into a single
    DataFrame (#790)
  - Appending length-0 DataFrame with new columns would not result in those new
    columns being part of the resulting concatenated DataFrame (#782)
  - Fixed groupby corner case when passing dictionary grouper and as_index is
    False (#819)
  - Fixed bug whereby bool array sometimes had object dtype (#820)
  - Fix exception thrown on np.diff (#816)
  - Fix to_records where columns are non-strings (#822)
  - Fix Index.intersection where indices have incomparable types (#811)
  - Fix ExcelFile throwing an exception for two-line file (#837)
  - Add clearer error message in csv parser (#835)
  - Fix loss of fractional seconds in HDFStore (#513)
  - Fix DataFrame join where columns have datetimes (#787)
  - Work around numpy performance issue in take (#817)
  - Improve comparison operations for NA-friendliness (#801)
  - Fix indexing operation for floating point values (#780, #798)
  - Fix groupby case resulting in malformed dataframe (#814)
  - Fix behavior of reindex of Series dropping name (#812)
  - Improve on redudant groupby computation (#775)
  - Catch possible NA assignment to int/bool series with exception (#839)

pandas 0.7.0
============

**Release date:** 2/9/2012

**New features / modules**

  - New ``merge`` function for efficiently performing full gamut of database /
    relational-algebra operations. Refactored existing join methods to use the
    new infrastructure, resulting in substantial performance gains (GH #220,
    #249, #267)
  - New ``concat`` function for concatenating DataFrame or Panel objects along
    an axis. Can form union or intersection of the other axes. Improves
    performance of ``DataFrame.append`` (#468, #479, #273)
  - Handle differently-indexed output values in ``DataFrame.apply`` (GH #498)
  - Can pass list of dicts (e.g., a list of shallow JSON objects) to DataFrame
    constructor (GH #526)
  - Add ``reorder_levels`` method to Series and DataFrame (PR #534)
  - Add dict-like ``get`` function to DataFrame and Panel (PR #521)
  - ``DataFrame.iterrows`` method for efficiently iterating through the rows of
    a DataFrame
  - Added ``DataFrame.to_panel`` with code adapted from ``LongPanel.to_long``
  - ``reindex_axis`` method added to DataFrame
  - Add ``level`` option to binary arithmetic functions on ``DataFrame`` and
    ``Series``
  - Add ``level`` option to the ``reindex`` and ``align`` methods on Series and
    DataFrame for broadcasting values across a level (GH #542, PR #552, others)
  - Add attribute-based item access to ``Panel`` and add IPython completion (PR
    #554)
  - Add ``logy`` option to ``Series.plot`` for log-scaling on the Y axis
  - Add ``index``, ``header``, and ``justify`` options to
    ``DataFrame.to_string``. Add option to   (GH #570, GH #571)
  - Can pass multiple DataFrames to ``DataFrame.join`` to join on index (GH #115)
  - Can pass multiple Panels to ``Panel.join`` (GH #115)
  - Can pass multiple DataFrames to `DataFrame.append` to concatenate (stack)
    and multiple Series to ``Series.append`` too
  - Added ``justify`` argument to ``DataFrame.to_string`` to allow different
    alignment of column headers
  - Add ``sort`` option to GroupBy to allow disabling sorting of the group keys
    for potential speedups (GH #595)
  - Can pass MaskedArray to Series constructor (PR #563)
  - Add Panel item access via attributes and IPython completion (GH #554)
  - Implement ``DataFrame.lookup``, fancy-indexing analogue for retrieving
    values given a sequence of row and column labels (GH #338)
  - Add ``verbose`` option to ``read_csv`` and ``read_table`` to show number of
    NA values inserted in non-numeric columns (GH #614)
  - Can pass a list of dicts or Series to ``DataFrame.append`` to concatenate
    multiple rows (GH #464)
  - Add ``level`` argument to ``DataFrame.xs`` for selecting data from other
    MultiIndex levels. Can take one or more levels with potentially a tuple of
    keys for flexible retrieval of data (GH #371, GH #629)
  - New ``crosstab`` function for easily computing frequency tables (GH #170)
  - Can pass a list of functions to aggregate with groupby on a DataFrame,
    yielding an aggregated result with hierarchical columns (GH #166)
  - Add integer-indexing functions ``iget`` in Series and ``irow`` / ``iget``
    in DataFrame (GH #628)
  - Add new ``Series.unique`` function, significantly faster than
    ``numpy.unique`` (GH #658)
  - Add new ``cummin`` and ``cummax`` instance methods to ``Series`` and
    ``DataFrame`` (GH #647)
  - Add new ``value_range`` function to return min/max of a dataframe (GH #288)
  - Add ``drop`` parameter to ``reset_index`` method of ``DataFrame`` and added
    method to ``Series`` as well (GH #699)
  - Add ``isin`` method to Index objects, works just like ``Series.isin`` (GH
    #657)
  - Implement array interface on Panel so that ufuncs work (re: #740)
  - Add ``sort`` option to ``DataFrame.join`` (GH #731)
  - Improved handling of NAs (propagation) in binary operations with
    dtype=object arrays (GH #737)
  - Add ``abs`` method to Pandas objects
  - Added ``algorithms`` module to start collecting central algos

**API Changes**

  - Label-indexing with integer indexes now raises KeyError if a label is not
    found instead of falling back on location-based indexing (GH #700)
  - Label-based slicing via ``ix`` or ``[]`` on Series will now only work if
    exact matches for the labels are found or if the index is monotonic (for
    range selections)
  - Label-based slicing and sequences of labels can be passed to ``[]`` on a
    Series for both getting and setting (GH #86)
  - `[]` operator (``__getitem__`` and ``__setitem__``) will raise KeyError
    with integer indexes when an index is not contained in the index. The prior
    behavior would fall back on position-based indexing if a key was not found
    in the index which would lead to subtle bugs. This is now consistent with
    the behavior of ``.ix`` on DataFrame and friends (GH #328)
  - Rename ``DataFrame.delevel`` to ``DataFrame.reset_index`` and add
    deprecation warning
  - `Series.sort` (an in-place operation) called on a Series which is a view on
    a larger array (e.g. a column in a DataFrame) will generate an Exception to
    prevent accidentally modifying the data source (GH #316)
  - Refactor to remove deprecated ``LongPanel`` class (PR #552)
  - Deprecated ``Panel.to_long``, renamed to ``to_frame``
  - Deprecated ``colSpace`` argument in ``DataFrame.to_string``, renamed to
    ``col_space``
  - Rename ``precision`` to ``accuracy`` in engineering float formatter (GH
    #395)
  - The default delimiter for ``read_csv`` is comma rather than letting
    ``csv.Sniffer`` infer it
  - Rename ``col_or_columns`` argument in ``DataFrame.drop_duplicates`` (GH
    #734)

**Improvements to existing features**

  - Better error message in DataFrame constructor when passed column labels
    don't match data (GH #497)
  - Substantially improve performance of multi-GroupBy aggregation when a
    Python function is passed, reuse ndarray object in Cython (GH #496)
  - Can store objects indexed by tuples and floats in HDFStore (GH #492)
  - Don't print length by default in Series.to_string, add `length` option (GH
    #489)
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
    also (GH #536)
  - Default name assignment when calling ``reset_index`` on DataFrame with a
    regular (non-hierarchical) index (GH #476)
  - Use Cythonized groupers when possible in Series/DataFrame stat ops with
    ``level`` parameter passed (GH #545)
  - Ported skiplist data structure to C to speed up ``rolling_median`` by about
    5-10x in most typical use cases (GH #374)
  - Some performance enhancements in constructing a Panel from a dict of
    DataFrame objects
  - Made ``Index._get_duplicates`` a public method by removing the underscore
  - Prettier printing of floats, and column spacing fix (GH #395, GH #571)
  - Add ``bold_rows`` option to DataFrame.to_html (GH #586)
  - Improve the performance of ``DataFrame.sort_index`` by up to 5x or more
    when sorting by multiple columns
  - Substantially improve performance of DataFrame and Series constructors when
    passed a nested dict or dict, respectively (GH #540, GH #621)
  - Modified setup.py so that pip / setuptools will install dependencies (GH
    #507, various pull requests)
  - Unstack called on DataFrame with non-MultiIndex will return Series (GH
    #477)
  - Improve DataFrame.to_string and console formatting to be more consistent in
    the number of displayed digits (GH #395)
  - Use bottleneck if available for performing NaN-friendly statistical
    operations that it implemented (GH #91)
  - Monkey-patch context to traceback in ``DataFrame.apply`` to indicate which
    row/column the function application failed on (GH #614)
  - Improved ability of read_table and read_clipboard to parse
    console-formatted DataFrames (can read the row of index names, etc.)
  - Can pass list of group labels (without having to convert to an ndarray
    yourself) to ``groupby`` in some cases (GH #659)
  - Use ``kind`` argument to Series.order for selecting different sort kinds
    (GH #668)
  - Add option to Series.to_csv to omit the index (PR #684)
  - Add ``delimiter`` as an alternative to ``sep`` in ``read_csv`` and other
    parsing functions
  - Substantially improved performance of groupby on DataFrames with many
    columns by aggregating blocks of columns all at once (GH #745)
  - Can pass a file handle or StringIO to Series/DataFrame.to_csv (GH #765)
  - Can pass sequence of integers to DataFrame.irow(icol) and Series.iget, (GH
    #654)
  - Prototypes for some vectorized string functions
  - Add float64 hash table to solve the Series.unique problem with NAs (GH #714)
  - Memoize objects when reading from file to reduce memory footprint
  - Can get and set a column of a DataFrame with hierarchical columns
    containing "empty" ('') lower levels without passing the empty levels (PR
    #768)

**Bug fixes**

  - Raise exception in out-of-bounds indexing of Series instead of
    seg-faulting, regression from earlier releases (GH #495)
  - Fix error when joining DataFrames of different dtypes within the same
    typeclass (e.g. float32 and float64) (GH #486)
  - Fix bug in Series.min/Series.max on objects like datetime.datetime (GH
    #487)
  - Preserve index names in Index.union (GH #501)
  - Fix bug in Index joining causing subclass information (like DateRange type)
    to be lost in some cases (GH #500)
  - Accept empty list as input to DataFrame constructor, regression from 0.6.0
    (GH #491)
  - Can output DataFrame and Series with ndarray objects in a dtype=object
    array (GH #490)
  - Return empty string from Series.to_string when called on empty Series (GH
    #488)
  - Fix exception passing empty list to DataFrame.from_records
  - Fix Index.format bug (excluding name field) with datetimes with time info
  - Fix scalar value access in Series to always return NumPy scalars,
    regression from prior versions (GH #510)
  - Handle rows skipped at beginning of file in read_* functions (GH #505)
  - Handle improper dtype casting in ``set_value`` methods
  - Unary '-' / __neg__ operator on DataFrame was returning integer values
  - Unbox 0-dim ndarrays from certain operators like all, any in Series
  - Fix handling of missing columns (was combine_first-specific) in
    DataFrame.combine for general case (GH #529)
  - Fix type inference logic with boolean lists and arrays in DataFrame indexing
  - Use centered sum of squares in R-square computation if entity_effects=True
    in panel regression
  - Handle all NA case in Series.{corr, cov}, was raising exception (GH #548)
  - Aggregating by multiple levels with ``level`` argument to DataFrame, Series
    stat method, was broken (GH #545)
  - Fix Cython buf when converter passed to read_csv produced a numeric array
    (buffer dtype mismatch when passed to Cython type inference function) (GH
    #546)
  - Fix exception when setting scalar value using .ix on a DataFrame with a
    MultiIndex (GH #551)
  - Fix outer join between two DateRanges with different offsets that returned
    an invalid DateRange
  - Cleanup DataFrame.from_records failure where index argument is an integer
  - Fix Data.from_records failure when passed a dictionary
  - Fix NA handling in {Series, DataFrame}.rank with non-floating point dtypes
  - Fix bug related to integer type-checking in .ix-based indexing
  - Handle non-string index name passed to DataFrame.from_records
  - DataFrame.insert caused the columns name(s) field to be discarded (GH #527)
  - Fix erroneous in monotonic many-to-one left joins
  - Fix DataFrame.to_string to remove extra column white space (GH #571)
  - Format floats to default to same number of digits (GH #395)
  - Added decorator to copy docstring from one function to another (GH #449)
  - Fix error in monotonic many-to-one left joins
  - Fix __eq__ comparison between DateOffsets with different relativedelta
    keywords passed
  - Fix exception caused by parser converter returning strings (GH #583)
  - Fix MultiIndex formatting bug with integer names (GH #601)
  - Fix bug in handling of non-numeric aggregates in Series.groupby (GH #612)
  - Fix TypeError with tuple subclasses (e.g. namedtuple) in
    DataFrame.from_records (GH #611)
  - Catch misreported console size when running IPython within Emacs
  - Fix minor bug in pivot table margins, loss of index names and length-1
    'All' tuple in row labels
  - Add support for legacy WidePanel objects to be read from HDFStore
  - Fix out-of-bounds segfault in pad_object and backfill_object methods when
    either source or target array are empty
  - Could not create a new column in a DataFrame from a list of tuples
  - Fix bugs preventing SparseDataFrame and SparseSeries working with groupby
    (GH #666)
  - Use sort kind in Series.sort / argsort (GH #668)
  - Fix DataFrame operations on non-scalar, non-pandas objects (GH #672)
  - Don't convert DataFrame column to integer type when passing integer to
    __setitem__ (GH #669)
  - Fix downstream bug in pivot_table caused by integer level names in
    MultiIndex (GH #678)
  - Fix SparseSeries.combine_first when passed a dense Series (GH #687)
  - Fix performance regression in HDFStore loading when DataFrame or Panel
    stored in table format with datetimes
  - Raise Exception in DateRange when offset with n=0 is passed (GH #683)
  - Fix get/set inconsistency with .ix property and integer location but
    non-integer index (GH #707)
  - Use right dropna function for SparseSeries. Return dense Series for NA fill
    value (GH #730)
  - Fix Index.format bug causing incorrectly string-formatted Series with
    datetime indexes (# 726, 758)
  - Fix errors caused by object dtype arrays passed to ols (GH #759)
  - Fix error where column names lost when passing list of labels to
    DataFrame.__getitem__, (GH #662)
  - Fix error whereby top-level week iterator overwrote week instance
  - Fix circular reference causing memory leak in sparse array / series /
    frame, (GH #663)
  - Fix integer-slicing from integers-as-floats (GH #670)
  - Fix zero division errors in nanops from object dtype arrays in all NA case
    (GH #676)
  - Fix csv encoding when using unicode (GH #705, #717, #738)
  - Fix assumption that each object contains every unique block type in concat,
    (GH #708)
  - Fix sortedness check of multiindex in to_panel (GH #719, 720)
  - Fix that None was not treated as NA in PyObjectHashtable
  - Fix hashing dtype because of endianness confusion (GH #747, #748)
  - Fix SparseSeries.dropna to return dense Series in case of NA fill value (GH
    #730)
  - Use map_infer instead of np.vectorize. handle NA sentinels if converter
    yields numeric array, (GH #753)
  - Fixes and improvements to DataFrame.rank (GH #742)
  - Fix catching AttributeError instead of NameError for bottleneck
  - Try to cast non-MultiIndex to better dtype when calling reset_index (GH #726
    #440)
  - Fix #1.QNAN0' float bug on 2.6/win64
  - Allow subclasses of dicts in DataFrame constructor, with tests
  - Fix problem whereby set_index destroys column multiindex (GH #764)
  - Hack around bug in generating DateRange from naive DateOffset (GH #770)
  - Fix bug in DateRange.intersection causing incorrect results with some
    overlapping ranges (GH #771)

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

pandas 0.6.1
============

**Release date:** 12/13/2011

**API Changes**

  - Rename `names` argument in DataFrame.from_records to `columns`. Add
    deprecation warning
  - Boolean get/set operations on Series with boolean Series will reindex
    instead of requiring that the indexes be exactly equal (GH #429)

**New features / modules**

  - Can pass Series to DataFrame.append with ignore_index=True for appending a
    single row (GH #430)
  - Add Spearman and Kendall correlation options to Series.corr and
    DataFrame.corr (GH #428)
  - Add new `get_value` and `set_value` methods to Series, DataFrame, and Panel
    to very low-overhead access to scalar elements. df.get_value(row, column)
    is about 3x faster than df[column][row] by handling fewer cases (GH #437,
    #438). Add similar methods to sparse data structures for compatibility
  - Add Qt table widget to sandbox (PR #435)
  - DataFrame.align can accept Series arguments, add axis keyword (GH #461)
  - Implement new SparseList and SparseArray data structures. SparseSeries now
    derives from SparseArray (GH #463)
  - max_columns / max_rows options in set_printoptions (PR #453)
  - Implement Series.rank and DataFrame.rank, fast versions of
    scipy.stats.rankdata (GH #428)
  - Implement DataFrame.from_items alternate constructor (GH #444)
  - DataFrame.convert_objects method for inferring better dtypes for object
    columns (GH #302)
  - Add rolling_corr_pairwise function for computing Panel of correlation
    matrices (GH #189)
  - Add `margins` option to `pivot_table` for computing subgroup aggregates (GH
    #114)
  - Add `Series.from_csv` function (PR #482)

**Improvements to existing features**

  - Improve memory usage of `DataFrame.describe` (do not copy data
    unnecessarily) (PR #425)
  - Use same formatting function for outputting floating point Series to console
    as in DataFrame (PR #420)
  - DataFrame.delevel will try to infer better dtype for new columns (GH #440)
  - Exclude non-numeric types in DataFrame.{corr, cov}
  - Override Index.astype to enable dtype casting (GH #412)
  - Use same float formatting function for Series.__repr__ (PR #420)
  - Use available console width to output DataFrame columns (PR #453)
  - Accept ndarrays when setting items in Panel (GH #452)
  - Infer console width when printing __repr__ of DataFrame to console (PR
    #453)
  - Optimize scalar value lookups in the general case by 25% or more in Series
    and DataFrame
  - Can pass DataFrame/DataFrame and DataFrame/Series to
    rolling_corr/rolling_cov (GH #462)
  - Fix performance regression in cross-sectional count in DataFrame, affecting
    DataFrame.dropna speed
  - Column deletion in DataFrame copies no data (computes views on blocks) (GH
    #158)
  - MultiIndex.get_level_values can take the level name
  - More helpful error message when DataFrame.plot fails on one of the columns
    (GH #478)
  - Improve performance of DataFrame.{index, columns} attribute lookup

**Bug fixes**

  - Fix O(K^2) memory leak caused by inserting many columns without
    consolidating, had been present since 0.4.0 (GH #467)
  - `DataFrame.count` should return Series with zero instead of NA with length-0
    axis (GH #423)
  - Fix Yahoo! Finance API usage in pandas.io.data (GH #419, PR #427)
  - Fix upstream bug causing failure in Series.align with empty Series (GH #434)
  - Function passed to DataFrame.apply can return a list, as long as it's the
    right length. Regression from 0.4 (GH #432)
  - Don't "accidentally" upcast scalar values when indexing using .ix (GH #431)
  - Fix groupby exception raised with as_index=False and single column selected
    (GH #421)
  - Implement DateOffset.__ne__ causing downstream bug (GH #456)
  - Fix __doc__-related issue when converting py -> pyo with py2exe
  - Bug fix in left join Cython code with duplicate monotonic labels
  - Fix bug when unstacking multiple levels described in #451
  - Exclude NA values in dtype=object arrays, regression from 0.5.0 (GH #469)
  - Use Cython map_infer function in DataFrame.applymap to properly infer
    output type, handle tuple return values and other things that were breaking
    (GH #465)
  - Handle floating point index values in HDFStore (GH #454)
  - Fixed stale column reference bug (cached Series object) caused by type
    change / item deletion in DataFrame (GH #473)
  - Index.get_loc should always raise Exception when there are duplicates
  - Handle differently-indexed Series input to DataFrame constructor (GH #475)
  - Omit nuisance columns in multi-groupby with Python function
  - Buglet in handling of single grouping in general apply
  - Handle type inference properly when passing list of lists or tuples to
    DataFrame constructor (GH #484)
  - Preserve Index / MultiIndex names in GroupBy.apply concatenation step (GH
    #481)

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

pandas 0.6.0
============

**Release date:** 11/25/2011

**API Changes**

  - Arithmetic methods like `sum` will attempt to sum dtype=object values by
    default instead of excluding them (GH #382)

**New features / modules**

  - Add `melt` function to `pandas.core.reshape`
  - Add `level` parameter to group by level in Series and DataFrame
    descriptive statistics (PR #313)
  - Add `head` and `tail` methods to Series, analogous to to DataFrame (PR
    #296)
  - Add `Series.isin` function which checks if each value is contained in a
    passed sequence (GH #289)
  - Add `float_format` option to `Series.to_string`
  - Add `skip_footer` (GH #291) and `converters` (GH #343) options to
    `read_csv` and `read_table`
  - Add proper, tested weighted least squares to standard and panel OLS (GH
    #303)
  - Add `drop_duplicates` and `duplicated` functions for removing duplicate
    DataFrame rows and checking for duplicate rows, respectively (GH #319)
  - Implement logical (boolean) operators &, |, ^ on DataFrame (GH #347)
  - Add `Series.mad`, mean absolute deviation, matching DataFrame
  - Add `QuarterEnd` DateOffset (PR #321)
  - Add matrix multiplication function `dot` to DataFrame (GH #65)
  - Add `orient` option to `Panel.from_dict` to ease creation of mixed-type
    Panels (GH #359, #301)
  - Add `DataFrame.from_dict` with similar `orient` option
  - Can now pass list of tuples or list of lists to `DataFrame.from_records`
    for fast conversion to DataFrame (GH #357)
  - Can pass multiple levels to groupby, e.g. `df.groupby(level=[0, 1])` (GH
    #103)
  - Can sort by multiple columns in `DataFrame.sort_index` (GH #92, PR #362)
  - Add fast `get_value` and `put_value` methods to DataFrame and
    micro-performance tweaks (GH #360)
  - Add `cov` instance methods to Series and DataFrame (GH #194, PR #362)
  - Add bar plot option to `DataFrame.plot` (PR #348)
  - Add `idxmin` and `idxmax` functions to Series and DataFrame for computing
    index labels achieving maximum and minimum values (PR #286)
  - Add `read_clipboard` function for parsing DataFrame from OS clipboard,
    should work across platforms (GH #300)
  - Add `nunique` function to Series for counting unique elements (GH #297)
  - DataFrame constructor will use Series name if no columns passed (GH #373)
  - Support regular expressions and longer delimiters in read_table/read_csv,
    but does not handle quoted strings yet (GH #364)
  - Add `DataFrame.to_html` for formatting DataFrame to HTML (PR #387)
  - MaskedArray can be passed to DataFrame constructor and masked values will be
    converted to NaN (PR #396)
  - Add `DataFrame.boxplot` function (GH #368, others)
  - Can pass extra args, kwds to DataFrame.apply (GH #376)

**Improvements to existing features**

  - Raise more helpful exception if date parsing fails in DateRange (GH #298)
  - Vastly improved performance of GroupBy on axes with a MultiIndex (GH #299)
  - Print level names in hierarchical index in Series repr (GH #305)
  - Return DataFrame when performing GroupBy on selected column and
    as_index=False (GH #308)
  - Can pass vector to `on` argument in `DataFrame.join` (GH #312)
  - Don't show Series name if it's None in the repr, also omit length for short
    Series (GH #317)
  - Show legend by default in `DataFrame.plot`, add `legend` boolean flag (GH
    #324)
  - Significantly improved performance of `Series.order`, which also makes
    np.unique called on a Series faster (GH #327)
  - Faster cythonized count by level in Series and DataFrame (GH #341)
  - Raise exception if dateutil 2.0 installed on Python 2.x runtime (GH #346)
  - Significant GroupBy performance enhancement with multiple keys with many
    "empty" combinations
  - New Cython vectorized function `map_infer` speeds up `Series.apply` and
    `Series.map` significantly when passed elementwise Python function,
    motivated by PR #355
  - Cythonized `cache_readonly`, resulting in substantial micro-performance
    enhancements throughout the codebase (GH #361)
  - Special Cython matrix iterator for applying arbitrary reduction operations
    with 3-5x better performance than `np.apply_along_axis` (GH #309)
  - Add `raw` option to `DataFrame.apply` for getting better performance when
    the passed function only requires an ndarray (GH #309)
  - Improve performance of `MultiIndex.from_tuples`
  - Can pass multiple levels to `stack` and `unstack` (GH #370)
  - Can pass multiple values columns to `pivot_table` (GH #381)
  - Can call `DataFrame.delevel` with standard Index with name set (GH #393)
  - Use Series name in GroupBy for result index (GH #363)
  - Refactor Series/DataFrame stat methods to use common set of NaN-friendly
    function
  - Handle NumPy scalar integers at C level in Cython conversion routines

**Bug fixes**

  - Fix bug in `DataFrame.to_csv` when writing a DataFrame with an index
    name (GH #290)
  - DataFrame should clear its Series caches on consolidation, was causing
    "stale" Series to be returned in some corner cases (GH #304)
  - DataFrame constructor failed if a column had a list of tuples (GH #293)
  - Ensure that `Series.apply` always returns a Series and implement
    `Series.round` (GH #314)
  - Support boolean columns in Cythonized groupby functions (GH #315)
  - `DataFrame.describe` should not fail if there are no numeric columns,
    instead return categorical describe (GH #323)
  - Fixed bug which could cause columns to be printed in wrong order in
    `DataFrame.to_string` if specific list of columns passed (GH #325)
  - Fix legend plotting failure if DataFrame columns are integers (GH #326)
  - Shift start date back by one month for Yahoo! Finance API in pandas.io.data
    (GH #329)
  - Fix `DataFrame.join` failure on unconsolidated inputs (GH #331)
  - DataFrame.min/max will no longer fail on mixed-type DataFrame (GH #337)
  - Fix `read_csv` / `read_table` failure when passing list to index_col that is
    not in ascending order (GH #349)
  - Fix failure passing Int64Index to Index.union when both are monotonic
  - Fix error when passing SparseSeries to (dense) DataFrame constructor
  - Added missing bang at top of setup.py (GH #352)
  - Change `is_monotonic` on MultiIndex so it properly compares the tuples
  - Fix MultiIndex outer join logic (GH #351)
  - Set index name attribute with single-key groupby (GH #358)
  - Bug fix in reflexive binary addition in Series and DataFrame for
    non-commutative operations (like string concatenation) (GH #353)
  - setupegg.py will invoke Cython (GH #192)
  - Fix block consolidation bug after inserting column into MultiIndex (GH #366)
  - Fix bug in join operations between Index and Int64Index (GH #367)
  - Handle min_periods=0 case in moving window functions (GH #365)
  - Fixed corner cases in DataFrame.apply/pivot with empty DataFrame (GH #378)
  - Fixed repr exception when Series name is a tuple
  - Always return DateRange from `asfreq` (GH #390)
  - Pass level names to `swaplavel` (GH #379)
  - Don't lose index names in `MultiIndex.droplevel` (GH #394)
  - Infer more proper return type in `DataFrame.apply` when no columns or rows
    depending on whether the passed function is a reduction (GH #389)
  - Always return NA/NaN from Series.min/max and DataFrame.min/max when all of a
    row/column/values are NA (GH #384)
  - Enable partial setting with .ix / advanced indexing (GH #397)
  - Handle mixed-type DataFrames correctly in unstack, do not lose type
    information (GH #403)
  - Fix integer name formatting bug in Index.format and in Series.__repr__
  - Handle label types other than string passed to groupby (GH #405)
  - Fix bug in .ix-based indexing with partial retrieval when a label is not
    contained in a level
  - Index name was not being pickled (GH #408)
  - Level name should be passed to result index in GroupBy.apply (GH #416)

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
    #225)
  - Removed `weights` option in panel regression which was not doing anything
    principled (GH #155)
  - Changed `buffer` argument name in `Series.to_string` to `buf`
  - `Series.to_string` and `DataFrame.to_string` now return strings by default
    instead of printing to sys.stdout
  - Deprecated `nanRep` argument in various `to_string` and `to_csv` functions
    in favor of `na_rep`. Will be removed in 0.6 (GH #275)
  - Renamed `delimiter` to `sep` in `DataFrame.from_csv` for consistency
  - Changed order of `Series.clip` arguments to match those of `numpy.clip` and
    added (unimplemented) `out` argument so `numpy.clip` can be called on a
    Series (GH #272)
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
    lazily iterating through chunks of a flat file (GH #242)
  - Added ability to join on multiple columns in `DataFrame.join` (GH #214)
  - Added private `_get_duplicates` function to `Index` for identifying
    duplicate values more easily
  - Added column attribute access to DataFrame, e.g. df.A equivalent to df['A']
    if 'A' is a column in the DataFrame (PR #213)
  - Added IPython tab completion hook for DataFrame columns. (PR #233, GH #230)
  - Implement `Series.describe` for Series containing objects (PR #241)
  - Add inner join option to `DataFrame.join` when joining on key(s) (GH #248)
  - Can select set of DataFrame columns by passing a list to `__getitem__` (GH
    #253)
  - Can use & and | to intersection / union Index objects, respectively (GH
    #261)
  - Added `pivot_table` convenience function to pandas namespace (GH #234)
  - Implemented `Panel.rename_axis` function (GH #243)
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
    performance (GH #211)
  - Improved speed of `DataFrame.xs` on mixed-type DataFrame objects by about
    5x, regression from 0.3.0 (GH #215)
  - With new `DataFrame.align` method, speeding up binary operations between
    differently-indexed DataFrame objects by 10-25%.
  - Significantly sped up conversion of nested dict into DataFrame (GH #212)
  - Can pass hierarchical index level name to `groupby` instead of the level
    number if desired (GH #223)
  - Add support for different delimiters in `DataFrame.to_csv` (PR #244)
  - Add more helpful error message when importing pandas post-installation from
    the source directory (GH #250)
  - Significantly speed up DataFrame `__repr__` and `count` on large mixed-type
    DataFrame objects
  - Better handling of pyx file dependencies in Cython module build (GH #271)

**Bug fixes**

  - `read_csv` / `read_table` fixes
    - Be less aggressive about converting float->int in cases of floating point
      representations of integers like 1.0, 2.0, etc.
    - "True"/"False" will not get correctly converted to boolean
    - Index name attribute will get set when specifying an index column
    - Passing column names should force `header=None` (GH #257)
    - Don't modify passed column names when `index_col` is not
      None (GH #258)
    - Can sniff CSV separator in zip file (since seek is not supported, was
      failing before)
  - Worked around matplotlib "bug" in which series[:, np.newaxis] fails. Should
    be reported upstream to matplotlib (GH #224)
  - DataFrame.iteritems was not returning Series with the name attribute
    set. Also neither was DataFrame._series
  - Can store datetime.date objects in HDFStore (GH #231)
  - Index and Series names are now stored in HDFStore
  - Fixed problem in which data would get upcasted to object dtype in
    GroupBy.apply operations (GH #237)
  - Fixed outer join bug with empty DataFrame (GH #238)
  - Can create empty Panel (GH #239)
  - Fix join on single key when passing list with 1 entry (GH #246)
  - Don't raise Exception on plotting DataFrame with an all-NA column (GH #251,
    PR #254)
  - Bug min/max errors when called on integer DataFrames (PR #241)
  - `DataFrame.iteritems` and `DataFrame._series` not assigning name attribute
  - Panel.__repr__ raised exception on length-0 major/minor axes
  - `DataFrame.join` on key with empty DataFrame produced incorrect columns
  - Implemented `MultiIndex.diff` (GH #260)
  - `Int64Index.take` and `MultiIndex.take` lost name field, fix downstream
    issue GH #262
  - Can pass list of tuples to `Series` (GH #270)
  - Can pass level name to `DataFrame.stack`
  - Support set operations between MultiIndex and Index
  - Fix many corner cases in MultiIndex set operations
    - Fix MultiIndex-handling bug with GroupBy.apply when returned groups are not
    indexed the same
  - Fix corner case bugs in DataFrame.apply
  - Setting DataFrame index did not cause Series cache to get cleared
  - Various int32 -> int64 platform-specific issues
  - Don't be too aggressive converting to integer when parsing file with
    MultiIndex (GH #285)
  - Fix bug when slicing Series with negative indices before beginning

Thanks
------

- Thomas Kluyver
- Daniel Fortunov
- Aman Thakral
- Luca Beltrame
- Wouter Overmeire

pandas 0.4.3
============

Release notes
-------------

**Release date:** 10/9/2011

This is largely a bugfix release from 0.4.2 but also includes a handful of new
and enhanced features. Also, pandas can now be installed and used on Python 3
(thanks Thomas Kluyver!).

**New features / modules**

  - Python 3 support using 2to3 (PR #200, Thomas Kluyver)
  - Add `name` attribute to `Series` and added relevant logic and tests. Name
    now prints as part of `Series.__repr__`
  - Add `name` attribute to standard Index so that stacking / unstacking does
    not discard names and so that indexed DataFrame objects can be reliably
    round-tripped to flat files, pickle, HDF5, etc.
  - Add `isnull` and `notnull` as instance methods on Series (PR #209, GH #203)

**Improvements to existing features**

  - Skip xlrd-related unit tests if not installed
  - `Index.append` and `MultiIndex.append` can accept a list of Index objects to
    concatenate together
  - Altered binary operations on differently-indexed SparseSeries objects to use
    the integer-based (dense) alignment logic which is faster with a larger
    number of blocks (GH #205)
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
  - `MultiIndex.sortlevel` discarded the level names (GH #202)
  - Fix bugs in groupby, join, and append due to improper concatenation of
    `MultiIndex` objects (GH #201)
  - Fix regression from 0.4.1, `isnull` and `notnull` ceased to work on other
    kinds of Python scalar objects like `datetime.datetime`
  - Raise more helpful exception when attempting to write empty DataFrame or
    LongPanel to `HDFStore` (GH #204)
  - Use stdlib csv module to properly escape strings with commas in
    `DataFrame.to_csv` (PR #206, Thomas Kluyver)
  - Fix Python ndarray access in Cython code for sparse blocked index integrity
    check
  - Fix bug writing Series to CSV in Python 3 (PR #209)
  - Miscellaneous Python 3 bugfixes

Thanks
------

  - Thomas Kluyver
  - rsamson

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
    (GH #187)
  - Wrote templating / code generation script to auto-generate Cython code for
    various functions which need to be available for the 4 major data types
    used in pandas (float64, bool, object, int64)
  - Refactored code related to `DataFrame.join` so that intermediate aligned
    copies of the data in each `DataFrame` argument do not need to be
    created. Substantial performance increases result (GH #176)
  - Substantially improved performance of generic `Index.intersection` and
    `Index.union`
  - Improved performance of `DateRange.union` with overlapping ranges and
    non-cacheable offsets (like Minute). Implemented analogous fast
    `DateRange.intersection` for overlapping ranges.
  - Implemented `BlockManager.take` resulting in significantly faster `take`
    performance on mixed-type `DataFrame` objects (GH #104)
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
  - Throw exception when step specified in label-based slice (GH #185)
  - Fix isnull to correctly work with np.float32. Fix upstream bug described in
    GH #182
  - Finish implementation of as_index=False in groupby for DataFrame
    aggregation (GH #181)
  - Raise SkipTest for pre-epoch HDFStore failure. Real fix will be sorted out
    via datetime64 dtype

Thanks
------

- Uri Laserson
- Scott Sinclair

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
    objects has been implemented (fixes GH #135)
  - `read_csv` can read multiple columns into a `MultiIndex`. DataFrame's
    `to_csv` method will properly write out a `MultiIndex` which can be read
    back (PR #151, thanks to Skipper Seabold)
  - Wrote fast time series merging / joining methods in Cython. Will be
    integrated later into DataFrame.join and related functions
  - Added `ignore_index` option to `DataFrame.append` for combining unindexed
    records stored in a DataFrame

**Improvements to existing features**

  - Some speed enhancements with internal Index type-checking function
  - `DataFrame.rename` has a new `copy` parameter which can rename a DataFrame
    in place
  - Enable unstacking by level name (PR #142)
  - Enable sortlevel to work by level name (PR #141)
  - `read_csv` can automatically "sniff" other kinds of delimiters using
    `csv.Sniffer` (PR #146)
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
  - Fixed single-key groupby on DataFrame with as_index=False (GH #160)
  - `Series.shift` was failing on integer Series (GH #154)
  - `unstack` methods were producing incorrect output in the case of duplicate
    hierarchical labels. An exception will now be raised (GH #147)
  - Calling `count` with level argument caused reduceat failure or segfault in
    earlier NumPy (GH #169)
  - Fixed `DataFrame.corrwith` to automatically exclude non-numeric data (GH
    #144)
  - Unicode handling bug fixes in `DataFrame.to_string` (GH #138)
  - Excluding OLS degenerate unit test case that was causing platform specific
    failure (GH #149)
  - Skip blosc-dependent unit tests for PyTables < 2.2 (PR #137)
  - Calling `copy` on `DateRange` did not copy over attributes to the new object
    (GH #168)
  - Fix bug in `HDFStore` in which Panel data could be appended to a Table with
    different item order, thus resulting in an incorrect result read back

Thanks
------
- Yaroslav Halchenko
- Jeff Reback
- Skipper Seabold
- Dan Lovell
- Nick Pentreath

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
