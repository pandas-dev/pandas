.. _release:

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

pandas 0.13
===========

**Release date:** not-yet-released

**New features**

**Improvements to existing features**

**API Changes**

**Experimental Features**

**Bug Fixes**

pandas 0.12
===========

**Release date:** 2013-07-24

**New features**

  - ``pd.read_html()`` can now parse HTML strings, files or urls and returns a
    list of ``DataFrame`` s courtesy of @cpcloud. (:issue:`3477`,
    :issue:`3605`, :issue:`3606`)
  - Support for reading Amazon S3 files. (:issue:`3504`)
  - Added module for reading and writing JSON strings/files: pandas.io.json
    includes ``to_json`` DataFrame/Series method, and a ``read_json`` top-level reader
    various issues (:issue:`1226`, :issue:`3804`, :issue:`3876`, :issue:`3867`, :issue:`1305`)
  - Added module for reading and writing Stata files: pandas.io.stata (:issue:`1512`)
    includes ``to_stata`` DataFrame method, and a ``read_stata`` top-level reader
  - Added support for writing in ``to_csv`` and reading in ``read_csv``,
    multi-index columns. The ``header`` option in ``read_csv`` now accepts a
    list of the rows from which to read the index. Added the option,
    ``tupleize_cols`` to provide compatiblity for the pre 0.12 behavior of
    writing and reading multi-index columns via a list of tuples. The default in
    0.12 is to write lists of tuples and *not* interpret list of tuples as a
    multi-index column.
    Note: The default value will change in 0.12 to make the default *to* write and
    read multi-index columns in the new format. (:issue:`3571`, :issue:`1651`, :issue:`3141`)
  - Add iterator to ``Series.str`` (:issue:`3638`)
  - ``pd.set_option()`` now allows N option, value pairs (:issue:`3667`).
  - Added keyword parameters for different types of scatter_matrix subplots
  - A ``filter`` method on grouped Series or DataFrames returns a subset of
    the original (:issue:`3680`, :issue:`919`)
  - Access to historical Google Finance data in pandas.io.data (:issue:`3814`)
  - DataFrame plotting methods can sample column colors from a Matplotlib
    colormap via the ``colormap`` keyword. (:issue:`3860`)

**Improvements to existing features**

  - Fixed various issues with internal pprinting code, the repr() for various objects
    including TimeStamp and Index now produces valid python code strings and
    can be used to recreate the object, (:issue:`3038`, :issue:`3379`, :issue:`3251`, :issue:`3460`)
  - ``convert_objects`` now accepts a ``copy`` parameter (defaults to ``True``)
  - ``HDFStore``

    - will retain index attributes (freq,tz,name) on recreation (:issue:`3499`,:issue:`4098`)
    - will warn with a ``AttributeConflictWarning`` if you are attempting to append
      an index with a different frequency than the existing, or attempting
      to append an index with a different name than the existing
    - support datelike columns with a timezone as data_columns (:issue:`2852`)
    - table writing performance improvements.
    - support python3 (via ``PyTables 3.0.0``) (:issue:`3750`)
  - Add modulo operator to Series, DataFrame
  - Add ``date`` method to DatetimeIndex
  - Add ``dropna`` argument to pivot_table (:issue: `3820`)
  - Simplified the API and added a describe method to Categorical
  - ``melt`` now accepts the optional parameters ``var_name`` and ``value_name``
    to specify custom column names of the returned DataFrame (:issue:`3649`),
    thanks @hoechenberger. If ``var_name`` is not specified and ``dataframe.columns.name``
    is not None, then this will be used as the ``var_name`` (:issue:`4144`).
    Also support for MultiIndex columns.
  - clipboard functions use pyperclip (no dependencies on Windows, alternative
    dependencies offered for Linux) (:issue:`3837`).
  - Plotting functions now raise a ``TypeError`` before trying to plot anything
    if the associated objects have have a dtype of ``object`` (:issue:`1818`,
    :issue:`3572`, :issue:`3911`, :issue:`3912`), but they will try to convert object
    arrays to numeric arrays if possible so that you can still plot, for example, an
    object array with floats. This happens before any drawing takes place which
    elimnates any spurious plots from showing up.
  - Added Faq section on repr display options, to help users customize their setup.
  - ``where`` operations that result in block splitting are much faster (:issue:`3733`)
  - Series and DataFrame hist methods now take a ``figsize`` argument (:issue:`3834`)
  - DatetimeIndexes no longer try to convert mixed-integer indexes during join
    operations (:issue:`3877`)
  - Add ``unit`` keyword to ``Timestamp`` and ``to_datetime`` to enable passing of
    integers or floats that are in an epoch unit of ``D, s, ms, us, ns``, thanks @mtkini (:issue:`3969`)
    (e.g. unix timestamps or epoch ``s``, with fracional seconds allowed) (:issue:`3540`)
  - DataFrame corr method (spearman) is now cythonized.
  - Improved ``network`` test decorator to catch ``IOError`` (and therefore
    ``URLError`` as well). Added ``with_connectivity_check`` decorator to allow
    explicitly checking a website as a proxy for seeing if there is network
    connectivity. Plus, new ``optional_args`` decorator factory for decorators.
    (:issue:`3910`, :issue:`3914`)
  - ``read_csv`` will now throw a more informative error message when a file
    contains no columns, e.g., all newline characters
  - Added ``layout`` keyword to DataFrame.hist() for more customizable layout (:issue:`4050`)
  - Timestamp.min and Timestamp.max now represent valid Timestamp instances instead
    of the default datetime.min and datetime.max (respectively), thanks @SleepingPills
  - ``read_html`` now raises when no tables are found and BeautifulSoup==4.2.0
    is detected (:issue:`4214`)

**API Changes**

  - ``HDFStore``

    - When removing an object, ``remove(key)`` raises
      ``KeyError`` if the key is not a valid store object.
    - raise a ``TypeError`` on passing ``where`` or ``columns``
      to select with a Storer; these are invalid parameters at this time (:issue:`4189`)
    - can now specify an ``encoding`` option to ``append/put``
      to enable alternate encodings (:issue:`3750`)
    - enable support for ``iterator/chunksize`` with ``read_hdf``
  - The repr() for (Multi)Index now obeys display.max_seq_items rather
    then numpy threshold print options. (:issue:`3426`, :issue:`3466`)
  - Added mangle_dupe_cols option to read_table/csv, allowing users
    to control legacy behaviour re dupe cols (A, A.1, A.2 vs A, A ) (:issue:`3468`)
    Note: The default value will change in 0.12 to the "no mangle" behaviour,
    If your code relies on this behaviour, explicitly specify mangle_dupe_cols=True
    in your calls.
  - Do not allow astypes on ``datetime64[ns]`` except to ``object``, and
    ``timedelta64[ns]`` to ``object/int`` (:issue:`3425`)
  - The behavior of ``datetime64`` dtypes has changed with respect to certain
    so-called reduction operations (:issue:`3726`). The following operations now
    raise a ``TypeError`` when perfomed on a ``Series`` and return an *empty*
    ``Series`` when performed on a ``DataFrame`` similar to performing these
    operations on, for example, a ``DataFrame`` of ``slice`` objects:
    - sum, prod, mean, std, var, skew, kurt, corr, and cov
  - Do not allow datetimelike/timedeltalike creation except with valid types
    (e.g. cannot pass ``datetime64[ms]``) (:issue:`3423`)
  - Add ``squeeze`` keyword to ``groupby`` to allow reduction from
    DataFrame -> Series if groups are unique. Regression from 0.10.1,
    partial revert on (:issue:`2893`) with (:issue:`3596`)
  - Raise on ``iloc`` when boolean indexing with a label based indexer mask
    e.g. a boolean Series, even with integer labels, will raise. Since ``iloc``
    is purely positional based, the labels on the Series are not alignable (:issue:`3631`)
  - The ``raise_on_error`` option to plotting methods is obviated by :issue:`3572`,
    so it is removed. Plots now always raise when data cannot be plotted or the
    object being plotted has a dtype of ``object``.
  - ``DataFrame.interpolate()`` is now deprecated. Please use
    ``DataFrame.fillna()`` and ``DataFrame.replace()`` instead (:issue:`3582`,
    :issue:`3675`, :issue:`3676`).
  - the ``method`` and ``axis`` arguments of ``DataFrame.replace()`` are
    deprecated
  - ``DataFrame.replace`` 's ``infer_types`` parameter is removed and now
    performs conversion by default. (:issue:`3907`)
  - Deprecated display.height, display.width is now only a formatting option
    does not control triggering of summary, similar to < 0.11.0.
  - Add the keyword ``allow_duplicates`` to ``DataFrame.insert`` to allow a duplicate column
    to be inserted if ``True``, default is ``False`` (same as prior to 0.12) (:issue:`3679`)
  - io API changes

    - added ``pandas.io.api`` for i/o imports
    - removed ``Excel`` support to ``pandas.io.excel``
    - added top-level ``pd.read_sql`` and ``to_sql`` DataFrame methods
    - removed ``clipboard`` support to ``pandas.io.clipboard``
    - replace top-level and instance methods ``save`` and ``load`` with
      top-level ``read_pickle`` and ``to_pickle`` instance method, ``save`` and
      ``load`` will give deprecation warning.
  - the ``method`` and ``axis`` arguments of ``DataFrame.replace()`` are
    deprecated
  - set FutureWarning to require data_source, and to replace year/month with
    expiry date in pandas.io options. This is in preparation to add options
    data from google (:issue:`3822`)
  - the ``method`` and ``axis`` arguments of ``DataFrame.replace()`` are
    deprecated
  - Implement ``__nonzero__`` for ``NDFrame`` objects (:issue:`3691`, :issue:`3696`)
  - ``as_matrix`` with mixed signed and unsigned dtypes will result in 2 x the lcd of the unsigned
    as an int, maxing with ``int64``, to avoid precision issues (:issue:`3733`)
  - ``na_values`` in a list provided to ``read_csv/read_excel`` will match string and numeric versions
    e.g. ``na_values=['99']`` will match 99 whether the column ends up being int, float, or string (:issue:`3611`)
  - ``read_html`` now defaults to ``None`` when reading, and falls back on
    ``bs4`` + ``html5lib`` when lxml fails to parse. a list of parsers to try
    until success is also valid
  - more consistency in the to_datetime return types (give string/array of string inputs) (:issue:`3888`)
  - The internal ``pandas`` class hierarchy has changed (slightly). The
    previous ``PandasObject`` now is called ``PandasContainer`` and a new
    ``PandasObject`` has become the baseclass for ``PandasContainer`` as well
    as ``Index``, ``Categorical``, ``GroupBy``, ``SparseList``, and
    ``SparseArray`` (+ their base classes). Currently, ``PandasObject``
    provides string methods (from ``StringMixin``). (:issue:`4090`, :issue:`4092`)
  - New ``StringMixin`` that, given a ``__unicode__`` method, gets python 2 and
    python 3 compatible string methods (``__str__``, ``__bytes__``, and
    ``__repr__``). Plus string safety throughout. Now employed in many places
    throughout the pandas library. (:issue:`4090`, :issue:`4092`)

**Experimental Features**

  - Added experimental ``CustomBusinessDay`` class to support ``DateOffsets``
    with custom holiday calendars and custom weekmasks. (:issue:`2301`)

**Bug Fixes**

  - Fixed an esoteric excel reading bug, xlrd>= 0.9.0 now required for excel
    support. Should provide python3 support (for reading) which has been
    lacking. (:issue:`3164`)
  - Disallow Series constructor called with MultiIndex which caused segfault (:issue:`4187`)
  - Allow unioning of date ranges sharing a timezone (:issue:`3491`)
  - Fix to_csv issue when having a large number of rows and ``NaT`` in some
    columns (:issue:`3437`)
  - ``.loc`` was not raising when passed an integer list (:issue:`3449`)
  - Unordered time series selection was misbehaving when using label slicing (:issue:`3448`)
  - Fix sorting in a frame with a list of columns which contains datetime64[ns] dtypes (:issue:`3461`)
  - DataFrames fetched via FRED now handle '.' as a NaN. (:issue:`3469`)
  - Fix regression in a DataFrame apply with axis=1, objects were not being converted back
    to base dtypes correctly (:issue:`3480`)
  - Fix issue when storing uint dtypes in an HDFStore. (:issue:`3493`)
  - Non-unique index support clarified (:issue:`3468`)

    - Addressed handling of dupe columns in df.to_csv new and old (:issue:`3454`, :issue:`3457`)
    - Fix assigning a new index to a duplicate index in a DataFrame would fail (:issue:`3468`)
    - Fix construction of a DataFrame with a duplicate index
    - ref_locs support to allow duplicative indices across dtypes,
      allows iget support to always find the index (even across dtypes) (:issue:`2194`)
    - applymap on a DataFrame with a non-unique index now works
      (removed warning) (:issue:`2786`), and fix (:issue:`3230`)
    - Fix to_csv to handle non-unique columns (:issue:`3495`)
    - Duplicate indexes with getitem will return items in the correct order (:issue:`3455`, :issue:`3457`)
      and handle missing elements like unique indices (:issue:`3561`)
    - Duplicate indexes with and empty DataFrame.from_records will return a correct frame (:issue:`3562`)
    - Concat to produce a non-unique columns when duplicates are across dtypes is fixed (:issue:`3602`)
    - Non-unique indexing with a slice via ``loc`` and friends fixed (:issue:`3659`)
    - Allow insert/delete to non-unique columns (:issue:`3679`)
    - Extend ``reindex`` to correctly deal with non-unique indices (:issue:`3679`)
    - ``DataFrame.itertuples()`` now works with frames with duplicate column
      names (:issue:`3873`)
    - Bug in non-unique indexing via ``iloc`` (:issue:`4017`); added ``takeable`` argument to
      ``reindex`` for location-based taking
    - Allow non-unique indexing in series via ``.ix/.loc`` and ``__getitem__`` (:issue:`4246`)
    - Fixed non-unique indexing memory allocation issue with ``.ix/.loc`` (:issue:`4280`)

  - Fixed bug in groupby with empty series referencing a variable before assignment. (:issue:`3510`)
  - Allow index name to be used in groupby for non MultiIndex (:issue:`4014`)
  - Fixed bug in mixed-frame assignment with aligned series (:issue:`3492`)
  - Fixed bug in selecting month/quarter/year from a series would not select the time element
    on the last day (:issue:`3546`)
  - Fixed a couple of MultiIndex rendering bugs in df.to_html() (:issue:`3547`, :issue:`3553`)
  - Properly convert np.datetime64 objects in a Series (:issue:`3416`)
  - Raise a ``TypeError`` on invalid datetime/timedelta operations
    e.g. add datetimes, multiple timedelta x datetime
  - Fix ``.diff`` on datelike and timedelta operations (:issue:`3100`)
  - ``combine_first`` not returning the same dtype in cases where it can (:issue:`3552`)
  - Fixed bug with ``Panel.transpose`` argument aliases (:issue:`3556`)
  - Fixed platform bug in ``PeriodIndex.take`` (:issue:`3579`)
  - Fixed bud in incorrect conversion of datetime64[ns] in ``combine_first`` (:issue:`3593`)
  - Fixed bug in reset_index with ``NaN`` in a multi-index (:issue:`3586`)
  - ``fillna`` methods now raise a ``TypeError`` when the ``value`` parameter
    is a ``list`` or ``tuple``.
  - Fixed bug where a time-series was being selected in preference to an actual column name
    in a frame (:issue:`3594`)
  - Make secondary_y work properly for bar plots (:issue:`3598`)
  - Fix modulo and integer division on Series,DataFrames to act similary to ``float`` dtypes to return
    ``np.nan`` or ``np.inf`` as appropriate (:issue:`3590`)
  - Fix incorrect dtype on groupby with ``as_index=False`` (:issue:`3610`)
  - Fix ``read_csv/read_excel`` to correctly encode identical na_values, e.g. ``na_values=[-999.0,-999]``
    was failing (:issue:`3611`)
  - Disable HTML output in qtconsole again. (:issue:`3657`)
  - Reworked the new repr display logic, which users found confusing. (:issue:`3663`)
  - Fix indexing issue in ndim >= 3 with ``iloc`` (:issue:`3617`)
  - Correctly parse date columns with embedded (nan/NaT) into datetime64[ns] dtype in ``read_csv``
    when ``parse_dates`` is specified (:issue:`3062`)
  - Fix not consolidating before to_csv (:issue:`3624`)
  - Fix alignment issue when setitem in a DataFrame with a piece of a DataFrame (:issue:`3626`) or
    a mixed DataFrame and a Series (:issue:`3668`)
  - Fix plotting of unordered DatetimeIndex (:issue:`3601`)
  - ``sql.write_frame`` failing when writing a single column to sqlite (:issue:`3628`),
    thanks to @stonebig
  - Fix pivoting with ``nan`` in the index (:issue:`3558`)
  - Fix running of bs4 tests when it is not installed (:issue:`3605`)
  - Fix parsing of html table (:issue:`3606`)
  - ``read_html()`` now only allows a single backend: ``html5lib`` (:issue:`3616`)
  - ``convert_objects`` with ``convert_dates='coerce'`` was parsing some single-letter strings into today's date
  - ``DataFrame.from_records`` did not accept empty recarrays (:issue:`3682`)
  - ``DataFrame.to_csv`` will succeed with the deprecated option ``nanRep``, @tdsmith
  - ``DataFrame.to_html`` and ``DataFrame.to_latex`` now accept a path for
    their first argument (:issue:`3702`)
  - Fix file tokenization error with \r delimiter and quoted fields (:issue:`3453`)
  - Groupby transform with item-by-item not upcasting correctly (:issue:`3740`)
  - Incorrectly read a HDFStore multi-index Frame witha column specification (:issue:`3748`)
  - ``read_html`` now correctly skips tests (:issue:`3741`)
  - PandasObjects raise TypeError when trying to hash (:issue:`3882`)
  - Fix incorrect arguments passed to concat that are not list-like (e.g. concat(df1,df2)) (:issue:`3481`)
  - Correctly parse when passed the ``dtype=str`` (or other variable-len string dtypes)
    in ``read_csv`` (:issue:`3795`)
  - Fix index name not propogating when using ``loc/ix`` (:issue:`3880`)
  - Fix groupby when applying a custom function resulting in a returned DataFrame was
    not converting dtypes (:issue:`3911`)
  - Fixed a bug where ``DataFrame.replace`` with a compiled regular expression
    in the ``to_replace`` argument wasn't working (:issue:`3907`)
  - Fixed ``__truediv__`` in Python 2.7 with ``numexpr`` installed to actually do true division when dividing
    two integer arrays with at least 10000 cells total (:issue:`3764`)
  - Indexing with a string with seconds resolution not selecting from a time index (:issue:`3925`)
  - csv parsers would loop infinitely if ``iterator=True`` but no ``chunksize`` was
    specified (:issue:`3967`), python parser failing with ``chunksize=1``
  - Fix index name not propogating when using ``shift``
  - Fixed dropna=False being ignored with multi-index stack (:issue:`3997`)
  - Fixed flattening of columns when renaming MultiIndex columns DataFrame (:issue:`4004`)
  - Fix ``Series.clip`` for datetime series. NA/NaN threshold values will now throw ValueError (:issue:`3996`)
  - Fixed insertion issue into DataFrame, after rename (:issue:`4032`)
  - Fixed testing issue where too many sockets where open thus leading to a
    connection reset issue (:issue:`3982`, :issue:`3985`, :issue:`4028`,
    :issue:`4054`)
  - Fixed failing tests in test_yahoo, test_google where symbols were not
    retrieved but were being accessed (:issue:`3982`, :issue:`3985`,
    :issue:`4028`, :issue:`4054`)
  - ``Series.hist`` will now take the figure from the current environment if
    one is not passed
  - Fixed bug where a 1xN DataFrame would barf on a 1xN mask (:issue:`4071`)
  - Fixed running of ``tox`` under python3 where the pickle import was getting
    rewritten in an incompatible way (:issue:`4062`, :issue:`4063`)
  - Fixed bug where sharex and sharey were not being passed to grouped_hist
    (:issue:`4089`)
  - Fix bug where ``HDFStore`` will fail to append because of a different block
    ordering on-disk (:issue:`4096`)
  - Better error messages on inserting incompatible columns to a frame (:issue:`4107`)
  - Fixed bug in ``DataFrame.replace`` where a nested dict wasn't being
    iterated over when regex=False (:issue:`4115`)
  - Fixed bug in ``convert_objects(convert_numeric=True)`` where a mixed numeric and
    object Series/Frame was not converting properly (:issue:`4119`)
  - Fixed bugs in multi-index selection with column multi-index and duplicates
    (:issue:`4145`, :issue:`4146`)
  - Fixed bug in the parsing of microseconds when using the ``format``
    argument in ``to_datetime`` (:issue:`4152`)
  - Fixed bug in ``PandasAutoDateLocator`` where ``invert_xaxis`` triggered
    incorrectly ``MilliSecondLocator``  (:issue:`3990`)
  - Fixed bug in ``Series.where`` where broadcasting a single element input vector
    to the length of the series resulted in multiplying the value
    inside the input (:issue:`4192`)
  - Fixed bug in plotting that wasn't raising on invalid colormap for
    matplotlib 1.1.1 (:issue:`4215`)
  - Fixed the legend displaying in ``DataFrame.plot(kind='kde')`` (:issue:`4216`)
  - Fixed bug where Index slices weren't carrying the name attribute
    (:issue:`4226`)
  - Fixed bug in initializing ``DatetimeIndex`` with an array of strings
    in a certain time zone (:issue:`4229`)
  - Fixed bug where html5lib wasn't being properly skipped (:issue:`4265`)
  - Fixed bug where get_data_famafrench wasn't using the correct file edges
    (:issue:`4281`)

pandas 0.11.0
=============

**Release date:** 2013-04-22

**New features**

  - New documentation section, ``10 Minutes to Pandas``
  - New documentation section, ``Cookbook``
  - Allow mixed dtypes (e.g ``float32/float64/int32/int16/int8``) to coexist in
    DataFrames and propogate in operations
  - Add function to pandas.io.data for retrieving stock index components from
    Yahoo! finance (:issue:`2795`)
  - Support slicing with time objects (:issue:`2681`)
  - Added ``.iloc`` attribute, to support strict integer based indexing,
    analogous to ``.ix`` (:issue:`2922`)
  - Added ``.loc`` attribute, to support strict label based indexing, analagous
    to ``.ix`` (:issue:`3053`)
  - Added ``.iat`` attribute, to support fast scalar access via integers
    (replaces ``iget_value/iset_value``)
  - Added ``.at`` attribute, to support fast scalar access via labels (replaces
    ``get_value/set_value``)
  - Moved functionaility from ``irow,icol,iget_value/iset_value`` to ``.iloc`` indexer
    (via ``_ixs`` methods in each object)
  - Added support for expression evaluation using the ``numexpr`` library
  - Added ``convert=boolean`` to ``take`` routines to translate negative
    indices to positive, defaults to True
  - Added to_series() method to indices, to facilitate the creation of indexeres
    (:issue:`3275`)

**Improvements to existing features**

  - Improved performance of df.to_csv() by up to 10x in some cases. (:issue:`3059`)
  - added ``blocks`` attribute to DataFrames, to return a dict of dtypes to
    homogeneously dtyped DataFrames
  - added keyword ``convert_numeric`` to ``convert_objects()`` to try to
    convert object dtypes to numeric types (default is False)
  - ``convert_dates`` in ``convert_objects`` can now be ``coerce`` which will
    return a datetime64[ns] dtype with non-convertibles set as ``NaT``; will
    preserve an all-nan object (e.g. strings), default is True (to perform
    soft-conversion
  - Series print output now includes the dtype by default
  - Optimize internal reindexing routines (:issue:`2819`, :issue:`2867`)
  - ``describe_option()`` now reports the default and current value of options.
  - Add ``format`` option to ``pandas.to_datetime`` with faster conversion of
    strings that can be parsed with datetime.strptime
  - Add ``axes`` property to ``Series`` for compatibility
  - Add ``xs`` function to ``Series`` for compatibility
  - Allow setitem in a frame where only mixed numerics are present (e.g. int
    and float), (:issue:`3037`)
  - ``HDFStore``

    - Provide dotted attribute access to ``get`` from stores
      (e.g. store.df == store['df'])
    - New keywords ``iterator=boolean``, and ``chunksize=number_in_a_chunk``
      are provided to support iteration on ``select`` and
      ``select_as_multiple`` (:issue:`3076`)
    - support ``read_hdf/to_hdf`` API similar to ``read_csv/to_csv`` (:issue:`3222`)

  - Add ``squeeze`` method to possibly remove length 1 dimensions from an
    object.

    .. ipython:: python

       p = Panel(randn(3,4,4),items=['ItemA','ItemB','ItemC'],
                          major_axis=date_range('20010102',periods=4),
                          minor_axis=['A','B','C','D'])
       p
       p.reindex(items=['ItemA']).squeeze()
       p.reindex(items=['ItemA'],minor=['B']).squeeze()

  - Improvement to Yahoo API access in ``pd.io.data.Options`` (:issue:`2758`)
  - added option `display.max_seq_items` to control the number of
    elements printed per sequence pprinting it. (:issue:`2979`)
  - added option `display.chop_threshold` to control display of small numerical
    values. (:issue:`2739`)
  - added option `display.max_info_rows` to prevent verbose_info from being
    calculated for frames above 1M rows (configurable). (:issue:`2807`, :issue:`2918`)
  - value_counts() now accepts a "normalize" argument, for normalized
    histograms. (:issue:`2710`).
  - DataFrame.from_records now accepts not only dicts but any instance of
    the collections.Mapping ABC.
  - Allow selection semantics via a string with a datelike index to work in both
    Series and DataFrames (:issue:`3070`)

    .. ipython:: python

        idx = date_range("2001-10-1", periods=5, freq='M')
        ts = Series(np.random.rand(len(idx)),index=idx)
        ts['2001']

        df = DataFrame(dict(A = ts))
        df['2001']

  - added option `display.mpl_style` providing a sleeker visual style
    for plots. Based on https://gist.github.com/huyng/816622 (:issue:`3075`).


  - Improved performance across several core functions by taking memory
    ordering of arrays into account. Courtesy of @stephenwlin (:issue:`3130`)
  - Improved performance of groupby transform method (:issue:`2121`)
  - Handle "ragged" CSV files missing trailing delimiters in rows with missing
    fields when also providing explicit list of column names (so the parser
    knows how many columns to expect in the result) (:issue:`2981`)
  - On a mixed DataFrame, allow setting with indexers with ndarray/DataFrame
    on rhs (:issue:`3216`)
  - Treat boolean values as integers (values 1 and 0) for numeric
    operations. (:issue:`2641`)
  - Add ``time`` method to DatetimeIndex (:issue:`3180`)
  - Return NA when using Series.str[...] for values that are not long enough
    (:issue:`3223`)
  - Display cursor coordinate information in time-series plots (:issue:`1670`)
  - to_html() now accepts an optional "escape" argument to control reserved
    HTML character escaping (enabled by default) and escapes ``&``, in addition
    to ``<`` and ``>``.  (:issue:`2919`)

**API Changes**

  - Do not automatically upcast numeric specified dtypes to ``int64`` or
    ``float64`` (:issue:`622` and :issue:`797`)
  - DataFrame construction of lists and scalars, with no dtype present, will
    result in casting to ``int64`` or ``float64``, regardless of platform.
    This is not an apparent change in the API, but noting it.
  - Guarantee that ``convert_objects()`` for Series/DataFrame always returns a
    copy
  - groupby operations will respect dtypes for numeric float operations
    (float32/float64); other types will be operated on, and will try to cast
    back to the input dtype (e.g. if an int is passed, as long as the output
    doesn't have nans, then an int will be returned)
  - backfill/pad/take/diff/ohlc will now support ``float32/int16/int8``
    operations
  - Block types will upcast as needed in where/masking operations (:issue:`2793`)
  - Series now automatically will try to set the correct dtype based on passed
    datetimelike objects (datetime/Timestamp)

    - timedelta64 are returned in appropriate cases (e.g. Series - Series,
      when both are datetime64)
    - mixed datetimes and objects (:issue:`2751`) in a constructor will be cast
      correctly
    - astype on datetimes to object are now handled (as well as NaT
      conversions to np.nan)
    - all timedelta like objects will be correctly assigned to ``timedelta64``
      with mixed ``NaN`` and/or ``NaT`` allowed

  - arguments to DataFrame.clip were inconsistent to numpy and Series clipping
    (:issue:`2747`)
  - util.testing.assert_frame_equal now checks the column and index names (:issue:`2964`)
  - Constructors will now return a more informative ValueError on failures
    when invalid shapes are passed
  - Don't suppress TypeError in GroupBy.agg (:issue:`3238`)
  - Methods return None when inplace=True (:issue:`1893`)
  - ``HDFStore``

     - added the method ``select_column`` to select a single column from a table as a Series.
     - deprecated the ``unique`` method, can be replicated by ``select_column(key,column).unique()``
     - ``min_itemsize`` parameter will now automatically create data_columns for passed keys

  - Downcast on pivot if possible (:issue:`3283`), adds argument ``downcast`` to ``fillna``
  - Introduced options `display.height/width` for explicitly specifying terminal
    height/width in characters. Deprecated display.line_width, now replaced by display.width.
    These defaults are in effect for scripts as well, so unless disabled, previously
    very wide output will now be output as "expand_repr" style wrapped output.
  - Various defaults for options (including display.max_rows) have been revised,
    after a brief survey concluded they were wrong for everyone. Now at w=80,h=60.
  - HTML repr output in IPython qtconsole is once again controlled by the option
    `display.notebook_repr_html`, and on by default.

**Bug Fixes**

  - Fix seg fault on empty data frame when fillna with ``pad`` or ``backfill``
    (:issue:`2778`)
  - Single element ndarrays of datetimelike objects are handled
    (e.g. np.array(datetime(2001,1,1,0,0))), w/o dtype being passed
  - 0-dim ndarrays with a passed dtype are handled correctly
    (e.g. np.array(0.,dtype='float32'))
  - Fix some boolean indexing inconsistencies in Series.__getitem__/__setitem__
    (:issue:`2776`)
  - Fix issues with DataFrame and Series constructor with integers that
    overflow ``int64`` and some mixed typed type lists (:issue:`2845`)

  - ``HDFStore``

    - Fix weird PyTables error when using too many selectors in a where
      also correctly filter on any number of values in a Term expression
      (so not using numexpr filtering, but isin filtering)
    - Internally, change all variables to be private-like (now have leading
      underscore)
    - Fixes for query parsing to correctly interpret boolean and != (:issue:`2849`, :issue:`2973`)
    - Fixes for pathological case on SparseSeries with 0-len array and
      compression (:issue:`2931`)
    - Fixes bug with writing rows if part of a block was all-nan (:issue:`3012`)
    - Exceptions are now ValueError or TypeError as needed
    - A table will now raise if min_itemsize contains fields which are not queryables

  - Bug showing up in applymap where some object type columns are converted (:issue:`2909`)
    had an incorrect default in convert_objects

  - TimeDeltas

    - Series ops with a Timestamp on the rhs was throwing an exception (:issue:`2898`)
      added tests for Series ops with datetimes,timedeltas,Timestamps, and datelike
      Series on both lhs and rhs
    - Fixed subtle timedelta64 inference issue on py3 & numpy 1.7.0 (:issue:`3094`)
    - Fixed some formatting issues on timedelta when negative
    - Support null checking on timedelta64, representing (and formatting) with NaT
    - Support setitem with np.nan value, converts to NaT
    - Support min/max ops in a Dataframe (abs not working, nor do we error on non-supported ops)
    - Support idxmin/idxmax/abs/max/min in a Series (:issue:`2989`, :issue:`2982`)

  - Bug on in-place putmasking on an ``integer`` series that needs to be converted to
    ``float`` (:issue:`2746`)
  - Bug in argsort of ``datetime64[ns]`` Series with ``NaT`` (:issue:`2967`)
  - Bug in value_counts of ``datetime64[ns]`` Series (:issue:`3002`)
  - Fixed printing of ``NaT`` in an index
  - Bug in idxmin/idxmax of ``datetime64[ns]`` Series with ``NaT`` (:issue:`2982`)
  - Bug in ``icol, take`` with negative indicies was producing incorrect return
    values (see :issue:`2922`, :issue:`2892`), also check for out-of-bounds indices (:issue:`3029`)
  - Bug in DataFrame column insertion when the column creation fails, existing frame is left in
    an irrecoverable state (:issue:`3010`)
  - Bug in DataFrame update, combine_first where non-specified values could cause
    dtype changes (:issue:`3016`, :issue:`3041`)
  - Bug in groupby with first/last where dtypes could change (:issue:`3041`, :issue:`2763`)
  - Formatting of an index that has ``nan`` was inconsistent or wrong (would fill from
    other values), (:issue:`2850`)
  - Unstack of a frame with no nans would always cause dtype upcasting (:issue:`2929`)
  - Fix scalar datetime.datetime parsing bug in read_csv (:issue:`3071`)
  - Fixed slow printing of large Dataframes, due to inefficient dtype
    reporting (:issue:`2807`)
  - Fixed a segfault when using a function as grouper in groupby (:issue:`3035`)
  - Fix pretty-printing of infinite data structures (closes :issue:`2978`)
  - Fixed exception when plotting timeseries bearing a timezone (closes :issue:`2877`)
  - str.contains ignored na argument (:issue:`2806`)
  - Substitute warning for segfault when grouping with categorical grouper
    of mismatched length (:issue:`3011`)
  - Fix exception in SparseSeries.density (:issue:`2083`)
  - Fix upsampling bug with closed='left' and daily to daily data (:issue:`3020`)
  - Fixed missing tick bars on scatter_matrix plot (:issue:`3063`)
  - Fixed bug in Timestamp(d,tz=foo) when d is date() rather then datetime() (:issue:`2993`)
  - series.plot(kind='bar') now respects pylab color schem (:issue:`3115`)
  - Fixed bug in reshape if not passed correct input, now raises TypeError (:issue:`2719`)
  - Fixed a bug where Series ctor did not respect ordering if OrderedDict passed in (:issue:`3282`)
  - Fix NameError issue on RESO_US (:issue:`2787`)
  - Allow selection in an *unordered* timeseries to work similary
    to an *ordered* timeseries (:issue:`2437`).
  - Fix implemented ``.xs`` when called with ``axes=1`` and a level parameter (:issue:`2903`)
  - Timestamp now supports the class method fromordinal similar to datetimes (:issue:`3042`)
  - Fix issue with indexing a series with a boolean key and specifiying a 1-len list on the rhs (:issue:`2745`)
    or a list on the rhs (:issue:`3235`)
  - Fixed bug in groupby apply when kernel generate list of arrays having unequal len (:issue:`1738`)
  - fixed handling of rolling_corr with center=True which could produce corr>1 (:issue:`3155`)
  - Fixed issues where indices can be passed as 'index/column' in addition to 0/1 for the axis parameter
  - PeriodIndex.tolist now boxes to Period (:issue:`3178`)
  - PeriodIndex.get_loc KeyError now reports Period instead of ordinal (:issue:`3179`)
  - df.to_records bug when handling MultiIndex (GH3189)
  - Fix Series.__getitem__ segfault when index less than -length (:issue:`3168`)
  - Fix bug when using Timestamp as a date parser (:issue:`2932`)
  - Fix bug creating date range from Timestamp with time zone and passing same
    time zone (:issue:`2926`)
  - Add comparison operators to Period object (:issue:`2781`)
  - Fix bug when concatenating two Series into a DataFrame when they have the
    same name (:issue:`2797`)
  - Fix automatic color cycling when plotting consecutive timeseries
    without color arguments (:issue:`2816`)
  - fixed bug in the pickling of PeriodIndex (:issue:`2891`)
  - Upcast/split blocks when needed in a mixed DataFrame when setitem
    with an indexer (:issue:`3216`)
  - Invoking df.applymap on a dataframe with dupe cols now raises a ValueError (:issue:`2786`)
  - Apply with invalid returned indices raise correct Exception (:issue:`2808`)
  - Fixed a bug in plotting log-scale bar plots (:issue:`3247`)
  - df.plot() grid on/off now obeys the mpl default style, just like
    series.plot(). (:issue:`3233`)
  - Fixed a bug in the legend of plotting.andrews_curves() (:issue:`3278`)
  - Produce a series on apply if we only generate a singular series and have
    a simple index (:issue:`2893`)
  - Fix Python ascii file parsing when integer falls outside of floating point
    spacing (:issue:`3258`)
  - fixed pretty priniting of sets (:issue:`3294`)
  - Panel() and Panel.from_dict() now respects ordering when give OrderedDict (:issue:`3303`)
  - DataFrame where with a datetimelike incorrectly selecting (:issue:`3311`)
  - Ensure index casts work even in Int64Index
  - Fix set_index segfault when passing MultiIndex (:issue:`3308`)
  - Ensure pickles created in py2 can be read in py3
  - Insert ellipsis in MultiIndex summary repr (:issue:`3348`)
  - Groupby will handle mutation among an input groups columns (and fallback
    to non-fast apply) (:issue:`3380`)
  - Eliminated unicode errors on FreeBSD when using MPL GTK backend (:issue:`3360`)
  - Period.strftime should return unicode strings always (:issue:`3363`)
  - Respect passed read_* chunksize in get_chunk function (:issue:`3406`)


pandas 0.10.1
=============

**Release date:** 2013-01-22

**New features**

  - Add data inferface to World Bank WDI pandas.io.wb (:issue:`2592`)

**API Changes**

  - Restored inplace=True behavior returning self (same object) with
    deprecation warning until 0.11 (:issue:`1893`)
  - ``HDFStore``

    - refactored HFDStore to deal with non-table stores as objects, will allow future enhancements
    - removed keyword ``compression`` from ``put`` (replaced by keyword
      ``complib`` to be consistent across library)
    - warn `PerformanceWarning` if you are attempting to store types that will be pickled by PyTables

**Improvements to existing features**

  - ``HDFStore``

    - enables storing of multi-index dataframes (closes :issue:`1277`)
    - support data column indexing and selection, via ``data_columns`` keyword
      in append
    - support write chunking to reduce memory footprint, via ``chunksize``
      keyword to append
    - support automagic indexing via ``index`` keyword to append
    - support ``expectedrows`` keyword in append to inform ``PyTables`` about
      the expected tablesize
    - support ``start`` and ``stop`` keywords in select to limit the row
      selection space
    - added ``get_store`` context manager to automatically import with pandas
    - added column filtering via ``columns`` keyword in select
    - added methods append_to_multiple/select_as_multiple/select_as_coordinates
      to do multiple-table append/selection
    - added support for datetime64 in columns
    - added method ``unique`` to select the unique values in an indexable or
      data column
    - added method ``copy`` to copy an existing store (and possibly upgrade)
    - show the shape of the data on disk for non-table stores when printing the
      store
    - added ability to read PyTables flavor tables (allows compatiblity to
      other HDF5 systems)

  - Add ``logx`` option to DataFrame/Series.plot (:issue:`2327`, :issue:`2565`)
  - Support reading gzipped data from file-like object
  - ``pivot_table`` aggfunc can be anything used in GroupBy.aggregate (:issue:`2643`)
  - Implement DataFrame merges in case where set cardinalities might overflow
    64-bit integer (:issue:`2690`)
  - Raise exception in C file parser if integer dtype specified and have NA
    values. (:issue:`2631`)
  - Attempt to parse ISO8601 format dates when parse_dates=True in read_csv for
    major performance boost in such cases (:issue:`2698`)
  - Add methods ``neg`` and ``inv`` to Series
  - Implement ``kind`` option in ``ExcelFile`` to indicate whether it's an XLS
    or XLSX file (:issue:`2613`)

**Bug fixes**

  - Fix read_csv/read_table multithreading issues (:issue:`2608`)
  - ``HDFStore``

    - correctly handle ``nan`` elements in string columns; serialize via the
      ``nan_rep`` keyword to append
    - raise correctly on non-implemented column types (unicode/date)
    - handle correctly ``Term`` passed types (e.g. ``index<1000``, when index
      is ``Int64``), (closes :issue:`512`)
    - handle Timestamp correctly in data_columns (closes :issue:`2637`)
    - contains correctly matches on non-natural names
    - correctly store ``float32`` dtypes in tables (if not other float types in
      the same table)

  - Fix DataFrame.info bug with UTF8-encoded columns. (:issue:`2576`)
  - Fix DatetimeIndex handling of FixedOffset tz (:issue:`2604`)
  - More robust detection of being in IPython session for wide DataFrame
    console formatting (:issue:`2585`)
  - Fix platform issues with ``file:///`` in unit test (:issue:`2564`)
  - Fix bug and possible segfault when grouping by hierarchical level that
    contains NA values (:issue:`2616`)
  - Ensure that MultiIndex tuples can be constructed with NAs (:issue:`2616`)
  - Fix int64 overflow issue when unstacking MultiIndex with many levels
    (:issue:`2616`)
  - Exclude non-numeric data from DataFrame.quantile by default (:issue:`2625`)
  - Fix a Cython C int64 boxing issue causing read_csv to return incorrect
    results (:issue:`2599`)
  - Fix groupby summing performance issue on boolean data (:issue:`2692`)
  - Don't bork Series containing datetime64 values with to_datetime (:issue:`2699`)
  - Fix DataFrame.from_records corner case when passed columns, index column,
    but empty record list (:issue:`2633`)
  - Fix C parser-tokenizer bug with trailing fields. (:issue:`2668`)
  - Don't exclude non-numeric data from GroupBy.max/min (:issue:`2700`)
  - Don't lose time zone when calling DatetimeIndex.drop (:issue:`2621`)
  - Fix setitem on a Series with a boolean key and a non-scalar as value
    (:issue:`2686`)
  - Box datetime64 values in Series.apply/map (:issue:`2627`, :issue:`2689`)
  - Upconvert datetime + datetime64 values when concatenating frames (:issue:`2624`)
  - Raise a more helpful error message in merge operations when one DataFrame
    has duplicate columns (:issue:`2649`)
  - Fix partial date parsing issue occuring only when code is run at EOM
    (:issue:`2618`)
  - Prevent MemoryError when using counting sort in sortlevel with
    high-cardinality MultiIndex objects (:issue:`2684`)
  - Fix Period resampling bug when all values fall into a single bin (:issue:`2070`)
  - Fix buggy interaction with usecols argument in read_csv when there is an
    implicit first index column (:issue:`2654`)


pandas 0.10.0
=============

**Release date:** 2012-12-17

**New features**

  - Brand new high-performance delimited file parsing engine written in C and
    Cython. 50% or better performance in many standard use cases with a
    fraction as much memory usage. (:issue:`407`, :issue:`821`)
  - Many new file parser (read_csv, read_table) features:

    - Support for on-the-fly gzip or bz2 decompression (`compression` option)
    - Ability to get back numpy.recarray instead of DataFrame
      (`as_recarray=True`)
    - `dtype` option: explicit column dtypes
    - `usecols` option: specify list of columns to be read from a file. Good
      for reading very wide files with many irrelevant columns (:issue:`1216` :issue:`926`, :issue:`2465`)
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
    - Easy of European (and other) decimal formats (`decimal` option) (:issue:`584`, :issue:`2466`)
    - Custom line terminators (e.g. lineterminator='~') (:issue:`2457`)
    - Handling of no trailing commas in CSV files (:issue:`2333`)
    - Ability to handle fractional seconds in date_converters (:issue:`2209`)
    - read_csv allow scalar arg to na_values (:issue:`1944`)
    - Explicit column dtype specification in read_* functions (:issue:`1858`)
    - Easier CSV dialect specification (:issue:`1743`)
    - Improve parser performance when handling special characters (:issue:`1204`)

  - Google Analytics API integration with easy oauth2 workflow (:issue:`2283`)
  - Add error handling to Series.str.encode/decode (:issue:`2276`)
  - Add ``where`` and ``mask`` to Series (:issue:`2337`)
  - Grouped histogram via `by` keyword in Series/DataFrame.hist (:issue:`2186`)
  - Support optional ``min_periods`` keyword in ``corr`` and ``cov``
    for both Series and DataFrame (:issue:`2002`)
  - Add ``duplicated`` and ``drop_duplicates`` functions to Series (:issue:`1923`)
  - Add docs for ``HDFStore table`` format
  - 'density' property in `SparseSeries` (:issue:`2384`)
  - Add ``ffill`` and ``bfill`` convenience functions for forward- and
    backfilling time series data (:issue:`2284`)
  - New option configuration system and functions `set_option`, `get_option`,
    `describe_option`, and `reset_option`. Deprecate `set_printoptions` and
    `reset_printoptions` (:issue:`2393`).
    You can also access options as attributes via ``pandas.options.X``
  - Wide DataFrames can be viewed more easily in the console with new
    `expand_frame_repr` and `line_width` configuration options. This is on by
    default now (:issue:`2436`)
  - Scikits.timeseries-like moving window functions via ``rolling_window`` (:issue:`1270`)

**Experimental Features**

  - Add support for Panel4D, a named 4 Dimensional stucture
  - Add support for ndpanel factory functions, to create custom,
    domain-specific N-Dimensional containers

**API Changes**

  - The default binning/labeling behavior for ``resample`` has been changed to
    `closed='left', label='left'` for daily and lower frequencies. This had
    been a large source of confusion for users. See "what's new" page for more
    on this. (:issue:`2410`)
  - Methods with ``inplace`` option now return None instead of the calling
    (modified) object (:issue:`1893`)
  - The special case DataFrame - TimeSeries doing column-by-column broadcasting
    has been deprecated. Users should explicitly do e.g. df.sub(ts, axis=0)
    instead. This is a legacy hack and can lead to subtle bugs.
  - inf/-inf are no longer considered as NA by isnull/notnull. To be clear, this
    is legacy cruft from early pandas. This behavior can be globally re-enabled
    using the new option ``mode.use_inf_as_null`` (:issue:`2050`, :issue:`1919`)
  - ``pandas.merge`` will now default to ``sort=False``. For many use cases
    sorting the join keys is not necessary, and doing it by default is wasteful
  - Specify ``header=0`` explicitly to replace existing column names in file in
    read_* functions.
  - Default column names for header-less parsed files (yielded by read_csv,
    etc.) are now the integers 0, 1, .... A new argument `prefix` has been
    added; to get the v0.9.x behavior specify ``prefix='X'`` (:issue:`2034`). This API
    change was made to make the default column names more consistent with the
    DataFrame constructor's default column names when none are specified.
  - DataFrame selection using a boolean frame now preserves input shape
  - If function passed to Series.apply yields a Series, result will be a
    DataFrame (:issue:`2316`)
  - Values like YES/NO/yes/no will not be considered as boolean by default any
    longer in the file parsers. This can be customized using the new
    ``true_values`` and ``false_values`` options (:issue:`2360`)
  - `obj.fillna()` is no longer valid; make `method='pad'` no longer the
    default option, to be more explicit about what kind of filling to
    perform. Add `ffill/bfill` convenience functions per above (:issue:`2284`)
  - `HDFStore.keys()` now returns an absolute path-name for each key
  - `to_string()` now always returns a unicode string. (:issue:`2224`)
  - File parsers will not handle NA sentinel values arising from passed
    converter functions

**Improvements to existing features**

  - Add ``nrows`` option to DataFrame.from_records for iterators (:issue:`1794`)
  - Unstack/reshape algorithm rewrite to avoid high memory use in cases where
    the number of observed key-tuples is much smaller than the total possible
    number that could occur (:issue:`2278`). Also improves performance in most cases.
  - Support duplicate columns in DataFrame.from_records (:issue:`2179`)
  - Add ``normalize`` option to Series/DataFrame.asfreq (:issue:`2137`)
  - SparseSeries and SparseDataFrame construction from empty and scalar
    values now no longer create dense ndarrays unnecessarily (:issue:`2322`)
  - ``HDFStore`` now supports hierarchial keys (:issue:`2397`)
  - Support multiple query selection formats for ``HDFStore tables`` (:issue:`1996`)
  - Support ``del store['df']`` syntax to delete HDFStores
  - Add multi-dtype support for ``HDFStore tables``
  - ``min_itemsize`` parameter can be specified in ``HDFStore table`` creation
  - Indexing support in ``HDFStore tables`` (:issue:`698`)
  - Add `line_terminator` option to DataFrame.to_csv (:issue:`2383`)
  - added implementation of str(x)/unicode(x)/bytes(x) to major pandas data
    structures, which should do the right thing on both py2.x and py3.x. (:issue:`2224`)
  - Reduce groupby.apply overhead substantially by low-level manipulation of
    internal NumPy arrays in DataFrames (:issue:`535`)
  - Implement ``value_vars`` in ``melt`` and add ``melt`` to pandas namespace
    (:issue:`2412`)
  - Added boolean comparison operators to Panel
  - Enable ``Series.str.strip/lstrip/rstrip`` methods to take an argument (:issue:`2411`)
  - The DataFrame ctor now respects column ordering when given
    an OrderedDict (:issue:`2455`)
  - Assigning DatetimeIndex to Series changes the class to TimeSeries (:issue:`2139`)
  - Improve performance of .value_counts method on non-integer data (:issue:`2480`)
  - ``get_level_values`` method for MultiIndex return Index instead of ndarray (:issue:`2449`)
  - ``convert_to_r_dataframe`` conversion for datetime values (:issue:`2351`)
  - Allow ``DataFrame.to_csv`` to represent inf and nan differently (:issue:`2026`)
  - Add ``min_i`` argument to ``nancorr`` to specify minimum required observations (:issue:`2002`)
  - Add ``inplace`` option to ``sortlevel`` / ``sort`` functions on DataFrame (:issue:`1873`)
  - Enable DataFrame to accept scalar constructor values like Series (:issue:`1856`)
  - DataFrame.from_records now takes optional ``size`` parameter (:issue:`1794`)
  - include iris dataset (:issue:`1709`)
  - No datetime64 DataFrame column conversion of datetime.datetime with tzinfo (:issue:`1581`)
  - Micro-optimizations in DataFrame for tracking state of internal consolidation (:issue:`217`)
  - Format parameter in DataFrame.to_csv (:issue:`1525`)
  - Partial string slicing for ``DatetimeIndex`` for daily and higher frequencies (:issue:`2306`)
  - Implement ``col_space`` parameter in ``to_html`` and ``to_string`` in DataFrame (:issue:`1000`)
  - Override ``Series.tolist`` and box datetime64 types (:issue:`2447`)
  - Optimize ``unstack`` memory usage by compressing indices (:issue:`2278`)
  - Fix HTML repr in IPython qtconsole if opening window is small (:issue:`2275`)
  - Escape more special characters in console output (:issue:`2492`)
  - df.select now invokes bool on the result of crit(x) (:issue:`2487`)

**Bug fixes**

  - Fix major performance regression in DataFrame.iteritems (:issue:`2273`)
  - Fixes bug when negative period passed to Series/DataFrame.diff (:issue:`2266`)
  - Escape tabs in console output to avoid alignment issues (:issue:`2038`)
  - Properly box datetime64 values when retrieving cross-section from
    mixed-dtype DataFrame (:issue:`2272`)
  - Fix concatenation bug leading to :issue:`2057`, :issue:`2257`
  - Fix regression in Index console formatting (:issue:`2319`)
  - Box Period data when assigning PeriodIndex to frame column (:issue:`2243`, :issue:`2281`)
  - Raise exception on calling reset_index on Series with inplace=True (:issue:`2277`)
  - Enable setting multiple columns in DataFrame with hierarchical columns
    (:issue:`2295`)
  - Respect dtype=object in DataFrame constructor (:issue:`2291`)
  - Fix DatetimeIndex.join bug with tz-aware indexes and how='outer' (:issue:`2317`)
  - pop(...) and del works with DataFrame with duplicate columns (:issue:`2349`)
  - Treat empty strings as NA in date parsing (rather than let dateutil do
    something weird) (:issue:`2263`)
  - Prevent uint64 -> int64 overflows (:issue:`2355`)
  - Enable joins between MultiIndex and regular Index (:issue:`2024`)
  - Fix time zone metadata issue when unioning non-overlapping DatetimeIndex
    objects (:issue:`2367`)
  - Raise/handle int64 overflows in parsers (:issue:`2247`)
  - Deleting of consecutive rows in ``HDFStore tables``` is much faster than before
  - Appending on a HDFStore would fail if the table was not first created via ``put``
  - Use `col_space` argument as minimum column width in DataFrame.to_html (:issue:`2328`)
  - Fix tz-aware DatetimeIndex.to_period (:issue:`2232`)
  - Fix DataFrame row indexing case with MultiIndex (:issue:`2314`)
  - Fix to_excel exporting issues with Timestamp objects in index (:issue:`2294`)
  - Fixes assigning scalars and array to hierarchical column chunk (:issue:`1803`)
  - Fixed a UnicdeDecodeError with series tidy_repr (:issue:`2225`)
  - Fixed issued with duplicate keys in an index (:issue:`2347`, :issue:`2380`)
  - Fixed issues re: Hash randomization, default on starting w/ py3.3 (:issue:`2331`)
  - Fixed issue with missing attributes after loading a pickled dataframe (:issue:`2431`)
  - Fix Timestamp formatting with tzoffset time zone in dateutil 2.1 (:issue:`2443`)
  - Fix GroupBy.apply issue when using BinGrouper to do ts binning (:issue:`2300`)
  - Fix issues resulting from datetime.datetime columns being converted to
    datetime64 when calling DataFrame.apply. (:issue:`2374`)
  - Raise exception when calling to_panel on non uniquely-indexed frame (:issue:`2441`)
  - Improved detection of console encoding on IPython zmq frontends (:issue:`2458`)
  - Preserve time zone when .append-ing two time series (:issue:`2260`)
  - Box timestamps when calling reset_index on time-zone-aware index rather
    than creating a tz-less datetime64 column (:issue:`2262`)
  - Enable searching non-string columns in DataFrame.filter(like=...) (:issue:`2467`)
  - Fixed issue with losing nanosecond precision upon conversion to DatetimeIndex(:issue:`2252`)
  - Handle timezones in Datetime.normalize (:issue:`2338`)
  - Fix test case where dtype specification with endianness causes
    failures on big endian machines (:issue:`2318`)
  - Fix plotting bug where upsampling causes data to appear shifted in time (:issue:`2448`)
  - Fix ``read_csv`` failure for UTF-16 with BOM and skiprows(:issue:`2298`)
  - read_csv with names arg not implicitly setting header=None(:issue:`2459`)
  - Unrecognized compression mode causes segfault in read_csv(:issue:`2474`)
  - In read_csv, header=0 and passed names should discard first row(:issue:`2269`)
  - Correctly route to stdout/stderr in read_table (:issue:`2071`)
  - Fix exception when Timestamp.to_datetime is called on a Timestamp with tzoffset (:issue:`2471`)
  - Fixed unintentional conversion of datetime64 to long in groupby.first() (:issue:`2133`)
  - Union of empty DataFrames now return empty with concatenated index (:issue:`2307`)
  - DataFrame.sort_index raises more helpful exception if sorting by column
    with duplicates (:issue:`2488`)
  - DataFrame.to_string formatters can be list, too (:issue:`2520`)
  - DataFrame.combine_first will always result in the union of the index and
    columns, even if one DataFrame is length-zero (:issue:`2525`)
  - Fix several DataFrame.icol/irow with duplicate indices issues (:issue:`2228`, :issue:`2259`)
  - Use Series names for column names when using concat with axis=1 (:issue:`2489`)
  - Raise Exception if start, end, periods all passed to date_range (:issue:`2538`)
  - Fix Panel resampling issue (:issue:`2537`)



pandas 0.9.1
============

**Release date:** 2012-11-14

**New features**

  - Can specify multiple sort orders in DataFrame/Series.sort/sort_index (:issue:`928`)
  - New `top` and `bottom` options for handling NAs in rank (:issue:`1508`, :issue:`2159`)
  - Add `where` and `mask` functions to DataFrame (:issue:`2109`, :issue:`2151`)
  - Add `at_time` and `between_time` functions to DataFrame (:issue:`2149`)
  - Add flexible `pow` and `rpow` methods to DataFrame (:issue:`2190`)

**API Changes**

  - Upsampling period index "spans" intervals. Example: annual periods
    upsampled to monthly will span all months in each year
  - Period.end_time will yield timestamp at last nanosecond in the interval
    (:issue:`2124`, :issue:`2125`, :issue:`1764`)
  - File parsers no longer coerce to float or bool for columns that have custom
    converters specified (:issue:`2184`)

**Improvements to existing features**

  - Time rule inference for week-of-month (e.g. WOM-2FRI) rules (:issue:`2140`)
  - Improve performance of datetime + business day offset with large number of
    offset periods
  - Improve HTML display of DataFrame objects with hierarchical columns
  - Enable referencing of Excel columns by their column names (:issue:`1936`)
  - DataFrame.dot can accept ndarrays (:issue:`2042`)
  - Support negative periods in Panel.shift (:issue:`2164`)
  - Make .drop(...) work with non-unique indexes (:issue:`2101`)
  - Improve performance of Series/DataFrame.diff (re: :issue:`2087`)
  - Support unary ~ (__invert__) in DataFrame (:issue:`2110`)
  - Turn off pandas-style tick locators and formatters (:issue:`2205`)
  - DataFrame[DataFrame] uses DataFrame.where to compute masked frame (:issue:`2230`)

**Bug fixes**

  - Fix some duplicate-column DataFrame constructor issues (:issue:`2079`)
  - Fix bar plot color cycle issues (:issue:`2082`)
  - Fix off-center grid for stacked bar plots (:issue:`2157`)
  - Fix plotting bug if inferred frequency is offset with N > 1 (:issue:`2126`)
  - Implement comparisons on date offsets with fixed delta (:issue:`2078`)
  - Handle inf/-inf correctly in read_* parser functions (:issue:`2041`)
  - Fix matplotlib unicode interaction bug
  - Make WLS r-squared match statsmodels 0.5.0 fixed value
  - Fix zero-trimming DataFrame formatting bug
  - Correctly compute/box datetime64 min/max values from Series.min/max (:issue:`2083`)
  - Fix unstacking edge case with unrepresented groups (:issue:`2100`)
  - Fix Series.str failures when using pipe pattern '|' (:issue:`2119`)
  - Fix pretty-printing of dict entries in Series, DataFrame (:issue:`2144`)
  - Cast other datetime64 values to nanoseconds in DataFrame ctor (:issue:`2095`)
  - Alias Timestamp.astimezone to tz_convert, so will yield Timestamp (:issue:`2060`)
  - Fix timedelta64 formatting from Series (:issue:`2165`, :issue:`2146`)
  - Handle None values gracefully in dict passed to Panel constructor (:issue:`2075`)
  - Box datetime64 values as Timestamp objects in Series/DataFrame.iget (:issue:`2148`)
  - Fix Timestamp indexing bug in DatetimeIndex.insert (:issue:`2155`)
  - Use index name(s) (if any) in DataFrame.to_records (:issue:`2161`)
  - Don't lose index names in Panel.to_frame/DataFrame.to_panel (:issue:`2163`)
  - Work around length-0 boolean indexing NumPy bug (:issue:`2096`)
  - Fix partial integer indexing bug in DataFrame.xs (:issue:`2107`)
  - Fix variety of cut/qcut string-bin formatting bugs (:issue:`1978`, :issue:`1979`)
  - Raise Exception when xs view not possible of MultiIndex'd DataFrame (:issue:`2117`)
  - Fix groupby(...).first() issue with datetime64 (:issue:`2133`)
  - Better floating point error robustness in some rolling_* functions
    (:issue:`2114`, :issue:`2527`)
  - Fix ewma NA handling in the middle of Series (:issue:`2128`)
  - Fix numerical precision issues in diff with integer data (:issue:`2087`)
  - Fix bug in MultiIndex.__getitem__ with NA values (:issue:`2008`)
  - Fix DataFrame.from_records dict-arg bug when passing columns (:issue:`2179`)
  - Fix Series and DataFrame.diff for integer dtypes (:issue:`2087`, :issue:`2174`)
  - Fix bug when taking intersection of DatetimeIndex with empty index (:issue:`2129`)
  - Pass through timezone information when calling DataFrame.align (:issue:`2127`)
  - Properly sort when joining on datetime64 values (:issue:`2196`)
  - Fix indexing bug in which False/True were being coerced to 0/1 (:issue:`2199`)
  - Many unicode formatting fixes (:issue:`2201`)
  - Fix improper MultiIndex conversion issue when assigning
    e.g. DataFrame.index (:issue:`2200`)
  - Fix conversion of mixed-type DataFrame to ndarray with dup columns (:issue:`2236`)
  - Fix duplicate columns issue (:issue:`2218`, :issue:`2219`)
  - Fix SparseSeries.__pow__ issue with NA input (:issue:`2220`)
  - Fix icol with integer sequence failure (:issue:`2228`)
  - Fixed resampling tz-aware time series issue (:issue:`2245`)
  - SparseDataFrame.icol was not returning SparseSeries (:issue:`2227`, :issue:`2229`)
  - Enable ExcelWriter to handle PeriodIndex (:issue:`2240`)
  - Fix issue constructing DataFrame from empty Series with name (:issue:`2234`)
  - Use console-width detection in interactive sessions only (:issue:`1610`)
  - Fix parallel_coordinates legend bug with mpl 1.2.0 (:issue:`2237`)
  - Make tz_localize work in corner case of empty Series (:issue:`2248`)



pandas 0.9.0
============

**Release date:** 10/7/2012

**New features**

  - Add ``str.encode`` and ``str.decode`` to Series (:issue:`1706`)
  - Add `to_latex` method to DataFrame (:issue:`1735`)
  - Add convenient expanding window equivalents of all rolling_* ops (:issue:`1785`)
  - Add Options class to pandas.io.data for fetching options data from Yahoo!
    Finance (:issue:`1748`, :issue:`1739`)
  - Recognize and convert more boolean values in file parsing (Yes, No, TRUE,
    FALSE, variants thereof) (:issue:`1691`, :issue:`1295`)
  - Add Panel.update method, analogous to DataFrame.update (:issue:`1999`, :issue:`1988`)

**Improvements to existing features**

  - Proper handling of NA values in merge operations (:issue:`1990`)
  - Add ``flags`` option for ``re.compile`` in some Series.str methods (:issue:`1659`)
  - Parsing of UTC date strings in read_* functions (:issue:`1693`)
  - Handle generator input to Series (:issue:`1679`)
  - Add `na_action='ignore'` to Series.map to quietly propagate NAs (:issue:`1661`)
  - Add args/kwds options to Series.apply (:issue:`1829`)
  - Add inplace option to Series/DataFrame.reset_index (:issue:`1797`)
  - Add ``level`` parameter to ``Series.reset_index``
  - Add quoting option for DataFrame.to_csv (:issue:`1902`)
  - Indicate long column value truncation in DataFrame output with ... (:issue:`1854`)
  - DataFrame.dot will not do data alignment, and also work with Series (:issue:`1915`)
  - Add ``na`` option for missing data handling in some vectorized string
    methods (:issue:`1689`)
  - If index_label=False in DataFrame.to_csv, do not print fields/commas in the
    text output. Results in easier importing into R (:issue:`1583`)
  - Can pass tuple/list of axes to DataFrame.dropna to simplify repeated calls
    (dropping both columns and rows) (:issue:`924`)
  - Improve DataFrame.to_html output for hierarchically-indexed rows (do not
    repeat levels) (:issue:`1929`)
  - TimeSeries.between_time can now select times across midnight (:issue:`1871`)
  - Enable `skip_footer` parameter in `ExcelFile.parse` (:issue:`1843`)

**API Changes**

  - Change default header names in read_* functions to more Pythonic X0, X1,
    etc. instead of X.1, X.2. (:issue:`2000`)
  - Deprecated ``day_of_year`` API removed from PeriodIndex, use ``dayofyear``
    (:issue:`1723`)
  - Don't modify NumPy suppress printoption at import time
  - The internal HDF5 data arrangement for DataFrames has been
    transposed. Legacy files will still be readable by HDFStore (:issue:`1834`, :issue:`1824`)
  - Legacy cruft removed: pandas.stats.misc.quantileTS
  - Use ISO8601 format for Period repr: monthly, daily, and on down (:issue:`1776`)
  - Empty DataFrame columns are now created as object dtype. This will prevent
    a class of TypeErrors that was occurring in code where the dtype of a
    column would depend on the presence of data or not (e.g. a SQL query having
    results) (:issue:`1783`)
  - Setting parts of DataFrame/Panel using ix now aligns input Series/DataFrame
    (:issue:`1630`)
  - `first` and `last` methods in `GroupBy` no longer drop non-numeric columns
    (:issue:`1809`)
  - Resolved inconsistencies in specifying custom NA values in text parser.
    `na_values` of type dict no longer override default NAs unless
    `keep_default_na` is set to false explicitly (:issue:`1657`)
  - Enable `skipfooter` parameter in text parsers as an alias for `skip_footer`

**Bug fixes**

  - Perform arithmetic column-by-column in mixed-type DataFrame to avoid type
    upcasting issues. Caused downstream DataFrame.diff bug (:issue:`1896`)
  - Fix matplotlib auto-color assignment when no custom spectrum passed. Also
    respect passed color keyword argument (:issue:`1711`)
  - Fix resampling logical error with closed='left' (:issue:`1726`)
  - Fix critical DatetimeIndex.union bugs (:issue:`1730`, :issue:`1719`, :issue:`1745`, :issue:`1702`, :issue:`1753`)
  - Fix critical DatetimeIndex.intersection bug with unanchored offsets (:issue:`1708`)
  - Fix MM-YYYY time series indexing case (:issue:`1672`)
  - Fix case where Categorical group key was not being passed into index in
    GroupBy result (:issue:`1701`)
  - Handle Ellipsis in Series.__getitem__/__setitem__ (:issue:`1721`)
  - Fix some bugs with handling datetime64 scalars of other units in NumPy 1.6
    and 1.7 (:issue:`1717`)
  - Fix performance issue in MultiIndex.format (:issue:`1746`)
  - Fixed GroupBy bugs interacting with DatetimeIndex asof / map methods (:issue:`1677`)
  - Handle factors with NAs in pandas.rpy (:issue:`1615`)
  - Fix statsmodels import in pandas.stats.var (:issue:`1734`)
  - Fix DataFrame repr/info summary with non-unique columns (:issue:`1700`)
  - Fix Series.iget_value for non-unique indexes (:issue:`1694`)
  - Don't lose tzinfo when passing DatetimeIndex as DataFrame column (:issue:`1682`)
  - Fix tz conversion with time zones that haven't had any DST transitions since
    first date in the array (:issue:`1673`)
  - Fix field access with  UTC->local conversion on unsorted arrays (:issue:`1756`)
  - Fix isnull handling of array-like (list) inputs (:issue:`1755`)
  - Fix regression in handling of Series in Series constructor (:issue:`1671`)
  - Fix comparison of Int64Index with DatetimeIndex (:issue:`1681`)
  - Fix min_periods handling in new rolling_max/min at array start (:issue:`1695`)
  - Fix errors with how='median' and generic NumPy resampling in some cases
    caused by SeriesBinGrouper (:issue:`1648`, :issue:`1688`)
  - When grouping by level, exclude unobserved levels (:issue:`1697`)
  - Don't lose tzinfo in DatetimeIndex when shifting by different offset (:issue:`1683`)
  - Hack to support storing data with a zero-length axis in HDFStore (:issue:`1707`)
  - Fix DatetimeIndex tz-aware range generation issue (:issue:`1674`)
  - Fix method='time' interpolation with intraday data (:issue:`1698`)
  - Don't plot all-NA DataFrame columns as zeros (:issue:`1696`)
  - Fix bug in scatter_plot with by option (:issue:`1716`)
  - Fix performance problem in infer_freq with lots of non-unique stamps (:issue:`1686`)
  - Fix handling of PeriodIndex as argument to create MultiIndex (:issue:`1705`)
  - Fix re: unicode MultiIndex level names in Series/DataFrame repr (:issue:`1736`)
  - Handle PeriodIndex in to_datetime instance method (:issue:`1703`)
  - Support StaticTzInfo in DatetimeIndex infrastructure (:issue:`1692`)
  - Allow MultiIndex setops with length-0 other type indexes (:issue:`1727`)
  - Fix handling of DatetimeIndex in DataFrame.to_records (:issue:`1720`)
  - Fix handling of general objects in isnull on which bool(...) fails (:issue:`1749`)
  - Fix .ix indexing with MultiIndex ambiguity (:issue:`1678`)
  - Fix .ix setting logic error with non-unique MultiIndex (:issue:`1750`)
  - Basic indexing now works on MultiIndex with > 1000000 elements, regression
    from earlier version of pandas (:issue:`1757`)
  - Handle non-float64 dtypes in fast DataFrame.corr/cov code paths (:issue:`1761`)
  - Fix DatetimeIndex.isin to function properly (:issue:`1763`)
  - Fix conversion of array of tz-aware datetime.datetime to DatetimeIndex with
    right time zone (:issue:`1777`)
  - Fix DST issues with generating ancxhored date ranges (:issue:`1778`)
  - Fix issue calling sort on result of Series.unique (:issue:`1807`)
  - Fix numerical issue leading to square root of negative number in
    rolling_std (:issue:`1840`)
  - Let Series.str.split accept no arguments (like str.split) (:issue:`1859`)
  - Allow user to have dateutil 2.1 installed on a Python 2 system (:issue:`1851`)
  - Catch ImportError less aggressively in pandas/__init__.py (:issue:`1845`)
  - Fix pip source installation bug when installing from GitHub (:issue:`1805`)
  - Fix error when window size > array size in rolling_apply (:issue:`1850`)
  - Fix pip source installation issues via SSH from GitHub
  - Fix OLS.summary when column is a tuple (:issue:`1837`)
  - Fix bug in __doc__ patching when -OO passed to interpreter
    (:issue:`1792` :issue:`1741` :issue:`1774`)
  - Fix unicode console encoding issue in IPython notebook (:issue:`1782`, :issue:`1768`)
  - Fix unicode formatting issue with Series.name (:issue:`1782`)
  - Fix bug in DataFrame.duplicated with datetime64 columns (:issue:`1833`)
  - Fix bug in Panel internals resulting in error when doing fillna after
    truncate not changing size of panel (:issue:`1823`)
  - Prevent segfault due to MultiIndex not being supported in HDFStore table
    format (:issue:`1848`)
  - Fix UnboundLocalError in Panel.__setitem__ and add better error (:issue:`1826`)
  - Fix to_csv issues with list of string entries. Isnull works on list of
    strings now too (:issue:`1791`)
  - Fix Timestamp comparisons with datetime values outside the nanosecond range
    (1677-2262)
  - Revert to prior behavior of normalize_date with datetime.date objects
    (return datetime)
  - Fix broken interaction between np.nansum and Series.any/all
  - Fix bug with multiple column date parsers (:issue:`1866`)
  - DatetimeIndex.union(Int64Index) was broken
  - Make plot x vs y interface consistent with integer indexing (:issue:`1842`)
  - set_index inplace modified data even if unique check fails (:issue:`1831`)
  - Only use Q-OCT/NOV/DEC in quarterly frequency inference (:issue:`1789`)
  - Upcast to dtype=object when unstacking boolean DataFrame (:issue:`1820`)
  - Fix float64/float32 merging bug (:issue:`1849`)
  - Fixes to Period.start_time for non-daily frequencies (:issue:`1857`)
  - Fix failure when converter used on index_col in read_csv (:issue:`1835`)
  - Implement PeriodIndex.append so that pandas.concat works correctly (:issue:`1815`)
  - Avoid Cython out-of-bounds access causing segfault sometimes in pad_2d,
    backfill_2d
  - Fix resampling error with intraday times and anchored target time (like
    AS-DEC) (:issue:`1772`)
  - Fix .ix indexing bugs with mixed-integer indexes (:issue:`1799`)
  - Respect passed color keyword argument in Series.plot (:issue:`1890`)
  - Fix rolling_min/max when the window is larger than the size of the input
    array. Check other malformed inputs (:issue:`1899`, :issue:`1897`)
  - Rolling variance / standard deviation with only a single observation in
    window (:issue:`1884`)
  - Fix unicode sheet name failure in to_excel (:issue:`1828`)
  - Override DatetimeIndex.min/max to return Timestamp objects (:issue:`1895`)
  - Fix column name formatting issue in length-truncated column (:issue:`1906`)
  - Fix broken handling of copying Index metadata to new instances created by
    view(...) calls inside the NumPy infrastructure
  - Support datetime.date again in DateOffset.rollback/rollforward
  - Raise Exception if set passed to Series constructor (:issue:`1913`)
  - Add TypeError when appending HDFStore table w/ wrong index type (:issue:`1881`)
  - Don't raise exception on empty inputs in EW functions (e.g. ewma) (:issue:`1900`)
  - Make asof work correctly with PeriodIndex (:issue:`1883`)
  - Fix extlinks in doc build
  - Fill boolean DataFrame with NaN when calling shift (:issue:`1814`)
  - Fix setuptools bug causing pip not to Cythonize .pyx files sometimes
  - Fix negative integer indexing regression in .ix from 0.7.x (:issue:`1888`)
  - Fix error while retrieving timezone and utc offset from subclasses of
    datetime.tzinfo without .zone and ._utcoffset attributes (:issue:`1922`)
  - Fix DataFrame formatting of small, non-zero FP numbers (:issue:`1911`)
  - Various fixes by upcasting of date -> datetime (:issue:`1395`)
  - Raise better exception when passing multiple functions with the same name,
    such as lambdas, to GroupBy.aggregate
  - Fix DataFrame.apply with axis=1 on a non-unique index (:issue:`1878`)
  - Proper handling of Index subclasses in pandas.unique (:issue:`1759`)
  - Set index names in DataFrame.from_records (:issue:`1744`)
  - Fix time series indexing error with duplicates, under and over hash table
    size cutoff (:issue:`1821`)
  - Handle list keys in addition to tuples in DataFrame.xs when
    partial-indexing a hierarchically-indexed DataFrame (:issue:`1796`)
  - Support multiple column selection in DataFrame.__getitem__ with duplicate
    columns (:issue:`1943`)
  - Fix time zone localization bug causing improper fields (e.g. hours) in time
    zones that have not had a UTC transition in a long time (:issue:`1946`)
  - Fix errors when parsing and working with with fixed offset timezones
    (:issue:`1922`, :issue:`1928`)
  - Fix text parser bug when handling UTC datetime objects generated by
    dateutil (:issue:`1693`)
  - Fix plotting bug when 'B' is the inferred frequency but index actually
    contains weekends (:issue:`1668`, :issue:`1669`)
  - Fix plot styling bugs (:issue:`1666`, :issue:`1665`, :issue:`1658`)
  - Fix plotting bug with index/columns with unicode (:issue:`1685`)
  - Fix DataFrame constructor bug when passed Series with datetime64 dtype
    in a dict (:issue:`1680`)
  - Fixed regression in generating DatetimeIndex using timezone aware
    datetime.datetime (:issue:`1676`)
  - Fix DataFrame bug when printing concatenated DataFrames with duplicated
    columns (:issue:`1675`)
  - Fixed bug when plotting time series with multiple intraday frequencies
    (:issue:`1732`)
  - Fix bug in DataFrame.duplicated to enable iterables other than list-types
    as input argument (:issue:`1773`)
  - Fix resample bug when passed list of lambdas as `how` argument (:issue:`1808`)
  - Repr fix for MultiIndex level with all NAs (:issue:`1971`)
  - Fix PeriodIndex slicing bug when slice start/end are out-of-bounds (:issue:`1977`)
  - Fix read_table bug when parsing unicode (:issue:`1975`)
  - Fix BlockManager.iget bug when dealing with non-unique MultiIndex as columns
    (:issue:`1970`)
  - Fix reset_index bug if both drop and level are specified (:issue:`1957`)
  - Work around unsafe NumPy object->int casting with Cython function (:issue:`1987`)
  - Fix datetime64 formatting bug in DataFrame.to_csv (:issue:`1993`)
  - Default start date in pandas.io.data to 1/1/2000 as the docs say (:issue:`2011`)




pandas 0.8.1
============

**Release date:** July 22, 2012

**New features**

  - Add vectorized, NA-friendly string methods to Series (:issue:`1621`, :issue:`620`)
  - Can pass dict of per-column line styles to DataFrame.plot (:issue:`1559`)
  - Selective plotting to secondary y-axis on same subplot (:issue:`1640`)
  - Add new ``bootstrap_plot`` plot function
  - Add new ``parallel_coordinates`` plot function (:issue:`1488`)
  - Add ``radviz`` plot function (:issue:`1566`)
  - Add ``multi_sparse`` option to ``set_printoptions`` to modify display of
    hierarchical indexes (:issue:`1538`)
  - Add ``dropna`` method to Panel (:issue:`171`)

**Improvements to existing features**

  - Use moving min/max algorithms from Bottleneck in rolling_min/rolling_max
    for > 100x speedup. (:issue:`1504`, :issue:`50`)
  - Add Cython group median method for >15x speedup (:issue:`1358`)
  - Drastically improve ``to_datetime`` performance on ISO8601 datetime strings
    (with no time zones) (:issue:`1571`)
  - Improve single-key groupby performance on large data sets, accelerate use of
    groupby with a Categorical variable
  - Add ability to append hierarchical index levels with ``set_index`` and to
    drop single levels with ``reset_index`` (:issue:`1569`, :issue:`1577`)
  - Always apply passed functions in ``resample``, even if upsampling (:issue:`1596`)
  - Avoid unnecessary copies in DataFrame constructor with explicit dtype (:issue:`1572`)
  - Cleaner DatetimeIndex string representation with 1 or 2 elements (:issue:`1611`)
  - Improve performance of array-of-Period to PeriodIndex, convert such arrays
    to PeriodIndex inside Index (:issue:`1215`)
  - More informative string representation for weekly Period objects (:issue:`1503`)
  - Accelerate 3-axis multi data selection from homogeneous Panel (:issue:`979`)
  - Add ``adjust`` option to ewma to disable adjustment factor (:issue:`1584`)
  - Add new matplotlib converters for high frequency time series plotting (:issue:`1599`)
  - Handling of tz-aware datetime.datetime objects in to_datetime; raise
    Exception unless utc=True given (:issue:`1581`)

**Bug fixes**

  - Fix NA handling in DataFrame.to_panel (:issue:`1582`)
  - Handle TypeError issues inside PyObject_RichCompareBool calls in khash
    (:issue:`1318`)
  - Fix resampling bug to lower case daily frequency (:issue:`1588`)
  - Fix kendall/spearman DataFrame.corr bug with no overlap (:issue:`1595`)
  - Fix bug in DataFrame.set_index (:issue:`1592`)
  - Don't ignore axes in boxplot if by specified (:issue:`1565`)
  - Fix Panel .ix indexing with integers bug (:issue:`1603`)
  - Fix Partial indexing bugs (years, months, ...) with PeriodIndex (:issue:`1601`)
  - Fix MultiIndex console formatting issue (:issue:`1606`)
  - Unordered index with duplicates doesn't yield scalar location for single
    entry (:issue:`1586`)
  - Fix resampling of tz-aware time series with "anchored" freq (:issue:`1591`)
  - Fix DataFrame.rank error on integer data (:issue:`1589`)
  - Selection of multiple SparseDataFrame columns by list in __getitem__ (:issue:`1585`)
  - Override Index.tolist for compatibility with MultiIndex (:issue:`1576`)
  - Fix hierarchical summing bug with MultiIndex of length 1 (:issue:`1568`)
  - Work around numpy.concatenate use/bug in Series.set_value (:issue:`1561`)
  - Ensure Series/DataFrame are sorted before resampling (:issue:`1580`)
  - Fix unhandled IndexError when indexing very large time series (:issue:`1562`)
  - Fix DatetimeIndex intersection logic error with irregular indexes (:issue:`1551`)
  - Fix unit test errors on Python 3 (:issue:`1550`)
  - Fix .ix indexing bugs in duplicate DataFrame index (:issue:`1201`)
  - Better handle errors with non-existing objects in HDFStore (:issue:`1254`)
  - Don't copy int64 array data in DatetimeIndex when copy=False (:issue:`1624`)
  - Fix resampling of conforming periods quarterly to annual (:issue:`1622`)
  - Don't lose index name on resampling (:issue:`1631`)
  - Support python-dateutil version 2.1 (:issue:`1637`)
  - Fix broken scatter_matrix axis labeling, esp. with time series (:issue:`1625`)
  - Fix cases where extra keywords weren't being passed on to matplotlib from
    Series.plot (:issue:`1636`)
  - Fix BusinessMonthBegin logic for dates before 1st bday of month (:issue:`1645`)
  - Ensure string alias converted (valid in DatetimeIndex.get_loc) in
    DataFrame.xs / __getitem__ (:issue:`1644`)
  - Fix use of string alias timestamps with tz-aware time series (:issue:`1647`)
  - Fix Series.max/min and Series.describe on len-0 series (:issue:`1650`)
  - Handle None values in dict passed to concat (:issue:`1649`)
  - Fix Series.interpolate with method='values' and DatetimeIndex (:issue:`1646`)
  - Fix IndexError in left merges on a DataFrame with 0-length (:issue:`1628`)
  - Fix DataFrame column width display with UTF-8 encoded characters (:issue:`1620`)
  - Handle case in pandas.io.data.get_data_yahoo where Yahoo! returns duplicate
    dates for most recent business day
  - Avoid downsampling when plotting mixed frequencies on the same subplot (:issue:`1619`)
  - Fix read_csv bug when reading a single line (:issue:`1553`)
  - Fix bug in C code causing monthly periods prior to December 1969 to be off (:issue:`1570`)



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
  - Time series string indexing shorthand (:issue:`222`)
  - Add week, dayofyear array and other timestamp array-valued field accessor
    functions to DatetimeIndex
  - Add GroupBy.prod optimized aggregation function and 'prod' fast time series
    conversion method (:issue:`1018`)
  - Implement robust frequency inference function and `inferred_freq` attribute
    on DatetimeIndex (:issue:`391`)
  - New ``tz_convert`` and ``tz_localize`` methods in Series / DataFrame
  - Convert DatetimeIndexes to UTC if time zones are different in join/setops
    (:issue:`864`)
  - Add limit argument for forward/backward filling to reindex, fillna,
    etc. (:issue:`825` and others)
  - Add support for indexes (dates or otherwise) with duplicates and common
    sense indexing/selection functionality
  - Series/DataFrame.update methods, in-place variant of combine_first (:issue:`961`)
  - Add ``match`` function to API (:issue:`502`)
  - Add Cython-optimized first, last, min, max, prod functions to GroupBy (:issue:`994`,
    :issue:`1043`)
  - Dates can be split across multiple columns (:issue:`1227`, :issue:`1186`)
  - Add experimental support for converting pandas DataFrame to R data.frame
    via rpy2 (:issue:`350`, :issue:`1212`)
  - Can pass list of (name, function) to GroupBy.aggregate to get aggregates in
    a particular order (:issue:`610`)
  - Can pass dicts with lists of functions or dicts to GroupBy aggregate to do
    much more flexible multiple function aggregation (:issue:`642`, :issue:`610`)
  - New ordered_merge functions for merging DataFrames with ordered
    data. Also supports group-wise merging for panel data (:issue:`813`)
  - Add keys() method to DataFrame
  - Add flexible replace method for replacing potentially values to Series and
    DataFrame (:issue:`929`, :issue:`1241`)
  - Add 'kde' plot kind for Series/DataFrame.plot (:issue:`1059`)
  - More flexible multiple function aggregation with GroupBy
  - Add pct_change function to Series/DataFrame
  - Add option to interpolate by Index values in Series.interpolate (:issue:`1206`)
  - Add ``max_colwidth`` option for DataFrame, defaulting to 50
  - Conversion of DataFrame through rpy2 to R data.frame (:issue:`1282`, )
  - Add keys() method on DataFrame (:issue:`1240`)
  - Add new ``match`` function to API (similar to R) (:issue:`502`)
  - Add dayfirst option to parsers (:issue:`854`)
  - Add ``method`` argument to ``align`` method for forward/backward fillin
    (:issue:`216`)
  - Add Panel.transpose method for rearranging axes (:issue:`695`)
  - Add new ``cut`` function (patterned after R) for discretizing data into
    equal range-length bins or arbitrary breaks of your choosing (:issue:`415`)
  - Add new ``qcut`` for cutting with quantiles (:issue:`1378`)
  - Add ``value_counts`` top level array method (:issue:`1392`)
  - Added Andrews curves plot tupe (:issue:`1325`)
  - Add lag plot (:issue:`1440`)
  - Add autocorrelation_plot (:issue:`1425`)
  - Add support for tox and Travis CI (:issue:`1382`)
  - Add support for Categorical use in GroupBy (:issue:`292`)
  - Add ``any`` and ``all`` methods to DataFrame (:issue:`1416`)
  - Add ``secondary_y`` option to Series.plot
  - Add experimental ``lreshape`` function for reshaping wide to long

**Improvements to existing features**

  - Switch to klib/khash-based hash tables in Index classes for better
    performance in many cases and lower memory footprint
  - Shipping some functions from scipy.stats to reduce dependency,
    e.g. Series.describe and DataFrame.describe (:issue:`1092`)
  - Can create MultiIndex by passing list of lists or list of arrays to Series,
    DataFrame constructor, etc. (:issue:`831`)
  - Can pass arrays in addition to column names to DataFrame.set_index (:issue:`402`)
  - Improve the speed of "square" reindexing of homogeneous DataFrame objects
    by significant margin (:issue:`836`)
  - Handle more dtypes when passed MaskedArrays in DataFrame constructor (:issue:`406`)
  - Improved performance of join operations on integer keys (:issue:`682`)
  - Can pass multiple columns to GroupBy object, e.g. grouped[[col1, col2]] to
    only aggregate a subset of the value columns (:issue:`383`)
  - Add histogram / kde plot options for scatter_matrix diagonals (:issue:`1237`)
  - Add inplace option to Series/DataFrame.rename and sort_index,
    DataFrame.drop_duplicates (:issue:`805`, :issue:`207`)
  - More helpful error message when nothing passed to Series.reindex (:issue:`1267`)
  - Can mix array and scalars as dict-value inputs to DataFrame ctor (:issue:`1329`)
  - Use DataFrame columns' name for legend title in plots
  - Preserve frequency in DatetimeIndex when possible in boolean indexing
    operations
  - Promote datetime.date values in data alignment operations (:issue:`867`)
  - Add ``order`` method to Index classes (:issue:`1028`)
  - Avoid hash table creation in large monotonic hash table indexes (:issue:`1160`)
  - Store time zones in HDFStore (:issue:`1232`)
  - Enable storage of sparse data structures in HDFStore (:issue:`85`)
  - Enable Series.asof to work with arrays of timestamp inputs
  - Cython implementation of DataFrame.corr speeds up by > 100x (:issue:`1349`, :issue:`1354`)
  - Exclude "nuisance" columns automatically in GroupBy.transform (:issue:`1364`)
  - Support functions-as-strings in GroupBy.transform (:issue:`1362`)
  - Use index name as xlabel/ylabel in plots (:issue:`1415`)
  - Add ``convert_dtype`` option to Series.apply to be able to leave data as
    dtype=object (:issue:`1414`)
  - Can specify all index level names in concat (:issue:`1419`)
  - Add ``dialect`` keyword to parsers for quoting conventions (:issue:`1363`)
  - Enable DataFrame[bool_DataFrame] += value (:issue:`1366`)
  - Add ``retries`` argument to ``get_data_yahoo`` to try to prevent Yahoo! API
    404s (:issue:`826`)
  - Improve performance of reshaping by using O(N) categorical sorting
  - Series names will be used for index of DataFrame if no index passed (:issue:`1494`)
  - Header argument in DataFrame.to_csv can accept a list of column names to
    use instead of the object's columns (:issue:`921`)
  - Add ``raise_conflict`` argument to DataFrame.update (:issue:`1526`)
  - Support file-like objects in ExcelFile (:issue:`1529`)

**API Changes**

  - Rename `pandas._tseries` to `pandas.lib`
  - Rename Factor to Categorical and add improvements. Numerous Categorical bug
    fixes
  - Frequency name overhaul, WEEKDAY/EOM and rules with @
    deprecated. get_legacy_offset_name backwards compatibility function added
  - Raise ValueError in DataFrame.__nonzero__, so "if df" no longer works
    (:issue:`1073`)
  - Change BDay (business day) to not normalize dates by default (:issue:`506`)
  - Remove deprecated DataMatrix name
  - Default merge suffixes for overlap now have underscores instead of periods
    to facilitate tab completion, etc. (:issue:`1239`)
  - Deprecation of offset, time_rule timeRule parameters throughout codebase
  - Series.append and DataFrame.append no longer check for duplicate indexes
    by default, add verify_integrity parameter (:issue:`1394`)
  - Refactor Factor class, old constructor moved to Factor.from_array
  - Modified internals of MultiIndex to use less memory (no longer represented
    as array of tuples) internally, speed up construction time and many methods
    which construct intermediate hierarchical indexes (:issue:`1467`)

**Bug fixes**

  - Fix OverflowError from storing pre-1970 dates in HDFStore by switching to
    datetime64 (:issue:`179`)
  - Fix logical error with February leap year end in YearEnd offset
  - Series([False, nan]) was getting casted to float64 (:issue:`1074`)
  - Fix binary operations between boolean Series and object Series with
    booleans and NAs (:issue:`1074`, :issue:`1079`)
  - Couldn't assign whole array to column in mixed-type DataFrame via .ix
    (:issue:`1142`)
  - Fix label slicing issues with float index values (:issue:`1167`)
  - Fix segfault caused by empty groups passed to groupby (:issue:`1048`)
  - Fix occasionally misbehaved reindexing in the presence of NaN labels (:issue:`522`)
  - Fix imprecise logic causing weird Series results from .apply (:issue:`1183`)
  - Unstack multiple levels in one shot, avoiding empty columns in some
    cases. Fix pivot table bug (:issue:`1181`)
  - Fix formatting of MultiIndex on Series/DataFrame when index name coincides
    with label (:issue:`1217`)
  - Handle Excel 2003 #N/A as NaN from xlrd (:issue:`1213`, :issue:`1225`)
  - Fix timestamp locale-related deserialization issues with HDFStore by moving
    to datetime64 representation (:issue:`1081`, :issue:`809`)
  - Fix DataFrame.duplicated/drop_duplicates NA value handling (:issue:`557`)
  - Actually raise exceptions in fast reducer (:issue:`1243`)
  - Fix various timezone-handling bugs from 0.7.3 (:issue:`969`)
  - GroupBy on level=0 discarded index name (:issue:`1313`)
  - Better error message with unmergeable DataFrames (:issue:`1307`)
  - Series.__repr__ alignment fix with unicode index values (:issue:`1279`)
  - Better error message if nothing passed to reindex (:issue:`1267`)
  - More robust NA handling in DataFrame.drop_duplicates (:issue:`557`)
  - Resolve locale-based and pre-epoch HDF5 timestamp deserialization issues
    (:issue:`973`, :issue:`1081`, :issue:`179`)
  - Implement Series.repeat (:issue:`1229`)
  - Fix indexing with namedtuple and other tuple subclasses (:issue:`1026`)
  - Fix float64 slicing bug (:issue:`1167`)
  - Parsing integers with commas (:issue:`796`)
  - Fix groupby improper data type when group consists of one value (:issue:`1065`)
  - Fix negative variance possibility in nanvar resulting from floating point
    error (:issue:`1090`)
  - Consistently set name on groupby pieces (:issue:`184`)
  - Treat dict return values as Series in GroupBy.apply (:issue:`823`)
  - Respect column selection for DataFrame in in GroupBy.transform (:issue:`1365`)
  - Fix MultiIndex partial indexing bug (:issue:`1352`)
  - Enable assignment of rows in mixed-type DataFrame via .ix (:issue:`1432`)
  - Reset index mapping when grouping Series in Cython (:issue:`1423`)
  - Fix outer/inner DataFrame.join with non-unique indexes (:issue:`1421`)
  - Fix MultiIndex groupby bugs with empty lower levels (:issue:`1401`)
  - Calling fillna with a Series will have same behavior as with dict (:issue:`1486`)
  - SparseSeries reduction bug (:issue:`1375`)
  - Fix unicode serialization issue in HDFStore (:issue:`1361`)
  - Pass keywords to pyplot.boxplot in DataFrame.boxplot (:issue:`1493`)
  - Bug fixes in MonthBegin (:issue:`1483`)
  - Preserve MultiIndex names in drop (:issue:`1513`)
  - Fix Panel DataFrame slice-assignment bug (:issue:`1533`)
  - Don't use locals() in read_* functions (:issue:`1547`)



pandas 0.7.3
============

**Release date:** April 12, 2012

**New features / modules**

  - Support for non-unique indexes: indexing and selection, many-to-one and
    many-to-many joins (:issue:`1306`)
  - Added fixed-width file reader, read_fwf (:issue:`952`)
  - Add group_keys argument to groupby to not add group names to MultiIndex in
    result of apply (:issue:`938`)
  - DataFrame can now accept non-integer label slicing (:issue:`946`). Previously
    only DataFrame.ix was able to do so.
  - DataFrame.apply now retains name attributes on Series objects (:issue:`983`)
  - Numeric DataFrame comparisons with non-numeric values now raises proper
    TypeError (:issue:`943`). Previously raise "PandasError: DataFrame constructor
    not properly called!"
  - Add ``kurt`` methods to Series and DataFrame (:issue:`964`)
  - Can pass dict of column -> list/set NA values for text parsers (:issue:`754`)
  - Allows users specified NA values in text parsers (:issue:`754`)
  - Parsers checks for openpyxl dependency and raises ImportError if not found
    (:issue:`1007`)
  - New factory function to create HDFStore objects that can be used in a with
    statement so users do not have to explicitly call HDFStore.close (:issue:`1005`)
  - pivot_table is now more flexible with same parameters as groupby (:issue:`941`)
  - Added stacked bar plots (:issue:`987`)
  - scatter_matrix method in pandas/tools/plotting.py (:issue:`935`)
  - DataFrame.boxplot returns plot results for ex-post styling (:issue:`985`)
  - Short version number accessible as pandas.version.short_version (:issue:`930`)
  - Additional documentation in panel.to_frame (:issue:`942`)
  - More informative Series.apply docstring regarding element-wise apply
    (:issue:`977`)
  - Notes on rpy2 installation (:issue:`1006`)
  - Add rotation and font size options to hist method (:issue:`1012`)
  - Use exogenous / X variable index in result of OLS.y_predict. Add
    OLS.predict method (:issue:`1027`, :issue:`1008`)

**API Changes**

  - Calling apply on grouped Series, e.g. describe(), will no longer yield
    DataFrame by default. Will have to call unstack() to get prior behavior
  - NA handling in non-numeric comparisons has been tightened up (:issue:`933`, :issue:`953`)
  - No longer assign dummy names key_0, key_1, etc. to groupby index (:issue:`1291`)

**Bug fixes**

  - Fix logic error when selecting part of a row in a DataFrame with a
    MultiIndex index (:issue:`1013`)
  - Series comparison with Series of differing length causes crash (:issue:`1016`).
  - Fix bug in indexing when selecting section of hierarchically-indexed row
    (:issue:`1013`)
  - DataFrame.plot(logy=True) has no effect (:issue:`1011`).
  - Broken arithmetic operations between SparsePanel-Panel (:issue:`1015`)
  - Unicode repr issues in MultiIndex with non-ascii characters (:issue:`1010`)
  - DataFrame.lookup() returns inconsistent results if exact match not present
    (:issue:`1001`)
  - DataFrame arithmetic operations not treating None as NA (:issue:`992`)
  - DataFrameGroupBy.apply returns incorrect result (:issue:`991`)
  - Series.reshape returns incorrect result for multiple dimensions (:issue:`989`)
  - Series.std and Series.var ignores ddof parameter (:issue:`934`)
  - DataFrame.append loses index names (:issue:`980`)
  - DataFrame.plot(kind='bar') ignores color argument (:issue:`958`)
  - Inconsistent Index comparison results (:issue:`948`)
  - Improper int dtype DataFrame construction from data with NaN (:issue:`846`)
  - Removes default 'result' name in grouby results (:issue:`995`)
  - DataFrame.from_records no longer mutate input columns (:issue:`975`)
  - Use Index name when grouping by it (:issue:`1313`)



pandas 0.7.2
============

**Release date:** March 16, 2012

**New features / modules**

  - Add additional tie-breaking methods in DataFrame.rank (:issue:`874`)
  - Add ascending parameter to rank in Series, DataFrame (:issue:`875`)
  - Add sort_columns parameter to allow unsorted plots (:issue:`918`)
  - IPython tab completion on GroupBy objects

**API Changes**

  - Series.sum returns 0 instead of NA when called on an empty
    series. Analogously for a DataFrame whose rows or columns are length 0
    (:issue:`844`)

**Improvements to existing features**

  - Don't use groups dict in Grouper.size (:issue:`860`)
  - Use khash for Series.value_counts, add raw function to algorithms.py (:issue:`861`)
  - Enable column access via attributes on GroupBy (:issue:`882`)
  - Enable setting existing columns (only) via attributes on DataFrame, Panel
    (:issue:`883`)
  - Intercept __builtin__.sum in groupby (:issue:`885`)
  - Can pass dict to DataFrame.fillna to use different values per column (:issue:`661`)
  - Can select multiple hierarchical groups by passing list of values in .ix
    (:issue:`134`)
  - Add level keyword to ``drop`` for dropping values from a level (:issue:`159`)
  - Add ``coerce_float`` option on DataFrame.from_records (:issue:`893`)
  - Raise exception if passed date_parser fails in ``read_csv``
  - Add ``axis`` option to DataFrame.fillna (:issue:`174`)
  - Fixes to Panel to make it easier to subclass (:issue:`888`)

**Bug fixes**

  - Fix overflow-related bugs in groupby (:issue:`850`, :issue:`851`)
  - Fix unhelpful error message in parsers (:issue:`856`)
  - Better err msg for failed boolean slicing of dataframe (:issue:`859`)
  - Series.count cannot accept a string (level name) in the level argument (:issue:`869`)
  - Group index platform int check (:issue:`870`)
  - concat on axis=1 and ignore_index=True raises TypeError (:issue:`871`)
  - Further unicode handling issues resolved (:issue:`795`)
  - Fix failure in multiindex-based access in Panel (:issue:`880`)
  - Fix DataFrame boolean slice assignment failure (:issue:`881`)
  - Fix combineAdd NotImplementedError for SparseDataFrame (:issue:`887`)
  - Fix DataFrame.to_html encoding and columns (:issue:`890`, :issue:`891`, :issue:`909`)
  - Fix na-filling handling in mixed-type DataFrame (:issue:`910`)
  - Fix to DataFrame.set_value with non-existant row/col (:issue:`911`)
  - Fix malformed block in groupby when excluding nuisance columns (:issue:`916`)
  - Fix inconsistant NA handling in dtype=object arrays (:issue:`925`)
  - Fix missing center-of-mass computation in ewmcov (:issue:`862`)
  - Don't raise exception when opening read-only HDF5 file (:issue:`847`)
  - Fix possible out-of-bounds memory access in 0-length Series (:issue:`917`)



pandas 0.7.1
============

**Release date:** February 29, 2012

**New features / modules**

  - Add ``to_clipboard`` function to pandas namespace for writing objects to
    the system clipboard (:issue:`774`)
  - Add ``itertuples`` method to DataFrame for iterating through the rows of a
    dataframe as tuples (:issue:`818`)
  - Add ability to pass fill_value and method to DataFrame and Series align
    method (:issue:`806`, :issue:`807`)
  - Add fill_value option to reindex, align methods (:issue:`784`)
  - Enable concat to produce DataFrame from Series (:issue:`787`)
  - Add ``between`` method to Series (:issue:`802`)
  - Add HTML representation hook to DataFrame for the IPython HTML notebook
    (:issue:`773`)
  - Support for reading Excel 2007 XML documents using openpyxl

**Improvements to existing features**

  - Improve performance and memory usage of fillna on DataFrame
  - Can concatenate a list of Series along axis=1 to obtain a DataFrame (:issue:`787`)

**Bug fixes**

  - Fix memory leak when inserting large number of columns into a single
    DataFrame (:issue:`790`)
  - Appending length-0 DataFrame with new columns would not result in those new
    columns being part of the resulting concatenated DataFrame (:issue:`782`)
  - Fixed groupby corner case when passing dictionary grouper and as_index is
    False (:issue:`819`)
  - Fixed bug whereby bool array sometimes had object dtype (:issue:`820`)
  - Fix exception thrown on np.diff (:issue:`816`)
  - Fix to_records where columns are non-strings (:issue:`822`)
  - Fix Index.intersection where indices have incomparable types (:issue:`811`)
  - Fix ExcelFile throwing an exception for two-line file (:issue:`837`)
  - Add clearer error message in csv parser (:issue:`835`)
  - Fix loss of fractional seconds in HDFStore (:issue:`513`)
  - Fix DataFrame join where columns have datetimes (:issue:`787`)
  - Work around numpy performance issue in take (:issue:`817`)
  - Improve comparison operations for NA-friendliness (:issue:`801`)
  - Fix indexing operation for floating point values (:issue:`780`, :issue:`798`)
  - Fix groupby case resulting in malformed dataframe (:issue:`814`)
  - Fix behavior of reindex of Series dropping name (:issue:`812`)
  - Improve on redudant groupby computation (:issue:`775`)
  - Catch possible NA assignment to int/bool series with exception (:issue:`839`)



pandas 0.7.0
============

**Release date:** 2/9/2012

**New features / modules**

  - New ``merge`` function for efficiently performing full gamut of database /
    relational-algebra operations. Refactored existing join methods to use the
    new infrastructure, resulting in substantial performance gains (:issue:`220`,
    :issue:`249`, :issue:`267`)
  - New ``concat`` function for concatenating DataFrame or Panel objects along
    an axis. Can form union or intersection of the other axes. Improves
    performance of ``DataFrame.append`` (:issue:`468`, :issue:`479`, :issue:`273`)
  - Handle differently-indexed output values in ``DataFrame.apply`` (:issue:`498`)
  - Can pass list of dicts (e.g., a list of shallow JSON objects) to DataFrame
    constructor (:issue:`526`)
  - Add ``reorder_levels`` method to Series and DataFrame (:issue:`534`)
  - Add dict-like ``get`` function to DataFrame and Panel (:issue:`521`)
  - ``DataFrame.iterrows`` method for efficiently iterating through the rows of
    a DataFrame
  - Added ``DataFrame.to_panel`` with code adapted from ``LongPanel.to_long``
  - ``reindex_axis`` method added to DataFrame
  - Add ``level`` option to binary arithmetic functions on ``DataFrame`` and
    ``Series``
  - Add ``level`` option to the ``reindex`` and ``align`` methods on Series and
    DataFrame for broadcasting values across a level (:issue:`542`, :issue:`552`, others)
  - Add attribute-based item access to ``Panel`` and add IPython completion (PR
    :issue:`554`)
  - Add ``logy`` option to ``Series.plot`` for log-scaling on the Y axis
  - Add ``index``, ``header``, and ``justify`` options to
    ``DataFrame.to_string``. Add option to   (:issue:`570`, :issue:`571`)
  - Can pass multiple DataFrames to ``DataFrame.join`` to join on index (:issue:`115`)
  - Can pass multiple Panels to ``Panel.join`` (:issue:`115`)
  - Can pass multiple DataFrames to `DataFrame.append` to concatenate (stack)
    and multiple Series to ``Series.append`` too
  - Added ``justify`` argument to ``DataFrame.to_string`` to allow different
    alignment of column headers
  - Add ``sort`` option to GroupBy to allow disabling sorting of the group keys
    for potential speedups (:issue:`595`)
  - Can pass MaskedArray to Series constructor (:issue:`563`)
  - Add Panel item access via attributes and IPython completion (:issue:`554`)
  - Implement ``DataFrame.lookup``, fancy-indexing analogue for retrieving
    values given a sequence of row and column labels (:issue:`338`)
  - Add ``verbose`` option to ``read_csv`` and ``read_table`` to show number of
    NA values inserted in non-numeric columns (:issue:`614`)
  - Can pass a list of dicts or Series to ``DataFrame.append`` to concatenate
    multiple rows (:issue:`464`)
  - Add ``level`` argument to ``DataFrame.xs`` for selecting data from other
    MultiIndex levels. Can take one or more levels with potentially a tuple of
    keys for flexible retrieval of data (:issue:`371`, :issue:`629`)
  - New ``crosstab`` function for easily computing frequency tables (:issue:`170`)
  - Can pass a list of functions to aggregate with groupby on a DataFrame,
    yielding an aggregated result with hierarchical columns (:issue:`166`)
  - Add integer-indexing functions ``iget`` in Series and ``irow`` / ``iget``
    in DataFrame (:issue:`628`)
  - Add new ``Series.unique`` function, significantly faster than
    ``numpy.unique`` (:issue:`658`)
  - Add new ``cummin`` and ``cummax`` instance methods to ``Series`` and
    ``DataFrame`` (:issue:`647`)
  - Add new ``value_range`` function to return min/max of a dataframe (:issue:`288`)
  - Add ``drop`` parameter to ``reset_index`` method of ``DataFrame`` and added
    method to ``Series`` as well (:issue:`699`)
  - Add ``isin`` method to Index objects, works just like ``Series.isin`` (GH
    :issue:`657`)
  - Implement array interface on Panel so that ufuncs work (re: :issue:`740`)
  - Add ``sort`` option to ``DataFrame.join`` (:issue:`731`)
  - Improved handling of NAs (propagation) in binary operations with
    dtype=object arrays (:issue:`737`)
  - Add ``abs`` method to Pandas objects
  - Added ``algorithms`` module to start collecting central algos

**API Changes**

  - Label-indexing with integer indexes now raises KeyError if a label is not
    found instead of falling back on location-based indexing (:issue:`700`)
  - Label-based slicing via ``ix`` or ``[]`` on Series will now only work if
    exact matches for the labels are found or if the index is monotonic (for
    range selections)
  - Label-based slicing and sequences of labels can be passed to ``[]`` on a
    Series for both getting and setting (:issue:`86`)
  - `[]` operator (``__getitem__`` and ``__setitem__``) will raise KeyError
    with integer indexes when an index is not contained in the index. The prior
    behavior would fall back on position-based indexing if a key was not found
    in the index which would lead to subtle bugs. This is now consistent with
    the behavior of ``.ix`` on DataFrame and friends (:issue:`328`)
  - Rename ``DataFrame.delevel`` to ``DataFrame.reset_index`` and add
    deprecation warning
  - `Series.sort` (an in-place operation) called on a Series which is a view on
    a larger array (e.g. a column in a DataFrame) will generate an Exception to
    prevent accidentally modifying the data source (:issue:`316`)
  - Refactor to remove deprecated ``LongPanel`` class (:issue:`552`)
  - Deprecated ``Panel.to_long``, renamed to ``to_frame``
  - Deprecated ``colSpace`` argument in ``DataFrame.to_string``, renamed to
    ``col_space``
  - Rename ``precision`` to ``accuracy`` in engineering float formatter (GH
    :issue:`395`)
  - The default delimiter for ``read_csv`` is comma rather than letting
    ``csv.Sniffer`` infer it
  - Rename ``col_or_columns`` argument in ``DataFrame.drop_duplicates`` (GH
    :issue:`734`)

**Improvements to existing features**

  - Better error message in DataFrame constructor when passed column labels
    don't match data (:issue:`497`)
  - Substantially improve performance of multi-GroupBy aggregation when a
    Python function is passed, reuse ndarray object in Cython (:issue:`496`)
  - Can store objects indexed by tuples and floats in HDFStore (:issue:`492`)
  - Don't print length by default in Series.to_string, add `length` option (GH
    :issue:`489`)
  - Improve Cython code for multi-groupby to aggregate without having to sort
    the data (:issue:`93`)
  - Improve MultiIndex reindexing speed by storing tuples in the MultiIndex,
    test for backwards unpickling compatibility
  - Improve column reindexing performance by using specialized Cython take
    function
  - Further performance tweaking of Series.__getitem__ for standard use cases
  - Avoid Index dict creation in some cases (i.e. when getting slices, etc.),
    regression from prior versions
  - Friendlier error message in setup.py if NumPy not installed
  - Use common set of NA-handling operations (sum, mean, etc.) in Panel class
    also (:issue:`536`)
  - Default name assignment when calling ``reset_index`` on DataFrame with a
    regular (non-hierarchical) index (:issue:`476`)
  - Use Cythonized groupers when possible in Series/DataFrame stat ops with
    ``level`` parameter passed (:issue:`545`)
  - Ported skiplist data structure to C to speed up ``rolling_median`` by about
    5-10x in most typical use cases (:issue:`374`)
  - Some performance enhancements in constructing a Panel from a dict of
    DataFrame objects
  - Made ``Index._get_duplicates`` a public method by removing the underscore
  - Prettier printing of floats, and column spacing fix (:issue:`395`, :issue:`571`)
  - Add ``bold_rows`` option to DataFrame.to_html (:issue:`586`)
  - Improve the performance of ``DataFrame.sort_index`` by up to 5x or more
    when sorting by multiple columns
  - Substantially improve performance of DataFrame and Series constructors when
    passed a nested dict or dict, respectively (:issue:`540`, :issue:`621`)
  - Modified setup.py so that pip / setuptools will install dependencies (GH
    :issue:`507`, various pull requests)
  - Unstack called on DataFrame with non-MultiIndex will return Series (GH
    :issue:`477`)
  - Improve DataFrame.to_string and console formatting to be more consistent in
    the number of displayed digits (:issue:`395`)
  - Use bottleneck if available for performing NaN-friendly statistical
    operations that it implemented (:issue:`91`)
  - Monkey-patch context to traceback in ``DataFrame.apply`` to indicate which
    row/column the function application failed on (:issue:`614`)
  - Improved ability of read_table and read_clipboard to parse
    console-formatted DataFrames (can read the row of index names, etc.)
  - Can pass list of group labels (without having to convert to an ndarray
    yourself) to ``groupby`` in some cases (:issue:`659`)
  - Use ``kind`` argument to Series.order for selecting different sort kinds
    (:issue:`668`)
  - Add option to Series.to_csv to omit the index (:issue:`684`)
  - Add ``delimiter`` as an alternative to ``sep`` in ``read_csv`` and other
    parsing functions
  - Substantially improved performance of groupby on DataFrames with many
    columns by aggregating blocks of columns all at once (:issue:`745`)
  - Can pass a file handle or StringIO to Series/DataFrame.to_csv (:issue:`765`)
  - Can pass sequence of integers to DataFrame.irow(icol) and Series.iget, (GH
    :issue:`654`)
  - Prototypes for some vectorized string functions
  - Add float64 hash table to solve the Series.unique problem with NAs (:issue:`714`)
  - Memoize objects when reading from file to reduce memory footprint
  - Can get and set a column of a DataFrame with hierarchical columns
    containing "empty" ('') lower levels without passing the empty levels (PR
    :issue:`768`)

**Bug fixes**

  - Raise exception in out-of-bounds indexing of Series instead of
    seg-faulting, regression from earlier releases (:issue:`495`)
  - Fix error when joining DataFrames of different dtypes within the same
    typeclass (e.g. float32 and float64) (:issue:`486`)
  - Fix bug in Series.min/Series.max on objects like datetime.datetime (GH
    :issue:`487`)
  - Preserve index names in Index.union (:issue:`501`)
  - Fix bug in Index joining causing subclass information (like DateRange type)
    to be lost in some cases (:issue:`500`)
  - Accept empty list as input to DataFrame constructor, regression from 0.6.0
    (:issue:`491`)
  - Can output DataFrame and Series with ndarray objects in a dtype=object
    array (:issue:`490`)
  - Return empty string from Series.to_string when called on empty Series (GH
    :issue:`488`)
  - Fix exception passing empty list to DataFrame.from_records
  - Fix Index.format bug (excluding name field) with datetimes with time info
  - Fix scalar value access in Series to always return NumPy scalars,
    regression from prior versions (:issue:`510`)
  - Handle rows skipped at beginning of file in read_* functions (:issue:`505`)
  - Handle improper dtype casting in ``set_value`` methods
  - Unary '-' / __neg__ operator on DataFrame was returning integer values
  - Unbox 0-dim ndarrays from certain operators like all, any in Series
  - Fix handling of missing columns (was combine_first-specific) in
    DataFrame.combine for general case (:issue:`529`)
  - Fix type inference logic with boolean lists and arrays in DataFrame indexing
  - Use centered sum of squares in R-square computation if entity_effects=True
    in panel regression
  - Handle all NA case in Series.{corr, cov}, was raising exception (:issue:`548`)
  - Aggregating by multiple levels with ``level`` argument to DataFrame, Series
    stat method, was broken (:issue:`545`)
  - Fix Cython buf when converter passed to read_csv produced a numeric array
    (buffer dtype mismatch when passed to Cython type inference function) (GH
    :issue:`546`)
  - Fix exception when setting scalar value using .ix on a DataFrame with a
    MultiIndex (:issue:`551`)
  - Fix outer join between two DateRanges with different offsets that returned
    an invalid DateRange
  - Cleanup DataFrame.from_records failure where index argument is an integer
  - Fix Data.from_records failure when passed a dictionary
  - Fix NA handling in {Series, DataFrame}.rank with non-floating point dtypes
  - Fix bug related to integer type-checking in .ix-based indexing
  - Handle non-string index name passed to DataFrame.from_records
  - DataFrame.insert caused the columns name(s) field to be discarded (:issue:`527`)
  - Fix erroneous in monotonic many-to-one left joins
  - Fix DataFrame.to_string to remove extra column white space (:issue:`571`)
  - Format floats to default to same number of digits (:issue:`395`)
  - Added decorator to copy docstring from one function to another (:issue:`449`)
  - Fix error in monotonic many-to-one left joins
  - Fix __eq__ comparison between DateOffsets with different relativedelta
    keywords passed
  - Fix exception caused by parser converter returning strings (:issue:`583`)
  - Fix MultiIndex formatting bug with integer names (:issue:`601`)
  - Fix bug in handling of non-numeric aggregates in Series.groupby (:issue:`612`)
  - Fix TypeError with tuple subclasses (e.g. namedtuple) in
    DataFrame.from_records (:issue:`611`)
  - Catch misreported console size when running IPython within Emacs
  - Fix minor bug in pivot table margins, loss of index names and length-1
    'All' tuple in row labels
  - Add support for legacy WidePanel objects to be read from HDFStore
  - Fix out-of-bounds segfault in pad_object and backfill_object methods when
    either source or target array are empty
  - Could not create a new column in a DataFrame from a list of tuples
  - Fix bugs preventing SparseDataFrame and SparseSeries working with groupby
    (:issue:`666`)
  - Use sort kind in Series.sort / argsort (:issue:`668`)
  - Fix DataFrame operations on non-scalar, non-pandas objects (:issue:`672`)
  - Don't convert DataFrame column to integer type when passing integer to
    __setitem__ (:issue:`669`)
  - Fix downstream bug in pivot_table caused by integer level names in
    MultiIndex (:issue:`678`)
  - Fix SparseSeries.combine_first when passed a dense Series (:issue:`687`)
  - Fix performance regression in HDFStore loading when DataFrame or Panel
    stored in table format with datetimes
  - Raise Exception in DateRange when offset with n=0 is passed (:issue:`683`)
  - Fix get/set inconsistency with .ix property and integer location but
    non-integer index (:issue:`707`)
  - Use right dropna function for SparseSeries. Return dense Series for NA fill
    value (:issue:`730`)
  - Fix Index.format bug causing incorrectly string-formatted Series with
    datetime indexes (:issue:`726`, :issue:`758`)
  - Fix errors caused by object dtype arrays passed to ols (:issue:`759`)
  - Fix error where column names lost when passing list of labels to
    DataFrame.__getitem__, (:issue:`662`)
  - Fix error whereby top-level week iterator overwrote week instance
  - Fix circular reference causing memory leak in sparse array / series /
    frame, (:issue:`663`)
  - Fix integer-slicing from integers-as-floats (:issue:`670`)
  - Fix zero division errors in nanops from object dtype arrays in all NA case
    (:issue:`676`)
  - Fix csv encoding when using unicode (:issue:`705`, :issue:`717`, :issue:`738`)
  - Fix assumption that each object contains every unique block type in concat,
    (:issue:`708`)
  - Fix sortedness check of multiindex in to_panel (:issue:`719`, 720)
  - Fix that None was not treated as NA in PyObjectHashtable
  - Fix hashing dtype because of endianness confusion (:issue:`747`, :issue:`748`)
  - Fix SparseSeries.dropna to return dense Series in case of NA fill value (GH
    :issue:`730`)
  - Use map_infer instead of np.vectorize. handle NA sentinels if converter
    yields numeric array, (:issue:`753`)
  - Fixes and improvements to DataFrame.rank (:issue:`742`)
  - Fix catching AttributeError instead of NameError for bottleneck
  - Try to cast non-MultiIndex to better dtype when calling reset_index (:issue:`726`
    :issue:`440`)
  - Fix #1.QNAN0' float bug on 2.6/win64
  - Allow subclasses of dicts in DataFrame constructor, with tests
  - Fix problem whereby set_index destroys column multiindex (:issue:`764`)
  - Hack around bug in generating DateRange from naive DateOffset (:issue:`770`)
  - Fix bug in DateRange.intersection causing incorrect results with some
    overlapping ranges (:issue:`771`)

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
    instead of requiring that the indexes be exactly equal (:issue:`429`)

**New features / modules**

  - Can pass Series to DataFrame.append with ignore_index=True for appending a
    single row (:issue:`430`)
  - Add Spearman and Kendall correlation options to Series.corr and
    DataFrame.corr (:issue:`428`)
  - Add new `get_value` and `set_value` methods to Series, DataFrame, and Panel
    to very low-overhead access to scalar elements. df.get_value(row, column)
    is about 3x faster than df[column][row] by handling fewer cases (:issue:`437`,
    :issue:`438`). Add similar methods to sparse data structures for compatibility
  - Add Qt table widget to sandbox (:issue:`435`)
  - DataFrame.align can accept Series arguments, add axis keyword (:issue:`461`)
  - Implement new SparseList and SparseArray data structures. SparseSeries now
    derives from SparseArray (:issue:`463`)
  - max_columns / max_rows options in set_printoptions (:issue:`453`)
  - Implement Series.rank and DataFrame.rank, fast versions of
    scipy.stats.rankdata (:issue:`428`)
  - Implement DataFrame.from_items alternate constructor (:issue:`444`)
  - DataFrame.convert_objects method for inferring better dtypes for object
    columns (:issue:`302`)
  - Add rolling_corr_pairwise function for computing Panel of correlation
    matrices (:issue:`189`)
  - Add `margins` option to `pivot_table` for computing subgroup aggregates (GH
    :issue:`114`)
  - Add `Series.from_csv` function (:issue:`482`)

**Improvements to existing features**

  - Improve memory usage of `DataFrame.describe` (do not copy data
    unnecessarily) (:issue:`425`)
  - Use same formatting function for outputting floating point Series to console
    as in DataFrame (:issue:`420`)
  - DataFrame.delevel will try to infer better dtype for new columns (:issue:`440`)
  - Exclude non-numeric types in DataFrame.{corr, cov}
  - Override Index.astype to enable dtype casting (:issue:`412`)
  - Use same float formatting function for Series.__repr__ (:issue:`420`)
  - Use available console width to output DataFrame columns (:issue:`453`)
  - Accept ndarrays when setting items in Panel (:issue:`452`)
  - Infer console width when printing __repr__ of DataFrame to console (PR
    :issue:`453`)
  - Optimize scalar value lookups in the general case by 25% or more in Series
    and DataFrame
  - Can pass DataFrame/DataFrame and DataFrame/Series to
    rolling_corr/rolling_cov (:issue:`462`)
  - Fix performance regression in cross-sectional count in DataFrame, affecting
    DataFrame.dropna speed
  - Column deletion in DataFrame copies no data (computes views on blocks) (GH
    :issue:`158`)
  - MultiIndex.get_level_values can take the level name
  - More helpful error message when DataFrame.plot fails on one of the columns
    (:issue:`478`)
  - Improve performance of DataFrame.{index, columns} attribute lookup

**Bug fixes**

  - Fix O(K^2) memory leak caused by inserting many columns without
    consolidating, had been present since 0.4.0 (:issue:`467`)
  - `DataFrame.count` should return Series with zero instead of NA with length-0
    axis (:issue:`423`)
  - Fix Yahoo! Finance API usage in pandas.io.data (:issue:`419`, :issue:`427`)
  - Fix upstream bug causing failure in Series.align with empty Series (:issue:`434`)
  - Function passed to DataFrame.apply can return a list, as long as it's the
    right length. Regression from 0.4 (:issue:`432`)
  - Don't "accidentally" upcast scalar values when indexing using .ix (:issue:`431`)
  - Fix groupby exception raised with as_index=False and single column selected
    (:issue:`421`)
  - Implement DateOffset.__ne__ causing downstream bug (:issue:`456`)
  - Fix __doc__-related issue when converting py -> pyo with py2exe
  - Bug fix in left join Cython code with duplicate monotonic labels
  - Fix bug when unstacking multiple levels described in :issue:`451`
  - Exclude NA values in dtype=object arrays, regression from 0.5.0 (:issue:`469`)
  - Use Cython map_infer function in DataFrame.applymap to properly infer
    output type, handle tuple return values and other things that were breaking
    (:issue:`465`)
  - Handle floating point index values in HDFStore (:issue:`454`)
  - Fixed stale column reference bug (cached Series object) caused by type
    change / item deletion in DataFrame (:issue:`473`)
  - Index.get_loc should always raise Exception when there are duplicates
  - Handle differently-indexed Series input to DataFrame constructor (:issue:`475`)
  - Omit nuisance columns in multi-groupby with Python function
  - Buglet in handling of single grouping in general apply
  - Handle type inference properly when passing list of lists or tuples to
    DataFrame constructor (:issue:`484`)
  - Preserve Index / MultiIndex names in GroupBy.apply concatenation step (GH
    :issue:`481`)

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
    default instead of excluding them (:issue:`382`)

**New features / modules**

  - Add `melt` function to `pandas.core.reshape`
  - Add `level` parameter to group by level in Series and DataFrame
    descriptive statistics (:issue:`313`)
  - Add `head` and `tail` methods to Series, analogous to to DataFrame (PR
    :issue:`296`)
  - Add `Series.isin` function which checks if each value is contained in a
    passed sequence (:issue:`289`)
  - Add `float_format` option to `Series.to_string`
  - Add `skip_footer` (:issue:`291`) and `converters` (:issue:`343`) options to
    `read_csv` and `read_table`
  - Add proper, tested weighted least squares to standard and panel OLS (GH
    :issue:`303`)
  - Add `drop_duplicates` and `duplicated` functions for removing duplicate
    DataFrame rows and checking for duplicate rows, respectively (:issue:`319`)
  - Implement logical (boolean) operators ``&``, ``|``, ``^`` on DataFrame
    (:issue:`347`)
  - Add `Series.mad`, mean absolute deviation, matching DataFrame
  - Add `QuarterEnd` DateOffset (:issue:`321`)
  - Add matrix multiplication function `dot` to DataFrame (:issue:`65`)
  - Add `orient` option to `Panel.from_dict` to ease creation of mixed-type
    Panels (:issue:`359`, :issue:`301`)
  - Add `DataFrame.from_dict` with similar `orient` option
  - Can now pass list of tuples or list of lists to `DataFrame.from_records`
    for fast conversion to DataFrame (:issue:`357`)
  - Can pass multiple levels to groupby, e.g. `df.groupby(level=[0, 1])` (GH
    :issue:`103`)
  - Can sort by multiple columns in `DataFrame.sort_index` (:issue:`92`, :issue:`362`)
  - Add fast `get_value` and `put_value` methods to DataFrame and
    micro-performance tweaks (:issue:`360`)
  - Add `cov` instance methods to Series and DataFrame (:issue:`194`, :issue:`362`)
  - Add bar plot option to `DataFrame.plot` (:issue:`348`)
  - Add `idxmin` and `idxmax` functions to Series and DataFrame for computing
    index labels achieving maximum and minimum values (:issue:`286`)
  - Add `read_clipboard` function for parsing DataFrame from OS clipboard,
    should work across platforms (:issue:`300`)
  - Add `nunique` function to Series for counting unique elements (:issue:`297`)
  - DataFrame constructor will use Series name if no columns passed (:issue:`373`)
  - Support regular expressions and longer delimiters in read_table/read_csv,
    but does not handle quoted strings yet (:issue:`364`)
  - Add `DataFrame.to_html` for formatting DataFrame to HTML (:issue:`387`)
  - MaskedArray can be passed to DataFrame constructor and masked values will be
    converted to NaN (:issue:`396`)
  - Add `DataFrame.boxplot` function (:issue:`368`, others)
  - Can pass extra args, kwds to DataFrame.apply (:issue:`376`)

**Improvements to existing features**

  - Raise more helpful exception if date parsing fails in DateRange (:issue:`298`)
  - Vastly improved performance of GroupBy on axes with a MultiIndex (:issue:`299`)
  - Print level names in hierarchical index in Series repr (:issue:`305`)
  - Return DataFrame when performing GroupBy on selected column and
    as_index=False (:issue:`308`)
  - Can pass vector to `on` argument in `DataFrame.join` (:issue:`312`)
  - Don't show Series name if it's None in the repr, also omit length for short
    Series (:issue:`317`)
  - Show legend by default in `DataFrame.plot`, add `legend` boolean flag (GH
    :issue:`324`)
  - Significantly improved performance of `Series.order`, which also makes
    np.unique called on a Series faster (:issue:`327`)
  - Faster cythonized count by level in Series and DataFrame (:issue:`341`)
  - Raise exception if dateutil 2.0 installed on Python 2.x runtime (:issue:`346`)
  - Significant GroupBy performance enhancement with multiple keys with many
    "empty" combinations
  - New Cython vectorized function `map_infer` speeds up `Series.apply` and
    `Series.map` significantly when passed elementwise Python function,
    motivated by :issue:`355`
  - Cythonized `cache_readonly`, resulting in substantial micro-performance
    enhancements throughout the codebase (:issue:`361`)
  - Special Cython matrix iterator for applying arbitrary reduction operations
    with 3-5x better performance than `np.apply_along_axis` (:issue:`309`)
  - Add `raw` option to `DataFrame.apply` for getting better performance when
    the passed function only requires an ndarray (:issue:`309`)
  - Improve performance of `MultiIndex.from_tuples`
  - Can pass multiple levels to `stack` and `unstack` (:issue:`370`)
  - Can pass multiple values columns to `pivot_table` (:issue:`381`)
  - Can call `DataFrame.delevel` with standard Index with name set (:issue:`393`)
  - Use Series name in GroupBy for result index (:issue:`363`)
  - Refactor Series/DataFrame stat methods to use common set of NaN-friendly
    function
  - Handle NumPy scalar integers at C level in Cython conversion routines

**Bug fixes**

  - Fix bug in `DataFrame.to_csv` when writing a DataFrame with an index
    name (:issue:`290`)
  - DataFrame should clear its Series caches on consolidation, was causing
    "stale" Series to be returned in some corner cases (:issue:`304`)
  - DataFrame constructor failed if a column had a list of tuples (:issue:`293`)
  - Ensure that `Series.apply` always returns a Series and implement
    `Series.round` (:issue:`314`)
  - Support boolean columns in Cythonized groupby functions (:issue:`315`)
  - `DataFrame.describe` should not fail if there are no numeric columns,
    instead return categorical describe (:issue:`323`)
  - Fixed bug which could cause columns to be printed in wrong order in
    `DataFrame.to_string` if specific list of columns passed (:issue:`325`)
  - Fix legend plotting failure if DataFrame columns are integers (:issue:`326`)
  - Shift start date back by one month for Yahoo! Finance API in pandas.io.data
    (:issue:`329`)
  - Fix `DataFrame.join` failure on unconsolidated inputs (:issue:`331`)
  - DataFrame.min/max will no longer fail on mixed-type DataFrame (:issue:`337`)
  - Fix `read_csv` / `read_table` failure when passing list to index_col that is
    not in ascending order (:issue:`349`)
  - Fix failure passing Int64Index to Index.union when both are monotonic
  - Fix error when passing SparseSeries to (dense) DataFrame constructor
  - Added missing bang at top of setup.py (:issue:`352`)
  - Change `is_monotonic` on MultiIndex so it properly compares the tuples
  - Fix MultiIndex outer join logic (:issue:`351`)
  - Set index name attribute with single-key groupby (:issue:`358`)
  - Bug fix in reflexive binary addition in Series and DataFrame for
    non-commutative operations (like string concatenation) (:issue:`353`)
  - setupegg.py will invoke Cython (:issue:`192`)
  - Fix block consolidation bug after inserting column into MultiIndex (:issue:`366`)
  - Fix bug in join operations between Index and Int64Index (:issue:`367`)
  - Handle min_periods=0 case in moving window functions (:issue:`365`)
  - Fixed corner cases in DataFrame.apply/pivot with empty DataFrame (:issue:`378`)
  - Fixed repr exception when Series name is a tuple
  - Always return DateRange from `asfreq` (:issue:`390`)
  - Pass level names to `swaplavel` (:issue:`379`)
  - Don't lose index names in `MultiIndex.droplevel` (:issue:`394`)
  - Infer more proper return type in `DataFrame.apply` when no columns or rows
    depending on whether the passed function is a reduction (:issue:`389`)
  - Always return NA/NaN from Series.min/max and DataFrame.min/max when all of a
    row/column/values are NA (:issue:`384`)
  - Enable partial setting with .ix / advanced indexing (:issue:`397`)
  - Handle mixed-type DataFrames correctly in unstack, do not lose type
    information (:issue:`403`)
  - Fix integer name formatting bug in Index.format and in Series.__repr__
  - Handle label types other than string passed to groupby (:issue:`405`)
  - Fix bug in .ix-based indexing with partial retrieval when a label is not
    contained in a level
  - Index name was not being pickled (:issue:`408`)
  - Level name should be passed to result index in GroupBy.apply (:issue:`416`)

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
    :issue:`225`)
  - Removed `weights` option in panel regression which was not doing anything
    principled (:issue:`155`)
  - Changed `buffer` argument name in `Series.to_string` to `buf`
  - `Series.to_string` and `DataFrame.to_string` now return strings by default
    instead of printing to sys.stdout
  - Deprecated `nanRep` argument in various `to_string` and `to_csv` functions
    in favor of `na_rep`. Will be removed in 0.6 (:issue:`275`)
  - Renamed `delimiter` to `sep` in `DataFrame.from_csv` for consistency
  - Changed order of `Series.clip` arguments to match those of `numpy.clip` and
    added (unimplemented) `out` argument so `numpy.clip` can be called on a
    Series (:issue:`272`)
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
    lazily iterating through chunks of a flat file (:issue:`242`)
  - Added ability to join on multiple columns in `DataFrame.join` (:issue:`214`)
  - Added private `_get_duplicates` function to `Index` for identifying
    duplicate values more easily
  - Added column attribute access to DataFrame, e.g. df.A equivalent to df['A']
    if 'A' is a column in the DataFrame (:issue:`213`)
  - Added IPython tab completion hook for DataFrame columns. (:issue:`233`, :issue:`230`)
  - Implement `Series.describe` for Series containing objects (:issue:`241`)
  - Add inner join option to `DataFrame.join` when joining on key(s) (:issue:`248`)
  - Can select set of DataFrame columns by passing a list to `__getitem__` (GH
    :issue:`253`)
  - Can use & and | to intersection / union Index objects, respectively (GH
    :issue:`261`)
  - Added `pivot_table` convenience function to pandas namespace (:issue:`234`)
  - Implemented `Panel.rename_axis` function (:issue:`243`)
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
    performance (:issue:`211`)
  - Improved speed of `DataFrame.xs` on mixed-type DataFrame objects by about
    5x, regression from 0.3.0 (:issue:`215`)
  - With new `DataFrame.align` method, speeding up binary operations between
    differently-indexed DataFrame objects by 10-25%.
  - Significantly sped up conversion of nested dict into DataFrame (:issue:`212`)
  - Can pass hierarchical index level name to `groupby` instead of the level
    number if desired (:issue:`223`)
  - Add support for different delimiters in `DataFrame.to_csv` (:issue:`244`)
  - Add more helpful error message when importing pandas post-installation from
    the source directory (:issue:`250`)
  - Significantly speed up DataFrame `__repr__` and `count` on large mixed-type
    DataFrame objects
  - Better handling of pyx file dependencies in Cython module build (:issue:`271`)

**Bug fixes**

  - `read_csv` / `read_table` fixes

    - Be less aggressive about converting float->int in cases of floating point
      representations of integers like 1.0, 2.0, etc.
    - "True"/"False" will not get correctly converted to boolean
    - Index name attribute will get set when specifying an index column
    - Passing column names should force `header=None` (:issue:`257`)
    - Don't modify passed column names when `index_col` is not None
      (:issue:`258`)
    - Can sniff CSV separator in zip file (since seek is not supported, was
      failing before)

  - Worked around matplotlib "bug" in which series[:, np.newaxis] fails. Should
    be reported upstream to matplotlib (:issue:`224`)
  - DataFrame.iteritems was not returning Series with the name attribute
    set. Also neither was DataFrame._series
  - Can store datetime.date objects in HDFStore (:issue:`231`)
  - Index and Series names are now stored in HDFStore
  - Fixed problem in which data would get upcasted to object dtype in
    GroupBy.apply operations (:issue:`237`)
  - Fixed outer join bug with empty DataFrame (:issue:`238`)
  - Can create empty Panel (:issue:`239`)
  - Fix join on single key when passing list with 1 entry (:issue:`246`)
  - Don't raise Exception on plotting DataFrame with an all-NA column (:issue:`251`,
    :issue:`254`)
  - Bug min/max errors when called on integer DataFrames (:issue:`241`)
  - `DataFrame.iteritems` and `DataFrame._series` not assigning name attribute
  - Panel.__repr__ raised exception on length-0 major/minor axes
  - `DataFrame.join` on key with empty DataFrame produced incorrect columns
  - Implemented `MultiIndex.diff` (:issue:`260`)
  - `Int64Index.take` and `MultiIndex.take` lost name field, fix downstream
    issue :issue:`262`
  - Can pass list of tuples to `Series` (:issue:`270`)
  - Can pass level name to `DataFrame.stack`
  - Support set operations between MultiIndex and Index
  - Fix many corner cases in MultiIndex set operations
    - Fix MultiIndex-handling bug with GroupBy.apply when returned groups are not
    indexed the same
  - Fix corner case bugs in DataFrame.apply
  - Setting DataFrame index did not cause Series cache to get cleared
  - Various int32 -> int64 platform-specific issues
  - Don't be too aggressive converting to integer when parsing file with
    MultiIndex (:issue:`285`)
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

  - Python 3 support using 2to3 (:issue:`200`, Thomas Kluyver)
  - Add `name` attribute to `Series` and added relevant logic and tests. Name
    now prints as part of `Series.__repr__`
  - Add `name` attribute to standard Index so that stacking / unstacking does
    not discard names and so that indexed DataFrame objects can be reliably
    round-tripped to flat files, pickle, HDF5, etc.
  - Add `isnull` and `notnull` as instance methods on Series (:issue:`209`, :issue:`203`)

**Improvements to existing features**

  - Skip xlrd-related unit tests if not installed
  - `Index.append` and `MultiIndex.append` can accept a list of Index objects to
    concatenate together
  - Altered binary operations on differently-indexed SparseSeries objects to use
    the integer-based (dense) alignment logic which is faster with a larger
    number of blocks (:issue:`205`)
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
  - `MultiIndex.sortlevel` discarded the level names (:issue:`202`)
  - Fix bugs in groupby, join, and append due to improper concatenation of
    `MultiIndex` objects (:issue:`201`)
  - Fix regression from 0.4.1, `isnull` and `notnull` ceased to work on other
    kinds of Python scalar objects like `datetime.datetime`
  - Raise more helpful exception when attempting to write empty DataFrame or
    LongPanel to `HDFStore` (:issue:`204`)
  - Use stdlib csv module to properly escape strings with commas in
    `DataFrame.to_csv` (:issue:`206`, Thomas Kluyver)
  - Fix Python ndarray access in Cython code for sparse blocked index integrity
    check
  - Fix bug writing Series to CSV in Python 3 (:issue:`209`)
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
    (:issue:`187`)
  - Wrote templating / code generation script to auto-generate Cython code for
    various functions which need to be available for the 4 major data types
    used in pandas (float64, bool, object, int64)
  - Refactored code related to `DataFrame.join` so that intermediate aligned
    copies of the data in each `DataFrame` argument do not need to be
    created. Substantial performance increases result (:issue:`176`)
  - Substantially improved performance of generic `Index.intersection` and
    `Index.union`
  - Improved performance of `DateRange.union` with overlapping ranges and
    non-cacheable offsets (like Minute). Implemented analogous fast
    `DateRange.intersection` for overlapping ranges.
  - Implemented `BlockManager.take` resulting in significantly faster `take`
    performance on mixed-type `DataFrame` objects (:issue:`104`)
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
  - Throw exception when step specified in label-based slice (:issue:`185`)
  - Fix isnull to correctly work with np.float32. Fix upstream bug described in
    :issue:`182`
  - Finish implementation of as_index=False in groupby for DataFrame
    aggregation (:issue:`181`)
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
    objects has been implemented (fixes :issue:`135`)
  - `read_csv` can read multiple columns into a `MultiIndex`. DataFrame's
    `to_csv` method will properly write out a `MultiIndex` which can be read
    back (:issue:`151`, thanks to Skipper Seabold)
  - Wrote fast time series merging / joining methods in Cython. Will be
    integrated later into DataFrame.join and related functions
  - Added `ignore_index` option to `DataFrame.append` for combining unindexed
    records stored in a DataFrame

**Improvements to existing features**

  - Some speed enhancements with internal Index type-checking function
  - `DataFrame.rename` has a new `copy` parameter which can rename a DataFrame
    in place
  - Enable unstacking by level name (:issue:`142`)
  - Enable sortlevel to work by level name (:issue:`141`)
  - `read_csv` can automatically "sniff" other kinds of delimiters using
    `csv.Sniffer` (:issue:`146`)
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
  - Fixed single-key groupby on DataFrame with as_index=False (:issue:`160`)
  - `Series.shift` was failing on integer Series (:issue:`154`)
  - `unstack` methods were producing incorrect output in the case of duplicate
    hierarchical labels. An exception will now be raised (:issue:`147`)
  - Calling `count` with level argument caused reduceat failure or segfault in
    earlier NumPy (:issue:`169`)
  - Fixed `DataFrame.corrwith` to automatically exclude non-numeric data (GH
    :issue:`144`)
  - Unicode handling bug fixes in `DataFrame.to_string` (:issue:`138`)
  - Excluding OLS degenerate unit test case that was causing platform specific
    failure (:issue:`149`)
  - Skip blosc-dependent unit tests for PyTables < 2.2 (:issue:`137`)
  - Calling `copy` on `DateRange` did not copy over attributes to the new object
    (:issue:`168`)
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

  - Fixed bug in IndexableSkiplist Cython code that was breaking rolling_max
    function
  - Numerous numpy.int64-related indexing fixes
  - Several NumPy 1.4.0 NaN-handling fixes
  - Bug fixes to pandas.io.parsers.parseCSV
  - Fixed `DateRange` caching issue with unusual date offsets
  - Fixed bug in `DateRange.union`
  - Fixed corner case in `IndexableSkiplist` implementation
