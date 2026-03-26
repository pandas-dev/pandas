{{ header }}

.. _api.typing.aliases:

{{ header }}

======================================
pandas typing aliases
======================================

**************
Typing aliases
**************

.. module:: pandas.api.typing.aliases

.. raw:: html

    <style>
    td > dl.py.type { margin-bottom: 0; }
    td > dl.py.type .descclassname { display: none; }
    </style>

The typing declarations in ``pandas/_typing.py`` are considered private, and used
by pandas developers for type checking of the pandas code base.  For users, it is
highly recommended to use the ``pandas-stubs`` package that represents the officially
supported type declarations for users of pandas.
They are documented here for users who wish to use these declarations in their
own python code that calls pandas or expects certain results.

.. warning::

    Note that the definitions and use cases of these aliases are subject to change without notice in any major, minor, or patch release of pandas.

Each of these aliases listed in the table below can be found by importing them from :py:mod:`pandas.api.typing.aliases`.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Alias
     - Meaning
   * - .. type:: AggFuncType
     - Type of functions that can be passed to :meth:`DataFrame's <pandas.DataFrame.agg>`, :meth:`Series' <pandas.Series.agg>`, and :meth:`DataFrameGroupBy's <pandas.api.typing.DataFrameGroupBy.aggregate>` ``aggregate()`` methods
   * - .. type:: AlignJoin
     - Argument type for ``join`` in :meth:`DataFrame's <pandas.DataFrame.align>` and :meth:`Series' <pandas.Series.align>` ``align()`` methods
   * - .. type:: AnyAll
     - Argument type for ``how`` in :meth:`DataFrame's <pandas.DataFrame.dropna>` and :meth:`Series' <pandas.Series.dropna>` ``dropna()`` methods
   * - .. type:: AnyArrayLike
     - Used to represent :class:`~pandas.api.extensions.ExtensionArray`, ``numpy`` arrays, :class:`~pandas.Index` and :class:`~pandas.Series`
   * - .. type:: ArrayLike
     - Used to represent :class:`~pandas.api.extensions.ExtensionArray`, ``numpy`` arrays
   * - .. type:: AstypeArg
     - Argument type in :meth:`DataFrame's <pandas.DataFrame.astype>` and :meth:`Series' <pandas.Series.astype>` ``astype()`` methods
   * - .. type:: Axes
     - :py:type:`AnyArrayLike` plus sequences (not strings) and ``range``
   * - .. type:: Axis
     - Argument type for ``axis`` in many methods
   * - .. type:: CSVEngine
     - Argument type for ``engine`` in :meth:`pandas.read_csv`
   * - .. type:: ColspaceArgType
     - Argument type for ``colspace`` in :meth:`pandas.DataFrame.to_html`
   * - .. type:: CompressionOptions
     - Argument type for ``compression`` in all I/O output methods except :meth:`pandas.DataFrame.to_parquet`
   * - .. type:: CorrelationMethod
     - Argument type for ``correlation`` in :meth:`DataFrame's <pandas.DataFrame.corr>` and :meth:`Series' <pandas.Series.corr>` ``corr()`` methods
   * - .. type:: DropKeep
     - Argument type for ``keep`` in :meth:`DataFrame's <pandas.DataFrame.drop_duplicates>` and :meth:`Series' <pandas.Series.drop_duplicates>` ``drop_duplicates()`` methods
   * - .. type:: Dtype
     - Types as objects that can be used to specify dtypes
   * - .. type:: DtypeArg
     - Argument type for ``dtype`` in various methods
   * - .. type:: DtypeBackend
     - Argument type for ``dtype_backend`` in various methods
   * - .. type:: DtypeObj
     - Numpy dtypes and Extension dtypes
   * - .. type:: ExcelWriterIfSheetExists
     - Argument type for ``if_sheet_exists`` in :class:`~pandas.ExcelWriter`
   * - .. type:: ExcelWriterMergeCells
     - Argument type for ``merge_cells`` in :meth:`DataFrame's <pandas.DataFrame.to_excel>` and :meth:`Series' <pandas.Series.to_excel>` ``to_excel()`` methods
   * - .. type:: FilePath
     - Type of paths for files for I/O methods
   * - .. type:: FillnaOptions
     - Argument type for ``method`` in various methods where NA values are filled
   * - .. type:: FloatFormatType
     - Argument type for ``float_format`` in :meth:`DataFrame's <pandas.DataFrame.to_string>` and :meth:`Series' <pandas.Series.to_string>` ``to_string()`` methods
   * - .. type:: FormattersType
     - Argument type for ``formatters`` in :meth:`DataFrame's <pandas.DataFrame.to_string>` and :meth:`Series' <pandas.Series.to_string>` ``to_string()`` methods
   * - .. type:: FromDictOrient
     - Argument type for ``orient`` in :meth:`DataFrame.from_dict() <pandas.DataFrame.from_dict>`
   * - .. type:: HTMLFlavors
     - Argument type for ``flavor`` in :meth:`pandas.read_html`
   * - .. type:: IgnoreRaise
     - Argument type for ``errors`` in multiple methods
   * - .. type:: IndexLabel
     - Argument type for ``level`` in multiple methods
   * - .. type:: InterpolateOptions
     - Argument type for ``interpolate`` in :meth:`DataFrame's <pandas.DataFrame.interpolate>` and :meth:`Series' <pandas.Series.interpolate>` ``interpolate()`` methods
   * - .. type:: IntervalClosedType
     - Argument type for ``closed`` in :class:`~pandas.Interval`, :class:`~pandas.IntervalIndex`, and ``inclusive`` in various methods
   * - .. type:: IntervalLeftRight
     - Restriction for ``closed`` to be ``left`` or ``right`` in :class:`~pandas.Interval`, :class:`~pandas.IntervalIndex`, and ``inclusive`` in various methods
   * - .. type:: JSONEngine
     - Argument type for ``engine`` in :meth:`pandas.read_json`
   * - .. type:: JSONSerializable
     - Argument type for the return type of a callable for argument ``default_handler`` in :meth:`DataFrame's <pandas.DataFrame.to_json>` and :meth:`Series' <pandas.Series.to_json>` ``to_json()`` methods
   * - .. type:: JoinHow
     - Argument type for ``how`` in :meth:`pandas.merge_ordered` and for ``join`` in :meth:`Series.align() <pandas.Series.align>`
   * - .. type:: JoinValidate
     - Argument type for ``validate`` in :meth:`DataFrame.join() <pandas.DataFrame.join>`
   * - .. type:: ListLike
     - Argument type for arguments that can be either a single value or a list of values in various methods
   * - .. type:: MergeHow
     - Argument type for ``how`` in :meth:`pandas.merge`
   * - .. type:: MergeValidate
     - Argument type for ``validate`` in :meth:`pandas.merge`
   * - .. type:: NaPosition
     - Argument type for ``na_position`` in :meth:`DataFrame's <pandas.DataFrame.sort_values>` and :meth:`Series' <pandas.Series.sort_values>` ``sort_values()`` methods
   * - .. type:: NsmallestNlargestKeep
     - Argument type for ``keep`` in :meth:`DataFrame's <pandas.DataFrame.nlargest>` and :meth:`Series' <pandas.Series.nlargest>` ``nlargest()``, :meth:`DataFrame's <pandas.DataFrame.nsmallest>` and :meth:`Series' <pandas.Series.nsmallest>` ``nsmallest()``, and :meth:`SeriesGroupBy's <pandas.api.typing.SeriesGroupBy.nlargest>` ``nlargest()`` methods
   * - .. type:: OpenFileErrors
     - Argument type for ``errors`` in :meth:`DataFrame's <pandas.DataFrame.to_hdf>`, :meth:`Series' <pandas.Series.to_hdf>` ``to_hdf()`` methods, and :meth:`DataFrame's <pandas.DataFrame.to_csv>` and :meth:`Series' <pandas.Series.to_csv>` ``to_csv()`` methods
   * - .. type:: Ordered
     - Return type for :attr:`ordered <pandas.CategoricalDtype.ordered>` in :class:`pandas.CategoricalDtype` and :class:`pandas.Categorical`
   * - .. type:: ParquetCompressionOptions
     - Argument type for ``compression`` in :meth:`DataFrame.to_parquet() <pandas.DataFrame.to_parquet>`
   * - .. type:: QuantileInterpolation
     - Argument type for ``interpolation`` in :meth:`DataFrame's <pandas.DataFrame.quantile>` and :meth:`Series' <pandas.Series.quantile>` ``quantile()`` methods
   * - .. type:: ReadBuffer
     - Additional argument type corresponding to buffers for various file reading methods
   * - .. type:: ReadCsvBuffer
     - Additional argument type corresponding to buffers for :meth:`pandas.read_csv`
   * - .. type:: ReadPickleBuffer
     - Additional argument type corresponding to buffers for :meth:`pandas.read_pickle`
   * - .. type:: ReindexMethod
     - Argument type for ``reindex`` in :meth:`DataFrame's <pandas.DataFrame.reindex>` and :meth:`Series' <pandas.Series.reindex>` ``reindex()`` methods
   * - .. type:: Scalar
     - Types that can be stored in :class:`~pandas.Series` with non-object dtype
   * - .. type:: ScalarIndexer
     - Argument type used for scalar indexing operations, such as the ``key`` argument in ``__getitem__()`` methods
   * - .. type:: SequenceIndexer
     - Argument type used for sequence indexing operations, such as the ``key`` argument in ``__getitem__()`` methods
   * - .. type:: SequenceNotStr
     - Used for arguments that require sequences, but not plain strings
   * - .. type:: SliceType
     - Argument types for ``start`` and ``end`` in :meth:`Index.slice_locs() <pandas.Index.slice_locs>`
   * - .. type:: SortKind
     - Argument type for ``kind`` in :meth:`DataFrame's <pandas.DataFrame.sort_values>` and :meth:`Series' <pandas.Series.sort_values>` ``sort_values()`` methods
   * - .. type:: StorageOptions
     - Argument type for ``storage_options`` in various file output methods
   * - .. type:: Suffixes
     - Argument type for ``suffixes`` in :meth:`pandas.merge`, :meth:`pandas.merge_ordered`,  and :meth:`DataFrame's <pandas.DataFrame.compare>` and :meth:`Series' <pandas.Series.compare>` ``compare()`` methods
   * - .. type:: TakeIndexer
     - Argument type for ``indexer`` and ``indices`` in :meth:`DataFrame's <pandas.DataFrame.take>` and :meth:`Series' <pandas.Series.take>` ``take()`` methods
   * - .. type:: TimeAmbiguous
     - Argument type for ``ambiguous`` in time operations
   * - .. type:: TimeGrouperOrigin
     - Argument type for ``origin`` in :meth:`DataFrame's <pandas.DataFrame.resample>`, :meth:`Series' <pandas.Series.resample>` ``resample()`` methods and for :class:`~pandas.Grouper`
   * - .. type:: TimeNonexistent
     - Argument type for ``nonexistent`` in time operations
   * - .. type:: TimeUnit
     - Time unit argument and return type for :attr:`pandas.Timedelta.unit`, arguments ``unit`` and ``date_unit``
   * - .. type:: TimedeltaConvertibleTypes
     - Argument type for ``offset`` in various methods, such as :meth:`DataFrame's <pandas.DataFrame.resample>` and :meth:`Series' <pandas.Series.resample>` ``resample()``, ``halflife`` in :meth:`DataFrame's <pandas.DataFrame.ewm>`, :meth:`DataFrameGroupBy's <pandas.api.typing.DataFrameGroupBy.ewm>`, and :meth:`Series' <pandas.Series.ewm>` ``ewm()``, and ``start`` and ``end`` in :meth:`pandas.timedelta_range`
   * - .. type:: TimestampConvertibleTypes
     - Argument type for ``origin`` in :meth:`DataFrame's <pandas.DataFrame.resample>` and :meth:`Series' <pandas.Series.resample>` ``resample()``, and in :meth:`pandas.to_datetime`
   * - .. type:: ToStataByteorder
     - Argument type for ``byteorder`` in :meth:`DataFrame.to_stata() <pandas.DataFrame.to_stata>`
   * - .. type:: ToTimestampHow
     - Argument type for ``how`` in :meth:`DataFrame's <pandas.DataFrame.to_timestamp>` and :meth:`Series' <pandas.Series.to_timestamp>` ``to_timestamp()`` methods, and ``convention`` in :meth:`DataFrame's <pandas.DataFrame.resample>` and :meth:`Series' <pandas.Series.resample>` ``resample()`` methods
   * - .. type:: UpdateJoin
     - Argument type for ``join`` in :meth:`DataFrame.update() <pandas.DataFrame.update>`
   * - .. type:: UsecolsArgType
     - Argument type for ``usecols`` in :meth:`pandas.read_clipboard`, :meth:`pandas.read_csv` and :meth:`pandas.read_excel`
   * - .. type:: WindowingRankType
     - Argument type for ``method`` in :meth:`Rolling's <pandas.api.typing.Rolling.rank>` and :meth:`Expanding's <pandas.api.typing.Expanding.rank>` ``rank()`` methods, applicable in rolling and expanding window operations
   * - .. type:: WriteBuffer
     - Additional argument type corresponding to buffers for various file output methods
   * - .. type:: WriteExcelBuffer
     - Additional argument type corresponding to buffers for :meth:`DataFrame's <pandas.DataFrame.to_excel>`, :meth:`Series' <pandas.Series.to_excel>`  and :meth:`Styler's <pandas.io.formats.style.Styler.to_excel>` ``to_excel()`` methods
   * - .. type:: XMLParsers
     - Argument type for ``parser`` in :meth:`DataFrame.to_xml() <pandas.DataFrame.to_xml>` and :meth:`pandas.read_xml`
