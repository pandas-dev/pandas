{{ header }}

.. _api.typing.aliases:

======================================
pandas typing aliases
======================================

**************
Typing aliases
**************

.. module:: pandas.api.typing.aliases

..  raw:: html

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

==================================== ================================================================
Alias                                Meaning
==================================== ================================================================
.. type:: AggFuncType                Type of functions that can be passed to :meth:`DataFrame.agg <pandas.DataFrame.agg>`, :meth:`Series.agg <pandas.Series.agg>`, and :meth:`DataFrameGroupBy.aggregate <pandas.core.groupby.DataFrameGroupBy.aggregate>`
.. type:: AlignJoin                  Argument type for ``join`` in :meth:`DataFrame.align <pandas.DataFrame.align>` and :meth:`Series.align <pandas.Series.align>`
.. type:: AnyAll                     Argument type for ``how`` in :meth:`~pandas.DataFrame.dropna`
.. type:: AnyArrayLike               Used to represent :class:`~pandas.api.extensions.ExtensionArray`, ``numpy`` arrays, :class:`Index` and :class:`Series`
.. type:: ArrayLike                  Used to represent :class:`~pandas.api.extensions.ExtensionArray`, ``numpy`` arrays
.. type:: AstypeArg                  Argument type in :meth:`astype`
.. type:: Axes                       :py:type:`AnyArrayLike` plus sequences (not strings) and ``range``
.. type:: Axis                       Argument type for ``axis`` in many methods
.. type:: CSVEngine                  Argument type for ``engine`` in :meth:`pandas.DataFrame.read_csv`
.. type:: ColspaceArgType            Argument type for ``colspace`` in :meth:`pandas.DataFrame.to_html`
.. type:: CompressionOptions         Argument type for ``compression`` in all I/O output methods except :meth:`~pandas.DataFrame.to_parquet`
.. type:: CorrelationMethod          Argument type for ``correlation`` in :meth:`DataFrame.corr <pandas.DataFrame.corr>` and :meth:`Series.corr <pandas.Series.corr>`
.. type:: DropKeep                   Argument type for ``keep`` in :meth:`DataFrame.drop_duplicates <pandas.DataFrame.drop_duplicates>` and :meth:`Series.drop_duplicates <pandas.Series.drop_duplicates>`
.. type:: Dtype                      Types as objects that can be used to specify dtypes
.. type:: DtypeArg                   Argument type for ``dtype`` in various methods
.. type:: DtypeBackend               Argument type for ``dtype_backend`` in various methods
.. type:: DtypeObj                   Numpy dtypes and Extension dtypes
.. type:: ExcelWriterIfSheetExists   Argument type for ``if_sheet_exists`` in :class:`~pandas.ExcelWriter`
.. type:: ExcelWriterMergeCells      Argument type for ``merge_cells`` in :meth:`DataFrame.to_excel <pandas.DataFrame.to_excel>` and :meth:`Series.to_excel <pandas.Series.to_excel>`
.. type:: FilePath                   Type of paths for files for I/O methods
.. type:: FillnaOptions              Argument type for ``method`` in various methods where NA values are filled
.. type:: FloatFormatType            Argument type for ``float_format`` in :meth:`DataFrame.to_string <pandas.DataFrame.to_string>` and :meth:`Series.to_string <pandas.Series.to_string>`
.. type:: FormattersType             Argument type for ``formatters`` in :meth:`DataFrame.to_string <pandas.DataFrame.to_string>` and :meth:`Series.to_string <pandas.Series.to_string>`
.. type:: FromDictOrient             Argument type for ``orient`` in :meth:`~pandas.DataFrame.from_dict`
.. type:: HTMLFlavors                Argument type for ``flavor`` in :meth:`pandas.read_html`
.. type:: IgnoreRaise                Argument type for ``errors`` in multiple methods
.. type:: IndexLabel                 Argument type for ``level`` in multiple methods
.. type:: InterpolateOptions         Argument type for ``interpolate`` in :meth:`DataFrame.interpolate <pandas.DataFrame.interpolate>` and :meth:`Series.interpolate <pandas.Series.interpolate>`
.. type:: JSONEngine                 Argument type for ``engine`` in :meth:`read_json`
.. type:: JSONSerializable           Argument type for the return type of a callable for argument ``default_handler`` in :meth:`to_json`
.. type:: JoinHow                    Argument type for ``how`` in :meth:`pandas.merge_ordered` and for ``join`` in :meth:`Series.align`
.. type:: JoinValidate               Argument type for ``validate`` in :meth:`~pandas.DataFrame.join`
.. type:: MergeHow                   Argument type for ``how`` in :meth:`merge`
.. type:: MergeValidate              Argument type for ``validate`` in :meth:`merge`
.. type:: NaPosition                 Argument type for ``na_position`` in :meth:`DataFrame.sort_values <pandas.DataFrame.sort_values>` and :meth:`Series.sort_values <pandas.Series.sort_values>`
.. type:: NsmallestNlargestKeep      Argument type for ``keep`` in :meth:`DataFrame.nlargest <pandas.DataFrame.nlargest>`, :meth:`DataFrame.nsmallest <pandas.DataFrame.nsmallest>`, :meth:`Series.nlargest <pandas.Series.nlargest>`, :meth:`Series.nsmallest <pandas.Series.nsmallest>`, and :meth:`SeriesGroupBy.nlargest <pandas.core.groupby.SeriesGroupBy.nlargest>`
.. type:: OpenFileErrors             Argument type for ``errors`` in :meth:`to_hdf` and :meth:`to_csv`
.. type:: Ordered                    Return type for :py:attr:`ordered` in :class:`CategoricalDtype` and :class:`Categorical`
.. type:: ParquetCompressionOptions  Argument type for ``compression`` in :meth:`~pandas.DataFrame.to_parquet`
.. type:: QuantileInterpolation      Argument type for ``interpolation`` in :meth:`DataFrame.quantile <pandas.DataFrame.quantile>` and :meth:`Series.quantile <pandas.Series.quantile>`
.. type:: ReadBuffer                 Additional argument type corresponding to buffers for various file reading methods
.. type:: ReadCsvBuffer              Additional argument type corresponding to buffers for :meth:`pandas.read_csv`
.. type:: ReadPickleBuffer           Additional argument type corresponding to buffers for :meth:`pandas.read_pickle`
.. type:: ReindexMethod              Argument type for ``reindex`` in :meth:`reindex`
.. type:: Scalar                     Types that can be stored in :class:`Series` with non-object dtype
.. type:: SequenceNotStr             Used for arguments that require sequences, but not plain strings
.. type:: SliceType                  Argument types for ``start`` and ``end`` in :meth:`Index.slice_locs`
.. type:: SortKind                   Argument type for ``kind`` in :meth:`DataFrame.sort_values <pandas.DataFrame.sort_values>` and :meth:`Series.sort_values <pandas.Series.sort_values>`
.. type:: StorageOptions             Argument type for ``storage_options`` in various file output methods
.. type:: Suffixes                   Argument type for ``suffixes`` in :meth:`merge`, :meth:`compare` and :meth:`merge_ordered`
.. type:: TakeIndexer                Argument type for ``indexer`` and ``indices`` in :meth:`take`
.. type:: TimeAmbiguous              Argument type for ``ambiguous`` in time operations
.. type:: TimeGrouperOrigin          Argument type for ``origin`` in :meth:`resample` and :class:`TimeGrouper`
.. type:: TimeNonexistent            Argument type for ``nonexistent`` in time operations
.. type:: TimeUnit                   Time unit argument and return type for :py:attr:`unit`, arguments ``unit`` and ``date_unit``
.. type:: TimedeltaConvertibleTypes  Argument type for ``offset`` in various methods, such as :meth:`DataFrame’s <pandas.DataFrame.resample>` and :meth:`Series’ <pandas.Series.resample>` ``resample()``, ``halflife`` in :meth:`DataFrame’s <pandas.DataFrame.ewm>`, :meth:`DataFrameGroupBy’s <.DataFrameGroupBy.ewm>`, and :meth:`Series’ <pandas.Series.ewm>` ``ewm()``, and ``start`` and ``end`` in :meth:`pandas.timedelta_range`
.. type:: TimestampConvertibleTypes  Argument type for ``origin`` in :meth:`DataFrame’s <pandas.DataFrame.resample>` and :meth:`Series’ <pandas.Series.resample>` ``resample()``, and in :meth:`pandas.to_datetime`
.. type:: ToStataByteorder           Argument type for ``byteorder`` in :meth:`~pandas.DataFrame.to_stata`
.. type:: ToTimestampHow             Argument type for ``how`` in :meth:`DataFrame.to_timestamp <pandas.DataFrame.to_timestamp>` and :meth:`Series.to_timestamp <pandas.Series.to_timestamp>`
.. type:: UpdateJoin                 Argument type for ``join`` in :meth:`DataFrame.update <pandas.DataFrame.update>` and :meth:`Series.update <pandas.Series.update>`
.. type:: UsecolsArgType             Argument type for ``usecols`` in :meth:`pandas.read_clipboard`, :meth:`pandas.read_csv` and :meth:`pandas.read_excel`
.. type:: WindowingRankType          Argument type for ``method`` in :meth:`DataFrame’s <pandas.DataFrame.rank>` and :meth:`Series’ <pandas.Series.rank>` ``rank()`` in rolling and expanding window operations
.. type:: WriteBuffer                Additional argument type corresponding to buffers for various file output methods
.. type:: WriteExcelBuffer           Additional argument type corresponding to buffers for :meth:`to_excel`
.. type:: XMLParsers                 Argument type for ``parser`` in :meth:`~pandas.DataFrame.to_xml` and :meth:`pandas.read_xml`

==================================== ================================================================
