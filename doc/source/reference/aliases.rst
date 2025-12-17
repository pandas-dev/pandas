{{ header }}

.. _api.typing.aliases:

======================================
pandas typing aliases
======================================

**************
Typing aliases
**************

.. currentmodule:: pandas.api.typing.aliases

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
:py:type:`AggFuncType`               Type of functions that can be passed to :meth:`agg` methods
:py:type:`AlignJoin`                 Argument type for ``join`` in :meth:`DataFrame.join`
:py:type:`AnyAll`                    Argument type for ``how`` in :meth:`dropna`
:py:type:`AnyArrayLike`              Used to represent :class:`ExtensionArray`, ``numpy`` arrays, :class:`Index` and :class:`Series`
:py:type:`ArrayLike`                 Used to represent :class:`ExtensionArray`, ``numpy`` arrays
:py:type:`AstypeArg`                 Argument type in :meth:`astype`
:py:type:`Axes`                      :py:type:`AnyArrayLike` plus sequences (not strings) and ``range``
:py:type:`Axis`                      Argument type for ``axis`` in many methods
:py:type:`CSVEngine`                 Argument type for ``engine`` in :meth:`DataFrame.read_csv`
:py:type:`ColspaceArgType`           Argument type for ``colspace`` in :meth:`DataFrame.to_html`
:py:type:`CompressionOptions`        Argument type for ``compression`` in all I/O output methods except :meth:`DataFrame.to_parquet`
:py:type:`CorrelationMethod`         Argument type for ``correlation`` in :meth:`corr`
:py:type:`DropKeep`                  Argument type for ``keep`` in :meth:`drop_duplicates`
:py:type:`Dtype`                     Types as objects that can be used to specify dtypes
:py:type:`DtypeArg`                  Argument type for ``dtype`` in various methods
:py:type:`DtypeBackend`              Argument type for ``dtype_backend`` in various methods
:py:type:`DtypeObj`                  Numpy dtypes and Extension dtypes
:py:type:`ExcelWriterIfSheetExists`  Argument type for ``if_sheet_exists`` in :class:`ExcelWriter`
:py:type:`ExcelWriterMergeCells`     Argument type for ``merge_cells`` in :meth:`to_excel`
:py:type:`FilePath`                  Type of paths for files for I/O methods
:py:type:`FillnaOptions`             Argument type for ``method`` in various methods where NA values are filled
:py:type:`FloatFormatType`           Argument type for ``float_format`` in :meth:`to_string`
:py:type:`FormattersType`            Argument type for ``formatters`` in :meth:`to_string`
:py:type:`FromDictOrient`            Argument type for ``orient`` in :meth:`DataFrame.from_dict`
:py:type:`HTMLFlavors`               Argument type for ``flavor`` in :meth:`pandas.read_html`
:py:type:`IgnoreRaise`               Argument type for ``errors`` in multiple methods
:py:type:`IndexLabel`                Argument type for ``level`` in multiple methods
:py:type:`InterpolateOptions`        Argument type for ``interpolate`` in :meth:`interpolate`
:py:type:`JSONEngine`                Argument type for ``engine`` in :meth:`read_json`
:py:type:`JSONSerializable`          Argument type for the return type of a callable for argument ``default_handler`` in :meth:`to_json`
:py:type:`JoinHow`                   Argument type for ``how`` in :meth:`pandas.merge_ordered` and for ``join`` in :meth:`Series.align`
:py:type:`JoinValidate`              Argument type for ``validate`` in :meth:`DataFrame.join`
:py:type:`MergeHow`                  Argument type for ``how`` in :meth:`merge`
:py:type:`MergeValidate`             Argument type for ``validate`` in :meth:`merge`
:py:type:`NaPosition`                Argument type for ``na_position`` in :meth:`sort_index` and :meth:`sort_values`
:py:type:`NsmallestNlargestKeep`     Argument type for ``keep`` in :meth:`nlargest` and :meth:`nsmallest`
:py:type:`OpenFileErrors`            Argument type for ``errors`` in :meth:`to_hdf` and :meth:`to_csv`
:py:type:`Ordered`                   Return type for :py:attr:`ordered` in :class:`CategoricalDtype` and :class:`Categorical`
:py:type:`ParquetCompressionOptions` Argument type for ``compression`` in :meth:`DataFrame.to_parquet`
:py:type:`QuantileInterpolation`     Argument type for ``interpolation`` in :meth:`quantile`
:py:type:`ReadBuffer`                Additional argument type corresponding to buffers for various file reading methods
:py:type:`ReadCsvBuffer`             Additional argument type corresponding to buffers for :meth:`pandas.read_csv`
:py:type:`ReadPickleBuffer`          Additional argument type corresponding to buffers for :meth:`pandas.read_pickle`
:py:type:`ReindexMethod`             Argument type for ``reindex`` in :meth:`reindex`
:py:type:`Scalar`                    Types that can be stored in :class:`Series` with non-object dtype
:py:type:`SequenceNotStr`            Used for arguments that require sequences, but not plain strings
:py:type:`SliceType`                 Argument types for ``start`` and ``end`` in :meth:`Index.slice_locs`
:py:type:`SortKind`                  Argument type for ``kind`` in :meth:`sort_index` and :meth:`sort_values`
:py:type:`StorageOptions`            Argument type for ``storage_options`` in various file output methods
:py:type:`Suffixes`                  Argument type for ``suffixes`` in :meth:`merge`, :meth:`compare` and :meth:`merge_ordered`
:py:type:`TakeIndexer`               Argument type for ``indexer`` and ``indices`` in :meth:`take`
:py:type:`TimeAmbiguous`             Argument type for ``ambiguous`` in time operations
:py:type:`TimeGrouperOrigin`         Argument type for ``origin`` in :meth:`resample` and :class:`TimeGrouper`
:py:type:`TimeNonexistent`           Argument type for ``nonexistent`` in time operations
:py:type:`TimeUnit`                  Time unit argument and return type for :py:attr:`unit`, arguments ``unit`` and ``date_unit``
:py:type:`TimedeltaConvertibleTypes` Argument type for ``offset`` in :meth:`resample`, ``halflife`` in :meth:`ewm` and ``start`` and ``end`` in :meth:`pandas.timedelta_range`
:py:type:`TimestampConvertibleTypes` Argument type for ``origin`` in :meth:`resample` and :meth:`pandas.to_datetime`
:py:type:`ToStataByteorder`          Argument type for ``byteorder`` in :meth:`DataFrame.to_stata`
:py:type:`ToTimestampHow`            Argument type for ``how`` in :meth:`to_timestamp` and ``convention`` in :meth:`resample`
:py:type:`UpdateJoin`                Argument type for ``join`` in :meth:`DataFrame.update`
:py:type:`UsecolsArgType`            Argument type for ``usecols`` in :meth:`pandas.read_clipboard`, :meth:`pandas.read_csv` and :meth:`pandas.read_excel`
:py:type:`WindowingRankType`         Argument type for ``method`` in :meth:`rank` in rolling and expanding window operations
:py:type:`WriteBuffer`               Additional argument type corresponding to buffers for various file output methods
:py:type:`WriteExcelBuffer`          Additional argument type corresponding to buffers for :meth:`to_excel`
:py:type:`XMLParsers`                Argument type for ``parser`` in :meth:`DataFrame.to_xml` and :meth:`pandas.read_xml`
==================================== ================================================================
