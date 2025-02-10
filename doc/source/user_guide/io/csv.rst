.. _io.read_csv_table:

================
CSV & text files
================

The workhorse function for reading text files (a.k.a. flat files) is
:func:`read_csv`. See the :ref:`cookbook<cookbook.csv>` for some advanced strategies.

Parsing options
'''''''''''''''

:func:`read_csv` accepts the following common arguments:

Basic
+++++

filepath_or_buffer : various
  Either a path to a file (a :class:`python:str`, :class:`python:pathlib.Path`)
  URL (including http, ftp, and S3
  locations), or any object with a ``read()`` method (such as an open file or
  :class:`~python:io.StringIO`).
sep : str, defaults to ``','`` for :func:`read_csv`, ``\t`` for :func:`read_table`
  Delimiter to use. If sep is ``None``, the C engine cannot automatically detect
  the separator, but the Python parsing engine can, meaning the latter will be
  used and automatically detect the separator by Python's builtin sniffer tool,
  :class:`python:csv.Sniffer`. In addition, separators longer than 1 character and
  different from ``'\s+'`` will be interpreted as regular expressions and
  will also force the use of the Python parsing engine. Note that regex
  delimiters are prone to ignoring quoted data. Regex example: ``'\\r\\t'``.
delimiter : str, default ``None``
  Alternative argument name for sep.

Column and index locations and names
++++++++++++++++++++++++++++++++++++

header : int or list of ints, default ``'infer'``
  Row number(s) to use as the column names, and the start of the
  data. Default behavior is to infer the column names: if no names are
  passed the behavior is identical to ``header=0`` and column names
  are inferred from the first line of the file, if column names are
  passed explicitly then the behavior is identical to
  ``header=None``. Explicitly pass ``header=0`` to be able to replace
  existing names.

  The header can be a list of ints that specify row locations
  for a MultiIndex on the columns e.g. ``[0,1,3]``. Intervening rows
  that are not specified will be skipped (e.g. 2 in this example is
  skipped). Note that this parameter ignores commented lines and empty
  lines if ``skip_blank_lines=True``, so header=0 denotes the first
  line of data rather than the first line of the file.
names : array-like, default ``None``
  List of column names to use. If file contains no header row, then you should
  explicitly pass ``header=None``. Duplicates in this list are not allowed.
index_col : int, str, sequence of int / str, or False, optional, default ``None``
  Column(s) to use as the row labels of the ``DataFrame``, either given as
  string name or column index. If a sequence of int / str is given, a
  MultiIndex is used.

  .. note::
     ``index_col=False`` can be used to force pandas to *not* use the first
     column as the index, e.g. when you have a malformed file with delimiters at
     the end of each line.

  The default value of ``None`` instructs pandas to guess. If the number of
  fields in the column header row is equal to the number of fields in the body
  of the data file, then a default index is used.  If it is larger, then
  the first columns are used as index so that the remaining number of fields in
  the body are equal to the number of fields in the header.

  The first row after the header is used to determine the number of columns,
  which will go into the index. If the subsequent rows contain less columns
  than the first row, they are filled with ``NaN``.

  This can be avoided through ``usecols``. This ensures that the columns are
  taken as is and the trailing data are ignored.
usecols : list-like or callable, default ``None``
  Return a subset of the columns. If list-like, all elements must either
  be positional (i.e. integer indices into the document columns) or strings
  that correspond to column names provided either by the user in ``names`` or
  inferred from the document header row(s). If ``names`` are given, the document
  header row(s) are not taken into account. For example, a valid list-like
  ``usecols`` parameter would be ``[0, 1, 2]`` or ``['foo', 'bar', 'baz']``.

  Element order is ignored, so ``usecols=[0, 1]`` is the same as ``[1, 0]``. To
  instantiate a DataFrame from ``data`` with element order preserved use
  ``pd.read_csv(data, usecols=['foo', 'bar'])[['foo', 'bar']]`` for columns
  in ``['foo', 'bar']`` order or
  ``pd.read_csv(data, usecols=['foo', 'bar'])[['bar', 'foo']]`` for
  ``['bar', 'foo']`` order.

  If callable, the callable function will be evaluated against the column names,
  returning names where the callable function evaluates to True:

  .. ipython:: python

     import pandas as pd
     from io import StringIO

     data = "col1,col2,col3\na,b,1\na,b,2\nc,d,3"
     pd.read_csv(StringIO(data))
     pd.read_csv(StringIO(data), usecols=lambda x: x.upper() in ["COL1", "COL3"])

  Using this parameter results in much faster parsing time and lower memory usage
  when using the c engine. The Python engine loads the data first before deciding
  which columns to drop.

General parsing configuration
+++++++++++++++++++++++++++++

dtype : Type name or dict of column -> type, default ``None``
  Data type for data or columns. E.g. ``{'a': np.float64, 'b': np.int32, 'c': 'Int64'}``
  Use ``str`` or ``object`` together with suitable ``na_values`` settings to preserve
  and not interpret dtype. If converters are specified, they will be applied INSTEAD
  of dtype conversion.

  .. versionadded:: 1.5.0

     Support for defaultdict was added. Specify a defaultdict as input where
     the default determines the dtype of the columns which are not explicitly
     listed.

dtype_backend : {"numpy_nullable", "pyarrow"}, defaults to NumPy backed DataFrames
  Which dtype_backend to use, e.g. whether a DataFrame should have NumPy
  arrays, nullable dtypes are used for all dtypes that have a nullable
  implementation when "numpy_nullable" is set, pyarrow is used for all
  dtypes if "pyarrow" is set.

  The dtype_backends are still experimental.

  .. versionadded:: 2.0

engine : {``'c'``, ``'python'``, ``'pyarrow'``}
  Parser engine to use. The C and pyarrow engines are faster, while the python engine
  is currently more feature-complete. Multithreading is currently only supported by
  the pyarrow engine.

  .. versionadded:: 1.4.0

     The "pyarrow" engine was added as an *experimental* engine, and some features
     are unsupported, or may not work correctly, with this engine.
converters : dict, default ``None``
  Dict of functions for converting values in certain columns. Keys can either be
  integers or column labels.
true_values : list, default ``None``
  Values to consider as ``True``.
false_values : list, default ``None``
  Values to consider as ``False``.
skipinitialspace : boolean, default ``False``
  Skip spaces after delimiter.
skiprows : list-like or integer, default ``None``
  Line numbers to skip (0-indexed) or number of lines to skip (int) at the start
  of the file.

  If callable, the callable function will be evaluated against the row
  indices, returning True if the row should be skipped and False otherwise:

  .. ipython:: python

     from io import StringIO

     data = "col1,col2,col3\na,b,1\na,b,2\nc,d,3"
     pd.read_csv(StringIO(data))
     pd.read_csv(StringIO(data), skiprows=lambda x: x % 2 != 0)

skipfooter : int, default ``0``
  Number of lines at bottom of file to skip (unsupported with engine='c').

nrows : int, default ``None``
  Number of rows of file to read. Useful for reading pieces of large files.
low_memory : boolean, default ``True``
  Internally process the file in chunks, resulting in lower memory use
  while parsing, but possibly mixed type inference.  To ensure no mixed
  types either set ``False``, or specify the type with the ``dtype`` parameter.
  Note that the entire file is read into a single ``DataFrame`` regardless,
  use the ``chunksize`` or ``iterator`` parameter to return the data in chunks.
  (Only valid with C parser)
memory_map : boolean, default False
  If a filepath is provided for ``filepath_or_buffer``, map the file object
  directly onto memory and access the data directly from there. Using this
  option can improve performance because there is no longer any I/O overhead.

NA and missing data handling
++++++++++++++++++++++++++++

na_values : scalar, str, list-like, or dict, default ``None``
  Additional strings to recognize as NA/NaN. If dict passed, specific per-column
  NA values. See :ref:`na values const <io.navaluesconst>` below
  for a list of the values interpreted as NaN by default.

keep_default_na : boolean, default ``True``
  Whether or not to include the default NaN values when parsing the data.
  Depending on whether ``na_values`` is passed in, the behavior is as follows:

  * If ``keep_default_na`` is ``True``, and ``na_values`` are specified, ``na_values``
    is appended to the default NaN values used for parsing.
  * If ``keep_default_na`` is ``True``, and ``na_values`` are not specified, only
    the default NaN values are used for parsing.
  * If ``keep_default_na`` is ``False``, and ``na_values`` are specified, only
    the NaN values specified ``na_values`` are used for parsing.
  * If ``keep_default_na`` is ``False``, and ``na_values`` are not specified, no
    strings will be parsed as NaN.

  Note that if ``na_filter`` is passed in as ``False``, the ``keep_default_na`` and
  ``na_values`` parameters will be ignored.
na_filter : boolean, default ``True``
  Detect missing value markers (empty strings and the value of na_values). In
  data without any NAs, passing ``na_filter=False`` can improve the performance
  of reading a large file.
verbose : boolean, default ``False``
  Indicate number of NA values placed in non-numeric columns.
skip_blank_lines : boolean, default ``True``
  If ``True``, skip over blank lines rather than interpreting as NaN values.

.. _io.read_csv_table.datetime:

Datetime handling
+++++++++++++++++

parse_dates : boolean or list of ints or names or list of lists or dict, default ``False``.
  * If ``True`` -> try parsing the index.
  * If ``[1, 2, 3]`` ->  try parsing columns 1, 2, 3 each as a separate date
    column.

  .. note::
     A fast-path exists for iso8601-formatted dates.
date_format : str or dict of column -> format, default ``None``
   If used in conjunction with ``parse_dates``, will parse dates according to this
   format. For anything more complex,
   please read in as ``object`` and then apply :func:`to_datetime` as-needed.

   .. versionadded:: 2.0.0
dayfirst : boolean, default ``False``
  DD/MM format dates, international and European format.
cache_dates : boolean, default True
  If True, use a cache of unique, converted dates to apply the datetime
  conversion. May produce significant speed-up when parsing duplicate
  date strings, especially ones with timezone offsets.

Iteration
+++++++++

iterator : boolean, default ``False``
  Return ``TextFileReader`` object for iteration or getting chunks with
  ``get_chunk()``.
chunksize : int, default ``None``
  Return ``TextFileReader`` object for iteration. See :ref:`iterating and chunking
  <io.chunking>` below.

Quoting, compression, and file format
+++++++++++++++++++++++++++++++++++++

compression : {``'infer'``, ``'gzip'``, ``'bz2'``, ``'zip'``, ``'xz'``, ``'zstd'``, ``None``, ``dict``}, default ``'infer'``
  For on-the-fly decompression of on-disk data. If 'infer', then use gzip,
  bz2, zip, xz, or zstandard if ``filepath_or_buffer`` is path-like ending in '.gz', '.bz2',
  '.zip', '.xz', '.zst', respectively, and no decompression otherwise. If using 'zip',
  the ZIP file must contain only one data file to be read in.
  Set to ``None`` for no decompression. Can also be a dict with key ``'method'``
  set to one of {``'zip'``, ``'gzip'``, ``'bz2'``, ``'zstd'``} and other key-value pairs are
  forwarded to ``zipfile.ZipFile``, ``gzip.GzipFile``, ``bz2.BZ2File``, or ``zstandard.ZstdDecompressor``.
  As an example, the following could be passed for faster compression and to
  create a reproducible gzip archive:
  ``compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}``.

  .. versionchanged:: 1.2.0 Previous versions forwarded dict entries for 'gzip' to ``gzip.open``.
thousands : str, default ``None``
  Thousands separator.
decimal : str, default ``'.'``
  Character to recognize as decimal point. E.g. use ``','`` for European data.
float_precision : string, default None
  Specifies which converter the C engine should use for floating-point values.
  The options are ``None`` for the ordinary converter, ``high`` for the
  high-precision converter, and ``round_trip`` for the round-trip converter.
lineterminator : str (length 1), default ``None``
  Character to break file into lines. Only valid with C parser.
quotechar : str (length 1)
  The character used to denote the start and end of a quoted item. Quoted items
  can include the delimiter and it will be ignored.
quoting : int or ``csv.QUOTE_*`` instance, default ``0``
  Control field quoting behavior per ``csv.QUOTE_*`` constants. Use one of
  ``QUOTE_MINIMAL`` (0), ``QUOTE_ALL`` (1), ``QUOTE_NONNUMERIC`` (2) or
  ``QUOTE_NONE`` (3).
doublequote : boolean, default ``True``
   When ``quotechar`` is specified and ``quoting`` is not ``QUOTE_NONE``,
   indicate whether or not to interpret two consecutive ``quotechar`` elements
   **inside** a field as a single ``quotechar`` element.
escapechar : str (length 1), default ``None``
  One-character string used to escape delimiter when quoting is ``QUOTE_NONE``.
comment : str, default ``None``
  Indicates remainder of line should not be parsed. If found at the beginning of
  a line, the line will be ignored altogether. This parameter must be a single
  character. Like empty lines (as long as ``skip_blank_lines=True``), fully
  commented lines are ignored by the parameter ``header`` but not by ``skiprows``.
  For example, if ``comment='#'``, parsing '#empty\\na,b,c\\n1,2,3' with
  ``header=0`` will result in 'a,b,c' being treated as the header.
encoding : str, default ``None``
  Encoding to use for UTF when reading/writing (e.g. ``'utf-8'``). `List of
  Python standard encodings
  <https://docs.python.org/3/library/codecs.html#standard-encodings>`_.
dialect : str or :class:`python:csv.Dialect` instance, default ``None``
  If provided, this parameter will override values (default or not) for the
  following parameters: ``delimiter``, ``doublequote``, ``escapechar``,
  ``skipinitialspace``, ``quotechar``, and ``quoting``. If it is necessary to
  override values, a ParserWarning will be issued. See :class:`python:csv.Dialect`
  documentation for more details.

Error handling
++++++++++++++

on_bad_lines : {{'error', 'warn', 'skip'}}, default 'error'
    Specifies what to do upon encountering a bad line (a line with too many fields).
    Allowed values are :

    - 'error', raise an ParserError when a bad line is encountered.
    - 'warn', print a warning when a bad line is encountered and skip that line.
    - 'skip', skip bad lines without raising or warning when they are encountered.

    .. versionadded:: 1.3.0

.. _io.dtypes:

Specifying column data types
''''''''''''''''''''''''''''

You can indicate the data type for the whole ``DataFrame`` or individual
columns:

.. ipython:: python

    import numpy as np
    from io import StringIO

    data = "a,b,c,d\n1,2,3,4\n5,6,7,8\n9,10,11"
    print(data)

    df = pd.read_csv(StringIO(data), dtype=object)
    df
    df["a"][0]
    df = pd.read_csv(StringIO(data), dtype={"b": object, "c": np.float64, "d": "Int64"})
    df.dtypes

Fortunately, pandas offers more than one way to ensure that your column(s)
contain only one ``dtype``. If you're unfamiliar with these concepts, you can
see :ref:`here<basics.dtypes>` to learn more about dtypes, and
:ref:`here<basics.object_conversion>` to learn more about ``object`` conversion in
pandas.


For instance, you can use the ``converters`` argument
of :func:`~pandas.read_csv`:

.. ipython:: python

    from io import StringIO

    data = "col_1\n1\n2\n'A'\n4.22"
    df = pd.read_csv(StringIO(data), converters={"col_1": str})
    df
    df["col_1"].apply(type).value_counts()

Or you can use the :func:`~pandas.to_numeric` function to coerce the
dtypes after reading in the data,

.. ipython:: python

    from io import StringIO

    df2 = pd.read_csv(StringIO(data))
    df2["col_1"] = pd.to_numeric(df2["col_1"], errors="coerce")
    df2
    df2["col_1"].apply(type).value_counts()

which will convert all valid parsing to floats, leaving the invalid parsing
as ``NaN``.

Ultimately, how you deal with reading in columns containing mixed dtypes
depends on your specific needs. In the case above, if you wanted to ``NaN`` out
the data anomalies, then :func:`~pandas.to_numeric` is probably your best option.
However, if you wanted for all the data to be coerced, no matter the type, then
using the ``converters`` argument of :func:`~pandas.read_csv` would certainly be
worth trying.

.. note::
   In some cases, reading in abnormal data with columns containing mixed dtypes
   will result in an inconsistent dataset. If you rely on pandas to infer the
   dtypes of your columns, the parsing engine will go and infer the dtypes for
   different chunks of the data, rather than the whole dataset at once. Consequently,
   you can end up with column(s) with mixed dtypes. For example,

   .. ipython:: python
        :okwarning:

        col_1 = list(range(500000)) + ["a", "b"] + list(range(500000))
        df = pd.DataFrame({"col_1": col_1})
        df.to_csv("foo.csv")
        mixed_df = pd.read_csv("foo.csv")
        mixed_df["col_1"].apply(type).value_counts()
        mixed_df["col_1"].dtype

   will result with ``mixed_df`` containing an ``int`` dtype for certain chunks
   of the column, and ``str`` for others due to the mixed dtypes from the
   data that was read in. It is important to note that the overall column will be
   marked with a ``dtype`` of ``object``, which is used for columns with mixed dtypes.

.. ipython:: python
   :suppress:

   import os

   os.remove("foo.csv")

Setting ``dtype_backend="numpy_nullable"`` will result in nullable dtypes for every column.

.. ipython:: python

   from io import StringIO

   data = """a,b,c,d,e,f,g,h,i,j
   1,2.5,True,a,,,,,12-31-2019,
   3,4.5,False,b,6,7.5,True,a,12-31-2019,
   """

   df = pd.read_csv(StringIO(data), dtype_backend="numpy_nullable", parse_dates=["i"])
   df
   df.dtypes

.. _io.categorical:

Specifying categorical dtype
''''''''''''''''''''''''''''

``Categorical`` columns can be parsed directly by specifying ``dtype='category'`` or
``dtype=CategoricalDtype(categories, ordered)``.

.. ipython:: python

   from io import StringIO

   data = "col1,col2,col3\na,b,1\na,b,2\nc,d,3"

   pd.read_csv(StringIO(data))
   pd.read_csv(StringIO(data)).dtypes
   pd.read_csv(StringIO(data), dtype="category").dtypes

Individual columns can be parsed as a ``Categorical`` using a dict
specification:

.. ipython:: python

   from io import StringIO

   pd.read_csv(StringIO(data), dtype={"col1": "category"}).dtypes

Specifying ``dtype='category'`` will result in an unordered ``Categorical``
whose ``categories`` are the unique values observed in the data. For more
control on the categories and order, create a
:class:`~pandas.api.types.CategoricalDtype` ahead of time, and pass that for
that column's ``dtype``.

.. ipython:: python

   from pandas.api.types import CategoricalDtype
   from io import StringIO

   dtype = CategoricalDtype(["d", "c", "b", "a"], ordered=True)
   pd.read_csv(StringIO(data), dtype={"col1": dtype}).dtypes

When using ``dtype=CategoricalDtype``, "unexpected" values outside of
``dtype.categories`` are treated as missing values.

.. ipython:: python

   from io import StringIO

   dtype = CategoricalDtype(["a", "b", "d"])  # No 'c'
   pd.read_csv(StringIO(data), dtype={"col1": dtype}).col1

This matches the behavior of :meth:`Categorical.set_categories`.

.. note::

   With ``dtype='category'``, the resulting categories will always be parsed
   as strings (object dtype). If the categories are numeric they can be
   converted using the :func:`to_numeric` function, or as appropriate, another
   converter such as :func:`to_datetime`.

   When ``dtype`` is a ``CategoricalDtype`` with homogeneous ``categories`` (
   all numeric, all datetimes, etc.), the conversion is done automatically.

   .. ipython:: python

      from io import StringIO

      df = pd.read_csv(StringIO(data), dtype="category")
      df.dtypes
      df["col3"]
      new_categories = pd.to_numeric(df["col3"].cat.categories)
      df["col3"] = df["col3"].cat.rename_categories(new_categories)
      df["col3"]


Naming and using columns
''''''''''''''''''''''''

.. _io.headers:

Handling column names
+++++++++++++++++++++

A file may or may not have a header row. pandas assumes the first row should be
used as the column names:

.. ipython:: python

    from io import StringIO

    data = "a,b,c\n1,2,3\n4,5,6\n7,8,9"
    print(data)
    pd.read_csv(StringIO(data))

By specifying the ``names`` argument in conjunction with ``header`` you can
indicate other names to use and whether or not to throw away the header row (if
any):

.. ipython:: python

    from io import StringIO

    print(data)
    pd.read_csv(StringIO(data), names=["foo", "bar", "baz"], header=0)
    pd.read_csv(StringIO(data), names=["foo", "bar", "baz"], header=None)

If the header is in a row other than the first, pass the row number to
``header``. This will skip the preceding rows:

.. ipython:: python

    from io import StringIO

    data = "skip this skip it\na,b,c\n1,2,3\n4,5,6\n7,8,9"
    pd.read_csv(StringIO(data), header=1)

.. note::

  Default behavior is to infer the column names: if no names are
  passed the behavior is identical to ``header=0`` and column names
  are inferred from the first non-blank line of the file, if column
  names are passed explicitly then the behavior is identical to
  ``header=None``.

.. _io.dupe_names:

Duplicate names parsing
'''''''''''''''''''''''

If the file or header contains duplicate names, pandas will by default
distinguish between them so as to prevent overwriting data:

.. ipython:: python

   from io import StringIO

   data = "a,b,a\n0,1,2\n3,4,5"
   pd.read_csv(StringIO(data))

There is no more duplicate data because duplicate columns 'X', ..., 'X' become
'X', 'X.1', ..., 'X.N'.

.. _io.usecols:

Filtering columns (``usecols``)
+++++++++++++++++++++++++++++++

The ``usecols`` argument allows you to select any subset of the columns in a
file, either using the column names, position numbers or a callable:

.. ipython:: python

    from io import StringIO

    data = "a,b,c,d\n1,2,3,foo\n4,5,6,bar\n7,8,9,baz"
    pd.read_csv(StringIO(data))
    pd.read_csv(StringIO(data), usecols=["b", "d"])
    pd.read_csv(StringIO(data), usecols=[0, 2, 3])
    pd.read_csv(StringIO(data), usecols=lambda x: x.upper() in ["A", "C"])

The ``usecols`` argument can also be used to specify which columns not to
use in the final result:

.. ipython:: python

   from io import StringIO

   pd.read_csv(StringIO(data), usecols=lambda x: x not in ["a", "c"])

In this case, the callable is specifying that we exclude the "a" and "c"
columns from the output.

Comments and empty lines
''''''''''''''''''''''''

.. _io.skiplines:

Ignoring line comments and empty lines
++++++++++++++++++++++++++++++++++++++

If the ``comment`` parameter is specified, then completely commented lines will
be ignored. By default, completely blank lines will be ignored as well.

.. ipython:: python

   from io import StringIO

   data = "\na,b,c\n  \n# commented line\n1,2,3\n\n4,5,6"
   print(data)
   pd.read_csv(StringIO(data), comment="#")

If ``skip_blank_lines=False``, then ``read_csv`` will not ignore blank lines:

.. ipython:: python

   from io import StringIO

   data = "a,b,c\n\n1,2,3\n\n\n4,5,6"
   pd.read_csv(StringIO(data), skip_blank_lines=False)

.. warning::

   The presence of ignored lines might create ambiguities involving line numbers;
   the parameter ``header`` uses row numbers (ignoring commented/empty
   lines), while ``skiprows`` uses line numbers (including commented/empty lines):

   .. ipython:: python

      from io import StringIO

      data = "#comment\na,b,c\nA,B,C\n1,2,3"
      pd.read_csv(StringIO(data), comment="#", header=1)
      data = "A,B,C\n#comment\na,b,c\n1,2,3"
      pd.read_csv(StringIO(data), comment="#", skiprows=2)

   If both ``header`` and ``skiprows`` are specified, ``header`` will be
   relative to the end of ``skiprows``. For example:

.. ipython:: python

   from io import StringIO

   data = (
       "# empty\n"
       "# second empty line\n"
       "# third emptyline\n"
       "X,Y,Z\n"
       "1,2,3\n"
       "A,B,C\n"
       "1,2.,4.\n"
       "5.,NaN,10.0\n"
   )
   print(data)
   pd.read_csv(StringIO(data), comment="#", skiprows=4, header=1)

.. _io.comments:

Comments
++++++++

Sometimes comments or meta data may be included in a file:

.. ipython:: python

   data = (
       "ID,level,category\n"
       "Patient1,123000,x # really unpleasant\n"
       "Patient2,23000,y # wouldn't take his medicine\n"
       "Patient3,1234018,z # awesome"
   )
   with open("tmp.csv", "w") as fh:
       fh.write(data)

   print(open("tmp.csv").read())

By default, the parser includes the comments in the output:

.. ipython:: python

   df = pd.read_csv("tmp.csv")
   df

We can suppress the comments using the ``comment`` keyword:

.. ipython:: python

   df = pd.read_csv("tmp.csv", comment="#")
   df

.. ipython:: python
   :suppress:

   os.remove("tmp.csv")

.. _io.unicode:

Dealing with Unicode data
'''''''''''''''''''''''''

The ``encoding`` argument should be used for encoded unicode data, which will
result in byte strings being decoded to unicode in the result:

.. ipython:: python

   from io import BytesIO

   data = b"word,length\n" b"Tr\xc3\xa4umen,7\n" b"Gr\xc3\xbc\xc3\x9fe,5"
   data = data.decode("utf8").encode("latin-1")
   df = pd.read_csv(BytesIO(data), encoding="latin-1")
   df
   df["word"][1]

Some formats which encode all characters as multiple bytes, like UTF-16, won't
parse correctly at all without specifying the encoding. `Full list of Python
standard encodings
<https://docs.python.org/3/library/codecs.html#standard-encodings>`_.

.. _io.index_col:

Index columns and trailing delimiters
'''''''''''''''''''''''''''''''''''''

If a file has one more column of data than the number of column names, the
first column will be used as the ``DataFrame``'s row names:

.. ipython:: python

    from io import StringIO

    data = "a,b,c\n4,apple,bat,5.7\n8,orange,cow,10"
    pd.read_csv(StringIO(data))

.. ipython:: python

    from io import StringIO

    data = "index,a,b,c\n4,apple,bat,5.7\n8,orange,cow,10"
    pd.read_csv(StringIO(data), index_col=0)

Ordinarily, you can achieve this behavior using the ``index_col`` option.

There are some exception cases when a file has been prepared with delimiters at
the end of each data line, confusing the parser. To explicitly disable the
index column inference and discard the last column, pass ``index_col=False``:

.. ipython:: python

    from io import StringIO

    data = "a,b,c\n4,apple,bat,\n8,orange,cow,"
    print(data)
    pd.read_csv(StringIO(data))
    pd.read_csv(StringIO(data), index_col=False)

If a subset of data is being parsed using the ``usecols`` option, the
``index_col`` specification is based on that subset, not the original data.

.. ipython:: python

    from io import StringIO

    data = "a,b,c\n4,apple,bat,\n8,orange,cow,"
    print(data)
    pd.read_csv(StringIO(data), usecols=["b", "c"])
    pd.read_csv(StringIO(data), usecols=["b", "c"], index_col=0)

.. _io.parse_dates:

Date Handling
'''''''''''''

Specifying date columns
+++++++++++++++++++++++

To better facilitate working with datetime data, :func:`read_csv`
uses the keyword arguments ``parse_dates`` and ``date_format``
to allow users to specify a variety of columns and date/time formats to turn the
input text data into ``datetime`` objects.

The simplest case is to just pass in ``parse_dates=True``:

.. ipython:: python

   with open("foo.csv", mode="w") as f:
       f.write("date,A,B,C\n20090101,a,1,2\n20090102,b,3,4\n20090103,c,4,5")

   # Use a column as an index, and parse it as dates.
   df = pd.read_csv("foo.csv", index_col=0, parse_dates=True)
   df

   # These are Python datetime objects
   df.index

It is often the case that we may want to store date and time data separately,
or store various date fields separately. the ``parse_dates`` keyword can be
used to specify columns to parse the dates and/or times.


.. note::
   If a column or index contains an unparsable date, the entire column or
   index will be returned unaltered as an object data type. For non-standard
   datetime parsing, use :func:`to_datetime` after ``pd.read_csv``.


.. note::
   read_csv has a fast_path for parsing datetime strings in iso8601 format,
   e.g "2000-01-01T00:01:02+00:00" and similar variations. If you can arrange
   for your data to store datetimes in this format, load times will be
   significantly faster, ~20x has been observed.


Date parsing functions
++++++++++++++++++++++

Finally, the parser allows you to specify a custom ``date_format``.
Performance-wise, you should try these methods of parsing dates in order:

1. If you know the format, use ``date_format``, e.g.:
   ``date_format="%d/%m/%Y"`` or ``date_format={column_name: "%d/%m/%Y"}``.

2. If you different formats for different columns, or want to pass any extra options (such
   as ``utc``) to ``to_datetime``, then you should read in your data as ``object`` dtype, and
   then use ``to_datetime``.


.. _io.csv.mixed_timezones:

Parsing a CSV with mixed timezones
++++++++++++++++++++++++++++++++++

pandas cannot natively represent a column or index with mixed timezones. If your CSV
file contains columns with a mixture of timezones, the default result will be
an object-dtype column with strings, even with ``parse_dates``.
To parse the mixed-timezone values as a datetime column, read in as ``object`` dtype and
then call :func:`to_datetime` with ``utc=True``.


.. ipython:: python

   from io import StringIO

   content = """\
   a
   2000-01-01T00:00:00+05:00
   2000-01-01T00:00:00+06:00"""
   df = pd.read_csv(StringIO(content))
   df["a"] = pd.to_datetime(df["a"], utc=True)
   df["a"]


.. _io.dayfirst:


Inferring datetime format
+++++++++++++++++++++++++

Here are some examples of datetime strings that can be guessed (all
representing December 30th, 2011 at 00:00:00):

* "20111230"
* "2011/12/30"
* "20111230 00:00:00"
* "12/30/2011 00:00:00"
* "30/Dec/2011 00:00:00"
* "30/December/2011 00:00:00"

Note that format inference is sensitive to ``dayfirst``.  With
``dayfirst=True``, it will guess "01/12/2011" to be December 1st. With
``dayfirst=False`` (default) it will guess "01/12/2011" to be January 12th.

If you try to parse a column of date strings, pandas will attempt to guess the format
from the first non-NaN element, and will then parse the rest of the column with that
format. If pandas fails to guess the format (for example if your first string is
``'01 December US/Pacific 2000'``), then a warning will be raised and each
row will be parsed individually by ``dateutil.parser.parse``. The safest
way to parse dates is to explicitly set ``format=``.

.. ipython:: python

   df = pd.read_csv(
       "foo.csv",
       index_col=0,
       parse_dates=True,
   )
   df

In the case that you have mixed datetime formats within the same column, you can
pass  ``format='mixed'``

.. ipython:: python

   from io import StringIO

   data = StringIO("date\n12 Jan 2000\n2000-01-13\n")
   df = pd.read_csv(data)
   df['date'] = pd.to_datetime(df['date'], format='mixed')
   df

or, if your datetime formats are all ISO8601 (possibly not identically-formatted):

.. ipython:: python

   from io import StringIO

   data = StringIO("date\n2020-01-01\n2020-01-01 03:00\n")
   df = pd.read_csv(data)
   df['date'] = pd.to_datetime(df['date'], format='ISO8601')
   df

.. ipython:: python
   :suppress:

   os.remove("foo.csv")

International date formats
++++++++++++++++++++++++++

While US date formats tend to be MM/DD/YYYY, many international formats use
DD/MM/YYYY instead. For convenience, a ``dayfirst`` keyword is provided:

.. ipython:: python

   data = "date,value,cat\n1/6/2000,5,a\n2/6/2000,10,b\n3/6/2000,15,c"
   print(data)
   with open("tmp.csv", "w") as fh:
       fh.write(data)

   pd.read_csv("tmp.csv", parse_dates=[0])
   pd.read_csv("tmp.csv", dayfirst=True, parse_dates=[0])

.. ipython:: python
   :suppress:

   os.remove("tmp.csv")

Writing CSVs to binary file objects
+++++++++++++++++++++++++++++++++++

.. versionadded:: 1.2.0

``df.to_csv(..., mode="wb")`` allows writing a CSV to a file object
opened binary mode. In most cases, it is not necessary to specify
``mode`` as pandas will auto-detect whether the file object is
opened in text or binary mode.

.. ipython:: python

   import io

   data = pd.DataFrame([0, 1, 2])
   buffer = io.BytesIO()
   data.to_csv(buffer, encoding="utf-8", compression="gzip")

.. _io.float_precision:

Specifying method for floating-point conversion
'''''''''''''''''''''''''''''''''''''''''''''''

The parameter ``float_precision`` can be specified in order to use
a specific floating-point converter during parsing with the C engine.
The options are the ordinary converter, the high-precision converter, and
the round-trip converter (which is guaranteed to round-trip values after
writing to a file). For example:

.. ipython:: python

   from io import StringIO

   val = "0.3066101993807095471566981359501369297504425048828125"
   data = "a,b,c\n1,2,{0}".format(val)
   abs(
       pd.read_csv(
           StringIO(data),
           engine="c",
           float_precision=None,
       )["c"][0] - float(val)
   )
   abs(
       pd.read_csv(
           StringIO(data),
           engine="c",
           float_precision="high",
       )["c"][0] - float(val)
   )
   abs(
       pd.read_csv(StringIO(data), engine="c", float_precision="round_trip")["c"][0]
       - float(val)
   )


.. _io.thousands:

Thousand separators
'''''''''''''''''''

For large numbers that have been written with a thousands separator, you can
set the ``thousands`` keyword to a string of length 1 so that integers will be parsed
correctly.

By default, numbers with a thousands separator will be parsed as strings:

.. ipython:: python

   data = (
       "ID|level|category\n"
       "Patient1|123,000|x\n"
       "Patient2|23,000|y\n"
       "Patient3|1,234,018|z"
   )

   with open("tmp.csv", "w") as fh:
       fh.write(data)

   df = pd.read_csv("tmp.csv", sep="|")
   df

   df.level.dtype

The ``thousands`` keyword allows integers to be parsed correctly:

.. ipython:: python

    df = pd.read_csv("tmp.csv", sep="|", thousands=",")
    df

    df.level.dtype

.. ipython:: python
   :suppress:

   os.remove("tmp.csv")

.. _io.na_values:

NA values
'''''''''

To control which values are parsed as missing values (which are signified by
``NaN``), specify a string in ``na_values``. If you specify a list of strings,
then all values in it are considered to be missing values. If you specify a
number (a ``float``, like ``5.0`` or an ``integer`` like ``5``), the
corresponding equivalent values will also imply a missing value (in this case
effectively ``[5.0, 5]`` are recognized as ``NaN``).

To completely override the default values that are recognized as missing, specify ``keep_default_na=False``.

.. _io.navaluesconst:

The default ``NaN`` recognized values are ``['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A',
'n/a', 'NA', '<NA>', '#NA', 'NULL', 'null', 'NaN', '-NaN', 'nan', '-nan', 'None', '']``.

Let us consider some examples:

.. code-block:: python

   pd.read_csv("path_to_file.csv", na_values=[5])

In the example above ``5`` and ``5.0`` will be recognized as ``NaN``, in
addition to the defaults. A string will first be interpreted as a numerical
``5``, then as a ``NaN``.

.. code-block:: python

   pd.read_csv("path_to_file.csv", keep_default_na=False, na_values=[""])

Above, only an empty field will be recognized as ``NaN``.

.. code-block:: python

   pd.read_csv("path_to_file.csv", keep_default_na=False, na_values=["NA", "0"])

Above, both ``NA`` and ``0`` as strings are ``NaN``.

.. code-block:: python

   pd.read_csv("path_to_file.csv", na_values=["Nope"])

The default values, in addition to the string ``"Nope"`` are recognized as
``NaN``.

.. _io.infinity:

Infinity
''''''''

``inf`` like values will be parsed as ``np.inf`` (positive infinity), and ``-inf`` as ``-np.inf`` (negative infinity).
These will ignore the case of the value, meaning ``Inf``, will also be parsed as ``np.inf``.

.. _io.boolean:

Boolean values
''''''''''''''

The common values ``True``, ``False``, ``TRUE``, and ``FALSE`` are all
recognized as boolean. Occasionally you might want to recognize other values
as being boolean. To do this, use the ``true_values`` and ``false_values``
options as follows:

.. ipython:: python

    from io import StringIO

    data = "a,b,c\n1,Yes,2\n3,No,4"
    print(data)
    pd.read_csv(StringIO(data))
    pd.read_csv(StringIO(data), true_values=["Yes"], false_values=["No"])

.. _io.bad_lines:

Handling "bad" lines
''''''''''''''''''''

Some files may have malformed lines with too few fields or too many. Lines with
too few fields will have NA values filled in the trailing fields. Lines with
too many fields will raise an error by default:

.. ipython:: python
    :okexcept:

    from io import StringIO

    data = "a,b,c\n1,2,3\n4,5,6,7\n8,9,10"
    pd.read_csv(StringIO(data))

You can elect to skip bad lines:

.. ipython:: python

    from io import StringIO

    data = "a,b,c\n1,2,3\n4,5,6,7\n8,9,10"
    pd.read_csv(StringIO(data), on_bad_lines="skip")

.. versionadded:: 1.4.0

Or pass a callable function to handle the bad line if ``engine="python"``.
The bad line will be a list of strings that was split by the ``sep``:

.. ipython:: python

    from io import StringIO

    external_list = []
    def bad_lines_func(line):
        external_list.append(line)
        return line[-3:]
    pd.read_csv(StringIO(data), on_bad_lines=bad_lines_func, engine="python")
    external_list

.. note::

   The callable function will handle only a line with too many fields.
   Bad lines caused by other errors will be silently skipped.

   .. ipython:: python

      from io import StringIO

      bad_lines_func = lambda line: print(line)

      data = 'name,type\nname a,a is of type a\nname b,"b\" is of type b"'
      data
      pd.read_csv(StringIO(data), on_bad_lines=bad_lines_func, engine="python")

   The line was not processed in this case, as a "bad line" here is caused by an escape character.

You can also use the ``usecols`` parameter to eliminate extraneous column
data that appear in some lines but not others:

.. ipython:: python
   :okexcept:

   from io import StringIO

   pd.read_csv(StringIO(data), usecols=[0, 1, 2])

In case you want to keep all data including the lines with too many fields, you can
specify a sufficient number of ``names``. This ensures that lines with not enough
fields are filled with ``NaN``.

.. ipython:: python

   from io import StringIO

   pd.read_csv(StringIO(data), names=['a', 'b', 'c', 'd'])

.. _io.dialect:

Dialect
'''''''

The ``dialect`` keyword gives greater flexibility in specifying the file format.
By default it uses the Excel dialect but you can specify either the dialect name
or a :class:`python:csv.Dialect` instance.

Suppose you had data with unenclosed quotes:

.. ipython:: python

   data = "label1,label2,label3\n" 'index1,"a,c,e\n' "index2,b,d,f"
   print(data)

By default, ``read_csv`` uses the Excel dialect and treats the double quote as
the quote character, which causes it to fail when it finds a newline before it
finds the closing double quote.

We can get around this using ``dialect``:

.. ipython:: python
   :okwarning:

   import csv
   from io import StringIO

   dia = csv.excel()
   dia.quoting = csv.QUOTE_NONE
   pd.read_csv(StringIO(data), dialect=dia)

All of the dialect options can be specified separately by keyword arguments:

.. ipython:: python

    from io import StringIO

    data = "a,b,c~1,2,3~4,5,6"
    pd.read_csv(StringIO(data), lineterminator="~")

Another common dialect option is ``skipinitialspace``, to skip any whitespace
after a delimiter:

.. ipython:: python

   from io import StringIO

   data = "a, b, c\n1, 2, 3\n4, 5, 6"
   print(data)
   pd.read_csv(StringIO(data), skipinitialspace=True)

The parsers make every attempt to "do the right thing" and not be fragile. Type
inference is a pretty big deal. If a column can be coerced to integer dtype
without altering the contents, the parser will do so. Any non-numeric
columns will come through as object dtype as with the rest of pandas objects.

.. _io.quoting:

Quoting and Escape Characters
'''''''''''''''''''''''''''''

Quotes (and other escape characters) in embedded fields can be handled in any
number of ways. One way is to use backslashes; to properly parse this data, you
should pass the ``escapechar`` option:

.. ipython:: python

   from io import StringIO

   data = 'a,b\n"hello, \\"Bob\\", nice to see you",5'
   print(data)
   pd.read_csv(StringIO(data), escapechar="\\")

.. _io.fwf_reader:
.. _io.fwf:

Files with fixed width columns
''''''''''''''''''''''''''''''

While :func:`read_csv` reads delimited data, the :func:`read_fwf` function works
with data files that have known and fixed column widths. The function parameters
to ``read_fwf`` are largely the same as ``read_csv`` with two extra parameters, and
a different usage of the ``delimiter`` parameter:

* ``colspecs``: A list of pairs (tuples) giving the extents of the
  fixed-width fields of each line as half-open intervals (i.e.,  [from, to[ ).
  String value 'infer' can be used to instruct the parser to try detecting
  the column specifications from the first 100 rows of the data. Default
  behavior, if not specified, is to infer.
* ``widths``: A list of field widths which can be used instead of 'colspecs'
  if the intervals are contiguous.
* ``delimiter``: Characters to consider as filler characters in the fixed-width file.
  Can be used to specify the filler character of the fields
  if it is not spaces (e.g., '~').

Consider a typical fixed-width data file:

.. ipython:: python

   data1 = (
       "id8141    360.242940   149.910199   11950.7\n"
       "id1594    444.953632   166.985655   11788.4\n"
       "id1849    364.136849   183.628767   11806.2\n"
       "id1230    413.836124   184.375703   11916.8\n"
       "id1948    502.953953   173.237159   12468.3"
   )
   with open("bar.csv", "w") as f:
       f.write(data1)

In order to parse this file into a ``DataFrame``, we simply need to supply the
column specifications to the ``read_fwf`` function along with the file name:

.. ipython:: python

   # Column specifications are a list of half-intervals
   colspecs = [(0, 6), (8, 20), (21, 33), (34, 43)]
   df = pd.read_fwf("bar.csv", colspecs=colspecs, header=None, index_col=0)
   df

Note how the parser automatically picks column names X.<column number> when
``header=None`` argument is specified. Alternatively, you can supply just the
column widths for contiguous columns:

.. ipython:: python

   # Widths are a list of integers
   widths = [6, 14, 13, 10]
   df = pd.read_fwf("bar.csv", widths=widths, header=None)
   df

The parser will take care of extra white spaces around the columns
so it's ok to have extra separation between the columns in the file.

By default, ``read_fwf`` will try to infer the file's ``colspecs`` by using the
first 100 rows of the file. It can do it only in cases when the columns are
aligned and correctly separated by the provided ``delimiter`` (default delimiter
is whitespace).

.. ipython:: python

   df = pd.read_fwf("bar.csv", header=None, index_col=0)
   df

``read_fwf`` supports the ``dtype`` parameter for specifying the types of
parsed columns to be different from the inferred type.

.. ipython:: python

   pd.read_fwf("bar.csv", header=None, index_col=0).dtypes
   pd.read_fwf("bar.csv", header=None, dtype={2: "object"}).dtypes

.. ipython:: python
   :suppress:

   os.remove("bar.csv")


Indexes
'''''''

Files with an "implicit" index column
+++++++++++++++++++++++++++++++++++++

Consider a file with one less entry in the header than the number of data
column:

.. ipython:: python

   data = "A,B,C\n20090101,a,1,2\n20090102,b,3,4\n20090103,c,4,5"
   print(data)
   with open("foo.csv", "w") as f:
       f.write(data)

In this special case, ``read_csv`` assumes that the first column is to be used
as the index of the ``DataFrame``:

.. ipython:: python

   pd.read_csv("foo.csv")

Note that the dates weren't automatically parsed. In that case you would need
to do as before:

.. ipython:: python

   df = pd.read_csv("foo.csv", parse_dates=True)
   df.index

.. ipython:: python
   :suppress:

   os.remove("foo.csv")


Reading an index with a ``MultiIndex``
++++++++++++++++++++++++++++++++++++++

.. _io.csv_multiindex:

Suppose you have data indexed by two columns:

.. ipython:: python

   data = 'year,indiv,zit,xit\n1977,"A",1.2,.6\n1977,"B",1.5,.5'
   print(data)
   with open("mindex_ex.csv", mode="w") as f:
       f.write(data)

The ``index_col`` argument to ``read_csv`` can take a list of
column numbers to turn multiple columns into a ``MultiIndex`` for the index of the
returned object:

.. ipython:: python

   df = pd.read_csv("mindex_ex.csv", index_col=[0, 1])
   df
   df.loc[1977]

.. ipython:: python
   :suppress:

   os.remove("mindex_ex.csv")

.. _io.multi_index_columns:

Reading columns with a ``MultiIndex``
+++++++++++++++++++++++++++++++++++++

By specifying list of row locations for the ``header`` argument, you
can read in a ``MultiIndex`` for the columns. Specifying non-consecutive
rows will skip the intervening rows.

.. ipython:: python

   mi_idx = pd.MultiIndex.from_arrays([[1, 2, 3, 4], list("abcd")], names=list("ab"))
   mi_col = pd.MultiIndex.from_arrays([[1, 2], list("ab")], names=list("cd"))
   df = pd.DataFrame(np.ones((4, 2)), index=mi_idx, columns=mi_col)
   df.to_csv("mi.csv")
   print(open("mi.csv").read())
   pd.read_csv("mi.csv", header=[0, 1, 2, 3], index_col=[0, 1])

``read_csv`` is also able to interpret a more common format
of multi-columns indices.

.. ipython:: python

   data = ",a,a,a,b,c,c\n,q,r,s,t,u,v\none,1,2,3,4,5,6\ntwo,7,8,9,10,11,12"
   print(data)
   with open("mi2.csv", "w") as fh:
       fh.write(data)

   pd.read_csv("mi2.csv", header=[0, 1], index_col=0)

.. note::
   If an ``index_col`` is not specified (e.g. you don't have an index, or wrote it
   with ``df.to_csv(..., index=False)``, then any ``names`` on the columns index will
   be *lost*.

.. ipython:: python
   :suppress:

   os.remove("mi.csv")
   os.remove("mi2.csv")

.. _io.sniff:

Automatically "sniffing" the delimiter
''''''''''''''''''''''''''''''''''''''

``read_csv`` is capable of inferring delimited (not necessarily
comma-separated) files, as pandas uses the :class:`python:csv.Sniffer`
class of the csv module. For this, you have to specify ``sep=None``.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 4))
   df.to_csv("tmp2.csv", sep=":", index=False)
   pd.read_csv("tmp2.csv", sep=None, engine="python")

.. ipython:: python
   :suppress:

   os.remove("tmp2.csv")

.. _io.multiple_files:

Reading multiple files to create a single DataFrame
'''''''''''''''''''''''''''''''''''''''''''''''''''

It's best to use :func:`~pandas.concat` to combine multiple files.
See the :ref:`cookbook<cookbook.csv.multiple_files>` for an example.

.. _io.chunking:

Iterating through files chunk by chunk
''''''''''''''''''''''''''''''''''''''

Suppose you wish to iterate through a (potentially very large) file lazily
rather than reading the entire file into memory, such as the following:


.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 4))
   df.to_csv("tmp.csv", index=False)
   table = pd.read_csv("tmp.csv")
   table


By specifying a ``chunksize`` to ``read_csv``, the return
value will be an iterable object of type ``TextFileReader``:

.. ipython:: python

   with pd.read_csv("tmp.csv", chunksize=4) as reader:
       print(reader)
       for chunk in reader:
           print(chunk)

.. versionchanged:: 1.2

  ``read_csv/json/sas`` return a context-manager when iterating through a file.

Specifying ``iterator=True`` will also return the ``TextFileReader`` object:

.. ipython:: python

   with pd.read_csv("tmp.csv", iterator=True) as reader:
       print(reader.get_chunk(5))

.. ipython:: python
   :suppress:

   os.remove("tmp.csv")

Specifying the parser engine
''''''''''''''''''''''''''''

pandas currently supports three engines, the C engine, the python engine, and an experimental
pyarrow engine (requires the ``pyarrow`` package). In general, the pyarrow engine is fastest
on larger workloads and is equivalent in speed to the C engine on most other workloads.
The python engine tends to be slower than the pyarrow and C engines on most workloads. However,
the pyarrow engine is much less robust than the C engine, which lacks a few features compared to the
Python engine.

Where possible, pandas uses the C parser (specified as ``engine='c'``), but it may fall
back to Python if C-unsupported options are specified.

Currently, options unsupported by the C and pyarrow engines include:

* ``sep`` other than a single character (e.g. regex separators)
* ``skipfooter``

Specifying any of the above options will produce a ``ParserWarning`` unless the
python engine is selected explicitly using ``engine='python'``.

Options that are unsupported by the pyarrow engine which are not covered by the list above include:

* ``float_precision``
* ``chunksize``
* ``comment``
* ``nrows``
* ``thousands``
* ``memory_map``
* ``dialect``
* ``on_bad_lines``
* ``quoting``
* ``lineterminator``
* ``converters``
* ``decimal``
* ``iterator``
* ``dayfirst``
* ``verbose``
* ``skipinitialspace``
* ``low_memory``

Specifying these options with ``engine='pyarrow'`` will raise a ``ValueError``.

.. _io.remote:

Reading/writing remote files
''''''''''''''''''''''''''''

You can pass in a URL to read or write remote files to many of pandas' IO
functions - the following example shows reading a CSV file:

.. code-block:: python

   df = pd.read_csv("https://download.bls.gov/pub/time.series/cu/cu.item", sep="\t")

.. versionadded:: 1.3.0

A custom header can be sent alongside HTTP(s) requests by passing a dictionary
of header key value mappings to the ``storage_options`` keyword argument as shown below:

.. code-block:: python

   headers = {"User-Agent": "pandas"}
   df = pd.read_csv(
       "https://download.bls.gov/pub/time.series/cu/cu.item",
       sep="\t",
       storage_options=headers
   )

All URLs which are not local files or HTTP(s) are handled by
`fsspec`_, if installed, and its various filesystem implementations
(including Amazon S3, Google Cloud, SSH, FTP, webHDFS...).
Some of these implementations will require additional packages to be
installed, for example
S3 URLs require the `s3fs
<https://pypi.org/project/s3fs/>`_ library:

.. code-block:: python

   df = pd.read_json("s3://pandas-test/adatafile.json")

When dealing with remote storage systems, you might need
extra configuration with environment variables or config files in
special locations. For example, to access data in your S3 bucket,
you will need to define credentials in one of the several ways listed in
the `S3Fs documentation
<https://s3fs.readthedocs.io/en/latest/#credentials>`_. The same is true
for several of the storage backends, and you should follow the links
at `fsimpl1`_ for implementations built into ``fsspec`` and `fsimpl2`_
for those not included in the main ``fsspec``
distribution.

You can also pass parameters directly to the backend driver. Since ``fsspec`` does not
utilize the ``AWS_S3_HOST`` environment variable, we can directly define a
dictionary containing the endpoint_url and pass the object into the storage
option parameter:

.. code-block:: python

   storage_options = {"client_kwargs": {"endpoint_url": "http://127.0.0.1:5555"}}
   df = pd.read_json("s3://pandas-test/test-1", storage_options=storage_options)

More sample configurations and documentation can be found at `S3Fs documentation
<https://s3fs.readthedocs.io/en/latest/index.html?highlight=host#s3-compatible-storage>`__.

If you do *not* have S3 credentials, you can still access public
data by specifying an anonymous connection, such as

.. versionadded:: 1.2.0

.. code-block:: python

   pd.read_csv(
       "s3://ncei-wcsd-archive/data/processed/SH1305/18kHz/SaKe2013"
       "-D20130523-T080854_to_SaKe2013-D20130523-T085643.csv",
       storage_options={"anon": True},
   )

``fsspec`` also allows complex URLs, for accessing data in compressed
archives, local caching of files, and more. To locally cache the above
example, you would modify the call to

.. code-block:: python

   pd.read_csv(
       "simplecache::s3://ncei-wcsd-archive/data/processed/SH1305/18kHz/"
       "SaKe2013-D20130523-T080854_to_SaKe2013-D20130523-T085643.csv",
       storage_options={"s3": {"anon": True}},
   )

where we specify that the "anon" parameter is meant for the "s3" part of
the implementation, not to the caching implementation. Note that this caches to a temporary
directory for the duration of the session only, but you can also specify
a permanent store.

.. _fsspec: https://filesystem-spec.readthedocs.io/en/latest/
.. _fsimpl1: https://filesystem-spec.readthedocs.io/en/latest/api.html#built-in-implementations
.. _fsimpl2: https://filesystem-spec.readthedocs.io/en/latest/api.html#other-known-implementations

Writing out data
''''''''''''''''

.. _io.store_in_csv:

Writing to CSV format
+++++++++++++++++++++

The ``Series`` and ``DataFrame`` objects have an instance method ``to_csv`` which
allows storing the contents of the object as a comma-separated-values file. The
function takes a number of arguments. Only the first is required.

* ``path_or_buf``: A string path to the file to write or a file object.  If a file object it must be opened with ``newline=''``
* ``sep`` : Field delimiter for the output file (default ",")
* ``na_rep``: A string representation of a missing value (default '')
* ``float_format``: Format string for floating point numbers
* ``columns``: Columns to write (default None)
* ``header``: Whether to write out the column names (default True)
* ``index``: whether to write row (index) names (default True)
* ``index_label``: Column label(s) for index column(s) if desired. If None
  (default), and ``header`` and ``index`` are True, then the index names are
  used. (A sequence should be given if the ``DataFrame`` uses MultiIndex).
* ``mode`` : Python write mode, default 'w'
* ``encoding``: a string representing the encoding to use if the contents are
  non-ASCII, for Python versions prior to 3
* ``lineterminator``: Character sequence denoting line end (default ``os.linesep``)
* ``quoting``: Set quoting rules as in csv module (default csv.QUOTE_MINIMAL). Note that if you have set a ``float_format`` then floats are converted to strings and csv.QUOTE_NONNUMERIC will treat them as non-numeric
* ``quotechar``: Character used to quote fields (default '"')
* ``doublequote``: Control quoting of ``quotechar`` in fields (default True)
* ``escapechar``: Character used to escape ``sep`` and ``quotechar`` when
  appropriate (default None)
* ``chunksize``: Number of rows to write at a time
* ``date_format``: Format string for datetime objects

Writing a formatted string
++++++++++++++++++++++++++

.. _io.formatting:

The ``DataFrame`` object has an instance method ``to_string`` which allows control
over the string representation of the object. All arguments are optional:

* ``buf`` default None, for example a StringIO object
* ``columns`` default None, which columns to write
* ``col_space`` default None, minimum width of each column.
* ``na_rep`` default ``NaN``, representation of NA value
* ``formatters`` default None, a dictionary (by column) of functions each of
  which takes a single argument and returns a formatted string
* ``float_format`` default None, a function which takes a single (float)
  argument and returns a formatted string; to be applied to floats in the
  ``DataFrame``.
* ``sparsify`` default True, set to False for a ``DataFrame`` with a hierarchical
  index to print every MultiIndex key at each row.
* ``index_names`` default True, will print the names of the indices
* ``index`` default True, will print the index (ie, row labels)
* ``header`` default True, will print the column labels
* ``justify`` default ``left``, will print column headers left- or
  right-justified

The ``Series`` object also has a ``to_string`` method, but with only the ``buf``,
``na_rep``, ``float_format`` arguments. There is also a ``length`` argument
which, if set to ``True``, will additionally output the length of the Series.
