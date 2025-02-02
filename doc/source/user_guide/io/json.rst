.. _io.json:

====
JSON
====

Read and write ``JSON`` format files and strings.

.. _io.json_writer:

Writing JSON
''''''''''''

A ``Series`` or ``DataFrame`` can be converted to a valid JSON string. Use ``to_json``
with optional parameters:

* ``path_or_buf`` : the pathname or buffer to write the output.
  This can be ``None`` in which case a JSON string is returned.
* ``orient`` :

  ``Series``:
      * default is ``index``
      * allowed values are {``split``, ``records``, ``index``}

  ``DataFrame``:
      * default is ``columns``
      * allowed values are {``split``, ``records``, ``index``, ``columns``, ``values``, ``table``}

  The format of the JSON string

  .. csv-table::
     :widths: 20, 150

     ``split``, dict like {index -> [index]; columns -> [columns]; data -> [values]}
     ``records``, list like [{column -> value}; ... ]
     ``index``, dict like {index -> {column -> value}}
     ``columns``, dict like {column -> {index -> value}}
     ``values``, just the values array
     ``table``, adhering to the JSON `Table Schema`_

* ``date_format`` : string, type of date conversion, 'epoch' for timestamp, 'iso' for ISO8601.
* ``double_precision`` : The number of decimal places to use when encoding floating point values, default 10.
* ``force_ascii`` : force encoded string to be ASCII, default True.
* ``date_unit`` : The time unit to encode to, governs timestamp and ISO8601 precision. One of 's', 'ms', 'us' or 'ns' for seconds, milliseconds, microseconds and nanoseconds respectively. Default 'ms'.
* ``default_handler`` : The handler to call if an object cannot otherwise be converted to a suitable format for JSON. Takes a single argument, which is the object to convert, and returns a serializable object.
* ``lines`` : If ``records`` orient, then will write each record per line as json.
* ``mode`` : string, writer mode when writing to path. 'w' for write, 'a' for append. Default 'w'

Note ``NaN``'s, ``NaT``'s and ``None`` will be converted to ``null`` and ``datetime`` objects will be converted based on the ``date_format`` and ``date_unit`` parameters.

.. ipython:: python

   dfj = pd.DataFrame(np.random.randn(5, 2), columns=list("AB"))
   json = dfj.to_json()
   json

Orient options
++++++++++++++

There are a number of different options for the format of the resulting JSON
file / string. Consider the following ``DataFrame`` and ``Series``:

.. ipython:: python

  dfjo = pd.DataFrame(
      dict(A=range(1, 4), B=range(4, 7), C=range(7, 10)),
      columns=list("ABC"),
      index=list("xyz"),
  )
  dfjo
  sjo = pd.Series(dict(x=15, y=16, z=17), name="D")
  sjo

**Column oriented** (the default for ``DataFrame``) serializes the data as
nested JSON objects with column labels acting as the primary index:

.. ipython:: python

  dfjo.to_json(orient="columns")
  # Not available for Series

**Index oriented** (the default for ``Series``) similar to column oriented
but the index labels are now primary:

.. ipython:: python

  dfjo.to_json(orient="index")
  sjo.to_json(orient="index")

**Record oriented** serializes the data to a JSON array of column -> value records,
index labels are not included. This is useful for passing ``DataFrame`` data to plotting
libraries, for example the JavaScript library ``d3.js``:

.. ipython:: python

  dfjo.to_json(orient="records")
  sjo.to_json(orient="records")

**Value oriented** is a bare-bones option which serializes to nested JSON arrays of
values only, column and index labels are not included:

.. ipython:: python

  dfjo.to_json(orient="values")
  # Not available for Series

**Split oriented** serializes to a JSON object containing separate entries for
values, index and columns. Name is also included for ``Series``:

.. ipython:: python

  dfjo.to_json(orient="split")
  sjo.to_json(orient="split")

**Table oriented** serializes to the JSON `Table Schema`_, allowing for the
preservation of metadata including but not limited to dtypes and index names.

.. note::

  Any orient option that encodes to a JSON object will not preserve the ordering of
  index and column labels during round-trip serialization. If you wish to preserve
  label ordering use the ``split`` option as it uses ordered containers.

Date handling
+++++++++++++

Writing in ISO date format:

.. ipython:: python

   dfd = pd.DataFrame(np.random.randn(5, 2), columns=list("AB"))
   dfd["date"] = pd.Timestamp("20130101")
   dfd = dfd.sort_index(axis=1, ascending=False)
   json = dfd.to_json(date_format="iso")
   json

Writing in ISO date format, with microseconds:

.. ipython:: python

   json = dfd.to_json(date_format="iso", date_unit="us")
   json

Writing to a file, with a date index and a date column:

.. ipython:: python

   dfj2 = dfj.copy()
   dfj2["date"] = pd.Timestamp("20130101")
   dfj2["ints"] = list(range(5))
   dfj2["bools"] = True
   dfj2.index = pd.date_range("20130101", periods=5)
   dfj2.to_json("test.json", date_format="iso")

   with open("test.json") as fh:
       print(fh.read())

Fallback behavior
+++++++++++++++++

If the JSON serializer cannot handle the container contents directly it will
fall back in the following manner:

* if the dtype is unsupported (e.g. ``np.complex_``) then the ``default_handler``, if provided, will be called
  for each value, otherwise an exception is raised.

* if an object is unsupported it will attempt the following:


    - check if the object has defined a ``toDict`` method and call it.
      A ``toDict`` method should return a ``dict`` which will then be JSON serialized.

    - invoke the ``default_handler`` if one was provided.

    - convert the object to a ``dict`` by traversing its contents. However this will often fail
      with an ``OverflowError`` or give unexpected results.

In general the best approach for unsupported objects or dtypes is to provide a ``default_handler``.
For example:

.. code-block:: python

  >>> DataFrame([1.0, 2.0, complex(1.0, 2.0)]).to_json()  # raises
  RuntimeError: Unhandled numpy dtype 15

can be dealt with by specifying a simple ``default_handler``:

.. ipython:: python

   pd.DataFrame([1.0, 2.0, complex(1.0, 2.0)]).to_json(default_handler=str)

.. _io.json_reader:

Reading JSON
''''''''''''

Reading a JSON string to pandas object can take a number of parameters.
The parser will try to parse a ``DataFrame`` if ``typ`` is not supplied or
is ``None``. To explicitly force ``Series`` parsing, pass ``typ=series``

* ``filepath_or_buffer`` : a **VALID** JSON string or file handle / StringIO. The string could be
  a URL. Valid URL schemes include http, ftp, S3, and file. For file URLs, a host
  is expected. For instance, a local file could be
  file ://localhost/path/to/table.json
* ``typ``    : type of object to recover (series or frame), default 'frame'
* ``orient`` :

  Series :
      * default is ``index``
      * allowed values are {``split``, ``records``, ``index``}

  DataFrame
      * default is ``columns``
      * allowed values are {``split``, ``records``, ``index``, ``columns``, ``values``, ``table``}

  The format of the JSON string

  .. csv-table::
     :widths: 20, 150

     ``split``, dict like {index -> [index]; columns -> [columns]; data -> [values]}
     ``records``, list like [{column -> value} ...]
     ``index``, dict like {index -> {column -> value}}
     ``columns``, dict like {column -> {index -> value}}
     ``values``, just the values array
     ``table``, adhering to the JSON `Table Schema`_


* ``dtype`` : if True, infer dtypes, if a dict of column to dtype, then use those, if ``False``, then don't infer dtypes at all, default is True, apply only to the data.
* ``convert_axes`` : boolean, try to convert the axes to the proper dtypes, default is ``True``
* ``convert_dates`` : a list of columns to parse for dates; If ``True``, then try to parse date-like columns, default is ``True``.
* ``keep_default_dates`` : boolean, default ``True``. If parsing dates, then parse the default date-like columns.
* ``precise_float`` : boolean, default ``False``. Set to enable usage of higher precision (strtod) function when decoding string to double values. Default (``False``) is to use fast but less precise builtin functionality.
* ``date_unit`` : string, the timestamp unit to detect if converting dates. Default
  None. By default the timestamp precision will be detected, if this is not desired
  then pass one of 's', 'ms', 'us' or 'ns' to force timestamp precision to
  seconds, milliseconds, microseconds or nanoseconds respectively.
* ``lines`` : reads file as one json object per line.
* ``encoding`` : The encoding to use to decode py3 bytes.
* ``chunksize`` : when used in combination with ``lines=True``, return a ``pandas.api.typing.JsonReader`` which reads in ``chunksize`` lines per iteration.
* ``engine``: Either ``"ujson"``, the built-in JSON parser, or ``"pyarrow"`` which dispatches to pyarrow's ``pyarrow.json.read_json``.
  The ``"pyarrow"`` is only available when ``lines=True``

The parser will raise one of ``ValueError/TypeError/AssertionError`` if the JSON is not parseable.

If a non-default ``orient`` was used when encoding to JSON be sure to pass the same
option here so that decoding produces sensible results, see `Orient Options`_ for an
overview.

Data conversion
+++++++++++++++

The default of ``convert_axes=True``, ``dtype=True``, and ``convert_dates=True``
will try to parse the axes, and all of the data into appropriate types,
including dates. If you need to override specific dtypes, pass a dict to
``dtype``. ``convert_axes`` should only be set to ``False`` if you need to
preserve string-like numbers (e.g. '1', '2') in an axes.

.. note::

  Large integer values may be converted to dates if ``convert_dates=True`` and the data and / or column labels appear 'date-like'. The exact threshold depends on the ``date_unit`` specified. 'date-like' means that the column label meets one of the following criteria:

  * it ends with ``'_at'``
  * it ends with ``'_time'``
  * it begins with ``'timestamp'``
  * it is ``'modified'``
  * it is ``'date'``

.. warning::

   When reading JSON data, automatic coercing into dtypes has some quirks:

   * an index can be reconstructed in a different order from serialization, that is, the returned order is not guaranteed to be the same as before serialization
   * a column that was ``float`` data will be converted to ``integer`` if it can be done safely, e.g. a column of ``1.``
   * bool columns will be converted to ``integer`` on reconstruction

   Thus there are times where you may want to specify specific dtypes via the ``dtype`` keyword argument.

Reading from a JSON string:

.. ipython:: python

   from io import StringIO
   pd.read_json(StringIO(json))

Reading from a file:

.. ipython:: python

   pd.read_json("test.json")

Don't convert any data (but still convert axes and dates):

.. ipython:: python

   pd.read_json("test.json", dtype=object).dtypes

Specify dtypes for conversion:

.. ipython:: python

   pd.read_json("test.json", dtype={"A": "float32", "bools": "int8"}).dtypes

Preserve string indices:

.. ipython:: python

   from io import StringIO
   si = pd.DataFrame(
       np.zeros((4, 4)), columns=list(range(4)), index=[str(i) for i in range(4)]
   )
   si
   si.index
   si.columns
   json = si.to_json()

   sij = pd.read_json(StringIO(json), convert_axes=False)
   sij
   sij.index
   sij.columns

Dates written in nanoseconds need to be read back in nanoseconds:

.. ipython:: python

   from io import StringIO
   json = dfj2.to_json(date_format="iso", date_unit="ns")

   # Try to parse timestamps as milliseconds -> Won't Work
   dfju = pd.read_json(StringIO(json), date_unit="ms")
   dfju

   # Let pandas detect the correct precision
   dfju = pd.read_json(StringIO(json))
   dfju

   # Or specify that all timestamps are in nanoseconds
   dfju = pd.read_json(StringIO(json), date_unit="ns")
   dfju

By setting the ``dtype_backend`` argument you can control the default dtypes used for the resulting DataFrame.

.. ipython:: python

    from io import StringIO

    data = (
     '{"a":{"0":1,"1":3},"b":{"0":2.5,"1":4.5},"c":{"0":true,"1":false},"d":{"0":"a","1":"b"},'
     '"e":{"0":null,"1":6.0},"f":{"0":null,"1":7.5},"g":{"0":null,"1":true},"h":{"0":null,"1":"a"},'
     '"i":{"0":"12-31-2019","1":"12-31-2019"},"j":{"0":null,"1":null}}'
    )
    df = pd.read_json(StringIO(data), dtype_backend="pyarrow")
    df
    df.dtypes

.. _io.json_normalize:

Normalization
'''''''''''''

pandas provides a utility function to take a dict or list of dicts and *normalize* this semi-structured data
into a flat table.

.. ipython:: python

   data = [
       {"id": 1, "name": {"first": "Coleen", "last": "Volk"}},
       {"name": {"given": "Mark", "family": "Regner"}},
       {"id": 2, "name": "Faye Raker"},
   ]
   pd.json_normalize(data)

.. ipython:: python

   data = [
       {
           "state": "Florida",
           "shortname": "FL",
           "info": {"governor": "Rick Scott"},
           "county": [
               {"name": "Dade", "population": 12345},
               {"name": "Broward", "population": 40000},
               {"name": "Palm Beach", "population": 60000},
           ],
       },
       {
           "state": "Ohio",
           "shortname": "OH",
           "info": {"governor": "John Kasich"},
           "county": [
               {"name": "Summit", "population": 1234},
               {"name": "Cuyahoga", "population": 1337},
           ],
       },
   ]

   pd.json_normalize(data, "county", ["state", "shortname", ["info", "governor"]])

The max_level parameter provides more control over which level to end normalization.
With max_level=1 the following snippet normalizes until 1st nesting level of the provided dict.

.. ipython:: python

    data = [
        {
            "CreatedBy": {"Name": "User001"},
            "Lookup": {
                "TextField": "Some text",
                "UserField": {"Id": "ID001", "Name": "Name001"},
            },
            "Image": {"a": "b"},
        }
    ]
    pd.json_normalize(data, max_level=1)

.. _io.jsonl:

Line delimited json
'''''''''''''''''''

pandas is able to read and write line-delimited json files that are common in data processing pipelines
using Hadoop or Spark.

For line-delimited json files, pandas can also return an iterator which reads in ``chunksize`` lines at a time. This can be useful for large files or to read from a stream.

.. ipython:: python

  from io import StringIO
  jsonl = """
      {"a": 1, "b": 2}
      {"a": 3, "b": 4}
  """
  df = pd.read_json(StringIO(jsonl), lines=True)
  df
  df.to_json(orient="records", lines=True)

  # reader is an iterator that returns ``chunksize`` lines each iteration
  with pd.read_json(StringIO(jsonl), lines=True, chunksize=1) as reader:
      reader
      for chunk in reader:
          print(chunk)

Line-limited json can also be read using the pyarrow reader by specifying ``engine="pyarrow"``.

.. ipython:: python

   from io import BytesIO
   df = pd.read_json(BytesIO(jsonl.encode()), lines=True, engine="pyarrow")
   df

.. versionadded:: 2.0.0

.. _io.table_schema:

Table schema
''''''''''''

`Table Schema`_ is a spec for describing tabular datasets as a JSON
object. The JSON includes information on the field names, types, and
other attributes. You can use the orient ``table`` to build
a JSON string with two fields, ``schema`` and ``data``.

.. ipython:: python

   df = pd.DataFrame(
       {
           "A": [1, 2, 3],
           "B": ["a", "b", "c"],
           "C": pd.date_range("2016-01-01", freq="D", periods=3),
       },
       index=pd.Index(range(3), name="idx"),
   )
   df
   df.to_json(orient="table", date_format="iso")

The ``schema`` field contains the ``fields`` key, which itself contains
a list of column name to type pairs, including the ``Index`` or ``MultiIndex``
(see below for a list of types).
The ``schema`` field also contains a ``primaryKey`` field if the (Multi)index
is unique.

The second field, ``data``, contains the serialized data with the ``records``
orient.
The index is included, and any datetimes are ISO 8601 formatted, as required
by the Table Schema spec.

The full list of types supported are described in the Table Schema
spec. This table shows the mapping from pandas types:

=============== =================
pandas type     Table Schema type
=============== =================
int64           integer
float64         number
bool            boolean
datetime64[ns]  datetime
timedelta64[ns] duration
categorical     any
object          str
=============== =================

A few notes on the generated table schema:

* The ``schema`` object contains a ``pandas_version`` field. This contains
  the version of pandas' dialect of the schema, and will be incremented
  with each revision.
* All dates are converted to UTC when serializing. Even timezone naive values,
  which are treated as UTC with an offset of 0.

  .. ipython:: python

     from pandas.io.json import build_table_schema

     s = pd.Series(pd.date_range("2016", periods=4))
     build_table_schema(s)

* datetimes with a timezone (before serializing), include an additional field
  ``tz`` with the time zone name (e.g. ``'US/Central'``).

  .. ipython:: python

     s_tz = pd.Series(pd.date_range("2016", periods=12, tz="US/Central"))
     build_table_schema(s_tz)

* Periods are converted to timestamps before serialization, and so have the
  same behavior of being converted to UTC. In addition, periods will contain
  and additional field ``freq`` with the period's frequency, e.g. ``'A-DEC'``.

  .. ipython:: python

     s_per = pd.Series(1, index=pd.period_range("2016", freq="Y-DEC", periods=4))
     build_table_schema(s_per)

* Categoricals use the ``any`` type and an ``enum`` constraint listing
  the set of possible values. Additionally, an ``ordered`` field is included:

  .. ipython:: python

     s_cat = pd.Series(pd.Categorical(["a", "b", "a"]))
     build_table_schema(s_cat)

* A ``primaryKey`` field, containing an array of labels, is included
  *if the index is unique*:

  .. ipython:: python

     s_dupe = pd.Series([1, 2], index=[1, 1])
     build_table_schema(s_dupe)

* The ``primaryKey`` behavior is the same with MultiIndexes, but in this
  case the ``primaryKey`` is an array:

  .. ipython:: python

     s_multi = pd.Series(1, index=pd.MultiIndex.from_product([("a", "b"), (0, 1)]))
     build_table_schema(s_multi)

* The default naming roughly follows these rules:

    - For series, the ``object.name`` is used. If that's none, then the
      name is ``values``
    - For ``DataFrames``, the stringified version of the column name is used
    - For ``Index`` (not ``MultiIndex``), ``index.name`` is used, with a
      fallback to ``index`` if that is None.
    - For ``MultiIndex``, ``mi.names`` is used. If any level has no name,
      then ``level_<i>`` is used.

``read_json`` also accepts ``orient='table'`` as an argument. This allows for
the preservation of metadata such as dtypes and index names in a
round-trippable manner.

.. ipython:: python

   df = pd.DataFrame(
       {
           "foo": [1, 2, 3, 4],
           "bar": ["a", "b", "c", "d"],
           "baz": pd.date_range("2018-01-01", freq="D", periods=4),
           "qux": pd.Categorical(["a", "b", "c", "c"]),
       },
       index=pd.Index(range(4), name="idx"),
   )
   df
   df.dtypes

   df.to_json("test.json", orient="table")
   new_df = pd.read_json("test.json", orient="table")
   new_df
   new_df.dtypes

Please note that the literal string 'index' as the name of an :class:`Index`
is not round-trippable, nor are any names beginning with ``'level_'`` within a
:class:`MultiIndex`. These are used by default in :func:`DataFrame.to_json` to
indicate missing values and the subsequent read cannot distinguish the intent.

.. ipython:: python
   :okwarning:

   df.index.name = "index"
   df.to_json("test.json", orient="table")
   new_df = pd.read_json("test.json", orient="table")
   print(new_df.index.name)

.. ipython:: python
   :suppress:

   os.remove("test.json")

When using ``orient='table'`` along with user-defined ``ExtensionArray``,
the generated schema will contain an additional ``extDtype`` key in the respective
``fields`` element. This extra key is not standard but does enable JSON roundtrips
for extension types (e.g. ``read_json(df.to_json(orient="table"), orient="table")``).

The ``extDtype`` key carries the name of the extension, if you have properly registered
the ``ExtensionDtype``, pandas will use said name to perform a lookup into the registry
and re-convert the serialized data into your custom dtype.

.. _Table Schema: https://specs.frictionlessdata.io/table-schema/
