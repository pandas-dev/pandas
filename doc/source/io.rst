.. _io:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import os
   import csv
   from StringIO import StringIO
   import pandas as pd
   ExcelWriter = pd.ExcelWriter

   import numpy as np
   np.random.seed(123456)
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

   import matplotlib.pyplot as plt
   plt.close('all')

   from pandas import *
   import pandas.util.testing as tm
   clipdf = DataFrame({'A':[1,2,3],'B':[4,5,6],'C':['p','q','r']},
                      index=['x','y','z'])

*******************************
IO Tools (Text, CSV, HDF5, ...)
*******************************

The Pandas I/O api is a set of top level ``reader`` functions accessed like ``pd.read_csv()`` that generally return a ``pandas``
object.

    * ``read_csv``
    * ``read_excel``
    * ``read_hdf``
    * ``read_sql``
    * ``read_json``
    * ``read_html``
    * ``read_stata``
    * ``read_clipboard``
    * ``read_pickle``

The corresponding ``writer`` functions are object methods that are accessed like ``df.to_csv()``

    * ``to_csv``
    * ``to_excel``
    * ``to_hdf``
    * ``to_sql``
    * ``to_json``
    * ``to_html``
    * ``to_stata``
    * ``to_clipboard``
    * ``to_pickle``

.. _io.read_csv_table:

CSV & Text files
----------------

The two workhorse functions for reading text files (a.k.a. flat files) are
:func:`~pandas.io.parsers.read_csv` and :func:`~pandas.io.parsers.read_table`.
They both use the same parsing code to intelligently convert tabular
data into a DataFrame object. See the :ref:`cookbook<cookbook.csv>`
for some advanced strategies

They can take a number of arguments:

  - ``filepath_or_buffer``: Either a string path to a file, url
    (including http, ftp, and s3 locations), or any object with a ``read``
    method (such as an open file or ``StringIO``).
  - ``sep`` or ``delimiter``: A delimiter / separator to split fields
    on. `read_csv` is capable of inferring the delimiter automatically in some
    cases by "sniffing." The separator may be specified as a regular
    expression; for instance you may use '\|\\s*' to indicate a pipe plus
    arbitrary whitespace.
  - ``delim_whitespace``: Parse whitespace-delimited (spaces or tabs) file
    (much faster than using a regular expression)
  - ``compression``: decompress ``'gzip'`` and ``'bz2'`` formats on the fly.
  - ``dialect``: string or :class:`python:csv.Dialect` instance to expose more
    ways to specify the file format
  - ``dtype``: A data type name or a dict of column name to data type. If not
    specified, data types will be inferred.
  - ``header``: row number to use as the column names, and the start of the
    data.  Defaults to 0 if no ``names`` passed, otherwise ``None``. Explicitly
    pass ``header=0`` to be able to replace existing names. The header can be
    a list of integers that specify row locations for a multi-index on the columns
    E.g. [0,1,3]. Interveaning rows that are not specified will be skipped.
    (E.g. 2 in this example are skipped)
  - ``skiprows``: A collection of numbers for rows in the file to skip. Can
    also be an integer to skip the first ``n`` rows
  - ``index_col``: column number, column name, or list of column numbers/names,
    to use as the ``index`` (row labels) of the resulting DataFrame. By default,
    it will number the rows without using any column, unless there is one more
    data column than there are headers, in which case the first column is taken
    as the index.
  - ``names``: List of column names to use as column names. To replace header
    existing in file, explicitly pass ``header=0``.
  - ``na_values``: optional list of strings to recognize as NaN (missing
    values), either in addition to or in lieu of the default set.
  - ``true_values``: list of strings to recognize as ``True``
  - ``false_values``: list of strings to recognize as ``False``
  - ``keep_default_na``: whether to include the default set of missing values
    in addition to the ones specified in ``na_values``
  - ``parse_dates``: if True then index will be parsed as dates
    (False by default). You can specify more complicated options to parse
    a subset of columns or a combination of columns into a single date column
    (list of ints or names, list of lists, or dict)
    [1, 2, 3] -> try parsing columns 1, 2, 3 each as a separate date column
    [[1, 3]] -> combine columns 1 and 3 and parse as a single date column
    {'foo' : [1, 3]} -> parse columns 1, 3 as date and call result 'foo'
  - ``keep_date_col``: if True, then date component columns passed into
    ``parse_dates`` will be retained in the output (False by default).
  - ``date_parser``: function to use to parse strings into datetime
    objects. If ``parse_dates`` is True, it defaults to the very robust
    ``dateutil.parser``. Specifying this implicitly sets ``parse_dates`` as True.
    You can also use functions from community supported date converters from
    date_converters.py
  - ``dayfirst``: if True then uses the DD/MM international/European date format
    (This is False by default)
  - ``thousands``: sepcifies the thousands separator. If not None, then parser
    will try to look for it in the output and parse relevant data to integers.
    Because it has to essentially scan through the data again, this causes a
    significant performance hit so only use if necessary.
  - ``lineterminator`` : string (length 1), default ``None``, Character to break file into lines. Only valid with C parser
  - ``quotechar`` : string, The character to used to denote the start and end of a quoted item.
    Quoted items can include the delimiter and it will be ignored.
  - ``quoting`` : int,
    Controls whether quotes should be recognized. Values are taken from `csv.QUOTE_*` values.
    Acceptable values are 0, 1, 2, and 3 for QUOTE_MINIMAL, QUOTE_ALL, QUOTE_NONE, and QUOTE_NONNUMERIC, respectively.
  - ``skipinitialspace`` : boolean, default ``False``, Skip spaces after delimiter
  - ``escapechar`` : string, to specify how to escape quoted data
  - ``comment``: denotes the start of a comment and ignores the rest of the line.
    Currently line commenting is not supported.
  - ``nrows``: Number of rows to read out of the file. Useful to only read a
    small portion of a large file
  - ``iterator``: If True, return a ``TextFileReader`` to enable reading a file
    into memory piece by piece
  - ``chunksize``: An number of rows to be used to "chunk" a file into
    pieces. Will cause an ``TextFileReader`` object to be returned. More on this
    below in the section on :ref:`iterating and chunking <io.chunking>`
  - ``skip_footer``: number of lines to skip at bottom of file (default 0)
  - ``converters``: a dictionary of functions for converting values in certain
    columns, where keys are either integers or column labels
  - ``encoding``: a string representing the encoding to use for decoding
    unicode data, e.g. ``'utf-8``` or ``'latin-1'``.
  - ``verbose``: show number of NA values inserted in non-numeric columns
  - ``squeeze``: if True then output with only one column is turned into Series
  - ``error_bad_lines``: if False then any lines causing an error will be skipped :ref:`bad lines <io.bad_lines>`
  - ``usecols``: a subset of columns to return, results in much faster parsing
    time and lower memory usage.
  - ``mangle_dupe_cols``: boolean, default True, then duplicate columns will be specified
    as 'X.0'...'X.N', rather than 'X'...'X'
  - ``tupleize_cols``: boolean, default True, if False, convert a list of tuples
    to a multi-index of columns, otherwise, leave the column index as a list of tuples

.. ipython:: python
   :suppress:

   f = open('foo.csv','w')
   f.write('date,A,B,C\n20090101,a,1,2\n20090102,b,3,4\n20090103,c,4,5')
   f.close()

Consider a typical CSV file containing, in this case, some time series data:

.. ipython:: python

   print open('foo.csv').read()

The default for `read_csv` is to create a DataFrame with simple numbered rows:

.. ipython:: python

   pd.read_csv('foo.csv')

In the case of indexed data, you can pass the column number or column name you
wish to use as the index:

.. ipython:: python

   pd.read_csv('foo.csv', index_col=0)

.. ipython:: python

   pd.read_csv('foo.csv', index_col='date')

You can also use a list of columns to create a hierarchical index:

.. ipython:: python

   pd.read_csv('foo.csv', index_col=[0, 'A'])

.. _io.dialect:

The ``dialect`` keyword gives greater flexibility in specifying the file format.
By default it uses the Excel dialect but you can specify either the dialect name
or a :class:`python:csv.Dialect` instance.

.. ipython:: python
   :suppress:

   data = ('label1,label2,label3\n'
           'index1,"a,c,e\n'
           'index2,b,d,f')

Suppose you had data with unenclosed quotes:

.. ipython:: python

   print data

By default, ``read_csv`` uses the Excel dialect and treats the double quote as
the quote character, which causes it to fail when it finds a newline before it
finds the closing double quote.

We can get around this using ``dialect``

.. ipython:: python

   dia = csv.excel()
   dia.quoting = csv.QUOTE_NONE
   pd.read_csv(StringIO(data), dialect=dia)

All of the dialect options can be specified separately by keyword arguments:

.. ipython:: python

    data = 'a,b,c~1,2,3~4,5,6'
    pd.read_csv(StringIO(data), lineterminator='~')

Another common dialect option is ``skipinitialspace``, to skip any whitespace
after a delimiter:

.. ipython:: python

   data = 'a, b, c\n1, 2, 3\n4, 5, 6'
   print data
   pd.read_csv(StringIO(data), skipinitialspace=True)

The parsers make every attempt to "do the right thing" and not be very
fragile. Type inference is a pretty big deal. So if a column can be coerced to
integer dtype without altering the contents, it will do so. Any non-numeric
columns will come through as object dtype as with the rest of pandas objects.

.. _io.dtypes:

Specifying column data types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Starting with v0.10, you can indicate the data type for the whole DataFrame or
individual columns:

.. ipython:: python

    data = 'a,b,c\n1,2,3\n4,5,6\n7,8,9'
    print data

    df = pd.read_csv(StringIO(data), dtype=object)
    df
    df['a'][0]
    df = pd.read_csv(StringIO(data), dtype={'b': object, 'c': np.float64})
    df.dtypes

.. _io.headers:

Handling column names
~~~~~~~~~~~~~~~~~~~~~

A file may or may not have a header row. pandas assumes the first row should be
used as the column names:

.. ipython:: python

    from StringIO import StringIO
    data = 'a,b,c\n1,2,3\n4,5,6\n7,8,9'
    print data
    pd.read_csv(StringIO(data))

By specifying the ``names`` argument in conjunction with ``header`` you can
indicate other names to use and whether or not to throw away the header row (if
any):

.. ipython:: python

    print data
    pd.read_csv(StringIO(data), names=['foo', 'bar', 'baz'], header=0)
    pd.read_csv(StringIO(data), names=['foo', 'bar', 'baz'], header=None)

If the header is in a row other than the first, pass the row number to
``header``. This will skip the preceding rows:

.. ipython:: python

    data = 'skip this skip it\na,b,c\n1,2,3\n4,5,6\n7,8,9'
    pd.read_csv(StringIO(data), header=1)

.. _io.usecols:

Filtering columns (``usecols``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``usecols`` argument allows you to select any subset of the columns in a
file, either using the column names or position numbers:

.. ipython:: python

    data = 'a,b,c,d\n1,2,3,foo\n4,5,6,bar\n7,8,9,baz'
    pd.read_csv(StringIO(data))
    pd.read_csv(StringIO(data), usecols=['b', 'd'])
    pd.read_csv(StringIO(data), usecols=[0, 2, 3])

.. _io.unicode:

Dealing with Unicode Data
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``encoding`` argument should be used for encoded unicode data, which will
result in byte strings being decoded to unicode in the result:

.. ipython:: python

   data = 'word,length\nTr\xe4umen,7\nGr\xfc\xdfe,5'
   df = pd.read_csv(StringIO(data), encoding='latin-1')
   df
   df['word'][1]

Some formats which encode all characters as multiple bytes, like UTF-16, won't
parse correctly at all without specifying the encoding.

.. _io.index_col:

Index columns and trailing delimiters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a file has one more column of data than the number of column names, the
first column will be used as the DataFrame's row names:

.. ipython:: python

    data = 'a,b,c\n4,apple,bat,5.7\n8,orange,cow,10'
    pd.read_csv(StringIO(data))

.. ipython:: python

    data = 'index,a,b,c\n4,apple,bat,5.7\n8,orange,cow,10'
    pd.read_csv(StringIO(data), index_col=0)

Ordinarily, you can achieve this behavior using the ``index_col`` option.

There are some exception cases when a file has been prepared with delimiters at
the end of each data line, confusing the parser. To explicitly disable the
index column inference and discard the last column, pass ``index_col=False``:

.. ipython:: python

    data = 'a,b,c\n4,apple,bat,\n8,orange,cow,'
    print data
    pd.read_csv(StringIO(data))
    pd.read_csv(StringIO(data), index_col=False)

.. _io.parse_dates:

Specifying Date Columns
~~~~~~~~~~~~~~~~~~~~~~~

To better facilitate working with datetime data,
:func:`~pandas.io.parsers.read_csv` and :func:`~pandas.io.parsers.read_table`
uses the keyword arguments ``parse_dates`` and ``date_parser`` to allow users
to specify a variety of columns and date/time formats to turn the input text
data into ``datetime`` objects.

The simplest case is to just pass in ``parse_dates=True``:

.. ipython:: python

   # Use a column as an index, and parse it as dates.
   df = pd.read_csv('foo.csv', index_col=0, parse_dates=True)
   df

   # These are python datetime objects
   df.index

.. ipython:: python
   :suppress:

   os.remove('foo.csv')

It is often the case that we may want to store date and time data separately,
or store various date fields separately. the ``parse_dates`` keyword can be
used to specify a combination of columns to parse the dates and/or times from.

You can specify a list of column lists to ``parse_dates``, the resulting date
columns will be prepended to the output (so as to not affect the existing column
order) and the new column names will be the concatenation of the component
column names:

.. ipython:: python
   :suppress:

   data =  ("KORD,19990127, 19:00:00, 18:56:00, 0.8100\n"
            "KORD,19990127, 20:00:00, 19:56:00, 0.0100\n"
            "KORD,19990127, 21:00:00, 20:56:00, -0.5900\n"
            "KORD,19990127, 21:00:00, 21:18:00, -0.9900\n"
            "KORD,19990127, 22:00:00, 21:56:00, -0.5900\n"
            "KORD,19990127, 23:00:00, 22:56:00, -0.5900")

   with open('tmp.csv', 'w') as fh:
       fh.write(data)

.. ipython:: python

    print open('tmp.csv').read()
    df = pd.read_csv('tmp.csv', header=None, parse_dates=[[1, 2], [1, 3]])
    df

By default the parser removes the component date columns, but you can choose
to retain them via the ``keep_date_col`` keyword:

.. ipython:: python

   df = pd.read_csv('tmp.csv', header=None, parse_dates=[[1, 2], [1, 3]],
                    keep_date_col=True)
   df

Note that if you wish to combine multiple columns into a single date column, a
nested list must be used. In other words, ``parse_dates=[1, 2]`` indicates that
the second and third columns should each be parsed as separate date columns
while ``parse_dates=[[1, 2]]`` means the two columns should be parsed into a
single column.

You can also use a dict to specify custom name columns:

.. ipython:: python

   date_spec = {'nominal': [1, 2], 'actual': [1, 3]}
   df = pd.read_csv('tmp.csv', header=None, parse_dates=date_spec)
   df

It is important to remember that if multiple text columns are to be parsed into
a single date column, then a new column is prepended to the data. The `index_col`
specification is based off of this new set of columns rather than the original
data columns:


.. ipython:: python

   date_spec = {'nominal': [1, 2], 'actual': [1, 3]}
   df = pd.read_csv('tmp.csv', header=None, parse_dates=date_spec,
                    index_col=0) #index is the nominal column
   df

**Note**: When passing a dict as the `parse_dates` argument, the order of
the columns prepended is not guaranteed, because `dict` objects do not impose
an ordering on their keys. On Python 2.7+ you may use `collections.OrderedDict`
instead of a regular `dict` if this matters to you. Because of this, when using a
dict for 'parse_dates' in conjunction with the `index_col` argument, it's best to
specify `index_col` as a column label rather then as an index on the resulting frame.

Date Parsing Functions
~~~~~~~~~~~~~~~~~~~~~~
Finally, the parser allows you can specify a custom ``date_parser`` function to
take full advantage of the flexiblity of the date parsing API:

.. ipython:: python

   import pandas.io.date_converters as conv
   df = pd.read_csv('tmp.csv', header=None, parse_dates=date_spec,
                    date_parser=conv.parse_date_time)
   df

You can explore the date parsing functionality in ``date_converters.py`` and
add your own. We would love to turn this module into a community supported set
of date/time parsers. To get you started, ``date_converters.py`` contains
functions to parse dual date and time columns, year/month/day columns,
and year/month/day/hour/minute/second columns. It also contains a
``generic_parser`` function so you can curry it with a function that deals with
a single date rather than the entire array.

.. ipython:: python
   :suppress:

   os.remove('tmp.csv')

.. _io.dayfirst:

International Date Formats
~~~~~~~~~~~~~~~~~~~~~~~~~~
While US date formats tend to be MM/DD/YYYY, many international formats use
DD/MM/YYYY instead. For convenience, a ``dayfirst`` keyword is provided:

.. ipython:: python
   :suppress:

   data = "date,value,cat\n1/6/2000,5,a\n2/6/2000,10,b\n3/6/2000,15,c"
   with open('tmp.csv', 'w') as fh:
        fh.write(data)

.. ipython:: python

   print open('tmp.csv').read()

   pd.read_csv('tmp.csv', parse_dates=[0])
   pd.read_csv('tmp.csv', dayfirst=True, parse_dates=[0])

.. _io.thousands:

Thousand Separators
~~~~~~~~~~~~~~~~~~~
For large integers that have been written with a thousands separator, you can
set the ``thousands`` keyword to ``True`` so that integers will be parsed
correctly:

.. ipython:: python
   :suppress:

   data =  ("ID|level|category\n"
            "Patient1|123,000|x\n"
            "Patient2|23,000|y\n"
            "Patient3|1,234,018|z")

   with open('tmp.csv', 'w') as fh:
       fh.write(data)

By default, integers with a thousands separator will be parsed as strings

.. ipython:: python

    print open('tmp.csv').read()
    df = pd.read_csv('tmp.csv', sep='|')
    df

    df.level.dtype

The ``thousands`` keyword allows integers to be parsed correctly

.. ipython:: python

    print open('tmp.csv').read()
    df = pd.read_csv('tmp.csv', sep='|', thousands=',')
    df

    df.level.dtype

.. ipython:: python
   :suppress:

   os.remove('tmp.csv')

.. _io.comments:

Comments
~~~~~~~~
Sometimes comments or meta data may be included in a file:

.. ipython:: python
   :suppress:

   data =  ("ID,level,category\n"
            "Patient1,123000,x # really unpleasant\n"
            "Patient2,23000,y # wouldn't take his medicine\n"
            "Patient3,1234018,z # awesome")

   with open('tmp.csv', 'w') as fh:
       fh.write(data)

.. ipython:: python

   print open('tmp.csv').read()

By default, the parse includes the comments in the output:

.. ipython:: python

   df = pd.read_csv('tmp.csv')
   df

We can suppress the comments using the ``comment`` keyword:

.. ipython:: python

   df = pd.read_csv('tmp.csv', comment='#')
   df

.. ipython:: python
   :suppress:

   os.remove('tmp.csv')

Returning Series
~~~~~~~~~~~~~~~~

Using the ``squeeze`` keyword, the parser will return output with a single column
as a ``Series``:

.. ipython:: python
   :suppress:

   data =  ("level\n"
            "Patient1,123000\n"
            "Patient2,23000\n"
            "Patient3,1234018")

   with open('tmp.csv', 'w') as fh:
       fh.write(data)

.. ipython:: python

   print open('tmp.csv').read()

   output =  pd.read_csv('tmp.csv', squeeze=True)
   output

   type(output)

.. ipython:: python
   :suppress:

   os.remove('tmp.csv')

.. _io.boolean:

Boolean values
~~~~~~~~~~~~~~

The common values ``True``, ``False``, ``TRUE``, and ``FALSE`` are all
recognized as boolean. Sometime you would want to recognize some other values
as being boolean. To do this use the ``true_values`` and ``false_values``
options:

.. ipython:: python

    data= 'a,b,c\n1,Yes,2\n3,No,4'
    print data
    pd.read_csv(StringIO(data))
    pd.read_csv(StringIO(data), true_values=['Yes'], false_values=['No'])

.. _io.bad_lines:

Handling "bad" lines
~~~~~~~~~~~~~~~~~~~~

Some files may have malformed lines with too few fields or too many. Lines with
too few fields will have NA values filled in the trailing fields. Lines with
too many will cause an error by default:

.. ipython:: python
   :suppress:

    data = 'a,b,c\n1,2,3\n4,5,6,7\n8,9,10'

.. code-block:: ipython

    In [27]: data = 'a,b,c\n1,2,3\n4,5,6,7\n8,9,10'

    In [28]: pd.read_csv(StringIO(data))
    ---------------------------------------------------------------------------
    CParserError                              Traceback (most recent call last)
    CParserError: Error tokenizing data. C error: Expected 3 fields in line 3, saw 4

You can elect to skip bad lines:

.. code-block:: ipython

    In [29]: pd.read_csv(StringIO(data), error_bad_lines=False)
    Skipping line 3: expected 3 fields, saw 4

    Out[29]:
       a  b   c
    0  1  2   3
    1  8  9  10

.. _io.quoting:

Quoting and Escape Characters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quotes (and other escape characters) in embedded fields can be handled in any
number of ways. One way is to use backslashes; to properly parse this data, you
should pass the ``escapechar`` option:

.. ipython:: python

   data = 'a,b\n"hello, \\"Bob\\", nice to see you",5'
   print data
   pd.read_csv(StringIO(data), escapechar='\\')

.. _io.fwf:

Files with Fixed Width Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
While ``read_csv`` reads delimited data, the :func:`~pandas.io.parsers.read_fwf`
function works with data files that have known and fixed column widths.
The function parameters to ``read_fwf`` are largely the same as `read_csv` with
two extra parameters:

  - ``colspecs``: a list of pairs (tuples), giving the extents of the
    fixed-width fields of each line as half-open intervals [from, to[
  - ``widths``: a list of field widths, which can be used instead of
    ``colspecs`` if the intervals are contiguous

.. ipython:: python
   :suppress:

   f = open('bar.csv', 'w')
   data1 = ("id8141    360.242940   149.910199   11950.7\n"
            "id1594    444.953632   166.985655   11788.4\n"
            "id1849    364.136849   183.628767   11806.2\n"
            "id1230    413.836124   184.375703   11916.8\n"
            "id1948    502.953953   173.237159   12468.3")
   f.write(data1)
   f.close()

Consider a typical fixed-width data file:

.. ipython:: python

   print open('bar.csv').read()

In order to parse this file into a DataFrame, we simply need to supply the
column specifications to the `read_fwf` function along with the file name:

.. ipython:: python

   #Column specifications are a list of half-intervals
   colspecs = [(0, 6), (8, 20), (21, 33), (34, 43)]
   df = pd.read_fwf('bar.csv', colspecs=colspecs, header=None, index_col=0)
   df

Note how the parser automatically picks column names X.<column number> when
``header=None`` argument is specified. Alternatively, you can supply just the
column widths for contiguous columns:

.. ipython:: python

   #Widths are a list of integers
   widths = [6, 14, 13, 10]
   df = pd.read_fwf('bar.csv', widths=widths, header=None)
   df

The parser will take care of extra white spaces around the columns
so it's ok to have extra separation between the columns in the file.

.. ipython:: python
   :suppress:

   os.remove('bar.csv')

Files with an "implicit" index column
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   :suppress:

   f = open('foo.csv', 'w')
   f.write('A,B,C\n20090101,a,1,2\n20090102,b,3,4\n20090103,c,4,5')
   f.close()

Consider a file with one less entry in the header than the number of data
column:

.. ipython:: python

   print open('foo.csv').read()

In this special case, ``read_csv`` assumes that the first column is to be used
as the index of the DataFrame:

.. ipython:: python

   pd.read_csv('foo.csv')

Note that the dates weren't automatically parsed. In that case you would need
to do as before:

.. ipython:: python

   df = pd.read_csv('foo.csv', parse_dates=True)
   df.index

.. ipython:: python
   :suppress:

   os.remove('foo.csv')


Reading an index with a ``MultiIndex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _io.csv_multiindex:

Suppose you have data indexed by two columns:

.. ipython:: python

   print open('data/mindex_ex.csv').read()

The ``index_col`` argument to ``read_csv`` and ``read_table`` can take a list of
column numbers to turn multiple columns into a ``MultiIndex`` for the index of the
returned object:

.. ipython:: python

   df = pd.read_csv("data/mindex_ex.csv", index_col=[0,1])
   df
   df.ix[1978]

.. _io.multi_index_columns:

Reading columns with a ``MultiIndex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By specifying list of row locations for the ``header`` argument, you
can read in a ``MultiIndex`` for the columns. Specifying non-consecutive
rows will skip the interveaning rows.

.. ipython:: python

   from pandas.util.testing import makeCustomDataframe as mkdf
   df = mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
   df.to_csv('mi.csv',tupleize_cols=False)
   print open('mi.csv').read()
   pd.read_csv('mi.csv',header=[0,1,2,3],index_col=[0,1],tupleize_cols=False)

Note: The default behavior in 0.12 remains unchanged (``tupleize_cols=True``) from prior versions,
but starting with 0.13, the default *to* write and read multi-index columns will be in the new
format (``tupleize_cols=False``)

Note: If an ``index_col`` is not specified (e.g. you don't have an index, or wrote it
with ``df.to_csv(..., index=False``), then any ``names`` on the columns index will be *lost*.

.. ipython:: python
   :suppress:

   import os
   os.remove('mi.csv')

.. _io.sniff:

Automatically "sniffing" the delimiter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``read_csv`` is capable of inferring delimited (not necessarily
comma-separated) files. YMMV, as pandas uses the :class:`python:csv.Sniffer`
class of the csv module.

.. ipython:: python
   :suppress:

   df = DataFrame(np.random.randn(10, 4))
   df.to_csv('tmp.sv', sep='|')
   df.to_csv('tmp2.sv', sep=':')

.. ipython:: python

    print open('tmp2.sv').read()
    pd.read_csv('tmp2.sv')

.. _io.chunking:

Iterating through files chunk by chunk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you wish to iterate through a (potentially very large) file lazily
rather than reading the entire file into memory, such as the following:


.. ipython:: python

   print open('tmp.sv').read()
   table = pd.read_table('tmp.sv', sep='|')
   table


By specifiying a ``chunksize`` to ``read_csv`` or ``read_table``, the return
value will be an iterable object of type ``TextFileReader``:

.. ipython:: python

   reader = pd.read_table('tmp.sv', sep='|', chunksize=4)
   reader

   for chunk in reader:
       print chunk


Specifying ``iterator=True`` will also return the ``TextFileReader`` object:

.. ipython:: python

   reader = pd.read_table('tmp.sv', sep='|', iterator=True)
   reader.get_chunk(5)

.. ipython:: python
   :suppress:

   os.remove('tmp.sv')
   os.remove('tmp2.sv')

Writing to CSV format
~~~~~~~~~~~~~~~~~~~~~

.. _io.store_in_csv:

The Series and DataFrame objects have an instance method ``to_csv`` which
allows storing the contents of the object as a comma-separated-values file. The
function takes a number of arguments. Only the first is required.

  - ``path``: A string path to the file to write
  - ``na_rep``: A string representation of a missing value (default '')
  - ``cols``: Columns to write (default None)
  - ``header``: Whether to write out the column names (default True)
  - ``index``: whether to write row (index) names (default True)
  - ``index_label``: Column label(s) for index column(s) if desired. If None
    (default), and `header` and `index` are True, then the index names are
    used. (A sequence should be given if the DataFrame uses MultiIndex).
  - ``mode`` : Python write mode, default 'w'
  - ``sep`` : Field delimiter for the output file (default ",")
  - ``encoding``: a string representing the encoding to use if the contents are
    non-ascii, for python versions prior to 3
  - ``tupleize_cols``: boolean, default True, if False, write as a list of tuples,
    otherwise write in an expanded line format suitable for ``read_csv``

Writing a formatted string
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _io.formatting:

The DataFrame object has an instance method ``to_string`` which allows control
over the string representation of the object. All arguments are optional:

  - ``buf`` default None, for example a StringIO object
  - ``columns`` default None, which columns to write
  - ``col_space`` default None, minimum width of each column.
  - ``na_rep`` default ``NaN``, representation of NA value
  - ``formatters`` default None, a dictionary (by column) of functions each of
    which takes a single argument and returns a formatted string
  - ``float_format`` default None, a function which takes a single (float)
    argument and returns a formatted string; to be applied to floats in the
    DataFrame.
  - ``sparsify`` default True, set to False for a DataFrame with a hierarchical
    index to print every multiindex key at each row.
  - ``index_names`` default True, will print the names of the indices
  - ``index`` default True, will print the index (ie, row labels)
  - ``header`` default True, will print the column labels
  - ``justify`` default ``left``, will print column headers left- or
    right-justified

The Series object also has a ``to_string`` method, but with only the ``buf``,
``na_rep``, ``float_format`` arguments. There is also a ``length`` argument
which, if set to ``True``, will additionally output the length of the Series.


JSON
----

Read and write ``JSON`` format files.

.. _io.json:

Writing JSON
~~~~~~~~~~~~

A ``Series`` or ``DataFrame`` can be converted to a valid JSON string. Use ``to_json``
with optional parameters:

- ``path_or_buf`` : the pathname or buffer to write the output
  This can be ``None`` in which case a JSON string is returned
- ``orient`` :

  Series :
      - default is ``index``
      - allowed values are {``split``, ``records``, ``index``}

  DataFrame
      - default is ``columns``
      - allowed values are {``split``, ``records``, ``index``, ``columns``, ``values``}

  The format of the JSON string

  .. csv-table::
     :widths: 20, 150
     :delim: ;

     ``split``; dict like {index -> [index], columns -> [columns], data -> [values]}
     ``records``; list like [{column -> value}, ... , {column -> value}]
     ``index``; dict like {index -> {column -> value}}
     ``columns``; dict like {column -> {index -> value}}
     ``values``; just the values array

- ``date_format`` : type of date conversion (epoch = epoch milliseconds, iso = ISO8601), default is epoch
- ``double_precision`` : The number of decimal places to use when encoding floating point values, default 10.
- ``force_ascii`` : force encoded string to be ASCII, default True.

Note NaN's and None will be converted to null and datetime objects will be converted based on the date_format parameter

.. ipython:: python

   dfj = DataFrame(randn(5, 2), columns=list('AB'))
   json = dfj.to_json()
   json

Writing in iso date format

.. ipython:: python

   dfd = DataFrame(randn(5, 2), columns=list('AB'))
   dfd['date'] = Timestamp('20130101')
   json = dfd.to_json(date_format='iso')
   json

Writing to a file, with a date index and a date column

.. ipython:: python

   dfj2 = dfj.copy()
   dfj2['date'] = Timestamp('20130101')
   dfj2['ints'] = range(5)
   dfj2['bools'] = True
   dfj2.index = date_range('20130101',periods=5)
   dfj2.to_json('test.json')
   open('test.json').read()

Reading JSON
~~~~~~~~~~~~

Reading a JSON string to pandas object can take a number of parameters.
The parser will try to parse a ``DataFrame`` if ``typ`` is not supplied or
is ``None``. To explicity force ``Series`` parsing, pass ``typ=series``

- ``filepath_or_buffer`` : a **VALID** JSON string or file handle / StringIO. The string could be
  a URL. Valid URL schemes include http, ftp, s3, and file. For file URLs, a host
  is expected. For instance, a local file could be
  file ://localhost/path/to/table.json
- ``typ``    : type of object to recover (series or frame), default 'frame'
- ``orient`` :

  Series :
      - default is ``index``
      - allowed values are {``split``, ``records``, ``index``}

  DataFrame
      - default is ``columns``
      - allowed values are {``split``, ``records``, ``index``, ``columns``, ``values``}

  The format of the JSON string

  .. csv-table::
     :widths: 20, 150
     :delim: ;

     ``split``; dict like {index -> [index], columns -> [columns], data -> [values]}
     ``records``; list like [{column -> value}, ... , {column -> value}]
     ``index``; dict like {index -> {column -> value}}
     ``columns``; dict like {column -> {index -> value}}
     ``values``; just the values array

- ``dtype`` : if True, infer dtypes, if a dict of column to dtype, then use those, if False, then don't infer dtypes at all, default is True, apply only to the data
- ``convert_axes`` : boolean, try to convert the axes to the proper dtypes, default is True
- ``convert_dates`` : a list of columns to parse for dates; If True, then try to parse datelike columns, default is True
- ``keep_default_dates`` : boolean, default True. If parsing dates, then parse the default datelike columns
- ``numpy`` : direct decoding to numpy arrays. default is False;
  Note that the JSON ordering **MUST** be the same for each term if ``numpy=True``
- ``precise_float`` : boolean, default ``False``. Set to enable usage of higher precision (strtod) function
  when decoding string to double values. Default (``False``) is to use fast but less precise builtin functionality

The parser will raise one of ``ValueError/TypeError/AssertionError`` if the JSON is
not parsable.

The default of ``convert_axes=True``, ``dtype=True``, and ``convert_dates=True`` will try to parse the axes, and all of the data
into appropriate types, including dates. If you need to override specific dtypes, pass a dict to ``dtype``. ``convert_axes`` should only
be set to ``False`` if you need to preserve string-like numbers (e.g. '1', '2') in an axes.

.. warning::

   When reading JSON data, automatic coercing into dtypes has some quirks:

     * an index can be reconstructed in a different order from serialization, that is, the returned order is not guaranteed to be the same as before serialization
     * a column that was ``float`` data will be converted to ``integer`` if it can be done safely, e.g. a column of ``1.``
     * bool columns will be converted to ``integer`` on reconstruction

   Thus there are times where you may want to specify specific dtypes via the ``dtype`` keyword argument.

Reading from a JSON string

.. ipython:: python

   pd.read_json(json)

Reading from a file

.. ipython:: python

   pd.read_json('test.json')

Don't convert any data (but still convert axes and dates)

.. ipython:: python

   pd.read_json('test.json',dtype=object).dtypes

Specify how I want to convert data

.. ipython:: python

   pd.read_json('test.json',dtype={'A' : 'float32', 'bools' : 'int8'}).dtypes

I like my string indicies

.. ipython:: python

   si = DataFrame(np.zeros((4, 4)),
            columns=range(4),
            index=[str(i) for i in range(4)])
   si
   si.index
   si.columns
   json = si.to_json()

   sij = pd.read_json(json,convert_axes=False)
   sij
   sij.index
   sij.columns

.. ipython:: python
   :suppress:

   import os
   os.remove('test.json')

HTML
----

Reading HTML Content
~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   We **highly encourage** you to read the :ref:`HTML parsing gotchas
   <html-gotchas>` regarding the issues surrounding the
   BeautifulSoup4/html5lib/lxml parsers.

.. _io.read_html:

.. versionadded:: 0.12

The top-level :func:`~pandas.io.html.read_html` function can accept an HTML
string/file/url and will parse HTML tables into list of pandas DataFrames.
Let's look at a few examples.

.. note::

   ``read_html`` returns a ``list`` of ``DataFrame`` objects, even if there is
   only a single table contained in the HTML content

Read a URL with no options

.. ipython:: python

   url = 'http://www.fdic.gov/bank/individual/failed/banklist.html'
   dfs = read_html(url)
   dfs

.. note::

   The data from the above URL changes every Monday so the resulting data above
   and the data below may be slightly different.

Read in the content of the file from the above URL and pass it to ``read_html``
as a string

.. ipython:: python
   :suppress:

   import os
   file_path = os.path.abspath(os.path.join('source', '_static', 'banklist.html'))

.. ipython:: python

   with open(file_path, 'r') as f:
       dfs = read_html(f.read())
   dfs

You can even pass in an instance of ``StringIO`` if you so desire

.. ipython:: python

   from cStringIO import StringIO

   with open(file_path, 'r') as f:
       sio = StringIO(f.read())

   dfs = read_html(sio)
   dfs

.. note::

   The following examples are not run by the IPython evaluator due to the fact
   that having so many network-accessing functions slows down the documentation
   build. If you spot an error or an example that doesn't run, please do not
   hesitate to report it over on `pandas GitHub issues page
   <http://www.github.com/pydata/pandas/issues>`__.


Read a URL and match a table that contains specific text

.. code-block:: python

   match = 'Metcalf Bank'
   df_list = read_html(url, match=match)

Specify a header row (by default ``<th>`` elements are used to form the column
index); if specified, the header row is taken from the data minus the parsed
header elements (``<th>`` elements).

.. code-block:: python

   dfs = read_html(url, header=0)

Specify an index column

.. code-block:: python

   dfs = read_html(url, index_col=0)

Specify a number of rows to skip

.. code-block:: python

   dfs = read_html(url, skiprows=0)

Specify a number of rows to skip using a list (``xrange`` (Python 2 only) works
as well)

.. code-block:: python

   dfs = read_html(url, skiprows=range(2))

Don't infer numeric and date types

.. code-block:: python

   dfs = read_html(url, infer_types=False)

Specify an HTML attribute

.. code-block:: python

   dfs1 = read_html(url, attrs={'id': 'table'})
   dfs2 = read_html(url, attrs={'class': 'sortable'})
   print np.array_equal(dfs1[0], dfs2[0])  # Should be True

Use some combination of the above

.. code-block:: python

   dfs = read_html(url, match='Metcalf Bank', index_col=0)

Read in pandas ``to_html`` output (with some loss of floating point precision)

.. code-block:: python

   df = DataFrame(randn(2, 2))
   s = df.to_html(float_format='{0:.40g}'.format)
   dfin = read_html(s, index_col=0)

The ``lxml`` backend will raise an error on a failed parse if that is the only
parser you provide (if you only have a single parser you can provide just a
string, but it is considered good practice to pass a list with one string if,
for example, the function expects a sequence of strings)

.. code-block:: python

   dfs = read_html(url, 'Metcalf Bank', index_col=0, flavor=['lxml'])

or

.. code-block:: python

   dfs = read_html(url, 'Metcalf Bank', index_col=0, flavor='lxml')

However, if you have bs4 and html5lib installed and pass ``None`` or ``['lxml',
'bs4']`` then the parse will most likely succeed. Note that *as soon as a parse
succeeds, the function will return*.

.. code-block:: python

   dfs = read_html(url, 'Metcalf Bank', index_col=0, flavor=['lxml', 'bs4'])


Writing to HTML files
~~~~~~~~~~~~~~~~~~~~~~

.. _io.html:

``DataFrame`` objects have an instance method ``to_html`` which renders the
contents of the ``DataFrame`` as an HTML table. The function arguments are as
in the method ``to_string`` described above.

.. note::

   Not all of the possible options for ``DataFrame.to_html`` are shown here for
   brevity's sake. See :func:`~pandas.core.frame.DataFrame.to_html` for the
   full set of options.

.. ipython:: python
   :suppress:

   def write_html(df, filename, *args, **kwargs):
       static = os.path.abspath(os.path.join('source', '_static'))
       with open(os.path.join(static, filename + '.html'), 'w') as f:
           df.to_html(f, *args, **kwargs)

.. ipython:: python

   df = DataFrame(randn(2, 2))
   df
   print df.to_html()  # raw html

.. ipython:: python
   :suppress:

   write_html(df, 'basic')

HTML:

.. raw:: html
   :file: _static/basic.html

The ``columns`` argument will limit the columns shown

.. ipython:: python

   print df.to_html(columns=[0])

.. ipython:: python
   :suppress:

   write_html(df, 'columns', columns=[0])

HTML:

.. raw:: html
   :file: _static/columns.html

``float_format`` takes a Python callable to control the precision of floating
point values

.. ipython:: python

   print df.to_html(float_format='{0:.10f}'.format)

.. ipython:: python
   :suppress:

   write_html(df, 'float_format', float_format='{0:.10f}'.format)

HTML:

.. raw:: html
   :file: _static/float_format.html

``bold_rows`` will make the row labels bold by default, but you can turn that
off

.. ipython:: python

   print df.to_html(bold_rows=False)

.. ipython:: python
   :suppress:

   write_html(df, 'nobold', bold_rows=False)

.. raw:: html
   :file: _static/nobold.html

The ``classes`` argument provides the ability to give the resulting HTML
table CSS classes. Note that these classes are *appended* to the existing
``'dataframe'`` class.

.. ipython:: python

   print df.to_html(classes=['awesome_table_class', 'even_more_awesome_class'])

Finally, the ``escape`` argument allows you to control whether the
"<", ">" and "&" characters escaped in the resulting HTML (by default it is
``True``). So to get the HTML without escaped characters pass ``escape=False``

.. ipython:: python

   df = DataFrame({'a': list('&<>'), 'b': randn(3)})


.. ipython:: python
   :suppress:

   write_html(df, 'escape')
   write_html(df, 'noescape', escape=False)

Escaped:

.. ipython:: python

   print df.to_html()

.. raw:: html
   :file: _static/escape.html

Not escaped:

.. ipython:: python

   print df.to_html(escape=False)

.. raw:: html
   :file: _static/noescape.html

.. note::

   Some browsers may not show a difference in the rendering of the previous two
   HTML tables.


Clipboard
---------

.. _io.clipboard:

A handy way to grab data is to use the ``read_clipboard`` method, which takes
the contents of the clipboard buffer and passes them to the ``read_table``
method. For instance, you can copy the following
text to the clipboard (CTRL-C on many operating systems):

.. code-block:: python

     A B C
   x 1 4 p
   y 2 5 q
   z 3 6 r

And then import the data directly to a DataFrame by calling:

.. code-block:: python

   clipdf = pd.read_clipboard()

.. ipython:: python

   clipdf

The ``to_clipboard`` method can be used to write the contents of a DataFrame to
the clipboard. Following which you can paste the clipboard contents into other
applications (CTRL-V on many operating systems). Here we illustrate writing a
DataFrame into clipboard and reading it back.

.. ipython:: python

    df=pd.DataFrame(randn(5,3))
    df
    df.to_clipboard()
    pd.read_clipboard()

We can see that we got the same content back, which we had earlier written to the clipboard.

.. note::

   You may need to install xclip or xsel (with gtk or PyQt4 modules) on Linux to use these methods.


.. _io.serialize:

Pickling and serialization
--------------------------

All pandas objects are equipped with ``to_pickle`` methods which use Python's
``cPickle`` module to save data structures to disk using the pickle format.

.. ipython:: python

   df
   df.to_pickle('foo.pkl')

The ``read_pickle`` function in the ``pandas`` namespace can be used to load
any pickled pandas object (or any other pickled object) from file:


.. ipython:: python

   read_pickle('foo.pkl')

.. ipython:: python
   :suppress:

   import os
   os.remove('foo.pkl')

.. warning::

   Loading pickled data received from untrusted sources can be unsafe.

   See: http://docs.python.org/2.7/library/pickle.html

.. note::

    These methods were previously ``save`` and ``load``, now deprecated.

.. _io.excel:

Excel files
-----------

The ``read_excel`` method can read Excel 2003 (``.xls``) and
Excel 2007 (``.xlsx``) files using the ``xlrd`` Python
module and use the same parsing code as the above to convert tabular data into
a DataFrame. See the :ref:`cookbook<cookbook.excel>` for some
advanced strategies

.. note::

   The prior method of accessing Excel is now deprecated as of 0.12,
   this will work but will be removed in a future version.

      .. code-block:: python

         from pandas.io.parsers import ExcelFile
         xls = ExcelFile('path_to_file.xls')
         xls.parse('Sheet1', index_col=None, na_values=['NA'])

   Replaced by

     .. code-block:: python

        read_excel('path_to_file.xls', 'Sheet1', index_col=None, na_values=['NA'])

It is often the case that users will insert columns to do temporary computations
in Excel and you may not want to read in those columns. `read_excel` takes
a `parse_cols` keyword to allow you to specify a subset of columns to parse.

If `parse_cols` is an integer, then it is assumed to indicate the last column
to be parsed.

.. code-block:: python

   read_excel('path_to_file.xls', 'Sheet1', parse_cols=2, index_col=None, na_values=['NA'])

If `parse_cols` is a list of integers, then it is assumed to be the file column
indices to be parsed.

.. code-block:: python

   read_excel('path_to_file.xls', Sheet1', parse_cols=[0, 2, 3], index_col=None, na_values=['NA'])

To write a DataFrame object to a sheet of an Excel file, you can use the
``to_excel`` instance method.  The arguments are largely the same as ``to_csv``
described above, the first argument being the name of the excel file, and the
optional second argument the name of the sheet to which the DataFrame should be
written.  For example:

.. code-block:: python

   df.to_excel('path_to_file.xlsx', sheet_name='sheet1')

Files with a ``.xls`` extension will be written using ``xlwt`` and those with
a ``.xlsx`` extension will be written using ``openpyxl``.
The Panel class also has a ``to_excel`` instance method,
which writes each DataFrame in the Panel to a separate sheet.

In order to write separate DataFrames to separate sheets in a single Excel file,
one can use the ExcelWriter class, as in the following example:

.. code-block:: python

   writer = ExcelWriter('path_to_file.xlsx')
   df1.to_excel(writer, sheet_name='sheet1')
   df2.to_excel(writer, sheet_name='sheet2')
   writer.save()

.. _io.hdf5:

HDF5 (PyTables)
---------------

``HDFStore`` is a dict-like object which reads and writes pandas using
the high performance HDF5 format using the excellent `PyTables
<http://www.pytables.org/>`__ library. See the :ref:`cookbook<cookbook.hdf>`
for some advanced strategies

.. note::

   ``PyTables`` 3.0.0 was recently released to enables support for Python 3.
   Pandas should be fully compatible (and previously written stores should be
   backwards compatible) with all ``PyTables`` >= 2.3

.. ipython:: python
   :suppress:
   :okexcept:

   os.remove('store.h5')

.. ipython:: python

   store = HDFStore('store.h5')
   print store

Objects can be written to the file just like adding key-value pairs to a
dict:

.. ipython:: python

   index = date_range('1/1/2000', periods=8)
   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   df = DataFrame(randn(8, 3), index=index,
                  columns=['A', 'B', 'C'])
   wp = Panel(randn(2, 5, 4), items=['Item1', 'Item2'],
              major_axis=date_range('1/1/2000', periods=5),
              minor_axis=['A', 'B', 'C', 'D'])

   # store.put('s', s) is an equivalent method
   store['s'] = s

   store['df'] = df

   store['wp'] = wp

   # the type of stored data
   store.root.wp._v_attrs.pandas_type

   store

In a current or later Python session, you can retrieve stored objects:

.. ipython:: python

   # store.get('df') is an equivalent method
   store['df']

   # dotted (attribute) access provides get as well
   store.df

Deletion of the object specified by the key

.. ipython:: python

   # store.remove('wp') is an equivalent method
   del store['wp']

   store

Closing a Store, Context Manager

.. ipython:: python

   # closing a store
   store.close()

   # Working with, and automatically closing the store with the context
   # manager
   with get_store('store.h5') as store:
        store.keys()

.. ipython:: python
   :suppress:

   store.close()
   import os
   os.remove('store.h5')

Read/Write API
~~~~~~~~~~~~~~

``HDFStore`` supports an top-level API using  ``read_hdf`` for reading and ``to_hdf`` for writing,
similar to how ``read_csv`` and ``to_csv`` work. (new in 0.11.0)

.. ipython:: python

   df_tl = DataFrame(dict(A=range(5), B=range(5)))
   df_tl.to_hdf('store_tl.h5','table',append=True)
   read_hdf('store_tl.h5', 'table', where = ['index>2'])

.. ipython:: python
   :suppress:
   :okexcept:

   os.remove('store_tl.h5')

.. _io.hdf5-storer:

Storer Format
~~~~~~~~~~~~~

The examples above show storing using ``put``, which write the HDF5 to ``PyTables`` in a fixed array format, called
the ``storer`` format. These types of stores are are **not** appendable once written (though you can simply
remove them and rewrite). Nor are they **queryable**; they must be
retrieved in their entirety. These offer very fast writing and slightly faster reading than ``table`` stores.

.. warning::

   A ``storer`` format will raise a ``TypeError`` if you try to retrieve using a ``where`` .

   .. code-block:: python

       DataFrame(randn(10,2)).to_hdf('test_storer.h5','df')

       pd.read_hdf('test_storer.h5','df',where='index>5')
       TypeError: cannot pass a where specification when reading a non-table
                  this store must be selected in its entirety


.. _io.hdf5-table:

Table Format
~~~~~~~~~~~~

``HDFStore`` supports another ``PyTables`` format on disk, the ``table``
format. Conceptually a ``table`` is shaped very much like a DataFrame,
with rows and columns. A ``table`` may be appended to in the same or
other sessions.  In addition, delete & query type operations are
supported.

.. ipython:: python
   :suppress:
   :okexcept:

   os.remove('store.h5')

.. ipython:: python

   store = HDFStore('store.h5')
   df1 = df[0:4]
   df2 = df[4:]

   # append data (creates a table automatically)
   store.append('df', df1)
   store.append('df', df2)
   store

   # select the entire object
   store.select('df')

   # the type of stored data
   store.root.df._v_attrs.pandas_type

.. note::

   You can also create a ``table`` by passing ``table=True`` to a ``put`` operation.

.. _io.hdf5-keys:

Hierarchical Keys
~~~~~~~~~~~~~~~~~

Keys to a store can be specified as a string. These can be in a
hierarchical path-name like format (e.g. ``foo/bar/bah``), which will
generate a hierarchy of sub-stores (or ``Groups`` in PyTables
parlance). Keys can be specified with out the leading '/' and are ALWAYS
absolute (e.g. 'foo' refers to '/foo'). Removal operations can remove
everying in the sub-store and BELOW, so be *careful*.

.. ipython:: python

   store.put('foo/bar/bah', df)
   store.append('food/orange', df)
   store.append('food/apple',  df)
   store

   # a list of keys are returned
   store.keys()

   # remove all nodes under this level
   store.remove('food')
   store

.. _io.hdf5-types:

Storing Mixed Types in a Table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Storing mixed-dtype data is supported. Strings are stored as a
fixed-width using the maximum size of the appended column. Subsequent
appends will truncate strings at this length.

Passing ``min_itemsize={`values`: size}`` as a parameter to append
will set a larger minimum for the string columns. Storing ``floats,
strings, ints, bools, datetime64`` are currently supported. For string
columns, passing ``nan_rep = 'nan'`` to append will change the default
nan representation on disk (which converts to/from `np.nan`), this
defaults to `nan`.

.. ipython:: python

    df_mixed = DataFrame({ 'A' : randn(8),
                           'B' : randn(8),
                           'C' : np.array(randn(8),dtype='float32'),
                           'string' :'string',
                           'int' : 1,
                           'bool' : True,
                           'datetime64' : Timestamp('20010102')},
                         index=range(8))
    df_mixed.ix[3:5,['A', 'B', 'string', 'datetime64']] = np.nan

    store.append('df_mixed', df_mixed, min_itemsize = {'values': 50})
    df_mixed1 = store.select('df_mixed')
    df_mixed1
    df_mixed1.get_dtype_counts()

    # we have provided a minimum string column size
    store.root.df_mixed.table

Storing Multi-Index DataFrames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Storing multi-index dataframes as tables is very similar to
storing/selecting from homogeneous index DataFrames.

.. ipython:: python

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['foo', 'bar'])
        df_mi = DataFrame(np.random.randn(10, 3), index=index,
                          columns=['A', 'B', 'C'])
        df_mi

        store.append('df_mi',df_mi)
        store.select('df_mi')

        # the levels are automatically included as data columns
        store.select('df_mi', Term('foo=bar'))


.. _io.hdf5-query:

Querying a Table
~~~~~~~~~~~~~~~~

``select`` and ``delete`` operations have an optional criterion that can
be specified to select/delete only a subset of the data. This allows one
to have a very large on-disk table and retrieve only a portion of the
data.

A query is specified using the ``Term`` class under the hood.

   - 'index' and 'columns' are supported indexers of a DataFrame
   - 'major_axis', 'minor_axis', and 'items' are supported indexers of
     the Panel

Valid terms can be created from ``dict, list, tuple, or
string``. Objects can be embeded as values. Allowed operations are: ``<,
<=, >, >=, =, !=``. ``=`` will be inferred as an implicit set operation
(e.g. if 2 or more values are provided). The following are all valid
terms.

       - ``dict(field = 'index', op = '>', value = '20121114')``
       - ``('index', '>', '20121114')``
       - ``'index > 20121114'``
       - ``('index', '>', datetime(2012, 11, 14))``
       - ``('index', ['20121114', '20121115'])``
       - ``('major_axis', '=', Timestamp('2012/11/14'))``
       - ``('minor_axis', ['A', 'B'])``

Queries are built up using a list of ``Terms`` (currently only
**anding** of terms is supported). An example query for a panel might be
specified as follows.  ``['major_axis>20000102', ('minor_axis', '=',
['A', 'B']) ]``. This is roughly translated to: `major_axis must be
greater than the date 20000102 and the minor_axis must be A or B`

.. ipython:: python

   store.append('wp',wp)
   store
   store.select('wp', [ Term('major_axis>20000102'), Term('minor_axis', '=', ['A', 'B']) ])

The ``columns`` keyword can be supplied to select a list of columns to be returned,
this is equivalent to passing a ``Term('columns', list_of_columns_to_filter)``:

.. ipython:: python

   store.select('df', columns=['A', 'B'])

``start`` and ``stop`` parameters can be specified to limit the total search
space. These are in terms of the total number of rows in a table.

.. ipython:: python

   # this is effectively what the storage of a Panel looks like
   wp.to_frame()

   # limiting the search
   store.select('wp',[ Term('major_axis>20000102'),
                       Term('minor_axis', '=', ['A','B']) ],
                start=0, stop=10)


Indexing
~~~~~~~~

You can create/modify an index for a table with ``create_table_index``
after data is already in the table (after and ``append/put``
operation). Creating a table index is **highly** encouraged. This will
speed your queries a great deal when you use a ``select`` with the
indexed dimension as the ``where``. **Indexes are automagically created
(starting 0.10.1)** on the indexables and any data columns you
specify. This behavior can be turned off by passing ``index=False`` to
``append``.

.. ipython:: python

   # we have automagically already created an index (in the first section)
   i = store.root.df.table.cols.index.index
   i.optlevel, i.kind

   # change an index by passing new parameters
   store.create_table_index('df', optlevel=9, kind='full')
   i = store.root.df.table.cols.index.index
   i.optlevel, i.kind


Query via Data Columns
~~~~~~~~~~~~~~~~~~~~~~

You can designate (and index) certain columns that you want to be able
to perform queries (other than the `indexable` columns, which you can
always query). For instance say you want to perform this common
operation, on-disk, and return just the frame that matches this
query. You can specify ``data_columns = True`` to force all columns to
be data_columns

.. ipython:: python

   df_dc = df.copy()
   df_dc['string'] = 'foo'
   df_dc.ix[4:6,'string'] = np.nan
   df_dc.ix[7:9,'string'] = 'bar'
   df_dc['string2'] = 'cool'
   df_dc

   # on-disk operations
   store.append('df_dc', df_dc, data_columns = ['B', 'C', 'string', 'string2'])
   store.select('df_dc', [ Term('B>0') ])

   # getting creative
   store.select('df_dc', ['B > 0', 'C > 0', 'string == foo'])

   # this is in-memory version of this type of selection
   df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (df_dc.string == 'foo')]

   # we have automagically created this index and the B/C/string/string2
   # columns are stored separately as ``PyTables`` columns
   store.root.df_dc.table

There is some performance degredation by making lots of columns into
`data columns`, so it is up to the user to designate these. In addition,
you cannot change data columns (nor indexables) after the first
append/put operation (Of course you can simply read in the data and
create a new table!)

Iterator
~~~~~~~~

Starting in 0.11, you can pass, ``iterator=True`` or ``chunksize=number_in_a_chunk``
to ``select`` and ``select_as_multiple`` to return an iterator on the results.
The default is 50,000 rows returned in a chunk.

.. ipython:: python

   for df in store.select('df', chunksize=3):
      print df

.. note::

   .. versionadded:: 0.12

   You can also use the iterator with ``read_hdf`` which will open, then
   automatically close the store when finished iterating.

   .. code-block:: python

      for df in read_hdf('store.h5','df', chunsize=3):
          print df

Note, that the chunksize keyword applies to the **returned** rows. So if you
are doing a query, then that set will be subdivided and returned in the
iterator. Keep in mind that if you do not pass a ``where`` selection criteria
then the ``nrows`` of the table are considered.

Advanced Queries
~~~~~~~~~~~~~~~~

**Select a Single Column**

To retrieve a single indexable or data column, use the
method ``select_column``. This will, for example, enable you to get the index
very quickly. These return a ``Series`` of the result, indexed by the row number.
These do not currently accept the ``where`` selector (coming soon)

.. ipython:: python

   store.select_column('df_dc', 'index')
   store.select_column('df_dc', 'string')

**Replicating or**

``not`` and ``or`` conditions are unsupported at this time; however,
``or`` operations are easy to replicate, by repeatedly applying the
criteria to the table, and then ``concat`` the results.

.. ipython:: python

   crit1 = [ Term('B>0'), Term('C>0'), Term('string=foo') ]
   crit2 = [ Term('B<0'), Term('C>0'), Term('string=foo') ]

   concat([store.select('df_dc',c) for c in [crit1, crit2]])

**Storer Object**

If you want to inspect the stored object, retrieve via
``get_storer``. You could use this programmatically to say get the number
of rows in an object.

.. ipython:: python

   store.get_storer('df_dc').nrows


Multiple Table Queries
~~~~~~~~~~~~~~~~~~~~~~

New in 0.10.1 are the methods ``append_to_multiple`` and
``select_as_multiple``, that can perform appending/selecting from
multiple tables at once. The idea is to have one table (call it the
selector table) that you index most/all of the columns, and perform your
queries. The other table(s) are data tables with an index matching the
selector table's index. You can then perform a very fast query
on the selector table, yet get lots of data back. This method works similar to
having a very wide table, but is more efficient in terms of queries.

Note, **THE USER IS RESPONSIBLE FOR SYNCHRONIZING THE TABLES**. This
means, append to the tables in the same order; ``append_to_multiple``
splits a single object to multiple tables, given a specification (as a
dictionary). This dictionary is a mapping of the table names to the
'columns' you want included in that table. Pass a `None` for a single
table (optional) to let it have the remaining columns. The argument
``selector`` defines which table is the selector table.

.. ipython:: python

   df_mt = DataFrame(randn(8, 6), index=date_range('1/1/2000', periods=8),
                                  columns=['A', 'B', 'C', 'D', 'E', 'F'])
   df_mt['foo'] = 'bar'

   # you can also create the tables individually
   store.append_to_multiple({'df1_mt': ['A', 'B'], 'df2_mt': None },
                             df_mt, selector='df1_mt')
   store

   # indiviual tables were created
   store.select('df1_mt')
   store.select('df2_mt')

   # as a multiple
   store.select_as_multiple(['df1_mt', 'df2_mt'], where=['A>0', 'B>0'],
                             selector = 'df1_mt')

.. _io.hdf5-delete:

Delete from a Table
~~~~~~~~~~~~~~~~~~~

You can delete from a table selectively by specifying a ``where``. In
deleting rows, it is important to understand the ``PyTables`` deletes
rows by erasing the rows, then **moving** the following data. Thus
deleting can potentially be a very expensive operation depending on the
orientation of your data. This is especially true in higher dimensional
objects (``Panel`` and ``Panel4D``). To get optimal performance, it's
worthwhile to have the dimension you are deleting be the first of the
``indexables``.

Data is ordered (on the disk) in terms of the ``indexables``. Here's a
simple use case. You store panel-type data, with dates in the
``major_axis`` and ids in the ``minor_axis``. The data is then
interleaved like this:

   - date_1
        - id_1
        - id_2
        -  .
        - id_n
   - date_2
        - id_1
        -  .
        - id_n

It should be clear that a delete operation on the ``major_axis`` will be
fairly quick, as one chunk is removed, then the following data moved. On
the other hand a delete operation on the ``minor_axis`` will be very
expensive. In this case it would almost certainly be faster to rewrite
the table using a ``where`` that selects all but the missing data.

.. ipython:: python

   # returns the number of rows deleted
   store.remove('wp', 'major_axis>20000102' )
   store.select('wp')

Please note that HDF5 **DOES NOT RECLAIM SPACE** in the h5 files
automatically. Thus, repeatedly deleting (or removing nodes) and adding
again **WILL TEND TO INCREASE THE FILE SIZE**. To *clean* the file, use
``ptrepack`` (see below).

Compression
~~~~~~~~~~~

``PyTables`` allows the stored data to be compressed. Tthis applies to
all kinds of stores, not just tables.

   - Pass ``complevel=int`` for a compression level (1-9, with 0 being no
     compression, and the default)
   - Pass ``complib=lib`` where lib is any of ``zlib, bzip2, lzo, blosc`` for
     whichever compression library you prefer.

``HDFStore`` will use the file based compression scheme if no overriding
``complib`` or ``complevel`` options are provided. ``blosc`` offers very
fast compression, and is my most used. Note that ``lzo`` and ``bzip2``
may not be installed (by Python) by default.

Compression for all objects within the file

   - ``store_compressed = HDFStore('store_compressed.h5', complevel=9, complib='blosc')``

Or on-the-fly compression (this only applies to tables). You can turn
off file compression for a specific table by passing ``complevel=0``

   - ``store.append('df', df, complib='zlib', complevel=5)``

**ptrepack**

``PyTables`` offers better write performance when tables are compressed after
they are written, as opposed to turning on compression at the very
beginning. You can use the supplied ``PyTables`` utility
``ptrepack``. In addition, ``ptrepack`` can change compression levels
after the fact.

   - ``ptrepack --chunkshape=auto --propindexes --complevel=9 --complib=blosc in.h5 out.h5``

Furthermore ``ptrepack in.h5 out.h5`` will *repack* the file to allow
you to reuse previously deleted space. Aalternatively, one can simply
remove the file and write again, or use the ``copy`` method.

Notes & Caveats
~~~~~~~~~~~~~~~

   - Once a ``table`` is created its items (Panel) / columns (DataFrame)
     are fixed; only exactly the same columns can be appended
   - If a row has ``np.nan`` for **EVERY COLUMN** (having a ``nan``
     in a string, or a ``NaT`` in a datetime-like column counts as having
     a value), then those rows **WILL BE DROPPED IMPLICITLY**. This limitation
     *may* be addressed in the future.
   - ``HDFStore`` is **not-threadsafe for writing**. The underlying
     ``PyTables`` only supports concurrent reads (via threading or
     processes). If you need reading and writing *at the same time*, you
     need to serialize these operations in a single thread in a single
     process. You will corrupt your data otherwise. See the issue
     <https://github.com/pydata/pandas/issues/2397> for more
     information.
   - ``PyTables`` only supports fixed-width string columns in
     ``tables``. The sizes of a string based indexing column
     (e.g. *columns* or *minor_axis*) are determined as the maximum size
     of the elements in that axis or by passing the parameter

DataTypes
~~~~~~~~~

``HDFStore`` will map an object dtype to the ``PyTables`` underlying
dtype. This means the following types are known to work:

    - floating : ``float64, float32, float16`` *(using* ``np.nan`` *to
      represent invalid values)*
    - integer : ``int64, int32, int8, uint64, uint32, uint8``
    - bool
    - datetime64[ns] *(using* ``NaT`` *to represent invalid values)*
    - object : ``strings`` *(using* ``np.nan`` *to represent invalid
      values)*

Currently, ``unicode`` and ``datetime`` columns (represented with a
dtype of ``object``), **WILL FAIL**. In addition, even though a column
may look like a ``datetime64[ns]``, if it contains ``np.nan``, this
**WILL FAIL**. You can try to convert datetimelike columns to proper
``datetime64[ns]`` columns, that possibily contain ``NaT`` to represent
invalid values. (Some of these issues have been addressed and these
conversion may not be necessary in future versions of pandas)

    .. ipython:: python

       import datetime
       df = DataFrame(dict(datelike=Series([datetime.datetime(2001, 1, 1),
                                            datetime.datetime(2001, 1, 2), np.nan])))
       df
       df.dtypes

       # to convert
       df['datelike'] = Series(df['datelike'].values, dtype='M8[ns]')
       df
       df.dtypes

String Columns
~~~~~~~~~~~~~~

The underlying implementation of ``HDFStore`` uses a fixed column width (itemsize) for string columns. A string column itemsize is calculated as the maximum of the
length of data (for that column) that is passed to the ``HDFStore``, **in the first append**. Subsequent appends, may introduce a string for a column **larger** than the column can hold, an Exception will be raised (otherwise you could have a silent truncation of these columns, leading to loss of information). In the future we may relax this and allow a user-specified truncation to occur.

Pass ``min_itemsize`` on the first table creation to a-priori specifiy the minimum length of a particular string column. ``min_itemsize`` can be an integer, or a dict mapping a column name to an integer. You can pass ``values`` as a key to allow all *indexables* or *data_columns* to have this min_itemsize.

Starting in 0.11, passing a ``min_itemsize`` dict will cause all passed columns to be created as *data_columns* automatically.

.. note::

   If you are not passing any *data_columns*, then the min_itemsize will be the maximum of the length of any string passed

.. ipython:: python

   dfs = DataFrame(dict(A = 'foo', B = 'bar'),index=range(5))
   dfs

   # A and B have a size of 30
   store.append('dfs', dfs, min_itemsize = 30)
   store.get_storer('dfs').table

   # A is created as a data_column with a size of 30
   # B is size is calculated
   store.append('dfs2', dfs, min_itemsize = { 'A' : 30 })
   store.get_storer('dfs2').table

External Compatibility
~~~~~~~~~~~~~~~~~~~~~~

``HDFStore`` write storer objects in specific formats suitable for
producing loss-less roundtrips to pandas objects. For external
compatibility, ``HDFStore`` can read native ``PyTables`` format
tables. It is possible to write an ``HDFStore`` object that can easily
be imported into ``R`` using the ``rhdf5`` library. Create a table
format store like this:

     .. ipython:: python

        store_export = HDFStore('export.h5')
	    store_export.append('df_dc', df_dc, data_columns=df_dc.columns)
	    store_export

     .. ipython:: python
        :suppress:

        store_export.close()
        import os
        os.remove('export.h5')

Backwards Compatibility
~~~~~~~~~~~~~~~~~~~~~~~

0.10.1 of ``HDFStore`` can read tables created in a prior version of pandas,
however query terms using the
prior (undocumented) methodology are unsupported. ``HDFStore`` will
issue a warning if you try to use a legacy-format file. You must
read in the entire file and write it out using the new format, using the
method ``copy`` to take advantage of the updates. The group attribute
``pandas_version`` contains the version information. ``copy`` takes a
number of options, please see the docstring.


     .. ipython:: python
        :suppress:

        import os
        legacy_file_path = os.path.abspath('source/_static/legacy_0.10.h5')

     .. ipython:: python

        # a legacy store
        legacy_store = HDFStore(legacy_file_path,'r')
        legacy_store

        # copy (and return the new handle)
	    new_store = legacy_store.copy('store_new.h5')
	    new_store
        new_store.close()

     .. ipython:: python
        :suppress:

        legacy_store.close()
        import os
        os.remove('store_new.h5')


Performance
~~~~~~~~~~~

   - ``Tables`` come with a writing performance penalty as compared to
     regular stores. The benefit is the ability to append/delete and
     query (potentially very large amounts of data).  Write times are
     generally longer as compared with regular stores. Query times can
     be quite fast, especially on an indexed axis.
   - You can pass ``chunksize=<int>`` to ``append``, specifying the
     write chunksize (default is 50000). This will signficantly lower
     your memory usage on writing.
   - You can pass ``expectedrows=<int>`` to the first ``append``,
     to set the TOTAL number of expected rows that ``PyTables`` will
     expected. This will optimize read/write performance.
   - Duplicate rows can be written to tables, but are filtered out in
     selection (with the last items being selected; thus a table is
     unique on major, minor pairs)
   - A ``PerformanceWarning`` will be raised if you are attempting to
     store types that will be pickled by PyTables (rather than stored as
     endemic types). See
     <http://stackoverflow.com/questions/14355151/how-to-make-pandas-hdfstore-put-operation-faster/14370190#14370190>
     for more information and some solutions.

Experimental
~~~~~~~~~~~~

HDFStore supports ``Panel4D`` storage.

.. ipython:: python

   p4d = Panel4D({ 'l1' : wp })
   p4d
   store.append('p4d', p4d)
   store

These, by default, index the three axes ``items, major_axis,
minor_axis``. On an ``AppendableTable`` it is possible to setup with the
first append a different indexing scheme, depending on how you want to
store your data. Pass the ``axes`` keyword with a list of dimensions
(currently must by exactly 1 less than the total dimensions of the
object). This cannot be changed after table creation.

.. ipython:: python

   store.append('p4d2', p4d, axes=['labels', 'major_axis', 'minor_axis'])
   store
   store.select('p4d2', [ Term('labels=l1'), Term('items=Item1'), Term('minor_axis=A_big_strings') ])

.. ipython:: python
   :suppress:

   store.close()
   import os
   os.remove('store.h5')


.. _io.sql:

SQL Queries
-----------

The :mod:`pandas.io.sql` module provides a collection of query wrappers to both
facilitate data retrieval and to reduce dependency on DB-specific API. These
wrappers only support the Python database adapters which respect the `Python
DB-API <http://www.python.org/dev/peps/pep-0249/>`_. See some
:ref:`cookbook examples <cookbook.sql>` for some advanced strategies

For example, suppose you want to query some data with different types from a
table such as:


+-----+------------+-------+-------+-------+
| id  |    Date    | Col_1 | Col_2 | Col_3 |
+=====+============+=======+=======+=======+
| 26  | 2012-10-18 |   X   |  25.7 | True  |
+-----+------------+-------+-------+-------+
| 42  | 2012-10-19 |   Y   | -12.4 | False |
+-----+------------+-------+-------+-------+
| 63  | 2012-10-20 |   Z   |  5.73 | True  |
+-----+------------+-------+-------+-------+


Functions from :mod:`pandas.io.sql` can extract some data into a DataFrame. In
the following example, we use the `SQlite <http://www.sqlite.org/>`_ SQL database
engine. You can use a temporary SQLite database where data are stored in
"memory". Just do:

.. code-block:: python

   import sqlite3
   from pandas.io import sql
   # Create your connection.
   cnx = sqlite3.connect(':memory:')

.. ipython:: python
   :suppress:

   import sqlite3
   from pandas.io import sql
   cnx = sqlite3.connect(':memory:')

.. ipython:: python
   :suppress:

   cu = cnx.cursor()
   # Create a table named 'data'.
   cu.execute("""CREATE TABLE data(id integer,
                                   date date,
                                   Col_1 string,
                                   Col_2 float,
                                   Col_3 bool);""")
   cu.executemany('INSERT INTO data VALUES (?,?,?,?,?)',
                  [(26, datetime.datetime(2010,10,18), 'X', 27.5, True),
                   (42, datetime.datetime(2010,10,19), 'Y', -12.5, False),
                   (63, datetime.datetime(2010,10,20), 'Z', 5.73, True)])


Let ``data`` be the name of your SQL table. With a query and your database
connection, just use the :func:`~pandas.io.sql.read_frame` function to get the
query results into a DataFrame:

.. ipython:: python

   sql.read_frame("SELECT * FROM data;", cnx)

You can also specify the name of the column as the DataFrame index:

.. ipython:: python

   sql.read_frame("SELECT * FROM data;", cnx, index_col='id')
   sql.read_frame("SELECT * FROM data;", cnx, index_col='date')

Of course, you can specify a more "complex" query.

.. ipython:: python

   sql.read_frame("SELECT id, Col_1, Col_2 FROM data WHERE id = 42;", cnx)

.. ipython:: python
   :suppress:

   cu.close()
   cnx.close()


There are a few other available functions:

  - ``tquery`` returns a list of tuples corresponding to each row.
  - ``uquery`` does the same thing as tquery, but instead of returning results
    it returns the number of related rows.
  - ``write_frame`` writes records stored in a DataFrame into the SQL table.
  - ``has_table`` checks if a given SQLite table exists.

.. note::

   For now, writing your DataFrame into a database works only with
   **SQLite**. Moreover, the **index** will currently be **dropped**.


STATA Format
------------

.. _io.stata:

Writing to STATA format
~~~~~~~~~~~~~~~~~~~~~~~

.. _io.stata_writer:

The method :func:`~pandas.core.frame.DataFrame.to_stata` will write a DataFrame
into a .dta file. The format version of this file is always the latest one, 115.

.. ipython:: python

   df = DataFrame(randn(10, 2), columns=list('AB'))
   df.to_stata('stata.dta')

Reading from STATA format
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _io.stata_reader:

.. versionadded:: 0.12

The top-level function ``read_stata`` will read a dta format file
and return a DataFrame:
The class :class:`~pandas.io.stata.StataReader` will read the header of the
given dta file at initialization. Its method
:func:`~pandas.io.stata.StataReader.data` will read the observations,
converting them to a DataFrame which is returned:

.. ipython:: python

   pd.read_stata('stata.dta')

Currently the ``index`` is retrieved as a column on read back.

The parameter ``convert_categoricals`` indicates wheter value labels should be
read and used to create a ``Categorical`` variable from them. Value labels can
also be retrieved by the function ``variable_labels``, which requires data to be
called before (see ``pandas.io.stata.StataReader``).

The StataReader supports .dta Formats 104, 105, 108, 113-115.
Alternatively, the function :func:`~pandas.io.stata.read_stata` can be used

.. ipython:: python
   :suppress:

   import os
   os.remove('stata.dta')

Data Reader
-----------

.. _io.data_reader:

Functions from :mod:`pandas.io.data` extract data from various Internet
sources into a DataFrame. Currently the following sources are supported:

    - Yahoo! Finance
    - Google Finance
    - St. Louis FED (FRED)
    - Kenneth French's data library

It should be noted, that various sources support different kinds of data, so not all sources implement the same methods and the data elements returned might also differ.

Yahoo! Finance
~~~~~~~~~~~~~~

.. ipython:: python

    import pandas.io.data as web
    start = datetime.datetime(2010, 1, 1)
    end = datetime.datetime(2013, 01, 27)
    f=web.DataReader("F", 'yahoo', start, end)
    f.ix['2010-01-04']

Google Finance
~~~~~~~~~~~~~~

.. ipython:: python

    import pandas.io.data as web
    start = datetime.datetime(2010, 1, 1)
    end = datetime.datetime(2013, 01, 27)
    f=web.DataReader("F", 'google', start, end)
    f.ix['2010-01-04']

FRED
~~~~

.. ipython:: python

    import pandas.io.data as web
    start = datetime.datetime(2010, 1, 1)
    end = datetime.datetime(2013, 01, 27)
    gdp=web.DataReader("GDP", "fred", start, end)
    gdp.ix['2013-01-01']


Fama/French
~~~~~~~~~~~

Tthe dataset names are listed at `Fama/French Data Library
<http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>`_)

.. ipython:: python

    import pandas.io.data as web
    ip=web.DataReader("5_Industry_Portfolios", "famafrench")
    ip[4].ix[192607]
