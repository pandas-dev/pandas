.. _io:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import os
   import csv
   from StringIO import StringIO
   import pandas as pd

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

.. _io.read_csv_table:

CSV & Text files
----------------

The two workhorse functions for reading text files (a.k.a. flat files) are
:func:`~pandas.io.parsers.read_csv` and :func:`~pandas.io.parsers.read_table`.
They both use the same parsing code to intelligently convert tabular
data into a DataFrame object. They can take a number of arguments:

See some :ref:`cookbook examples <cookbook.csv>` for some advanced strategies
See some :ref:`cookbook examples <cookbook.csv>` for some advanced strategies

  - ``filepath_or_buffer``: Either a string path to a file, or any object with a
    ``read`` method (such as an open file or ``StringIO``).
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
    pass ``header=0`` to be able to replace existing names.
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
  - ``comment``: denotes the start of a comment and ignores the rest of the line.
    Currently line commenting is not supported.
  - ``nrows``: Number of rows to read out of the file. Useful to only read a
    small portion of a large file
  - ``iterator``: If True, return a ``TextParser`` to enable reading a file
    into memory piece by piece
  - ``chunksize``: An number of rows to be used to "chunk" a file into
    pieces. Will cause an ``TextParser`` object to be returned. More on this
    below in the section on :ref:`iterating and chunking <io.chunking>`
  - ``skip_footer``: number of lines to skip at bottom of file (default 0)
  - ``converters``: a dictionary of functions for converting values in certain
    columns, where keys are either integers or column labels
  - ``encoding``: a string representing the encoding to use for decoding
    unicode data, e.g. ``'utf-8``` or ``'latin-1'``.
  - ``verbose``: show number of NA values inserted in non-numeric columns
  - ``squeeze``: if True then output with only one column is turned into Series
  - ``error_bad_lines``: if False then any lines causing an error will be skipped :ref:`bad lines <io.bad_lines>`

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
While `read_csv` reads delimited data, the :func:`~pandas.io.parsers.read_fwf`
function works with data files that have known and fixed column widths.
The function parameters to `read_fwf` are largely the same as `read_csv` with
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


Reading DataFrame objects with ``MultiIndex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _io.csv_multiindex:

Suppose you have data indexed by two columns:

.. ipython:: python

   print open('data/mindex_ex.csv').read()

The ``index_col`` argument to ``read_csv`` and ``read_table`` can take a list of
column numbers to turn multiple columns into a ``MultiIndex``:

.. ipython:: python

   df = pd.read_csv("data/mindex_ex.csv", index_col=[0,1])
   df
   df.ix[1978]

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
value will be an iterable object of type ``TextParser``:

.. ipython:: python

   reader = pd.read_table('tmp.sv', sep='|', chunksize=4)
   reader

   for chunk in reader:
       print chunk


Specifying ``iterator=True`` will also return the ``TextParser`` object:

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
  - ``nanRep``: A string representation of a missing value (default '')
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


Writing to HTML format
~~~~~~~~~~~~~~~~~~~~~~

.. _io.html:

DataFrame object has an instance method ``to_html`` which renders the contents
of the DataFrame as an html table. The function arguments are as in the method
``to_string`` described above.

Clipboard
---------

.. _io.clipboard:

A handy way to grab data is to use the ``read_clipboard`` method, which takes
the contents of the clipboard buffer and passes them to the ``read_table``
method described in the next section. For instance, you can copy the following
text to the clipboard (CTRL-C on many operating systems):

.. code-block:: python

     A B C
   x 1 4 p
   y 2 5 q
   z 3 6 r

And then import the data directly to a DataFrame by calling:

.. code-block:: python

   clipdf = pd.read_clipboard(delim_whitespace=True)

.. ipython:: python

   clipdf


.. _io.excel:

Excel files
-----------

The ``ExcelFile`` class can read an Excel 2003 file using the ``xlrd`` Python
module and use the same parsing code as the above to convert tabular data into
a DataFrame. To use it, create the ``ExcelFile`` object:

.. code-block:: python

   xls = ExcelFile('path_to_file.xls')

Then use the ``parse`` instance method with a sheetname, then use the same
additional arguments as the parsers above:

.. code-block:: python

   xls.parse('Sheet1', index_col=None, na_values=['NA'])

To read sheets from an Excel 2007 file, you can pass a filename with a ``.xlsx``
extension, in which case the ``openpyxl`` module will be used to read the file.

It is often the case that users will insert columns to do temporary computations
in Excel and you may not want to read in those columns. `ExcelFile.parse` takes
a `parse_cols` keyword to allow you to specify a subset of columns to parse.

If `parse_cols` is an integer, then it is assumed to indicate the last column
to be parsed.

.. code-block:: python

   xls.parse('Sheet1', parse_cols=2, index_col=None, na_values=['NA'])

If `parse_cols` is a list of integers, then it is assumed to be the file column
indices to be parsed.

.. code-block:: python

   xls.parse('Sheet1', parse_cols=[0, 2, 3], index_col=None, na_values=['NA'])

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
<http://www.pytables.org/>`__ library.

See some :ref:`cookbook examples <cookbook.hdf>` for some advanced strategies

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
Closing a Store

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


These stores are **not** appendable once written (though you can simply
remove them and rewrite). Nor are they **queryable**; they must be
retrieved in their entirety.

.. _io.hdf5-table:

Storing in Table format
~~~~~~~~~~~~~~~~~~~~~~~

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

The ``columns`` keyword can be supplied to select to filter a list of
the return columns, this is equivalent to passing a
``Term('columns', list_of_columns_to_filter)``

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

   # we have automagically created this index and that the B/C/string/string2
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

Note, that the chunksize keyword applies to the **returned** rows. So if you
are doing a query, then that set will be subdivided and returned in the
iterator. Keep in mind that if you do not pass a ``where`` selection criteria
then the ``nrows`` of the table are considered.

Advanced Queries
~~~~~~~~~~~~~~~~

**Unique**

To retrieve the *unique* values of an indexable or data column, use the
method ``unique``. This will, for example, enable you to get the index
very quickly. Note ``nan`` are excluded from the result set.

.. ipython:: python

   store.unique('df_dc', 'index')
   store.unique('df_dc', 'string')

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

New in 0.10.1 are the methods ``append_to_multple`` and
``select_as_multiple``, that can perform appending/selecting from
multiple tables at once. The idea is to have one table (call it the
selector table) that you index most/all of the columns, and perform your
queries. The other table(s) are data tables that are indexed the same as
the selector table. You can then perform a very fast query on the
selector table, yet get lots of data back. This method works similar to
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
objects (``Panel`` and ``Panel4D``). To get optimal deletion speed, it
pays to have the dimension you are deleting be the first of the
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

``PyTables`` offer better write performance when compressed after
writing them, as opposed to turning on compression at the very
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
   - You can not append/select/delete to a non-table (table creation is
     determined on the first append, or by passing ``table=True`` in a
     put operation)
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
     ``min_itemsize`` on the first table creation (``min_itemsize`` can
     be an integer or a dict of column name to an integer). If
     subsequent appends introduce elements in the indexing axis that are
     larger than the supported indexer, an Exception will be raised
     (otherwise you could have a silent truncation of these indexers,
     leading to loss of information). Just to be clear, this fixed-width
     restriction applies to **indexables** (the indexing columns) and
     **string values** in a mixed_type table.

     .. ipython:: python

       store.append('wp_big_strings', wp, min_itemsize = { 'minor_axis' : 30 })
       wp = wp.rename_axis(lambda x: x + '_big_strings', axis=2)
       store.append('wp_big_strings', wp)
       store.select('wp_big_strings')

       # we have provided a minimum minor_axis indexable size
       store.root.wp_big_strings.table

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

0.10.1 of ``HDFStore`` is backwards compatible for reading tables
created in a prior version of pandas however, query terms using the
prior (undocumented) methodology are unsupported. ``HDFStore`` will
issue a warning if you try to use a prior-version format file. You must
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
   - You can pass ``chunksize=an integer`` to ``append``, to change the
     writing chunksize (default is 50000). This will signficantly lower
     your memory usage on writing.
   - You can pass ``expectedrows=an integer`` to the first ``append``,
     to set the TOTAL number of expectedrows that ``PyTables`` will
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
store your data. Pass the ``axes`` keyword with a list of dimension
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
facilitate data retrieval and to reduce dependency on DB-specific API. There
wrappers only support the Python database adapters which respect the `Python
DB-API <http://www.python.org/dev/peps/pep-0249/>`_.

See some :ref:`cookbook examples <cookbook.sql>` for some advanced strategies

Suppose you want to query some data with different types from a table such as:

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
the following example, we use `SQlite <http://www.sqlite.org/>`_ SQL database
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

Of course, you can specify more "complex" query.

.. ipython:: python

   sql.read_frame("SELECT id, Col_1, Col_2 FROM data WHERE id = 42;", cnx)

.. ipython:: python
   :suppress:

   cu.close()
   cnx.close()


There are a few other available functions:

  - ``tquery`` returns list of tuples corresponding to each row.
  - ``uquery`` does the same thing as tquery, but instead of returning results,
    it returns the number of related rows.
  - ``write_frame`` writes records stored in a DataFrame into the SQL table.
  - ``has_table`` checks if a given SQLite table exists.

.. note::

   For now, writing your DataFrame into a database works only with
   **SQLite**. Moreover, the **index** will currently be **dropped**.
