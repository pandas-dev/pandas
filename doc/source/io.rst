.. _io:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import os
   np.random.seed(123456)
   from pandas import *
   from StringIO import StringIO
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')
   clipdf = DataFrame({'A':[1,2,3],'B':[4,5,6],'C':['p','q','r']},
                      index=['x','y','z'])

*******************************
IO Tools (Text, CSV, HDF5, ...)
*******************************

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

   clipdf = read_clipboard(sep='\s*')

.. ipython:: python

   clipdf

CSV & Text files
----------------

.. _io.parse_dates:

The two workhorse functions for reading text files (a.k.a. flat files) are
:func:`~pandas.io.parsers.read_csv` and :func:`~pandas.io.parsers.read_table`.
They both use the same parsing code to intelligently convert tabular
data into a DataFrame object. They can take a number of arguments:

  - ``path_or_buffer``: Either a string path to a file, or any object with a
    ``read`` method (such as an open file or ``StringIO``).
  - ``sep`` or ``delimiter``: A delimiter / separator to split fields
    on. `read_csv` is capable of inferring the delimiter automatically in some
    cases by "sniffing." The separator may be specified as a regular
    expression; for instance you may use '\s*' to indicate arbitrary
    whitespace.
  - ``header``: row number to use as the column names, and the start of the data.
    Defaults to 0 (first row); specify None if there is no header row.
  - ``names``: List of column names to use. If passed, header will be
    implicitly set to None.
  - ``skiprows``: A collection of numbers for rows in the file to skip. Can
    also be an integer to skip the first ``n`` rows
  - ``index_col``: column number, or list of column numbers, to use as the
    ``index`` (row labels) of the resulting DataFrame. By default, it will number
    the rows without using any column, unless there is one more data column than
    there are headers, in which case the first column is taken as the index.
  - ``parse_dates``: If True, attempt to parse the index column as dates. False
    by default.
  - ``date_parser``: function to use to parse strings into datetime
    objects. If ``parse_dates`` is True, it defaults to the very robust
    ``dateutil.parser``. Specifying this implicitly sets ``parse_dates`` as True.
  - ``na_values``: optional list of strings to recognize as NaN (missing values),
    in addition to a default set.
  - ``nrows``: Number of rows to read out of the file. Useful to only read a
    small portion of a large file
  - ``chunksize``: An number of rows to be used to "chunk" a file into
    pieces. Will cause an ``TextParser`` object to be returned. More on this
    below in the section on :ref:`iterating and chunking <io.chunking>`
  - ``iterator``: If True, return a ``TextParser`` to enable reading a file
    into memory piece by piece
  - ``skip_footer``: number of lines to skip at bottom of file (default 0)
  - ``converters``: a dictionary of functions for converting values in certain
    columns, where keys are either integers or column labels
  - ``encoding``: a string representing the encoding to use if the contents are
    non-ascii, for python versions prior to 3 
  - ``verbose`` : show number of NA values inserted in non-numeric columns

.. ipython:: python
   :suppress:

   f = open('foo.csv', 'w')
   f.write('date,A,B,C\n20090101,a,1,2\n20090102,b,3,4\n20090103,c,4,5')
   f.close()

Consider a typical CSV file containing, in this case, some time series data:

.. ipython:: python

   print open('foo.csv').read()

The default for `read_csv` is to create a DataFrame with simple numbered rows:

.. ipython:: python

   read_csv('foo.csv')

In the case of indexed data, you can pass the column number (or a list of
column numbers, for a hierarchical index) you wish to use as the index. If the
index values are dates and you want them to be converted to ``datetime``
objects, pass ``parse_dates=True``:

.. ipython:: python

   # Use a column as an index, and parse it as dates.
   df = read_csv('foo.csv', index_col=0, parse_dates=True)
   df
   # These are python datetime objects
   df.index

.. ipython:: python
   :suppress:

   os.remove('foo.csv')

The parsers make every attempt to "do the right thing" and not be very
fragile. Type inference is a pretty big deal. So if a column can be coerced to
integer dtype without altering the contents, it will do so. Any non-numeric
columns will come through as object dtype as with the rest of pandas objects.

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

   read_csv('foo.csv')

Note that the dates weren't automatically parsed. In that case you would need
to do as before:

.. ipython:: python

   df = read_csv('foo.csv', parse_dates=True)
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

   df = read_csv("data/mindex_ex.csv", index_col=[0,1])
   df
   df.ix[1978]

.. _io.sniff:

Automatically "sniffing" the delimiter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``read_csv`` is capable of inferring delimited (not necessarily
comma-separated) files. YMMV, as pandas uses the Sniffer_ class of the csv
module.

.. ipython:: python
   :suppress:

   df[:7].to_csv('tmp.sv', sep='|')
   df[:7].to_csv('tmp2.sv', sep=':')

.. ipython:: python

    print open('tmp2.sv').read()
    read_csv('tmp2.sv')

.. _Sniffer: http://docs.python.org/library/csv.html#csv.Sniffer

.. _io.chunking:

Iterating through files chunk by chunk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you wish to iterate through a (potentially very large) file lazily
rather than reading the entire file into memory, such as the following:


.. ipython:: python

   print open('tmp.sv').read()
   table = read_table('tmp.sv', sep='|')
   table


By specifiying a ``chunksize`` to ``read_csv`` or ``read_table``, the return
value will be an iterable object of type ``TextParser``:

.. ipython::

   In [1]: reader = read_table('tmp.sv', sep='|', chunksize=4)

   In [1]: reader

   In [2]: for chunk in reader:
      ...:     print chunk
      ...:

Specifying ``iterator=True`` will also return the ``TextParser`` object:

.. ipython:: python

   reader = read_table('tmp.sv', sep='|', iterator=True)
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
    ``nanRep``: A string representation of a missing value (default '')
  - ``cols``: Columns to write (default None)
  - ``header``: Whether to write out the column names (default True)
  - ``index``: whether to write row (index) names (default True)
  - ``index_label``: Column label(s) for index column(s) if desired. If None
    (default), and `header` and `index` are True, then the index names are
    used. (A sequence should be given if the DataFrame uses MultiIndex).
  - ``mode`` : Python write mode, default 'w'
  - ``sep`` : Field delimiter for the output file (default "'")
  - ``encoding``: a string representing the encoding to use if the contents are
    non-ascii, for python versions prior to 3 

Writing a formatted string
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _io.formatting:

The DataFrame object has an instance method ``to_string`` which allows control
over the string representation of the object. All arguments are optional:

  - ``buf`` default None, for example a StringIO object
  - ``columns`` default None, which columns to write
  - ``col_space`` default None, number of spaces to write between columns
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


Excel 2003 files
----------------

The ``ExcelFile`` class can read an Excel 2003 file using the ``xlrd`` Python
module and use the same parsing code as the above to convert tabular data into
a DataFrame. To use it, create the ``ExcelFile`` object:

.. code-block:: python

   xls = ExcelFile('path_to_file.xls')

Then use the ``parse`` instance method with a sheetname, then use the same
additional arguments as the parsers above:

.. code-block:: python

   xls.parse('Sheet1', index_col=None, na_values=['NA'])

HDF5 (PyTables)
---------------

``HDFStore`` is a dict-like object which reads and writes pandas to the high
performance HDF5 format using the excellent `PyTables
<http://www.pytables.org/>`__ library.

.. ipython:: python
   :suppress:
   :okexcept:

   os.remove('store.h5')

.. ipython:: python

   store = HDFStore('store.h5')
   print store

Objects can be written to the file just like adding key-value pairs to a dict:

.. ipython:: python

   index = DateRange('1/1/2000', periods=8)
   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   df = DataFrame(randn(8, 3), index=index,
                  columns=['A', 'B', 'C'])
   wp = Panel(randn(2, 5, 4), items=['Item1', 'Item2'],
              major_axis=DateRange('1/1/2000', periods=5),
              minor_axis=['A', 'B', 'C', 'D'])

   store['s'] = s
   store['df'] = df
   store['wp'] = wp
   store

In a current or later Python session, you can retrieve stored objects:

.. ipython:: python

   store['df']

.. ipython:: python
   :suppress:

   store.close()
   import os
   os.remove('store.h5')


.. Storing in Table format
.. ~~~~~~~~~~~~~~~~~~~~~~~

.. Querying objects stored in Table format
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
