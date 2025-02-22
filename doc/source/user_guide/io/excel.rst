.. _io.excel:

===========
Excel files
===========

The :func:`~pandas.read_excel` method can read Excel 2007+ (``.xlsx``) files
using the ``openpyxl`` Python module. Excel 2003 (``.xls``) files
can be read using ``xlrd``. Binary Excel (``.xlsb``)
files can be read using ``pyxlsb``. All formats can be read
using :ref:`calamine<io.calamine>` engine.
The :meth:`~DataFrame.to_excel` instance method is used for
saving a ``DataFrame`` to Excel.  Generally the semantics are
similar to working with :ref:`csv<io.read_csv_table>` data.
See the :ref:`cookbook<cookbook.excel>` for some advanced strategies.

.. note::

   When ``engine=None``, the following logic will be used to determine the engine:

   - If ``path_or_buffer`` is an OpenDocument format (.odf, .ods, .odt),
     then `odf <https://pypi.org/project/odfpy/>`_ will be used.
   - Otherwise if ``path_or_buffer`` is an xls format, ``xlrd`` will be used.
   - Otherwise if ``path_or_buffer`` is in xlsb format, ``pyxlsb`` will be used.
   - Otherwise ``openpyxl`` will be used.

.. _io.excel_reader:

Reading Excel files
'''''''''''''''''''

In the most basic use-case, ``read_excel`` takes a path to an Excel
file, and the ``sheet_name`` indicating which sheet to parse.

When using the ``engine_kwargs`` parameter, pandas will pass these arguments to the
engine. For this, it is important to know which function pandas is
using internally.

* For the engine openpyxl, pandas is using :func:`openpyxl.load_workbook` to read in (``.xlsx``) and (``.xlsm``) files.

* For the engine xlrd, pandas is using :func:`xlrd.open_workbook` to read in (``.xls``) files.

* For the engine pyxlsb, pandas is using :func:`pyxlsb.open_workbook` to read in (``.xlsb``) files.

* For the engine odf, pandas is using :func:`odf.opendocument.load` to read in (``.ods``) files.

* For the engine calamine, pandas is using :func:`python_calamine.load_workbook`
  to read in (``.xlsx``), (``.xlsm``), (``.xls``), (``.xlsb``), (``.ods``) files.

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.xls", sheet_name="Sheet1")


.. _io.excel.excelfile_class:

``ExcelFile`` class
+++++++++++++++++++

To facilitate working with multiple sheets from the same file, the ``ExcelFile``
class can be used to wrap the file and can be passed into ``read_excel``
There will be a performance benefit for reading multiple sheets as the file is
read into memory only once.

.. code-block:: python

   xlsx = pd.ExcelFile("path_to_file.xls")
   df = pd.read_excel(xlsx, "Sheet1")

The ``ExcelFile`` class can also be used as a context manager.

.. code-block:: python

   with pd.ExcelFile("path_to_file.xls") as xls:
       df1 = pd.read_excel(xls, "Sheet1")
       df2 = pd.read_excel(xls, "Sheet2")

The ``sheet_names`` property will generate
a list of the sheet names in the file.

The primary use-case for an ``ExcelFile`` is parsing multiple sheets with
different parameters:

.. code-block:: python

    data = {}
    # For when Sheet1's format differs from Sheet2
    with pd.ExcelFile("path_to_file.xls") as xls:
        data["Sheet1"] = pd.read_excel(xls, "Sheet1", index_col=None, na_values=["NA"])
        data["Sheet2"] = pd.read_excel(xls, "Sheet2", index_col=1)

Note that if the same parsing parameters are used for all sheets, a list
of sheet names can simply be passed to ``read_excel`` with no loss in performance.

.. code-block:: python

    # using the ExcelFile class
    data = {}
    with pd.ExcelFile("path_to_file.xls") as xls:
        data["Sheet1"] = pd.read_excel(xls, "Sheet1", index_col=None, na_values=["NA"])
        data["Sheet2"] = pd.read_excel(xls, "Sheet2", index_col=None, na_values=["NA"])

    # equivalent using the read_excel function
    data = pd.read_excel(
        "path_to_file.xls", ["Sheet1", "Sheet2"], index_col=None, na_values=["NA"]
    )

``ExcelFile`` can also be called with a ``xlrd.book.Book`` object
as a parameter. This allows the user to control how the excel file is read.
For example, sheets can be loaded on demand by calling ``xlrd.open_workbook()``
with ``on_demand=True``.

.. code-block:: python

    import xlrd

    xlrd_book = xlrd.open_workbook("path_to_file.xls", on_demand=True)
    with pd.ExcelFile(xlrd_book) as xls:
        df1 = pd.read_excel(xls, "Sheet1")
        df2 = pd.read_excel(xls, "Sheet2")

.. _io.excel.specifying_sheets:

Specifying sheets
+++++++++++++++++

.. note:: The second argument is ``sheet_name``, not to be confused with ``ExcelFile.sheet_names``.

.. note:: An ExcelFile's attribute ``sheet_names`` provides access to a list of sheets.

* The arguments ``sheet_name`` allows specifying the sheet or sheets to read.
* The default value for ``sheet_name`` is 0, indicating to read the first sheet
* Pass a string to refer to the name of a particular sheet in the workbook.
* Pass an integer to refer to the index of a sheet. Indices follow Python
  convention, beginning at 0.
* Pass a list of either strings or integers, to return a dictionary of specified sheets.
* Pass a ``None`` to return a dictionary of all available sheets.

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.xls", "Sheet1", index_col=None, na_values=["NA"])

Using the sheet index:

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.xls", 0, index_col=None, na_values=["NA"])

Using all default values:

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.xls")

Using None to get all sheets:

.. code-block:: python

   # Returns a dictionary of DataFrames
   pd.read_excel("path_to_file.xls", sheet_name=None)

Using a list to get multiple sheets:

.. code-block:: python

   # Returns the 1st and 4th sheet, as a dictionary of DataFrames.
   pd.read_excel("path_to_file.xls", sheet_name=["Sheet1", 3])

``read_excel`` can read more than one sheet, by setting ``sheet_name`` to either
a list of sheet names, a list of sheet positions, or ``None`` to read all sheets.
Sheets can be specified by sheet index or sheet name, using an integer or string,
respectively.

.. _io.excel.reading_multiindex:

Reading a ``MultiIndex``
++++++++++++++++++++++++

``read_excel`` can read a ``MultiIndex`` index, by passing a list of columns to ``index_col``
and a ``MultiIndex`` column by passing a list of rows to ``header``.  If either the ``index``
or ``columns`` have serialized level names those will be read in as well by specifying
the rows/columns that make up the levels.

For example, to read in a ``MultiIndex`` index without names:

.. ipython:: python

   df = pd.DataFrame(
       {"a": [1, 2, 3, 4], "b": [5, 6, 7, 8]},
       index=pd.MultiIndex.from_product([["a", "b"], ["c", "d"]]),
   )
   df.to_excel("path_to_file.xlsx")
   df = pd.read_excel("path_to_file.xlsx", index_col=[0, 1])
   df

If the index has level names, they will be parsed as well, using the same
parameters.

.. ipython:: python

   df.index = df.index.set_names(["lvl1", "lvl2"])
   df.to_excel("path_to_file.xlsx")
   df = pd.read_excel("path_to_file.xlsx", index_col=[0, 1])
   df


If the source file has both ``MultiIndex`` index and columns, lists specifying each
should be passed to ``index_col`` and ``header``:

.. ipython:: python

   df.columns = pd.MultiIndex.from_product([["a"], ["b", "d"]], names=["c1", "c2"])
   df.to_excel("path_to_file.xlsx")
   df = pd.read_excel("path_to_file.xlsx", index_col=[0, 1], header=[0, 1])
   df

.. ipython:: python
   :suppress:

   import os
   os.remove("path_to_file.xlsx")

Missing values in columns specified in ``index_col`` will be forward filled to
allow roundtripping with ``to_excel`` for ``merged_cells=True``. To avoid forward
filling the missing values use ``set_index`` after reading the data instead of
``index_col``.

Parsing specific columns
++++++++++++++++++++++++

It is often the case that users will insert columns to do temporary computations
in Excel and you may not want to read in those columns. ``read_excel`` takes
a ``usecols`` keyword to allow you to specify a subset of columns to parse.

You can specify a comma-delimited set of Excel columns and ranges as a string:

.. code-block:: python

   pd.read_excel("path_to_file.xls", "Sheet1", usecols="A,C:E")

If ``usecols`` is a list of integers, then it is assumed to be the file column
indices to be parsed.

.. code-block:: python

   pd.read_excel("path_to_file.xls", "Sheet1", usecols=[0, 2, 3])

Element order is ignored, so ``usecols=[0, 1]`` is the same as ``[1, 0]``.

If ``usecols`` is a list of strings, it is assumed that each string corresponds
to a column name provided either by the user in ``names`` or inferred from the
document header row(s). Those strings define which columns will be parsed:

.. code-block:: python

    pd.read_excel("path_to_file.xls", "Sheet1", usecols=["foo", "bar"])

Element order is ignored, so ``usecols=['baz', 'joe']`` is the same as ``['joe', 'baz']``.

If ``usecols`` is callable, the callable function will be evaluated against
the column names, returning names where the callable function evaluates to ``True``.

.. code-block:: python

    pd.read_excel("path_to_file.xls", "Sheet1", usecols=lambda x: x.isalpha())

Parsing dates
+++++++++++++

Datetime-like values are normally automatically converted to the appropriate
dtype when reading the excel file. But if you have a column of strings that
*look* like dates (but are not actually formatted as dates in excel), you can
use the ``parse_dates`` keyword to parse those strings to datetimes:

.. code-block:: python

   pd.read_excel("path_to_file.xls", "Sheet1", parse_dates=["date_strings"])


Cell converters
+++++++++++++++

It is possible to transform the contents of Excel cells via the ``converters``
option. For instance, to convert a column to boolean:

.. code-block:: python

   pd.read_excel("path_to_file.xls", "Sheet1", converters={"MyBools": bool})

This options handles missing values and treats exceptions in the converters
as missing data. Transformations are applied cell by cell rather than to the
column as a whole, so the array dtype is not guaranteed. For instance, a
column of integers with missing values cannot be transformed to an array
with integer dtype, because NaN is strictly a float. You can manually mask
missing data to recover integer dtype:

.. code-block:: python

   def cfun(x):
       return int(x) if x else -1


   pd.read_excel("path_to_file.xls", "Sheet1", converters={"MyInts": cfun})

Dtype specifications
++++++++++++++++++++

As an alternative to converters, the type for an entire column can
be specified using the ``dtype`` keyword, which takes a dictionary
mapping column names to types.  To interpret data with
no type inference, use the type ``str`` or ``object``.

.. code-block:: python

   pd.read_excel("path_to_file.xls", dtype={"MyInts": "int64", "MyText": str})

.. _io.excel_writer:

Writing Excel files
'''''''''''''''''''

Writing Excel files to disk
+++++++++++++++++++++++++++

To write a ``DataFrame`` object to a sheet of an Excel file, you can use the
``to_excel`` instance method.  The arguments are largely the same as ``to_csv``
described above, the first argument being the name of the excel file, and the
optional second argument the name of the sheet to which the ``DataFrame`` should be
written. For example:

.. code-block:: python

   df.to_excel("path_to_file.xlsx", sheet_name="Sheet1")

Files with a
``.xlsx`` extension will be written using ``xlsxwriter`` (if available) or
``openpyxl``.

The ``DataFrame`` will be written in a way that tries to mimic the REPL output.
The ``index_label`` will be placed in the second
row instead of the first. You can place it in the first row by setting the
``merge_cells`` option in ``to_excel()`` to ``False``:

.. code-block:: python

   df.to_excel("path_to_file.xlsx", index_label="label", merge_cells=False)

In order to write separate ``DataFrames`` to separate sheets in a single Excel file,
one can pass an :class:`~pandas.io.excel.ExcelWriter`.

.. code-block:: python

   with pd.ExcelWriter("path_to_file.xlsx") as writer:
       df1.to_excel(writer, sheet_name="Sheet1")
       df2.to_excel(writer, sheet_name="Sheet2")

.. _io.excel_writing_buffer:

When using the ``engine_kwargs`` parameter, pandas will pass these arguments to the
engine. For this, it is important to know which function pandas is using internally.

* For the engine openpyxl, pandas is using :func:`openpyxl.Workbook` to create a new sheet and :func:`openpyxl.load_workbook` to append data to an existing sheet. The openpyxl engine writes to (``.xlsx``) and (``.xlsm``) files.

* For the engine xlsxwriter, pandas is using :func:`xlsxwriter.Workbook` to write to (``.xlsx``) files.

* For the engine odf, pandas is using :func:`odf.opendocument.OpenDocumentSpreadsheet` to write to (``.ods``) files.

Writing Excel files to memory
+++++++++++++++++++++++++++++

pandas supports writing Excel files to buffer-like objects such as ``StringIO`` or
``BytesIO`` using :class:`~pandas.io.excel.ExcelWriter`.

.. code-block:: python

   from io import BytesIO

   bio = BytesIO()

   # By setting the 'engine' in the ExcelWriter constructor.
   writer = pd.ExcelWriter(bio, engine="xlsxwriter")
   df.to_excel(writer, sheet_name="Sheet1")

   # Save the workbook
   writer.save()

   # Seek to the beginning and read to copy the workbook to a variable in memory
   bio.seek(0)
   workbook = bio.read()

.. note::

    ``engine`` is optional but recommended.  Setting the engine determines
    the version of workbook produced. Setting ``engine='xlrd'`` will produce an
    Excel 2003-format workbook (xls).  Using either ``'openpyxl'`` or
    ``'xlsxwriter'`` will produce an Excel 2007-format workbook (xlsx). If
    omitted, an Excel 2007-formatted workbook is produced.


.. _io.excel.writers:

Excel writer engines
''''''''''''''''''''

pandas chooses an Excel writer via two methods:

1. the ``engine`` keyword argument
2. the filename extension (via the default specified in config options)

By default, pandas uses the `XlsxWriter`_  for ``.xlsx``, `openpyxl`_
for ``.xlsm``. If you have multiple
engines installed, you can set the default engine through :ref:`setting the
config options <options>` ``io.excel.xlsx.writer`` and
``io.excel.xls.writer``. pandas will fall back on `openpyxl`_ for ``.xlsx``
files if `Xlsxwriter`_ is not available.

.. _XlsxWriter: https://xlsxwriter.readthedocs.io
.. _openpyxl: https://openpyxl.readthedocs.io/

To specify which writer you want to use, you can pass an engine keyword
argument to ``to_excel`` and to ``ExcelWriter``. The built-in engines are:

* ``openpyxl``: version 2.4 or higher is required
* ``xlsxwriter``

.. code-block:: python

   # By setting the 'engine' in the DataFrame 'to_excel()' methods.
   df.to_excel("path_to_file.xlsx", sheet_name="Sheet1", engine="xlsxwriter")

   # By setting the 'engine' in the ExcelWriter constructor.
   writer = pd.ExcelWriter("path_to_file.xlsx", engine="xlsxwriter")

   # Or via pandas configuration.
   from pandas import options  # noqa: E402

   options.io.excel.xlsx.writer = "xlsxwriter"

   df.to_excel("path_to_file.xlsx", sheet_name="Sheet1")

.. _io.excel.style:

Style and formatting
''''''''''''''''''''

The look and feel of Excel worksheets created from pandas can be modified using the following parameters on the ``DataFrame``'s ``to_excel`` method.

* ``float_format`` : Format string for floating point numbers (default ``None``).
* ``freeze_panes`` : A tuple of two integers representing the bottommost row and rightmost column to freeze. Each of these parameters is one-based, so (1, 1) will freeze the first row and first column (default ``None``).

.. note::

    As of pandas 3.0, by default spreadsheets created with the ``to_excel`` method
    will not contain any styling. Users wishing to bold text, add bordered styles,
    etc in a worksheet output by ``to_excel`` can do so by using :meth:`Styler.to_excel`
    to create styled excel files. For documentation on styling spreadsheets, see
    `here <https://pandas.pydata.org/docs/user_guide/style.html#Export-to-Excel>`__.


.. code-block:: python

    css = "border: 1px solid black; font-weight: bold;"
    df.style.map_index(lambda x: css).map_index(lambda x: css, axis=1).to_excel("myfile.xlsx")

Using the `Xlsxwriter`_ engine provides many options for controlling the
format of an Excel worksheet created with the ``to_excel`` method.  Excellent examples can be found in the
`Xlsxwriter`_ documentation here: https://xlsxwriter.readthedocs.io/working_with_pandas.html

.. _io.ods:

OpenDocument Spreadsheets
'''''''''''''''''''''''''

The io methods for `Excel files`_ also support reading and writing OpenDocument spreadsheets
using the `odfpy <https://pypi.org/project/odfpy/>`__ module. The semantics and features for reading and writing
OpenDocument spreadsheets match what can be done for `Excel files`_ using
``engine='odf'``. The optional dependency 'odfpy' needs to be installed.

The :func:`~pandas.read_excel` method can read OpenDocument spreadsheets

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.ods", engine="odf")

Similarly, the :func:`~pandas.to_excel` method can write OpenDocument spreadsheets

.. code-block:: python

   # Writes DataFrame to a .ods file
   df.to_excel("path_to_file.ods", engine="odf")

.. _io.xlsb:

Binary Excel (.xlsb) files
''''''''''''''''''''''''''

The :func:`~pandas.read_excel` method can also read binary Excel files
using the ``pyxlsb`` module. The semantics and features for reading
binary Excel files mostly match what can be done for `Excel files`_ using
``engine='pyxlsb'``. ``pyxlsb`` does not recognize datetime types
in files and will return floats instead (you can use :ref:`calamine<io.calamine>`
if you need recognize datetime types).

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.xlsb", engine="pyxlsb")

.. note::

   Currently pandas only supports *reading* binary Excel files. Writing
   is not implemented.

.. _io.calamine:

Calamine (Excel and ODS files)
''''''''''''''''''''''''''''''

The :func:`~pandas.read_excel` method can read Excel file (``.xlsx``, ``.xlsm``, ``.xls``, ``.xlsb``)
and OpenDocument spreadsheets (``.ods``) using the ``python-calamine`` module.
This module is a binding for Rust library `calamine <https://crates.io/crates/calamine>`__
and is faster than other engines in most cases. The optional dependency 'python-calamine' needs to be installed.

.. code-block:: python

   # Returns a DataFrame
   pd.read_excel("path_to_file.xlsb", engine="calamine")
