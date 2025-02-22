.. _io.stata:

============
STATA format
============

.. _io.stata_writer:

Writing to stata format
'''''''''''''''''''''''

The method :func:`.DataFrame.to_stata` will write a DataFrame
into a .dta file. The format version of this file is always 115 (Stata 12).

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 2), columns=list("AB"))
   df.to_stata("stata.dta")

*Stata* data files have limited data type support; only strings with
244 or fewer characters, ``int8``, ``int16``, ``int32``, ``float32``
and ``float64`` can be stored in ``.dta`` files.  Additionally,
*Stata* reserves certain values to represent missing data. Exporting a
non-missing value that is outside of the permitted range in Stata for
a particular data type will retype the variable to the next larger
size.  For example, ``int8`` values are restricted to lie between -127
and 100 in Stata, and so variables with values above 100 will trigger
a conversion to ``int16``. ``nan`` values in floating points data
types are stored as the basic missing data type (``.`` in *Stata*).

.. note::

    It is not possible to export missing data values for integer data types.


The *Stata* writer gracefully handles other data types including ``int64``,
``bool``, ``uint8``, ``uint16``, ``uint32`` by casting to
the smallest supported type that can represent the data.  For example, data
with a type of ``uint8`` will be cast to ``int8`` if all values are less than
100 (the upper bound for non-missing ``int8`` data in *Stata*), or, if values are
outside of this range, the variable is cast to ``int16``.


.. warning::

   Conversion from ``int64`` to ``float64`` may result in a loss of precision
   if ``int64`` values are larger than 2**53.

.. warning::

  :class:`~pandas.io.stata.StataWriter` and
  :func:`.DataFrame.to_stata` only support fixed width
  strings containing up to 244 characters, a limitation imposed by the version
  115 dta file format. Attempting to write *Stata* dta files with strings
  longer than 244 characters raises a ``ValueError``.

.. _io.stata_reader:

Reading from Stata format
'''''''''''''''''''''''''

The top-level function ``read_stata`` will read a dta file and return
either a ``DataFrame`` or a :class:`pandas.api.typing.StataReader` that can
be used to read the file incrementally.

.. ipython:: python

   pd.read_stata("stata.dta")

Specifying a ``chunksize`` yields a
:class:`pandas.api.typing.StataReader` instance that can be used to
read ``chunksize`` lines from the file at a time.  The ``StataReader``
object can be used as an iterator.

.. ipython:: python

  with pd.read_stata("stata.dta", chunksize=3) as reader:
      for df in reader:
          print(df.shape)

For more fine-grained control, use ``iterator=True`` and specify
``chunksize`` with each call to
:func:`~pandas.io.stata.StataReader.read`.

.. ipython:: python

  with pd.read_stata("stata.dta", iterator=True) as reader:
      chunk1 = reader.read(5)
      chunk2 = reader.read(5)

Currently the ``index`` is retrieved as a column.

The parameter ``convert_categoricals`` indicates whether value labels should be
read and used to create a ``Categorical`` variable from them. Value labels can
also be retrieved by the function ``value_labels``, which requires :func:`~pandas.io.stata.StataReader.read`
to be called before use.

The parameter ``convert_missing`` indicates whether missing value
representations in Stata should be preserved.  If ``False`` (the default),
missing values are represented as ``np.nan``.  If ``True``, missing values are
represented using ``StataMissingValue`` objects, and columns containing missing
values will have ``object`` data type.

.. note::

   :func:`~pandas.read_stata` and
   :class:`~pandas.io.stata.StataReader` support .dta formats 113-115
   (Stata 10-12), 117 (Stata 13), and 118 (Stata 14).

.. note::

   Setting ``preserve_dtypes=False`` will upcast to the standard pandas data types:
   ``int64`` for all integer types and ``float64`` for floating point data.  By default,
   the Stata data types are preserved when importing.

.. note::

   All :class:`~pandas.io.stata.StataReader` objects, whether created by :func:`~pandas.read_stata`
   (when using ``iterator=True`` or ``chunksize``) or instantiated by hand, must be used as context
   managers (e.g. the ``with`` statement).
   While the :meth:`~pandas.io.stata.StataReader.close` method is available, its use is unsupported.
   It is not part of the public API and will be removed in with future without warning.

.. ipython:: python
   :suppress:

   import os
   os.remove("stata.dta")

.. _io.stata-categorical:

Categorical data
++++++++++++++++

``Categorical`` data can be exported to *Stata* data files as value labeled data.
The exported data consists of the underlying category codes as integer data values
and the categories as value labels.  *Stata* does not have an explicit equivalent
to a ``Categorical`` and information about *whether* the variable is ordered
is lost when exporting.

.. warning::

    *Stata* only supports string value labels, and so ``str`` is called on the
    categories when exporting data.  Exporting ``Categorical`` variables with
    non-string categories produces a warning, and can result a loss of
    information if the ``str`` representations of the categories are not unique.

Labeled data can similarly be imported from *Stata* data files as ``Categorical``
variables using the keyword argument ``convert_categoricals`` (``True`` by default).
The keyword argument ``order_categoricals`` (``True`` by default) determines
whether imported ``Categorical`` variables are ordered.

.. note::

    When importing categorical data, the values of the variables in the *Stata*
    data file are not preserved since ``Categorical`` variables always
    use integer data types between ``-1`` and ``n-1`` where ``n`` is the number
    of categories. If the original values in the *Stata* data file are required,
    these can be imported by setting ``convert_categoricals=False``, which will
    import original data (but not the variable labels). The original values can
    be matched to the imported categorical data since there is a simple mapping
    between the original *Stata* data values and the category codes of imported
    Categorical variables: missing values are assigned code ``-1``, and the
    smallest original value is assigned ``0``, the second smallest is assigned
    ``1`` and so on until the largest original value is assigned the code ``n-1``.

.. note::

    *Stata* supports partially labeled series. These series have value labels for
    some but not all data values. Importing a partially labeled series will produce
    a ``Categorical`` with string categories for the values that are labeled and
    numeric categories for values with no label.
