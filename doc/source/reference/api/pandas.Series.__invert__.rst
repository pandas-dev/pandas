.. currentmodule:: pandas

pandas.Series.__invert__
========================

Elementwise invert (``~``) for :class:`pandas.Series`.

**Signature**

``Series.__invert__(self) -> Series``

**Summary**

For boolean and nullable-boolean dtypes, ``~s`` toggles the mask
(``True`` ↔ ``False``) and propagates ``pd.NA``. For integer dtypes,
``~`` performs bitwise invert as in Python.

.. seealso::
   :ref:`indexing.boolean`

**Examples**

Boolean / nullable boolean::

    >>> s = pd.Series([True, pd.NA, False], dtype="boolean")
    >>> ~s
    0    False
    1     <NA>
    2     True
    dtype: boolean

Integer vs boolean::

    >>> s_int = pd.Series([0, 1, 2], dtype="int64")
    >>> ~s_int
    0    -1
    1    -2
    2    -3
    dtype: int64

Arrow-backed boolean (if pyarrow installed)::

    >>> s_arrow = pd.Series([True, pd.NA, False], dtype="boolean[pyarrow]")
    >>> ~s_arrow
    0    False
    1     <NA>
    2     True
    dtype: boolean[pyarrow]

**Notes**

- In Python’s stdlib, :func:`operator.__invert__` is bitwise invert on integers.
  In pandas, ``~`` on boolean arrays is elementwise invert.
