.. _sparse:

{{ header }}

**********************
Sparse data structures
**********************

.. note::

   ``SparseSeries`` and ``SparseDataFrame`` have been deprecated. Their purpose
   is served equally well by a :class:`Series` or :class:`DataFrame` with
   sparse values. See :ref:`sparse.migration` for tips on migrating.

Pandas provides data structures for efficiently storing sparse data.
These are not necessarily sparse in the typical "mostly 0". Rather, you can view these
objects as being "compressed" where any data matching a specific value (``NaN`` / missing value, though any value
can be chosen, including 0) is omitted. A special ``SparseIndex`` object tracks where data has been
"sparsified". For example,

.. ipython:: python

   arr = np.random.randn(10)
   arr[2:-2] = np.nan
   ts = pd.Series(pd.SparseArray(arr))
   ts

Notice the dtype, ``Sparse[float64, nan]``. The ``nan`` means that elements in the
array that are ``nan`` aren't actually stored, only the non-``nan`` elements are.
Those non-``nan`` elements have a ``float64`` dtype.

The sparse objects exist for memory efficiency reasons. Suppose you had a
large, mostly NA ``DataFrame``:

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10000, 4))
   df.iloc[:9998] = np.nan
   sdf = df.astype(pd.SparseDtype("float", np.nan))
   sdf
   sdf.sparse.density

As you can see, the density (% of values that have not been "compressed") is
extremely low. This sparse object takes up much less memory on disk (pickled)
and in the Python interpreter. Functionally, their behavior should be nearly
identical to their dense counterparts.

.. _sparse.array:

SparseArray
-----------

:class:`SparseArray` is a :class:`~pandas.api.extensions.ExtensionArray`
for storing an array of sparse values (see :ref:`basics.dtypes` for more
on extension arrays). It is a 1-dimensional ndarray-like object storing
only values distinct from the ``fill_value``:

.. ipython:: python

   arr = np.random.randn(10)
   arr[2:5] = np.nan
   arr[7:8] = np.nan
   sparr = pd.SparseArray(arr)
   sparr

A sparse array can be converted to a regular (dense) ndarray with :meth:`numpy.asarray`

.. ipython:: python

   np.asarray(sparr)

The :attr:`SparseArray.dtype` property stores two pieces of information

1. The dtype of the non-sparse values
2. The scalar fill value

A :class:`SparseDtype` may be constructed by passing each of these

.. ipython:: python

   pd.SparseDtype(np.dtype('datetime64[ns]'))

The default fill value for a given NumPy dtype is the "missing" value for that dtype,
though it may be overridden.

.. ipython:: python

   pd.SparseDtype(np.dtype('datetime64[ns]'),
                  fill_value=pd.Timestamp('2017-01-01'))

Finally, the string alias ``'Sparse[dtype]'`` may be used to specify a sparse dtype
in many places

.. ipython:: python

   pd.array([1, 0, 0, 2], dtype='Sparse[int]')

.. _sparse.accessor:

Sparse Accessor
---------------

.. versionadded:: 0.24.0

Pandas provides a ``.sparse`` accessor, similar to ``.str`` for string data, ``.cat``
for categorical data, and ``.dt`` for datetime-like data. This namespace provides
attributes and methods that are specific to sparse data.

.. ipython:: python

   s = pd.Series([0, 0, 1, 2], dtype="Sparse[int]")
   s.sparse.density
   s.sparse.fill_value

This accessor is available only on data with ``SparseDtype``, and on the :class:`Series`
class itself for creating a Series with sparse data from a scipy COO matrix with.


.. versionadded:: 0.25.0

A ``.sparse`` accessor has been added for :class:`DataFrame` as well.
See :ref:`api.dataframe.sparse` for more.

SparseIndex objects
-------------------

Two kinds of ``SparseIndex`` are implemented, ``block`` and ``integer``. We
recommend using ``block`` as it's more memory efficient. The ``integer`` format
keeps an arrays of all of the locations where the data are not equal to the
fill value. The ``block`` format tracks only the locations and sizes of blocks
of data.

.. _sparse.calculation:

Sparse Calculation
------------------

You can apply NumPy *ufuncs* to ``SparseArray`` and get a ``SparseArray`` as a result.

.. ipython:: python

   arr = pd.SparseArray([1., np.nan, np.nan, -2., np.nan])
   np.abs(arr)


The *ufunc* is also applied to ``fill_value``. This is needed to get
the correct dense result.

.. ipython:: python

   arr = pd.SparseArray([1., -1, -1, -2., -1], fill_value=-1)
   np.abs(arr)
   np.abs(arr).to_dense()

.. _sparse.migration:

Migrating
---------

In older versions of pandas, the ``SparseSeries`` and ``SparseDataFrame`` classes (documented below)
were the preferred way to work with sparse data. With the advent of extension arrays, these subclasses
are no longer needed. Their purpose is better served by using a regular Series or DataFrame with
sparse values instead.

**There's no performance or memory penalty to using a Series or DataFrame with sparse values,
rather than a SparseSeries or SparseDataFrame**.

This section provides some guidance on migrating your code to the new style. As a reminder, you can
use the python warnings module to control warnings. If you wish to ignore the warnings,

.. code-block:: python

   >>> import warnings

   >>> warnings.filterwarnings('ignore', 'Sparse', FutureWarning
   >>> pd.SparseSeries()  # No warning message
   Series([], dtype: Sparse[float64, nan])
   BlockIndex
   Block locations: array([], dtype=int32)
   Block lengths: array([], dtype=int32)

But we recommend modifying your code, rather than ignoring the warning.l

**Construction**

From an array-like, use the regular :class:`Series` or
:class:`DataFrame` constructors with :class:`SparseArray` values.

.. code-block:: python

   # Old way
   >>> pd.SparseDataFrame({"A": [0, 1]})

.. ipython:: python

   # New way
   pd.DataFrame({"A": pd.SparseArray([0, 1])})

From a SciPy sparse matrix, use :meth:`DataFrame.sparse.from_spmatrix`,

.. code-block:: python

   # Old way
   df = pd.SparseDataFrame(sp_matrix, columns=['A', 'B', 'C'])

.. ipython:: python

   # New way
   from scipy import sparse
   mat = sparse.eye(3)
   df = pd.DataFrame.sparse.from_spmatrix(mat, columns=['A', 'B', 'C'])
   df

**Conversion**

From sparse to dense, use the ``.sparse`` accessors

.. ipython:: python

   df.sparse.to_dense()
   df.sparse.to_coo()
   df['A']

From dense to sparse, use :meth:`DataFrame.astype` with a :class:`SparseDtype`.

.. ipython:: python

   dense = pd.DataFrame({"A": [1, 0, 0, 1]})
   dtype = pd.SparseDtype(int, fill_value=0)
   dense.astype(dtype)['A

**Sparse Properties**

Sparse-specific properties, like ``density``, are available on the ``.sparse`` accssor.

.. ipython:: python

   df.sparse.density

The ``SparseDataFrame.default_kind`` and ``SparseDataFrame.default_fill_value`` attributes
have no replacement.

.. _sparse.scipysparse:

Interaction with scipy.sparse
-----------------------------

SparseDataFrame
~~~~~~~~~~~~~~~

.. versionadded:: 0.20.0

Pandas supports creating sparse dataframes directly from ``scipy.sparse`` matrices.

.. ipython:: python
   :okwarning:

   from scipy.sparse import csr_matrix

   arr = np.random.random(size=(1000, 5))
   arr[arr < .9] = 0

   sp_arr = csr_matrix(arr)
   sp_arr

   sdf = pd.SparseDataFrame(sp_arr)
   sdf

All sparse formats are supported, but matrices that are not in :mod:`COOrdinate <scipy.sparse>` format will be converted, copying data as needed.
To convert a ``SparseDataFrame`` back to sparse SciPy matrix in COO format, you can use the :meth:`SparseDataFrame.to_coo` method:

.. ipython:: python

   sdf.to_coo()

SparseSeries
~~~~~~~~~~~~

A :meth:`SparseSeries.to_coo` method is implemented for transforming a ``SparseSeries`` indexed by a ``MultiIndex`` to a ``scipy.sparse.coo_matrix``.

The method requires a ``MultiIndex`` with two or more levels.

.. ipython:: python

   s = pd.Series([3.0, np.nan, 1.0, 3.0, np.nan, np.nan])
   s.index = pd.MultiIndex.from_tuples([(1, 2, 'a', 0),
                                        (1, 2, 'a', 1),
                                        (1, 1, 'b', 0),
                                        (1, 1, 'b', 1),
                                        (2, 1, 'b', 0),
                                        (2, 1, 'b', 1)],
                                       names=['A', 'B', 'C', 'D'])

   s
   # SparseSeries
   ss = s.to_sparse()
   ss

In the example below, we transform the ``SparseSeries`` to a sparse representation of a 2-d array by specifying that the first and second ``MultiIndex`` levels define labels for the rows and the third and fourth levels define labels for the columns. We also specify that the column and row labels should be sorted in the final sparse representation.

.. ipython:: python

   A, rows, columns = ss.to_coo(row_levels=['A', 'B'],
                                column_levels=['C', 'D'],
                                sort_labels=True)

   A
   A.todense()
   rows
   columns

Specifying different row and column labels (and not sorting them) yields a different sparse matrix:

.. ipython:: python

   A, rows, columns = ss.to_coo(row_levels=['A', 'B', 'C'],
                                column_levels=['D'],
                                sort_labels=False)

   A
   A.todense()
   rows
   columns

A convenience method :meth:`SparseSeries.from_coo` is implemented for creating a ``SparseSeries`` from a ``scipy.sparse.coo_matrix``.

.. ipython:: python

   from scipy import sparse
   A = sparse.coo_matrix(([3.0, 1.0, 2.0], ([1, 0, 0], [0, 2, 3])),
                         shape=(3, 4))
   A
   A.todense()

The default behaviour (with ``dense_index=False``) simply returns a ``SparseSeries`` containing
only the non-null entries.

.. ipython:: python

   ss = pd.SparseSeries.from_coo(A)
   ss

Specifying ``dense_index=True`` will result in an index that is the Cartesian product of the
row and columns coordinates of the matrix. Note that this will consume a significant amount of memory
(relative to ``dense_index=False``) if the sparse matrix is large (and sparse) enough.

.. ipython:: python

   ss_dense = pd.SparseSeries.from_coo(A, dense_index=True)
   ss_dense


