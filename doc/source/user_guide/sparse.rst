.. _sparse:

{{ header }}

**********************
Sparse data structures
**********************

Pandas provides data structures for efficiently storing sparse data.
These are not necessarily sparse in the typical "mostly 0". Rather, you can view these
objects as being "compressed" where any data matching a specific value (``NaN`` / missing value, though any value
can be chosen, including 0) is omitted. The compressed values are not actually stored in the array.

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
   sdf.head()
   sdf.dtypes
   sdf.sparse.density

As you can see, the density (% of values that have not been "compressed") is
extremely low. This sparse object takes up much less memory on disk (pickled)
and in the Python interpreter.

.. ipython:: python

   'dense : {:0.2f} bytes'.format(df.memory_usage().sum() / 1e3)
   'sparse: {:0.2f} bytes'.format(sdf.memory_usage().sum() / 1e3)

Functionally, their behavior should be nearly
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


.. _sparse.dtype:

SparseDtype
-----------

The :attr:`SparseArray.dtype` property stores two pieces of information

1. The dtype of the non-sparse values
2. The scalar fill value


.. ipython:: python

   sparr.dtype


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

Sparse accessor
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
See :ref:`api.frame.sparse` for more.

.. _sparse.calculation:

Sparse calculation
------------------

You can apply NumPy `ufuncs <https://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_
to ``SparseArray`` and get a ``SparseArray`` as a result.

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

.. note::

   ``SparseSeries`` and ``SparseDataFrame`` were removed in pandas 1.0.0. This migration
   guide is present to aid in migrating from previous versions.

In older versions of pandas, the ``SparseSeries`` and ``SparseDataFrame`` classes (documented below)
were the preferred way to work with sparse data. With the advent of extension arrays, these subclasses
are no longer needed. Their purpose is better served by using a regular Series or DataFrame with
sparse values instead.

.. note::

  There's no performance or memory penalty to using a Series or DataFrame with sparse values,
  rather than a SparseSeries or SparseDataFrame.

This section provides some guidance on migrating your code to the new style. As a reminder,
you can use the python warnings module to control warnings. But we recommend modifying
your code, rather than ignoring the warning.

**Construction**

From an array-like, use the regular :class:`Series` or
:class:`DataFrame` constructors with :class:`SparseArray` values.

.. code-block:: python

   # Previous way
   >>> pd.SparseDataFrame({"A": [0, 1]})

.. ipython:: python

   # New way
   pd.DataFrame({"A": pd.SparseArray([0, 1])})

From a SciPy sparse matrix, use :meth:`DataFrame.sparse.from_spmatrix`,

.. code-block:: python

   # Previous way
   >>> from scipy import sparse
   >>> mat = sparse.eye(3)
   >>> df = pd.SparseDataFrame(mat, columns=['A', 'B', 'C'])

.. ipython:: python

   # New way
   from scipy import sparse
   mat = sparse.eye(3)
   df = pd.DataFrame.sparse.from_spmatrix(mat, columns=['A', 'B', 'C'])
   df.dtypes

**Conversion**

From sparse to dense, use the ``.sparse`` accessors

.. ipython:: python

   df.sparse.to_dense()
   df.sparse.to_coo()

From dense to sparse, use :meth:`DataFrame.astype` with a :class:`SparseDtype`.

.. ipython:: python

   dense = pd.DataFrame({"A": [1, 0, 0, 1]})
   dtype = pd.SparseDtype(int, fill_value=0)
   dense.astype(dtype)

**Sparse Properties**

Sparse-specific properties, like ``density``, are available on the ``.sparse`` accessor.

.. ipython:: python

   df.sparse.density

**General differences**

In a ``SparseDataFrame``, *all* columns were sparse. A :class:`DataFrame` can have a mixture of
sparse and dense columns. As a consequence, assigning new columns to a ``DataFrame`` with sparse
values will not automatically convert the input to be sparse.

.. code-block:: python

   # Previous Way
   >>> df = pd.SparseDataFrame({"A": [0, 1]})
   >>> df['B'] = [0, 0]  # implicitly becomes Sparse
   >>> df['B'].dtype
   Sparse[int64, nan]

Instead, you'll need to ensure that the values being assigned are sparse

.. ipython:: python

   df = pd.DataFrame({"A": pd.SparseArray([0, 1])})
   df['B'] = [0, 0]  # remains dense
   df['B'].dtype
   df['B'] = pd.SparseArray([0, 0])
   df['B'].dtype

The ``SparseDataFrame.default_kind`` and ``SparseDataFrame.default_fill_value`` attributes
have no replacement.

.. _sparse.scipysparse:

Interaction with scipy.sparse
-----------------------------

Use :meth:`DataFrame.sparse.from_spmatrix` to create a ``DataFrame`` with sparse values from a sparse matrix.

.. versionadded:: 0.25.0

.. ipython:: python

   from scipy.sparse import csr_matrix

   arr = np.random.random(size=(1000, 5))
   arr[arr < .9] = 0

   sp_arr = csr_matrix(arr)
   sp_arr

   sdf = pd.DataFrame.sparse.from_spmatrix(sp_arr)
   sdf.head()
   sdf.dtypes

All sparse formats are supported, but matrices that are not in :mod:`COOrdinate <scipy.sparse>` format will be converted, copying data as needed.
To convert back to sparse SciPy matrix in COO format, you can use the :meth:`DataFrame.sparse.to_coo` method:

.. ipython:: python

   sdf.sparse.to_coo()

meth:`Series.sparse.to_coo` is implemented for transforming a ``Series`` with sparse values indexed by a :class:`MultiIndex` to a :class:`scipy.sparse.coo_matrix`.

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
   ss = s.astype('Sparse')
   ss

In the example below, we transform the ``Series`` to a sparse representation of a 2-d array by specifying that the first and second ``MultiIndex`` levels define labels for the rows and the third and fourth levels define labels for the columns. We also specify that the column and row labels should be sorted in the final sparse representation.

.. ipython:: python

   A, rows, columns = ss.sparse.to_coo(row_levels=['A', 'B'],
                                       column_levels=['C', 'D'],
                                       sort_labels=True)

   A
   A.todense()
   rows
   columns

Specifying different row and column labels (and not sorting them) yields a different sparse matrix:

.. ipython:: python

   A, rows, columns = ss.sparse.to_coo(row_levels=['A', 'B', 'C'],
                                       column_levels=['D'],
                                       sort_labels=False)

   A
   A.todense()
   rows
   columns

A convenience method :meth:`Series.sparse.from_coo` is implemented for creating a ``Series`` with sparse values from a ``scipy.sparse.coo_matrix``.

.. ipython:: python

   from scipy import sparse
   A = sparse.coo_matrix(([3.0, 1.0, 2.0], ([1, 0, 0], [0, 2, 3])),
                         shape=(3, 4))
   A
   A.todense()

The default behaviour (with ``dense_index=False``) simply returns a ``Series`` containing
only the non-null entries.

.. ipython:: python

   ss = pd.Series.sparse.from_coo(A)
   ss

Specifying ``dense_index=True`` will result in an index that is the Cartesian product of the
row and columns coordinates of the matrix. Note that this will consume a significant amount of memory
(relative to ``dense_index=False``) if the sparse matrix is large (and sparse) enough.

.. ipython:: python

   ss_dense = pd.Series.sparse.from_coo(A, dense_index=True)
   ss_dense
