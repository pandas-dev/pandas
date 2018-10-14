.. currentmodule:: pandas
.. _sparse:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   import pandas as pd
   import pandas.util.testing as tm
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows = 15

**********************
Sparse data structures
**********************

.. note:: The ``SparsePanel`` class has been removed in 0.19.0

We have implemented "sparse" versions of ``Series`` and ``DataFrame``. These are not sparse
in the typical "mostly 0". Rather, you can view these objects as being "compressed"
where any data matching a specific value (``NaN`` / missing value, though any value
can be chosen) is omitted. A special ``SparseIndex`` object tracks where data has been
"sparsified". This will make much more sense with an example. All of the standard pandas
data structures have a ``to_sparse`` method:

.. ipython:: python

   ts = pd.Series(randn(10))
   ts[2:-2] = np.nan
   sts = ts.to_sparse()
   sts

The ``to_sparse`` method takes a ``kind`` argument (for the sparse index, see
below) and a ``fill_value``. So if we had a mostly zero ``Series``, we could
convert it to sparse with ``fill_value=0``:

.. ipython:: python

   ts.fillna(0).to_sparse(fill_value=0)

The sparse objects exist for memory efficiency reasons. Suppose you had a
large, mostly NA ``DataFrame``:

.. ipython:: python

   df = pd.DataFrame(randn(10000, 4))
   df.iloc[:9998] = np.nan
   sdf = df.to_sparse()
   sdf
   sdf.density

As you can see, the density (% of values that have not been "compressed") is
extremely low. This sparse object takes up much less memory on disk (pickled)
and in the Python interpreter. Functionally, their behavior should be nearly
identical to their dense counterparts.

Any sparse object can be converted back to the standard dense form by calling
``to_dense``:

.. ipython:: python

   sts.to_dense()

.. _sparse.accessor:

Sparse Accessor
---------------

Pandas provides a ``.sparse`` accessor, similar to ``.str`` for string data, ``.cat``
for categorical data, and ``.dt`` for datetime-like data. This namespace provides
attributes and methods that are specific to sparse data.

.. ipython:: python

   s = pd.Series([0, 0, 1, 2], dtype="Sparse[int]")
   s.sparse.density
   s.sparse.fill_value

This accessor is available only on data with ``SparseDtype``, and on the :class:`Series`
class itself for creating a Series with sparse data from a scipy COO matrix with.

.. _sparse.array:

SparseArray
-----------

``SparseArray`` is the base layer for all of the sparse indexed data
structures. It is a 1-dimensional ndarray-like object storing only values
distinct from the ``fill_value``:

.. ipython:: python

   arr = np.random.randn(10)
   arr[2:5] = np.nan; arr[7:8] = np.nan
   sparr = pd.SparseArray(arr)
   sparr

Like the indexed objects (SparseSeries, SparseDataFrame), a ``SparseArray``
can be converted back to a regular ndarray by calling ``to_dense``:

.. ipython:: python

   sparr.to_dense()


SparseIndex objects
-------------------

Two kinds of ``SparseIndex`` are implemented, ``block`` and ``integer``. We
recommend using ``block`` as it's more memory efficient. The ``integer`` format
keeps an arrays of all of the locations where the data are not equal to the
fill value. The ``block`` format tracks only the locations and sizes of blocks
of data.

.. _sparse.dtype:

Sparse Dtypes
-------------

Sparse data should have the same dtype as its dense representation. Currently,
``float64``, ``int64`` and ``bool`` dtypes are supported. Depending on the original
dtype, ``fill_value`` default changes:

* ``float64``: ``np.nan``
* ``int64``: ``0``
* ``bool``: ``False``

.. ipython:: python

   s = pd.Series([1, np.nan, np.nan])
   s
   s.to_sparse()

   s = pd.Series([1, 0, 0])
   s
   s.to_sparse()

   s = pd.Series([True, False, True])
   s
   s.to_sparse()

You can change the dtype using ``.astype()``, the result is also sparse. Note that
``.astype()`` also affects to the ``fill_value`` to keep its dense representation.


.. ipython:: python

   s = pd.Series([1, 0, 0, 0, 0])
   s
   ss = s.to_sparse()
   ss
   ss.astype(np.float64)

It raises if any value cannot be coerced to specified dtype.

.. code-block:: ipython

   In [1]: ss = pd.Series([1, np.nan, np.nan]).to_sparse()
   0    1.0
   1    NaN
   2    NaN
   dtype: float64
   BlockIndex
   Block locations: array([0], dtype=int32)
   Block lengths: array([1], dtype=int32)

   In [2]: ss.astype(np.int64)
   ValueError: unable to coerce current fill_value nan to int64 dtype

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

.. _sparse.scipysparse:

Interaction with scipy.sparse
-----------------------------

SparseDataFrame
~~~~~~~~~~~~~~~

.. versionadded:: 0.20.0

Pandas supports creating sparse dataframes directly from ``scipy.sparse`` matrices.

.. ipython:: python

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
   :suppress:


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
   :suppress:

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
