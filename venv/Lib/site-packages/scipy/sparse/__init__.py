"""
===================================
Sparse arrays (:mod:`scipy.sparse`)
===================================

.. currentmodule:: scipy.sparse

.. toctree::
   :hidden:

   sparse.csgraph
   sparse.linalg
   sparse.migration_to_sparray
   sparse.spmatrix_api

SciPy 2-D sparse array package for numeric data.

.. warning::

   SciPy sparse is shifting from a sparse matrix interface to a sparse
   array interface. In the next few releases we expect to deprecate the
   sparse matrix interface. For documentation of the matrix
   interface, see the :ref:`spmatrix interface docs <spmatrix_api>`.
   For guidance on converting existing code to sparse arrays, see
   :ref:`Migration from spmatrix to sparray <migration_to_sparray>`.


Submodules
==========

.. autosummary::

   csgraph - Compressed sparse graph routines
   linalg - Sparse linear algebra routines


Sparse array classes
====================

.. autosummary::
   :toctree: generated/

   bsr_array - Block Sparse Row array
   coo_array - A sparse array in COOrdinate format
   csc_array - Compressed Sparse Column array
   csr_array - Compressed Sparse Row array
   dia_array - Sparse array with DIAgonal storage
   dok_array - Dictionary Of Keys based sparse array
   lil_array - Row-based list of lists sparse array
   sparray - Sparse array base class

.. _sparse-construction-functions:

Building sparse arrays
----------------------

.. autosummary::
   :toctree: generated/

   diags_array - Return a sparse array from diagonals
   eye_array - Sparse MxN array whose k-th diagonal is all ones
   random_array - Random values in a given shape array
   block_array - Build a sparse array from sub-blocks

.. _combining-arrays:

Combining arrays
----------------

.. autosummary::
   :toctree: generated/

   kron - Kronecker product of two sparse arrays
   kronsum - Kronecker sum of sparse arrays
   block_diag - Build a block diagonal sparse array
   tril - Lower triangular portion of a sparse array
   triu - Upper triangular portion of a sparse array
   hstack - Stack sparse arrays horizontally (column wise)
   vstack - Stack sparse arrays vertically (row wise)

Manipulating arrays
-------------------

.. autosummary::
   :toctree: generated/

   matrix_transpose - Transpose a matrix (or a batch of matrices)
   swapaxes - swap two axes of a sparse array
   expand_dims - add a new (trivial) axis to a sparse array
   permute_dims - reorder the axes of a sparse array

Sparse tools
------------

.. autosummary::
   :toctree: generated/

   save_npz - Save a sparse array to a file using ``.npz`` format.
   load_npz - Load a sparse array from a file using ``.npz`` format.
   find - Return the indices and values of the nonzero elements
   get_index_dtype - determine a good dtype for index arrays.
   safely_cast_index_arrays - cast index array dtype or raise if shape too big

Identifying sparse arrays
-------------------------

.. autosummary::
   :toctree: generated/

   issparse - Check if the argument is a sparse object (array or matrix).
   isspmatrix - Check if the argument is an old style sparse matrix.


Warnings
========

.. autosummary::
   :toctree: generated/

   SparseEfficiencyWarning
   SparseWarning


Usage information
=================

There are seven available sparse array types:

1. csc_array: Compressed Sparse Column format
2. csr_array: Compressed Sparse Row format
3. bsr_array: Block Sparse Row format
4. lil_array: List of Lists format
5. dok_array: Dictionary of Keys format
6. coo_array: COOrdinate format (aka IJV, triplet format)
7. dia_array: DIAgonal format

To construct an array efficiently, use any of `coo_array`,
`dok_array` or `lil_array`. They each support basic slicing
and fancy indexing with a similar syntax to NumPy arrays.
The COO format is recommended and other formats use it under the hood to
allow efficient construction using data values and coordinate arrays.

Despite their similarity to NumPy arrays, it is **strongly discouraged**
to use NumPy functions directly on these arrays because NumPy typically
treats them as generic Python objects rather than arrays, leading to
unexpected (and incorrect) results. If you are tempted to apply a NumPy
function to these arrays, check if SciPy has its own implementation
for the given sparse array class and, if not, **convert the sparse array
to a NumPy array** (e.g., using the `toarray` method of the class)
before applying the method.

All conversions among the CSR, CSC, and COO formats are efficient,
linear-time operations.

To perform manipulations such as multiplication or inversion, first
convert the array to either CSC or CSR format. The `lil_array`
format is row-based, so conversion to CSR is efficient, but
conversion to CSC is less so.

Matrix vector product
---------------------

To do a vector product between a 2D sparse array and a vector use
the matmul operator (i.e., ``@``) which performs a dot product (like the
``dot`` method):

>>> import numpy as np
>>> from scipy.sparse import csr_array
>>> A = csr_array([[1, 2, 0], [0, 0, 3], [4, 0, 5]])
>>> v = np.array([1, 0, -1])
>>> A @ v
array([ 1, -3, -1], dtype=int64)

The CSR format is especially suitable for fast matrix vector products.

Example 1
---------

Construct a 1000x1000 `lil_array` and add some values to it:

>>> from scipy.sparse import lil_array
>>> from scipy.sparse.linalg import spsolve
>>> from numpy.linalg import solve, norm
>>> from numpy.random import rand

>>> A = lil_array((1000, 1000))
>>> A[0, :100] = rand(100)
>>> A.setdiag(rand(1000))

Now convert it to CSR format and solve A x = b for x:

>>> A = A.tocsr()
>>> b = rand(1000)
>>> x = spsolve(A, b)

Convert it to a dense array and solve, and check that the result
is the same:

>>> x_ = solve(A.toarray(), b)

Now we can compute norm of the error with:

>>> err = norm(x-x_)
>>> err < 1e-9
True

It should be small :)


Example 2
---------

Construct an array in COO format:

>>> from scipy import sparse
>>> from numpy import array
>>> I = array([0,3,1,0])
>>> J = array([0,3,1,2])
>>> V = array([4,5,7,9])
>>> A = sparse.coo_array((V,(I,J)),shape=(4,4))

Notice that the indices do not need to be sorted.

Duplicate (i,j) entries are summed when converting to CSR or CSC.

>>> I = array([0,0,1,3,1,0,0])
>>> J = array([0,2,1,3,1,0,0])
>>> V = array([1,1,1,1,1,1,1])
>>> B = sparse.coo_array((V,(I,J)),shape=(4,4)).tocsr()

This is useful for constructing finite-element stiffness and mass matrices.

Further details
---------------

CSR column indices are not necessarily sorted. Likewise for CSC row
indices. And similarly for COO coordinates. Use the ``.sorted_indices()``
and ``.sort_indices()`` methods when sorted indices are required (e.g.,
when passing data to other libraries).

"""

# Original code by Travis Oliphant.
# Modified and extended by Ed Schofield, Robert Cimrman,
# Nathan Bell, and Jake Vanderplas.

import warnings as _warnings
import importlib as _importlib

from ._base import *
from ._csr import *
from ._csc import *
from ._lil import *
from ._dok import *
from ._coo import *
from ._dia import *
from ._bsr import *
from ._construct import *
from ._extract import *
from ._matrix import spmatrix
from ._matrix_io import *
from ._sputils import get_index_dtype, safely_cast_index_arrays


# Deprecated namespaces, to be removed in v2.0.0
from . import (
    base, bsr, compressed, construct, coo, csc, csr, data, dia, dok, extract,
    lil, sparsetools, sputils
)

_submodules = ["csgraph", "linalg"]

__all__ = [s for s in dir() if not s.startswith('_')] + _submodules

# Filter PendingDeprecationWarning for np.matrix introduced with numpy 1.15
msg = 'the matrix subclass is not the recommended way'
_warnings.filterwarnings('ignore', message=msg)

def __dir__():
   return __all__


def __getattr__(name):
    if name in _submodules:
        return _importlib.import_module(f'scipy.sparse.{name}')
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(
                f"Module 'scipy.sparse' has no attribute '{name}'"
            )

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
