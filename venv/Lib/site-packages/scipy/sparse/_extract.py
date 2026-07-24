"""Functions to extract parts of sparse matrices
"""

__docformat__ = "restructuredtext en"

__all__ = ['find', 'tril', 'triu']

import os
from warnings import warn

from ._coo import coo_array, coo_matrix
from ._base import sparray, spmatrix


def find(A):
    """Return the indices and values of the nonzero elements of a matrix.

    Parameters
    ----------
    A : dense or sparse array or matrix
        Matrix whose nonzero elements are desired.

    Returns
    -------
    (I,J,V) : tuple of arrays
        I,J, and V contain the row indices, column indices, and values
        of the nonzero entries.


    Examples
    --------
    >>> from scipy.sparse import csr_array, find
    >>> A = csr_array([[7.0, 8.0, 0],[0, 0, 9.0]])
    >>> find(A)
    (array([0, 0, 1], dtype=int32),
     array([0, 1, 2], dtype=int32),
     array([ 7.,  8.,  9.]))

    """

    A = coo_array(A, copy=True)
    A.sum_duplicates()
    # remove explicit zeros
    nz_mask = A.data != 0
    return A.row[nz_mask], A.col[nz_mask], A.data[nz_mask]


def tril(A, k=0, format=None):
    """Return the lower triangular portion of a sparse array or matrix.

    Returns the elements on or below the k-th diagonal of A.
        - k = 0 corresponds to the main diagonal
        - k > 0 is above the main diagonal
        - k < 0 is below the main diagonal

    .. warning::

        `tril` is switching to the sparse array interface.

        For the case where no input arrays are sparse, this function is
        switching to returning a sparse array instead of sparse matrix.
        Control the sparse return class by making at least one input sparse,
        e.g., ``tril(coo_matrix(A))``, or ``tril(coo_array(A))``.
        That removes any deprecation warnings as well.
        For more general information about sparrays, see
        :ref:`Migration from spmatrix to sparray <migration_to_sparray>`.
        Handling of this no sparse input case will change no earlier than v1.20.

    Parameters
    ----------
    A : dense or sparse array or matrix
        Matrix whose lower trianglar portion is desired.
    k : int : optional
        The top-most diagonal of the lower triangle.
    format : str
        Sparse format of the result, e.g. format="csr", etc.

    Returns
    -------
    L : sparse matrix
        Lower triangular portion of A in sparse format.

    See Also
    --------
    triu : upper triangle in sparse format

    Examples
    --------
    >>> from scipy.sparse import csr_array, tril
    >>> A = csr_array([[1, 2, 0, 0, 3], [4, 5, 0, 6, 7], [0, 0, 8, 9, 0]],
    ...               dtype='int32')
    >>> A.toarray()
    array([[1, 2, 0, 0, 3],
           [4, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]], dtype=int32)
    >>> tril(A).toarray()
    array([[1, 0, 0, 0, 0],
           [4, 5, 0, 0, 0],
           [0, 0, 8, 0, 0]], dtype=int32)
    >>> tril(A).nnz
    4
    >>> tril(A, k=1).toarray()
    array([[1, 2, 0, 0, 0],
           [4, 5, 0, 0, 0],
           [0, 0, 8, 9, 0]], dtype=int32)
    >>> tril(A, k=-1).toarray()
    array([[0, 0, 0, 0, 0],
           [4, 0, 0, 0, 0],
           [0, 0, 0, 0, 0]], dtype=int32)
    >>> tril(A, format='csc')
    <Compressed Sparse Column sparse array of dtype 'int32'
        with 4 stored elements and shape (3, 5)>

    """
    if isinstance(A, sparray):
        coo_sparse = coo_array
    elif isinstance(A, spmatrix):
        coo_sparse = coo_matrix
    else:  # dense
        msg = """`tril` is switching to the sparse array interface.

        For the case where input arrays are numpy arrays, this function is
        switching to returning a sparse array instead of sparse matrix.
        Recover the sparse matrix return value by making one input a sparse matrix.
        For example, tril(coo_matrix(A)).
        Avoid this message for sparse array output by using tril(coo_array(A)).
        For more information, see the spmatrix to sparray migration guide
        https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html

        This function will be changed no earlier than v1.20.
        """
        prefixes = (os.path.dirname(__file__),)
        warn(msg, category=DeprecationWarning, skip_file_prefixes=prefixes)
        # default when input is ndarray
        coo_sparse = coo_matrix

    # convert to COOrdinate format where things are easy
    A = coo_sparse(A, copy=False)
    mask = A.row + k >= A.col

    row = A.row[mask]
    col = A.col[mask]
    data = A.data[mask]
    new_coo = coo_sparse((data, (row, col)), shape=A.shape, dtype=A.dtype)
    return new_coo.asformat(format)


def triu(A, k=0, format=None):
    """Return the upper triangular portion of a sparse array or matrix.

    Returns the elements on or above the k-th diagonal of A.
        - k = 0 corresponds to the main diagonal
        - k > 0 is above the main diagonal
        - k < 0 is below the main diagonal

    .. warning::

        `triu` is switching to the sparse array interface.

        For the case where no input arrays are sparse, this function is
        switching to returning a sparse array instead of sparse matrix.
        Control the sparse return class by making at least one input sparse,
        e.g., ``triu(coo_matrix(A))``, or ``triu(coo_array(A))``.
        That removes any deprecation warnings as well.
        For more general information about sparrays, see
        :ref:`Migration from spmatrix to sparray <migration_to_sparray>`.
        Handling of this no sparse input case will change no earlier than v1.20.

    Parameters
    ----------
    A : dense or sparse array or matrix
        Matrix whose upper trianglar portion is desired.
    k : int : optional
        The bottom-most diagonal of the upper triangle.
    format : str
        Sparse format of the result, e.g. format="csr", etc.

    Returns
    -------
    L : sparse array or matrix
        Upper triangular portion of A in sparse format.
        Sparse array if A is a sparse array, otherwise matrix.

    See Also
    --------
    tril : lower triangle in sparse format

    Examples
    --------
    >>> from scipy.sparse import csr_array, triu
    >>> A = csr_array([[1, 2, 0, 0, 3], [4, 5, 0, 6, 7], [0, 0, 8, 9, 0]],
    ...                dtype='int32')
    >>> A.toarray()
    array([[1, 2, 0, 0, 3],
           [4, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]], dtype=int32)
    >>> triu(A).toarray()
    array([[1, 2, 0, 0, 3],
           [0, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]], dtype=int32)
    >>> triu(A).nnz
    8
    >>> triu(A, k=1).toarray()
    array([[0, 2, 0, 0, 3],
           [0, 0, 0, 6, 7],
           [0, 0, 0, 9, 0]], dtype=int32)
    >>> triu(A, k=-1).toarray()
    array([[1, 2, 0, 0, 3],
           [4, 5, 0, 6, 7],
           [0, 0, 8, 9, 0]], dtype=int32)
    >>> triu(A, format='csc')
    <Compressed Sparse Column sparse array of dtype 'int32'
        with 8 stored elements and shape (3, 5)>

    """
    if isinstance(A, sparray):
        coo_sparse = coo_array
    elif isinstance(A, spmatrix):
        coo_sparse = coo_matrix
    else:  # dense
        msg = """`triu` is switching to the sparse array interface.

        For the case where input arrays are numpy arrays, this function is
        switching to returning a sparse array instead of sparse matrix.
        Recover the sparse matrix return value by making one input a sparse matrix.
        For example, triu(coo_matrix(A)).
        Avoid this message for sparse array output by using triu(coo_array(A)).
        For more information, see the spmatrix to sparray migration guide
        https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html

        This function will be changed no earlier than v1.20.
        """
        prefixes = (os.path.dirname(__file__),)
        warn(msg, category=DeprecationWarning, skip_file_prefixes=prefixes)
        # default when input is ndarray
        coo_sparse = coo_matrix

    # convert to COOrdinate format where things are easy
    A = coo_sparse(A, copy=False)
    mask = A.row + k <= A.col

    row = A.row[mask]
    col = A.col[mask]
    data = A.data[mask]
    new_coo = coo_sparse((data, (row, col)), shape=A.shape, dtype=A.dtype)
    return new_coo.asformat(format)
