"""SVD decomposition functions."""
import numpy as np

from scipy._lib._util import _apply_over_batch, _deprecate_dtypes
from . import _batched_linalg

# Local imports.
from ._misc import LinAlgError, _datacopied
from .lapack import _normalize_lapack_dtype, _ensure_aligned_and_native, HAS_ILP64
from scipy.linalg.lapack import get_lapack_funcs   # noqa: F401  (backwards compat)
from ._decomp import _asarray_validated


__all__ = ['svd', 'svdvals', 'diagsvd', 'orth', 'subspace_angles', 'null_space']


def _format_emit_errors_warnings(err_lst, lapack_driver):
    """Format/emit errors/warnings from a lowlevel batched routine.
    """
    # NB the low-level routine currently stops processing a batch at the first error
    for entry in err_lst:
        info = entry["lapack_info"]
        num = entry["num"]
        if info != 0:
            if info > 0:
                raise LinAlgError(f"SVD did not converge for slice = {num}.")
            if info < 0:
                if lapack_driver == "gesdd" and info == -4:
                    msg = f"slice {num} has a NaN entry"
                    raise ValueError(msg)
                raise ValueError(
                    f'illegal value in {-info}th argument of internal {lapack_driver}'
                    f'  for slice {num}.'
                )


def svd(a, full_matrices=True, compute_uv=True, overwrite_a=False,
        check_finite=True, lapack_driver='gesdd'):
    """
    Singular Value Decomposition.

    Factorizes the matrix `a` into two unitary matrices ``U`` and ``Vh``, and
    a 1-D array ``s`` of singular values (real, non-negative) such that
    ``a == U @ S @ Vh``, where ``S`` is a suitably shaped matrix of zeros with
    main diagonal ``s``.

    Parameters
    ----------
    a : (..., M, N) array_like
        Matrix to decompose.
    full_matrices : bool, optional
        If True (default), `U` and `Vh` are of shape ``(M, M)``, ``(N, N)``.
        If False, the shapes are ``(M, K)`` and ``(K, N)``, where
        ``K = min(M, N)``.
    compute_uv : bool, optional
        Whether to compute also ``U`` and ``Vh`` in addition to ``s``.
        Default is True.
    overwrite_a : bool, optional
        Whether to overwrite data in `a` (may improve performance). Default is False.
        See :ref:`tutorial_linalg_overwrite` for details.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    lapack_driver : {'gesdd', 'gesvd'}, optional
        Whether to use the more efficient divide-and-conquer approach
        (``'gesdd'``) or general rectangular approach (``'gesvd'``)
        to compute the SVD. MATLAB and Octave use the ``'gesvd'`` approach.
        Default is ``'gesdd'``.

    Returns
    -------
    U : ndarray
        Unitary matrix having left singular vectors as columns.
        Of shape ``(M, M)`` or ``(M, K)``, depending on `full_matrices`.
        Only present when ``compute_uv=True``.
    s : ndarray
        The singular values, sorted in non-increasing order.
        Of shape (K,), with ``K = min(M, N)``.
    Vh : ndarray
        Unitary matrix having right singular vectors as rows.
        Of shape ``(N, N)`` or ``(K, N)`` depending on `full_matrices`.
        Only present when ``compute_uv=True``.

    Raises
    ------
    LinAlgError
        If SVD computation does not converge.

    See Also
    --------
    svdvals : Compute singular values of a matrix.
    diagsvd : Construct the Sigma matrix, given the vector s.

    Notes
    -----
    The array argument of this function, `a`, may have additional
    "batch" dimensions prepended to the core shape. In this case, the array is treated
    as a batch of lower-dimensional slices; see :ref:`linalg_batch` for details.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import linalg
    >>> rng = np.random.default_rng()
    >>> m, n = 9, 6
    >>> a = rng.standard_normal((m, n)) + 1.j*rng.standard_normal((m, n))
    >>> U, s, Vh = linalg.svd(a)
    >>> U.shape,  s.shape, Vh.shape
    ((9, 9), (6,), (6, 6))

    Reconstruct the original matrix from the decomposition:

    >>> sigma = np.zeros((m, n))
    >>> for i in range(min(m, n)):
    ...     sigma[i, i] = s[i]
    >>> a1 = U @ sigma @ Vh
    >>> np.allclose(a, a1)
    True

    Alternatively, use ``full_matrices=False`` (notice that the shape of
    ``U`` is then ``(m, n)`` instead of ``(m, m)``):

    >>> U, s, Vh = linalg.svd(a, full_matrices=False)
    >>> U.shape, s.shape, Vh.shape
    ((9, 6), (6,), (6, 6))
    >>> S = np.diag(s)
    >>> np.allclose(a, U @ S @ Vh)
    True

    >>> s2 = linalg.svd(a, compute_uv=False)
    >>> np.allclose(s, s2)
    True

    If the input matrix has more than two dimensions, it is interpreted as a batch of
    two-dimensional matrices:

    >>> aa = np.stack((a, 2*a))
    >>> linalg.svdvals(aa)[0] == linalg.svdvals(a)
    array([ True,  True,  True,  True,  True,  True])
    >>> linalg.svdvals(aa)[1] == 2 * linalg.svdvals(a)
    array([ True,  True,  True,  True,  True,  True])
    """
    if not isinstance(lapack_driver, str):
        raise TypeError('lapack_driver must be a string')
    if lapack_driver not in ('gesdd', 'gesvd'):
        message = f'lapack_driver must be "gesdd" or "gesvd", not "{lapack_driver}"'
        raise ValueError(message)

    # basic sanity checks of the input matrix
    a1 = _asarray_validated(a, check_finite=check_finite)
    _deprecate_dtypes("svd", a1)

    if a1.ndim < 2:
        raise ValueError(f"Expected at least ndim=2, got {a1.ndim=}")

    m, n = a1.shape[-2], a1.shape[-1]

    # Also check if dtype is LAPACK compatible
    a1, overwrite_a = _normalize_lapack_dtype(a1, overwrite_a)
    a1, overwrite_a = _ensure_aligned_and_native(a1, overwrite_a)

    overwrite_a = overwrite_a or (_datacopied(a1, a))
    overwrite_a = overwrite_a and (a1.ndim == 2) and (a1.flags["F_CONTIGUOUS"])

    # accommodate empty matrix
    if a1.size == 0:
        u0, s0, v0 = svd(np.eye(2, dtype=a1.dtype))

        batch_shape = a1.shape[:-2]
        s = np.empty_like(a1, shape=batch_shape + (0,), dtype=s0.dtype)
        if full_matrices:
            u = np.empty_like(a1, shape=batch_shape + (m, m), dtype=u0.dtype)
            u[...] = np.identity(m)
            v = np.empty_like(a1, shape=batch_shape + (n, n), dtype=v0.dtype)
            v[...] = np.identity(n)
        else:
            u = np.empty_like(a1, shape=batch_shape + (m, 0), dtype=u0.dtype)
            v = np.empty_like(a1, shape=batch_shape + (0, n), dtype=v0.dtype)
        if compute_uv:
            return u, s, v
        else:
            return s

    if compute_uv:
        max_mn, min_mn = (m, n) if m > n else (n, m)
        if full_matrices:
            if not HAS_ILP64 and max_mn*max_mn > np.iinfo(np.int32).max:
                raise ValueError(f"Indexing a matrix size {max_mn} x {max_mn} "
                                  "would incur integer overflow in LAPACK. "
                                  "Instead, either use using numpy.linalg.svd or build"
                                  "SciPy with ILP64 support.")
        else:
            sz = max(m * min_mn, n * min_mn)
            if not HAS_ILP64 and max(m * min_mn, n * min_mn) > np.iinfo(np.int32).max:
                raise ValueError(f"Indexing a matrix of {sz} elements would "
                                  "incur an in integer overflow in LAPACK. "
                                  "Instead, either use using numpy.linalg.svd or build"
                                  "SciPy with ILP64 support.")

    res = _batched_linalg._svd(
        a1, lapack_driver, compute_uv, full_matrices, overwrite_a
    )

    err_lst = res[-1]
    if err_lst:
        _format_emit_errors_warnings(err_lst, lapack_driver)

    if compute_uv:
        return res[:-1]    # u, s, v
    else:
        return res[0]   # s


def svdvals(a, overwrite_a=False, check_finite=True):
    """
    Compute singular values of a matrix.

    Parameters
    ----------
    a : (M, N) array_like
        Matrix to decompose.
    overwrite_a : bool, optional
        Whether to overwrite data in `a` (may improve performance). Default is False.
        See :ref:`tutorial_linalg_overwrite` for details.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    s : (min(M, N),) ndarray
        The singular values, sorted in decreasing order.

    Raises
    ------
    LinAlgError
        If SVD computation does not converge.

    See Also
    --------
    svd : Compute the full singular value decomposition of a matrix.
    diagsvd : Construct the Sigma matrix, given the vector s.

    Notes
    -----
    Array argument of this function, `a`, may have additional
    "batch" dimensions prepended to the core shape. In this case, the array is treated
    as a batch of lower-dimensional slices; see :ref:`linalg_batch` for details.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import svdvals
    >>> m = np.array([[1.0, 0.0],
    ...               [2.0, 3.0],
    ...               [1.0, 1.0],
    ...               [0.0, 2.0],
    ...               [1.0, 0.0]])
    >>> svdvals(m)
    array([ 4.28091555,  1.63516424])

    If the input matrix has more than two dimensions, it is interpreted as a batch of
    two-dimensional matrices:

    >>> mm = np.stack((m, 2*m))
    >>> svdvals(mm)
    array([[4.28091555, 1.63516424],
           [8.56183109, 3.27032847]])

    We can verify the maximum singular value of `m` by computing the maximum
    length of `m @ u` over all the unit vectors `u` in the (x,y) plane.
    We approximate "all" the unit vectors with a large sample. Because
    of linearity, we only need the unit vectors with angles in ``[0, pi]``.

    >>> t = np.linspace(0, np.pi, 2000)
    >>> u = np.array([np.cos(t), np.sin(t)])
    >>> np.linalg.norm(m @ u, axis=0).max()
    4.2809152422538475

    `p` is a projection matrix with rank 1. With exact arithmetic,
    its singular values would be ``[1, 0, 0, 0]``.

    >>> v = np.array([0.1, 0.3, 0.9, 0.3])
    >>> p = np.outer(v, v)
    >>> svdvals(p)
    array([  1.00000000e+00,   2.02021698e-17,   1.56692500e-17,
             8.15115104e-34])

    The singular values of an orthogonal matrix are all 1. Here, we
    create a random orthogonal matrix by using the ``rvs()`` method of
    `scipy.stats.ortho_group`.

    >>> from scipy.stats import ortho_group
    >>> orth = ortho_group.rvs(4)
    >>> svdvals(orth)
    array([ 1.,  1.,  1.,  1.])
    """
    return svd(a, compute_uv=0, overwrite_a=overwrite_a,
               check_finite=check_finite)


@_apply_over_batch(('s', 1))
def diagsvd(s, M, N):
    """
    Construct the sigma matrix in SVD from singular values and size M, N.

    Parameters
    ----------
    s : (M,) or (N,) array_like
        Singular values
    M : int
        Size of the matrix whose singular values are `s`.
    N : int
        Size of the matrix whose singular values are `s`.

    Returns
    -------
    S : (M, N) ndarray
        The S-matrix in the singular value decomposition

    See Also
    --------
    svd : Singular value decomposition of a matrix
    svdvals : Compute singular values of a matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import diagsvd
    >>> vals = np.array([1, 2, 3])  # The array representing the computed svd
    >>> diagsvd(vals, 3, 4)
    array([[1, 0, 0, 0],
           [0, 2, 0, 0],
           [0, 0, 3, 0]])
    >>> diagsvd(vals, 4, 3)
    array([[1, 0, 0],
           [0, 2, 0],
           [0, 0, 3],
           [0, 0, 0]])

    """
    part = np.diag(s)
    typ = part.dtype.char
    MorN = len(s)
    if MorN == M:
        return np.hstack((part, np.zeros((M, N - M), dtype=typ)))
    elif MorN == N:
        return np.r_[part, np.zeros((M - N, N), dtype=typ)]
    else:
        raise ValueError("Length of s must be M or N.")


# Orthonormal decomposition

@_apply_over_batch(('A', 2))
def orth(A, rcond=None):
    """
    Construct an orthonormal basis for the range of A using SVD.

    Parameters
    ----------
    A : (M, N) array_like
        Input array
    rcond : float, optional
        Relative condition number. Singular values ``s`` smaller than
        ``rcond * max(s)`` are considered zero.
        Default: floating point eps * max(M,N).

    Returns
    -------
    Q : (M, K) ndarray
        Orthonormal basis for the range of A.
        K = effective rank of A, as determined by rcond

    See Also
    --------
    svd : Singular value decomposition of a matrix
    null_space : Matrix null space

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import orth
    >>> A = np.array([[2, 0, 0], [0, 5, 0]])  # rank 2 array
    >>> orth(A)
    array([[0., 1.],
           [1., 0.]])
    >>> orth(A.T)
    array([[0., 1.],
           [1., 0.],
           [0., 0.]])

    """
    u, s, vh = svd(A, full_matrices=False)
    M, N = u.shape[0], vh.shape[1]
    if rcond is None:
        rcond = np.finfo(s.dtype).eps * max(M, N)
    tol = np.amax(s, initial=0.) * rcond
    num = np.sum(s > tol, dtype=int)
    Q = u[:, :num]
    return Q


@_apply_over_batch(('A', 2))
def null_space(A, rcond=None, *, overwrite_a=False, check_finite=True,
               lapack_driver='gesdd'):
    """
    Construct an orthonormal basis for the null space of A using SVD.

    Parameters
    ----------
    A : (M, N) array_like
        Input array
    rcond : float, optional
        Relative condition number. Singular values ``s`` smaller than
        ``rcond * max(s)`` are considered zero.
        Default: floating point eps * max(M,N).
    overwrite_a : bool, optional
        Whether to overwrite `a`; may improve performance. Default is False.
        See :ref:`tutorial_linalg_overwrite` for details.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    lapack_driver : {'gesdd', 'gesvd'}, optional
        Whether to use the more efficient divide-and-conquer approach
        (``'gesdd'``) or general rectangular approach (``'gesvd'``)
        to compute the SVD. MATLAB and Octave use the ``'gesvd'`` approach.
        Default is ``'gesdd'``.

    Returns
    -------
    Z : (N, K) ndarray
        Orthonormal basis for the null space of A.
        K = dimension of effective null space, as determined by rcond

    See Also
    --------
    svd : Singular value decomposition of a matrix
    orth : Matrix range

    Examples
    --------
    1-D null space:

    >>> import numpy as np
    >>> from scipy.linalg import null_space
    >>> A = np.array([[1, 1], [1, 1]])
    >>> ns = null_space(A)
    >>> ns * np.copysign(1, ns[0,0])  # Remove the sign ambiguity of the vector
    array([[ 0.70710678],
           [-0.70710678]])

    2-D null space:

    >>> from numpy.random import default_rng
    >>> rng = default_rng()
    >>> B = rng.random((3, 5))
    >>> Z = null_space(B)
    >>> Z.shape
    (5, 2)
    >>> np.allclose(B.dot(Z), 0)
    True

    The basis vectors are orthonormal (up to rounding error):

    >>> Z.T.dot(Z)
    array([[  1.00000000e+00,   6.92087741e-17],
           [  6.92087741e-17,   1.00000000e+00]])

    """
    u, s, vh = svd(A, full_matrices=True, overwrite_a=overwrite_a,
                   check_finite=check_finite, lapack_driver=lapack_driver)
    M, N = u.shape[0], vh.shape[1]
    if rcond is None:
        rcond = np.finfo(s.dtype).eps * max(M, N)
    tol = np.amax(s, initial=0.) * rcond
    num = np.sum(s > tol, dtype=int)
    Q = vh[num:,:].T.conj()
    return Q


@_apply_over_batch(('A', 2), ('B', 2))
def subspace_angles(A, B):
    r"""
    Compute the subspace angles between two matrices.

    Parameters
    ----------
    A : (M, N) array_like
        The first input array.
    B : (M, K) array_like
        The second input array.

    Returns
    -------
    angles : ndarray, shape (min(N, K),)
        The subspace angles between the column spaces of `A` and `B` in
        descending order.

    See Also
    --------
    orth
    svd

    Notes
    -----
    This computes the subspace angles according to the formula
    provided in [1]_. For equivalence with MATLAB and Octave behavior,
    use ``angles[0]``.

    .. versionadded:: 1.0

    References
    ----------
    .. [1] Knyazev A, Argentati M (2002) Principal Angles between Subspaces
           in an A-Based Scalar Product: Algorithms and Perturbation
           Estimates. SIAM J. Sci. Comput. 23:2008-2040.

    Examples
    --------
    A Hadamard matrix, which has orthogonal columns, so we expect that
    the subspace angle to be :math:`\frac{\pi}{2}`:

    >>> import numpy as np
    >>> from scipy.linalg import hadamard, subspace_angles
    >>> rng = np.random.default_rng()
    >>> H = hadamard(4)
    >>> print(H)
    [[ 1  1  1  1]
     [ 1 -1  1 -1]
     [ 1  1 -1 -1]
     [ 1 -1 -1  1]]
    >>> np.rad2deg(subspace_angles(H[:, :2], H[:, 2:]))
    array([ 90.,  90.])

    And the subspace angle of a matrix to itself should be zero:

    >>> subspace_angles(H[:, :2], H[:, :2]) <= 2 * np.finfo(float).eps
    array([ True,  True], dtype=bool)

    The angles between non-orthogonal subspaces are in between these extremes:

    >>> x = rng.standard_normal((4, 3))
    >>> np.rad2deg(subspace_angles(x[:, :2], x[:, [2]]))
    array([ 55.832])  # random
    """
    # Steps here omit the U and V calculation steps from the paper

    # 1. Compute orthonormal bases of column-spaces
    A = _asarray_validated(A, check_finite=True)
    if len(A.shape) != 2:
        raise ValueError(f'expected 2D array, got shape {A.shape}')
    QA = orth(A)
    del A

    B = _asarray_validated(B, check_finite=True)
    if len(B.shape) != 2:
        raise ValueError(f'expected 2D array, got shape {B.shape}')
    if len(B) != len(QA):
        raise ValueError('A and B must have the same number of rows, got '
                         f'{QA.shape[0]} and {B.shape[0]}')
    QB = orth(B)
    del B

    # 2. Compute SVD for cosine
    QA_H_QB = np.dot(QA.T.conj(), QB)
    sigma = svdvals(QA_H_QB)

    # 3. Compute matrix B
    if QA.shape[1] >= QB.shape[1]:
        B = QB - np.dot(QA, QA_H_QB)
    else:
        B = QA - np.dot(QB, QA_H_QB.T.conj())
    del QA, QB, QA_H_QB

    # 4. Compute SVD for sine
    mask = sigma ** 2 >= 0.5
    if mask.any():
        mu_arcsin = np.arcsin(np.clip(svdvals(B, overwrite_a=True), -1., 1.))
    else:
        mu_arcsin = 0.

    # 5. Compute the principal angles
    # with reverse ordering of sigma because smallest sigma belongs to largest
    # angle theta
    theta = np.where(mask, mu_arcsin, np.arccos(np.clip(sigma[::-1], -1., 1.)))
    return theta
