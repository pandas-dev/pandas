"""Abstract linear algebra library.

This module defines a class hierarchy that implements a kind of "lazy"
matrix representation, called the ``LinearOperator``. It can be used to do
linear algebra with extremely large sparse or structured matrices, without
representing those explicitly in memory. Such matrices can be added,
multiplied, transposed, etc.

As a motivating example, suppose you want have a matrix where almost all of
the elements have the value one. The standard sparse matrix representation
skips the storage of zeros, but not ones. By contrast, a LinearOperator is
able to represent such matrices efficiently. First, we need a compact way to
represent an all-ones matrix::

    >>> import numpy as np
    >>> from scipy.sparse.linalg._interface import LinearOperator
    >>> class Ones(LinearOperator):
    ...     def __init__(self, shape):
    ...         super().__init__(dtype=None, shape=shape)
    ...     def _matvec(self, x):
    ...         return np.repeat(x.sum(), self.shape[0])

Instances of this class emulate ``np.ones(shape)``, but using a constant
amount of storage, independent of ``shape``. The ``_matvec`` method specifies
how this linear operator multiplies with (operates on) a vector. We can now
add this operator to a sparse matrix that stores only offsets from one::

    >>> from scipy.sparse.linalg._interface import aslinearoperator
    >>> from scipy.sparse import csr_array
    >>> offsets = csr_array([[1, 0, 2], [0, -1, 0], [0, 0, 3]])
    >>> A = aslinearoperator(offsets) + Ones(offsets.shape)
    >>> A.dot([1, 2, 3])
    array([13,  4, 15])

The result is the same as that given by its dense, explicitly-stored
counterpart::

    >>> (np.ones(A.shape, A.dtype) + offsets.toarray()).dot([1, 2, 3])
    array([13,  4, 15])

Several algorithms in the ``scipy.sparse`` library are able to operate on
``LinearOperator`` instances.
"""

import os
import types
import warnings

import numpy as np

from scipy import sparse
from scipy._external import array_api_extra as xpx
from scipy._lib._array_api import (
    SCIPY_ARRAY_API,
    _asarray,
    array_namespace,
    is_array_api_obj,
    is_pydata_sparse_array,
    np_compat,
    xp_capabilities,
    xp_copy,
    xp_isscalar,
    xp_result_type,
)
from scipy.sparse import issparse
from scipy.sparse._sputils import asmatrix, is_pydata_spmatrix, isintlike, isshape

__all__ = ["LinearOperator", "aslinearoperator"]


@xp_capabilities()
class LinearOperator:
    """Common interface for performing matrix vector products.

    Many iterative methods (e.g. `cg`, `gmres`) do not need to know the
    individual entries of a matrix to solve a linear system ``A@x = b``.
    Such solvers only require the computation of matrix vector
    products, ``A@v``, where ``v`` is a dense vector.  This class serves as
    an abstract interface between iterative solvers and matrix-like
    objects.

    To construct a concrete `LinearOperator`, either pass appropriate
    callables to the constructor of this class, or subclass it.

    A subclass must implement either one of the methods ``_matvec``
    and ``_matmat``, and the attributes/properties ``shape`` (pair of
    integers, optionally with additional batch dimensions at the front)
    and ``dtype`` (may be None). It may call the ``__init__``
    on this class to have these attributes validated. Implementing
    ``_matvec`` automatically implements ``_matmat`` (using a naive
    algorithm) and vice-versa.

    Optionally, a subclass may implement ``_rmatvec`` or ``_adjoint``
    to implement the Hermitian adjoint (conjugate transpose). As with
    ``_matvec`` and ``_matmat``, implementing either ``_rmatvec`` or
    ``_adjoint`` implements the other automatically. Implementing
    ``_adjoint`` is preferable; ``_rmatvec`` is mostly there for
    backwards compatibility.

    The defined operator may have additional "batch" dimensions
    prepended to the core shape, to represent a batch of 2-D operators;
    see :ref:`linalg_batch` for details.

    Parameters
    ----------
    shape : tuple
        Matrix dimensions ``(..., M, N)``,
        where ``...`` represents any additional batch dimensions.
    matvec : callable f(v)
        Applies ``A`` to ``v``, where ``v`` is a dense vector
        with shape ``(..., N)``.
    rmatvec : callable f(v)
        Applies ``A^H`` to ``v``, where ``A^H`` is the conjugate transpose of ``A``,
        and ``v`` is a dense vector of shape ``(..., M)``.
    matmat : callable f(V)
        Returns ``A @ V``, where ``V`` is a dense matrix
        with dimensions ``(..., N, K)``.
    rmatmat : callable f(V)
        Returns ``A^H @ V``, where ``A^H`` is the conjugate transpose of ``A``,
        and where ``V`` is a dense matrix with dimensions ``(..., M, K)``.
    dtype : dtype
        Data type of the matrix or matrices.
    xp : array_namespace, optional
        A namespace compatible with the array API standard for use in array operations.
        Default: ``numpy``.

    Attributes
    ----------
    args : tuple
        For linear operators describing products etc. of other linear
        operators, the operands of the binary operation.
    ndim : int
        Number of dimensions (greater than 2 in the case of batch dimensions).
    T : LinearOperator
        Transpose.
    H : LinearOperator
        Hermitian adjoint.

    Methods
    -------
    matvec
    matmat
    adjoint
    transpose
    rmatvec
    rmatmat
    dot
    rdot
    __mul__
    __matmul__
    __call__
    __add__
    __truediv__
    __rmul__
    __rmatmul__

    See Also
    --------
    aslinearoperator : Construct a `LinearOperator`.

    Notes
    -----
    The user-defined `matvec` function must properly handle the case
    where ``v`` has shape ``(..., N)``.

    It is highly recommended to explicitly specify the `dtype`, otherwise
    it is determined automatically at the cost of a single matvec application
    on ``int8`` zero vector using the promoted `dtype` of the output.
    It is assumed that `matmat`, `rmatvec`, and `rmatmat` would result in
    the same dtype of the output given an ``int8`` input as `matvec`.

    `LinearOperator` instances can also be multiplied, added with each
    other, and raised to integral powers, all lazily: the result of these
    operations
    is always a new, composite `LinearOperator`, that defers linear
    operations to the original operators and combines the results.

    More details regarding how to subclass a `LinearOperator` and several
    examples of concrete `LinearOperator` instances can be found in the
    external project `PyLops <https://pylops.readthedocs.io>`_.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse.linalg import LinearOperator
    >>> def mv(v):
    ...     return np.array([2*v[0], 3*v[1]])
    ...
    >>> A = LinearOperator((2,2), matvec=mv)
    >>> A
    <2x2 _CustomLinearOperator with dtype=int8>
    >>> A.matvec(np.ones(2))
    array([ 2.,  3.])
    >>> A @ np.ones(2)
    array([ 2.,  3.])

    """  # numpydoc ignore=PR01,PR02

    # Necessary for right matmul with numpy arrays.
    __array_ufunc__ = None

    # generic type compatibility with scipy-stubs
    __class_getitem__: classmethod = classmethod(types.GenericAlias)

    ndim: int

    def __new__(cls, *args, **kwargs):
        if cls is LinearOperator:
            # Operate as _CustomLinearOperator factory.
            return super().__new__(_CustomLinearOperator)
        else:
            obj = super().__new__(cls)

            if (
                type(obj)._matvec == LinearOperator._matvec
                and type(obj)._matmat == LinearOperator._matmat
            ):
                warnings.warn(
                    "LinearOperator subclass should implement"
                    " at least one of _matvec and _matmat.",
                    category=RuntimeWarning,
                    stacklevel=2,
                )

            return obj

    def __init__(self, dtype, shape, xp=None):
        """Initialize this LinearOperator.

        To be called by subclasses. ``dtype`` may be None; ``shape`` should
        be convertible to a length >=2 tuple.
        """
        xp = np_compat if xp is None else xp
        if dtype is not None:
            dtype = xp.empty(0, dtype=dtype).dtype

        shape = tuple(shape)
        if len(shape) < 2:
            raise ValueError(f"invalid shape {shape!r} (must be at least 2-d)")
        if not isshape(shape, check_nd=False):
            raise ValueError(f"invalid shape {shape!r}")

        self.dtype = dtype
        self.shape = shape
        self.ndim = len(shape)
        self._xp = xp

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_xp"] = state["_xp"].empty(0)
        return state

    def __setstate__(self, state):
        self._xp = array_namespace(state.pop("_xp"))
        self.__dict__.update(state)

    def _init_dtype(self):
        """Determine the dtype by executing `matvec` on an `int8` test vector.

        In `np.promote_types` hierarchy, the type `int8` is the smallest,
        so we call `matvec` on `int8` and use the promoted dtype of the output
        to set the default `dtype` of the `LinearOperator`.
        We assume that `matmat`, `rmatvec`, and `rmatmat` would result in
        the same dtype of the output given an `int8` input as `matvec`.

        Called from subclasses at the end of the __init__ routine.
        """
        if self.dtype is None:
            xp = self._xp
            v = xp.zeros(self.shape[-1], dtype=xp.int8)
            try:
                matvec_v = xp.asarray(self.matvec(v))
            except (OverflowError, TypeError, RuntimeError):
                self.dtype = xpx.default_dtype(xp=xp, kind="integral")
            else:
                self.dtype = matvec_v.dtype

    def _matmat(self, X, /):
        """Default matrix-matrix multiplication handler.

        If ``self`` is a linear operator of shape ``(..., M, N)``,
        then this method will be called on a shape ``(..., N, K)`` array,
        and should return a shape ``(..., M, K)`` array.

        Falls back to `_matvec`, so defining that will
        define matrix multiplication too (though in a very suboptimal way).
        """
        xp = self._xp

        # NOTE: we can't use `_matvec` directly (for the unbatched case)
        # as we can't assume that user-defined `matvec` functions support batching.
        # TODO: determine whether user-defined `_matvec` supports batching,
        # and use it directly.
        return xp.stack(
            [self._matvec(X[..., :, i]) for i in range(X.shape[-1])], axis=-1
        )

    def _matvec(self, x):
        """Default matrix-vector multiplication handler.

        If ``self`` is a linear operator of shape ``(..., M, N)``,
        then this method will be called on a shape
        ``(..., N)`` array,
        and should return a shape ``(..., M)`` array.

        Falls back to `_matmat`, so defining that
        will define matrix-vector multiplication as well.
        """
        xp = self._xp
        return self._matmat(x[..., xp.newaxis])[..., 0]

    def _shared_matvec(self, x, adjoint: bool = False):
        xp = self._xp

        x = _asarray(x, subok=True, xp=xp)

        *self_broadcast_dims, M, N = self.shape
        inner_dim, outer_dim = (M, N) if adjoint else (N, M)

        # TODO: deprecate `np.matrix` support
        if isinstance(x, np.matrix):
            y = self._rmatvec(x) if adjoint else self._matvec(x)
            if x.shape == (inner_dim, 1):
                y = y.reshape(outer_dim, 1)
            return asmatrix(y)

        x_broadcast_dims: tuple[int, ...] = ()
        row_vector: bool = False
        if column_vector := x.shape == (inner_dim, 1):
            func_name = "`rmatvec`" if adjoint else "`matvec`"
            matmat_func_name = "`rmatmat`" if adjoint else "`matmat`"
            msg = (
                f"Calling {func_name} on 'column vectors' of shape "
                f"`({inner_dim}, 1)` was deprecated in SciPy 1.18.0 and will no "
                f"longer be possible in SciPy 1.20.0. "
                f"Please call {matmat_func_name} instead for identical behaviour."
            )
            warnings.warn(
                msg, FutureWarning, skip_file_prefixes=(os.path.dirname(__file__),)
            )
            x = xp.reshape(x, (inner_dim,))
        elif x.ndim >= 1 and (row_vector := x.shape[-1] == inner_dim):
            x_broadcast_dims = x.shape[:-1]

        if not (row_vector or column_vector):
            msg = (
                f"Dimension mismatch: `x` must have a shape ending in "
                f"`({inner_dim},)`, or shape `({inner_dim}, 1)`. "
                f"Given shape: {x.shape}"
            )
            raise ValueError(msg)

        y = self._rmatvec(x) if adjoint else self._matvec(x)

        broadcasted_dims = xpx.broadcast_shapes(self_broadcast_dims, x_broadcast_dims)
        if row_vector:
            y = xp.reshape(y, (*broadcasted_dims, outer_dim))
        elif column_vector:
            y = xp.reshape(y, (*broadcasted_dims, outer_dim, 1))

        return y

    def matvec(self, x):
        """Matrix-vector multiplication.

        Applies ``A`` to `x`, where ``A`` is an ``M`` x ``N``
        linear operator (or batch of linear operators)
        and `x` is a row vector (or batch of such vectors).

        Parameters
        ----------
        x : {matrix, ndarray}
            An array with shape ``(..., N)`` representing a row vector
            (or batch of row vectors).

            .. versionadded:: 1.18.0
                A ``FutureWarning`` is emitted for column vector input of shape
                ``(N, 1)``, for which an array with shape ``(M, 1)`` is returned.
                `matmat` can be called instead for identical behaviour on such input.

        Returns
        -------
        y : {matrix, ndarray}
            An array with shape ``(..., M)``.

        Notes
        -----
        This method wraps the user-specified ``matvec`` routine or overridden
        ``_matvec`` method to ensure that `y` has the correct shape and type.
        """
        return self._shared_matvec(x)

    def rmatvec(self, x):
        """Adjoint matrix-vector multiplication.

        Applies ``A^H`` to `x`, where ``A`` is an
        ``M`` x ``N`` linear operator (or batch of linear operators)
        and `x` is a row vector (or batch of such vectors).

        Parameters
        ----------
        x : {matrix, ndarray}
            An array with shape ``(..., M)`` representing a row vector
            (or batch of row vectors),
            or an array with shape ``(M, 1)`` representing a column vector.

            .. versionadded:: 1.18.0
                A ``FutureWarning`` is now emitted for column vector input of shape
                ``(M, 1)``, for which an array with shape ``(N, 1)`` is returned.
                `rmatmat` can be called instead for identical behaviour on such input.

        Returns
        -------
        y : {matrix, ndarray}
            An array with shape ``(..., N)``.

        Notes
        -----
        This method wraps the user-specified ``rmatvec`` routine or overridden
        ``_rmatvec`` method to ensure that `y` has the correct shape and type.
        """
        return self._shared_matvec(x, adjoint=True)

    def _rmatvec(self, x):
        """Default implementation of `_rmatvec`.
        Defers to `_rmatmat` or `adjoint`."""
        if type(self)._adjoint == LinearOperator._adjoint:
            # _adjoint not overridden, prevent infinite recursion
            if (
                hasattr(self, "_rmatmat")
                and type(self)._rmatmat != LinearOperator._rmatmat
            ):
                xp = self._xp
                # Try to use _rmatmat as a fallback
                return self._rmatmat(x[..., xp.newaxis])[..., 0]
            raise NotImplementedError
        else:
            return self.H._matvec(x)

    def _shared_matmat(self, X, adjoint: bool = False):
        if not (issparse(X) or is_pydata_spmatrix(X)):
            X = _asarray(X, subok=True, xp=self._xp)

        if X.ndim < 2:
            raise ValueError(f"Expected at least 2-d ndarray or matrix, not {X.ndim}-d")

        if X.shape[-2] != (self.shape[-2] if adjoint else self.shape[-1]):
            raise ValueError(f"Dimension mismatch: {self.shape}, {X.shape}")

        try:
            Y = self._rmatmat(X) if adjoint else self._matmat(X)
        except Exception as e:
            if issparse(X) or is_pydata_spmatrix(X):
                raise TypeError(
                    "Multipliying LinearOperator with a sparse matrix failed."
                    " Try wrapping the matrix with `aslinearoperator` first,"
                    " or ensuring the operator's `matmat` function"
                    " supports sparse input."
                ) from e
            raise

        if isinstance(Y, np.matrix):
            Y = asmatrix(Y)

        return Y

    def matmat(self, X):
        """Matrix-matrix multiplication.

        Performs the operation ``A @ X`` where ``A`` is an ``M`` x ``N``
        linear operator (or batch of linear operators)
        and `X` is a dense ``N`` x ``K`` matrix
        (or batch of dense matrices).

        Parameters
        ----------
        X : {matrix, ndarray}
            An array with shape ``(..., N, K)`` representing the dense matrix
            (or batch of dense matrices).

        Returns
        -------
        Y : {matrix, ndarray}
            An array with shape ``(..., M, K)``.

        Notes
        -----
        This method wraps any user-specified ``matmat`` routine or overridden
        ``_matmat`` method to ensure that `Y` has the correct type.
        """
        return self._shared_matmat(X)

    def rmatmat(self, X):
        """Adjoint matrix-matrix multiplication.

        Performs the operation ``A^H @ X`` where ``A`` is an ``M`` x ``N``
        linear operator (or batch of linear operators)
        and `X` is a dense ``M`` x ``K`` matrix
        (or batch of dense matrices).
        The default implementation defers to the adjoint.

        Parameters
        ----------
        X : {matrix, ndarray}
            An array with shape ``(..., M, K)`` representing the dense matrix
            (or batch of dense matrices).

        Returns
        -------
        Y : {matrix, ndarray}
            An array with shape ``(..., N, K)``.

        Notes
        -----
        This method wraps any user-specified ``rmatmat`` routine or overridden
        ``_rmatmat`` method to ensure that `Y` has the correct type.

        """
        return self._shared_matmat(X, adjoint=True)

    def _rmatmat(self, X, /):
        """Default implementation of `_rmatmat`; defers to `rmatvec` or `adjoint`."""
        if type(self)._adjoint == LinearOperator._adjoint:
            xp = self._xp

            # NOTE: we can't use `_rmatvec` directly as we can't assume that
            # user-defined `rmatvec` functions support batching.
            return xp.stack(
                [self._rmatvec(X[..., :, i]) for i in range(X.shape[-1])], axis=-1
            )
        else:
            return self.H._matmat(X)

    def __call__(self, x):
        """Apply this linear operator.

        Equivalent to `__matmul__`.
        """
        return self @ x

    def __mul__(self, x):
        """Multiplication.

        Used by the ``*`` operator. Equivalent to `dot`.
        """
        return self.dot(x)

    def __truediv__(self, other):
        """Scalar Division.

        Returns a lazily scaled linear operator.
        """
        if not xp_isscalar(other):
            raise ValueError("Can only divide a linear operator by a scalar.")

        return _ScaledLinearOperator(self, 1.0 / other, xp=self._xp)

    def _check_matching_namespace(self, x):
        xp_x = getattr(x, "_xp", None)
        if xp_x is None:
            xp_x = array_namespace(x, self._xp.empty(0), sparse_ok=True)
        if xp_x != self._xp:
            msg = (
                f"Mismatched array namespaces."
                f"Namespace for self is {self._xp}, namespace for x is {xp_x}"
            )
            raise TypeError(msg)

    def dot(self, x):
        """Multi-purpose multiplication method.

        Parameters
        ----------
        x : array_like or `LinearOperator` or scalar
            Array-like input will be interpreted as a 1-D row vector or
            2-D matrix (or batch of matrices)
            depending on its shape. See the Returns section for details.

        Returns
        -------
        Ax : array or `LinearOperator`
            - For `LinearOperator` input, operator composition is performed.

            - For scalar input, a lazily scaled operator is returned.

            - Otherwise, the input is expected to take the form of a dense
              1-D vector or 2-D matrix (or batch of matrices),
              interpreted as follows
              (where ``self`` is an ``M`` by ``N`` linear operator):

              - If `x` has shape ``(N,)``
                it is interpreted as a row vector
                and `matvec` is called.
              - If `x` has shape ``(..., N, K)`` for some
                integer ``K``, it is interpreted as a matrix
                (or batch of matrices if there are batch dimensions)
                and `matmat` is called.

        See Also
        --------
        __mul__ : Equivalent method used by the ``*`` operator.
        __matmul__ :
            Method used by the ``@`` operator which rejects scalar
            input before calling this method.

        Notes
        -----
        To perform matrix-vector multiplication on batches of vectors,
        use `matvec`.

        For clarity, it is recommended to use the `matvec` or
        `matmat` methods directly instead of this method
        when interacting with dense vectors and matrices.

        """
        self._check_matching_namespace(x)
        if isinstance(x, LinearOperator):
            return _ProductLinearOperator(self, x, xp=self._xp)
        elif xp_isscalar(x):
            return _ScaledLinearOperator(self, x, xp=self._xp)
        else:
            if not issparse(x) and not is_pydata_spmatrix(x):
                # scipy.sparse arrays shouldn't be converted to numpy arrays
                # pydata_sparse arrays do not require conversion
                # (and should avoid hitting https://github.com/pydata/sparse/pull/918)
                x = self._xp.asarray(x)

            N = self.shape[-1]

            # treat 1-D input as a vector and >1-D input as a matrix, if the shape fits
            vector = x.shape == (N,)
            matrix = x.ndim >= 2 and x.shape[-2] == N

            if not (vector or matrix):
                msg = (
                    f"Dimension mismatch: array input `x` must have shape `({N},)` "
                    f"or a shape ending in `({N}, K)` for some integer `K`. "
                    f"Given shape: {x.shape}"
                )
                raise ValueError(msg)

            if vector:
                return self.matvec(x)
            elif matrix:
                return self.matmat(x)

    def __matmul__(self, other):
        """Matrix Multiplication.

        Used by the ``@`` operator.
        Rejects scalar input.
        Otherwise, equivalent to `dot`.
        """
        if xp_isscalar(other):
            raise ValueError("Scalar operands are not allowed, use '*' instead")
        return self.__mul__(other)

    def __rmatmul__(self, other):
        """Matrix Multiplication from the right.

        Used by the ``@`` operator from the right.
        Rejects scalar input.
        Otherwise, equivalent to `rdot`.
        """
        if xp_isscalar(other):
            raise ValueError("Scalar operands are not allowed, use '*' instead")
        return self.__rmul__(other)

    def __rmul__(self, x):
        """Multiplication from the right.

        Used by the ``*`` operator from the right. Equivalent to `rdot`.
        """
        return self.rdot(x)

    def rdot(self, x):
        """Multi-purpose multiplication method from the right.

        .. note ::

            This method returns ``x A``.
            To perform adjoint multiplication instead, use one of
            `rmatvec` or `rmatmat`, or take the adjoint first,
            like ``self.H.rdot(x)`` or ``x * self.H``.

        Parameters
        ----------
        x : array_like or `LinearOperator` or scalar
            Array-like input will be interpreted as a 1-D row vector or
            2-D matrix (or batch of matrices)
            depending on its shape. See the Returns section for details.

        Returns
        -------
        xA : array or `LinearOperator`
            - For `LinearOperator` input, operator composition is performed.

            - For scalar input, a lazily scaled operator is returned.

            - Otherwise, the input is expected to take the form of a dense
              1-D vector or 2-D matrix (or batch of matrices),
              interpreted as follows
              (where ``self`` is an ``M`` by ``N`` linear operator):

              - If `x` has shape ``(M,)``
                it is interpreted as a row vector.
              - If `x` has shape ``(..., K, M)`` for some
                integer ``K``, it is interpreted as a matrix
                (or batch of matrices if there are batch dimensions).

        See Also
        --------
        dot : Multi-purpose multiplication method from the left.
        __rmul__ :
            Equivalent method, used by the ``*`` operator from the right.
        __rmatmul__ :
            Method used by the ``@`` operator from the right
            which rejects scalar input before calling this method.
        """
        self._check_matching_namespace(x)
        if isinstance(x, LinearOperator):
            return _ProductLinearOperator(x, self, xp=self._xp)
        elif xp_isscalar(x):
            return _ScaledLinearOperator(self, x, xp=self._xp)
        else:
            if not issparse(x) and not is_pydata_spmatrix(x):
                # scipy.sparse arrays shouldn't be converted to numpy arrays
                # pydata_sparse arrays do not require conversion
                # (and should avoid hitting https://github.com/pydata/sparse/pull/918)
                x = self._xp.asarray(x)

            M = self.shape[-2]

            # treat 1-D input as a vector and >1-D input as a matrix, if the shape fits
            vector = x.shape == (M,)
            matrix = x.ndim >= 2 and x.shape[-1] == M

            if not (vector or matrix):
                msg = (
                    f"Dimension mismatch: `x` must have shape `({M},)` "
                    f"or a shape ending in `(K, {M})` for some integer `K`. "
                    f"Given shape: {x.shape}."
                )
                raise ValueError(msg)

            # We use transpose instead of rmatvec/rmatmat to avoid
            # unnecessary complex conjugation if possible.
            # scipy/scipy#24157
            def mT(x):
                match x.ndim:
                    case 0 | 1:
                        return x
                    case 2:
                        return x.T
                    case _:
                        return self._xp.moveaxis(x, -2, -1)
            if vector:
                return self.T.matvec(mT(x))
            elif matrix:
                return mT(self.T.matmat(mT(x)))

    def __pow__(self, p):
        self._check_matching_namespace(p)
        if xp_isscalar(p):
            return _PowerLinearOperator(self, p, xp=self._xp)
        else:
            return NotImplemented

    def __add__(self, x):
        """Linear operator addition.

        The input must be a `LinearOperator`.
        A lazily summed linear operator is returned.
        """
        self._check_matching_namespace(x)
        if isinstance(x, LinearOperator):
            return _SumLinearOperator(self, x, xp=self._xp)
        else:
            return NotImplemented

    def __neg__(self):
        return _ScaledLinearOperator(self, -1, xp=self._xp)

    def __sub__(self, x):
        return self.__add__(-x)

    def __repr__(self):
        if self.dtype is None:
            dt = "unspecified dtype"
        else:
            dt = "dtype=" + str(self.dtype)

        shape = "x".join(str(dim) for dim in self.shape)
        return f"<{shape} {self.__class__.__name__} with {dt}>"

    def adjoint(self):
        """Hermitian adjoint.

        Returns the Hermitian adjoint of this linear operator,
        also known as the Hermitian
        conjugate or Hermitian transpose. For a complex matrix, the
        Hermitian adjoint is equal to the conjugate transpose.

        Returns
        -------
        `LinearOperator`
            Hermitian adjoint of self.

        See Also
        --------
        :attr:`~scipy.sparse.linalg.LinearOperator.H` : Equivalent attribute.
        """
        return self._adjoint()

    @property
    def H(self):
        """Hermitian adjoint.

        See Also
        --------
        scipy.sparse.linalg.LinearOperator.adjoint : Equivalent method.
        """
        return self.adjoint()

    def transpose(self):
        """Transpose.

        Returns
        -------
        `LinearOperator`
            Transpose of the linear operator.

        See Also
        --------
        :attr:`~scipy.sparse.linalg.LinearOperator.T` : Equivalent attribute.
        """
        return self._transpose()

    @property
    def T(self):
        """Transpose.

        See Also
        --------
        scipy.sparse.linalg.LinearOperator.transpose : Equivalent method.
        """
        return self.transpose()

    def _adjoint(self):
        """Default implementation of `_adjoint`.
        Defers to adjoint functions, e.g. `_rmatvec` for `_matvec`."""
        return _AdjointLinearOperator(self, xp=self._xp)

    def _transpose(self):
        """Default implementation of `_transpose`.
        For `_matvec`, defers to `_rmatvec` + `np.conj`."""
        return _TransposedLinearOperator(self, xp=self._xp)


class _CustomLinearOperator(LinearOperator):
    """Linear operator defined in terms of user-specified operations."""

    def __init__(
        self,
        shape,
        matvec,
        rmatvec=None,
        matmat=None,
        dtype=None,
        rmatmat=None,
        xp=None,
    ):
        super().__init__(dtype, shape, xp)

        self.args = ()

        self.__matvec_impl = matvec
        self.__rmatvec_impl = rmatvec
        self.__rmatmat_impl = rmatmat
        self.__matmat_impl = matmat

        self._init_dtype()

    def _matmat(self, X):
        if self.__matmat_impl is not None:
            return self.__matmat_impl(X)
        else:
            return super()._matmat(X)

    def _matvec(self, x):
        return self.__matvec_impl(x)

    def _rmatvec(self, x):
        func = self.__rmatvec_impl
        if func is None:
            raise NotImplementedError("rmatvec is not defined")
        return self.__rmatvec_impl(x)

    def _rmatmat(self, X):
        if self.__rmatmat_impl is not None:
            return self.__rmatmat_impl(X)
        else:
            return super()._rmatmat(X)

    def _adjoint(self):
        return _CustomLinearOperator(
            shape=(*self.shape[:-2], self.shape[-1], self.shape[-2]),
            matvec=self.__rmatvec_impl,
            rmatvec=self.__matvec_impl,
            matmat=self.__rmatmat_impl,
            rmatmat=self.__matmat_impl,
            dtype=self.dtype,
            xp=self._xp,
        )


class _AdjointLinearOperator(LinearOperator):
    """Adjoint of arbitrary Linear Operator"""

    def __init__(self, A, xp=None):
        shape = (*A.shape[:-2], A.shape[-1], A.shape[-2])
        super().__init__(A.dtype, shape, xp)
        self.A = A
        self.args = (A,)

    def _matvec(self, x):
        return self.A._rmatvec(x)

    def _rmatvec(self, x):
        return self.A._matvec(x)

    def _matmat(self, x):
        return self.A._rmatmat(x)

    def _rmatmat(self, x):
        return self.A._matmat(x)


class _TransposedLinearOperator(LinearOperator):
    """Transposition of arbitrary Linear Operator"""

    def __init__(self, A, xp=None):
        shape = (*A.shape[:-2], A.shape[-1], A.shape[-2])
        super().__init__(A.dtype, shape, xp)
        self.A = A
        self.args = (A,)

    def _matvec(self, x):
        return self._xp.conj(self.A._rmatvec(self._xp.conj(x)))

    def _rmatvec(self, x):
        return self._xp.conj(self.A._matvec(self._xp.conj(x)))

    def _matmat(self, x):
        return self._xp.conj(self.A._rmatmat(self._xp.conj(x)))

    def _rmatmat(self, x):
        return self._xp.conj(self.A._matmat(self._xp.conj(x)))


def _get_dtype(operators, dtypes=None, xp=None):
    """Returns the promoted dtype from input dtypes and operators."""
    xp = np_compat if xp is None else xp
    if dtypes is None:
        dtypes = []
    for obj in operators:
        if obj is not None and hasattr(obj, "dtype"):
            dtypes.append(obj.dtype)
    return xp_result_type(*dtypes, xp=xp)


class _SumLinearOperator(LinearOperator):
    """Representing ``A + B``"""

    def __init__(self, A, B, xp=None):
        if not isinstance(A, LinearOperator) or not isinstance(B, LinearOperator):
            raise ValueError("both operands have to be a LinearOperator")
        *A_broadcast_dims, A_M, A_N = A.shape
        *B_broadcast_dims, B_M, B_N = B.shape
        if (A_M, A_N) != (B_M, B_N):
            raise ValueError(f"cannot add {A} and {B}: shape mismatch")
        broadcasted_dims = xpx.broadcast_shapes(A_broadcast_dims, B_broadcast_dims)
        self.args = (A, B)
        super().__init__(_get_dtype([A, B], xp=xp), (*broadcasted_dims, A_M, A_N), xp)

    def _matvec(self, x):
        return self.args[0].matvec(x) + self.args[1].matvec(x)

    def _rmatvec(self, x):
        return self.args[0].rmatvec(x) + self.args[1].rmatvec(x)

    def _rmatmat(self, x):
        return self.args[0].rmatmat(x) + self.args[1].rmatmat(x)

    def _matmat(self, x):
        return self.args[0].matmat(x) + self.args[1].matmat(x)

    def _adjoint(self):
        A, B = self.args
        return A.H + B.H


class _ProductLinearOperator(LinearOperator):
    """Representing ``A @ B``"""

    def __init__(self, A, B, xp=None):
        if not isinstance(A, LinearOperator) or not isinstance(B, LinearOperator):
            raise ValueError("both operands have to be a LinearOperator")
        *A_broadcast_dims, A_M, A_N = A.shape
        *B_broadcast_dims, B_M, B_N = B.shape
        if A_N != B_M:
            raise ValueError(f"cannot multiply {A} and {B}: shape mismatch")
        broadcasted_dims = np.broadcast_shapes(A_broadcast_dims, B_broadcast_dims)
        super().__init__(_get_dtype([A, B], xp=xp), (*broadcasted_dims, A_M, B_N), xp)
        self.args = (A, B)

    def _matvec(self, x):
        return self.args[0].matvec(self.args[1].matvec(x))

    def _rmatvec(self, x):
        return self.args[1].rmatvec(self.args[0].rmatvec(x))

    def _rmatmat(self, x):
        return self.args[1].rmatmat(self.args[0].rmatmat(x))

    def _matmat(self, x):
        return self.args[0].matmat(self.args[1].matmat(x))

    def _adjoint(self):
        A, B = self.args
        return B.H @ A.H


class _ScaledLinearOperator(LinearOperator):
    """Representing ``alpha * A``"""

    def __init__(self, A, alpha, xp=None):
        if not isinstance(A, LinearOperator):
            raise ValueError("LinearOperator expected as A")
        if not np.isscalar(alpha):
            raise ValueError("scalar expected as alpha")
        if isinstance(A, _ScaledLinearOperator):
            A, alpha_original = A.args
            # Avoid in-place multiplication so that we don't accidentally mutate
            # the original prefactor.
            alpha = alpha * alpha_original

        dtype = _get_dtype([A], [alpha], xp=xp)
        super().__init__(dtype, A.shape, xp)
        self.args = (A, alpha)
        # Note: args[1] is alpha (a scalar), so use `*` below, not `@`

    def _matvec(self, x):
        return self.args[1] * self.args[0].matvec(x)

    def _rmatvec(self, x):
        return self.args[1].conjugate() * self.args[0].rmatvec(x)

    def _rmatmat(self, x):
        return self.args[1].conjugate() * self.args[0].rmatmat(x)

    def _matmat(self, x):
        return self.args[1] * self.args[0].matmat(x)

    def _adjoint(self):
        A, alpha = self.args
        return A.H * alpha.conjugate()


class _PowerLinearOperator(LinearOperator):
    """Representing ``A ** p``"""

    def __init__(self, A, p, xp=None):
        if not isinstance(A, LinearOperator):
            raise ValueError("LinearOperator expected as A")
        if A.shape[-2] != A.shape[-1]:
            msg = f"square core-dimensions of LinearOperator expected, got {A!r}"
            raise ValueError(msg)
        if not isintlike(p) or p < 0:
            raise ValueError("non-negative integer expected as p")

        super().__init__(_get_dtype([A], xp=xp), A.shape, xp)
        self.args = (A, p)

    def _power(self, fun, x):
        res = xp_copy(x)
        for i in range(self.args[1]):
            res = fun(res)
        return res

    def _matvec(self, x):
        return self._power(self.args[0].matvec, x)

    def _rmatvec(self, x):
        return self._power(self.args[0].rmatvec, x)

    def _rmatmat(self, x):
        return self._power(self.args[0].rmatmat, x)

    def _matmat(self, x):
        return self._power(self.args[0].matmat, x)

    def _adjoint(self):
        A, p = self.args
        return A.H**p


class MatrixLinearOperator(LinearOperator):
    """Operator defined by a matrix `A` which implements ``@``."""

    def __init__(self, A, xp=None):
        super().__init__(A.dtype, A.shape, xp)
        self.A = A
        self.__adj = None
        self.args = (A,)

    def _matmat(self, X):
        return self.A @ X

    def _adjoint(self):
        if self.__adj is None:
            self.__adj = _AdjointMatrixOperator(self.A, xp=self._xp)
        return self.__adj


class _AdjointMatrixOperator(MatrixLinearOperator):
    """Representing ``A.H``, for `MatrixLinearOperator` `A`."""

    def __init__(self, A, xp=None):
        xp = np_compat if xp is None else xp
        if A.ndim > 2:
            if issparse(A):
                A_T = sparse.swapaxes(A, -1, -2)
            else:
                A_T = A.mT
        else:
            A_T = A.T
        super().__init__(xp.conj(A_T), xp=xp)
        self.args = (A,)  # override to ensure `self.args[0]` is accurate

    def _adjoint(self):
        return MatrixLinearOperator(self.args[0], xp=self._xp)


class IdentityOperator(LinearOperator):
    def __init__(self, shape, dtype=None, xp=None):
        super().__init__(dtype, shape, xp)

    def _matvec(self, x):
        return x

    def _rmatvec(self, x):
        return x

    def _rmatmat(self, x):
        return x

    def _matmat(self, x):
        return x

    def _adjoint(self):
        return self


@xp_capabilities()
def aslinearoperator(A):
    """Return `A` as a `LinearOperator`.

    See the `LinearOperator` documentation for additional information.

    Parameters
    ----------
    A : object
        Object to convert to a `LinearOperator`. May be any one of the following types:

        - `numpy.ndarray`
        - `numpy.matrix`
        - `scipy.sparse` array
          (e.g. `~scipy.sparse.csr_array`, `~scipy.sparse.lil_array`, etc.)
        - `LinearOperator`
        - An object with ``.shape`` and ``.matvec`` attributes

    Returns
    -------
    B : LinearOperator
        A `LinearOperator` corresponding with `A`

    Notes
    -----
    If `A` has no ``.dtype`` attribute, the data type is determined by calling
    :func:`LinearOperator.matvec()` - set the ``.dtype`` attribute to prevent this
    call upon the linear operator creation.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse.linalg import aslinearoperator
    >>> M = np.array([[1,2,3],[4,5,6]], dtype=np.int32)
    >>> aslinearoperator(M)
    <2x3 MatrixLinearOperator with dtype=int32>
    """
    if isinstance(A, LinearOperator):
        return A

    # https://github.com/data-apis/array-api-compat/issues/388
    if is_array_api_obj(A) and not isinstance(A, np.matrix):
        xp = array_namespace(A)
        if is_pydata_sparse_array(A) and not SCIPY_ARRAY_API:
            pass  # skip `atleast_nd`, since `xp=np` does not work with `sparse`
        else:
            A = xpx.atleast_nd(A, ndim=2, xp=xp)
        return MatrixLinearOperator(A, xp=xp)

    if issparse(A):
        return MatrixLinearOperator(A)

    if isinstance(A, np.matrix):
        A = np.atleast_2d(np.asarray(A))
        return MatrixLinearOperator(A)

    if hasattr(A, "shape") and hasattr(A, "matvec"):
        rmatvec = None
        rmatmat = None
        dtype = None

        if hasattr(A, "rmatvec"):
            rmatvec = A.rmatvec
        if hasattr(A, "rmatmat"):
            rmatmat = A.rmatmat
        if hasattr(A, "dtype"):
            dtype = A.dtype
        return LinearOperator(
            A.shape, A.matvec, rmatvec=rmatvec, rmatmat=rmatmat, dtype=dtype
        )

    raise TypeError("type not understood")
