import collections
import numbers
import numpy as np

from ._input_validation import _nonneg_int_or_fail

from ._special_ufuncs import (legendre_p, assoc_legendre_p,
                              sph_legendre_p, sph_harm_y)
from ._gufuncs import (legendre_p_all, assoc_legendre_p_all,
                       sph_legendre_p_all, sph_harm_y_all)

__all__ = [
    "assoc_legendre_p",
    "assoc_legendre_p_all",
    "legendre_p",
    "legendre_p_all",
    "sph_harm_y",
    "sph_harm_y_all",
    "sph_legendre_p",
    "sph_legendre_p_all",
]


class MultiUFunc:
    def __init__(self, ufunc_or_ufuncs, name=None, doc=None, *,
                 force_complex_output=False, **default_kwargs):
        if not isinstance(ufunc_or_ufuncs, np.ufunc):
            if isinstance(ufunc_or_ufuncs, collections.abc.Mapping):
                ufuncs_iter = ufunc_or_ufuncs.values()
            elif isinstance(ufunc_or_ufuncs, collections.abc.Iterable):
                ufuncs_iter = ufunc_or_ufuncs
            else:
                raise ValueError("ufunc_or_ufuncs should be a ufunc or a"
                                 " ufunc collection")

            # Perform input validation to ensure all ufuncs in ufuncs are
            # actually ufuncs and all take the same input types.
            seen_input_types = set()
            for ufunc in ufuncs_iter:
                if not isinstance(ufunc, np.ufunc):
                    raise ValueError("All ufuncs must have type `numpy.ufunc`."
                                     f" Received {ufunc_or_ufuncs}")
                seen_input_types.add(frozenset(x.split("->")[0] for x in ufunc.types))
            if len(seen_input_types) > 1:
                raise ValueError("All ufuncs must take the same input types.")

        self.__name__ = name
        self._ufunc_or_ufuncs = ufunc_or_ufuncs
        self.__doc = doc
        self.__force_complex_output = force_complex_output
        self._default_kwargs = default_kwargs
        self._resolve_out_shapes = None
        self._finalize_out = None
        self._key = None
        self._ufunc_default_args = lambda *args, **kwargs: ()
        self._ufunc_default_kwargs = lambda *args, **kwargs: {}
        self.__text_signature__ = doc.split("\n\n")[0] if doc else None

    @property
    def __doc__(self):  # pyrefly:ignore[bad-override]
        return self.__doc

    def _override_key(self, func):
        """Set `key` method by decorating a function.
        """
        self._key = func

    def _override_ufunc_default_args(self, func):
        self._ufunc_default_args = func

    def _override_ufunc_default_kwargs(self, func):
        self._ufunc_default_kwargs = func

    def _override_resolve_out_shapes(self, func):
        """Set `resolve_out_shapes` method by decorating a function."""
        if func.__doc__ is None:
            func.__doc__ = \
                """Resolve to output shapes based on relevant inputs."""
        func.__name__ = "resolve_out_shapes"
        self._resolve_out_shapes = func

    def _override_finalize_out(self, func):
        self._finalize_out = func

    def _resolve_ufunc(self, **kwargs):
        """Resolve to a ufunc based on keyword arguments."""

        if isinstance(self._ufunc_or_ufuncs, np.ufunc):
            return self._ufunc_or_ufuncs

        ufunc_key = self._key(**kwargs)
        return self._ufunc_or_ufuncs[ufunc_key]

    def __call__(self, *args, **kwargs):
        kwargs = self._default_kwargs | kwargs

        args += self._ufunc_default_args(**kwargs)

        ufunc = self._resolve_ufunc(**kwargs)

        # array arguments to be passed to the ufunc
        ufunc_args = [np.asarray(arg) for arg in args[-ufunc.nin:]]

        ufunc_kwargs = self._ufunc_default_kwargs(**kwargs)

        if (self._resolve_out_shapes is not None):
            ufunc_arg_shapes = tuple(np.shape(ufunc_arg) for ufunc_arg in ufunc_args)
            ufunc_out_shapes = self._resolve_out_shapes(*args[:-ufunc.nin],
                                                        *ufunc_arg_shapes, ufunc.nout,
                                                        **kwargs)

            ufunc_arg_dtypes = tuple(ufunc_arg.dtype if hasattr(ufunc_arg, 'dtype')
                                     else np.dtype(type(ufunc_arg))
                                     for ufunc_arg in ufunc_args)

            if hasattr(ufunc, 'resolve_dtypes'):
                ufunc_dtypes = ufunc_arg_dtypes + ufunc.nout * (None,)
                ufunc_dtypes = ufunc.resolve_dtypes(ufunc_dtypes)
                ufunc_out_dtypes = ufunc_dtypes[-ufunc.nout:]
            else:
                ufunc_out_dtype = np.result_type(*ufunc_arg_dtypes)
                if (not np.issubdtype(ufunc_out_dtype, np.inexact)):
                    ufunc_out_dtype = np.float64

                ufunc_out_dtypes = ufunc.nout * (ufunc_out_dtype,)

            if self.__force_complex_output:
                ufunc_out_dtypes = tuple(np.result_type(1j, ufunc_out_dtype)
                                         for ufunc_out_dtype in ufunc_out_dtypes)

            out = tuple(np.empty(ufunc_out_shape, dtype=ufunc_out_dtype)
                        for ufunc_out_shape, ufunc_out_dtype
                        in zip(ufunc_out_shapes, ufunc_out_dtypes))

            ufunc_kwargs['out'] = out

        out = ufunc(*ufunc_args, **ufunc_kwargs)
        if (self._finalize_out is not None):
            out = self._finalize_out(out)

        return out


sph_legendre_p = MultiUFunc(
    sph_legendre_p,
    "sph_legendre_p",
    r"""sph_legendre_p(n, m, theta, *, diff_n=0)

    Spherical Legendre polynomial of the first kind.

    Parameters
    ----------
    n : array_like of ints
        Degree of the spherical Legendre polynomial. Must have ``n >= 0``.
    m : array_like of ints
        Order of the spherical Legendre polynomial.
    theta : array_like
        Input value.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order `diff_n`. Default is 0.

    Returns
    -------
    ndarray or tuple of ndarray
        Spherical Legendre polynomial with ``diff_n`` derivatives.

    Notes
    -----
    The spherical counterpart of an (unnormalized) associated Legendre polynomial has
    the additional factor

    .. math::

        \sqrt{\frac{(2 n + 1) (n - m)!}{4 \pi (n + m)!}}

    It is the same as the spherical harmonic :math:`Y_{n}^{m}(\theta, \phi)`
    with :math:`\phi = 0`.
    """, diff_n=0
)


@sph_legendre_p._override_key
def _(diff_n):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return diff_n


@sph_legendre_p._override_finalize_out
def _(out):
    return np.moveaxis(out, -1, 0)


sph_legendre_p_all = MultiUFunc(
    sph_legendre_p_all,
    "sph_legendre_p_all",
    """sph_legendre_p_all(n, m, theta, *, diff_n=0)

    All spherical Legendre polynomials of the first kind up to the
    specified degree `n`, order `m`, and all derivatives up
    to order `diff_n`.

    Parameters
    ----------
    n : int
        Degree of the spherical Legendre polynomials. Must have ``n >= 0``.
    m : int
        Order of the spherical Legendre polynomials.
    theta : array_like
        Input value.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order `diff_n`. Default is 0.

    Returns
    -------
    ndarray
        Output shape is ``(diff_n + 1, n + 1, 2 * m + 1, ...)``. The entry at
        ``(i, j, k)`` corresponds to the ``i``-th derivative, degree ``j``, and
        order ``k`` for all ``0 <= i <= diff_n``, ``0 <= j <= n``, and
        ``-m <= k <= m``.

    See Also
    --------
    sph_legendre_p
    """, diff_n=0
)


@sph_legendre_p_all._override_key
def _(diff_n):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return diff_n


@sph_legendre_p_all._override_ufunc_default_kwargs
def _(diff_n):
    return {'axes': [()] + [(0, 1, -1)]}


@sph_legendre_p_all._override_resolve_out_shapes
def _(n, m, theta_shape, nout, diff_n):
    if not isinstance(n, numbers.Integral) or (n < 0):
        raise ValueError("n must be a non-negative integer.")

    return ((n + 1, 2 * abs(m) + 1) + theta_shape + (diff_n + 1,),)


@sph_legendre_p_all._override_finalize_out
def _(out):
    return np.moveaxis(out, -1, 0)


assoc_legendre_p = MultiUFunc(
    assoc_legendre_p,
    "assoc_legendre_p",
    r"""assoc_legendre_p(n, m, z, *, branch_cut=2, norm=False, diff_n=0)

    Associated Legendre polynomial of the first kind.

    Defined as

    .. math::

        P_n^m(z) = (-1)^m (1 - z^2)^{m/2}
            \frac{d^m}{dz^m} P_n(z)

    where :math:`P_n` is the Legendre polynomial.

    This definition includes the Condon-Shortley phase :math:`(-1)^m`.

    Parameters
    ----------
    n : array_like of ints
        Degree of the associated Legendre polynomial. Must have ``n >= 0``.
    m : array_like of ints
        Order of the associated Legendre polynomial.
    z : array_like
        Input value.
    branch_cut : array_like of ints, optional
        Selects branch cut. Must be 2 (default) or 3.
        2: cut on the real axis ``|z| > 1``
        3: cut on the real axis ``-1 < z < 1``
    norm : bool, optional
        If ``True``, compute the normalized associated Legendre polynomial.
        Default is ``False``.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``. Default is 0.

    Returns
    -------
    ndarray or tuple of ndarray
        Associated Legendre polynomial with ``diff_n`` derivatives.

    Notes
    -----
    The normalized counterpart of an (unnormalized) associated Legendre
    polynomial has the additional factor

    .. math::

        \sqrt{\frac{(2 n + 1) (n - m)!}{2 (n + m)!}}

    Examples
    --------
    Evaluate :math:`P_2^1(z)` on :math:`[-1, 1]`:

    >>> import numpy as np
    >>> from scipy.special import assoc_legendre_p
    >>> z = np.linspace(-1, 1, 11)
    >>> expected = -3 * z * np.sqrt(1 - z**2)
    >>> np.allclose(assoc_legendre_p(2, 1, z), expected)
    True

    Compute the normalized associated Legendre polynomial:

    >>> scale = np.sqrt(5/12)
    >>> np.allclose(assoc_legendre_p(2, 1, z, norm=True), scale * expected)
    True
    """, branch_cut=2, norm=False, diff_n=0
)


@assoc_legendre_p._override_key
def _(branch_cut, norm, diff_n):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return norm, diff_n


@assoc_legendre_p._override_ufunc_default_args
def _(branch_cut, norm, diff_n):
    return branch_cut,


@assoc_legendre_p._override_finalize_out
def _(out):
    return np.moveaxis(out, -1, 0)


assoc_legendre_p_all = MultiUFunc(
    assoc_legendre_p_all,
    "assoc_legendre_p_all",
    """assoc_legendre_p_all(n, m, z, *, branch_cut=2, norm=False, diff_n=0)

    All associated Legendre polynomials of the first kind up to the
    specified degree `n`, order `m`, and all derivatives up
    to order `diff_n`.

    Parameters
    ----------
    n : int
        Degree of the associated Legendre polynomials. Must have ``n >= 0``.
    m : int
        Order of the associated Legendre polynomials.
    z : array_like
        Input value.
    branch_cut : array_like of ints, optional
        Selects branch cut. Must be 2 (default) or 3.
        2: cut on the real axis ``|z| > 1``
        3: cut on the real axis ``-1 < z < 1``
    norm : bool, optional
        If ``True``, compute the normalized associated Legendre polynomials.
        Default is ``False``.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``. Default is 0.

    Returns
    -------
    ndarray
        Output shape is ``(diff_n + 1, n + 1, 2 * m + 1, ...)``. The entry at
        ``(i, j, k)`` corresponds to the ``i``-th derivative, degree ``j``, and
        order ``k`` for all ``0 <= i <= diff_n``, ``0 <= j <= n``, and
        ``-m <= k <= m``.

    See Also
    --------
    assoc_legendre_p
    """, branch_cut=2, norm=False, diff_n=0
)


@assoc_legendre_p_all._override_key
def _(branch_cut, norm, diff_n):
    if not ((isinstance(diff_n, numbers.Integral))
            and diff_n >= 0):
        raise ValueError(
            f"diff_n must be a non-negative integer, received: {diff_n}."
        )
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return norm, diff_n


@assoc_legendre_p_all._override_ufunc_default_args
def _(branch_cut, norm, diff_n):
    return branch_cut,


@assoc_legendre_p_all._override_ufunc_default_kwargs
def _(branch_cut, norm, diff_n):
    return {'axes': [(), ()] + [(0, 1, -1)]}


@assoc_legendre_p_all._override_resolve_out_shapes
def _(n, m, z_shape, branch_cut_shape, nout, **kwargs):
    diff_n = kwargs['diff_n']

    if not isinstance(n, numbers.Integral) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if not isinstance(m, numbers.Integral) or (m < 0):
        raise ValueError("m must be a non-negative integer.")

    return ((n + 1, 2 * abs(m) + 1) +
        np.broadcast_shapes(z_shape, branch_cut_shape) + (diff_n + 1,),)


@assoc_legendre_p_all._override_finalize_out
def _(out):
    return np.moveaxis(out, -1, 0)


legendre_p = MultiUFunc(
    legendre_p,
    "legendre_p",
    """legendre_p(n, z, *, diff_n=0)

    Legendre polynomial of the first kind.

    Parameters
    ----------
    n : array_like of ints
        Degree of the Legendre polynomial. Must have ``n >= 0``.
    z : array_like
        Input value.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``. Default is 0.

    Returns
    -------
    p : ndarray or tuple of ndarray
        Legendre polynomial with ``diff_n`` derivatives.

    See Also
    --------
    legendre

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.html

    Examples
    --------
    Evaluate the Legendre polynomial :math:`P_3` at :math:`z = 0.5`:

    >>> import numpy as np
    >>> from scipy.special import legendre_p
    >>> np.allclose(legendre_p(3, 0.5), -0.4375)
    True

    Compute the value and first derivative with respect to ``z``:

    >>> p, dp = legendre_p(3, 0.5, diff_n=1)
    >>> np.allclose([p, dp], [-0.4375, 0.375])
    True
    """, diff_n=0
)


@legendre_p._override_key
def _(diff_n):
    if (not isinstance(diff_n, numbers.Integral)) or (diff_n < 0):
        raise ValueError(
            f"diff_n must be a non-negative integer, received: {diff_n}."
        )
    if not 0 <= diff_n <= 2:
        raise NotImplementedError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return diff_n


@legendre_p._override_finalize_out
def _(out):
    return np.moveaxis(out, -1, 0)


legendre_p_all = MultiUFunc(
    legendre_p_all,
    "legendre_p_all",
    """legendre_p_all(n, z, *, diff_n=0)

    All Legendre polynomials of the first kind up to the specified degree
    `n` and all derivatives up to order `diff_n`.

    Parameters
    ----------
    n : int
        Degree of the Legendre polynomials. Must have ``n >= 0``.
    z : array_like
        Input value.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order ``diff_n``. Default is 0.

    Returns
    -------
    ndarray
        Output shape is ``(diff_n + 1, n + 1, ...)``. The entry at ``(i, j)``
        corresponds to the ``i``-th derivative and degree ``j`` for all
        ``0 <= i <= diff_n`` and ``0 <= j <= n``.

    See Also
    --------
    legendre_p
    """, diff_n=0
)


@legendre_p_all._override_key
def _(diff_n):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return diff_n


@legendre_p_all._override_ufunc_default_kwargs
def _(diff_n):
    return {'axes': [(), (0, -1)]}


@legendre_p_all._override_resolve_out_shapes
def _(n, z_shape, nout, diff_n):
    n = _nonneg_int_or_fail(n, 'n', strict=False)

    return nout * ((n + 1,) + z_shape + (diff_n + 1,),)


@legendre_p_all._override_finalize_out
def _(out):
    return np.moveaxis(out, -1, 0)


sph_harm_y = MultiUFunc(
    sph_harm_y,
    "sph_harm_y",
    r"""sph_harm_y(n, m, theta, phi, *, diff_n=0)

    Spherical harmonics.

    They are defined as

    .. math::

        Y_n^m(\theta,\phi) = \sqrt{\frac{2 n + 1}{4 \pi} \frac{(n - m)!}{(n + m)!}}
            P_n^m(\cos(\theta)) e^{i m \phi}

    where :math:`P_n^m` are the (unnormalized) associated Legendre polynomials.

    Parameters
    ----------
    n : array_like of ints
        Degree of the harmonic. Must have ``n >= 0``. This is
        often denoted by ``l`` (lower case L) in descriptions of
        spherical harmonics.
    m : array_like of ints
        Order of the harmonic.
    theta : array_like
        Polar (colatitudinal) coordinate; must be in ``[0, pi]``.
    phi : array_like
        Azimuthal (longitudinal) coordinate; must be in ``[0, 2*pi]``.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order `diff_n`. Default is 0.

    Returns
    -------
    ndarray or tuple of ndarray
       Spherical harmonics with `diff_n` derivatives.

    Notes
    -----
    There are different conventions for the meanings of the input
    arguments ``theta`` and ``phi``. In SciPy ``theta`` is the
    polar angle and ``phi`` is the azimuthal angle. It is common to
    see the opposite convention, that is, ``theta`` as the azimuthal angle
    and ``phi`` as the polar angle.

    Note that SciPy's spherical harmonics include the Condon-Shortley
    phase [2]_ because it is part of `sph_legendre_p`.

    With SciPy's conventions, the first several spherical harmonics
    are

    .. math::

        Y_0^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{1}{\pi}} \\
        Y_1^{-1}(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                    e^{-i\phi} \sin(\theta) \\
        Y_1^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{\pi}}
                                 \cos(\theta) \\
        Y_1^1(\theta, \phi) &= -\frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                 e^{i\phi} \sin(\theta).

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 14.30.
           https://dlmf.nist.gov/14.30
    .. [2] https://en.wikipedia.org/wiki/Spherical_harmonics#Condon.E2.80.93Shortley_phase
    """, force_complex_output=True, diff_n=0
)


@sph_harm_y._override_key
def _(diff_n):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 0, 1, and 2,"
            f" received: {diff_n}."
        )
    return diff_n


@sph_harm_y._override_finalize_out
def _(out):
    if (out.shape[-1] == 1):
        return out[..., 0, 0]

    if (out.shape[-1] == 2):
        return out[..., 0, 0], out[..., [1, 0], [0, 1]]

    if (out.shape[-1] == 3):
        return (out[..., 0, 0], out[..., [1, 0], [0, 1]],
            out[..., [[2, 1], [1, 0]], [[0, 1], [1, 2]]])


sph_harm_y_all = MultiUFunc(
    sph_harm_y_all,
    "sph_harm_y_all",
    """sph_harm_y_all(n, m, theta, phi, *, diff_n=0)

    All spherical harmonics up to the specified degree `n`, order `m`,
    and all derivatives up to order `diff_n`.

    Parameters
    ----------
    n : int
        Degree of the harmonics. Must have ``n >= 0``. This is
        often denoted by ``l`` (lower case L) in descriptions of
        spherical harmonics.
    m : int
        Order of the harmonics.
    theta : array_like
        Polar (colatitudinal) coordinate; must be in ``[0, pi]``.
    phi : array_like
        Azimuthal (longitudinal) coordinate; must be in ``[0, 2*pi]``.
    diff_n : int, optional
        A non-negative integer. Compute and return all derivatives up
        to order `diff_n`. Default is 0.

    Returns
    -------
    ndarray or tuple of ndarray
        Returns a tuple of length ``diff_n + 1`` (if ``diff_n > 0``). The first
        entry corresponds to the spherical harmonics, the second entry
        (if ``diff_n >= 1``) to the gradient, and the third entry
        (if ``diff_n >= 2``)  to the Hessian matrix. Each entry is an array of
        shape ``(n + 1, 2 * m + 1, ...)``, where the entry at ``(i, j)``
        corresponds to degree ``i`` and order ``j`` for all ``0 <= i <= n``
        and ``-m <= j <= m``.

    See Also
    --------
    sph_harm_y
    """, force_complex_output=True, diff_n=0
)


@sph_harm_y_all._override_key
def _(diff_n):
    diff_n = _nonneg_int_or_fail(diff_n, "diff_n", strict=False)
    if not 0 <= diff_n <= 2:
        raise ValueError(
            "diff_n is currently only implemented for orders 2,"
            f" received: {diff_n}."
        )
    return diff_n


@sph_harm_y_all._override_ufunc_default_kwargs
def _(diff_n):
    return {'axes': [(), ()] + [(0, 1, -2, -1)]}


@sph_harm_y_all._override_resolve_out_shapes
def _(n, m, theta_shape, phi_shape, nout, **kwargs):
    diff_n = kwargs['diff_n']

    if not isinstance(n, numbers.Integral) or (n < 0):
        raise ValueError("n must be a non-negative integer.")

    return ((n + 1, 2 * abs(m) + 1) + np.broadcast_shapes(theta_shape, phi_shape) +
        (diff_n + 1, diff_n + 1),)


@sph_harm_y_all._override_finalize_out
def _(out):
    if (out.shape[-1] == 1):
        return out[..., 0, 0]

    if (out.shape[-1] == 2):
        return out[..., 0, 0], out[..., [1, 0], [0, 1]]

    if (out.shape[-1] == 3):
        return (out[..., 0, 0], out[..., [1, 0], [0, 1]],
            out[..., [[2, 1], [1, 0]], [[0, 1], [1, 2]]])
