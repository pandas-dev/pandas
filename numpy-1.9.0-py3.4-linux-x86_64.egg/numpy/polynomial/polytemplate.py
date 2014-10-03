"""
Template for the Chebyshev and Polynomial classes.

This module houses a Python string module Template object (see, e.g.,
http://docs.python.org/library/string.html#template-strings) used by
the `polynomial` and `chebyshev` modules to implement their respective
`Polynomial` and `Chebyshev` classes.  It provides a mechanism for easily
creating additional specific polynomial classes (e.g., Legendre, Jacobi,
etc.) in the future, such that all these classes will have a common API.

"""
from __future__ import division, absolute_import, print_function

import string
import sys
import warnings
from number import Number

from numpy import ModuleDeprecationWarning

warnings.warn("The polytemplate module will be removed in Numpy 1.10.0.",
              ModuleDeprecationWarning)

polytemplate = string.Template('''
from __future__ import division, absolute_import, print_function
import numpy as np
import warnings
from . import polyutils as pu

class $name(pu.PolyBase) :
    """A $name series class.

    $name instances provide the standard Python numerical methods '+',
    '-', '*', '//', '%', 'divmod', '**', and '()' as well as the listed
    methods.

    Parameters
    ----------
    coef : array_like
        $name coefficients, in increasing order.  For example,
        ``(1, 2, 3)`` implies ``P_0 + 2P_1 + 3P_2`` where the
        ``P_i`` are a graded polynomial basis.
    domain : (2,) array_like, optional
        Domain to use. The interval ``[domain[0], domain[1]]`` is mapped to
        the interval ``[window[0], window[1]]`` by shifting and scaling.
        The default value is $domain.
    window : (2,) array_like, optional
        Window, see ``domain`` for its use. The default value is $domain.
        .. versionadded:: 1.6.0

    Attributes
    ----------
    coef : (N,) ndarray
        $name coefficients, from low to high.
    domain : (2,) ndarray
        Domain that is mapped to ``window``.
    window : (2,) ndarray
        Window that ``domain`` is mapped to.

    Class Attributes
    ----------------
    maxpower : int
        Maximum power allowed, i.e., the largest number ``n`` such that
        ``p(x)**n`` is allowed. This is to limit runaway polynomial size.
    domain : (2,) ndarray
        Default domain of the class.
    window : (2,) ndarray
        Default window of the class.

    Notes
    -----
    It is important to specify the domain in many cases, for instance in
    fitting data, because many of the important properties of the
    polynomial basis only hold in a specified interval and consequently
    the data must be mapped into that interval in order to benefit.

    Examples
    --------

    """
    # Limit runaway size. T_n^m has degree n*2^m
    maxpower = 16
    # Default domain
    domain = np.array($domain)
    # Default window
    window = np.array($domain)
    # Don't let participate in array operations. Value doesn't matter.
    __array_priority__ = 1000
    # Not hashable
    __hash__ = None

    def has_samecoef(self, other):
        """Check if coefficients match.

        Parameters
        ----------
        other : class instance
            The other class must have the ``coef`` attribute.

        Returns
        -------
        bool : boolean
            True if the coefficients are the same, False otherwise.

        Notes
        -----
        .. versionadded:: 1.6.0

        """
        if len(self.coef) != len(other.coef):
            return False
        elif not np.all(self.coef == other.coef):
            return False
        else:
            return True

    def has_samedomain(self, other):
        """Check if domains match.

        Parameters
        ----------
        other : class instance
            The other class must have the ``domain`` attribute.

        Returns
        -------
        bool : boolean
            True if the domains are the same, False otherwise.

        Notes
        -----
        .. versionadded:: 1.6.0

        """
        return np.all(self.domain == other.domain)

    def has_samewindow(self, other):
        """Check if windows match.

        Parameters
        ----------
        other : class instance
            The other class must have the ``window`` attribute.

        Returns
        -------
        bool : boolean
            True if the windows are the same, False otherwise.

        Notes
        -----
        .. versionadded:: 1.6.0

        """
        return np.all(self.window == other.window)

    def has_sametype(self, other):
        """Check if types match.

        Parameters
        ----------
        other : object
            Class instance.

        Returns
        -------
        bool : boolean
            True if other is same class as self

        Notes
        -----
        .. versionadded:: 1.7.0

        """
        return isinstance(other, self.__class__)

    def __init__(self, coef, domain=$domain, window=$domain) :
        [coef, dom, win] = pu.as_series([coef, domain, window], trim=False)
        if len(dom) != 2 :
            raise ValueError("Domain has wrong number of elements.")
        if len(win) != 2 :
            raise ValueError("Window has wrong number of elements.")
        self.coef = coef
        self.domain = dom
        self.window = win

    def __repr__(self):
        format = "%s(%s, %s, %s)"
        coef = repr(self.coef)[6:-1]
        domain = repr(self.domain)[6:-1]
        window = repr(self.window)[6:-1]
        return format % ('$name', coef, domain, window)

    def __str__(self) :
        format = "%s(%s)"
        coef = str(self.coef)
        return format % ('$nick', coef)

    # Pickle and copy

    def __getstate__(self) :
        ret = self.__dict__.copy()
        ret['coef'] = self.coef.copy()
        ret['domain'] = self.domain.copy()
        ret['window'] = self.window.copy()
        return ret

    def __setstate__(self, dict) :
        self.__dict__ = dict

    # Call

    def __call__(self, arg) :
        off, scl = pu.mapparms(self.domain, self.window)
        arg = off + scl*arg
        return ${nick}val(arg, self.coef)

    def __iter__(self) :
        return iter(self.coef)

    def __len__(self) :
        return len(self.coef)

    # Numeric properties.

    def __neg__(self) :
        return self.__class__(-self.coef, self.domain, self.window)

    def __pos__(self) :
        return self

    def __add__(self, other) :
        """Returns sum"""
        if isinstance(other, pu.PolyBase):
            if not self.has_sametype(other):
                raise TypeError("Polynomial types differ")
            elif not self.has_samedomain(other):
                raise TypeError("Domains differ")
            elif not self.has_samewindow(other):
                raise TypeError("Windows differ")
            else:
                coef = ${nick}add(self.coef, other.coef)
        else :
            try :
                coef = ${nick}add(self.coef, other)
            except :
                return NotImplemented
        return self.__class__(coef, self.domain, self.window)

    def __sub__(self, other) :
        """Returns difference"""
        if isinstance(other, pu.PolyBase):
            if not self.has_sametype(other):
                raise TypeError("Polynomial types differ")
            elif not self.has_samedomain(other):
                raise TypeError("Domains differ")
            elif not self.has_samewindow(other):
                raise TypeError("Windows differ")
            else:
                coef = ${nick}sub(self.coef, other.coef)
        else :
            try :
                coef = ${nick}sub(self.coef, other)
            except :
                return NotImplemented
        return self.__class__(coef, self.domain, self.window)

    def __mul__(self, other) :
        """Returns product"""
        if isinstance(other, pu.PolyBase):
            if not self.has_sametype(other):
                raise TypeError("Polynomial types differ")
            elif not self.has_samedomain(other):
                raise TypeError("Domains differ")
            elif not self.has_samewindow(other):
                raise TypeError("Windows differ")
            else:
                coef = ${nick}mul(self.coef, other.coef)
        else :
            try :
                coef = ${nick}mul(self.coef, other)
            except :
                return NotImplemented
        return self.__class__(coef, self.domain, self.window)

    def __div__(self, other):
        # set to __floordiv__,  /, for now.
        return self.__floordiv__(other)

    def __truediv__(self, other) :
        # there is no true divide if the rhs is not a Number, although it
        # could return the first n elements of an infinite series.
        # It is hard to see where n would come from, though.
        if not isinstance(other, Number) or isinstance(other, bool):
            form = "unsupported types for true division: '%s', '%s'"
            raise TypeError(form % (type(self), type(other)))
        return self.__floordiv__(other)

    def __floordiv__(self, other) :
        """Returns the quotient."""
        if isinstance(other, pu.PolyBase):
            if not self.has_sametype(other):
                raise TypeError("Polynomial types differ")
            elif not self.has_samedomain(other):
                raise TypeError("Domains differ")
            elif not self.has_samewindow(other):
                raise TypeError("Windows differ")
            else:
                quo, rem = ${nick}div(self.coef, other.coef)
        else :
            try :
                quo, rem = ${nick}div(self.coef, other)
            except :
                return NotImplemented
        return self.__class__(quo, self.domain, self.window)

    def __mod__(self, other) :
        """Returns the remainder."""
        if isinstance(other, pu.PolyBase):
            if not self.has_sametype(other):
                raise TypeError("Polynomial types differ")
            elif not self.has_samedomain(other):
                raise TypeError("Domains differ")
            elif not self.has_samewindow(other):
                raise TypeError("Windows differ")
            else:
                quo, rem = ${nick}div(self.coef, other.coef)
        else :
            try :
                quo, rem = ${nick}div(self.coef, other)
            except :
                return NotImplemented
        return self.__class__(rem, self.domain, self.window)

    def __divmod__(self, other) :
        """Returns quo, remainder"""
        if isinstance(other, self.__class__) :
            if not self.has_samedomain(other):
                raise TypeError("Domains are not equal")
            elif not self.has_samewindow(other):
                raise TypeError("Windows are not equal")
            else:
                quo, rem = ${nick}div(self.coef, other.coef)
        else :
            try :
                quo, rem = ${nick}div(self.coef, other)
            except :
                return NotImplemented
        quo = self.__class__(quo, self.domain, self.window)
        rem = self.__class__(rem, self.domain, self.window)
        return quo, rem

    def __pow__(self, other) :
        try :
            coef = ${nick}pow(self.coef, other, maxpower = self.maxpower)
        except :
            raise
        return self.__class__(coef, self.domain, self.window)

    def __radd__(self, other) :
        try :
            coef = ${nick}add(other, self.coef)
        except :
            return NotImplemented
        return self.__class__(coef, self.domain, self.window)

    def __rsub__(self, other):
        try :
            coef = ${nick}sub(other, self.coef)
        except :
            return NotImplemented
        return self.__class__(coef, self.domain, self.window)

    def __rmul__(self, other) :
        try :
            coef = ${nick}mul(other, self.coef)
        except :
            return NotImplemented
        return self.__class__(coef, self.domain, self.window)

    def __rdiv__(self, other):
        # set to __floordiv__ /.
        return self.__rfloordiv__(other)

    def __rtruediv__(self, other) :
        # An instance of PolyBase is not considered a
        # Number.
        return NotImplemented

    def __rfloordiv__(self, other) :
        try :
            quo, rem = ${nick}div(other, self.coef)
        except:
            return NotImplemented
        return self.__class__(quo, self.domain, self.window)

    def __rmod__(self, other) :
        try :
            quo, rem = ${nick}div(other, self.coef)
        except :
            return NotImplemented
        return self.__class__(rem, self.domain, self.window)

    def __rdivmod__(self, other) :
        try :
            quo, rem = ${nick}div(other, self.coef)
        except :
            return NotImplemented
        quo = self.__class__(quo, self.domain, self.window)
        rem = self.__class__(rem, self.domain, self.window)
        return quo, rem

    # Enhance me
    # some augmented arithmetic operations could be added here

    def __eq__(self, other) :
        res = isinstance(other, self.__class__) \
                and self.has_samecoef(other) \
                and self.has_samedomain(other) \
                and self.has_samewindow(other)
        return res

    def __ne__(self, other) :
        return not self.__eq__(other)

    #
    # Extra methods.
    #

    def copy(self) :
        """Return a copy.

        Return a copy of the current $name instance.

        Returns
        -------
        new_instance : $name
            Copy of current instance.

        """
        return self.__class__(self.coef, self.domain, self.window)

    def degree(self) :
        """The degree of the series.

        Notes
        -----
        .. versionadded:: 1.5.0

        """
        return len(self) - 1

    def cutdeg(self, deg) :
        """Truncate series to the given degree.

        Reduce the degree of the $name series to `deg` by discarding the
        high order terms. If `deg` is greater than the current degree a
        copy of the current series is returned. This can be useful in least
        squares where the coefficients of the high degree terms may be very
        small.

        Parameters
        ----------
        deg : non-negative int
            The series is reduced to degree `deg` by discarding the high
            order terms. The value of `deg` must be a non-negative integer.

        Returns
        -------
        new_instance : $name
            New instance of $name with reduced degree.

        Notes
        -----
        .. versionadded:: 1.5.0

        """
        return self.truncate(deg + 1)

    def trim(self, tol=0) :
        """Remove small leading coefficients

        Remove leading coefficients until a coefficient is reached whose
        absolute value greater than `tol` or the beginning of the series is
        reached. If all the coefficients would be removed the series is set to
        ``[0]``. A new $name instance is returned with the new coefficients.
        The current instance remains unchanged.

        Parameters
        ----------
        tol : non-negative number.
            All trailing coefficients less than `tol` will be removed.

        Returns
        -------
        new_instance : $name
            Contains the new set of coefficients.

        """
        coef = pu.trimcoef(self.coef, tol)
        return self.__class__(coef, self.domain, self.window)

    def truncate(self, size) :
        """Truncate series to length `size`.

        Reduce the $name series to length `size` by discarding the high
        degree terms. The value of `size` must be a positive integer. This
        can be useful in least squares where the coefficients of the
        high degree terms may be very small.

        Parameters
        ----------
        size : positive int
            The series is reduced to length `size` by discarding the high
            degree terms. The value of `size` must be a positive integer.

        Returns
        -------
        new_instance : $name
            New instance of $name with truncated coefficients.

        """
        isize = int(size)
        if isize != size or isize < 1 :
            raise ValueError("size must be a positive integer")
        if isize >= len(self.coef) :
            coef = self.coef
        else :
            coef = self.coef[:isize]
        return self.__class__(coef, self.domain, self.window)

    def convert(self, domain=None, kind=None, window=None) :
        """Convert to different class and/or domain.

        Parameters
        ----------
        domain : array_like, optional
            The domain of the converted series. If the value is None,
            the default domain of `kind` is used.
        kind : class, optional
            The polynomial series type class to which the current instance
            should be converted. If kind is None, then the class of the
            current instance is used.
        window : array_like, optional
            The window of the converted series. If the value is None,
            the default window of `kind` is used.

        Returns
        -------
        new_series_instance : `kind`
            The returned class can be of different type than the current
            instance and/or have a different domain.

        Notes
        -----
        Conversion between domains and class types can result in
        numerically ill defined series.

        Examples
        --------

        """
        if kind is None:
            kind = $name
        if domain is None:
            domain = kind.domain
        if window is None:
            window = kind.window
        return self(kind.identity(domain, window=window))

    def mapparms(self) :
        """Return the mapping parameters.

        The returned values define a linear map ``off + scl*x`` that is
        applied to the input arguments before the series is evaluated. The
        map depends on the ``domain`` and ``window``; if the current
        ``domain`` is equal to the ``window`` the resulting map is the
        identity.  If the coefficients of the ``$name`` instance are to be
        used by themselves outside this class, then the linear function
        must be substituted for the ``x`` in the standard representation of
        the base polynomials.

        Returns
        -------
        off, scl : floats or complex
            The mapping function is defined by ``off + scl*x``.

        Notes
        -----
        If the current domain is the interval ``[l_1, r_1]`` and the window
        is ``[l_2, r_2]``, then the linear mapping function ``L`` is
        defined by the equations::

            L(l_1) = l_2
            L(r_1) = r_2

        """
        return pu.mapparms(self.domain, self.window)

    def integ(self, m=1, k=[], lbnd=None) :
        """Integrate.

        Return an instance of $name that is the definite integral of the
        current series. Refer to `${nick}int` for full documentation.

        Parameters
        ----------
        m : non-negative int
            The number of integrations to perform.
        k : array_like
            Integration constants. The first constant is applied to the
            first integration, the second to the second, and so on. The
            list of values must less than or equal to `m` in length and any
            missing values are set to zero.
        lbnd : Scalar
            The lower bound of the definite integral.

        Returns
        -------
        integral : $name
            The integral of the series using the same domain.

        See Also
        --------
        ${nick}int : similar function.
        ${nick}der : similar function for derivative.

        """
        off, scl = self.mapparms()
        if lbnd is None :
            lbnd = 0
        else :
            lbnd = off + scl*lbnd
        coef = ${nick}int(self.coef, m, k, lbnd, 1./scl)
        return self.__class__(coef, self.domain, self.window)

    def deriv(self, m=1):
        """Differentiate.

        Return an instance of $name that is the derivative of the current
        series.  Refer to `${nick}der` for full documentation.

        Parameters
        ----------
        m : non-negative int
            The number of integrations to perform.

        Returns
        -------
        derivative : $name
            The derivative of the series using the same domain.

        See Also
        --------
        ${nick}der : similar function.
        ${nick}int : similar function for integration.

        """
        off, scl = self.mapparms()
        coef = ${nick}der(self.coef, m, scl)
        return self.__class__(coef, self.domain, self.window)

    def roots(self) :
        """Return list of roots.

        Return ndarray of roots for this series. See `${nick}roots` for
        full documentation. Note that the accuracy of the roots is likely to
        decrease the further outside the domain they lie.

        See Also
        --------
        ${nick}roots : similar function
        ${nick}fromroots : function to go generate series from roots.

        """
        roots = ${nick}roots(self.coef)
        return pu.mapdomain(roots, self.window, self.domain)

    def linspace(self, n=100, domain=None):
        """Return x,y values at equally spaced points in domain.

        Returns x, y values at `n` linearly spaced points across domain.
        Here y is the value of the polynomial at the points x. By default
        the domain is the same as that of the $name instance.  This method
        is intended mostly as a plotting aid.

        Parameters
        ----------
        n : int, optional
            Number of point pairs to return. The default value is 100.
        domain : {None, array_like}
            If not None, the specified domain is used instead of that of
            the calling instance. It should be of the form ``[beg,end]``.
            The default is None.

        Returns
        -------
        x, y : ndarrays
            ``x`` is equal to linspace(self.domain[0], self.domain[1], n)
            ``y`` is the polynomial evaluated at ``x``.

        .. versionadded:: 1.5.0

        """
        if domain is None:
            domain = self.domain
        x = np.linspace(domain[0], domain[1], n)
        y = self(x)
        return x, y



    @staticmethod
    def fit(x, y, deg, domain=None, rcond=None, full=False, w=None,
        window=$domain):
        """Least squares fit to data.

        Return a `$name` instance that is the least squares fit to the data
        `y` sampled at `x`. Unlike `${nick}fit`, the domain of the returned
        instance can be specified and this will often result in a superior
        fit with less chance of ill conditioning. Support for NA was added
        in version 1.7.0. See `${nick}fit` for full documentation of the
        implementation.

        Parameters
        ----------
        x : array_like, shape (M,)
            x-coordinates of the M sample points ``(x[i], y[i])``.
        y : array_like, shape (M,) or (M, K)
            y-coordinates of the sample points. Several data sets of sample
            points sharing the same x-coordinates can be fitted at once by
            passing in a 2D-array that contains one dataset per column.
        deg : int
            Degree of the fitting polynomial.
        domain : {None, [beg, end], []}, optional
            Domain to use for the returned $name instance. If ``None``,
            then a minimal domain that covers the points `x` is chosen.  If
            ``[]`` the default domain ``$domain`` is used. The default
            value is $domain in numpy 1.4.x and ``None`` in later versions.
            The ``[]`` value was added in numpy 1.5.0.
        rcond : float, optional
            Relative condition number of the fit. Singular values smaller
            than this relative to the largest singular value will be
            ignored. The default value is len(x)*eps, where eps is the
            relative precision of the float type, about 2e-16 in most
            cases.
        full : bool, optional
            Switch determining nature of return value. When it is False
            (the default) just the coefficients are returned, when True
            diagnostic information from the singular value decomposition is
            also returned.
        w : array_like, shape (M,), optional
            Weights. If not None the contribution of each point
            ``(x[i],y[i])`` to the fit is weighted by `w[i]`. Ideally the
            weights are chosen so that the errors of the products
            ``w[i]*y[i]`` all have the same variance.  The default value is
            None.
            .. versionadded:: 1.5.0
        window : {[beg, end]}, optional
            Window to use for the returned $name instance. The default
            value is ``$domain``
            .. versionadded:: 1.6.0

        Returns
        -------
        least_squares_fit : instance of $name
            The $name instance is the least squares fit to the data and
            has the domain specified in the call.

        [residuals, rank, singular_values, rcond] : only if `full` = True
            Residuals of the least squares fit, the effective rank of the
            scaled Vandermonde matrix and its singular values, and the
            specified value of `rcond`. For more details, see
            `linalg.lstsq`.

        See Also
        --------
        ${nick}fit : similar function

        """
        if domain is None:
            domain = pu.getdomain(x)
        elif type(domain) is list and len(domain) == 0:
            domain = $domain

        if type(window) is list and len(window) == 0:
            window = $domain

        xnew = pu.mapdomain(x, domain, window)
        res = ${nick}fit(xnew, y, deg, w=w, rcond=rcond, full=full)
        if full :
            [coef, status] = res
            return $name(coef, domain=domain, window=window), status
        else :
            coef = res
            return $name(coef, domain=domain, window=window)

    @staticmethod
    def fromroots(roots, domain=$domain, window=$domain) :
        """Return $name instance with specified roots.

        Returns an instance of $name representing the product
        ``(x - r[0])*(x - r[1])*...*(x - r[n-1])``, where ``r`` is the
        list of roots.

        Parameters
        ----------
        roots : array_like
            List of roots.
        domain : {array_like, None}, optional
            Domain for the resulting instance of $name. If none the domain
            is the interval from the smallest root to the largest. The
            default is $domain.
        window : array_like, optional
            Window for the resulting instance of $name. The default value
            is $domain.

        Returns
        -------
        object : $name instance
            Series with the specified roots.

        See Also
        --------
        ${nick}fromroots : equivalent function

        """
        [roots] = pu.as_series([roots], trim=False)
        if domain is None :
            domain = pu.getdomain(roots)
        deg = len(roots)
        off, scl = pu.mapparms(domain, window)
        rnew = off + scl*roots
        coef = ${nick}fromroots(rnew) / scl**deg
        return $name(coef, domain=domain, window=window)

    @staticmethod
    def identity(domain=$domain, window=$domain) :
        """Identity function.

        If ``p`` is the returned $name object, then ``p(x) == x`` for all
        values of x.

        Parameters
        ----------
        domain : array_like
            The resulting array must be of the form ``[beg, end]``, where
            ``beg`` and ``end`` are the endpoints of the domain.
        window : array_like
            The resulting array must be if the form ``[beg, end]``, where
            ``beg`` and ``end`` are the endpoints of the window.

        Returns
        -------
        identity : $name instance

        """
        off, scl = pu.mapparms(window, domain)
        coef = ${nick}line(off, scl)
        return $name(coef, domain, window)

    @staticmethod
    def basis(deg, domain=$domain, window=$domain):
        """$name polynomial of degree `deg`.

        Returns an instance of the $name polynomial of degree `d`.

        Parameters
        ----------
        deg : int
            Degree of the $name polynomial. Must be >= 0.
        domain : array_like
            The resulting array must be of the form ``[beg, end]``, where
            ``beg`` and ``end`` are the endpoints of the domain.
        window : array_like
            The resulting array must be if the form ``[beg, end]``, where
            ``beg`` and ``end`` are the endpoints of the window.

        Returns
        p : $name instance

        Notes
        -----
        .. versionadded:: 1.7.0

        """
        ideg = int(deg)
        if ideg != deg or ideg < 0:
            raise ValueError("deg must be non-negative integer")
        return $name([0]*ideg + [1], domain, window)

    @staticmethod
    def cast(series, domain=$domain, window=$domain):
        """Convert instance to equivalent $name series.

        The `series` is expected to be an instance of some polynomial
        series of one of the types supported by by the numpy.polynomial
        module, but could be some other class that supports the convert
        method.

        Parameters
        ----------
        series : series
            The instance series to be converted.
        domain : array_like
            The resulting array must be of the form ``[beg, end]``, where
            ``beg`` and ``end`` are the endpoints of the domain.
        window : array_like
            The resulting array must be if the form ``[beg, end]``, where
            ``beg`` and ``end`` are the endpoints of the window.

        Returns
        p : $name instance
            A $name series equal to the `poly` series.

        See Also
        --------
        convert -- similar instance method

        Notes
        -----
        .. versionadded:: 1.7.0

        """
        return series.convert(domain, $name, window)

''')
