"""Machine limits for Float32 and Float64 and (long double) if available...

"""
from __future__ import division, absolute_import, print_function

__all__ = ['finfo', 'iinfo']

from .machar import MachAr
from . import numeric
from . import numerictypes as ntypes
from .numeric import array

def _frz(a):
    """fix rank-0 --> rank-1"""
    if a.ndim == 0: a.shape = (1,)
    return a

_convert_to_float = {
    ntypes.csingle: ntypes.single,
    ntypes.complex_: ntypes.float_,
    ntypes.clongfloat: ntypes.longfloat
    }

class finfo(object):
    """
    finfo(dtype)

    Machine limits for floating point types.

    Attributes
    ----------
    eps : float
        The smallest representable positive number such that
        ``1.0 + eps != 1.0``.  Type of `eps` is an appropriate floating
        point type.
    epsneg : floating point number of the appropriate type
        The smallest representable positive number such that
        ``1.0 - epsneg != 1.0``.
    iexp : int
        The number of bits in the exponent portion of the floating point
        representation.
    machar : MachAr
        The object which calculated these parameters and holds more
        detailed information.
    machep : int
        The exponent that yields `eps`.
    max : floating point number of the appropriate type
        The largest representable number.
    maxexp : int
        The smallest positive power of the base (2) that causes overflow.
    min : floating point number of the appropriate type
        The smallest representable number, typically ``-max``.
    minexp : int
        The most negative power of the base (2) consistent with there
        being no leading 0's in the mantissa.
    negep : int
        The exponent that yields `epsneg`.
    nexp : int
        The number of bits in the exponent including its sign and bias.
    nmant : int
        The number of bits in the mantissa.
    precision : int
        The approximate number of decimal digits to which this kind of
        float is precise.
    resolution : floating point number of the appropriate type
        The approximate decimal resolution of this type, i.e.,
        ``10**-precision``.
    tiny : float
        The smallest positive usable number.  Type of `tiny` is an
        appropriate floating point type.

    Parameters
    ----------
    dtype : float, dtype, or instance
        Kind of floating point data-type about which to get information.

    See Also
    --------
    MachAr : The implementation of the tests that produce this information.
    iinfo : The equivalent for integer data types.

    Notes
    -----
    For developers of NumPy: do not instantiate this at the module level.
    The initial calculation of these parameters is expensive and negatively
    impacts import times.  These objects are cached, so calling ``finfo()``
    repeatedly inside your functions is not a problem.

    """

    _finfo_cache = {}

    def __new__(cls, dtype):
        try:
            dtype = numeric.dtype(dtype)
        except TypeError:
            # In case a float instance was given
            dtype = numeric.dtype(type(dtype))

        obj = cls._finfo_cache.get(dtype, None)
        if obj is not None:
            return obj
        dtypes = [dtype]
        newdtype = numeric.obj2sctype(dtype)
        if newdtype is not dtype:
            dtypes.append(newdtype)
            dtype = newdtype
        if not issubclass(dtype, numeric.inexact):
            raise ValueError("data type %r not inexact" % (dtype))
        obj = cls._finfo_cache.get(dtype, None)
        if obj is not None:
            return obj
        if not issubclass(dtype, numeric.floating):
            newdtype = _convert_to_float[dtype]
            if newdtype is not dtype:
                dtypes.append(newdtype)
                dtype = newdtype
        obj = cls._finfo_cache.get(dtype, None)
        if obj is not None:
            return obj
        obj = object.__new__(cls)._init(dtype)
        for dt in dtypes:
            cls._finfo_cache[dt] = obj
        return obj

    def _init(self, dtype):
        self.dtype = numeric.dtype(dtype)
        if dtype is ntypes.double:
            itype = ntypes.int64
            fmt = '%24.16e'
            precname = 'double'
        elif dtype is ntypes.single:
            itype = ntypes.int32
            fmt = '%15.7e'
            precname = 'single'
        elif dtype is ntypes.longdouble:
            itype = ntypes.longlong
            fmt = '%s'
            precname = 'long double'
        elif dtype is ntypes.half:
            itype = ntypes.int16
            fmt = '%12.5e'
            precname = 'half'
        else:
            raise ValueError(repr(dtype))

        machar = MachAr(lambda v:array([v], dtype),
                        lambda v:_frz(v.astype(itype))[0],
                        lambda v:array(_frz(v)[0], dtype),
                        lambda v: fmt % array(_frz(v)[0], dtype),
                        'numpy %s precision floating point number' % precname)

        for word in ['precision', 'iexp',
                     'maxexp', 'minexp', 'negep',
                     'machep']:
            setattr(self, word, getattr(machar, word))
        for word in ['tiny', 'resolution', 'epsneg']:
            setattr(self, word, getattr(machar, word).flat[0])
        self.max = machar.huge.flat[0]
        self.min = -self.max
        self.eps = machar.eps.flat[0]
        self.nexp = machar.iexp
        self.nmant = machar.it
        self.machar = machar
        self._str_tiny = machar._str_xmin.strip()
        self._str_max = machar._str_xmax.strip()
        self._str_epsneg = machar._str_epsneg.strip()
        self._str_eps = machar._str_eps.strip()
        self._str_resolution = machar._str_resolution.strip()
        return self

    def __str__(self):
        return '''\
Machine parameters for %(dtype)s
---------------------------------------------------------------------
precision=%(precision)3s   resolution= %(_str_resolution)s
machep=%(machep)6s   eps=        %(_str_eps)s
negep =%(negep)6s   epsneg=     %(_str_epsneg)s
minexp=%(minexp)6s   tiny=       %(_str_tiny)s
maxexp=%(maxexp)6s   max=        %(_str_max)s
nexp  =%(nexp)6s   min=        -max
---------------------------------------------------------------------
''' % self.__dict__

    def __repr__(self):
        c = self.__class__.__name__
        d = self.__dict__.copy()
        d['klass'] = c
        return ("%(klass)s(resolution=%(resolution)s, min=-%(_str_max)s," \
               + " max=%(_str_max)s, dtype=%(dtype)s)") \
                % d


class iinfo(object):
    """
    iinfo(type)

    Machine limits for integer types.

    Attributes
    ----------
    min : int
        The smallest integer expressible by the type.
    max : int
        The largest integer expressible by the type.

    Parameters
    ----------
    type : integer type, dtype, or instance
        The kind of integer data type to get information about.

    See Also
    --------
    finfo : The equivalent for floating point data types.

    Examples
    --------
    With types:

    >>> ii16 = np.iinfo(np.int16)
    >>> ii16.min
    -32768
    >>> ii16.max
    32767
    >>> ii32 = np.iinfo(np.int32)
    >>> ii32.min
    -2147483648
    >>> ii32.max
    2147483647

    With instances:

    >>> ii32 = np.iinfo(np.int32(10))
    >>> ii32.min
    -2147483648
    >>> ii32.max
    2147483647

    """

    _min_vals = {}
    _max_vals = {}

    def __init__(self, int_type):
        try:
            self.dtype = numeric.dtype(int_type)
        except TypeError:
            self.dtype = numeric.dtype(type(int_type))
        self.kind = self.dtype.kind
        self.bits = self.dtype.itemsize * 8
        self.key = "%s%d" % (self.kind, self.bits)
        if not self.kind in 'iu':
            raise ValueError("Invalid integer data type.")

    def min(self):
        """Minimum value of given dtype."""
        if self.kind == 'u':
            return 0
        else:
            try:
                val = iinfo._min_vals[self.key]
            except KeyError:
                val = int(-(1 << (self.bits-1)))
                iinfo._min_vals[self.key] = val
            return val

    min = property(min)

    def max(self):
        """Maximum value of given dtype."""
        try:
            val = iinfo._max_vals[self.key]
        except KeyError:
            if self.kind == 'u':
                val = int((1 << self.bits) - 1)
            else:
                val = int((1 << (self.bits-1)) - 1)
            iinfo._max_vals[self.key] = val
        return val

    max = property(max)

    def __str__(self):
        """String representation."""
        return '''\
Machine parameters for %(dtype)s
---------------------------------------------------------------------
min = %(min)s
max = %(max)s
---------------------------------------------------------------------
''' % {'dtype': self.dtype, 'min': self.min, 'max': self.max}

    def __repr__(self):
        return "%s(min=%s, max=%s, dtype=%s)" % (self.__class__.__name__,
                                    self.min, self.max, self.dtype)

if __name__ == '__main__':
    f = finfo(ntypes.single)
    print('single epsilon:', f.eps)
    print('single tiny:', f.tiny)
    f = finfo(ntypes.float)
    print('float epsilon:', f.eps)
    print('float tiny:', f.tiny)
    f = finfo(ntypes.longfloat)
    print('longfloat epsilon:', f.eps)
    print('longfloat tiny:', f.tiny)
