# doubledouble.pyx
"""
Double-double precision arithmetic mentioned at start of referenced paper.
Provides ~106 bits of precision via error-free transformations.

Built referencing paper and fork of QD library

Reference:
    Hida, Li, Bailey
    https://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf

    Official QD repository:
    https://github.com/BL-highprecision/QD

#? Note for future

We could implement fast_two_sum in the future.
However, that requires |a| >= |b| so we'll keep it
concise for now.
"""

from libc.math cimport fma


cdef inline void two_sum(
    double a,
    double b,
    double* out_sum,
    double* out_err,
) noexcept:
    cdef:
        double s = a + b
        double v = s - a

    out_sum[0] = s
    out_err[0] = (a - (s - v)) + (b - v)

cdef inline void two_prod(
    double a,
    double b,
    double* out_prod,
    double* out_err,
) noexcept:

    out_prod[0] = a * b
    out_err[0] = fma(a, b, -out_prod[0])

# ---------------------------------------------------------------------------
# DoubleDouble Class
# ---------------------------------------------------------------------------

# Theoretically we could support float128 on all platforms with this
# I do not wish to open that pandora's box just yet
# !Therefore, a decision has to be made here before merge
# ?Actually, could cause awkward usage, tbd
cdef class DoubleDouble:
    """
    Double-double floating point number.

    Represents a number as an unevaluated sum of two doubles (hi + lo),
    providing approximately 106 bits of precision.

    Parameters
    ----------
    hi : float, default 0.0
        High-order component.
    lo : float, default 0.0
        Low-order component.
    """
    cdef:
        public double hi
        public double lo

    def __init__(self, double hi = 0.0, double lo=0.0):
        self.hi = hi
        self.lo = lo

    def __repr__(self) -> str:
        return f"DoubleDouble({self.hi}, {self.lo})"

    def collapse(self) -> float:
        """
        Collapse to a single float64.

        Returns
        -------
        float
            The sum hi + lo as a standard double.
        """
        return self.hi + self.lo

    def __float__(self):
        return self.hi + self.lo

    def __str__(self):
        return str(self.hi + self.lo)

    def __add__(self, other):
        cdef:
            double a_hi, a_lo, b_hi, b_lo
            double s, e, v, res_hi, res_lo

        a_hi = self.hi
        a_lo = self.lo
        if isinstance(other, DoubleDouble):
            # !I dont know if we should allow this due to logistical complications
            b_hi = (<DoubleDouble>other).hi
            b_lo = (<DoubleDouble>other).lo
        elif isinstance(other, (int, float)):
            b_hi = <double>other
            b_lo = 0.0
        else:
            return NotImplemented

        two_sum(a_hi, b_hi, &s, &e)
        # This is the 'sloppy' version due to this line
        # Could use two_sum here for more precision in future
        v = a_lo + b_lo + e
        two_sum(s, v, &res_hi, &res_lo)

        return DoubleDouble(res_hi, res_lo)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        cdef:
            double a_hi, a_lo, b_hi, b_lo
            double s, e, v, res_hi, res_lo

        a_hi = self.hi
        a_lo = self.lo
        if isinstance(other, DoubleDouble):
            # !I dont know if we should allow this due to logistical complications
            b_hi = -(<DoubleDouble>other).hi
            b_lo = -(<DoubleDouble>other).lo
        elif isinstance(other, (int, float)):
            b_hi = -<double>other
            b_lo = 0.0
        else:
            return NotImplemented

        two_sum(a_hi, b_hi, &s, &e)
        # This is the 'sloppy' version due to this line
        # Could use two_sum here for more precision in future
        v = a_lo + b_lo + e
        two_sum(s, v, &res_hi, &res_lo)

        return DoubleDouble(res_hi, res_lo)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        cdef:
            double a_hi, a_lo, b_hi, b_lo
            double p, e, rh, rl

        a_hi = self.hi
        a_lo = self.lo
        if isinstance(other, DoubleDouble):
            b_hi = (<DoubleDouble>other).hi
            b_lo = (<DoubleDouble>other).lo
        elif isinstance(other, (int, float)):
            b_hi = <double>other
            b_lo = 0.0
        else:
            return NotImplemented

        two_prod(a_hi, b_hi, &p, &e)
        e = fma(a_hi, b_lo, e)
        e = fma(a_lo, b_hi, e)
        two_sum(p, e, &rh, &rl)

        return DoubleDouble(rh, rl)

    def __rmul__(self, other):
        return self.__mul__(other)
