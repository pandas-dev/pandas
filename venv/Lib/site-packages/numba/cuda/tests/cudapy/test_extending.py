from numba.cuda.testing import skip_on_cudasim, unittest, CUDATestCase

import numpy as np
from numba import config, cuda, njit, types


class Interval:
    """
    A half-open interval on the real number line.
    """
    def __init__(self, lo, hi):
        self.lo = lo
        self.hi = hi

    def __repr__(self):
        return 'Interval(%f, %f)' % (self.lo, self.hi)

    @property
    def width(self):
        return self.hi - self.lo


@njit
def interval_width(interval):
    return interval.width


@njit
def sum_intervals(i, j):
    return Interval(i.lo + j.lo, i.hi + j.hi)


if not config.ENABLE_CUDASIM:
    from numba.core import cgutils
    from numba.core.extending import (lower_builtin, make_attribute_wrapper,
                                      models, register_model, type_callable,
                                      typeof_impl)
    from numba.core.typing.templates import AttributeTemplate
    from numba.cuda.cudadecl import registry as cuda_registry
    from numba.cuda.cudaimpl import lower_attr as cuda_lower_attr

    class IntervalType(types.Type):
        def __init__(self):
            super().__init__(name='Interval')

    interval_type = IntervalType()

    @typeof_impl.register(Interval)
    def typeof_interval(val, c):
        return interval_type

    @type_callable(Interval)
    def type_interval(context):
        def typer(lo, hi):
            if isinstance(lo, types.Float) and isinstance(hi, types.Float):
                return interval_type
        return typer

    @register_model(IntervalType)
    class IntervalModel(models.StructModel):
        def __init__(self, dmm, fe_type):
            members = [
                ('lo', types.float64),
                ('hi', types.float64),
            ]
            models.StructModel.__init__(self, dmm, fe_type, members)

    make_attribute_wrapper(IntervalType, 'lo', 'lo')
    make_attribute_wrapper(IntervalType, 'hi', 'hi')

    @lower_builtin(Interval, types.Float, types.Float)
    def impl_interval(context, builder, sig, args):
        typ = sig.return_type
        lo, hi = args
        interval = cgutils.create_struct_proxy(typ)(context, builder)
        interval.lo = lo
        interval.hi = hi
        return interval._getvalue()

    @cuda_registry.register_attr
    class Interval_attrs(AttributeTemplate):
        key = IntervalType

        def resolve_width(self, mod):
            return types.float64

    @cuda_lower_attr(IntervalType, 'width')
    def cuda_Interval_width(context, builder, sig, arg):
        lo = builder.extract_value(arg, 0)
        hi = builder.extract_value(arg, 1)
        return builder.fsub(hi, lo)


@skip_on_cudasim('Extensions not supported in the simulator')
class TestExtending(CUDATestCase):
    def test_attributes(self):
        @cuda.jit
        def f(r, x):
            iv = Interval(x[0], x[1])
            r[0] = iv.lo
            r[1] = iv.hi

        x = np.asarray((1.5, 2.5))
        r = np.zeros_like(x)

        f[1, 1](r, x)

        np.testing.assert_equal(r, x)

    def test_property(self):
        @cuda.jit
        def f(r, x):
            iv = Interval(x[0], x[1])
            r[0] = iv.width

        x = np.asarray((1.5, 2.5))
        r = np.zeros(1)

        f[1, 1](r, x)

        np.testing.assert_allclose(r[0], x[1] - x[0])

    def test_extension_type_as_arg(self):
        @cuda.jit
        def f(r, x):
            iv = Interval(x[0], x[1])
            r[0] = interval_width(iv)

        x = np.asarray((1.5, 2.5))
        r = np.zeros(1)

        f[1, 1](r, x)

        np.testing.assert_allclose(r[0], x[1] - x[0])

    def test_extension_type_as_retvalue(self):
        @cuda.jit
        def f(r, x):
            iv1 = Interval(x[0], x[1])
            iv2 = Interval(x[2], x[3])
            iv_sum = sum_intervals(iv1, iv2)
            r[0] = iv_sum.lo
            r[1] = iv_sum.hi

        x = np.asarray((1.5, 2.5, 3.0, 4.0))
        r = np.zeros(2)

        f[1, 1](r, x)

        expected = np.asarray((x[0] + x[2], x[1] + x[3]))
        np.testing.assert_allclose(r, expected)


if __name__ == '__main__':
    unittest.main()
