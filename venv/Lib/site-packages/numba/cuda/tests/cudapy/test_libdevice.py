import numpy as np
from numba.core import types
from numba.cuda.testing import skip_on_cudasim, unittest, CUDATestCase
from numba import cuda
from numba.cuda import libdevice, compile_ptx
from numba.cuda.libdevicefuncs import functions, create_signature


def use_sincos(s, c, x):
    i = cuda.grid(1)

    if i < len(x):
        sr, cr = libdevice.sincos(x[i])
        s[i] = sr
        c[i] = cr


def use_frexp(frac, exp, x):
    i = cuda.grid(1)

    if i < len(x):
        fracr, expr = libdevice.frexp(x[i])
        frac[i] = fracr
        exp[i] = expr


def use_sad(r, x, y, z):
    i = cuda.grid(1)

    if i < len(x):
        r[i] = libdevice.sad(x[i], y[i], z[i])


@skip_on_cudasim('Libdevice functions are not supported on cudasim')
class TestLibdevice(CUDATestCase):
    """
    Some tests of libdevice function wrappers that check the returned values.

    These are mainly to check that the generation of the implementations
    results in correct typing and lowering for each type of function return
    (e.g. scalar return, UniTuple return, Tuple return, etc.).
    """

    def test_sincos(self):
        # Tests return of a UniTuple from a libdevice function
        arr = np.arange(100, dtype=np.float64)
        sres = np.zeros_like(arr)
        cres = np.zeros_like(arr)

        cufunc = cuda.jit(use_sincos)
        cufunc[4, 32](sres, cres, arr)

        np.testing.assert_allclose(np.cos(arr), cres)
        np.testing.assert_allclose(np.sin(arr), sres)

    def test_frexp(self):
        # Test return of a Tuple from a libdevice function
        arr = np.linspace(start=1.0, stop=10.0, num=100, dtype=np.float64)
        fracres = np.zeros_like(arr)
        expres = np.zeros(shape=arr.shape, dtype=np.int32)

        cufunc = cuda.jit(use_frexp)
        cufunc[4, 32](fracres, expres, arr)

        frac_expect, exp_expect = np.frexp(arr)

        np.testing.assert_array_equal(frac_expect, fracres)
        np.testing.assert_array_equal(exp_expect, expres)

    def test_sad(self):
        # Test return of a scalar from a libdevice function
        x = np.arange(0, 200, 2)
        y = np.arange(50, 150)
        z = np.arange(15, 115)
        r = np.zeros_like(x)

        cufunc = cuda.jit(use_sad)
        cufunc[4, 32](r, x, y, z)

        np.testing.assert_array_equal(np.abs(x - y) + z, r)


# A template for generating tests of compiling calls to libdevice functions.
# The purpose of the call and assignment of the return variables is to ensure
# the actual function implementations are not thrown away resulting in a PTX
# implementation that only contains the ret instruction - this may hide certain
# errors.
function_template = """\
from numba.cuda import libdevice

def pyfunc(%(pyargs)s):
    ret = libdevice.%(func)s(%(funcargs)s)
    %(retvars)s = ret
"""


def make_test_call(libname):
    """
    Generates a test function for each libdevice function.
    """

    def _test_call_functions(self):
        # Strip off '__nv_' from libdevice name to get Python name
        apiname = libname[5:]
        apifunc = getattr(libdevice, apiname)
        retty, args = functions[libname]
        sig = create_signature(retty, args)

        # Construct arguments to the libdevice function. These are all
        # non-pointer arguments to the underlying bitcode function.
        funcargs = ", ".join(['a%d' % i for i, arg in enumerate(args) if not
                              arg.is_ptr])

        # Arguments to the Python function (`pyfunc` in the template above) are
        # the arguments to the libdevice function, plus as many extra arguments
        # as there are in the return type of the libdevice function - one for
        # scalar-valued returns, or the length of the tuple for tuple-valued
        # returns.
        if isinstance(sig.return_type, (types.Tuple, types.UniTuple)):
            # Start with the parameters for the return values
            pyargs = ", ".join(['r%d' % i for i in
                                range(len(sig.return_type))])
            # Add the parameters for the argument values
            pyargs += ", " + funcargs
            # Generate the unpacking of the return value from the libdevice
            # function into the Python function return values (`r0`, `r1`,
            # etc.).
            retvars = ", ".join(['r%d[0]' % i for i in
                                 range(len(sig.return_type))])
        else:
            # Scalar return is a more straightforward case
            pyargs = "r0, " + funcargs
            retvars = "r0[0]"

        # Create the string containing the function to compile
        d = { 'func': apiname,
              'pyargs': pyargs,
              'funcargs': funcargs,
              'retvars': retvars }
        code = function_template % d

        # Convert the string to a Python function
        locals = {}
        exec(code, globals(), locals)
        pyfunc = locals['pyfunc']

        # Compute the signature for compilation. This mirrors the creation of
        # arguments to the Python function above.
        pyargs = [ arg.ty for arg in args if not arg.is_ptr ]
        if isinstance(sig.return_type, (types.Tuple, types.UniTuple)):
            pyreturns = [ret[::1] for ret in sig.return_type]
            pyargs = pyreturns + pyargs
        else:
            pyargs.insert(0, sig.return_type[::1])

        pyargs = tuple(pyargs)
        ptx, resty = compile_ptx(pyfunc, pyargs)

        # If the function body was discarded by optimization (therefore making
        # the test a bit weak), there won't be any loading of parameters -
        # ensure that a load from parameters occurs somewhere in the PTX
        self.assertIn('ld.param', ptx)

        # Returning the result (through a passed-in array) should also require
        # a store to global memory, so check for at least one of those too.
        self.assertIn('st.global', ptx)

    return _test_call_functions


@skip_on_cudasim('Compilation to PTX is not supported on cudasim')
class TestLibdeviceCompilation(unittest.TestCase):
    """
    Class for holding all tests of compiling calls to libdevice functions. We
    generate the actual tests in this class (as opposed to using subTest and
    one test within this class) because there are a lot of tests, and it makes
    the test suite appear frozen to test them all as subTests in one test.
    """


for libname in functions:
    setattr(TestLibdeviceCompilation, 'test_%s' % libname,
            make_test_call(libname))


if __name__ == '__main__':
    unittest.main()
