import re
from io import StringIO

import numba
from numba.core import types
from numba import jit, njit
from numba.tests.support import override_config, TestCase
import unittest

try:
    import jinja2
except ImportError:
    jinja2 = None

try:
    import pygments
except ImportError:
    pygments = None


@unittest.skipIf(jinja2 is None, "please install the 'jinja2' package")
class TestAnnotation(TestCase):

    @TestCase.run_test_in_subprocess # annotations compound per module
    def test_exercise_code_path(self):
        """
        Ensures template.html is available
        """

        def foo(n, a):
            s = a
            for i in range(n):
                s += i
            return s

        cfunc = njit((types.int32, types.int32))(foo)
        cres = cfunc.overloads[cfunc.signatures[0]]
        ta = cres.type_annotation

        buf = StringIO()
        ta.html_annotate(buf)
        output = buf.getvalue()
        buf.close()
        self.assertIn("foo", output)

    @TestCase.run_test_in_subprocess # annotations compound per module
    def test_exercise_code_path_with_lifted_loop(self):
        """
        Ensures that lifted loops are handled correctly in obj mode
        """
        # the functions to jit
        def bar(x):
            return x

        def foo(x):
            h = 0.
            for i in range(x): # py 38 needs two loops for one to lift?!
                h = h + i
            for k in range(x):
                h = h + k
            if x:
                h = h - bar(x)
            return h

        # compile into an isolated context
        cfunc = jit((types.intp,), forceobj=True, looplift=True)(foo)
        cres = cfunc.overloads[cfunc.signatures[0]]

        ta = cres.type_annotation

        buf = StringIO()
        ta.html_annotate(buf)
        output = buf.getvalue()
        buf.close()
        self.assertIn("bar", output)
        self.assertIn("foo", output)
        self.assertIn("LiftedLoop", output)

    @TestCase.run_test_in_subprocess # annotations compound per module
    def test_html_output_with_lifted_loop(self):
        """
        Test some format and behavior of the html annotation with lifted loop
        """
        @numba.jit(forceobj=True)
        def udt(x):
            object()  # to force object mode
            z = 0
            for i in range(x):  # this line is tagged
                z += i
            return z

        # Regex pattern to check for the "lifted_tag" in the line of the loop
        re_lifted_tag = re.compile(
            r'<td class="lifted_tag">\s*'
            r'\s*<details>'
            r'\s*<summary>'
            r'\s*<code>'
            r'\s*[0-9]+:'
            r'\s*[&nbsp;]+for i in range\(x\):  # this line is tagged\s*',
            re.MULTILINE)

        # Compile int64 version
        sig_i64 = (types.int64,)
        udt.compile(sig_i64)  # compile with lifted loop
        cres = udt.overloads[sig_i64]

        # Make html output
        buf = StringIO()
        cres.type_annotation.html_annotate(buf)
        output = buf.getvalue()
        buf.close()

        # There should be only one function output.
        self.assertEqual(output.count("Function name: udt"), 1)

        sigfmt = "with signature: {} -&gt; pyobject"
        self.assertEqual(output.count(sigfmt.format(sig_i64)), 1)
        # Ensure the loop is tagged
        self.assertEqual(len(re.findall(re_lifted_tag, output)), 1,
                         msg='%s not found in %s' % (re_lifted_tag, output))

        # Compile float64 version
        sig_f64 = (types.float64,)
        udt.compile(sig_f64)
        cres = udt.overloads[sig_f64]

        # Make html output
        buf = StringIO()
        cres.type_annotation.html_annotate(buf)
        output = buf.getvalue()
        buf.close()

        # There should be two function output
        self.assertEqual(output.count("Function name: udt"), 2)
        self.assertEqual(output.count(sigfmt.format(sig_i64)), 1)
        self.assertEqual(output.count(sigfmt.format(sig_f64)), 1)
        # Ensure the loop is tagged in both output
        self.assertEqual(len(re.findall(re_lifted_tag, output)), 2)

    @unittest.skipIf(pygments is None, "please install the 'pygments' package")
    def test_pretty_print(self):

        @numba.njit
        def foo(x, y):
            return x, y

        foo(1, 2)
        # Exercise the method
        foo.inspect_types(pretty=True)

        # Exercise but supply a not None file kwarg, this is invalid
        with self.assertRaises(ValueError) as raises:
            foo.inspect_types(pretty=True, file='should be None')
        self.assertIn('`file` must be None if `pretty=True`',
                      str(raises.exception))


class TestTypeAnnotation(unittest.TestCase):

    def findpatloc(self, lines, pat):
        for i, ln in enumerate(lines):
            if pat in ln:
                return i
        raise ValueError("can't find {!r}".format(pat))

    def getlines(self, func):
        strbuf = StringIO()
        func.inspect_types(strbuf)
        return strbuf.getvalue().splitlines()

    def test_delete(self):
        @numba.njit
        def foo(appleorange, berrycherry):
            return appleorange + berrycherry

        foo(1, 2)

        lines = self.getlines(foo)

        # Ensure deletion show up after their use
        sa = self.findpatloc(lines, 'appleorange = arg(0, name=appleorange)')
        sb = self.findpatloc(lines, 'berrycherry = arg(1, name=berrycherry)')

        ea = self.findpatloc(lines, 'del appleorange')
        eb = self.findpatloc(lines, 'del berrycherry')

        self.assertLess(sa, ea)
        self.assertLess(sb, eb)

    def _lifetimes_impl(self, extend):
        with override_config('EXTEND_VARIABLE_LIFETIMES', extend):
            @njit
            def foo(a):
                b = a
                return b
            x = 10
            b = foo(x)
            self.assertEqual(b, x)

        lines = self.getlines(foo)

        sa = self.findpatloc(lines, 'a = arg(0, name=a)')
        sb = self.findpatloc(lines, 'b = a')

        cast_ret = self.findpatloc(lines, 'cast(value=b)')

        dela = self.findpatloc(lines, 'del a')
        delb = self.findpatloc(lines, 'del b')

        return sa, sb, cast_ret, dela, delb

    def test_delete_standard_lifetimes(self):
        # without extended lifetimes, dels occur as soon as dead
        #
        # label 0
        #   a = arg(0, name=a)  :: int64
        #   b = a  :: int64
        #   del a
        #   $8return_value.2 = cast(value=b)  :: int64
        #   del b
        #   return $8return_value.2

        sa, sb, cast_ret, dela, delb = self._lifetimes_impl(extend=0)

        self.assertLess(sa, dela)
        self.assertLess(sb, delb)
        # del a is before cast and del b is after
        self.assertLess(dela, cast_ret)
        self.assertGreater(delb, cast_ret)

    def test_delete_extended_lifetimes(self):
        # with extended lifetimes, dels are last in block:
        #
        # label 0
        #   a = arg(0, name=a)  :: int64
        #   b = a  :: int64
        #   $8return_value.2 = cast(value=b)  :: int64
        #   del a
        #   del b
        #   return $8return_value.2

        sa, sb, cast_ret, dela, delb = self._lifetimes_impl(extend=1)

        self.assertLess(sa, dela)
        self.assertLess(sb, delb)
        # dels are after the cast
        self.assertGreater(dela, cast_ret)
        self.assertGreater(delb, cast_ret)


if __name__ == '__main__':
    unittest.main()
