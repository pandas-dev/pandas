# NOTE: This test is sensitive to line numbers as it checks breakpoints
from numba import njit
import numpy as np
from numba.tests.gdb_support import GdbMIDriver, needs_gdb_py3
from numba.tests.support import TestCase, needs_subprocess
from numba.misc.numba_gdbinfo import collect_gdbinfo
import unittest
import re


@needs_gdb_py3
@needs_subprocess
class Test(TestCase):

    def test(self):
        rdt_a = np.dtype([("x", np.int16), ("y", np.float64)], align=True)

        @njit(debug=True)
        def foo():
            a = 1.234
            b = (1, 2, 3)
            c = ('a', b, 4)
            d = np.arange(5.)
            e = np.array([[1, 3j], [2, 4j]])
            f = "Some string" + "           L-Padded string".lstrip()
            g = 11 + 22j
            h = np.arange(24).reshape((4, 6))[::2, ::3]
            i = np.zeros(2, dtype=rdt_a)
            return a, b, c, d, e, f, g, h, i

        foo()

        extension = collect_gdbinfo().extension_loc
        driver = GdbMIDriver(__file__, init_cmds=['-x', extension], debug=False)
        driver.set_breakpoint(line=29)
        driver.run()
        driver.check_hit_breakpoint(1)

        # Ideally the function would be run to get the string repr of locals
        # but not everything appears in DWARF e.g. string literals. Further,
        # str on NumPy arrays seems to vary a bit in output. Therefore a custom
        # match is used.

        driver.stack_list_variables(1)
        output = driver._captured.after.decode('UTF-8')
        done_str = output.splitlines()[0]
        pat = r'^\^done,variables=\[\{(.*)\}\]$'
        lcls_strs = re.match(pat, done_str).groups()[0].split('},{')
        lcls = {k: v for k, v in [re.match(r'name="(.*)",value="(.*)"',
                x).groups() for x in lcls_strs]}
        expected = dict()
        expected['a'] = r'1\.234'
        expected['b'] = r'\(1, 2, 3\)'
        expected['c'] = r'\(0x0, \(1, 2, 3\), 4\)'
        expected['d'] = r'\\n\[0. 1. 2. 3. 4.\]'
        expected['e'] = r'\\n\[\[1.\+0.j 0.\+3.j\]\\n \[2.\+0.j 0.\+4.j\]\]'
        expected['f'] = "'Some stringL-Padded string'"
        expected['g'] = r"11\+22j"
        expected['h'] = r'\\n\[\[ 0  3\]\\n \[12 15\]\]'
        expected['i'] = r'\\n\[\(0, 0.\) \(0, 0.\)\]'

        for k, v in expected.items():
            self.assertRegex(lcls[k], v)

        driver.quit()


if __name__ == '__main__':
    unittest.main()
