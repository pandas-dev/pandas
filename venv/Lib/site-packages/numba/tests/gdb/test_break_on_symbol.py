# NOTE: This test is sensitive to line numbers as it checks breakpoints
from numba import njit, types
from numba.tests.gdb_support import GdbMIDriver
from numba.tests.support import TestCase, needs_subprocess
import unittest


@njit(debug=True)
def foo(x):
    z = 7 + x
    return x, z


@needs_subprocess
class Test(TestCase):

    def test(self):
        foo(120)
        sz = types.intp.bitwidth
        driver = GdbMIDriver(__file__)
        driver.set_breakpoint(symbol="__main__::foo")
        driver.run() # will hit cpython symbol match
        driver.check_hit_breakpoint(number=1)
        driver.cont() # will hit njit symbol match
        driver.check_hit_breakpoint(number=1, line=10) # Ensure line number
        driver.stack_list_arguments(2)
        expect = ('[frame={level="0",args=[{name="x",type="int%s",'
                  'value="120"}]}]' % sz)
        driver.assert_output(expect)
        driver.quit()


if __name__ == '__main__':
    unittest.main()
