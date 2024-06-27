# NOTE: This test is sensitive to line numbers as it checks breakpoints
from numba import njit
from numba.tests.gdb_support import GdbMIDriver
from numba.tests.support import TestCase, needs_subprocess
import unittest


@needs_subprocess
class Test(TestCase):

    def test(self):

        @njit(debug=True)
        def foo(x, y):
            c = x + y  # break-here
            return c

        @njit(debug=True)
        def call_foo(a):
            acc = 0
            for i in range(10):
                acc += foo(i, a)
            return acc

        call_foo(10)

        driver = GdbMIDriver(__file__)
        driver.set_breakpoint(line=15, condition='x == 4')
        driver.run()
        driver.check_hit_breakpoint(1)
        driver.stack_list_arguments(1)
        expect = ('[frame={level="0",args=[{name="x",value="4"},'
                  '{name="y",value="10"}]}]')
        driver.assert_output(expect)
        driver.set_breakpoint(line=22, condition='i == 8')
        driver.cont()
        driver.check_hit_breakpoint(2)
        driver.stack_list_variables(1)
        # i should be 8
        driver.assert_output('{name="i",value="8"}')
        driver.quit()


if __name__ == '__main__':
    unittest.main()
