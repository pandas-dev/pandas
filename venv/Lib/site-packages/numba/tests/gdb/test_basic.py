# NOTE: This test is sensitive to line numbers as it checks breakpoints
from numba import njit, types
from numba.tests.gdb_support import GdbMIDriver
from numba.tests.support import TestCase, needs_subprocess
import unittest


@needs_subprocess
class Test(TestCase):

    def test(self):
        @njit(debug=True)
        def foo(x):
            z = 7 + x # break here
            return x, z

        foo(120)

        sz = types.intp.bitwidth
        driver = GdbMIDriver(__file__)
        driver.set_breakpoint(line=14)
        driver.run()
        driver.check_hit_breakpoint(1)
        driver.stack_list_arguments(2)
        expect = ('[frame={level="0",args=[{name="x",type="int%s",'
                  'value="120"}]}]' % sz)
        driver.assert_output(expect)
        driver.stack_list_variables(1)
        expect = '[{name="x",arg="1",value="120"},{name="z",value="0"}]'
        driver.assert_output(expect)
        driver.next()
        driver.stack_list_variables(1)
        expect = '[{name="x",arg="1",value="120"},{name="z",value="127"}]'
        driver.assert_output(expect)
        driver.quit()


if __name__ == '__main__':
    unittest.main()
