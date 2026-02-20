# NOTE: This test is sensitive to line numbers as it checks breakpoints
from numba import njit, types
import numpy as np
from numba.tests.gdb_support import GdbMIDriver
from numba.tests.support import TestCase, needs_subprocess
import unittest


@needs_subprocess
class Test(TestCase):

    def test(self):
        @njit(debug=True)
        def foo(x):
            z = np.ones_like(x) # break here
            return x, z

        tmp = np.ones(5)
        foo(tmp)

        driver = GdbMIDriver(__file__)
        driver.set_breakpoint(line=15)
        driver.run()
        driver.check_hit_breakpoint(1)
        driver.stack_list_arguments(2)
        llvm_intp = f"i{types.intp.bitwidth}"
        expect = (
            '[frame={level="0",args=[{name="x",type="array(float64, 1d, C) '
            f'({{i8*, i8*, {llvm_intp}, {llvm_intp}, double*, '
            f'[1 x {llvm_intp}], [1 x {llvm_intp}]}})"}}]}}]'
        )
        driver.assert_output(expect)
        driver.stack_list_variables(1)
        # 'z' should be zero-init
        expect = ('{name="z",value="{meminfo = 0x0, parent = 0x0, nitems = 0, '
                  'itemsize = 0, data = 0x0, shape = {0}, strides = {0}}"}')
        driver.assert_output(expect)
        driver.set_breakpoint(line=16)
        driver.cont()
        driver.check_hit_breakpoint(2)
        driver.stack_list_variables(1)
        # 'z' should be populated
        expect = (r'^.*\{name="z",value="\{meminfo = 0x[0-9a-f]+ .*, '
                  r'parent = 0x0, nitems = 5, itemsize = 8, '
                  r'data = 0x[0-9a-f]+, shape = \{5\}, strides = \{8\}\}.*$')
        driver.assert_regex_output(expect)
        driver.quit()


if __name__ == '__main__':
    unittest.main()
