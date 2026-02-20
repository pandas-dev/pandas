# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import captured_stdout, override_config


class DocsLLVMPassTimings(unittest.TestCase):

    def test_pass_timings(self):
        with override_config('LLVM_PASS_TIMINGS', True):
            with captured_stdout() as stdout:
                # magictoken.ex_llvm_pass_timings.begin
                import numba

                @numba.njit
                def foo(n):
                    c = 0
                    for i in range(n):
                        for j in range(i):
                            c += j
                    return c

                foo(10)
                md = foo.get_metadata(foo.signatures[0])
                print(md['llvm_pass_timings'])
                # magictoken.ex_llvm_pass_timings.end
            self.assertIn("Finalize object", stdout.getvalue())


if __name__ == '__main__':
    unittest.main()
