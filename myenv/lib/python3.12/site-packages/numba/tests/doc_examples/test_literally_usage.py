# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import captured_stdout


class DocsLiterallyUsageTest(unittest.TestCase):

    def test_literally_usage(self):
        with captured_stdout() as stdout:
            # magictoken.ex_literally_usage.begin
            import numba

            def power(x, n):
                raise NotImplementedError

            @numba.extending.overload(power)
            def ov_power(x, n):
                if isinstance(n, numba.types.Literal):
                    # only if `n` is a literal
                    if n.literal_value == 2:
                        # special case: square
                        print("square")
                        return lambda x, n: x * x
                    elif n.literal_value == 3:
                        # special case: cubic
                        print("cubic")
                        return lambda x, n: x * x * x
                else:
                    # If `n` is not literal, request literal dispatch
                    return lambda x, n: numba.literally(n)

                print("generic")
                return lambda x, n: x ** n

            @numba.njit
            def test_power(x, n):
                return power(x, n)

            # should print "square" and "9"
            print(test_power(3, 2))

            # should print "cubic" and "27"
            print(test_power(3, 3))

            # should print "generic" and "81"
            print(test_power(3, 4))

            # magictoken.ex_literally_usage.end
            assert test_power(3, 2) == 3 ** 2
            assert test_power(3, 3) == 3 ** 3
            assert test_power(3, 4) == 3 ** 4

        self.assertEqual('square\n9\ncubic\n27\ngeneric\n81\n',
                         stdout.getvalue())


if __name__ == '__main__':
    unittest.main()
