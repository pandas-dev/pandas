# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import captured_stdout
from numba import typed


class DocsLiteralContainerUsageTest(unittest.TestCase):

    def test_ex_literal_dict_compile_time_consts(self):
        with captured_stdout():
            # magictoken.test_ex_literal_dict_compile_time_consts.begin
            import numpy as np
            from numba import njit, types
            from numba.extending import overload

            # overload this function
            def specialize(x):
                pass

            @overload(specialize)
            def ol_specialize(x):
                ld = x.literal_value
                const_expr = []
                for k, v in ld.items():
                    if isinstance(v, types.Literal):
                        lv = v.literal_value
                        if lv == 'cat':
                            const_expr.append("Meow!")
                        elif lv == 'dog':
                            const_expr.append("Woof!")
                        elif isinstance(lv, int):
                            const_expr.append(k.literal_value * lv)
                    else: # it's an array
                        const_expr.append("Array(dim={dim}".format(dim=v.ndim))
                const_strings = tuple(const_expr)

                def impl(x):
                    return const_strings
                return impl

            @njit
            def foo():
                pets_ints_and_array = {'a': 1,
                                       'b': 2,
                                       'c': 'cat',
                                       'd': 'dog',
                                       'e': np.ones(5,)}
                return specialize(pets_ints_and_array)

            result = foo()
            print(result) # ('a', 'bb', 'Meow!', 'Woof!', 'Array(dim=1')
            # magictoken.test_ex_literal_dict_compile_time_consts.end

        self.assertEqual(result, ('a', 'bb', 'Meow!', 'Woof!', 'Array(dim=1'))

    def test_ex_initial_value_dict_compile_time_consts(self):
        with captured_stdout():
            # magictoken.test_ex_initial_value_dict_compile_time_consts.begin
            from numba import njit, literally
            from numba.extending import overload

            # overload this function
            def specialize(x):
                pass

            @overload(specialize)
            def ol_specialize(x):
                iv = x.initial_value
                if iv is None:
                    return lambda x: literally(x) # Force literal dispatch
                assert iv == {'a': 1, 'b': 2, 'c': 3} # INITIAL VALUE
                return lambda x: literally(x)

            @njit
            def foo():
                d = {'a': 1, 'b': 2, 'c': 3}
                d['c'] = 20 # no impact on .initial_value
                d['d'] = 30 # no impact on .initial_value
                return specialize(d)

            result = foo()
            print(result) # {a: 1, b: 2, c: 20, d: 30} # NOT INITIAL VALUE!
            # magictoken.test_ex_initial_value_dict_compile_time_consts.end

        expected = typed.Dict()
        for k, v in {'a': 1, 'b': 2, 'c': 20, 'd': 30}.items():
            expected[k] = v
        self.assertEqual(result, expected)

    def test_ex_literal_list(self):
        with captured_stdout():
            # magictoken.test_ex_literal_list.begin
            from numba import njit
            from numba.extending import overload

            # overload this function
            def specialize(x):
                pass

            @overload(specialize)
            def ol_specialize(x):
                l = x.literal_value
                const_expr = []
                for v in l:
                    const_expr.append(str(v))
                const_strings = tuple(const_expr)

                def impl(x):
                    return const_strings
                return impl

            @njit
            def foo():
                const_list = ['a', 10, 1j, ['another', 'list']]
                return specialize(const_list)

            result = foo()
            print(result) # ('Literal[str](a)', 'Literal[int](10)', 'complex128', 'list(unicode_type)') # noqa E501
            # magictoken.test_ex_literal_list.end

        expected = ('Literal[str](a)', 'Literal[int](10)', 'complex128',
                    "list(unicode_type)<iv=['another', 'list']>")
        self.assertEqual(result, expected)

    def test_ex_initial_value_list_compile_time_consts(self):
        with captured_stdout():
            # magictoken.test_ex_initial_value_list_compile_time_consts.begin
            from numba import njit, literally
            from numba.extending import overload

            # overload this function
            def specialize(x):
                pass

            @overload(specialize)
            def ol_specialize(x):
                iv = x.initial_value
                if iv is None:
                    return lambda x: literally(x) # Force literal dispatch
                assert iv == [1, 2, 3] # INITIAL VALUE
                return lambda x: x

            @njit
            def foo():
                l = [1, 2, 3]
                l[2] = 20 # no impact on .initial_value
                l.append(30) # no impact on .initial_value
                return specialize(l)

            result = foo()
            print(result) # [1, 2, 20, 30] # NOT INITIAL VALUE!
            # magictoken.test_ex_initial_value_list_compile_time_consts.end

        expected = [1, 2, 20, 30]
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
