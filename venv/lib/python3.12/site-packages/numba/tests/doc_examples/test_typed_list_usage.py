# Contents in this file are referenced from the sphinx-generated docs.
# "magictoken" is used for markers as beginning and ending of example text.

import unittest
from numba.tests.support import captured_stdout


class DocsTypedListUsageTest(unittest.TestCase):

    def test_ex_inferred_list_jit(self):
        with captured_stdout():

            # magictoken.ex_inferred_list_jit.begin
            from numba import njit
            from numba.typed import List

            @njit
            def foo():
                # Instantiate a typed-list
                l = List()
                # Append a value to it, this will set the type to int32/int64
                # (depending on platform)
                l.append(42)
                # The usual list operations, getitem, pop and length are
                # supported
                print(l[0])   # 42
                l[0] = 23
                print(l[0])   # 23
                print(len(l)) # 1
                l.pop()
                print(len(l)) # 0
                return l

            foo()

            # magictoken.ex_inferred_list_jit.end

    def test_ex_inferred_list(self):
        with captured_stdout():
            # magictoken.ex_inferred_list.begin
            from numba import njit
            from numba.typed import List

            @njit
            def foo(mylist):
                for i in range(10, 20):
                    mylist.append(i)
                return mylist

            # Instantiate a typed-list, outside of a jit context
            l = List()
            # Append a value to it, this will set the type to int32/int64
            # (depending on platform)
            l.append(42)
            # The usual list operations, getitem, pop and length are supported
            print(l[0])   # 42
            l[0] = 23
            print(l[0])   # 23
            print(len(l)) # 1
            l.pop()
            print(len(l)) # 0

            # And you can use the typed-list as an argument for a jit compiled
            # function
            l = foo(l)
            print(len(l)) # 10

            # You can also directly construct a typed-list from an existing
            # Python list
            py_list = [2, 3, 5]
            numba_list = List(py_list)
            print(len(numba_list)) # 3

            # magictoken.ex_inferred_list.end

    def test_ex_nested_list(self):
        with captured_stdout():
            # magictoken.ex_nested_list.begin
            from numba.typed import List

            # typed-lists can be nested in typed-lists
            mylist = List()
            for i in range(10):
                l = List()
                for i in range(10):
                    l.append(i)
                mylist.append(l)
            # mylist is now a list of 10 lists, each containing 10 integers
            print(mylist)

            # magictoken.ex_nested_list.end


if __name__ == '__main__':
    unittest.main()
