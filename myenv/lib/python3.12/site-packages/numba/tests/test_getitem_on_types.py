import unittest
from itertools import product
from numba import types, njit, typed, errors
from numba.tests.support import TestCase


class TestGetitemOnTypes(TestCase):
    # Tests getitem on the type types.

    def test_static_getitem_on_type(self):

        def gen(numba_type, index):
            def foo():
                ty = numba_type[index] # a static_getitem
                return typed.List.empty_list(ty)
            return foo

        # test a few types
        tys = (types.bool_, types.float64, types.uint8, types.complex128,)

        # and a few indexes of increasing complexity
        contig = slice(None, None, 1) # unit stride
        noncontig = slice(None, None, None)
        indexes = (contig, # 1d contig -> C order
                   noncontig, # 1d non-contig -> A order
                   (noncontig, contig), # 2d C order
                   (contig, noncontig), # 2d F order
                   (noncontig, noncontig), # 2d A order
                   (noncontig, noncontig, contig), # 3d C order
                   (contig, noncontig, noncontig), # 3d F order
                   (noncontig, noncontig, noncontig), # 3d A order
                   )

        for ty, idx in product(tys, indexes):
            compilable = njit(gen(ty, idx))
            # check the type of the typed list returned matches the type
            # as constructed in the interpreter
            expected = ty[idx]
            # check execution
            self.assertEqual(compilable()._dtype, expected)
            got = compilable.nopython_signatures[0].return_type.dtype
            # check sig
            self.assertEqual(got, expected)

    def test_shorthand_syntax(self):
        # tests a couple of shorthand syntax examples
        # (test_static_getitem_on_type is a more extensive test of the
        # functionality but it uses slices directly).

        @njit
        def foo1():
            ty = types.float32[::1, :] # 2d F order
            return typed.List.empty_list(ty)

        self.assertEqual(foo1()._dtype, types.float32[::1, :])

        @njit
        def foo2():
            ty = types.complex64[:, :, :] # 3d A order
            return typed.List.empty_list(ty)

        self.assertEqual(foo2()._dtype, types.complex64[:, :, :])

    def test_static_getitem_on_invalid_type(self):
        # check that an unsupported type cannot be instantiated in njit code

        # check this executes in the interpreter:
        types.void[:]

        # check the same fails in compilation as it's not supported
        # it'll fall back to a generic getitem
        with self.assertRaises(errors.TypingError) as raises:
            @njit
            def foo():
                types.void[:]

            foo()

        msg = ("No implementation",
               "getitem(typeref[none], slice<a:b>)")

        excstr = str(raises.exception)
        for m in msg:
            self.assertIn(m, excstr)

    def test_standard_getitem_on_type(self):
        # not supported at present, should be doable if the slice is a literal
        # though.

        # check using a non-static arg to the getitem raises
        with self.assertRaises(errors.TypingError) as raises:
            @njit
            def foo(not_static):
                types.float64[not_static]

            foo(slice(None, None, 1))

        msg = ("No implementation",
               "getitem(class(float64), slice<a:b>)")

        excstr = str(raises.exception)
        for m in msg:
            self.assertIn(m, excstr)


if __name__ == '__main__':
    unittest.main()
