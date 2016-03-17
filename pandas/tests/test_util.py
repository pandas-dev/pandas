# -*- coding: utf-8 -*-
import nose

from pandas.util._move import move_into_mutable_buffer, BadMove
from pandas.util.decorators import deprecate_kwarg
from pandas.util.validators import validate_args, validate_kwargs

import pandas.util.testing as tm


class TestDecorators(tm.TestCase):

    def setUp(self):
        @deprecate_kwarg('old', 'new')
        def _f1(new=False):
            return new

        @deprecate_kwarg('old', 'new', {'yes': True, 'no': False})
        def _f2(new=False):
            return new

        @deprecate_kwarg('old', 'new', lambda x: x + 1)
        def _f3(new=0):
            return new

        self.f1 = _f1
        self.f2 = _f2
        self.f3 = _f3

    def test_deprecate_kwarg(self):
        x = 78
        with tm.assert_produces_warning(FutureWarning):
            result = self.f1(old=x)
        self.assertIs(result, x)
        with tm.assert_produces_warning(None):
            self.f1(new=x)

    def test_dict_deprecate_kwarg(self):
        x = 'yes'
        with tm.assert_produces_warning(FutureWarning):
            result = self.f2(old=x)
        self.assertEqual(result, True)

    def test_missing_deprecate_kwarg(self):
        x = 'bogus'
        with tm.assert_produces_warning(FutureWarning):
            result = self.f2(old=x)
        self.assertEqual(result, 'bogus')

    def test_callable_deprecate_kwarg(self):
        x = 5
        with tm.assert_produces_warning(FutureWarning):
            result = self.f3(old=x)
        self.assertEqual(result, x + 1)
        with tm.assertRaises(TypeError):
            self.f3(old='hello')

    def test_bad_deprecate_kwarg(self):
        with tm.assertRaises(TypeError):
            @deprecate_kwarg('old', 'new', 0)
            def f4(new=None):
                pass


def test_rands():
    r = tm.rands(10)
    assert(len(r) == 10)


def test_rands_array():
    arr = tm.rands_array(5, size=10)
    assert(arr.shape == (10,))
    assert(len(arr[0]) == 5)

    arr = tm.rands_array(7, size=(10, 10))
    assert(arr.shape == (10, 10))
    assert(len(arr[1, 1]) == 7)


class TestValidateArgs(tm.TestCase):

    def test_bad_min_length(self):
        msg = "'min_length' must be non-negative"
        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args((None,), min_length=-1, max_length=5)

    def test_bad_arg_length_no_max(self):
        min_length = 5
        msg = "expected at least {min_length} arguments".format(
            min_length=min_length)

        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args((None,), min_length=min_length, max_length=None)

    def test_bad_arg_length_with_max(self):
        min_length = 5
        max_length = 10
        msg = ("expected between {min_length} and {max_length}"
               " arguments inclusive".format(min_length=min_length,
                                             max_length=max_length))

        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args((None,), min_length=min_length,
                          max_length=max_length)

    def test_bad_min_max_length(self):
        msg = "'min_length' > 'max_length'"
        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args((None,), min_length=5, max_length=2)

    def test_not_all_none(self):
        msg = "All arguments must be None"
        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args(('foo',), min_length=0,
                          max_length=1, msg=msg)

        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args(('foo', 'bar', 'baz'), min_length=2,
                          max_length=5, msg=msg)

        with tm.assertRaisesRegexp(ValueError, msg):
            validate_args((None, 'bar', None), min_length=2,
                          max_length=5, msg=msg)

    def test_validation(self):
        # No exceptions should be thrown
        validate_args((None,), min_length=0, max_length=1)
        validate_args((None, None), min_length=1, max_length=5)


class TestValidateKwargs(tm.TestCase):

    def test_bad_kwarg(self):
        goodarg = 'f'
        badarg = goodarg + 'o'

        kwargs = {goodarg: 'foo', badarg: 'bar'}
        compat_args = (goodarg, badarg + 'o')
        fname = 'func'

        msg = ("{fname}\(\) got an unexpected "
               "keyword argument '{arg}'".format(
                   fname=fname, arg=badarg))

        with tm.assertRaisesRegexp(TypeError, msg):
            validate_kwargs(fname, kwargs, *compat_args)

    def test_validation(self):
        # No exceptions should be thrown
        compat_args = ('f', 'b', 'ba')
        kwargs = {'f': 'foo', 'b': 'bar'}
        validate_kwargs('func', kwargs, *compat_args)


class TestMove(tm.TestCase):
    def test_more_than_one_ref(self):
        """Test case for when we try to use ``move_into_mutable_buffer`` when
        the object being moved has other references.
        """
        b = b'testing'

        with tm.assertRaises(BadMove) as e:
            def handle_success(type_, value, tb):
                self.assertIs(value.args[0], b)
                return type(e).handle_success(e, type_, value, tb)  # super

            e.handle_success = handle_success
            move_into_mutable_buffer(b)

    def test_exactly_one_ref(self):
        """Test case for when the object being moved has exactly one reference.
        """
        b = b'testing'

        # We need to pass an expression on the stack to ensure that there are
        # not extra references hanging around. We cannot rewrite this test as
        #   buf = b[:-3]
        #   as_stolen_buf = move_into_mutable_buffer(buf)
        # because then we would have more than one reference to buf.
        as_stolen_buf = move_into_mutable_buffer(b[:-3])

        # materialize as bytearray to show that it is mutable
        self.assertEqual(bytearray(as_stolen_buf), b'test')


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
