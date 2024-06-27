import numpy as np
import sys
import traceback

from numba import jit, njit
from numba.core import types, errors, utils
from numba.tests.support import (TestCase, expected_failure_py311,
                                 expected_failure_py312,
                                 )
import unittest


force_pyobj_flags = {'nopython': False, 'forceobj': True}
no_pyobj_flags = {'nopython': True, '_nrt': False}
no_pyobj_flags_w_nrt = {'nopython': True, '_nrt': True}
no_gil_flags = {'nopython': True, 'nogil': True, '_nrt': True}


class MyError(Exception):
    pass


class OtherError(Exception):
    pass


class UDEArgsToSuper(Exception):
    def __init__(self, arg, value0):
        super(UDEArgsToSuper, self).__init__(arg)
        self.value0 = value0

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        same = True
        same |= self.args == other.args
        same |= self.value0 == other.value0
        return same

    def __hash__(self):
        return hash((super(UDEArgsToSuper).__hash__(), self.value0))


class UDENoArgSuper(Exception):
    def __init__(self, arg, value0):
        super(UDENoArgSuper, self).__init__()
        self.deferarg = arg
        self.value0 = value0

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        same = True
        same |= self.args == other.args
        same |= self.deferarg == other.deferarg
        same |= self.value0 == other.value0
        return same

    def __hash__(self):
        return hash((super(UDENoArgSuper).__hash__(), self.deferarg,
                     self.value0))


def raise_class(exc):
    def raiser(i):
        if i == 1:
            raise exc
        elif i == 2:
            raise ValueError
        elif i == 3:
            # The exception type is looked up on a module (issue #1624)
            raise np.linalg.LinAlgError
        return i
    return raiser


def raise_instance(exc, arg):
    def raiser(i):
        if i == 1:
            raise exc(arg, 1)
        elif i == 2:
            raise ValueError(arg, 2)
        elif i == 3:
            raise np.linalg.LinAlgError(arg, 3)
        return i
    return raiser


def raise_instance_runtime_args(exc):
    def raiser(i, arg):
        if i == 1:
            raise exc(arg, 1)
        elif i == 2:
            raise ValueError(arg, 2)
        elif i == 3:
            raise np.linalg.LinAlgError(arg, 3)
        return i
    return raiser


def reraise():
    raise


def outer_function(inner):
    def outer(i):
        if i == 3:
            raise OtherError("bar", 3)
        return inner(i)
    return outer


def assert_usecase(i):
    assert i == 1, "bar"


def ude_bug_usecase():
    raise UDEArgsToSuper()  # oops user forgot args to exception ctor


def raise_runtime_value(arg):
    raise ValueError(arg)


class TestRaising(TestCase):

    def test_unituple_index_error(self):
        def pyfunc(a, i):
            return a.shape[i]

        cfunc = njit((types.Array(types.int32, 1, 'A'), types.int32),)(pyfunc)

        a = np.empty(2, dtype=np.int32)

        self.assertEqual(cfunc(a, 0), pyfunc(a, 0))

        with self.assertRaises(IndexError) as cm:
            cfunc(a, 2)
        self.assertEqual(str(cm.exception), "tuple index out of range")

    def check_against_python(self, exec_mode, pyfunc, cfunc,
                             expected_error_class, *args):

        assert exec_mode in (force_pyobj_flags, no_pyobj_flags,
                             no_pyobj_flags_w_nrt, no_gil_flags)

        # invariant of mode, check the error class and args are the same
        with self.assertRaises(expected_error_class) as pyerr:
            pyfunc(*args)
        with self.assertRaises(expected_error_class) as jiterr:
            cfunc(*args)
        self.assertEqual(pyerr.exception.args, jiterr.exception.args)

        # special equality check for UDEs
        if isinstance(pyerr.exception, (UDEArgsToSuper, UDENoArgSuper)):
            self.assertTrue(pyerr.exception == jiterr.exception)

        # in npm check bottom of traceback matches as frame injection with
        # location info should ensure this
        if exec_mode is no_pyobj_flags:

            # we only care about the bottom two frames, the error and the
            # location it was raised.
            try:
                pyfunc(*args)
            except Exception:
                py_frames = traceback.format_exception(*sys.exc_info())
                expected_frames = py_frames[-2:]

            try:
                cfunc(*args)
            except Exception:
                c_frames = traceback.format_exception(*sys.exc_info())
                got_frames = c_frames[-2:]

            # check exception and the injected frame are the same
            for expf, gotf in zip(expected_frames, got_frames):
                # Note use of assertIn not assertEqual, Py 3.11 has markers (^)
                # that point to the variable causing the problem, Numba doesn't
                # do this so only the start of the string will match.
                self.assertIn(gotf, expf)

    def check_raise_class(self, flags):
        pyfunc = raise_class(MyError)
        cfunc = jit((types.int32,), **flags)(pyfunc)
        self.assertEqual(cfunc(0), 0)
        self.check_against_python(flags, pyfunc, cfunc, MyError, 1)
        self.check_against_python(flags, pyfunc, cfunc, ValueError, 2)
        self.check_against_python(flags, pyfunc, cfunc,
                                  np.linalg.linalg.LinAlgError, 3)

    def test_raise_class_nopython(self):
        self.check_raise_class(flags=no_pyobj_flags)

    def test_raise_class_objmode(self):
        self.check_raise_class(flags=force_pyobj_flags)

    def check_raise_instance(self, flags):
        for clazz in [MyError, UDEArgsToSuper,
                      UDENoArgSuper]:
            pyfunc = raise_instance(clazz, "some message")
            cfunc = jit((types.int32,), **flags)(pyfunc)

            self.assertEqual(cfunc(0), 0)
            self.check_against_python(flags, pyfunc, cfunc, clazz, 1)
            self.check_against_python(flags, pyfunc, cfunc, ValueError, 2)
            self.check_against_python(flags, pyfunc, cfunc,
                                      np.linalg.linalg.LinAlgError, 3)

    def test_raise_instance_objmode(self):
        self.check_raise_instance(flags=force_pyobj_flags)

    def test_raise_instance_nopython(self):
        self.check_raise_instance(flags=no_pyobj_flags)

    def check_raise_nested(self, flags, **jit_args):
        """
        Check exception propagation from nested functions.
        """
        for clazz in [MyError, UDEArgsToSuper,
                      UDENoArgSuper]:
            inner_pyfunc = raise_instance(clazz, "some message")
            pyfunc = outer_function(inner_pyfunc)
            inner_cfunc = jit(**jit_args)(inner_pyfunc)
            cfunc = jit(**jit_args)(outer_function(inner_cfunc))

            self.check_against_python(flags, pyfunc, cfunc, clazz, 1)
            self.check_against_python(flags, pyfunc, cfunc, ValueError, 2)
            self.check_against_python(flags, pyfunc, cfunc, OtherError, 3)

    def test_raise_nested_objmode(self):
        self.check_raise_nested(force_pyobj_flags, forceobj=True)

    def test_raise_nested_nopython(self):
        self.check_raise_nested(no_pyobj_flags, nopython=True)

    def check_reraise(self, flags):
        def raise_exc(exc):
            raise exc
        pyfunc = reraise
        cfunc = jit((), **flags)(pyfunc)
        for op, err in [(lambda : raise_exc(ZeroDivisionError),
                         ZeroDivisionError),
                        (lambda : raise_exc(UDEArgsToSuper("msg", 1)),
                         UDEArgsToSuper),
                        (lambda : raise_exc(UDENoArgSuper("msg", 1)),
                         UDENoArgSuper)]:
            def gen_impl(fn):
                def impl():
                    try:
                        op()
                    except err:
                        fn()
                return impl
            pybased = gen_impl(pyfunc)
            cbased = gen_impl(cfunc)
            self.check_against_python(flags, pybased, cbased, err,)

    def test_reraise_objmode(self):
        self.check_reraise(flags=force_pyobj_flags)

    def test_reraise_nopython(self):
        self.check_reraise(flags=no_pyobj_flags)

    def check_raise_invalid_class(self, cls, flags):
        pyfunc = raise_class(cls)
        cfunc = jit((types.int32,), **flags)(pyfunc)
        with self.assertRaises(TypeError) as cm:
            cfunc(1)
        self.assertEqual(str(cm.exception),
                         "exceptions must derive from BaseException")

    def test_raise_invalid_class_objmode(self):
        self.check_raise_invalid_class(int, flags=force_pyobj_flags)
        self.check_raise_invalid_class(1, flags=force_pyobj_flags)

    def test_raise_invalid_class_nopython(self):
        msg = "Encountered unsupported constant type used for exception"
        with self.assertRaises(errors.UnsupportedError) as raises:
            self.check_raise_invalid_class(int, flags=no_pyobj_flags)
        self.assertIn(msg, str(raises.exception))
        with self.assertRaises(errors.UnsupportedError) as raises:
            self.check_raise_invalid_class(1, flags=no_pyobj_flags)
        self.assertIn(msg, str(raises.exception))

    def test_raise_bare_string_nopython(self):
        @njit
        def foo():
            raise "illegal"
        msg = ("Directly raising a string constant as an exception is not "
               "supported")
        with self.assertRaises(errors.UnsupportedError) as raises:
            foo()
        self.assertIn(msg, str(raises.exception))

    def check_assert_statement(self, flags):
        pyfunc = assert_usecase
        cfunc = jit((types.int32,), **flags)(pyfunc)
        cfunc(1)
        self.check_against_python(flags, pyfunc, cfunc, AssertionError, 2)

    def test_assert_statement_objmode(self):
        self.check_assert_statement(flags=force_pyobj_flags)

    def test_assert_statement_nopython(self):
        self.check_assert_statement(flags=no_pyobj_flags)

    def check_raise_from_exec_string(self, flags):
        # issue #3428
        simple_raise = "def f(a):\n  raise exc('msg', 10)"
        assert_raise = "def f(a):\n  assert a != 1"
        py312_pep695_raise = "def f[T: int](a: T) -> T:\n  assert a != 1"
        py312_pep695_raise_2 = "def f[T: int\n](a: T) -> T:\n  assert a != 1"
        test_cases = [
            (assert_raise, AssertionError),
            (simple_raise, UDEArgsToSuper),
            (simple_raise, UDENoArgSuper),
        ]
        if utils.PYVERSION >= (3, 12):
            # Added for https://github.com/numba/numba/issues/9443
            test_cases.append((py312_pep695_raise, AssertionError))
            test_cases.append((py312_pep695_raise_2, AssertionError))
        for f_text, exc in test_cases:
            loc = {}
            exec(f_text, {'exc': exc}, loc)
            pyfunc = loc['f']
            cfunc = jit((types.int32,), **flags)(pyfunc)
            self.check_against_python(flags, pyfunc, cfunc, exc, 1)

    def test_assert_from_exec_string_objmode(self):
        self.check_raise_from_exec_string(flags=force_pyobj_flags)

    def test_assert_from_exec_string_nopython(self):
        self.check_raise_from_exec_string(flags=no_pyobj_flags)

    def check_user_code_error_traceback(self, flags):
        # this test checks that if a user tries to compile code that contains
        # a bug in exception initialisation (e.g. missing arg) then this also
        # has a frame injected with the location information.
        pyfunc = ude_bug_usecase
        cfunc = jit((), **flags)(pyfunc)
        self.check_against_python(flags, pyfunc, cfunc, TypeError)

    def test_user_code_error_traceback_objmode(self):
        self.check_user_code_error_traceback(flags=force_pyobj_flags)

    def test_user_code_error_traceback_nopython(self):
        self.check_user_code_error_traceback(flags=no_pyobj_flags)

    def check_raise_runtime_value(self, flags):
        pyfunc = raise_runtime_value
        cfunc = jit((types.string,), **flags)(pyfunc)
        self.check_against_python(flags, pyfunc, cfunc, ValueError, 'hello')

    def test_raise_runtime_value_objmode(self):
        self.check_raise_runtime_value(flags=force_pyobj_flags)

    def test_raise_runtime_value_nopython(self):
        self.check_raise_runtime_value(flags=no_pyobj_flags_w_nrt)

    def test_raise_runtime_value_nogil(self):
        self.check_raise_runtime_value(flags=no_gil_flags)

    def check_raise_instance_with_runtime_args(self, flags):
        for clazz in [MyError, UDEArgsToSuper,
                      UDENoArgSuper]:
            pyfunc = raise_instance_runtime_args(clazz)
            cfunc = jit((types.int32, types.string), **flags)(pyfunc)

            self.assertEqual(cfunc(0, 'test'), 0)
            self.check_against_python(flags, pyfunc, cfunc, clazz, 1, 'hello')
            self.check_against_python(flags, pyfunc, cfunc, ValueError, 2,
                                      'world')
            self.check_against_python(flags, pyfunc, cfunc,
                                      np.linalg.linalg.LinAlgError, 3, 'linalg')

    def test_raise_instance_with_runtime_args_objmode(self):
        self.check_raise_instance_with_runtime_args(flags=force_pyobj_flags)

    def test_raise_instance_with_runtime_args_nopython(self):
        self.check_raise_instance_with_runtime_args(flags=no_pyobj_flags_w_nrt)

    def test_raise_instance_with_runtime_args_nogil(self):
        self.check_raise_instance_with_runtime_args(flags=no_gil_flags)

    def test_dynamic_raise_bad_args(self):
        def raise_literal_dict():
            raise ValueError({'a': 1, 'b': np.ones(4)})

        def raise_range():
            raise ValueError(range(3))

        def raise_rng(rng):
            raise ValueError(rng.bit_generator)

        funcs = [
            (raise_literal_dict, ()),
            (raise_range, ()),
            (raise_rng, (types.npy_rng,)),
        ]

        for pyfunc, argtypes in funcs:
            msg = '.*Cannot convert native .* to a Python object.*'
            with self.assertRaisesRegex(errors.TypingError, msg):
                njit(argtypes)(pyfunc)

    def test_dynamic_raise_dict(self):
        @njit
        def raise_literal_dict2():
            raise ValueError({'a': 1, 'b': 3})

        msg = "{a: 1, b: 3}"
        with self.assertRaisesRegex(ValueError, msg):
            raise_literal_dict2()

    def test_disable_nrt(self):
        @njit(_nrt=False)
        def raise_with_no_nrt(i):
            raise ValueError(i)

        msg = 'NRT required but not enabled'
        with self.assertRaisesRegex(errors.NumbaRuntimeError, msg):
            raise_with_no_nrt(123)

    def test_try_raise(self):

        @njit
        def raise_(a):
            raise ValueError(a)

        @njit
        def try_raise(a):
            try:
                raise_(a)
            except Exception:
                pass
            return a + 1

        self.assertEqual(try_raise.py_func(3), try_raise(3))

    @expected_failure_py311
    @expected_failure_py312
    def test_dynamic_raise(self):

        @njit
        def raise_(a):
            raise ValueError(a)

        @njit
        def try_raise_(a):
            try:
                raise_(a)
            except Exception:
                raise ValueError(a)

        args = [
            1,
            1.1,
            'hello',
            np.ones(3),
            [1, 2],
            (1, 2),
            set([1, 2]),
        ]
        for fn in (raise_, try_raise_):
            for arg in args:
                with self.assertRaises(ValueError) as e:
                    fn(arg)
                self.assertEqual((arg,), e.exception.args)


if __name__ == '__main__':
    unittest.main()
