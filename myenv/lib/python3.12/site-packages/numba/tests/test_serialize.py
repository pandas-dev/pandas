import contextlib
import gc
import pickle
import runpy
import subprocess
import sys
import unittest
from multiprocessing import get_context

import numba
from numba.core.errors import TypingError
from numba.tests.support import TestCase
from numba.core.target_extension import resolve_dispatcher_from_str
from numba.cloudpickle import dumps, loads


class TestDispatcherPickling(TestCase):

    def run_with_protocols(self, meth, *args, **kwargs):
        for proto in range(pickle.HIGHEST_PROTOCOL + 1):
            meth(proto, *args, **kwargs)

    @contextlib.contextmanager
    def simulate_fresh_target(self):
        hwstr = 'cpu'
        dispatcher_cls = resolve_dispatcher_from_str(hwstr)
        old_descr = dispatcher_cls.targetdescr
        # Simulate fresh targetdescr
        dispatcher_cls.targetdescr = type(dispatcher_cls.targetdescr)(hwstr)
        try:
            yield
        finally:
            # Be sure to reinstantiate old descriptor, otherwise other
            # objects may be out of sync.
            dispatcher_cls.targetdescr = old_descr

    def check_call(self, proto, func, expected_result, args):
        def check_result(func):
            if (isinstance(expected_result, type)
                and issubclass(expected_result, Exception)):
                self.assertRaises(expected_result, func, *args)
            else:
                self.assertPreciseEqual(func(*args), expected_result)

        # Control
        check_result(func)
        pickled = pickle.dumps(func, proto)
        with self.simulate_fresh_target():
            new_func = pickle.loads(pickled)
            check_result(new_func)

    def test_call_with_sig(self):
        from .serialize_usecases import add_with_sig
        self.run_with_protocols(self.check_call, add_with_sig, 5, (1, 4))
        # Compilation has been disabled => float inputs will be coerced to int
        self.run_with_protocols(self.check_call, add_with_sig, 5, (1.2, 4.2))

    def test_call_without_sig(self):
        from .serialize_usecases import add_without_sig
        self.run_with_protocols(self.check_call, add_without_sig, 5, (1, 4))
        self.run_with_protocols(self.check_call, add_without_sig, 5.5, (1.2, 4.3))
        # Object mode is enabled
        self.run_with_protocols(self.check_call, add_without_sig, "abc", ("a", "bc"))

    def test_call_nopython(self):
        from .serialize_usecases import add_nopython
        self.run_with_protocols(self.check_call, add_nopython, 5.5, (1.2, 4.3))
        # Object mode is disabled
        self.run_with_protocols(self.check_call, add_nopython, TypingError, (object(), object()))

    def test_call_nopython_fail(self):
        from .serialize_usecases import add_nopython_fail
        # Compilation fails
        self.run_with_protocols(self.check_call, add_nopython_fail, TypingError, (1, 2))

    def test_call_objmode_with_global(self):
        from .serialize_usecases import get_global_objmode
        self.run_with_protocols(self.check_call, get_global_objmode, 7.5, (2.5,))

    def test_call_closure(self):
        from .serialize_usecases import closure
        inner = closure(1)
        self.run_with_protocols(self.check_call, inner, 6, (2, 3))

    def check_call_closure_with_globals(self, **jit_args):
        from .serialize_usecases import closure_with_globals
        inner = closure_with_globals(3.0, **jit_args)
        self.run_with_protocols(self.check_call, inner, 7.0, (4.0,))

    def test_call_closure_with_globals_nopython(self):
        self.check_call_closure_with_globals(nopython=True)

    def test_call_closure_with_globals_objmode(self):
        self.check_call_closure_with_globals(forceobj=True)

    def test_call_closure_calling_other_function(self):
        from .serialize_usecases import closure_calling_other_function
        inner = closure_calling_other_function(3.0)
        self.run_with_protocols(self.check_call, inner, 11.0, (4.0, 6.0))

    def test_call_closure_calling_other_closure(self):
        from .serialize_usecases import closure_calling_other_closure
        inner = closure_calling_other_closure(3.0)
        self.run_with_protocols(self.check_call, inner, 8.0, (4.0,))

    def test_call_dyn_func(self):
        from .serialize_usecases import dyn_func
        # Check serializing a dynamically-created function
        self.run_with_protocols(self.check_call, dyn_func, 36, (6,))

    def test_call_dyn_func_objmode(self):
        from .serialize_usecases import dyn_func_objmode
        # Same with an object mode function
        self.run_with_protocols(self.check_call, dyn_func_objmode, 36, (6,))

    def test_renamed_module(self):
        from .serialize_usecases import get_renamed_module
        # Issue #1559: using a renamed module (e.g. `import numpy as np`)
        # should not fail serializing
        expected = get_renamed_module(0.0)
        self.run_with_protocols(self.check_call, get_renamed_module,
                                expected, (0.0,))

    def test_other_process(self):
        """
        Check that reconstructing doesn't depend on resources already
        instantiated in the original process.
        """
        from .serialize_usecases import closure_calling_other_closure
        func = closure_calling_other_closure(3.0)
        pickled = pickle.dumps(func)
        code = """if 1:
            import pickle

            data = {pickled!r}
            func = pickle.loads(data)
            res = func(4.0)
            assert res == 8.0, res
            """.format(**locals())
        subprocess.check_call([sys.executable, "-c", code])

    def test_reuse(self):
        """
        Check that deserializing the same function multiple times re-uses
        the same dispatcher object.

        Note that "same function" is intentionally under-specified.
        """
        from .serialize_usecases import closure
        func = closure(5)
        pickled = pickle.dumps(func)
        func2 = closure(6)
        pickled2 = pickle.dumps(func2)

        f = pickle.loads(pickled)
        g = pickle.loads(pickled)
        h = pickle.loads(pickled2)
        self.assertIs(f, g)
        self.assertEqual(f(2, 3), 10)
        g.disable_compile()
        self.assertEqual(g(2, 4), 11)

        self.assertIsNot(f, h)
        self.assertEqual(h(2, 3), 11)

        # Now make sure the original object doesn't exist when deserializing
        func = closure(7)
        func(42, 43)
        pickled = pickle.dumps(func)
        del func
        gc.collect()

        f = pickle.loads(pickled)
        g = pickle.loads(pickled)
        self.assertIs(f, g)
        self.assertEqual(f(2, 3), 12)
        g.disable_compile()
        self.assertEqual(g(2, 4), 13)

    def test_imp_deprecation(self):
        """
        The imp module was deprecated in v3.4 in favour of importlib
        """
        code = """if 1:
            import pickle
            import warnings
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always', DeprecationWarning)
                from numba import njit
                @njit
                def foo(x):
                    return x + 1
                foo(1)
                serialized_foo = pickle.dumps(foo)
            for x in w:
                if 'serialize.py' in x.filename:
                    assert "the imp module is deprecated" not in x.msg
        """
        subprocess.check_call([sys.executable, "-c", code])


class TestSerializationMisc(TestCase):
    def test_numba_unpickle(self):
        # Test that _numba_unpickle is memorizing its output
        from numba.core.serialize import _numba_unpickle

        random_obj = object()
        bytebuf = pickle.dumps(random_obj)
        hashed = hash(random_obj)

        got1 = _numba_unpickle(id(random_obj), bytebuf, hashed)
        # not the original object
        self.assertIsNot(got1, random_obj)
        got2 = _numba_unpickle(id(random_obj), bytebuf, hashed)
        # unpickled results are the same objects
        self.assertIs(got1, got2)


class TestCloudPickleIssues(TestCase):
    """This test case includes issues specific to the cloudpickle implementation.
    """
    _numba_parallel_test_ = False

    def test_dynamic_class_reset_on_unpickle(self):
        # a dynamic class
        class Klass:
            classvar = None

        def mutator():
            Klass.classvar = 100

        def check():
            self.assertEqual(Klass.classvar, 100)

        saved = dumps(Klass)
        mutator()
        check()
        loads(saved)
        # Without the patch, each `loads(saved)` will reset `Klass.classvar`
        check()
        loads(saved)
        check()

    @unittest.skipIf(__name__ == "__main__",
                     "Test cannot run as when module is __main__")
    def test_main_class_reset_on_unpickle(self):
        mp = get_context('spawn')
        proc = mp.Process(target=check_main_class_reset_on_unpickle)
        proc.start()
        proc.join(timeout=60)
        self.assertEqual(proc.exitcode, 0)

    def test_dynamic_class_reset_on_unpickle_new_proc(self):
        # a dynamic class
        class Klass:
            classvar = None

        # serialize Klass in this process
        saved = dumps(Klass)

        # Check the reset problem in a new process
        mp = get_context('spawn')
        proc = mp.Process(target=check_unpickle_dyn_class_new_proc, args=(saved,))
        proc.start()
        proc.join(timeout=60)
        self.assertEqual(proc.exitcode, 0)

    def test_dynamic_class_issue_7356(self):
        cfunc = numba.njit(issue_7356)
        self.assertEqual(cfunc(), (100, 100))


class DynClass(object):
    # For testing issue #7356
    a = None


def issue_7356():
    with numba.objmode(before="intp"):
        DynClass.a = 100
        before = DynClass.a
    with numba.objmode(after="intp"):
        after = DynClass.a
    return before, after


def check_main_class_reset_on_unpickle():
    # Load module and get its global dictionary
    glbs = runpy.run_module(
        "numba.tests.cloudpickle_main_class",
        run_name="__main__",
    )
    # Get the Klass and check it is from __main__
    Klass = glbs['Klass']
    assert Klass.__module__ == "__main__"
    assert Klass.classvar != 100
    saved = dumps(Klass)
    # mutate
    Klass.classvar = 100
    # check
    _check_dyn_class(Klass, saved)


def check_unpickle_dyn_class_new_proc(saved):
    Klass = loads(saved)
    assert Klass.classvar != 100
    # mutate
    Klass.classvar = 100
    # check
    _check_dyn_class(Klass, saved)


def _check_dyn_class(Klass, saved):
    def check():
        if Klass.classvar != 100:
            raise AssertionError("Check failed. Klass reset.")

    check()
    loaded = loads(saved)
    if loaded is not Klass:
        raise AssertionError("Expected reuse")
    # Without the patch, each `loads(saved)` will reset `Klass.classvar`
    check()
    loaded = loads(saved)
    if loaded is not Klass:
        raise AssertionError("Expected reuse")
    check()


if __name__ == '__main__':
    unittest.main()
