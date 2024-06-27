import sys
from unittest import mock

import types
import warnings
import unittest
import os
import subprocess
import threading

from numba import config, njit
from numba.tests.support import TestCase
from numba.testing.main import _TIMEOUT as _RUNNER_TIMEOUT

if config.PYVERSION < (3, 9):
    import importlib_metadata
else:
    from importlib import metadata as importlib_metadata

_TEST_TIMEOUT = _RUNNER_TIMEOUT - 60.


class _DummyClass(object):
    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return '_DummyClass(%f, %f)' % self.value


class TestEntrypoints(TestCase):
    """
    Test registration of init() functions from Numba extensions
    """

    def test_init_entrypoint(self):
        # loosely based on Pandas test from:
        #   https://github.com/pandas-dev/pandas/pull/27488

        mod = mock.Mock(__name__='_test_numba_extension')

        try:
            # will remove this module at the end of the test
            sys.modules[mod.__name__] = mod

            my_entrypoint = importlib_metadata.EntryPoint(
                'init', '_test_numba_extension:init_func', 'numba_extensions',
            )

            with mock.patch.object(
                importlib_metadata,
                'entry_points',
                return_value={'numba_extensions': (my_entrypoint,)},
            ):

                from numba.core import entrypoints

                # Allow reinitialization
                entrypoints._already_initialized = False

                entrypoints.init_all()

                # was our init function called?
                mod.init_func.assert_called_once()

                # ensure we do not initialize twice
                entrypoints.init_all()
                mod.init_func.assert_called_once()
        finally:
            # remove fake module
            if mod.__name__ in sys.modules:
                del sys.modules[mod.__name__]

    def test_entrypoint_tolerance(self):
        # loosely based on Pandas test from:
        #   https://github.com/pandas-dev/pandas/pull/27488

        mod = mock.Mock(__name__='_test_numba_bad_extension')
        mod.configure_mock(**{'init_func.side_effect': ValueError('broken')})

        try:
            # will remove this module at the end of the test
            sys.modules[mod.__name__] = mod

            my_entrypoint = importlib_metadata.EntryPoint(
                'init',
                '_test_numba_bad_extension:init_func',
                'numba_extensions',
            )

            with mock.patch.object(
                importlib_metadata,
                'entry_points',
                return_value={'numba_extensions': (my_entrypoint,)},
            ):

                from numba.core import entrypoints
                # Allow reinitialization
                entrypoints._already_initialized = False

                with warnings.catch_warnings(record=True) as w:
                    entrypoints.init_all()

                bad_str = "Numba extension module '_test_numba_bad_extension'"
                for x in w:
                    if bad_str in str(x):
                        break
                else:
                    raise ValueError("Expected warning message not found")

                # was our init function called?
                mod.init_func.assert_called_once()

        finally:
            # remove fake module
            if mod.__name__ in sys.modules:
                del sys.modules[mod.__name__]

    _EP_MAGIC_TOKEN = 'RUN_ENTRY'

    @unittest.skipIf(os.environ.get('_EP_MAGIC_TOKEN', None) != _EP_MAGIC_TOKEN,
                     "needs token")
    def test_entrypoint_handles_type_extensions(self):
        # loosely based on Pandas test from:
        #   https://github.com/pandas-dev/pandas/pull/27488
        import numba

        def init_function():
            # This init function would normally just call a module init via
            # import or similar, for the sake of testing, inline registration
            # of how to handle the global "_DummyClass".
            class DummyType(numba.types.Type):
                def __init__(self):
                    super(DummyType, self).__init__(name='DummyType')

            @numba.extending.typeof_impl.register(_DummyClass)
            def typer_DummyClass(val, c):
                return DummyType()

            @numba.extending.register_model(DummyType)
            class DummyModel(numba.extending.models.StructModel):
                def __init__(self, dmm, fe_type):
                    members = [
                        ('value', numba.types.float64), ]
                    super(DummyModel, self).__init__(dmm, fe_type, members)

            @numba.extending.unbox(DummyType)
            def unbox_dummy(typ, obj, c):
                value_obj = c.pyapi.object_getattr_string(obj, "value")
                dummy_struct_proxy = numba.core.cgutils.create_struct_proxy(typ)
                dummy_struct = dummy_struct_proxy(c.context, c.builder)
                dummy_struct.value = c.pyapi.float_as_double(value_obj)
                c.pyapi.decref(value_obj)
                err_flag = c.pyapi.err_occurred()
                is_error = numba.core.cgutils.is_not_null(c.builder, err_flag)
                return numba.extending.NativeValue(dummy_struct._getvalue(),
                                                   is_error=is_error)

            @numba.extending.box(DummyType)
            def box_dummy(typ, val, c):
                dummy_struct_proxy = numba.core.cgutils.create_struct_proxy(typ)
                dummy_struct = dummy_struct_proxy(c.context, c.builder)
                value_obj = c.pyapi.float_from_double(dummy_struct.value)
                serialized_clazz = c.pyapi.serialize_object(_DummyClass)
                class_obj = c.pyapi.unserialize(serialized_clazz)
                res = c.pyapi.call_function_objargs(class_obj, (value_obj,))
                c.pyapi.decref(value_obj)
                c.pyapi.decref(class_obj)
                return res

        mod = types.ModuleType("_test_numba_init_sequence")
        mod.init_func = init_function

        try:
            # will remove this module at the end of the test
            sys.modules[mod.__name__] = mod

            my_entrypoint = importlib_metadata.EntryPoint(
                'init',
                '_test_numba_init_sequence:init_func',
                'numba_extensions',
            )

            with mock.patch.object(
                importlib_metadata,
                'entry_points',
                return_value={'numba_extensions': (my_entrypoint,)},
            ):
                @njit
                def foo(x):
                    return x

                ival = _DummyClass(10)
                foo(ival)
        finally:
            # remove fake module
            if mod.__name__ in sys.modules:
                del sys.modules[mod.__name__]

    def run_cmd(self, cmdline, env):
        popen = subprocess.Popen(cmdline,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 env=env)
        # finish in _TEST_TIMEOUT seconds or kill it
        timeout = threading.Timer(_TEST_TIMEOUT, popen.kill)
        try:
            timeout.start()
            out, err = popen.communicate()
            if popen.returncode != 0:
                raise AssertionError(
                    "process failed with code %s: stderr follows\n%s\n" %
                    (popen.returncode, err.decode()))
            return out.decode(), err.decode()
        finally:
            timeout.cancel()
        return None, None

    def test_entrypoint_extension_sequence(self):
        env_copy = os.environ.copy()
        env_copy['_EP_MAGIC_TOKEN'] = str(self._EP_MAGIC_TOKEN)
        themod = self.__module__
        thecls = type(self).__name__
        methname = 'test_entrypoint_handles_type_extensions'
        injected_method = '%s.%s.%s' % (themod, thecls, methname)
        cmdline = [sys.executable, "-m", "numba.runtests", injected_method]
        out, err = self.run_cmd(cmdline, env_copy)
        _DEBUG = False
        if _DEBUG:
            print(out, err)
