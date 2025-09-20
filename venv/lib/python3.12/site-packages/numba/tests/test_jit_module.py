import os
import sys
import inspect
import contextlib
import numpy as np
import logging
from io import StringIO

import unittest
from numba.tests.support import SerialMixin, create_temp_module
from numba.core import dispatcher


@contextlib.contextmanager
def captured_logs(l):
    try:
        buffer = StringIO()
        handler = logging.StreamHandler(buffer)
        l.addHandler(handler)
        yield buffer
    finally:
        l.removeHandler(handler)


class TestJitModule(SerialMixin, unittest.TestCase):

    source_lines = """
from numba import jit_module

def inc(x):
    return x + 1

def add(x, y):
    return x + y

def inc_add(x):
    y = inc(x)
    return add(x, y)

import numpy as np
mean = np.mean

class Foo(object):
    pass

jit_module({jit_options})
"""

    def test_create_temp_jitted_module(self):
        sys_path_original = list(sys.path)
        sys_modules_original = dict(sys.modules)
        with create_temp_module(self.source_lines) as test_module:
            temp_module_dir = os.path.dirname(test_module.__file__)
            self.assertEqual(temp_module_dir, sys.path[0])
            self.assertEqual(sys.path[1:], sys_path_original)
            self.assertTrue(test_module.__name__ in sys.modules)
        # Test that modifications to sys.path / sys.modules are reverted
        self.assertEqual(sys.path, sys_path_original)
        self.assertEqual(sys.modules, sys_modules_original)

    def test_create_temp_jitted_module_with_exception(self):
        try:
            sys_path_original = list(sys.path)
            sys_modules_original = dict(sys.modules)
            with create_temp_module(self.source_lines):
                raise ValueError("Something went wrong!")
        except ValueError:
            # Test that modifications to sys.path / sys.modules are reverted
            self.assertEqual(sys.path, sys_path_original)
            self.assertEqual(sys.modules, sys_modules_original)

    def test_jit_module(self):
        with create_temp_module(self.source_lines) as test_module:
            self.assertIsInstance(test_module.inc, dispatcher.Dispatcher)
            self.assertIsInstance(test_module.add, dispatcher.Dispatcher)
            self.assertIsInstance(test_module.inc_add, dispatcher.Dispatcher)
            self.assertTrue(test_module.mean is np.mean)
            self.assertTrue(inspect.isclass(test_module.Foo))

            # Test output of jitted functions is as expected
            x, y = 1.7, 2.3
            self.assertEqual(test_module.inc(x),
                             test_module.inc.py_func(x))
            self.assertEqual(test_module.add(x, y),
                             test_module.add.py_func(x, y))
            self.assertEqual(test_module.inc_add(x),
                             test_module.inc_add.py_func(x))

    def test_jit_module_jit_options(self):
        jit_options = {"nopython": True,
                       "nogil": False,
                       "error_model": "numpy",
                       "boundscheck": False,
                       }
        with create_temp_module(self.source_lines,
                                **jit_options) as test_module:
            self.assertEqual(test_module.inc.targetoptions, jit_options)

    def test_jit_module_jit_options_override(self):
        source_lines = """
from numba import jit, jit_module

@jit(nogil=True, forceobj=True)
def inc(x):
    return x + 1

def add(x, y):
    return x + y

jit_module({jit_options})
"""
        jit_options = {"nopython": True,
                       "error_model": "numpy",
                       "boundscheck": False,
                       }
        with create_temp_module(source_lines=source_lines,
                                **jit_options) as test_module:
            self.assertEqual(test_module.add.targetoptions, jit_options)
            # Test that manual jit-wrapping overrides jit_module options,
            # `forceobj` will automatically apply `nopython=False`.
            self.assertEqual(test_module.inc.targetoptions,
                             {'nogil': True, 'forceobj': True,
                              'boundscheck': None, 'nopython': False})

    def test_jit_module_logging_output(self):
        logger = logging.getLogger('numba.core.decorators')
        logger.setLevel(logging.DEBUG)
        jit_options = {"nopython": True,
                       "error_model": "numpy",
                       }
        with captured_logs(logger) as logs:
            with create_temp_module(self.source_lines,
                                    **jit_options) as test_module:
                logs = logs.getvalue()
                expected = ["Auto decorating function",
                            "from module {}".format(test_module.__name__),
                            "with jit and options: {}".format(jit_options)]
                self.assertTrue(all(i in logs for i in expected))

    def test_jit_module_logging_level(self):
        logger = logging.getLogger('numba.core.decorators')
        # Test there's no logging for INFO level
        logger.setLevel(logging.INFO)
        with captured_logs(logger) as logs:
            with create_temp_module(self.source_lines):
                self.assertEqual(logs.getvalue(), '')
