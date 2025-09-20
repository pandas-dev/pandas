from io import StringIO
import logging

import unittest
from numba.core import tracing

logger = logging.getLogger('trace')
logger.setLevel(logging.INFO)

# Make sure tracing is enabled
orig_trace = tracing.trace
tracing.trace = tracing.dotrace

class CapturedTrace:
    """Capture the trace temporarily for validation."""

    def __init__(self):
        self.buffer = StringIO()
        self.handler = logging.StreamHandler(self.buffer)
    def __enter__(self):
        self._handlers = logger.handlers
        self.buffer = StringIO()
        logger.handlers = [logging.StreamHandler(self.buffer)]
    def __exit__(self, type, value, traceback):
        logger.handlers = self._handlers
    def getvalue(self):

        # Depending on how the tests are run, object names may be
        # qualified by their containing module.
        # Remove that to make the trace output independent from the testing mode.
        log = self.buffer.getvalue()
        log = log.replace(__name__ + '.','')
        return log

class Class(object):

    @tracing.trace
    @classmethod
    def class_method(cls):
        pass

    @tracing.trace
    @staticmethod
    def static_method():
        pass

    __test = None

    def _test_get(self):
        return self.__test

    def _test_set(self, value):
        self.__test = value

    test = tracing.trace(property(_test_get, _test_set))
        
    @tracing.trace
    def method(self, some, other='value', *args, **kwds):
        pass

    def __repr__(self):
        """Generate a deterministic string for testing."""
        return '<Class instance>'

class Class2(object):
    @classmethod
    def class_method(cls):
        pass

    @staticmethod
    def static_method():
        pass

    __test = None
    @property
    def test(self):
        return self.__test
    @test.setter
    def test(self, value):
        self.__test = value

    def method(self):
        pass

    def __str__(self):
        return 'Test(' + str(self.test) + ')'

    def __repr__(self):
        """Generate a deterministic string for testing."""
        return '<Class2 instance>'


@tracing.trace
def test(x, y, z = True):
    a = x + y
    b = x * y
    if z: return a
    else: return b

class TestTracing(unittest.TestCase):

    def __init__(self, *args):
        super(TestTracing, self).__init__(*args)

    def setUp(self):
        self.capture = CapturedTrace()

    def tearDown(self):
        del self.capture
        
    def test_method(self):

        with self.capture:
            Class().method('foo', bar='baz')
        self.assertEqual(self.capture.getvalue(),
                         ">> Class.method(self=<Class instance>, some='foo', other='value', bar='baz')\n" +
                         "<< Class.method\n")

    def test_class_method(self):

        with self.capture:
            Class.class_method()
        self.assertEqual(self.capture.getvalue(),
                         ">> Class.class_method(cls=<class 'Class'>)\n" +
                         "<< Class.class_method\n")

    def test_static_method(self):

        with self.capture:
            Class.static_method()
        self.assertEqual(self.capture.getvalue(),
                         ">> static_method()\n" +
                         "<< static_method\n")


    def test_property(self):

        with self.capture:
            test = Class()
            test.test = 1
            assert 1 == test.test
        self.assertEqual(self.capture.getvalue(),
                         ">> Class._test_set(self=<Class instance>, value=1)\n" +
                         "<< Class._test_set\n" +
                         ">> Class._test_get(self=<Class instance>)\n" +
                         "<< Class._test_get -> 1\n")

    def test_function(self):

        with self.capture:
            test(5, 5)
            test(5, 5, False)
        self.assertEqual(self.capture.getvalue(),
                         ">> test(x=5, y=5, z=True)\n" +
                         "<< test -> 10\n" +
                         ">> test(x=5, y=5, z=False)\n" +
                         "<< test -> 25\n")

    @unittest.skip("recursive decoration not yet implemented")
    def test_injected(self):

        with self.capture:
            tracing.trace(Class2, recursive=True)
            Class2.class_method()
            Class2.static_method()
            test = Class2()
            test.test = 1
            assert 1 == test.test
            test.method()

            self.assertEqual(self.capture.getvalue(),
                         ">> Class2.class_method(cls=<type 'Class2'>)\n" +
                         "<< Class2.class_method\n"
                         ">> static_method()\n"
                         "<< static_method\n")

            
# Reset tracing to its original value
tracing.trace = orig_trace

if __name__ == '__main__':
    unittest.main()
