import unittest
from numba.core.compiler_lock import (
    global_compiler_lock,
    require_global_compiler_lock,
)
from numba.tests.support import TestCase


class TestCompilerLock(TestCase):
    def test_gcl_as_context_manager(self):
        with global_compiler_lock:
            require_global_compiler_lock()

    def test_gcl_as_decorator(self):
        @global_compiler_lock
        def func():
            require_global_compiler_lock()

        func()


if __name__ == '__main__':
    unittest.main()
