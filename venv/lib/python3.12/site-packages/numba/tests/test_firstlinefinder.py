import unittest
import linecache
import inspect
import textwrap

from numba import njit
from numba.tests.support import TestCase
from numba.misc.firstlinefinder import get_func_body_first_lineno


class TestFirstLineFinder(TestCase):
    """
    The following methods contains tests that are sensitive to the source
    locations w.r.t. the beginning of each method.
    """

    def _get_grandparent_caller_code(self):
        frame = inspect.currentframe()
        caller_frame = inspect.getouterframes(frame)
        return caller_frame[2].frame.f_code

    def assert_line_location(self, expected, offset_from_caller):
        grandparent_co = self._get_grandparent_caller_code()
        lno = grandparent_co.co_firstlineno
        self.assertEqual(expected, lno + offset_from_caller)

    def test_decorated_odd_comment_indent(self):
        @njit
        def foo():
# NOTE: THIS COMMENT MUST START AT COLUMN 0 FOR THIS SAMPLE CODE TO BE VALID # noqa: E115, E501
            return 1

        first_def_line = get_func_body_first_lineno(foo)
        self.assert_line_location(first_def_line, 4)

    def test_undecorated_odd_comment_indent(self):
        def foo():
# NOTE: THIS COMMENT MUST START AT COLUMN 0 FOR THIS SAMPLE CODE TO BE VALID # noqa: E115, E501
            return 1

        first_def_line = get_func_body_first_lineno(njit(foo))
        self.assert_line_location(first_def_line, 3)

    def test_unnamed_lambda(self):
        foo = lambda: 1
        first_def_line = get_func_body_first_lineno(njit(foo))
        # Cannot determine first line of lambda function
        self.assertIsNone(first_def_line)

    def test_nested_function(self):
        def foo():
            @njit
            def foo():
                # This function is intentionally named the same as the parent
                return 1

            return foo

        inner = foo()
        first_def_line = get_func_body_first_lineno(inner)
        self.assert_line_location(first_def_line, 5)

    def test_pass_statement(self):
        @njit
        def foo():
            pass

        first_def_line = get_func_body_first_lineno(foo)
        self.assert_line_location(first_def_line, 3)

    def test_string_eval(self):
        source = """def foo():
            pass
        """

        globalns = {}
        exec(source, globalns)
        foo = globalns['foo']

        first_def_line = get_func_body_first_lineno(foo)
        # Cannot determine first line of string evaled functions
        self.assertIsNone(first_def_line)

    def _test_with_patched_linecache(self, filename, source,
                                     function_name, expected_first_line):
        # Modify the line cache in a similar manner to Jupyter, so that
        # get_func_body_first_lineno can find the code using the fallback to
        # inspect.getsourcelines()
        timestamp = None
        entry = (len(source), timestamp, source.splitlines(True), filename)
        linecache.cache[filename] = entry

        # We need to compile the code so we can give it the fake filename used
        # in the linecache
        code = compile(source, filename, "exec")

        globalns = {}
        exec(code, globalns)
        function = globalns[function_name]

        # We should be able to determine the first line number even though the
        # source does not exist on disk
        first_def_line = get_func_body_first_lineno(function)
        self.assertEqual(first_def_line, expected_first_line)

    def test_string_eval_linecache_basic(self):
        source = """def foo():
            pass
        """

        filename = "<foo-basic>"
        function_name = "foo"
        expected_first_line = 2
        self._test_with_patched_linecache(filename, source, function_name,
                                          expected_first_line)

    def test_string_eval_linecache_indent(self):
        source = """if True:
        # indent designed to test against potential indent error in ast.parse

        def foo():
            pass
        """

        filename = "<foo-indent>"
        function_name = "foo"
        expected_first_line = 5
        self._test_with_patched_linecache(filename, source, function_name,
                                          expected_first_line)

    def test_string_eval_linecache_closure(self):
        source = textwrap.dedent("""
        def foo_gen():
            def foo():
                pass
            return foo

        generated_foo = foo_gen()
        """)

        filename = "<foo-gen>"
        function_name = "generated_foo"
        expected_first_line = 4
        self._test_with_patched_linecache(filename, source, function_name,
                                          expected_first_line)

    def test_string_eval_linecache_stacked_decorators(self):
        source = textwrap.dedent("""
        def decorator(function):
            return function

        @decorator
        @decorator
        @decorator
        def decorated():
            pass
        """)

        filename = "<foo-stacked-decorator>"
        function_name = "decorated"
        expected_first_line = 9
        self._test_with_patched_linecache(filename, source, function_name,
                                          expected_first_line)

    def test_string_eval_linecache_all(self):
        # A test combining indented code, a closure, and stacked decorators
        source = """if 1:
        def decorator(function):
            return function

        def gen_decorated_foo():
            @decorator
            @decorator
            @decorator
            def _foo():
                pass

            return _foo

        foo_all = gen_decorated_foo()
        """

        filename = "<foo-all>"
        function_name = "foo_all"
        expected_first_line = 10
        self._test_with_patched_linecache(filename, source, function_name,
                                          expected_first_line)

    def test_single_line_function(self):
        @njit
        def foo(): pass   # noqa: E704

        first_def_line = get_func_body_first_lineno(foo)
        self.assert_line_location(first_def_line, 2)

    def test_docstring(self):
        @njit
        def foo():
            """Docstring
            """
            pass

        first_def_line = get_func_body_first_lineno(foo)
        self.assert_line_location(first_def_line, 5)

    def test_docstring_2(self):
        @njit
        def foo():
            """Docstring
            """
            """Not Docstring, but a bare string literal
            """
            pass
        # Variation on test_docstring but with a "fake" docstring following
        # the true docstring.
        first_def_line = get_func_body_first_lineno(foo)
        self.assert_line_location(first_def_line, 5)


if __name__ == "__main__":
    unittest.main()
