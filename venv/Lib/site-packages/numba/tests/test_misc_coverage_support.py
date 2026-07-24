import unittest
from unittest.mock import MagicMock, patch

from numba.tests.support import TestCase

from numba import njit
from numba.core import ir
from numba.misc.coverage_support import (
    NotifyLocBase,
    _the_registry,
    coverage_available,
)


class TestMiscCoverageSupport(TestCase):
    @TestCase.run_test_in_subprocess(envvars={"NUMBA_JIT_COVERAGE": "1"})
    def test_custom_loc_notifier(self):
        class MyNotify(NotifyLocBase):
            records = []

            def notify(self, loc):
                self.records.append(("NOTIFY", loc))

            def close(self):
                self.records.append(("CLOSE", None))

        # Patch to install registry for testing
        new_the_registry = _the_registry + [MyNotify]
        gv = "numba.misc.coverage_support._the_registry"
        with patch(gv, new_the_registry):

            @njit
            def foo():
                return 123

            res = foo()

        self.assertEqual(res, 123)

        # offset by +2 because:
        # +1 for the decorator
        # +1 for the `def` line
        first_offset = 2
        offset = foo.__code__.co_firstlineno + first_offset
        loc = ir.Loc(__file__, 1)
        self.assertIn(("NOTIFY", loc.with_lineno(offset)), MyNotify.records)
        self.assertIn(("CLOSE", None), MyNotify.records)

        # Test dead branch pruned
        with patch(gv, new_the_registry):
            cond = False

            @njit
            def foo():
                if cond:
                    return 321
                return 123

            res = foo()

        self.assertEqual(res, 123)

        # `if cond` line is compiled
        offset = foo.__code__.co_firstlineno + first_offset
        self.assertIn(("NOTIFY", loc.with_lineno(offset)), MyNotify.records)

        # `    return 321` line is not compiled
        self.assertNotIn(
            ("NOTIFY", loc.with_lineno(offset + 1)), MyNotify.records
        )

        # `    return 123` line is compiled
        self.assertIn(("NOTIFY", loc.with_lineno(offset + 2)), MyNotify.records)

        self.assertIn(("CLOSE", None), MyNotify.records)


class TestNumbaTracerToolId(TestCase):
    @unittest.skipUnless(coverage_available, "coverage not installed")
    def test_tool_id_in_tracer_kwargs_accepted(self):
        # Regression test for https://github.com/numba/numba/issues/10589.
        # coverage >= 7.13.1 injects tool_id into tracer_kwargs when using the
        # sysmon core; NumbaTracer must accept it without raising TypeError.
        from numba.misc.coverage_support import NotifyCompilerCoverage

        collector = MagicMock()
        collector.core.tracer_kwargs = {'tool_id': 1, 'packed_arcs': False}
        collector.branch = False
        collector.should_start_context = None
        collector.switch_context = None
        collector.tracers = []

        NotifyCompilerCoverage(collector)


if __name__ == "__main__":
    unittest.main()
