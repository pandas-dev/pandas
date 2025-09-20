"""
Implement code coverage support.

Currently contains logic to extend ``coverage`` with lines covered by the
compiler.
"""

from typing import Optional, Sequence, Callable, no_type_check
from collections.abc import Mapping
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

from numba.core import ir, config


try:
    import coverage
except ImportError:
    coverage_available = False
else:
    coverage_available = True


@no_type_check
def get_active_coverage():
    """Get active coverage instance or return None if not found.
    """
    cov = None
    if coverage_available:
        cov = coverage.Coverage.current()
    return cov


_the_registry: Callable[[], Optional["NotifyLocBase"]] = []


def get_registered_loc_notify() -> Sequence["NotifyLocBase"]:
    """
    Returns a list of the registered NotifyLocBase instances.
    """
    if not config.JIT_COVERAGE:
        # Coverage disabled.
        return []
    return list(filter(lambda x: x is not None,
                       (factory() for factory in _the_registry)))


class NotifyLocBase(ABC):
    """Interface for notifying visiting of a ``numba.core.ir.Loc``."""

    @abstractmethod
    def notify(self, loc: ir.Loc) -> None:
        pass

    @abstractmethod
    def close(self) -> None:
        pass


class NotifyCompilerCoverage(NotifyLocBase):
    """
    Use to notify ``coverage`` about compiled lines.

    The compiled lines are under the "numba_compiled" context in the coverage
    data.
    """

    def __init__(self, collector):
        self._collector = collector

        # see https://github.com/nedbat/coveragepy/blob/e7c05fe91ee36c0c94e144bb88d25db4fc3d02fd/coverage/collector.py#L261 # noqa E501
        tracer_kwargs = collector.core.tracer_kwargs.copy()
        tracer_kwargs.update(
            dict(
                data=collector.data,
                lock_data=collector.lock_data,
                unlock_data=collector.unlock_data,
                trace_arcs=collector.branch,
                should_trace=collector.should_trace,
                should_trace_cache=collector.should_trace_cache,
                warn=collector.warn,
                should_start_context=collector.should_start_context,
                switch_context=collector.switch_context,
                packed_arcs=collector.core.packed_arcs,
            )
        )
        self._tracer = NumbaTracer(**tracer_kwargs)
        collector.tracers.append(self._tracer)

    def notify(self, loc: ir.Loc):
        tracer = self._tracer
        if loc.filename.endswith(".py"):
            tracer.switch_context("numba_compiled")
            tracer.trace(loc)
            tracer.switch_context(None)

    def close(self):
        pass


@_the_registry.append
def _register_coverage_notifier():
    cov = get_active_coverage()
    if cov is not None:
        col = cov._collector
        # Is coverage started?
        if col.tracers:
            return NotifyCompilerCoverage(col)


if coverage_available:

    @dataclass(kw_only=True)
    class NumbaTracer(coverage.types.Tracer):
        """
        Not actually a tracer as in the coverage implementation, which will
        setup a Python trace function. This implementation pretends to trace
        but instead receives fake trace events for each line the compiler has
        visited.

        See coverage.PyTracer
        """

        data: coverage.types.TTraceData
        trace_arcs: bool
        should_trace: coverage.types.TShouldTraceFn
        should_trace_cache: Mapping[
            str, coverage.types.TFileDisposition | None
        ]
        should_start_context: coverage.types.TShouldStartContextFn | None
        switch_context: Callable[[str | None], None] | None
        lock_data: Callable[[], None]
        unlock_data: Callable[[], None]
        warn: coverage.types.TWarnFn
        packed_arcs: bool

        _activity: bool = field(default=False)

        def start(self) -> coverage.types.TTraceFn | None:
            """Start this tracer, return a trace function if based on
            sys.settrace."""
            return None

        def stop(self) -> None:
            """Stop this tracer."""
            return None

        def activity(self) -> bool:
            """Has there been any activity?"""
            return self._activity

        def reset_activity(self) -> None:
            """Reset the activity() flag."""
            self._activity = False

        def get_stats(self) -> dict[str, int] | None:
            """Return a dictionary of statistics, or None."""
            return None

        def trace(self, loc: ir.Loc) -> None:
            """Insert coverage data given source location.
            """
            # Check whether the file should be traced
            disp = self.should_trace_cache.get(loc.filename)
            if disp is None:
                disp = self.should_trace(loc.filename, None)
                self.should_trace_cache[loc.filename] = disp
            if not disp.trace:
                # Bail if not tracing the file
                return
            # Insert trace data
            tracename = disp.source_filename
            self.lock_data()
            cur_file_data = self.data.setdefault(tracename, set())
            if self.trace_arcs:
                if self.packed_arcs:
                    cur_file_data.add(_pack_arcs(loc.line, loc.line))
                else:
                    cur_file_data.add((loc.line, loc.line))
            else:
                cur_file_data.add(loc.line)
            self.unlock_data()
            # Mark activity for this tracer
            self._activity = True

    def _pack_arcs(l1: int, l2: int) -> int:
        """Pack arcs into a single integer for compatibility with .packed_arcs
        option.

        See
        https://github.com/nedbat/coveragepy/blob/e7c05fe91ee36c0c94e144bb88d25db4fc3d02fd/coverage/ctracer/tracer.c#L171
        """
        packed = 0
        if l1 < 0:
            packed |= 1 << 40
            l1 = -l1
        if l2 < 0:
            packed |= 1 << 41
            l2 = -l2
        packed |= (l2 << 20) + l1
        return packed
