"""
Implement code coverage support.

Currently contains logic to extend ``coverage`` with lines covered by the
compiler.
"""
from typing import Optional, Sequence, Callable, no_type_check
from collections import defaultdict
from abc import ABC, abstractmethod
import atexit
from functools import cache

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


@cache
def _get_coverage_data():
    """
    Make a singleton ``CoverageData``.
    Avoid writing to disk. Other processes can corrupt the file.
    """
    covdata = coverage.CoverageData(no_disk=True)
    cov = get_active_coverage()
    assert cov is not None, "no active Coverage instance"

    @atexit.register
    def _finalize():
        # Update active coverage
        cov.get_data().update(covdata)

    return covdata


class NotifyLocBase(ABC):
    """Interface for notifying visiting of a ``numba.core.ir.Loc``.
    """
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
    def __init__(self):
        self._arcs_data = defaultdict(set)

    def notify(self, loc: ir.Loc):
        if loc.filename.endswith(".py"):
            # The compiler doesn't actually know about arc.
            self._arcs_data[loc.filename].add((loc.line, loc.line))

    def close(self):
        covdata = _get_coverage_data()
        with covdata._lock:
            covdata.set_context("numba_compiled")
            covdata.add_arcs(self._arcs_data)


@_the_registry.append
def _register_coverage_notifier():
    if get_active_coverage() is not None:
        return NotifyCompilerCoverage()
