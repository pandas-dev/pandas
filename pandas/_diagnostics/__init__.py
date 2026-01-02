from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable
from dataclasses import dataclass
from typing import (
    Any,
    Generic,
    TypeVar,
)

from pandas._diagnostics.collector import DiagnosticsCollector
from pandas._diagnostics.introspect import infer_op_metadata
from pandas._diagnostics.report import DiagnosticsReport

T = TypeVar("T")


@dataclass(frozen=True)
class DiagnosticsResult(Generic[T]):
    """
    Wrapper so we don't mutate user/pandas objects to attach reports.
    """

    value: T
    report: DiagnosticsReport


def report(
    func: Callable[..., T],
    *args: Any,
    memory: str = "off",  # "off" | "peak" | "delta"
    phase_memory: bool = False,  # requires memory="delta"
    allocations: bool = False,  # requires memory="delta"
    profile: bool = False,  # cProfile hotspots + auto buckets
    bucket_phases: bool = True,  # derive auto phases from profile + memory buckets
    profile_limit: int = 50,
    allocations_limit: int = 50,
    include_io_names: bool = False,  # record CSV basename (not full path)
    include_diagnostics_phases: bool = False,
    **kwargs: Any,
) -> DiagnosticsResult[T]:
    """
    Run a callable with diagnostics around it (no monkeypatching, no internal plumbing).

    Wrapper phases:
      - diag.setup: metadata inference (diagnostics overhead)
      - diag.call: the callable itself
      - diag.finalize: profiler + tracemalloc post-processing (diagnostics overhead)

    By default, diag.* wrapper phases are filtered out of the returned report
    (include_diagnostics_phases=False). Set it True to see wrapper overhead phases.
    """
    if memory not in ("off", "peak", "delta"):
        raise ValueError("memory must be one of: 'off', 'peak', 'delta'")
    if phase_memory and memory != "delta":
        raise ValueError("phase_memory requires memory='delta'")
    if allocations and memory != "delta":
        raise ValueError("allocations requires memory='delta'")

    # Auto bucket phases only makes sense when profiling is enabled
    if not profile:
        bucket_phases = False

    collector = DiagnosticsCollector(
        memory_mode=memory,
        track_phase_memory=phase_memory,
        track_allocations=allocations,
        allocations_limit=allocations_limit,
        track_profile=profile,
        profile_limit=profile_limit,
        bucket_phases=bucket_phases,
    )

    with collector.activate():
        with collector.phase("diag.setup"):
            try:
                meta = infer_op_metadata(
                    func,
                    args,
                    kwargs,
                    include_io_names=include_io_names,
                )
                collector.set_meta(**meta)
            except Exception:
                pass

        with collector.phase("diag.call"):
            value = collector.run(func, *args, **kwargs)

        with collector.phase("diag.finalize"):
            collector.finalize()

        collector.close()

    rep = collector.to_report()

    # filter out wrapper phases unless user explicitly wants them
    if not include_diagnostics_phases:
        rep = rep.filtered(include_diagnostics_phases=False)

    return DiagnosticsResult(value=value, report=rep)


__all__ = ["DiagnosticsReport", "DiagnosticsResult", "report"]
