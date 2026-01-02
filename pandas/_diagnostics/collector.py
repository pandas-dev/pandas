from __future__ import annotations

from collections import defaultdict
import contextlib
import contextvars
from dataclasses import dataclass
import os
import time
from typing import (
    TYPE_CHECKING,
    Any,
    TypeVar,
)

from pandas._diagnostics.phase_rules import bucket_for
from pandas._diagnostics.report import (
    DiagnosticsReport,
    PhaseRecord,
)

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Iterator,
    )

T = TypeVar("T")

_ACTIVE_COLLECTOR: contextvars.ContextVar[DiagnosticsCollector | None] = (
    contextvars.ContextVar("pandas__diagnostics_active_collector", default=None)
)


def collector() -> DiagnosticsCollector | None:
    return _ACTIVE_COLLECTOR.get()


@contextlib.contextmanager
def phase(c: DiagnosticsCollector | None, name: str, **meta: Any) -> Iterator[None]:
    if c is None:
        yield
        return
    with c.phase(name, **meta):
        yield


@dataclass(frozen=True)
class _PhaseDone:
    name: str
    seconds: float
    mem_bytes_net: int | None = None
    mem_current_bytes: int | None = None
    mem_peak_bytes: int | None = None


class DiagnosticsCollector:
    """
    Pay-for-play collector.

    memory_mode:
      - "off": no tracemalloc
      - "peak": current/peak traced memory only (no snapshots)
      - "delta": snapshots to compute net bytes and optional per-phase deltas
    """

    def __init__(
        self,
        *,
        memory_mode: str = "off",
        track_phase_memory: bool = False,  # requires memory_mode="delta"
        track_allocations: bool = False,  # requires memory_mode="delta"
        allocations_limit: int = 50,
        track_profile: bool = False,
        profile_limit: int = 50,
        bucket_phases: bool = True,
    ) -> None:
        self._memory_mode = memory_mode
        self._track_phase_memory = track_phase_memory
        self._track_allocations = track_allocations
        self._allocations_limit = int(allocations_limit)

        self._track_profile = track_profile
        self._profile_limit = int(profile_limit)
        self._bucket_phases = bool(bucket_phases)

        self._t0: float | None = None
        self._total_s: float | None = None

        self._phase_stack: list[tuple[str, float]] = []
        self._phases_done: list[_PhaseDone] = []

        self._counters: dict[str, int] = {}
        self._meta: dict[str, Any] = {}

        self._mem_bytes_net: int | None = None
        self._mem_current: int | None = None
        self._mem_peak: int | None = None

        self._profile_top: list[dict[str, Any]] = []
        self._alloc_top: list[dict[str, Any]] = []
        self._auto_phases_done: list[_PhaseDone] = []

        self._finalized = False
        self._closed = False

        self._snap_before: Any | None = None
        self._stats_by_file: Any | None = None
        self._prof_stats: Any | None = None
        self._profile_obj: Any | None = None

        self._tracemalloc_was_on: bool = False
        self._tracemalloc_started_by_us: bool = False

    @contextlib.contextmanager
    def activate(self) -> Iterator[None]:
        token = _ACTIVE_COLLECTOR.set(self)
        try:
            yield
        finally:
            _ACTIVE_COLLECTOR.reset(token)

    def set_meta(self, **meta: Any) -> None:
        self._meta.update(meta)

    def inc(self, name: str, value: int = 1) -> None:
        self._counters[name] = self._counters.get(name, 0) + int(value)

    @contextlib.contextmanager
    def phase(self, name: str, **meta: Any) -> Iterator[None]:
        # meta currently unused; reserved for future
        snap_before = None
        if self._memory_mode == "delta" and self._track_phase_memory:
            import tracemalloc

            if tracemalloc.is_tracing():
                snap_before = tracemalloc.take_snapshot()

        start = time.perf_counter()
        self._phase_stack.append((name, start))
        try:
            yield
        finally:
            _n, s = self._phase_stack.pop()
            end = time.perf_counter()

            mem_delta = None
            mem_cur = mem_peak = None

            if snap_before is not None:
                try:
                    import tracemalloc

                    if tracemalloc.is_tracing():
                        snap_after = tracemalloc.take_snapshot()
                        stats = snap_after.compare_to(snap_before, "filename")
                        mem_delta = int(sum(st.size_diff for st in stats))
                except Exception:
                    mem_delta = None

            # Only touch tracemalloc if we are actually using it.
            if self._memory_mode != "off":
                try:
                    import tracemalloc

                    if tracemalloc.is_tracing():
                        mem_cur, mem_peak = tracemalloc.get_traced_memory()
                except Exception:
                    mem_cur = mem_peak = None

            self._phases_done.append(
                _PhaseDone(
                    name=name,
                    seconds=end - s,
                    mem_bytes_net=mem_delta,
                    mem_current_bytes=int(mem_cur) if mem_cur is not None else None,
                    mem_peak_bytes=int(mem_peak) if mem_peak is not None else None,
                )
            )

    def run(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        self._t0 = time.perf_counter()

        # Setup profiling/memory tracing here; heavy work deferred to finalize().
        if self._track_profile:
            import cProfile

            self._profile_obj = cProfile.Profile()
            self._profile_obj.enable()

        if self._memory_mode != "off":
            import tracemalloc

            self._tracemalloc_was_on = tracemalloc.is_tracing()
            if not self._tracemalloc_was_on:
                tracemalloc.start()
                self._tracemalloc_started_by_us = True
            if self._memory_mode == "delta":
                self._snap_before = tracemalloc.take_snapshot()

        try:
            return func(*args, **kwargs)
        finally:
            t0 = self._t0 if self._t0 is not None else time.perf_counter()
            self._total_s = time.perf_counter() - t0

    def finalize(self) -> None:
        """
        Collect profiler + tracemalloc results.
        """
        if self._finalized:
            return
        self._finalized = True

        # Stop profiler and capture stats
        if self._profile_obj is not None:
            import pstats

            self._profile_obj.disable()
            self._prof_stats = pstats.Stats(self._profile_obj)
            self._profile_top = self._build_profile_top(self._prof_stats)

        # Collect memory snapshots / peak
        if self._memory_mode != "off":
            import tracemalloc

            try:
                if tracemalloc.is_tracing():
                    cur, peak = tracemalloc.get_traced_memory()
                    self._mem_current = int(cur)
                    self._mem_peak = int(peak)
            except Exception:
                self._mem_current = self._mem_peak = None

            if self._memory_mode == "delta" and self._snap_before is not None:
                try:
                    if tracemalloc.is_tracing():
                        snap_after = tracemalloc.take_snapshot()
                        self._stats_by_file = snap_after.compare_to(
                            self._snap_before, "filename"
                        )
                        self._mem_bytes_net = int(
                            sum(st.size_diff for st in self._stats_by_file)
                        )

                        if self._track_allocations:
                            allocs = []
                            for st in self._stats_by_file:
                                try:
                                    filename = st.traceback[0].filename.replace(
                                        os.sep, "/"
                                    )
                                    allocs.append(
                                        (
                                            int(st.size_diff),
                                            int(st.count_diff),
                                            filename,
                                        )
                                    )
                                except Exception:
                                    continue
                            allocs.sort(key=lambda x: abs(x[0]), reverse=True)
                            self._alloc_top = [
                                {"bytes_net": b, "count_net": c, "filename": f}
                                for (b, c, f) in allocs[: self._allocations_limit]
                            ]
                except Exception:
                    self._stats_by_file = None

        # Auto phases: aggregate profiler exclusive time + memory net diffs into buckets
        if self._bucket_phases and self._prof_stats is not None:
            time_b = self._bucketize_profile_exclusive(self._prof_stats)
            mem_b: dict[str, int] = {}
            if self._stats_by_file is not None:
                mem_b = self._bucketize_memory_by_filename(self._stats_by_file)

            done: list[_PhaseDone] = []
            for name, secs in sorted(
                time_b.items(), key=lambda kv: kv[1], reverse=True
            ):
                done.append(
                    _PhaseDone(
                        name=f"auto:{name}",
                        seconds=float(secs),
                        mem_bytes_net=int(mem_b.get(name, 0)) if mem_b else None,
                    )
                )
            self._auto_phases_done = done

    def close(self) -> None:
        """
        Stop tracemalloc if we started it.
        This is separated from finalize() so phase snapshots remain valid.
        """
        if self._closed:
            return
        self._closed = True

        if self._tracemalloc_started_by_us:
            try:
                import tracemalloc

                if tracemalloc.is_tracing():
                    tracemalloc.stop()
            except Exception:
                pass

    def _build_profile_top(self, stats: Any) -> list[dict[str, Any]]:
        items: list[tuple[float, str, int | None, str, int]] = []
        for (filename, line, funcname), (
            _cc,
            nc,
            _tt,
            ct,
            _callers,
        ) in stats.stats.items():
            loc = f"{filename}:{funcname}".replace(os.sep, "/")
            if "pandas" not in loc.lower():
                continue
            items.append(
                (
                    float(ct),
                    filename.replace(os.sep, "/"),
                    int(line) if isinstance(line, int) else None,
                    funcname,
                    int(nc),
                )
            )
        items.sort(key=lambda x: x[0], reverse=True)
        out: list[dict[str, Any]] = []
        for ct, filename, line, funcname, nc in items[: self._profile_limit]:
            out.append(
                {
                    "cum_seconds": ct,
                    "calls": nc,
                    "filename": filename,
                    "line": line,
                    "function": funcname,
                }
            )
        return out

    def _bucketize_profile_exclusive(self, stats: Any) -> dict[str, float]:
        buckets: dict[str, float] = defaultdict(float)
        for (filename, _line, funcname), (
            _cc,
            _nc,
            tt,
            _ct,
            _callers,
        ) in stats.stats.items():
            loc = f"{filename}:{funcname}".replace(os.sep, "/")
            b = bucket_for(loc)
            if not b:
                continue
            buckets[b] += float(tt)  # exclusive time
        return dict(buckets)

    def _bucketize_memory_by_filename(self, stats_by_file: Any) -> dict[str, int]:
        buckets: dict[str, int] = defaultdict(int)
        for st in stats_by_file:
            try:
                filename = st.traceback[0].filename.replace(os.sep, "/")
                b = bucket_for(filename)
                if not b:
                    continue
                buckets[b] += int(st.size_diff)
            except Exception:
                continue
        return dict(buckets)

    def to_report(self) -> DiagnosticsReport:
        return DiagnosticsReport(
            total_seconds=float(self._total_s or 0.0),
            phases=[
                PhaseRecord(
                    p.name,
                    p.seconds,
                    mem_bytes_net=p.mem_bytes_net,
                    mem_current_bytes=p.mem_current_bytes,
                    mem_peak_bytes=p.mem_peak_bytes,
                )
                for p in self._phases_done
            ],
            auto_phases=[
                PhaseRecord(p.name, p.seconds, mem_bytes_net=p.mem_bytes_net)
                for p in self._auto_phases_done
            ],
            counters=dict(self._counters),
            metadata=dict(self._meta),
            memory_bytes_net=self._mem_bytes_net,
            memory_current_bytes=self._mem_current,
            memory_peak_bytes=self._mem_peak,
            profile_top=list(self._profile_top) if self._profile_top else [],
            allocations_top=list(self._alloc_top) if self._alloc_top else [],
        )


__all__ = ["DiagnosticsCollector", "collector", "phase"]
