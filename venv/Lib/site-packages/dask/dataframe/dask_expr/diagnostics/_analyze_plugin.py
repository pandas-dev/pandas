from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, ClassVar

from distributed import Scheduler, SchedulerPlugin, Worker, WorkerPlugin
from distributed.protocol.pickle import dumps

if TYPE_CHECKING:
    from crick import TDigest


class Digest:
    count: int
    total: float
    sketch: TDigest

    def __init__(self) -> None:
        from crick import TDigest

        self.count = 0
        self.total = 0.0
        self.sketch = TDigest()

    def add(self, sample: float) -> None:
        self.count = self.count + 1
        self.total = self.total + sample
        self.sketch.add(sample)

    def merge(self, other: Digest) -> None:
        self.count = self.count + other.count
        self.total = self.total + other.total
        self.sketch.merge(other.sketch)

    @property
    def mean(self):
        return self.total / self.count


class AnalyzePlugin(SchedulerPlugin):
    idempotent: ClassVar[bool] = True
    name: ClassVar[str] = "analyze"
    _scheduler: Scheduler | None

    def __init__(self) -> None:
        self._scheduler = None

    async def start(self, scheduler: Scheduler) -> None:
        self._scheduler = scheduler
        scheduler.handlers["analyze_get_statistics"] = self.get_statistics
        worker_plugin = _AnalyzeWorkerPlugin()
        await self._scheduler.register_worker_plugin(
            None,
            dumps(worker_plugin),
            name=worker_plugin.name,
            idempotent=True,
        )

    async def get_statistics(self, id: str):
        assert self._scheduler is not None
        worker_statistics = await self._scheduler.broadcast(
            msg={"op": "analyze_get_statistics", "id": id}
        )
        cluster_statistics = Statistics()
        for statistics in worker_statistics.values():
            cluster_statistics.merge(statistics)
        return cluster_statistics


class ExpressionStatistics:
    _metric_digests: defaultdict[str, Digest]

    def __init__(self) -> None:
        self._metric_digests = defaultdict(Digest)

    def add(self, metric: str, value: float) -> None:
        self._metric_digests[metric].add(value)

    def merge(self, other: ExpressionStatistics) -> None:
        for metric, digest in other._metric_digests.items():
            self._metric_digests[metric].merge(digest)


class Statistics:
    _expr_statistics: defaultdict[str, ExpressionStatistics]

    def __init__(self) -> None:
        self._expr_statistics = defaultdict(ExpressionStatistics)

    def add(self, expr: str, metric: str, value: float):
        self._expr_statistics[expr].add(metric, value)

    def merge(self, other: Statistics):
        for expr, statistics in other._expr_statistics.items():
            self._expr_statistics[expr].merge(statistics)


class _AnalyzeWorkerPlugin(WorkerPlugin):
    idempotent: ClassVar[bool] = True
    name: ClassVar[str] = "analyze"
    _statistics: defaultdict[str, Statistics]
    _worker: Worker | None

    def __init__(self) -> None:
        self._worker = None
        self._statistics = defaultdict(Statistics)

    def setup(self, worker: Worker) -> None:
        self._digests = defaultdict(lambda: defaultdict(lambda: defaultdict(Digest)))  # type: ignore
        self._worker = worker
        self._worker.handlers["analyze_get_statistics"] = self.get_statistics

    def add(self, id: str, expr: str, metric: str, value: float):
        self._statistics[id].add(expr, metric, value)

    def get_statistics(self, id: str) -> Statistics:
        return self._statistics.pop(id)


def get_worker_plugin() -> _AnalyzeWorkerPlugin:
    from distributed import get_worker

    try:
        worker = get_worker()
    except ValueError as e:
        raise RuntimeError(
            "``.analyze()`` requires Dask's distributed scheduler"
        ) from e

    try:
        return worker.plugins["analyze"]
    except KeyError as e:
        raise RuntimeError(
            f"The worker {worker.address} does not have an Analyze plugin."
        ) from e
