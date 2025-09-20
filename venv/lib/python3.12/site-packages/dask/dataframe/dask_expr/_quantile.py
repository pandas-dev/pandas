from __future__ import annotations

import functools

import numpy as np

from dask.dataframe.dask_expr._expr import DropnaSeries, Expr
from dask.dataframe.dispatch import make_meta, meta_nonempty
from dask.utils import import_required, is_series_like


def _finalize_scalar_result(cons, *args, **kwargs):
    return cons(*args, **kwargs)[0]


class SeriesQuantile(Expr):
    _parameters = ["frame", "q", "method"]
    _defaults = {"method": "default"}

    @functools.cached_property
    def q(self):
        q = np.array(self.operand("q"))
        if q.ndim > 0:
            assert len(q) > 0, f"must provide non-empty q={q}"
            q.sort(kind="mergesort")
            return q
        return np.asarray([self.operand("q")])

    @functools.cached_property
    def method(self):
        if self.operand("method") == "default":
            return "dask"
        else:
            return self.operand("method")

    @functools.cached_property
    def _meta(self):
        meta = self.frame._meta
        if not is_series_like(self.frame._meta):
            meta = meta.to_series()
        return make_meta(meta_nonempty(meta).quantile(self.operand("q")))

    def _divisions(self):
        if is_series_like(self._meta):
            return (np.min(self.q), np.max(self.q))
        return (None, None)

    @functools.cached_property
    def _constructor(self):
        meta = self.frame._meta
        if not is_series_like(self.frame._meta):
            meta = meta.to_series()
        return meta._constructor

    @functools.cached_property
    def _finalizer(self):
        if is_series_like(self._meta):
            return lambda tsk: (
                self._constructor,
                tsk,
                self.q,
                None,
                self.frame._meta.name,
            )
        else:
            return lambda tsk: (_finalize_scalar_result, self._constructor, tsk, [0])

    def _lower(self):
        frame = DropnaSeries(self.frame)
        if self.method == "tdigest":
            return SeriesQuantileTdigest(
                frame, self.operand("q"), self.operand("method")
            )
        else:
            return SeriesQuantileDask(frame, self.operand("q"), self.operand("method"))


class SeriesQuantileTdigest(SeriesQuantile):
    @functools.cached_property
    def _meta(self):
        import_required(
            "crick", "crick is a required dependency for using the tdigest method."
        )
        return super()._meta

    def _layer(self) -> dict:
        from dask.array.percentile import _percentiles_from_tdigest, _tdigest_chunk

        dsk = {}
        for i in range(self.frame.npartitions):
            dsk[("chunk-" + self._name, i)] = (
                _tdigest_chunk,
                (getattr, (self.frame._name, i), "values"),
            )

        dsk[(self._name, 0)] = self._finalizer(
            (_percentiles_from_tdigest, self.q * 100, sorted(dsk))
        )
        return dsk

    def _lower(self):
        return None


class SeriesQuantileDask(SeriesQuantile):
    def _layer(self) -> dict:
        from dask.array.dispatch import percentile_lookup as _percentile
        from dask.array.percentile import merge_percentiles

        dsk = {}
        # Add 0 and 100 during calculation for more robust behavior (hopefully)
        calc_qs = np.pad(self.q * 100, 1, mode="constant")
        calc_qs[-1] = 100

        for i in range(self.frame.npartitions):
            dsk[("chunk-" + self._name, i)] = (
                _percentile,
                (self.frame._name, i),
                calc_qs,
            )
        dsk[(self._name, 0)] = self._finalizer(
            (
                merge_percentiles,
                self.q * 100,
                [calc_qs] * self.frame.npartitions,
                sorted(dsk),
                "lower",
                None,
                False,
            )
        )
        return dsk

    def _lower(self):
        return None
