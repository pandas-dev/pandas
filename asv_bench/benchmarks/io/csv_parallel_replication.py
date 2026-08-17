"""
TEMPORARY benchmarks for GH-66152 -- revert before merging.

These exist so reviewers can reproduce the numbers in the PR description on
their own hardware.  They are an asv port of the ``core``/``threads`` packs of
the standalone suite the PR was developed against, and they are deliberately
heavier than a normal asv benchmark (668 MB of fixtures): the parallel
``read_csv`` path only engages on real files above 5 MB, so there is no cheap
way to measure it.

GH-66152 changes only the *default* worker count -- pure Python, no C or
Cython.  ``mode.max_threads`` is honoured on both main and this branch, so
both behaviours are reachable from a single build:

    workers=4         what main does by default (min(os.cpu_count(), 4))
    workers=None      what this branch does (physical_core_count(), clamped)

That means no second build and no ``asv continuous`` is needed to reproduce
the PR's headline ratio -- one binary, one page cache, one process tree.

Reviewer quickstart -- from the repo root, against an already-built pandas.
Must be ``-m``: running this file by path would put ``benchmarks/io/`` on
``sys.path``, where its ``csv.py``/``json.py``/``pickle.py`` shadow the stdlib
modules of those names and break numpy's import.

    # standalone, no asv environment required; prints the A/B table directly
    python -m asv_bench.benchmarks.io.csv_parallel_replication

    # the full worker-count curve, to find the knee on your machine
    python -m asv_bench.benchmarks.io.csv_parallel_replication --sweep

    # or via asv, also no rebuild
    cd asv_bench && asv run --python=same -b "io\\.csv_parallel_replication"

The standalone runner prints the path of the pandas it imported -- check it is
the build you meant, since an editable install elsewhere on ``sys.path`` will
quietly win.

Fixtures are written once to ``$PANDAS_CSV_BENCH_DIR`` (default: a
``pandas-csv-parallel-fixtures`` directory under the system temp dir) and
reused on later runs.  Delete that directory to reclaim the space.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
import zlib

import numpy as np

import pandas as pd
from pandas import (
    DataFrame,
    read_csv,
)

# GH-66152 adds pandas.compat._cpu.  Reviewers running `asv continuous main
# HEAD` build main too, where it does not exist, so the detection helpers have
# to be optional or benchmark collection fails on the main side.
try:
    from pandas.compat._cpu import (
        available_cpu_count,
        physical_core_count,
    )
except ImportError:  # pandas without GH-66152
    available_cpu_count = None
    physical_core_count = None

try:
    from pandas.io.parsers.readers import _PARALLEL_READ_MIN_BYTES
except ImportError:  # pandas without GH-64347
    _PARALLEL_READ_MIN_BYTES = 5 * 1024 * 1024


FIXTURE_DIR = os.environ.get("PANDAS_CSV_BENCH_DIR") or os.path.join(
    tempfile.gettempdir(), "pandas-csv-parallel-fixtures"
)

WORDS = [
    "alpha", "bravo", "charlie", "delta", "echo", "foxtrot", "golf", "hotel",
    "india", "juliett", "kilo", "lima", "mike", "november", "oscar", "papa",
]  # fmt: skip

# Shapes match the standalone suite so the ms figures here line up with the
# PR table.  The comment on each is the speedup that suite measured on an
# M3 Pro (12 physical cores) for 4 workers -> default workers.
FIXTURES = {
    "float": "1M x 10 float64                        (108 MB, 1.77x)",
    "int": "1M x 10 int64                          (99 MB, 1.78x)",
    "bool": "1M x 10 bool                           (55 MB, 1.94x)",
    "str": "1M x 10 low-cardinality words          (62 MB, 1.69x)",
    "str_highcard": "1M x 10 near-distinct 12-char strings   (130 MB, 1.30x)",
    "mixed": "1M x 10 int/float/str/bool             (85 MB, 1.82x)",
    "wide": "100K x 100 int/float                   (103 MB, 1.23x)",
    # 70K rows, not 60K: at 60K this landed at 5,123,859 bytes, *under* the
    # 5 MiB gate, so every worker count read it serially and the sweep row was
    # flat by construction rather than by measurement.
    "medium": "70K x 10 mixed, just over the 5 MiB cut (6 MB)",
    # The one case that gets SLOWER (0.93x on 2026-07-27).  Carried here
    # deliberately: the full suite found it, and a reviewer-facing benchmark
    # that showed only the wins would be cherry-picking.
    "str_mixed_width": "300K x 3 strings of width 3/12/140      (47 MB, 0.93x)",
}

# The two shapes worth sweeping worker count over: the representative
# multi-dtype file, and the smallest file that still takes the parallel path
# (where per-worker startup is the largest share of the work, so it is where
# too many workers would show up as a regression).
SWEEP_FIXTURES = ["mixed", "medium"]

WORKER_SWEEP = [1, 2, 4, 6, 8, 12, 16]


def _rng(name: str) -> np.random.Generator:
    """Per-fixture seeded generator, so any subset reproduces identically.

    crc32 rather than hash(): str hashing is salted per process, so hash()
    would give a different file every time one was regenerated.
    """
    return np.random.default_rng(zlib.crc32(name.encode()))


def _strings(rng, nrows, pool=WORDS):
    return np.asarray(pool, dtype=object)[rng.integers(0, len(pool), nrows)]


def _highcard(rng, nrows, width=12):
    letters = np.array(list("abcdefghijklmnopqrstuvwxyz"))
    raw = letters[rng.integers(0, 26, (nrows, width))]
    return np.array(["".join(row) for row in raw], dtype=object)


def _floats(rng, nrows):
    # Rounded to 6 decimals: a full float64 repr would write ~17 significant
    # digits, which no real CSV looks like and which inflates the file 2x.
    return np.round(rng.random(nrows) * 1000, 6)


def _mixed(rng, nrows, ncols):
    out = {}
    for i in range(ncols):
        kind = i % 4
        if kind == 0:
            out[f"c{i}"] = rng.integers(0, 10**9, nrows)
        elif kind == 1:
            out[f"c{i}"] = _floats(rng, nrows)
        elif kind == 2:
            out[f"c{i}"] = _strings(rng, nrows)
        else:
            out[f"c{i}"] = rng.random(nrows) < 0.5
    return out


def _build(name, rng):
    nrows, ncols = 1_000_000, 10
    cols = [f"c{i}" for i in range(ncols)]
    if name == "float":
        return {col: _floats(rng, nrows) for col in cols}
    if name == "int":
        return {col: rng.integers(0, 10**9, nrows) for col in cols}
    if name == "bool":
        return {col: rng.random(nrows) < 0.5 for col in cols}
    if name == "str":
        return {col: _strings(rng, nrows) for col in cols}
    if name == "str_highcard":
        return {col: _highcard(rng, nrows) for col in cols}
    if name == "mixed":
        return _mixed(rng, nrows, ncols)
    if name == "wide":
        return {
            f"c{i}": (
                rng.integers(0, 10**9, 100_000) if i % 2 else _floats(rng, 100_000)
            )
            for i in range(100)
        }
    if name == "medium":
        return _mixed(rng, 70_000, ncols)
    if name == "str_mixed_width":
        # Width varied at fixed (high) cardinality, so the effect is token
        # length rather than dictionary-hit rate.
        return {
            "tiny": _highcard(rng, 300_000, width=3),
            "mid": _highcard(rng, 300_000, width=12),
            "long": _highcard(rng, 300_000, width=140),
        }
    raise ValueError(name)


def fixture_path(name: str) -> str:
    return os.path.join(FIXTURE_DIR, f"{name}.csv")


def ensure_fixtures(names=None) -> None:
    """Write any missing fixtures.  Idempotent, so it is cheap to re-call."""
    os.makedirs(FIXTURE_DIR, exist_ok=True)
    for name in names or FIXTURES:
        path = fixture_path(name)
        if os.path.exists(path):
            _check_eligible(path)
            continue
        print(f"generating {name}.csv ...", file=sys.stderr, flush=True)
        # Write-then-rename: a to_csv interrupted partway through a 130 MB
        # fixture would otherwise leave a truncated file that every later run
        # silently reuses -- and a truncation below the gate flips the read to
        # serial with no signal at all.
        tmp = f"{path}.tmp"
        DataFrame(_build(name, _rng(name))).to_csv(tmp, index=False)
        os.replace(tmp, path)
        _check_eligible(path)


def _check_eligible(path: str) -> None:
    """
    Fail loudly if a fixture is too small to take the parallel path.

    Every fixture here exists to measure parallel reading, so one that slips
    under the gate does not fail -- it silently reports serial timings that
    look like a flat, well-behaved result.
    """
    size = os.path.getsize(path)
    if size < _PARALLEL_READ_MIN_BYTES:
        raise AssertionError(
            f"{os.path.basename(path)} is {size:,} bytes, under the "
            f"{_PARALLEL_READ_MIN_BYTES:,}-byte parallel-read gate: it would be "
            "read serially at every worker count"
        )


def set_workers(workers: int | None) -> None:
    try:
        pd.set_option("mode.max_threads", workers)
    except pd.errors.OptionError:
        # A pandas predating the option has only the serial path, so every
        # setting collapses to it; swallowing this is what lets an older
        # commit serve as the baseline at all.
        pass


class ParallelWorkerDefault:
    """
    Reproduce the PR's headline ratio: 4 workers vs the new default.

    ``workers=1`` is included as the serial reference, so the table also shows
    parallel efficiency rather than just the ratio the PR claims.
    """

    params = (list(FIXTURES), [1, 4, None])
    param_names = ["fixture", "workers"]

    # Fixture generation happens in setup on a cold cache; a 100 MB to_csv is
    # well past asv's 60 s default.
    timeout = 2400
    # One read per sample, and no asv warmup: setup does a single untimed read
    # instead, which keeps page-cache warming out of the timed region without
    # burning a warmup read per sample on the cheap cases.
    number = 1
    warmup_time = 0
    repeat = 7

    def setup(self, fixture, workers):
        ensure_fixtures()
        self.path = fixture_path(fixture)
        set_workers(workers)
        read_csv(self.path)

    def teardown(self, fixture, workers):
        set_workers(None)

    def time_read_csv(self, fixture, workers):
        read_csv(self.path)


class ParallelWorkerSweep:
    """
    Worker-count curve, to locate the knee on the reviewer's own machine.

    This is the measurement that actually decides GH-66152: the PR's premise is
    that throughput keeps climbing past 4 workers up to the physical core
    count.  Compare where this curve flattens against
    ``track_physical_core_count`` below.
    """

    params = (SWEEP_FIXTURES, WORKER_SWEEP)
    param_names = ["fixture", "workers"]

    timeout = 2400
    number = 1
    warmup_time = 0
    repeat = 7

    def setup(self, fixture, workers):
        ensure_fixtures(SWEEP_FIXTURES)
        self.path = fixture_path(fixture)
        set_workers(workers)
        read_csv(self.path)

    def teardown(self, fixture, workers):
        set_workers(None)

    def time_read_csv(self, fixture, workers):
        read_csv(self.path)


# There is deliberately no `track_` benchmark for the detected core count.
# One reported the worker count each side dispatches (an honest 4 -> 12), but
# asv scores every track_ metric as lower-is-better, so a *correct* 4 -> 12
# rendered as a 3.00x "regression" and made `asv continuous` print
# "PERFORMANCE DECREASED" and exit 1 on a run where every timing improved.
# Worker counts come from the standalone runner's header instead:
#     python -m asv_bench.benchmarks.io.csv_parallel_replication --fixtures medium


# ---------------------------------------------------------------------------
# Standalone runner: same measurements without asv, for reviewers who would
# rather not stand up an asv environment.  asv ignores everything below.
# ---------------------------------------------------------------------------


def _best_of(path, reps):
    best = float("inf")
    for _ in range(reps):
        start = time.perf_counter()
        read_csv(path)
        best = min(best, time.perf_counter() - start)
    return best * 1000


def _measure(names, worker_counts, rounds, reps):
    """Interleave worker counts within each round so drift hits them equally."""
    results = {
        (name, workers): float("inf") for name in names for workers in worker_counts
    }
    for name in names:
        set_workers(1)
        read_csv(fixture_path(name))  # warm the page cache, untimed
    for round_ in range(rounds):
        print(f"  round {round_ + 1}/{rounds}", file=sys.stderr, flush=True)
        for name in names:
            for workers in worker_counts:
                set_workers(workers)
                key = (name, workers)
                results[key] = min(results[key], _best_of(fixture_path(name), reps))
    set_workers(None)
    return results


def _print_environment():
    print(f"pandas      {pd.__version__}  ({os.path.dirname(pd.__file__)})")
    print(f"os.cpu_count()        {os.cpu_count()}")
    if physical_core_count is None:
        print("physical_core_count() unavailable (build predates GH-66152)")
    else:
        print(f"physical_core_count() {physical_core_count()}")
        constrained = available_cpu_count()
        print(
            f"available_cpu_count() {constrained if constrained else 'unconstrained'}"
        )
    try:
        from pandas.io.parsers.readers import _default_n_workers

        print(f"default workers       {_default_n_workers()}")
    except ImportError:
        print(f"default workers       {min(os.cpu_count() or 1, 4)} (main's formula)")
    print(f"fixtures    {FIXTURE_DIR}")
    print()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--sweep",
        action="store_true",
        help="measure the full worker-count curve instead of the 4-vs-default A/B",
    )
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument(
        "--reps", type=int, default=3, help="reads per round; best wins"
    )
    parser.add_argument(
        "--fixtures",
        help="comma-separated subset, e.g. mixed,bool (default: all)",
    )
    args = parser.parse_args(argv)

    _print_environment()

    if args.fixtures:
        names = [name.strip() for name in args.fixtures.split(",")]
        unknown = set(names) - set(FIXTURES)
        if unknown:
            parser.error(f"unknown fixture(s): {', '.join(sorted(unknown))}")
    elif args.sweep:
        names = SWEEP_FIXTURES
    else:
        names = list(FIXTURES)

    ensure_fixtures(names)
    worker_counts = WORKER_SWEEP if args.sweep else [1, 4, None]
    results = _measure(names, worker_counts, args.rounds, args.reps)

    headers = [
        "default" if workers is None else str(workers) for workers in worker_counts
    ]
    width = max(len(name) for name in names) + 2
    print()
    print("read_csv, best-of-N milliseconds by worker count")
    print()
    print("fixture".ljust(width) + "".join(head.rjust(11) for head in headers), end="")
    print("   vs 4 workers" if not args.sweep else "")
    for name in names:
        row = name.ljust(width)
        row += "".join(f"{results[name, workers]:11.1f}" for workers in worker_counts)
        if not args.sweep:
            row += f"{results[name, 4] / results[name, None]:14.3f}x"
        print(row)
    print()
    for name in names:
        print(f"  {name:<17}{FIXTURES[name]}")


if __name__ == "__main__":
    sys.exit(main())
