import sys

from ._aux import update_sys_path
from .discovery import disc_benchmarks


def _check(args):
    """
    Checks all the discovered benchmarks in the provided benchmark directory.

    #### Parameters
    **args** (`tuple`)
    : A tuple containing the benchmark directory.

    #### Notes
    This function updates the system path with the root directory of the
    benchmark suite. Then, it iterates over all benchmarks discovered in the
    root directory. For each benchmark, it calls the check method of the
    benchmark and updates the 'ok' flag.

    If all benchmarks pass the check, it exits with a status code 0. If any
    benchmark fails, it exits with a status code 1.
    """
    (benchmark_dir,) = args

    update_sys_path(benchmark_dir)

    ok = True
    for benchmark in disc_benchmarks(benchmark_dir):
        ok = ok and benchmark.check(benchmark_dir)

    sys.exit(0 if ok else 1)
