#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Short one-line summary

long summary
"""


def main():
    import shutil
    import tempfile
    import warnings

    from pandas import Series

    from vbench.api import BenchmarkRunner
    from suite import (REPO_PATH, BUILD, DB_PATH, PREPARE,
                       dependencies, benchmarks)

    from memory_profiler import memory_usage

    warnings.filterwarnings('ignore', category=FutureWarning)

    try:
        TMP_DIR = tempfile.mkdtemp()
        runner = BenchmarkRunner(
            benchmarks, REPO_PATH, REPO_PATH, BUILD, DB_PATH,
            TMP_DIR, PREPARE, always_clean=True,
            # run_option='eod', start_date=START_DATE,
            module_dependencies=dependencies)
        results = {}
        for b in runner.benchmarks:
            k = b.name
            try:
                vs = memory_usage((b.run,))
                v = max(vs)
                # print(k, v)
                results[k] = v
            except Exception as e:
                print("Exception caught in %s\n" % k)
                print(str(e))

        s = Series(results)
        s.sort()
        print((s))

    finally:
        shutil.rmtree(TMP_DIR)


if __name__ == "__main__":
    main()
