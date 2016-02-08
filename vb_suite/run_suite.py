#!/usr/bin/env python
from vbench.api import BenchmarkRunner
from suite import *


def run_process():
    runner = BenchmarkRunner(benchmarks, REPO_PATH, REPO_URL,
                             BUILD, DB_PATH, TMP_DIR, PREPARE,
                             always_clean=True,
                             run_option='eod', start_date=START_DATE,
                             module_dependencies=dependencies)
    runner.run()

if __name__ == '__main__':
    run_process()
