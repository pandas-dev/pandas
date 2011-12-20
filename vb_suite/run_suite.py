from vbench.api import Benchmark, GitRepo, BenchmarkRunner
from datetime import datetime

modules = ['groupby', 'indexing', 'reindex', 'binary_ops',
           'sparse']

all_benchmarks = []
for modname in modules:
    ref = __import__(modname)
    for k, v in ref.__dict__.iteritems():
        if isinstance(v, Benchmark):
            all_benchmarks.append(v)

REPO_PATH = '/home/wesm/code/pandas'
REPO_URL = 'git@github.com:wesm/pandas.git'
DB_PATH = '/home/wesm/code/pandas/vb_suite/benchmarks.db'
TMP_DIR = '/home/wesm/tmp/vb_pandas'
PREPARE = """
python setup.py clean
"""
BUILD = """
python setup.py build_ext --inplace
"""
START_DATE = datetime(2011, 3, 1)

repo = GitRepo(REPO_PATH)

to_consider = repo.shas.truncate(START_DATE)

runner = BenchmarkRunner(all_benchmarks, REPO_PATH, REPO_URL,
                         BUILD, DB_PATH, TMP_DIR, PREPARE,
                         run_option='eod', start_date=START_DATE)

runner.run()
