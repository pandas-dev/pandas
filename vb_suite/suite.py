from vbench.api import Benchmark, GitRepo, BenchmarkRunner
from datetime import datetime

import os

modules = ['groupby', 'indexing', 'reindex', 'binary_ops',
           'sparse', 'index_object']

benchmarks = []
for modname in modules:
    ref = __import__(modname)
    for k, v in ref.__dict__.iteritems():
        if isinstance(v, Benchmark):
            benchmarks.append(v)

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
dependencies = ['pandas_vb_common.py']

START_DATE = datetime(2011, 3, 1)

repo = GitRepo(REPO_PATH)

RST_BASE = '../doc/source'

# HACK!

timespan = [datetime(2011, 1, 1), datetime(2012, 1, 1)]

def generate_rst_files(benchmarks):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    vb_path = os.path.join(RST_BASE, 'vbench')
    fig_base_path = os.path.join(vb_path, 'figures')

    if not os.path.exists(vb_path):
        print 'creating %s' % vb_path
        os.makedirs(vb_path)

    if not os.path.exists(fig_base_path):
        print 'creating %s' % fig_base_path
        os.makedirs(fig_base_path)

    for bmk in benchmarks:
        print 'Generating rst file for %s' % bmk.name
        rst_path = os.path.join(RST_BASE, 'vbench/%s.rst' % bmk.name)

        fig_full_path = os.path.join(fig_base_path, '%s.png' % bmk.name)

        # make the figure
        plt.figure(figsize=(10, 6))
        ax = plt.gca()
        bmk.plot(DB_PATH, ax=ax)
        plt.xlim(timespan)
        plt.savefig(fig_full_path, bbox_inches='tight')
        plt.close('all')

        fig_rel_path = 'vbench/figures/%s.png' % bmk.name
        rst_text = bmk.to_rst(image_path=fig_rel_path)
        with open(rst_path, 'w') as f:
            f.write(rst_text)

    with open(os.path.join(RST_BASE, 'vbench.rst'), 'w') as f:
        print >> f, """
VBENCH
------
"""
        for bmk in benchmarks:
            print >> f, '.. include:: vbench/%s.rst' % bmk.name
