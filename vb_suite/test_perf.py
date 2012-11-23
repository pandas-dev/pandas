#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
What
----
vbench is a library which can be used to benchmark the performance
of a codebase over time.
Although vbench can collect data over many commites, generate plots
and other niceties, for Pull-Requests the important thing is the
performance of the HEAD commit against a known-good baseline.

This script tries to automate the process of comparing these
two commits, and is meant to run out of the box on a fresh
clone.

How
---
These are the steps taken:
1) create a temp directory into which vbench will clone the temporary repo.
2) parse the Git tree to obtain metadata, and determine the HEAD.
3) instantiate a vbench runner, using the local repo as the source repo.
4) If results for the BASELINE_COMMIT aren't already in the db, have vbench
do a run for it and store the results.
5) perform a vbench run for HEAD and store the results.
6) pull the results for both commits from the db. use pandas to align
everything and calculate a ration for the timing information.
7) print the results to the log file and to stdout.

Known Issues: vbench fails to locate a baseline if HEAD is not a descendent
"""
import sys
import shutil

from pandas import *
from vbench.api import BenchmarkRunner
from vbench.db import BenchmarkDB
from vbench.git import GitRepo
import tempfile

from suite import *

BASELINE_COMMIT = 'bdbca8e'  # v0.9,1 + regression fix
LOG_FILE = os.path.abspath(os.path.join(REPO_PATH, 'vb_suite.log'))

def get_results_df(db,rev):
    """Takes a git commit hash and returns a Dataframe of benchmark results
    """
    bench = DataFrame(db.get_benchmarks())
    results = DataFrame(db.get_rev_results(rev).values())

    # Sinch vbench.db._reg_rev_results returns an unlabeled dict,
    # we have to break encapsulation a bit.
    results.columns = db._results.c.keys()
    results = results.join(bench['name'], on='checksum').set_index("checksum")
    return results

def prprint(s):
    print("*** %s"%s)

def main():
    TMP_DIR =  tempfile.mkdtemp()
    prprint("TMP_DIR = %s" % TMP_DIR)
    prprint("LOG_FILE = %s\n" % LOG_FILE)

    try:
        logfile = open(LOG_FILE, 'w')

        prprint( "Processing Repo at '%s'..." % REPO_PATH)
        repo = GitRepo(REPO_PATH)

        # get hashes of baseline and current head
        h_head = repo.shas[-1]
        h_baseline = BASELINE_COMMIT

        prprint( "Opening DB at '%s'...\n" % DB_PATH)
        db = BenchmarkDB(DB_PATH)

        prprint( 'Comparing Head [%s] : %s ' % (h_head, repo.messages.get(h_head,"")))
        prprint( 'Against baseline [%s] : %s \n' % (h_baseline,
                repo.messages.get(h_baseline,"")))

        prprint("Initializing Runner...")
        runner = BenchmarkRunner(benchmarks, REPO_PATH, REPO_PATH, BUILD, DB_PATH,
                                 TMP_DIR, PREPARE, always_clean=True,
            #                             run_option='eod', start_date=START_DATE,
                                 module_dependencies=dependencies)

        prprint ("removing any previous measurements for the commits." )
        db.delete_rev_results(h_baseline)
        db.delete_rev_results(h_head)

        # TODO: we could skip this, but we need to make sure all
        # results are in the DB, which is a little tricky with
        # start dates and so on.
        prprint( "Running benchmarks for baseline commit '%s'" % h_baseline)
        runner._run_and_write_results(h_baseline)

        prprint ("Running benchmarks for current HEAD '%s'" % h_head)
        runner._run_and_write_results(h_head)

        prprint( 'Processing results...')

        head_res = get_results_df(db,h_head)
        baseline_res = get_results_df(db,h_baseline)
        ratio = head_res['timing']/baseline_res['timing']
        totals = DataFrame(dict(t_head=head_res['timing'],
                                t_baseline=baseline_res['timing'],
                                ratio=ratio,
                                name=baseline_res.name),columns=["t_head","t_baseline","ratio","name"])
        totals = totals.ix[totals.t_head > 1.0] # ignore sub 1ms
        totals = totals.dropna().sort("ratio").set_index('name') # sort in ascending order

        s = "\n\nResults:\n" + totals.to_string(float_format=lambda x: "%0.2f" %x) + "\n\n"
        s += "Columns: test_name | head_time [ms] | baseline_time [ms] | ratio\n\n"
        s += "- a Ratio of 1.30 means HEAD is 30% slower then the Baseline.\n\n"

        s += 'Head [%s] : %s\n' % (h_head, repo.messages.get(h_head,""))
        s += 'Baseline [%s] : %s\n\n' % (h_baseline,repo.messages.get(h_baseline,""))

        logfile.write(s)
        logfile.close()

        prprint(s )
        prprint("Results were also written to the logfile at '%s'\n" % LOG_FILE)

    finally:
        #        print("Disposing of TMP_DIR: %s" % TMP_DIR)
        shutil.rmtree(TMP_DIR)
        logfile.close()

if __name__ == '__main__':
    main()
