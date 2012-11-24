#!/usr/bin/env python
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

import shutil
import os
import argparse
import tempfile

from pandas import DataFrame

DEFAULT_MIN_DURATION = 0.01
BASELINE_COMMIT = 'bdbca8e3dc' # 9,1 + regression fix # TODO: detect upstream/master

parser = argparse.ArgumentParser(description='Use vbench to generate a report comparing performance between two commits.')
parser.add_argument('-a', '--auto',
                    help='Execute a run using the defaults for the base and target commits.',
                    action='store_true',
                    default=False)
parser.add_argument('-b','--base-commit',
                    help='The commit serving as performance baseline (default: %s).' % BASELINE_COMMIT,
                    type=str)
parser.add_argument('-t','--target-commit',
                    help='The commit to compare against the baseline (default: HEAD).',
                    type=str)
parser.add_argument('-m', '--min-duration',
                    help='Minimum duration (in ms) of baseline test for inclusion in report (default: %.3f).' % DEFAULT_MIN_DURATION,
                    type=float,
                    default=0.01)
parser.add_argument('-o', '--output',
                    metavar="<file>",
                    dest='log_file',
                    help='path of file in which to save the report (default: vb_suite.log).')
args = parser.parse_args()

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

    from vbench.api import BenchmarkRunner
    from vbench.db import BenchmarkDB
    from vbench.git import GitRepo
    from suite import REPO_PATH, BUILD, DB_PATH, PREPARE, dependencies, benchmarks

    if not args.base_commit:
        args.base_commit = BASELINE_COMMIT

    # GitRepo wants exactly 7 character hash?
    args.base_commit = args.base_commit[:7]
    if args.target_commit:
        args.target_commit = args.target_commit[:7]

    if not args.log_file:
        args.log_file = os.path.abspath(os.path.join(REPO_PATH, 'vb_suite.log'))

    TMP_DIR =  tempfile.mkdtemp()
    prprint("TMP_DIR = %s" % TMP_DIR)
    prprint("LOG_FILE = %s\n" % args.log_file)

    try:
        logfile = open(args.log_file, 'w')

        prprint( "Processing Repo at '%s'..." % REPO_PATH)
        repo = GitRepo(REPO_PATH)

        # get hashes of baseline and current head
        h_head = args.target_commit or repo.shas[-1]
        h_baseline = args.base_commit

        prprint( "Opening DB at '%s'...\n" % DB_PATH)
        db = BenchmarkDB(DB_PATH)

        prprint('Target [%s] : %s' % (h_head, repo.messages.get(h_head,"")))
        prprint('Baseline [%s] : %s\n' % (h_baseline,repo.messages.get(h_baseline,"")))

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
        prprint( "Running benchmarks for baseline [%s]" % h_baseline)
        runner._run_and_write_results(h_baseline)

        prprint ("Running benchmarks for target [%s]" % h_head)
        runner._run_and_write_results(h_head)

        prprint( 'Processing results...')

        head_res = get_results_df(db,h_head)
        baseline_res = get_results_df(db,h_baseline)
        ratio = head_res['timing']/baseline_res['timing']
        totals = DataFrame(dict(t_head=head_res['timing'],
                                t_baseline=baseline_res['timing'],
                                ratio=ratio,
                                name=baseline_res.name),columns=["t_head","t_baseline","ratio","name"])
        totals = totals.ix[totals.t_head > args.min_duration] # ignore below threshold
        totals = totals.dropna().sort("ratio").set_index('name') # sort in ascending order

        s = "\n\nResults:\n" + totals.to_string(float_format=lambda x: "%0.4f" %x) + "\n\n"
        s += "Columns: test_name | target_duration [ms] | baseline_duration [ms] | ratio\n\n"
        s += "- a Ratio of 1.30 means the target commit is 30% slower then the baseline.\n\n"

        s += 'Target [%s] : %s' % (h_head, repo.messages.get(h_head,""))
        s += 'Baseline [%s] : %s\n\n' % (h_baseline,repo.messages.get(h_baseline,""))

        logfile.write(s)
        logfile.close()

        prprint(s )
        prprint("Results were also written to the logfile at '%s'\n" % args.log_file)

    finally:
        #        print("Disposing of TMP_DIR: %s" % TMP_DIR)
        shutil.rmtree(TMP_DIR)
        logfile.close()

if __name__ == '__main__':
    if not args.auto and not args.base_commit and not args.target_commit:
        parser.print_help()
    else:
        main()
