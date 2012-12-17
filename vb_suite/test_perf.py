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
2) instantiate a vbench runner, using the local repo as the source repo.
3) perform a vbench run for the baseline commit, then the target commit.
4) pull the results for both commits from the db. use pandas to align
everything and calculate a ration for the timing information.
5) print the results to the log file and to stdout.

"""

import shutil
import os
import argparse
import tempfile
import time

DEFAULT_MIN_DURATION = 0.01
BASELINE_COMMIT = '2149c50' # 0.9.1 + regression fix + vb fixes # TODO: detect upstream/master

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

def get_results_df(db,rev):
    from pandas import DataFrame
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
    from pandas import DataFrame
    from vbench.api import BenchmarkRunner
    from vbench.db import BenchmarkDB
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

        prprint( "Opening DB at '%s'...\n" % DB_PATH)
        db = BenchmarkDB(DB_PATH)

        prprint("Initializing Runner...")
        runner = BenchmarkRunner(benchmarks, REPO_PATH, REPO_PATH, BUILD, DB_PATH,
                                 TMP_DIR, PREPARE, always_clean=True,
            #                             run_option='eod', start_date=START_DATE,
                                 module_dependencies=dependencies)

        repo = runner.repo #(steal the parsed git repo used by runner)

        # ARGH. reparse the repo, without discarding any commits,
        # then overwrite the previous parse results
        #prprint ("Slaughtering kittens..." )
        (repo.shas, repo.messages,
         repo.timestamps, repo.authors) = _parse_commit_log(REPO_PATH)

        h_head = args.target_commit or repo.shas[-1]
        h_baseline = args.base_commit

        prprint('Target [%s] : %s\n' % (h_head, repo.messages.get(h_head,"")))
        prprint('Baseline [%s] : %s\n' % (h_baseline,repo.messages.get(h_baseline,"")))

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

        s = "\n\nResults:\n"
        s += totals.to_string(float_format=lambda x: "{:4.4f}".format(x).rjust(10))
        s += "\n\n"
        s += "Columns: test_name | target_duration [ms] | baseline_duration [ms] | ratio\n\n"
        s += "- a Ratio of 1.30 means the target commit is 30% slower then the baseline.\n\n"

        s += 'Target [%s] : %s\n' % (h_head, repo.messages.get(h_head,""))
        s += 'Baseline [%s] : %s\n\n' % (h_baseline,repo.messages.get(h_baseline,""))

        logfile.write(s)
        logfile.close()

        prprint(s )
        prprint("Results were also written to the logfile at '%s'\n" % args.log_file)

    finally:
        #        print("Disposing of TMP_DIR: %s" % TMP_DIR)
        shutil.rmtree(TMP_DIR)
        logfile.close()


# hack , vbench.git ignores some commits, but we
# need to be able to reference any commit.
# modified from vbench.git
def _parse_commit_log(repo_path):
    from vbench.git import parser, _convert_timezones
    from pandas import Series
    git_cmd = 'git --git-dir=%s/.git --work-tree=%s ' % (repo_path, repo_path)
    githist = git_cmd + ('log --graph --pretty=format:'
                          '\"::%h::%cd::%s::%an\" > githist.txt')
    os.system(githist)
    githist = open('githist.txt').read()
    os.remove('githist.txt')

    shas = []
    timestamps = []
    messages = []
    authors = []
    for line in githist.split('\n'):
        if '*' not in line.split("::")[0]: # skip non-commit lines
            continue

        _, sha, stamp, message, author = line.split('::', 4)

        # parse timestamp into datetime object
        stamp = parser.parse(stamp)

        shas.append(sha)
        timestamps.append(stamp)
        messages.append(message)
        authors.append(author)

    # to UTC for now
    timestamps = _convert_timezones(timestamps)

    shas = Series(shas, timestamps)
    messages = Series(messages, shas)
    timestamps = Series(timestamps, shas)
    authors = Series(authors, shas)
    return shas[::-1], messages[::-1], timestamps[::-1], authors[::-1]


if __name__ == '__main__':
    args = parser.parse_args()
    if not args.auto and not args.base_commit and not args.target_commit:
        parser.print_help()
    else:
        main()
