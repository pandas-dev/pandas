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
import sys
import argparse
import tempfile
import time
import re

import random
import numpy as np

from pandas import DataFrame, Series

from suite import REPO_PATH

DEFAULT_MIN_DURATION = 0.01
HEAD_COL="head[ms]"
BASE_COL="base[ms]"


class RevParseAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        import subprocess
        cmd = 'git rev-parse {0}'.format(values)
        rev_parse = subprocess.check_output(cmd, shell=True)
        setattr(namespace, self.dest, rev_parse.strip())


parser = argparse.ArgumentParser(description='Use vbench to measure and compare the performance of commits.')
parser.add_argument('-H', '--head',
                    help='Execute vbenches using the currently checked out copy.',
                    dest='head',
                    action='store_true',
                    default=False)
parser.add_argument('-b', '--base-commit',
                    help='The commit serving as performance baseline ',
                    type=str, action=RevParseAction)
parser.add_argument('-t', '--target-commit',
                    help='The commit to compare against the baseline (default: HEAD).',
                    type=str, action=RevParseAction)
parser.add_argument('-m', '--min-duration',
                    help='Minimum duration (in ms) of baseline test for inclusion in report (default: %.3f).' % DEFAULT_MIN_DURATION,
                    type=float,
                    default=0.01)
parser.add_argument('-o', '--output',
                    metavar="<file>",
                    dest='log_file',
                    help='Path of file in which to save the textual report (default: vb_suite.log).')
parser.add_argument('-d', '--outdf',
                    metavar="FNAME",
                    dest='outdf',
                    default=None,
                    help='Name of file to df.save() the result table into. Will overwrite')
parser.add_argument('-r', '--regex',
                    metavar="REGEX",
                    dest='regex',
                    default="",
                    help='Regex pat, only tests whose name matches the regext will be run.')
parser.add_argument('-s', '--seed',
                    metavar="SEED",
                    dest='seed',
                    default=1234,
                    type=int,
                    help='Integer value to seed PRNG with')
parser.add_argument('-n', '--repeats',
                    metavar="N",
                    dest='repeats',
                    default=3,
                    type=int,
                    help='Number of times to run each vbench, result value is the best of')
parser.add_argument('-c', '--ncalls',
                    metavar="N",
                    dest='ncalls',
                    default=3,
                    type=int,
                    help='Number of calls to in each repetition of a vbench')
parser.add_argument('-N', '--hrepeats',
                    metavar="N",
                    dest='hrepeats',
                    default=1,
                    type=int,
                    help='implies -H, number of times to run the vbench suite on the head commit.\n'
                    'Each iteration will yield another column in the output' )
parser.add_argument('-a', '--affinity',
                    metavar="a",
                    dest='affinity',
                    default=1,
                    type=int,
                    help='set processor affinity of process by default bind to cpu/core #1 only. '
                         'Requires the "affinity" or "psutil" python module, will raise Warning otherwise')
parser.add_argument('-u', '--burnin',
                    metavar="u",
                    dest='burnin',
                    default=1,
                    type=int,
                    help='Number of extra iteration per benchmark to perform first, then throw away. '  )

parser.add_argument('-S', '--stats',
                    default=False,
                    action='store_true',
                    help='when specified with -N, prints the output of describe() per vbench results. '  )

parser.add_argument('-q', '--quiet',
                    default=False,
                    action='store_true',
                    help='Suppress report output to stdout. '  )

def get_results_df(db, rev):
    """Takes a git commit hash and returns a Dataframe of benchmark results
    """
    bench = DataFrame(db.get_benchmarks())
    results = DataFrame(map(list,db.get_rev_results(rev).values()))

    # Sinch vbench.db._reg_rev_results returns an unlabeled dict,
    # we have to break encapsulation a bit.
    results.columns = db._results.c.keys()
    results = results.join(bench['name'], on='checksum').set_index("checksum")
    return results


def prprint(s):
    print("*** %s" % s)

def pre_hook():
    import gc
    gc.disable()

def post_hook():
    import gc
    gc.enable()

def profile_comparative(benchmarks):

    from vbench.api import BenchmarkRunner
    from vbench.db import BenchmarkDB
    from vbench.git import GitRepo
    from suite import BUILD, DB_PATH, PREPARE, dependencies
    TMP_DIR = tempfile.mkdtemp()

    try:

        prprint("Opening DB at '%s'...\n" % DB_PATH)
        db = BenchmarkDB(DB_PATH)

        prprint("Initializing Runner...")

        # all in a good cause...
        GitRepo._parse_commit_log = _parse_wrapper(args.base_commit)

        runner = BenchmarkRunner(
            benchmarks, REPO_PATH, REPO_PATH, BUILD, DB_PATH,
            TMP_DIR, PREPARE, always_clean=True,
            # run_option='eod', start_date=START_DATE,
            module_dependencies=dependencies)

        repo = runner.repo  # (steal the parsed git repo used by runner)
        h_head = args.target_commit or repo.shas[-1]
        h_baseline = args.base_commit

        # ARGH. reparse the repo, without discarding any commits,
        # then overwrite the previous parse results
        # prprint ("Slaughtering kittens..." )
        (repo.shas, repo.messages,
         repo.timestamps, repo.authors) = _parse_commit_log(None,REPO_PATH,
                                                                args.base_commit)

        prprint('Target [%s] : %s\n' % (h_head, repo.messages.get(h_head, "")))
        prprint('Baseline [%s] : %s\n' % (h_baseline,
                repo.messages.get(h_baseline, "")))

        prprint("Removing any previous measurements for the commits.")
        db.delete_rev_results(h_baseline)
        db.delete_rev_results(h_head)

        # TODO: we could skip this, but we need to make sure all
        # results are in the DB, which is a little tricky with
        # start dates and so on.
        prprint("Running benchmarks for baseline [%s]" % h_baseline)
        runner._run_and_write_results(h_baseline)

        prprint("Running benchmarks for target [%s]" % h_head)
        runner._run_and_write_results(h_head)

        prprint('Processing results...')

        head_res = get_results_df(db, h_head)
        baseline_res = get_results_df(db, h_baseline)
        ratio = head_res['timing'] / baseline_res['timing']
        totals = DataFrame({HEAD_COL:head_res['timing'],
                                BASE_COL:baseline_res['timing'],
                                'ratio':ratio,
                                'name':baseline_res.name},
                                columns=[HEAD_COL, BASE_COL, "ratio", "name"])
        totals = totals.ix[totals[HEAD_COL] > args.min_duration]
            # ignore below threshold
        totals = totals.dropna(
        ).sort("ratio").set_index('name')  # sort in ascending order

        h_msg =   repo.messages.get(h_head, "")
        b_msg =   repo.messages.get(h_baseline, "")

        print_report(totals,h_head=h_head,h_msg=h_msg,
                     h_baseline=h_baseline,b_msg=b_msg)

        if args.outdf:
            prprint("The results DataFrame was written to '%s'\n" %  args.outdf)
            totals.save(args.outdf)
    finally:
        #        print("Disposing of TMP_DIR: %s" % TMP_DIR)
        shutil.rmtree(TMP_DIR)


def profile_head_single(benchmark):
    import gc
    results = []

    # just in case
    gc.collect()

    try:
        from ctypes import cdll, CDLL
        cdll.LoadLibrary("libc.so.6")
        libc = CDLL("libc.so.6")
        libc.malloc_trim(0)
    except:
        pass


    N =  args.hrepeats + args.burnin

    results = []
    try:
        for i in range(N):
            gc.disable()
            d=dict()

            try:
                d = benchmark.run()

            except KeyboardInterrupt:
                raise
            except Exception as e: # if a single vbench bursts into flames, don't die.
                err=""
                try:
                    err =  d.get("traceback")
                    if err is None:
                        err = str(e)
                except:
                    pass
                print("%s died with:\n%s\nSkipping...\n" % (benchmark.name, err))

            results.append(d.get('timing',np.nan))
            gc.enable()
            gc.collect()

    finally:
        gc.enable()

    if results:
        # throw away the burn_in
        results = results[args.burnin:]
    sys.stdout.write('.')
    sys.stdout.flush()
    return Series(results, name=benchmark.name)

    # df = DataFrame(results)
    # df.columns = ["name",HEAD_COL]
    # return df.set_index("name")[HEAD_COL]

def profile_head(benchmarks):
    print( "Performing %d benchmarks (%d runs each)" % ( len(benchmarks), args.hrepeats))

    ss= [profile_head_single(b) for b in benchmarks]
    print("\n")

    results = DataFrame(ss)
    results.columns=[ "#%d" %i for i in range(args.hrepeats)]
    # results.index = ["#%d" % i for i in range(len(ss))]
    # results = results.T

    shas, messages, _,_  = _parse_commit_log(None,REPO_PATH,base_commit="HEAD^")
    print_report(results,h_head=shas[-1],h_msg=messages[-1])


    if args.outdf:
        prprint("The results DataFrame was written to '%s'\n" %  args.outdf)
        DataFrame(results).save(args.outdf)

def print_report(df,h_head=None,h_msg="",h_baseline=None,b_msg=""):

    name_width=45
    col_width = 10

    hdr = ("{:%s}" % name_width).format("Test name")
    hdr += ("|{:^%d}"  % col_width)* len(df.columns)
    hdr += "|"
    hdr = hdr.format(*df.columns)
    hdr = "-"*len(hdr) + "\n" + hdr + "\n" + "-"*len(hdr) + "\n"
    ftr=hdr
    s = "\n"
    s+= "Invoked with :\n"
    s+= "--ncalls: %s\n" % (args.ncalls or 'Auto')
    s+= "--repeats: %s\n" % (args.repeats)
    s+= "\n\n"

    s += hdr
    # import ipdb
    # ipdb.set_trace()
    for i in range(len(df)):
        lfmt = ("{:%s}" % name_width)
        lfmt += ("| {:%d.4f} " % (col_width-2))* len(df.columns)
        lfmt += "|\n"
        s += lfmt.format(df.index[i],*list(df.irow(i).values))

    s+= ftr + "\n"

    s += "Ratio < 1.0 means the target commit is faster then the baseline.\n"
    s += "Seed used: %d\n\n" % args.seed

    if  h_head:
        s += 'Target [%s] : %s\n' % (h_head, h_msg)
    if  h_baseline:
        s += 'Base   [%s] : %s\n\n' % (
            h_baseline, b_msg)

    stats_footer = "\n"
    if args.stats :
        stats_footer += str(df.T.describe().T) + "\n\n"

    s+= stats_footer
    logfile = open(args.log_file, 'w')
    logfile.write(s)
    logfile.close()

    if not args.quiet:
        prprint(s)

    if args.stats and args.quiet:
        prprint(stats_footer)

    prprint("Results were also written to the logfile at '%s'" %
            args.log_file)



def main():
    from suite import benchmarks
    # GitRepo wants exactly 7 character hash?
    if args.base_commit:
        args.base_commit = args.base_commit[:7]
    if args.target_commit:
        args.target_commit = args.target_commit[:7]

    if not args.log_file:
        args.log_file = os.path.abspath(
            os.path.join(REPO_PATH, 'vb_suite.log'))

    saved_dir = os.path.curdir
    if args.outdf:
        # not bullet-proof but enough for us
        args.outdf = os.path.realpath(args.outdf)

    if args.log_file:
        # not bullet-proof but enough for us
        args.log_file = os.path.realpath(args.log_file)

    random.seed(args.seed)
    np.random.seed(args.seed)

    affinity_set = False

    # try psutil first since it is more commonly present and better
    # maintained.  Some people experienced problems with affinity package
    # (see https://code.google.com/p/psutil/issues/detail?id=238 for more references)
    try:
        import psutil
        if hasattr(psutil.Process, 'set_cpu_affinity'):
            psutil.Process(os.getpid()).set_cpu_affinity([args.affinity])
            affinity_set = True
    except ImportError:
        pass

    if not affinity_set:
        try:
            import affinity
            affinity.set_process_affinity_mask(0, args.affinity)
            assert affinity.get_process_affinity_mask(0) == args.affinity
            affinity_set = True
        except ImportError:
            pass

    if not affinity_set:
        import warnings
        warnings.warn("\n\n"
              "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
              "The 'affinity' or 'psutil' >= 0.5.0 modules are not available, results may be unreliable\n"
              "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
            )
        time.sleep(2)
    else:
        print("CPU affinity set to %d" % args.affinity)

    print("\n")
    prprint("LOG_FILE = %s" % args.log_file)
    if args.outdf:
        prprint("PICKE_FILE = %s" % args.outdf)

    print("\n")

    # move away from the pandas root dit, to avoid possible import
    # surprises
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    benchmarks = [x for x in benchmarks if re.search(args.regex,x.name)]

    for b in benchmarks:
        b.repeat = args.repeats
        if args.ncalls:
            b.ncalls = args.ncalls

    if benchmarks:
        if args.head:
            profile_head(benchmarks)
        else:
            profile_comparative(benchmarks)
    else:
        print( "No matching benchmarks")

    os.chdir(saved_dir)

# hack , vbench.git ignores some commits, but we
# need to be able to reference any commit.
# modified from vbench.git
def _parse_commit_log(this,repo_path,base_commit=None):
    from vbench.git import _convert_timezones
    from pandas import Series
    from dateutil import parser as dparser

    git_cmd = 'git --git-dir=%s/.git --work-tree=%s ' % (repo_path, repo_path)
    githist = git_cmd + ('log --graph --pretty=format:'+
                         '\"::%h::%cd::%s::%an\"'+
                         ('%s..' % base_commit)+
                         '> githist.txt')
    os.system(githist)
    githist = open('githist.txt').read()
    os.remove('githist.txt')

    shas = []
    timestamps = []
    messages = []
    authors = []
    for line in githist.split('\n'):
        if '*' not in line.split("::")[0]:  # skip non-commit lines
            continue

        _, sha, stamp, message, author = line.split('::', 4)

        # parse timestamp into datetime object
        stamp = dparser.parse(stamp)

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

# even worse, monkey patch vbench
def _parse_wrapper(base_commit):
    def inner(repo_path):
        return _parse_commit_log(repo_path,base_commit)
    return inner

if __name__ == '__main__':
    args = parser.parse_args()
    if not args.head and (not args.base_commit and not args.target_commit):
        parser.print_help()
    else:
        import warnings
        warnings.filterwarnings('ignore',category=FutureWarning)
        warnings.filterwarnings('ignore',category=DeprecationWarning)
        main()
