from dateutil import parser
import subprocess
import os
import re
import sys

import numpy as np

from pandas import *

repo_path = '/home/wesm/code/pandas'
githist = ('git log --pretty=format:\"%h %ad | %s%d [%an]\" --date=short ' +
           repo_path + ' > githist.txt')

def rungithist():
    os.system(githist)

def get_commit_history():
    # return TimeSeries

    rungithist()

    githist = open('githist.txt').read()
    os.remove('githist.txt')

    sha_date = []
    for line in githist.split('\n'):
        sha_date.append(line.split()[:2])

    shas, dates = zip(*sha_date)

    hists = dict(zip(shas, githist.split('\n')))

    dates = [parser.parse(d) for d in dates]

    return Series(dates, shas), hists

def get_commit_churn(sha, prev_sha):
    stdout = subprocess.Popen(['git', 'diff', sha, prev_sha, '--numstat'],
                              stdout=subprocess.PIPE).stdout

    stdout = stdout.read()

    insertions = {}
    deletions = {}

    for line in stdout.split('\n'):
        try:
            i, d, path = line.split('\t')
            insertions[path] = int(i)
            deletions[path] = int(d)
        except: # EAFP
            pass

    # statline = stdout.split('\n')[-2]

    # match = re.match('.*\s(.*)\sinsertions.*\s(.*)\sdeletions', statline)

    # insertions = int(match.group(1))
    # deletions = int(match.group(2))

    return insertions, deletions

def get_code_churn(commits):
    shas = commits.index[::-1]

    prev = shas[0]

    insertions = [np.nan]
    deletions = [np.nan]

    insertions = {}
    deletions = {}

    for cur in shas[1:]:
        i, d = get_commit_churn(cur, prev)

        insertions[cur] = i
        deletions[cur] = d

        # insertions.append(i)
        # deletions.append(d)

        prev = cur

    return Panel({'insertions' : DataFrame(insertions),
                  'deletions' : DataFrame(deletions)}, minor_axis=shas)


    # return DataFrame({'insertions' : insertions,
    #                   'deletions' : deletions}, index=shas)

if __name__ == '__main__':
    commits, hists = get_commit_history()
    churn = get_code_churn(commits)

    file_include = []
    for path in churn.major_axis:
        if path.endswith('.pyx') or path.endswith('.py'):
            file_include.append(path)
    commits_include = [sha for sha in churn.minor_axis
                       if 'LF' not in hists[sha]]
    commits_include.remove('dcf3490')

    clean_churn = churn.reindex(major=file_include, minor=commits_include)

    by_commit = clean_churn.sum('major').sum(1)

    by_date = by_commit.groupby(commits).sum()

    by_date = by_date.drop([datetime(2011, 6, 10)])

    # clean out days where I touched Cython

    by_date = by_date[by_date < 5000]
