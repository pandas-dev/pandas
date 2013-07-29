import subprocess
import os
import re
import sys

import numpy as np

from pandas import *


if __name__ == '__main__':
    from vbench.git import GitRepo
    repo = GitRepo('/Users/wesm/code/pandas')
    churn = repo.get_churn_by_file()

    file_include = []
    for path in churn.major_axis:
        if path.endswith('.pyx') or path.endswith('.py'):
            file_include.append(path)
    commits_include = [sha for sha in churn.minor_axis
                       if 'LF' not in repo.messages[sha]]
    commits_include.remove('dcf3490')

    clean_churn = churn.reindex(major=file_include, minor=commits_include)

    by_commit = clean_churn.sum('major').sum(1)

    by_date = by_commit.groupby(repo.commit_date).sum()

    by_date = by_date.drop([datetime(2011, 6, 10)])

    # clean out days where I touched Cython

    by_date = by_date[by_date < 5000]
