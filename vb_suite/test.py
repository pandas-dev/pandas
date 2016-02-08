from pandas import *
import matplotlib.pyplot as plt

import sqlite3

from vbench.git import GitRepo


REPO_PATH = '/home/adam/code/pandas'
repo = GitRepo(REPO_PATH)

con = sqlite3.connect('vb_suite/benchmarks.db')

bmk = '36900a889961162138c140ce4ae3c205'
# bmk = '9d7b8c04b532df6c2d55ef497039b0ce'
bmk = '4481aa4efa9926683002a673d2ed3dac'
bmk = '00593cd8c03d769669d7b46585161726'
bmk = '3725ab7cd0a0657d7ae70f171c877cea'
bmk = '3cd376d6d6ef802cdea49ac47a67be21'
bmk2 = '459225186023853494bc345fd180f395'
bmk = 'c22ca82e0cfba8dc42595103113c7da3'
bmk = 'e0e651a8e9fbf0270ab68137f8b9df5f'
bmk = '96bda4b9a60e17acf92a243580f2a0c3'


def get_results(bmk):
    results = con.execute(
        "select * from results where checksum='%s'" % bmk).fetchall()
    x = Series(dict((t[1], t[3]) for t in results))
    x.index = x.index.map(repo.timestamps.get)
    x = x.sort_index()
    return x

x = get_results(bmk)


def graph1():
    dm_getitem = get_results('459225186023853494bc345fd180f395')
    dm_getvalue = get_results('c22ca82e0cfba8dc42595103113c7da3')

    plt.figure()
    ax = plt.gca()

    dm_getitem.plot(label='df[col][idx]', ax=ax)
    dm_getvalue.plot(label='df.get_value(idx, col)', ax=ax)

    plt.ylabel('ms')
    plt.legend(loc='best')


def graph2():
    bm = get_results('96bda4b9a60e17acf92a243580f2a0c3')
    plt.figure()
    ax = plt.gca()

    bm.plot(ax=ax)
    plt.ylabel('ms')

bm = get_results('36900a889961162138c140ce4ae3c205')
fig = plt.figure()
ax = plt.gca()
bm.plot(ax=ax)
fig.autofmt_xdate()

plt.xlim([bm.dropna().index[0] - datetools.MonthEnd(),
          bm.dropna().index[-1] + datetools.MonthEnd()])
plt.ylabel('ms')
