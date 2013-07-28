from pandas import *
import numpy as np
import string
import pandas.compat as compat

g1 = np.array(list(string.letters))[:-1]
g2 = np.arange(510)
df_small = DataFrame({'group1': ["a", "b", "a", "a", "b", "c", "c", "c", "c",
                                 "c", "a", "a", "a", "b", "b", "b", "b"],
                      'group2': [1, 2, 3, 4, 1, 3, 5, 6, 5, 4, 1, 2, 3, 4, 3, 2, 1],
                      'value': ["apple", "pear", "orange", "apple",
                                "banana", "durian", "lemon", "lime",
                                "raspberry", "durian", "peach", "nectarine",
                                "banana", "lemon", "guava", "blackberry",
                                "grape"]})
value = df_small['value'].values.repeat(3)
df = DataFrame({'group1': g1.repeat(4000 * 5),
                'group2': np.tile(g2, 400 * 5),
                'value': value.repeat(4000 * 5)})


def random_sample():
    grouped = df.groupby(['group1', 'group2'])['value']
    from random import choice
    choose = lambda group: choice(group.index)
    indices = grouped.apply(choose)
    return df.reindex(indices)


def random_sample_v2():
    grouped = df.groupby(['group1', 'group2'])['value']
    from random import choice
    choose = lambda group: choice(group.index)
    indices = [choice(v) for k, v in compat.iteritems(grouped.groups)]
    return df.reindex(indices)


def do_shuffle(arr):
    from random import shuffle
    result = arr.copy().values
    shuffle(result)
    return result


def shuffle_uri(df, grouped):
    perm = np.r_[tuple([np.random.permutation(
        idxs) for idxs in compat.itervalues(grouped.groups)])]
    df['state_permuted'] = np.asarray(df.ix[perm]['value'])

df2 = df.copy()
grouped = df2.groupby('group1')
shuffle_uri(df2, grouped)

df2['state_perm'] = grouped['value'].transform(do_shuffle)
