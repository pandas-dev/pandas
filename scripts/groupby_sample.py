from pandas import *
import numpy as np
import string

g1 = np.array(list(string.letters))[:-1]
g2 = np.arange(51)
df_small = DataFrame({'group1' : ["a","b","a","a","b","c","c","c","c",
                                  "c","a","a","a","b","b","b","b"],
                      'group2' : [1,2,3,4,1,3,5,6,5,4,1,2,3,4,3,2,1],
                      'value'  : ["apple","pear","orange","apple",
                                  "banana","durian","lemon","lime",
                                  "raspberry","durian","peach","nectarine",
                                  "banana","lemon","guava","blackberry",
                                  "grape"]})
value = df_small['value'].values.repeat(3)
df = DataFrame({'group1' : g1.repeat(40000),
                'group2' : np.tile(g2, 40000),
                'value' : value.repeat(40000)})


def random_sample():
    grouped = df.groupby(['group1','group2'])['value']
    from random import choice
    choose = lambda group: choice(group.index)
    indices = grouped.apply(choose)
    return df.reindex(indices)
