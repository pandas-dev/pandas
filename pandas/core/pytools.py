"""A collection of tools for various purely Python operations"""
from random import Random
import base64
import os
import string

# In theory should be few to no imports outside perhaps stdlib here

def rands(n):
    """Generates a random alphanumeric string of length *n*"""
    return ''.join(Random().sample(string.letters+string.digits, n))

def adjoin(space, *lists):
    """
    Glues together two sets of strings using the amount of space requested.
    The idea is to prettify.
    """
    outLines = []
    newLists = []
    lengths = [max(map(len, x)) + space for x in lists]
    maxLen = max(map(len, lists))
    for i, lst in enumerate(lists):
        nl = [x.ljust(lengths[i]) for x in lst]
        nl.extend([' ' * lengths[i]] * (maxLen - len(lst)))
        newLists.append(nl)
    toJoin = zip(*newLists)
    for lines in toJoin:
        outLines.append(''.join(lines))
    return '\n'.join(outLines)

def indent(string, spaces=4):
    dent = ' ' * spaces
    return '\n'.join([dent + x for x in string.split('\n')])

def banner(message):
    bar = '=' * 80
    return '%s\n%s\n%s' % (bar, message, bar)


class groupby(dict):
    """
    A simple groupby different from the one in itertools.
    
    Does not require the sequence elements to be sorted by keys,
    however it is slower. 
    """
    def __init__(self, seq, key=lambda x:x):
        for value in seq:
            k = key(value)
            self.setdefault(k, []).append(value)
    __iter__ = dict.iteritems

    
def map_indices_py(arr):
    """
    Returns a dictionary with (element, index) pairs for each element in the 
    given array/list
    """
    return dict([(x, i) for i, x in enumerate(arr)])

