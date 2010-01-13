"""A collection of tools for various purely Python operations"""
from random import Random
import itertools
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
    lengths = [max(map(len, x)) + space for x in lists[:-1]]

    # not the last one
    lengths.append(max(map(len, lists[-1])))

    maxLen = max(map(len, lists))
    for i, lst in enumerate(lists):
        nl = [x.ljust(lengths[i]) for x in lst]
        nl.extend([' ' * lengths[i]] * (maxLen - len(lst)))
        newLists.append(nl)
    toJoin = zip(*newLists)
    for lines in toJoin:
        outLines.append(''.join(lines))
    return '\n'.join(outLines)


def iterpairs(seq):
    """
    Parameters
    ----------
    seq: sequence

    Returns
    -------
    iterator returning overlapping pairs of elements

    Example
    -------
    >>> iterpairs([1, 2, 3, 4])
    [(1, 2), (2, 3), (3, 4)
    """
    # input may not be sliceable
    seq_it = iter(seq)
    seq_it_next = iter(seq)
    _ = seq_it_next.next()

    return itertools.izip(seq_it, seq_it_next)

def indent(string, spaces=4):
    dent = ' ' * spaces
    return '\n'.join([dent + x for x in string.split('\n')])

def banner(message):
    """
    Return 80-char width message declaration with = bars on top and bottom.
    """
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

#===============================================================================
# Set operations
#===============================================================================

def union(*seqs):
    result = set([])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result |= seq
    return type(seqs[0])(list(result))

def difference(a, b):
    return type(a)(list(set(a) - set(b)))

def intersection(*seqs):
    result = set(seqs[0])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result &= seq
    return type(seqs[0])(list(result))
