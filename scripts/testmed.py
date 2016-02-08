## {{{ Recipe 576930 (r10): Efficient Running Median using an Indexable Skiplist

from random import random
from math import log, ceil
from pandas.compat import range
from numpy.random import randn
from pandas.lib.skiplist import rolling_median


class Node(object):
    __slots__ = 'value', 'next', 'width'

    def __init__(self, value, next, width):
        self.value, self.next, self.width = value, next, width


class End(object):
    'Sentinel object that always compares greater than another object'
    def __cmp__(self, other):
        return 1

NIL = Node(End(), [], [])               # Singleton terminator node


class IndexableSkiplist:
    'Sorted collection supporting O(lg n) insertion, removal, and lookup by rank.'

    def __init__(self, expected_size=100):
        self.size = 0
        self.maxlevels = int(1 + log(expected_size, 2))
        self.head = Node('HEAD', [NIL] * self.maxlevels, [1] * self.maxlevels)

    def __len__(self):
        return self.size

    def __getitem__(self, i):
        node = self.head
        i += 1
        for level in reversed(range(self.maxlevels)):
            while node.width[level] <= i:
                i -= node.width[level]
                node = node.next[level]
        return node.value

    def insert(self, value):
        # find first node on each level where node.next[levels].value > value
        chain = [None] * self.maxlevels
        steps_at_level = [0] * self.maxlevels
        node = self.head
        for level in reversed(range(self.maxlevels)):
            while node.next[level].value <= value:
                steps_at_level[level] += node.width[level]
                node = node.next[level]
            chain[level] = node

        # insert a link to the newnode at each level
        d = min(self.maxlevels, 1 - int(log(random(), 2.0)))
        newnode = Node(value, [None] * d, [None] * d)
        steps = 0
        for level in range(d):
            prevnode = chain[level]
            newnode.next[level] = prevnode.next[level]
            prevnode.next[level] = newnode
            newnode.width[level] = prevnode.width[level] - steps
            prevnode.width[level] = steps + 1
            steps += steps_at_level[level]
        for level in range(d, self.maxlevels):
            chain[level].width[level] += 1
        self.size += 1

    def remove(self, value):
        # find first node on each level where node.next[levels].value >= value
        chain = [None] * self.maxlevels
        node = self.head
        for level in reversed(range(self.maxlevels)):
            while node.next[level].value < value:
                node = node.next[level]
            chain[level] = node
        if value != chain[0].next[0].value:
            raise KeyError('Not Found')

        # remove one link at each level
        d = len(chain[0].next[0].next)
        for level in range(d):
            prevnode = chain[level]
            prevnode.width[level] += prevnode.next[level].width[level] - 1
            prevnode.next[level] = prevnode.next[level].next[level]
        for level in range(d, self.maxlevels):
            chain[level].width[level] -= 1
        self.size -= 1

    def __iter__(self):
        'Iterate over values in sorted order'
        node = self.head.next[0]
        while node is not NIL:
            yield node.value
            node = node.next[0]

from collections import deque
from itertools import islice


class RunningMedian:
    'Fast running median with O(lg n) updates where n is the window size'

    def __init__(self, n, iterable):
        from pandas.lib.skiplist import IndexableSkiplist as skiplist

        self.it = iter(iterable)
        self.queue = deque(islice(self.it, n))
        self.skiplist = IndexableSkiplist(n)
        for elem in self.queue:
            self.skiplist.insert(elem)

    def __iter__(self):
        queue = self.queue
        skiplist = self.skiplist
        midpoint = len(queue) // 2
        yield skiplist[midpoint]
        for newelem in self.it:
            oldelem = queue.popleft()
            skiplist.remove(oldelem)
            queue.append(newelem)
            skiplist.insert(newelem)
            yield skiplist[midpoint]

N = 100000
K = 10000

import time


def test():
    from numpy.random import randn

    arr = randn(N)

    def _test(arr, k):
        meds = RunningMedian(k, arr)
        return list(meds)

    _test(arr, K)



def test2():

    arr = randn(N)

    return rolling_median(arr, K)


def runmany(f, arr, arglist):
    timings = []

    for arg in arglist:
        tot = 0
        for i in range(5):
            tot += _time(f, arr, arg)
        timings.append(tot / 5)

    return timings


def _time(f, *args):
    _start = time.clock()
    result = f(*args)
    return time.clock() - _start

if __name__ == '__main__':
    test2()
