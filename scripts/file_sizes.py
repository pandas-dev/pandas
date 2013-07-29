from __future__ import print_function
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from pandas import DataFrame
from pandas.util.testing import set_trace
from pandas import compat

dirs = []
names = []
lengths = []

if len(sys.argv) > 1:
    loc = sys.argv[1]
else:
    loc = '.'
walked = os.walk(loc)


def _should_count_file(path):
    return path.endswith('.py') or path.endswith('.pyx')


def _is_def_line(line):
    """def/cdef/cpdef, but not `cdef class`"""
    return (line.endswith(':') and not 'class' in line.split() and
            (line.startswith('def ') or
             line.startswith('cdef ') or
             line.startswith('cpdef ') or
             ' def ' in line or ' cdef ' in line or ' cpdef ' in line))


class LengthCounter(object):
    """
    should add option for subtracting nested function lengths??
    """
    def __init__(self, lines):
        self.lines = lines
        self.pos = 0
        self.counts = []
        self.n = len(lines)

    def get_counts(self):
        self.pos = 0
        self.counts = []
        while self.pos < self.n:
            line = self.lines[self.pos]
            self.pos += 1
            if _is_def_line(line):
                level = _get_indent_level(line)
                self._count_function(indent_level=level)
        return self.counts

    def _count_function(self, indent_level=1):
        indent = '    ' * indent_level

        def _end_of_function(line):
            return (line != '' and
                    not line.startswith(indent) and
                    not line.startswith('#'))

        start_pos = self.pos
        while self.pos < self.n:
            line = self.lines[self.pos]
            if _end_of_function(line):
                self._push_count(start_pos)
                return

            self.pos += 1

            if _is_def_line(line):
                self._count_function(indent_level=indent_level + 1)

        # end of file
        self._push_count(start_pos)

    def _push_count(self, start_pos):
        func_lines = self.lines[start_pos:self.pos]

        if len(func_lines) > 300:
            set_trace()

        # remove blank lines at end
        while len(func_lines) > 0 and func_lines[-1] == '':
            func_lines = func_lines[:-1]

        # remove docstrings and comments
        clean_lines = []
        in_docstring = False
        for line in func_lines:
            line = line.strip()
            if in_docstring and _is_triplequote(line):
                in_docstring = False
                continue

            if line.startswith('#'):
                continue

            if _is_triplequote(line):
                in_docstring = True
                continue

        self.counts.append(len(func_lines))


def _get_indent_level(line):
    level = 0
    while line.startswith('    ' * level):
        level += 1
    return level


def _is_triplequote(line):
    return line.startswith('"""') or line.startswith("'''")


def _get_file_function_lengths(path):
    lines = [x.rstrip() for x in open(path).readlines()]
    counter = LengthCounter(lines)
    return counter.get_counts()

# def test_get_function_lengths():
text = """
class Foo:

def foo():
    def bar():
        a = 1

        b = 2

        c = 3

    foo = 'bar'

def x():
    a = 1

    b = 3

    c = 7

    pass
"""

expected = [5, 8, 7]

lines = [x.rstrip() for x in text.splitlines()]
counter = LengthCounter(lines)
result = counter.get_counts()
assert(result == expected)


def doit():
    for directory, _, files in walked:
        print(directory)
        for path in files:
            if not _should_count_file(path):
                continue

            full_path = os.path.join(directory, path)
            print(full_path)
            lines = len(open(full_path).readlines())

            dirs.append(directory)
            names.append(path)
            lengths.append(lines)

    result = DataFrame({'dirs': dirs, 'names': names,
                        'lengths': lengths})


def doit2():
    counts = {}
    for directory, _, files in walked:
        print(directory)
        for path in files:
            if not _should_count_file(path) or path.startswith('test_'):
                continue

            full_path = os.path.join(directory, path)
            counts[full_path] = _get_file_function_lengths(full_path)

    return counts

counts = doit2()

# counts = _get_file_function_lengths('pandas/tests/test_series.py')

all_counts = []
for k, v in compat.iteritems(counts):
    all_counts.extend(v)
all_counts = np.array(all_counts)

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.hist(all_counts, bins=100)
n = len(all_counts)
nmore = (all_counts > 50).sum()
ax.set_title('%s function lengths, n=%d' % ('pandas', n))
ax.set_ylabel('N functions')
ax.set_xlabel('Function length')
ax.text(100, 300, '%.3f%% with > 50 lines' % ((n - nmore) / float(n)),
        fontsize=18)
plt.show()
