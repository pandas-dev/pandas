#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script that compares the signature arguments with the ones in the docsting
and returns the differences in plain text or GitHub task list format.

Usage::
    $ ./find_undoc_args.py  (see arguments below)
"""
from __future__ import print_function
import sys
from collections import namedtuple
import types
import os
import re
import argparse
import inspect


parser = argparse.ArgumentParser(description='Program description.')
parser.add_argument('-p', '--path', metavar='PATH', type=str, required=False,
                    default=None, action='store',
                    help='full path relative to which paths wills be reported')
parser.add_argument('-m', '--module', metavar='MODULE', type=str,
                    required=True, action='store',
                    help='name of package to import and examine')
parser.add_argument('-G', '--github_repo', metavar='REPO', type=str,
                    required=False, default=None, action='store',
                    help='github project where the code lives, '
                    'e.g. "pandas-dev/pandas"')
args = parser.parse_args()

Entry = namedtuple('Entry',
                   'func path lnum undoc_names missing_args '
                   'nsig_names ndoc_names')


def entry_gen(root_ns, module_name):
    """Walk and yield all methods and functions in the module root_ns and
    submodules."""
    q = [root_ns]
    seen = set()
    while q:
        ns = q.pop()
        for x in dir(ns):
            cand = getattr(ns, x)
            if (isinstance(cand, types.ModuleType) and
                    cand.__name__ not in seen and
                    cand.__name__.startswith(module_name)):
                seen.add(cand.__name__)
                q.insert(0, cand)
            elif (isinstance(cand, (types.MethodType, types.FunctionType)) and
                  cand not in seen and cand.__doc__):
                seen.add(cand)
                yield cand


def cmp_docstring_sig(f):
    """Return an `Entry` object describing the differences between the
    arguments in the signature and the documented ones."""
    def build_loc(f):
        path = f.__code__.co_filename.split(args.path, 1)[-1][1:]
        return dict(path=path, lnum=f.__code__.co_firstlineno)

    sig_names = set(inspect.getargspec(f).args)
    # XXX numpydoc can be used to get the list of parameters
    doc = f.__doc__.lower()
    doc = re.split('^\s*parameters\s*', doc, 1, re.M)[-1]
    doc = re.split('^\s*returns*', doc, 1, re.M)[0]
    doc_names = {x.split(":")[0].strip() for x in doc.split('\n')
                 if re.match('\s+[\w_]+\s*:', x)}
    sig_names.discard('self')
    doc_names.discard('kwds')
    doc_names.discard('kwargs')
    doc_names.discard('args')
    return Entry(func=f, path=build_loc(f)['path'], lnum=build_loc(f)['lnum'],
                 undoc_names=sig_names.difference(doc_names),
                 missing_args=doc_names.difference(sig_names),
                 nsig_names=len(sig_names), ndoc_names=len(doc_names))


def format_id(i):
    return i


def format_item_as_github_task_list(i, item, repo):
    tmpl = ('- [ ] {id_}) [{fname}:{lnum} ({func_name}())]({link}) -  '
            '__Missing__[{nmissing}/{total_args}]: {undoc_names}')
    link_tmpl = "https://github.com/{repo}/blob/master/{file}#L{lnum}"
    link = link_tmpl.format(repo=repo, file=item.path, lnum=item.lnum)
    s = tmpl.format(id_=i, fname=item.path, lnum=item.lnum,
                    func_name=item.func.__name__, link=link,
                    nmissing=len(item.undoc_names),
                    total_args=item.nsig_names,
                    undoc_names=list(item.undoc_names))
    if item.missing_args:
        s += '    __Extra__(?): %s' % list(item.missing_args)
    return s


def format_item_as_plain(i, item):
    tmpl = ('+{lnum} {path} {func_name}(): '
            'Missing[{nmissing}/{total_args}]={undoc_names}')
    s = tmpl.format(path=item.path, lnum=item.lnum,
                    func_name=item.func.__name__,
                    nmissing=len(item.undoc_names),
                    total_args=item.nsig_names,
                    undoc_names=list(item.undoc_names))
    if item.missing_args:
        s += ' Extra(?)=%s' % list(item.missing_args)
    return s


def main():
    module = __import__(args.module)
    if not args.path:
        args.path = os.path.dirname(module.__file__)
    collect = [cmp_docstring_sig(e)
               for e in entry_gen(module, module.__name__)]
    # only include if there are missing arguments in the docstring
    # (fewer false positives) and there are at least some documented arguments
    collect = [e for e in collect
               if e.undoc_names and len(e.undoc_names) != e.nsig_names]
    collect.sort(key=lambda x: x.path)

    if args.github_repo:
        for i, item in enumerate(collect, 1):
            print(format_item_as_github_task_list(i, item, args.github_repo))
    else:
        for i, item in enumerate(collect, 1):
            print(format_item_as_plain(i, item))


if __name__ == '__main__':
    sys.exit(main())
