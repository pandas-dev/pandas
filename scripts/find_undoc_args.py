#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from collections import namedtuple
from itertools import islice
import types
import os
import re
import argparse
#http://docs.python.org/2/library/argparse.html
# arg name is positional is not prefixed with - or --

parser = argparse.ArgumentParser(description='Program description.')
parser.add_argument('-p', '--path', metavar='PATH', type=str,required=False,
                    default=None,
                   help='full path relative to which paths wills be reported',action='store')
parser.add_argument('-m', '--module', metavar='MODULE', type=str,required=True,
                   help='name of package to import and examine',action='store')

args = parser.parse_args()

Entry=namedtuple("Entry","func loc undoc_names missing_args nsig_names ndoc_names")

def entry_gen(root_ns,module_name):

    q=[root_ns]
    seen=set()
    while q:
        ns = q.pop()
        for x in dir(ns):
            cand = getattr(ns,x)
            if (isinstance(cand,types.ModuleType)
                and cand.__name__ not in seen
                and cand.__name__.startswith(module_name)):
                # print(cand.__name__)
                seen.add(cand.__name__)
                q.insert(0,cand)
            elif (isinstance(cand,(types.MethodType,types.FunctionType)) and
                  cand not in seen and cand.func_doc):
                seen.add(cand)
                yield cand

def cmp_docstring_sig(f):
    def build_loc(f):
        path=f.func_code.co_filename.split(args.path,1)[-1][1:]
        return "+{} {}".format(f.func_code.co_firstlineno,path)
    import inspect
    sig_names=set(inspect.getargspec(f).args)
    doc = f.func_doc.lower()
    doc = re.split("^\s*parameters\s*",doc,1,re.M)[-1]
    doc = re.split("^\s*returns*",doc,1,re.M)[0]
    doc_names={x.split(":")[0].strip() for x in doc.split("\n")
                if re.match("\s+[\w_]+\s*:",x)}
    sig_names.discard("self")
    doc_names.discard("kwds")
    doc_names.discard("kwargs")
    doc_names.discard("args")
    return Entry(func=f,loc=build_loc(f),undoc_names=sig_names.difference(doc_names),
                 missing_args=doc_names.difference(sig_names),nsig_names=len(sig_names),
                 ndoc_names=len(doc_names))

def main():
    module = __import__(args.module)
    if not args.path:
        args.path=os.path.dirname(module.__file__)
    collect=[cmp_docstring_sig(e) for e in entry_gen(module,module.__name__)]
    # only include if there are missing arguments in the docstring (less false positives)
    # and there are at least some documented arguments
    collect = [e for e in collect if e.undoc_names and len(e.undoc_names) != e.nsig_names]

    tmpl = "{}:[{}]\t missing[{}/{}]={}"
    for x in collect:
        s=  tmpl.format(x.loc,x.func.__name__,len(x.undoc_names),
                        x.nsig_names,list(x.undoc_names))
        if x.missing_args:
            s+= " extra(?)={}".format(list(x.missing_args))
        print(s)

if __name__ == "__main__":
    import sys
    sys.exit(main())
