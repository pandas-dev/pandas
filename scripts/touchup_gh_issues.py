#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from collections import OrderedDict
import sys
import re

"""
Reads in stdin, replace all occurences of '#num' or 'GH #num' with
links to github issue. dumps the issue anchors before the next
section header
"""

pat = "((?:\s*GH\s*)?)#(\d{3,4})([^_]|$)?"
rep_pat = r"\1GH\2_\3"
anchor_pat = ".. _GH{id}: https://github.com/pandas-dev/pandas/issues/{id}"
section_pat = "^pandas\s[\d\.]+\s*$"


def main():
    issues = OrderedDict()
    while True:

        line = sys.stdin.readline()
        if not line:
            break

        if re.search(section_pat, line):
            for id in issues:
                print(anchor_pat.format(id=id).rstrip())
            if issues:
                print("\n")
            issues = OrderedDict()

        for m in re.finditer(pat, line):
            id = m.group(2)
            if id not in issues:
                issues[id] = True
        print(re.sub(pat, rep_pat, line).rstrip())
    pass

if __name__ == "__main__":
    main()
