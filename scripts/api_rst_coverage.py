#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
Script to generate a report with the coverage of the API in the docs.

The output of this script shows the existing methods that are not
included in the API documentation, as well as the methods documented
that do not exist. Ideally, no method should be listed. Currently it
considers the methods of Series, DataFrame and Panel.

Deprecated methods are usually removed from the documentation, while
still available for three minor versions. They are listed with the
word deprecated and the version number next to them.

Usage::

    $ PYTHONPATH=.. ./api_rst_coverage.py

"""
import os
import re
import inspect
import pandas as pd


def main():
    # classes whose members to check
    classes = [pd.Series, pd.DataFrame, pd.Panel]

    def class_name_sort_key(x):
        if x.startswith('Series'):
            # make sure Series precedes DataFrame, and Panel.
            return ' ' + x
        else:
            return x

    def get_docstring(x):
        class_name, method = x.split('.')
        obj = getattr(getattr(pd, class_name), method)
        return obj.__doc__

    def deprecation_version(x):
        pattern = re.compile('\.\. deprecated:: ([0-9]+\.[0-9]+\.[0-9]+)')
        doc = get_docstring(x)
        match = pattern.search(doc)
        if match:
            return match.groups()[0]

    def add_notes(x):
        # Some methods are not documented in api.rst because they
        # have been deprecated. Adding a comment to detect them easier.
        doc = get_docstring(x)
        note = None
        if not doc:
            note = 'no docstring'
        else:
            version = deprecation_version(x)
            if version:
                note = 'deprecated in {}'.format(version)

        return '{} ({})'.format(x, note) if note else x

    # class members
    class_members = set()
    for cls in classes:
        for member in inspect.getmembers(cls):
            class_members.add('{cls}.{member}'.format(cls=cls.__name__,
                                                      member=member[0]))

    # class members referenced in api.rst
    api_rst_members = set()
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    api_rst_fname = os.path.join(base_path, 'doc', 'source', 'api.rst')
    class_names = (cls.__name__ for cls in classes)
    pattern = re.compile('({})\.(\w+)'.format('|'.join(class_names)))
    with open(api_rst_fname, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                api_rst_members.add(match.group(0))

    print()
    print("Documented members in api.rst that aren't actual class members:")
    for x in sorted(api_rst_members.difference(class_members),
                    key=class_name_sort_key):
        print(x)

    print()
    print("Class members (other than those beginning with '_') "
          "missing from api.rst:")
    for x in sorted(class_members.difference(api_rst_members),
                    key=class_name_sort_key):
        if '._' not in x:
            print(add_notes(x))


if __name__ == "__main__":
    main()
