#!/usr/bin/env python
"""
Analyze documentation to detect errors.

Requires the documentation to be build before usage.

If no argument is provided, it validates all pages in the source folder for
which a doctree has been generated.

If a page is provided like "index", "whatsnew/v1.0.0" a list of all errors
in that page is printed.

Usage::
    $ ./validate_documentation.py
    $ ./validate_documentation.py index
"""

import argparse
import docutils.nodes
from fnmatch import fnmatch
import json
import os
import pickle
import re
import sys

BASE_DIR = os.path.join(os.path.dirname(__file__), '..', 'doc')
DOCUMENTATION_SOURCE = os.path.join(BASE_DIR, 'build', 'doctrees', '')
DOCTREE_PATH = os.path.join(DOCUMENTATION_SOURCE, '{}.doctree')
RST_PATH = os.path.join(BASE_DIR, 'source', '{}.rst')

ERROR_MSGS = {
    'DT01': "Found bullet list within block quote. Use 0 spaces for top-level "
            "list and 2 spaces for sub-lists",
    'WS01': "Indentation uses tabulation",
    'WS02': "Trailing whitespaces",
    'WS03': "Whitespace in empty line",
}

STARTING_WHITESPACE_RE = re.compile(r'^(\s+).', re.MULTILINE)
TRAILING_WHITESPACE_RE = re.compile(r'.([\t ]+)$', re.MULTILINE)
EMPTY_LINE_WHITESPACE_RE = re.compile(r'^([\t ]+)$', re.MULTILINE)


class DocumentChecker(object):
    """
    Checker to validate one page in the documentation

    Attributes
    ----------
    page : str
        Path to page relative to documentation's source folder without file
        extension. (e.g. io)
    doctree : docutils.nodes.document
        Generated doctree for associated the page
    raw_lines : list of str
        Lines from rst-file for associated page
    raw_doc : str
        Joined lines from rst-file for associated page

    Notes
    -----
    Add a method starting with `check` to add additional checks.

    """

    def __init__(self, page, raw_lines, doctree):
        self.page = page
        self.doctree = doctree
        self.raw_lines = raw_lines
        self.raw_doc = ''.join(raw_lines)
        self.errs = None

    def error(self, code, line=None, **kwargs):
        """
        Parameters
        ----------
        code : str
            Error code.
        line : Tuple[int, str]
        **kwargs
            Values for the variables in the error messages
        """
        errs = self.errs.setdefault(code, [])
        errs.append((line, kwargs))

    def find_line(self, match):
        """
        Find rows in documentation that were matched

        Parameters
        ----------
        match : typing.Match

        Returns
        -------
        row_start : int
        row_end : int
        """
        if not match:
            return None

        row_start = self.raw_doc[:match.start(0)].count('\n')
        row_end = self.raw_doc[:match.end(0)].count('\n')
        return row_start, row_end + 1

    def check_bullet_list_in_block_quote(self):
        for node in self.doctree.traverse(docutils.nodes.block_quote):
            match = node.first_child_matching_class(docutils.nodes.bullet_list)
            if match is not None:
                self.error('DT01')

    def check_tabulator_as_indentation(self):
        matches = STARTING_WHITESPACE_RE.finditer(self.raw_doc)
        for match in matches:
            if '\t' in match.group(1):
                self.error('WS01', line=self.find_line(match))

    def check_line_ends_with_whitespace(self):
        matches = TRAILING_WHITESPACE_RE.finditer(self.raw_doc)
        for match in matches:
            self.error('WS02', line=self.find_line(match))

    def check_empty_line_contains_whitespace(self):
        matches = EMPTY_LINE_WHITESPACE_RE.finditer(self.raw_doc)
        for match in matches:
            self.error('WS03', line=self.find_line(match))

    def validate(self):
        """Execute methods starting with 'check'"""
        self.errs = {}
        for func in dir(self):
            if func.startswith('check'):
                self.__class__.__dict__[func](self)

        return self.errs

    def report(self, errors=None, output_format='default'):
        """
        Output errors to stdout

        Parameters
        ----------
        errors : list of str, optional
            If provided, filter output by these error codes.
        output_format : str, optional
            One of 'default', 'json', 'azure'

        Returns
        -------
        int
            A integer with number of found issues
        """
        n_errs = 0
        if output_format == 'json':
            output = json.dumps(self.errs)
        else:
            if output_format == 'default':
                output_format = '{path}:{row}:: {code} {text}\n'
            elif output_format == 'azure':
                output_format = ('##vso[task.logissue type=error;'
                                 'sourcepath={path};'
                                 'linenumber={row};'
                                 'code={code};'
                                 ']{text}\n')
            else:
                raise ValueError('Unknown output_format "{}"'.format(
                    output_format))

            output = ''

            for err_code, errs in self.errs.items():
                # The script would be faster if instead of filtering the
                # errors after validating them, it didn't validate them
                # initially. But that would complicate the code too much
                if errors and err_code not in errors:
                    continue
                for line, kwargs in errs:
                    n_errs += 1
                    row_start, row_end = line if line else (0, 0)
                    output += output_format.format(
                        name=self.page,
                        path='doc/source/{}.rst'.format(self.page),
                        row=row_start + 1 if line else '?',
                        code=err_code,
                        source=''.join(self.raw_lines[row_start:row_end]),
                        text=ERROR_MSGS[err_code].format(kwargs))

        sys.stdout.write(output)
        return n_errs


def validate_one(page):
    """
    Validate the page for the given page

    Parameters
    ----------
    page : str
        Path to page relative to documentation's source folder without file
        extension. (e.g. io)

    Returns
    -------
    dict
        A dictionary containing all the information obtained from
        validating the page.

    Notes
    -----
    The errors codes are defined as:
    - First two characters: Type of errors:
       * DT: Error with unwanted node constellations inside the doctree
       * WS: Issues with whitespace characters
    - Last two characters: Numeric error code
    """
    try:
        with open(DOCTREE_PATH.format(page), 'r+b') as file:
            doctree = pickle.load(file)

        with open(RST_PATH.format(page), 'r') as file:
            raw_doc = file.readlines()
    except FileNotFoundError:
        return None

    checker = DocumentChecker(page, raw_doc, doctree)
    checker.validate()
    return checker


def validate_all(exclude_patterns):
    """
    Execute the validation of all pages, and return a dict with the
    results.

    Parameters
    ----------
    exclude_patterns : List[str] or None
        If provided, the pages that match with one of these patterns
        will be ignored. If None, all pages will be validated.

    Returns
    -------
    dict
        A dictionary with an item for every page containing
        all the validation information.
    """

    checkers = {}
    for root, dirs, files in os.walk(DOCUMENTATION_SOURCE):
        _, base_dir = root.split(DOCUMENTATION_SOURCE)
        for file in files:
            docname, ext = os.path.splitext(file)
            if not ext == '.doctree':
                continue
            page = os.path.join(base_dir, docname)
            if exclude_patterns:
                for pattern in exclude_patterns:
                    if fnmatch(page, pattern):
                        continue

            checker = validate_one(page)
            if checker:
                checkers[page] = checker
    return checkers


def main(page, errors, output_format, exclude_patterns=None):
    if page:
        checkers = {page: validate_one(page)}
    else:
        checkers = validate_all(exclude_patterns=exclude_patterns)
    exit_code = 0
    for page, checker in checkers.items():
        exit_code += checker.report(errors=errors, output_format=output_format)
    return exit_code


if __name__ == '__main__':
    format_opts = 'default', 'json', 'azure'
    argparser = argparse.ArgumentParser(
        description='validate pandas documentation')
    add = argparser.add_argument
    add('page', nargs='?', default=None,
        help='page to validate (e.g. io) '
             'if not provided, all pages are validated')
    add('--format', default='default', choices=format_opts,
        help='format of the output when validating '
             'multiple documents (ignored when validating one).'
             'It can be {}'.format(str(format_opts)[1:-1]))
    add('--errors', default=None,
        help='comma separated '
             'list of error codes to validate. By default it '
             'validates all errors (ignored when validating '
             'a single document)')
    add('--exclude', default=None,
        help='comma separated '
             'patterns of pages to exclude. Utilises '
             '`Unix filename pattern matching`'
             '(ignored when validating a single document)')

    args = argparser.parse_args()
    sys.exit(main(args.page,
                  args.errors.split(',') if args.errors else None,
                  args.format,
                  args.exclude.split(',') if args.exclude else None))
