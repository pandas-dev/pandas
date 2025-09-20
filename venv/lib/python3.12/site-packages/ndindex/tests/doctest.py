"""
Custom script to run the doctests

This runs the doctests but ignores trailing ``` in Markdown documents.

Running this separately from pytest also allows us to not include the doctests
in the coverage. It also allows us to force a separate namespace for each
docstring's doctest, which the pytest doctest integration does not allow.

TODO: Make these tests also run with pytest, but still keeping them out of the
coverage.

WARNING: This file registers these flags globally in the standard library
doctest module and monkeypatches the DoctestRunner object. It should only be
run as a standalone script. Do not attempt to import it.

"""

import numpy

import sys
import unittest
import glob
import os
from contextlib import contextmanager
from doctest import (DocTestRunner, DocFileSuite, DocTestSuite,
                     NORMALIZE_WHITESPACE, ELLIPSIS, register_optionflag)
import doctest

SKIPNP1 = register_optionflag("SKIPNP1")
NP1 = numpy.__version__.startswith('1')
if NP1:
    SKIP_THIS_VERSION = SKIPNP1
else:
    SKIP_THIS_VERSION = 0

@contextmanager
def patch_doctest():
    """
    Context manager to patch the doctester

    The doctester must be patched
    """
    orig_run = DocTestRunner.run
    orig_indent = doctest._indent

    def run(self, test, **kwargs):
        filtered_examples = []

        for example in test.examples:
            if SKIP_THIS_VERSION not in example.options:
                filtered_examples.append(example)

            # Remove ```
            example.want = example.want.replace('```\n', '')
            example.exc_msg = example.exc_msg and example.exc_msg.replace('```\n', '')

        test.examples = filtered_examples
        return orig_run(self, test, **kwargs)

    # Doctest indents the output, which is annoying for copy-paste, so disable
    # it.
    def _indent(s, **kwargs):
        return s

    try:
        DocTestRunner.run = run
        doctest._indent = _indent
        yield
    finally:
        DocTestRunner.run = orig_run
        doctest._indent = orig_indent

DOCS = os.path.realpath(os.path.join(__file__, os.path.pardir, os.path.pardir,
                                     os.pardir, 'docs'))
MARKDOWN = glob.glob(os.path.join(DOCS, '**', '*.md'), recursive=True)
RST = glob.glob(os.path.join(DOCS, '**', '*.rst'), recursive=True)
README = os.path.realpath(os.path.join(__file__, os.path.pardir, os.path.pardir,
                                     os.pardir, 'README.md'))
def load_tests(loader, tests, ignore):
    for mod in sys.modules:
        if mod.startswith('ndindex'):
            # globs={} makes the doctests not include module names
            tests.addTests(DocTestSuite(sys.modules[mod], globs={},
                                        optionflags=NORMALIZE_WHITESPACE))
    tests.addTests(DocFileSuite(*MARKDOWN, *RST, README,
                                optionflags=NORMALIZE_WHITESPACE | ELLIPSIS,
                                module_relative=False))
    return tests

def run_doctests():
    numpy.seterr(all='ignore')
    with patch_doctest():
        return unittest.main(module='ndindex.tests.doctest', exit=False).result

if __name__ == '__main__':
    # TODO: Allow specifying which doctests to run at the command line
    run_doctests()
