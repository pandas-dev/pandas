#!/usr/bin/env python
"""
Python script for building documentation.

To build the docs you must have all optional dependencies for pandas
installed. See the installation instructions for a list of these.

Usage
-----
    $ python make.py clean
    $ python make.py html
    $ python make.py latex
"""
import importlib
import sys
import os
import shutil
# import subprocess
import argparse
from contextlib import contextmanager
import webbrowser


DOC_PATH = os.path.dirname(os.path.abspath(__file__))
SOURCE_PATH = os.path.join(DOC_PATH, 'source')
BUILD_PATH = os.path.join(DOC_PATH, 'build')
BUILD_DIRS = ['doctrees', 'html', 'latex', 'plots', '_static', '_templates']


class DocBuilder:
    """
    Class to wrap the different commands of this script.

    All public methods of this class can be called as parameters of the
    script.
    """
    def __init__(self, num_jobs=1, include_api=True, single_doc=None,
                 verbosity=0, warnings_are_errors=False):
        self.num_jobs = num_jobs
        self.verbosity = verbosity
        self.warnings_are_errors = warnings_are_errors

        if single_doc:
            single_doc = self._process_single_doc(single_doc)
            include_api = False
            os.environ['SPHINX_PATTERN'] = single_doc
        elif not include_api:
            os.environ['SPHINX_PATTERN'] = '-api'

        # if self.single_doc_type == 'docstring':
        #    self._run_os('sphinx-autogen', os.path.join(SOURCE_PATH, 'index.rst'))

        if single_doc:
            self.single_doc_html = os.path.splitext(single_doc)[0] + '.html'


    def _process_single_doc(self, single_doc):
        """
        Make sure the provided value for --single is a path to an existing
        .rst/.ipynb file, or a pandas object that can be imported.

        For example, categorial.rst or pandas.DataFrame.head. For the latter,
        return the corresponding file path
        (e.g. generated/pandas.DataFrame.head.rst).
        """
        base_name, extension = os.path.splitext(single_doc)
        if extension in ('.rst', '.ipynb'):
            if os.path.exists(os.path.join(SOURCE_PATH, single_doc)):
                return single_doc
            else:
                raise FileNotFoundError('File {} not found'.format(single_doc))

        elif single_doc.startswith('pandas.'):
            try:
                obj = pandas  # noqa: F821
                for name in single_doc.split('.'):
                    obj = getattr(obj, name)
            except AttributeError:
                raise ImportError('Could not import {}'.format(single_doc))
            else:
                return os.path.join('generated', '{}.rst'.format(single_doc))
        else:
            raise ValueError('--single value should be a valid path to a '
                             '.rst or .ipynb file, or a valid pandas object '
                             '(e.g. categorical.rst or pandas.DataFrame.head)')

    @staticmethod
    def _run_os(*args):
        """
        Execute a command as a OS terminal.

        Parameters
        ----------
        *args : list of str
            Command and parameters to be executed

        Examples
        --------
        >>> DocBuilder()._run_os('python', '--version')
        """
        # TODO check_call should be more safe, but it fails with
        # exclude patterns, needs investigation
        # subprocess.check_call(args, stderr=subprocess.STDOUT)
        exit_status = os.system(' '.join(args))
        if exit_status:
            msg = 'Command "{}" finished with exit code {}'
            raise RuntimeError(msg.format(' '.join(args), exit_status))
        #subprocess.check_call(args, stderr=subprocess.STDOUT)

    def _sphinx_build(self, kind):
        """
        Call sphinx to build documentation.

        Attribute `num_jobs` from the class is used.

        Parameters
        ----------
        kind : {'html', 'latex'}

        Examples
        --------
        >>> DocBuilder(num_jobs=4)._sphinx_build('html')
        """
        if kind not in ('html', 'latex'):
            raise ValueError('kind must be html or latex, '
                             'not {}'.format(kind))

        self._run_os('sphinx-build',
                     '-j{}'.format(self.num_jobs),
                     '-b{}'.format(kind),
                     '-W' if self.warnings_are_errors else '',
                     '-{}'.format(
                         'v' * self.verbosity) if self.verbosity else '',
                     '-d"{}"'.format(os.path.join(BUILD_PATH, 'doctrees')),
                     '"{}"'.format(SOURCE_PATH),
                     '"{}"'.format(os.path.join(BUILD_PATH, kind)))

    def _open_browser(self, single_doc_html):
        """
        Open a browser tab showing single
        """
        url = os.path.join('file://', DOC_PATH, 'build', 'html', single_doc_html)
        webbrowser.open(url, new=2)

    def html(self):
        """
        Build HTML documentation.
        """
        self._sphinx_build('html')
        zip_fname = os.path.join(BUILD_PATH, 'html', 'pandas.zip')
        if os.path.exists(zip_fname):
            os.remove(zip_fname)

        if self.single_doc_html is not None:
            self._open_browser(self.single_doc_html)

    def latex(self, force=False):
        """
        Build PDF documentation.
        """
        if sys.platform == 'win32':
            sys.stderr.write('latex build has not been tested on windows\n')
        else:
            self._sphinx_build('latex')
            os.chdir(os.path.join(BUILD_PATH, 'latex'))
            if force:
                for i in range(3):
                    self._run_os('pdflatex',
                                 '-interaction=nonstopmode',
                                 'pandas.tex')
                raise SystemExit('You should check the file '
                                 '"build/latex/pandas.pdf" for problems.')
            else:
                self._run_os('make')

    def latex_forced(self):
        """
        Build PDF documentation with retries to find missing references.
        """
        self.latex(force=True)

    @staticmethod
    def clean():
        """
        Clean documentation generated files.
        """
        shutil.rmtree(BUILD_PATH, ignore_errors=True)
        shutil.rmtree(os.path.join(SOURCE_PATH, 'generated'),
                      ignore_errors=True)

    def zip_html(self):
        """
        Compress HTML documentation into a zip file.
        """
        zip_fname = os.path.join(BUILD_PATH, 'html', 'pandas.zip')
        if os.path.exists(zip_fname):
            os.remove(zip_fname)
        dirname = os.path.join(BUILD_PATH, 'html')
        fnames = os.listdir(dirname)
        os.chdir(dirname)
        self._run_os('zip',
                     zip_fname,
                     '-r',
                     '-q',
                     *fnames)


def main():
    cmds = [method for method in dir(DocBuilder) if not method.startswith('_')]

    argparser = argparse.ArgumentParser(
        description='pandas documentation builder',
        epilog='Commands: {}'.format(','.join(cmds)))
    argparser.add_argument('command',
                           nargs='?',
                           default='html',
                           help='command to run: {}'.format(', '.join(cmds)))
    argparser.add_argument('--num-jobs',
                           type=int,
                           default=1,
                           help='number of jobs used by sphinx-build')
    argparser.add_argument('--no-api',
                           default=False,
                           help='ommit api and autosummary',
                           action='store_true')
    argparser.add_argument('--single',
                           metavar='FILENAME',
                           type=str,
                           default=None,
                           help=('filename of section or method name to '
                                 'compile, e.g. "indexing", "DataFrame.join"'))
    argparser.add_argument('--python-path',
                           type=str,
                           default=os.path.dirname(DOC_PATH),
                           help='path')
    argparser.add_argument('-v', action='count', dest='verbosity', default=0,
                           help=('increase verbosity (can be repeated), '
                                 'passed to the sphinx build command'))
    argparser.add_argument('--warnings-are-errors', '-W',
                           action='store_true',
                           help='fail if warnings are raised')
    args = argparser.parse_args()

    if args.command not in cmds:
        raise ValueError('Unknown command {}. Available options: {}'.format(
            args.command, ', '.join(cmds)))

    # Below we update both os.environ and sys.path. The former is used by
    # external libraries (namely Sphinx) to compile this module and resolve
    # the import of `python_path` correctly. The latter is used to resolve
    # the import within the module, injecting it into the global namespace
    os.environ['PYTHONPATH'] = args.python_path
    sys.path.append(args.python_path)
    globals()['pandas'] = importlib.import_module('pandas')

    # Set the matplotlib backend to the non-interactive Agg backend for all
    # child processes.
    os.environ['MPLBACKEND'] = 'module://matplotlib.backends.backend_agg'

    builder = DocBuilder(args.num_jobs, not args.no_api, args.single,
                         args.verbosity, args.warnings_are_errors)
    getattr(builder, args.command)()


if __name__ == '__main__':
    sys.exit(main())
