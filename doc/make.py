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
import re
import textwrap
import shutil
# import subprocess
import argparse
from collections import namedtuple
from contextlib import contextmanager
import webbrowser
import jinja2


DOC_PATH = os.path.dirname(os.path.abspath(__file__))
SOURCE_PATH = os.path.join(DOC_PATH, 'source')
BUILD_PATH = os.path.join(DOC_PATH, 'build')
BUILD_DIRS = ['doctrees', 'html', 'latex', 'plots', '_static', '_templates']


@contextmanager
def _maybe_exclude_notebooks():
    """Skip building the notebooks if pandoc is not installed.

    This assumes that nbsphinx is installed.

    Skip notebook conversion if:
    1. nbconvert isn't installed, or
    2. nbconvert is installed, but pandoc isn't
    """
    # TODO move to exclude_pattern
    base = os.path.dirname(__file__)
    notebooks = [os.path.join(base, 'source', nb)
                 for nb in ['style.ipynb']]
    contents = {}

    def _remove_notebooks():
        for nb in notebooks:
            with open(nb, 'rt') as f:
                contents[nb] = f.read()
            os.remove(nb)

    try:
        import nbconvert
    except ImportError:
        sys.stderr.write('Warning: nbconvert not installed. '
                         'Skipping notebooks.\n')
        _remove_notebooks()
    else:
        try:
            nbconvert.utils.pandoc.get_pandoc_version()
        except nbconvert.utils.pandoc.PandocMissing:
            sys.stderr.write('Warning: Pandoc is not installed. '
                             'Skipping notebooks.\n')
            _remove_notebooks()

    yield

    for nb, content in contents.items():
        with open(nb, 'wt') as f:
            f.write(content)


class DocBuilder:
    """Class to wrap the different commands of this script.

    All public methods of this class can be called as parameters of the
    script.
    """
    def __init__(self, num_jobs=1, include_api=True, single_doc=None,
                 verbosity=0, warnings_are_errors=False, log_file=None):
        self.num_jobs = num_jobs
        self.include_api = include_api
        self.verbosity = verbosity
        self.single_doc = None
        self.single_doc_type = None
        if single_doc is not None:
            self._process_single_doc(single_doc)
        self.exclude_patterns = self._exclude_patterns
        self.warnings_are_errors = warnings_are_errors
        self.log_file = log_file

        self._generate_index()
        if self.single_doc_type == 'docstring':
            self._run_os('sphinx-autogen', '-o',
                         'source/generated_single', 'source/index.rst')

    @property
    def _exclude_patterns(self):
        """Docs source files that will be excluded from building."""
        # TODO move maybe_exclude_notebooks here
        if self.single_doc is not None:
            rst_files = [f for f in os.listdir(SOURCE_PATH)
                         if ((f.endswith('.rst') or f.endswith('.ipynb'))
                             and (f != 'index.rst')
                             and (f != '{0}.rst'.format(self.single_doc)))]
            if self.single_doc_type != 'api':
                rst_files += ['generated/*.rst']
        elif not self.include_api:
            rst_files = ['api.rst', 'generated/*.rst']
        else:
            rst_files = ['generated_single/*.rst']

        exclude_patterns = ','.join(
            '{!r}'.format(i) for i in ['**.ipynb_checkpoints'] + rst_files)

        return exclude_patterns

    def _process_single_doc(self, single_doc):
        """Extract self.single_doc (base name) and self.single_doc_type from
        passed single_doc kwarg.

        """
        self.include_api = False

        if single_doc == 'api.rst' or single_doc == 'api':
            self.single_doc_type = 'api'
            self.single_doc = 'api'
        elif os.path.exists(os.path.join(SOURCE_PATH, single_doc)):
            self.single_doc_type = 'rst'
            self.single_doc = os.path.splitext(os.path.basename(single_doc))[0]
        elif os.path.exists(
                os.path.join(SOURCE_PATH, '{}.rst'.format(single_doc))):
            self.single_doc_type = 'rst'
            self.single_doc = single_doc
        elif single_doc is not None:
            try:
                obj = pandas  # noqa: F821
                for name in single_doc.split('.'):
                    try:
                        # for names not in the top-level namespace by default,
                        # e.g. pandas.io.formats.style.Styler
                        importlib.import_module('.'.join([obj.__name__, name]))
                    except ModuleNotFoundError:
                        pass
                    obj = getattr(obj, name)
            except AttributeError:
                raise ValueError('Single document not understood, it should '
                                 'be a file in doc/source/*.rst (e.g. '
                                 '"contributing.rst" or a pandas function or '
                                 'method (e.g. "pandas.DataFrame.head")')
            else:
                self.single_doc_type = 'docstring'
                if single_doc.startswith('pandas.'):
                    self.single_doc = single_doc[len('pandas.'):]
                else:
                    self.single_doc = single_doc

    def _copy_generated_docstring(self):
        """Copy existing generated (from api.rst) docstring page because
        this is more correct in certain cases (where a custom autodoc
        template is used).

        """
        fname = os.path.join(SOURCE_PATH, 'generated',
                             'pandas.{}.rst'.format(self.single_doc))
        temp_dir = os.path.join(SOURCE_PATH, 'generated_single')

        try:
            os.makedirs(temp_dir)
        except OSError:
            pass

        if os.path.exists(fname):
            try:
                # copying to make sure sphinx always thinks it is new
                # and needs to be re-generated (to pick source code changes)
                shutil.copy(fname, temp_dir)
            except:  # noqa
                pass

    def _generate_index(self):
        """Create index.rst file with the specified sections."""
        if self.single_doc_type == 'docstring':
            self._copy_generated_docstring()

        with open(os.path.join(SOURCE_PATH, 'index.rst.template')) as f:
            t = jinja2.Template(f.read())

        with open(os.path.join(SOURCE_PATH, 'index.rst'), 'w') as f:
            f.write(t.render(include_api=self.include_api,
                             single_doc=self.single_doc,
                             single_doc_type=self.single_doc_type))

    @staticmethod
    def _create_build_structure():
        """Create directories required to build documentation."""
        for dirname in BUILD_DIRS:
            try:
                os.makedirs(os.path.join(BUILD_PATH, dirname))
            except OSError:
                pass

    @staticmethod
    def _run_os(*args):
        """Execute a command as a OS terminal.

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
        os.system(' '.join(args))

    def _sphinx_build(self, kind):
        """Call sphinx to build documentation.

        Attribute `num_jobs` from the class is used.

        Parameters
        ----------
        kind : {'html', 'latex'}

        Examples
        --------
        >>> DocBuilder(num_jobs=4)._sphinx_build('html')
        """
        if kind not in ('html', 'latex', 'spelling'):
            raise ValueError('kind must be html, latex or '
                             'spelling, not {}'.format(kind))
        self._run_os('sphinx-build',
                     '-j{}'.format(self.num_jobs),
                     '-b{}'.format(kind),
                     '{}'.format("W" if self.warnings_are_errors else ""),
                     '-{}'.format(
                         'v' * self.verbosity) if self.verbosity else '',
                     '-d{}'.format(os.path.join(BUILD_PATH, 'doctrees')),
                     '-Dexclude_patterns={}'.format(self.exclude_patterns),
                     SOURCE_PATH,
                     os.path.join(BUILD_PATH, kind))

    def _open_browser(self):
        base_url = os.path.join('file://', DOC_PATH, 'build', 'html')
        if self.single_doc_type == 'docstring':
            url = os.path.join(
                base_url,
                'generated_single', 'pandas.{}.html'.format(self.single_doc))
        else:
            url = os.path.join(base_url, '{}.html'.format(self.single_doc))
        webbrowser.open(url, new=2)

    def html(self):
        """Build HTML documentation."""
        self._create_build_structure()
        with _maybe_exclude_notebooks():
            self._sphinx_build('html')
            zip_fname = os.path.join(BUILD_PATH, 'html', 'pandas.zip')
            if os.path.exists(zip_fname):
                os.remove(zip_fname)

        if self.single_doc is not None:
            self._open_browser()
            shutil.rmtree(os.path.join(SOURCE_PATH, 'generated_single'),
                          ignore_errors=True)

    def latex(self, force=False):
        """Build PDF documentation."""
        self._create_build_structure()
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
        """Build PDF documentation with retries to find missing references."""
        self.latex(force=True)

    @staticmethod
    def clean():
        """Clean documentation generated files."""
        shutil.rmtree(BUILD_PATH, ignore_errors=True)
        shutil.rmtree(os.path.join(SOURCE_PATH, 'generated'),
                      ignore_errors=True)

    def zip_html(self):
        """Compress HTML documentation into a zip file."""
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

    def spellcheck(self):
        """Spell check the documentation."""
        self._sphinx_build('spelling')
        output_location = os.path.join('build', 'spelling', 'output.txt')
        with open(output_location) as output:
            lines = output.readlines()
            if lines:
                raise SyntaxError(
                    'Found misspelled words.'
                    ' Check pandas/doc/build/spelling/output.txt'
                    ' for more details.')

    def lint_log(self):
        with open(self.log_file) as f:
            log = f.read()

        tokens = tokenize_log(log)
        failed = [tok for tok in tokens if tok.kind != 'OK']
        if failed:
            report_failures(failed)
            sys.exit(1)

# ------
# Linter
# ------

LinterToken = namedtuple("Token", ['kind', 'value'])
IPY_ERROR = r'(?P<IPY_ERROR>>>>-*\n.*?<<<-*\n)'
SPHINX_WARNING = r'(?P<SPHINX_WARNING>^[^\n]*?: WARNING:.*?$\n?)'
OK = r'(?P<OK>^.*?\n)'


def tokenize_log(log):
    master_pat = re.compile("|".join([IPY_ERROR, SPHINX_WARNING, OK]),
                            flags=re.MULTILINE | re.DOTALL)

    def generate_tokens(pat, text):
        scanner = pat.scanner(text)
        for m in iter(scanner.match, None):
            yield LinterToken(m.lastgroup, m.group(m.lastgroup))

    tok = list(generate_tokens(master_pat, log))
    return tok


def report_failures(failed):
    tpl = textwrap.dedent("""\
    {n} failure{s}

    {individual}
    """)
    joined = []
    for i, tok in enumerate(failed):
        line = "Failure [{}]: {}".format(i, tok.value.strip())
        joined.append(line)
    joined = '\n'.join(joined)

    print(tpl.format(n=len(failed),
                     s="s" if len(failed) != 1 else "",
                     individual=joined))


# ---
# CLI
# ---


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
    argparser.add_argument("--warnings-are-errors",
                           default=False,
                           action="store_true",
                           help="Whether to fail the build on warnings.")
    argparser.add_argument("--log-file",
                           default="doc-build.log",
                           help="Log file of the build to lint for warnings.")
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
                         args.verbosity, args.warnings_are_errors,
                         args.log_file)
    getattr(builder, args.command)()


if __name__ == '__main__':
    sys.exit(main())
