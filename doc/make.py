#!/usr/bin/env python3
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
import argparse
import csv
import importlib
import os
import shutil
import subprocess
import sys
import webbrowser

import docutils
import docutils.parsers.rst

DOC_PATH = os.path.dirname(os.path.abspath(__file__))
SOURCE_PATH = os.path.join(DOC_PATH, "source")
BUILD_PATH = os.path.join(DOC_PATH, "build")
REDIRECTS_FILE = os.path.join(DOC_PATH, "redirects.csv")


class DocBuilder:
    """
    Class to wrap the different commands of this script.

    All public methods of this class can be called as parameters of the
    script.
    """

    def __init__(
        self,
        num_jobs="auto",
        include_api=True,
        whatsnew=False,
        single_doc=None,
        verbosity=0,
        warnings_are_errors=False,
    ):
        self.num_jobs = num_jobs
        self.include_api = include_api
        self.whatsnew = whatsnew
        self.verbosity = verbosity
        self.warnings_are_errors = warnings_are_errors

        if single_doc:
            single_doc = self._process_single_doc(single_doc)
            os.environ["SPHINX_PATTERN"] = single_doc
        elif not include_api:
            os.environ["SPHINX_PATTERN"] = "-api"
        elif whatsnew:
            os.environ["SPHINX_PATTERN"] = "whatsnew"

        self.single_doc_html = None
        if single_doc and single_doc.endswith(".rst"):
            self.single_doc_html = os.path.splitext(single_doc)[0] + ".html"
        elif single_doc:
            self.single_doc_html = f"reference/api/pandas.{single_doc}.html"

    def _process_single_doc(self, single_doc):
        """
        Make sure the provided value for --single is a path to an existing
        .rst/.ipynb file, or a pandas object that can be imported.

        For example, categorial.rst or pandas.DataFrame.head. For the latter,
        return the corresponding file path
        (e.g. reference/api/pandas.DataFrame.head.rst).
        """
        base_name, extension = os.path.splitext(single_doc)
        if extension in (".rst", ".ipynb"):
            if os.path.exists(os.path.join(SOURCE_PATH, single_doc)):
                return single_doc
            else:
                raise FileNotFoundError(f"File {single_doc} not found")

        elif single_doc.startswith("pandas."):
            try:
                obj = pandas  # noqa: F821
                for name in single_doc.split("."):
                    obj = getattr(obj, name)
            except AttributeError as err:
                raise ImportError(f"Could not import {single_doc}") from err
            else:
                return single_doc[len("pandas.") :]
        else:
            raise ValueError(
                f"--single={single_doc} not understood. "
                "Value should be a valid path to a .rst or .ipynb file, "
                "or a valid pandas object "
                "(e.g. categorical.rst or pandas.DataFrame.head)"
            )

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
        subprocess.check_call(args, stdout=sys.stdout, stderr=sys.stderr)

    def _sphinx_build(self, kind: str):
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
        if kind not in ("html", "latex"):
            raise ValueError(f"kind must be html or latex, not {kind}")

        cmd = ["sphinx-build", "-b", kind]
        if self.num_jobs:
            cmd += ["-j", self.num_jobs]
        if self.warnings_are_errors:
            cmd += ["-W", "--keep-going"]
        if self.verbosity:
            cmd.append(f"-{'v' * self.verbosity}")
        cmd += [
            "-d",
            os.path.join(BUILD_PATH, "doctrees"),
            SOURCE_PATH,
            os.path.join(BUILD_PATH, kind),
        ]
        return subprocess.call(cmd)

    def _open_browser(self, single_doc_html):
        """
        Open a browser tab showing single
        """
        url = os.path.join("file://", DOC_PATH, "build", "html", single_doc_html)
        webbrowser.open(url, new=2)

    def _get_page_title(self, page):
        """
        Open the rst file `page` and extract its title.
        """
        fname = os.path.join(SOURCE_PATH, f"{page}.rst")
        option_parser = docutils.frontend.OptionParser(
            components=(docutils.parsers.rst.Parser,)
        )
        doc = docutils.utils.new_document("<doc>", option_parser.get_default_values())
        with open(fname) as f:
            data = f.read()

        parser = docutils.parsers.rst.Parser()
        # do not generate any warning when parsing the rst
        with open(os.devnull, "a") as f:
            doc.reporter.stream = f
            parser.parse(data, doc)

        section = next(
            node for node in doc.children if isinstance(node, docutils.nodes.section)
        )
        title = next(
            node for node in section.children if isinstance(node, docutils.nodes.title)
        )

        return title.astext()

    def _add_redirects(self):
        """
        Create in the build directory an html file with a redirect,
        for every row in REDIRECTS_FILE.
        """
        with open(REDIRECTS_FILE) as mapping_fd:
            reader = csv.reader(mapping_fd)
            for row in reader:
                if not row or row[0].strip().startswith("#"):
                    continue

                html_path = os.path.join(BUILD_PATH, "html")
                path = os.path.join(html_path, *row[0].split("/")) + ".html"

                if not self.include_api and (
                    os.path.join(html_path, "reference") in path
                    or os.path.join(html_path, "generated") in path
                ):
                    continue

                try:
                    title = self._get_page_title(row[1])
                except Exception:
                    # the file can be an ipynb and not an rst, or docutils
                    # may not be able to read the rst because it has some
                    # sphinx specific stuff
                    title = "this page"

                with open(path, "w") as moved_page_fd:
                    html = f"""\
<html>
    <head>
        <meta http-equiv="refresh" content="0;URL={row[1]}.html"/>
    </head>
    <body>
        <p>
            The page has been moved to <a href="{row[1]}.html">{title}</a>
        </p>
    </body>
<html>"""

                    moved_page_fd.write(html)

    def html(self):
        """
        Build HTML documentation.
        """
        ret_code = self._sphinx_build("html")
        zip_fname = os.path.join(BUILD_PATH, "html", "pandas.zip")
        if os.path.exists(zip_fname):
            os.remove(zip_fname)

        if ret_code == 0:
            if self.single_doc_html is not None:
                self._open_browser(self.single_doc_html)
            else:
                self._add_redirects()
                if self.whatsnew:
                    self._open_browser(os.path.join("whatsnew", "index.html"))

        return ret_code

    def latex(self, force=False):
        """
        Build PDF documentation.
        """
        if sys.platform == "win32":
            sys.stderr.write("latex build has not been tested on windows\n")
        else:
            ret_code = self._sphinx_build("latex")
            os.chdir(os.path.join(BUILD_PATH, "latex"))
            if force:
                for i in range(3):
                    self._run_os("pdflatex", "-interaction=nonstopmode", "pandas.tex")
                raise SystemExit(
                    "You should check the file "
                    '"build/latex/pandas.pdf" for problems.'
                )
            else:
                self._run_os("make")
            return ret_code

    def latex_forced(self):
        """
        Build PDF documentation with retries to find missing references.
        """
        return self.latex(force=True)

    @staticmethod
    def clean():
        """
        Clean documentation generated files.
        """
        shutil.rmtree(BUILD_PATH, ignore_errors=True)
        shutil.rmtree(os.path.join(SOURCE_PATH, "reference", "api"), ignore_errors=True)

    def zip_html(self):
        """
        Compress HTML documentation into a zip file.
        """
        zip_fname = os.path.join(BUILD_PATH, "html", "pandas.zip")
        if os.path.exists(zip_fname):
            os.remove(zip_fname)
        dirname = os.path.join(BUILD_PATH, "html")
        fnames = os.listdir(dirname)
        os.chdir(dirname)
        self._run_os("zip", zip_fname, "-r", "-q", *fnames)


def main():
    cmds = [method for method in dir(DocBuilder) if not method.startswith("_")]

    joined = ",".join(cmds)
    argparser = argparse.ArgumentParser(
        description="pandas documentation builder", epilog=f"Commands: {joined}"
    )

    joined = ", ".join(cmds)
    argparser.add_argument(
        "command", nargs="?", default="html", help=f"command to run: {joined}"
    )
    argparser.add_argument(
        "--num-jobs", default="auto", help="number of jobs used by sphinx-build"
    )
    argparser.add_argument(
        "--no-api", default=False, help="omit api and autosummary", action="store_true"
    )
    argparser.add_argument(
        "--whatsnew",
        default=False,
        help="only build whatsnew (and api for links)",
        action="store_true",
    )
    argparser.add_argument(
        "--single",
        metavar="FILENAME",
        type=str,
        default=None,
        help=(
            "filename (relative to the 'source' folder) of section or method name to "
            "compile, e.g. 'development/contributing.rst', "
            "'ecosystem.rst', 'pandas.DataFrame.join'"
        ),
    )
    argparser.add_argument(
        "--python-path", type=str, default=os.path.dirname(DOC_PATH), help="path"
    )
    argparser.add_argument(
        "-v",
        action="count",
        dest="verbosity",
        default=0,
        help=(
            "increase verbosity (can be repeated), "
            "passed to the sphinx build command"
        ),
    )
    argparser.add_argument(
        "--warnings-are-errors",
        "-W",
        action="store_true",
        help="fail if warnings are raised",
    )
    args = argparser.parse_args()

    if args.command not in cmds:
        joined = ", ".join(cmds)
        raise ValueError(f"Unknown command {args.command}. Available options: {joined}")

    # Below we update both os.environ and sys.path. The former is used by
    # external libraries (namely Sphinx) to compile this module and resolve
    # the import of `python_path` correctly. The latter is used to resolve
    # the import within the module, injecting it into the global namespace
    os.environ["PYTHONPATH"] = args.python_path
    sys.path.insert(0, args.python_path)
    globals()["pandas"] = importlib.import_module("pandas")

    # Set the matplotlib backend to the non-interactive Agg backend for all
    # child processes.
    os.environ["MPLBACKEND"] = "module://matplotlib.backends.backend_agg"

    builder = DocBuilder(
        args.num_jobs,
        not args.no_api,
        args.whatsnew,
        args.single,
        args.verbosity,
        args.warnings_are_errors,
    )
    return getattr(builder, args.command)()


if __name__ == "__main__":
    sys.exit(main())
