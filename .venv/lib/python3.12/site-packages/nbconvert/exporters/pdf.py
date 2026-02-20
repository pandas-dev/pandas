"""Export to PDF via latex"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from traitlets import Bool, Instance, Integer, List, Unicode, default

from nbconvert.utils import _contextlib_chdir

from .latex import LatexExporter


class LatexFailed(IOError):
    """Exception for failed latex run

    Captured latex output is in error.output.
    """

    def __init__(self, output):
        """Initialize the error."""
        self.output = output

    def __unicode__(self):
        """Unicode representation."""
        return "PDF creating failed, captured latex output:\n%s" % self.output

    def __str__(self):
        """String representation."""
        return self.__unicode__()


def prepend_to_env_search_path(varname, value, envdict):
    """Add value to the environment variable varname in envdict

    e.g. prepend_to_env_search_path('BIBINPUTS', '/home/sally/foo', os.environ)
    """
    if not value:
        return  # Nothing to add

    envdict[varname] = value + os.pathsep + envdict.get(varname, "")


class PDFExporter(LatexExporter):
    """Writer designed to write to PDF files.

    This inherits from `LatexExporter`. It creates a LaTeX file in
    a temporary directory using the template machinery, and then runs LaTeX
    to create a pdf.
    """

    export_from_notebook = "PDF via LaTeX"

    latex_count = Integer(3, help="How many times latex will be called.").tag(config=True)

    latex_command = List(
        ["xelatex", "{filename}", "-quiet"], help="Shell command used to compile latex."
    ).tag(config=True)

    bib_command = List(["bibtex", "{filename}"], help="Shell command used to run bibtex.").tag(
        config=True
    )

    verbose = Bool(False, help="Whether to display the output of latex commands.").tag(config=True)

    texinputs = Unicode(help="texinputs dir. A notebook's directory is added")
    writer = Instance("nbconvert.writers.FilesWriter", args=(), kw={"build_directory": "."})

    output_mimetype = "application/pdf"

    _captured_output = List(Unicode())

    @default("file_extension")
    def _file_extension_default(self):
        return ".pdf"

    @default("template_extension")
    def _template_extension_default(self):
        return ".tex.j2"

    def run_command(self, command_list, filename, count, log_function, raise_on_failure=None):
        """Run command_list count times.

        Parameters
        ----------
        command_list : list
            A list of args to provide to Popen. Each element of this
            list will be interpolated with the filename to convert.
        filename : unicode
            The name of the file to convert.
        count : int
            How many times to run the command.
        raise_on_failure: Exception class (default None)
            If provided, will raise the given exception for if an instead of
            returning False on command failure.

        Returns
        -------
        success : bool
            A boolean indicating if the command was successful (True)
            or failed (False).
        """
        command = [c.format(filename=filename) for c in command_list]

        # This will throw a clearer error if the command is not found
        cmd = shutil.which(command_list[0])
        if cmd is None:
            link = "https://nbconvert.readthedocs.io/en/latest/install.html#installing-tex"
            msg = (
                f"{command_list[0]} not found on PATH, if you have not installed "
                f"{command_list[0]} you may need to do so. Find further instructions "
                f"at {link}."
            )
            raise OSError(msg)

        times = "time" if count == 1 else "times"
        self.log.info("Running %s %i %s: %s", command_list[0], count, times, command)

        shell = sys.platform == "win32"
        if shell:
            command = subprocess.list2cmdline(command)  # type:ignore[assignment]
        env = os.environ.copy()
        prepend_to_env_search_path("TEXINPUTS", self.texinputs, env)
        prepend_to_env_search_path("BIBINPUTS", self.texinputs, env)
        prepend_to_env_search_path("BSTINPUTS", self.texinputs, env)

        with open(os.devnull, "rb") as null:
            stdout = subprocess.PIPE if not self.verbose else None
            for _ in range(count):
                p = subprocess.Popen(  # noqa: S603
                    command,
                    stdout=stdout,
                    stderr=subprocess.STDOUT,
                    stdin=null,
                    shell=shell,
                    env=env,
                )
                out, _ = p.communicate()
                if p.returncode:
                    if self.verbose:  # noqa: SIM108
                        # verbose means I didn't capture stdout with PIPE,
                        # so it's already been displayed and `out` is None.
                        out_str = ""
                    else:
                        out_str = out.decode("utf-8", "replace")
                    log_function(command, out)
                    self._captured_output.append(out_str)
                    if raise_on_failure:
                        msg = f'Failed to run "{command}" command:\n{out_str}'
                        raise raise_on_failure(msg)
                    return False  # failure
        return True  # success

    def run_latex(self, filename, raise_on_failure=LatexFailed):
        """Run xelatex self.latex_count times."""

        def log_error(command, out):
            self.log.critical("%s failed: %s\n%s", command[0], command, out)

        return self.run_command(
            self.latex_command, filename, self.latex_count, log_error, raise_on_failure
        )

    def run_bib(self, filename, raise_on_failure=False):
        """Run bibtex one time."""
        filename = os.path.splitext(filename)[0]

        def log_error(command, out):
            self.log.warning(
                "%s had problems, most likely because there were no citations", command[0]
            )
            self.log.debug("%s output: %s\n%s", command[0], command, out)

        return self.run_command(self.bib_command, filename, 1, log_error, raise_on_failure)

    def from_notebook_node(self, nb, resources=None, **kw):
        """Convert from notebook node."""
        latex, resources = super().from_notebook_node(nb, resources=resources, **kw)
        # set texinputs directory, so that local files will be found
        if resources and resources.get("metadata", {}).get("path"):
            self.texinputs = os.path.abspath(resources["metadata"]["path"])
        else:
            self.texinputs = os.getcwd()

        self._captured_outputs = []
        with TemporaryDirectory() as td, _contextlib_chdir.chdir(td):
            notebook_name = "notebook"
            resources["output_extension"] = ".tex"
            tex_file = self.writer.write(latex, resources, notebook_name=notebook_name)
            self.log.info("Building PDF")
            self.run_latex(tex_file)
            if self.run_bib(tex_file):
                self.run_latex(tex_file)

            pdf_file = notebook_name + ".pdf"
            if not os.path.isfile(pdf_file):
                raise LatexFailed("\n".join(self._captured_output))
            self.log.info("PDF successfully created")
            with open(pdf_file, "rb") as f:
                pdf_data = f.read()

        # convert output extension to pdf
        # the writer above required it to be tex
        resources["output_extension"] = ".pdf"
        # clear figure outputs and attachments, extracted by latex export,
        # so we don't claim to be a multi-file export.
        resources.pop("outputs", None)
        resources.pop("attachments", None)

        return pdf_data, resources
