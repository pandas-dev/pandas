"""Module containing a preprocessor that converts outputs in the notebook from
one format to another.
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import base64
import os
import subprocess
import sys
import warnings
from pathlib import Path
from shutil import which
from tempfile import TemporaryDirectory

from traitlets import List, Unicode, Union, default

from nbconvert.utils.io import FormatSafeDict

from .convertfigures import ConvertFiguresPreprocessor

# inkscape path for darwin (macOS)
INKSCAPE_APP = "/Applications/Inkscape.app/Contents/Resources/bin/inkscape"
# Recent versions of Inkscape (v1.0) moved the executable from
# Resources/bin/inkscape to MacOS/inkscape
INKSCAPE_APP_v1 = "/Applications/Inkscape.app/Contents/MacOS/inkscape"

if sys.platform == "win32":
    try:
        import winreg
    except ImportError:
        import _winreg as winreg


class SVG2PDFPreprocessor(ConvertFiguresPreprocessor):
    """
    Converts all of the outputs in a notebook from SVG to PDF.
    """

    @default("from_format")
    def _from_format_default(self):
        return "image/svg+xml"

    @default("to_format")
    def _to_format_default(self):
        return "application/pdf"

    inkscape_version = Unicode(
        help="""The version of inkscape being used.

        This affects how the conversion command is run.
        """
    ).tag(config=True)

    @default("inkscape_version")
    def _inkscape_version_default(self):
        p = subprocess.Popen(  # noqa:S603
            [self.inkscape, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output, _ = p.communicate()
        if p.returncode != 0:
            msg = "Unable to find inkscape executable --version"
            raise RuntimeError(msg)
        return output.decode("utf-8").split(" ")[1]

    # FIXME: Deprecate passing a string here
    command = Union(
        [Unicode(), List()],
        help="""
        The command to use for converting SVG to PDF

        This traitlet is a template, which will be formatted with the keys
        to_filename and from_filename.

        The conversion call must read the SVG from {from_filename},
        and write a PDF to {to_filename}.

        It could be a List (recommended) or a String. If string, it will
        be passed to a shell for execution.
        """,
    ).tag(config=True)

    @default("command")
    def _command_default(self):
        major_version = self.inkscape_version.split(".")[0]
        command = [self.inkscape]

        if int(major_version) < 1:
            # --without-gui is only needed for inkscape 0.x
            command.append("--without-gui")
            # --export-pdf is old name for --export-filename
            command.append("--export-pdf={to_filename}")
        else:
            command.append("--export-filename={to_filename}")

        command.append("{from_filename}")
        return command

    inkscape = Unicode(help="The path to Inkscape, if necessary").tag(config=True)

    @default("inkscape")
    def _inkscape_default(self):
        # Windows: Secure registry lookup FIRST (CVE-2025-53000 fix)
        if sys.platform == "win32":
            wr_handle = winreg.ConnectRegistry(None, winreg.HKEY_LOCAL_MACHINE)
            try:
                rkey = winreg.OpenKey(wr_handle, r"SOFTWARE\Classes\inkscape.svg\DefaultIcon")
                inkscape_full = winreg.QueryValueEx(rkey, "")[0].split(",")[0]  # Fix: remove ",0"
                if os.path.isfile(inkscape_full):
                    return inkscape_full
            except (FileNotFoundError, OSError, IndexError):
                pass  # Safe fallback

        # Block CWD in PATH search (CVE-2025-53000)
        os.environ["NODEFAULTCURRENTDIRECTORYINEXEPATH"] = "1"

        inkscape_path = which("inkscape")

        # Extra safety for Python < 3.12 on Windows:
        # If which() resolved to a path in CWD even though CWD is not on PATH,
        # warn and treat as "not found".
        if sys.platform == "win32" and inkscape_path and sys.version_info < (3, 12):
            try:
                cwd = Path.cwd().resolve()
                in_cwd = Path(inkscape_path).resolve().parent == cwd
                cwd_on_path = cwd in {
                    Path(p).resolve() for p in os.environ.get("PATH", os.defpath).split(os.pathsep)
                }

                if in_cwd and not cwd_on_path:
                    warnings.warn(
                        "shutil.which('inkscape') resolved to an executable in the current "
                        "working directory even though CWD is not on PATH. Ignoring this "
                        "result for security reasons (CVE-2025-53000).",
                        RuntimeWarning,
                        stacklevel=2,
                    )
                    inkscape_path = None
            except Exception:
                # If detection fails for any reason, prefer safety: ignore CWD result
                inkscape_path = None

        if inkscape_path is not None:
            return inkscape_path

        # macOS: EXACT original order preserved
        if sys.platform == "darwin":
            if os.path.isfile(INKSCAPE_APP_v1):
                return INKSCAPE_APP_v1
            # Order is important. If INKSCAPE_APP exists, prefer it over
            # the executable in the MacOS directory.
            if os.path.isfile(INKSCAPE_APP):
                return INKSCAPE_APP

        msg = "Inkscape executable not found in safe paths"
        raise FileNotFoundError(msg)

    def convert_figure(self, data_format, data):
        """
        Convert a single SVG figure to PDF.  Returns converted data.
        """

        # Work in a temporary directory
        with TemporaryDirectory() as tmpdir:
            # Write fig to temp file
            input_filename = os.path.join(tmpdir, "figure.svg")
            # SVG data is unicode text
            with open(input_filename, "w", encoding="utf8") as f:
                f.write(data)

            # Call conversion application
            output_filename = os.path.join(tmpdir, "figure.pdf")

            template_vars = {"from_filename": input_filename, "to_filename": output_filename}
            if isinstance(self.command, list):
                full_cmd = [s.format_map(FormatSafeDict(**template_vars)) for s in self.command]
            else:
                # For backwards compatibility with specifying strings
                # Okay-ish, since the string is trusted
                full_cmd = self.command.format(**template_vars)
            subprocess.call(full_cmd, shell=isinstance(full_cmd, str))  # noqa: S603

            # Read output from drive
            # return value expects a filename
            if os.path.isfile(output_filename):
                with open(output_filename, "rb") as f:
                    # PDF is a nb supported binary, data type, so base64 encode.
                    return base64.encodebytes(f.read()).decode("utf-8")
            else:
                msg = "Inkscape svg to pdf conversion failed"
                raise TypeError(msg)
