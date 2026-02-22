"""Read and write notebooks as regular .py files.

Authors:

* Brian Granger
"""

# -----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file LICENSE, distributed as part of this software.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
from __future__ import annotations

import re

from .nbbase import new_code_cell, new_notebook, new_text_cell, new_worksheet
from .rwbase import NotebookReader, NotebookWriter

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

_encoding_declaration_re = re.compile(r"^#.*coding[:=]\s*([-\w.]+)")


class PyReaderError(Exception):
    """An error raised by the PyReader."""


class PyReader(NotebookReader):
    """A Python notebook reader."""

    def reads(self, s, **kwargs):
        """Convert a string to a notebook."""
        return self.to_notebook(s, **kwargs)

    def to_notebook(self, s, **kwargs):
        """Convert a string to a notebook."""
        lines = s.splitlines()
        cells = []
        cell_lines: list[str] = []
        state = "codecell"
        for line in lines:
            if line.startswith("# <nbformat>") or _encoding_declaration_re.match(line):
                pass
            elif line.startswith("# <codecell>"):
                cell = self.new_cell(state, cell_lines)
                if cell is not None:
                    cells.append(cell)
                state = "codecell"
                cell_lines = []
            elif line.startswith("# <htmlcell>"):
                cell = self.new_cell(state, cell_lines)
                if cell is not None:
                    cells.append(cell)
                state = "htmlcell"
                cell_lines = []
            elif line.startswith("# <markdowncell>"):
                cell = self.new_cell(state, cell_lines)
                if cell is not None:
                    cells.append(cell)
                state = "markdowncell"
                cell_lines = []
            else:
                cell_lines.append(line)
        if cell_lines and state == "codecell":
            cell = self.new_cell(state, cell_lines)
            if cell is not None:
                cells.append(cell)
        ws = new_worksheet(cells=cells)
        return new_notebook(worksheets=[ws])

    def new_cell(self, state, lines):
        """Create a new cell."""
        if state == "codecell":
            input_ = "\n".join(lines)
            input_ = input_.strip("\n")
            if input_:
                return new_code_cell(input=input_)
        elif state == "htmlcell":
            text = self._remove_comments(lines)
            if text:
                return new_text_cell("html", source=text)
        elif state == "markdowncell":
            text = self._remove_comments(lines)
            if text:
                return new_text_cell("markdown", source=text)

    def _remove_comments(self, lines):
        new_lines = []
        for line in lines:
            if line.startswith("#"):
                new_lines.append(line[2:])
            else:
                new_lines.append(line)
        text = "\n".join(new_lines)
        text = text.strip("\n")
        return text  # noqa: RET504

    def split_lines_into_blocks(self, lines):
        """Split lines into code blocks."""
        if len(lines) == 1:
            yield lines[0]
            raise StopIteration()
        import ast

        source = "\n".join(lines)
        code = ast.parse(source)
        starts = [x.lineno - 1 for x in code.body]
        for i in range(len(starts) - 1):
            yield "\n".join(lines[starts[i] : starts[i + 1]]).strip("\n")
        yield "\n".join(lines[starts[-1] :]).strip("\n")


class PyWriter(NotebookWriter):
    """A Python notebook writer."""

    def writes(self, nb, **kwargs):
        """Convert a notebook object to a string."""
        lines = ["# -*- coding: utf-8 -*-"]
        lines.extend(["# <nbformat>2</nbformat>", ""])
        for ws in nb.worksheets:
            for cell in ws.cells:
                if cell.cell_type == "code":
                    input_ = cell.get("input")
                    if input_ is not None:
                        lines.extend(["# <codecell>", ""])
                        lines.extend(input_.splitlines())
                        lines.append("")
                elif cell.cell_type == "html":
                    input_ = cell.get("source")
                    if input_ is not None:
                        lines.extend(["# <htmlcell>", ""])
                        lines.extend(["# " + line for line in input_.splitlines()])
                        lines.append("")
                elif cell.cell_type == "markdown":
                    input_ = cell.get("source")
                    if input_ is not None:
                        lines.extend(["# <markdowncell>", ""])
                        lines.extend(["# " + line for line in input_.splitlines()])
                        lines.append("")
        lines.append("")
        return str("\n".join(lines))


_reader = PyReader()
_writer = PyWriter()

reads = _reader.reads
read = _reader.read
to_notebook = _reader.to_notebook
write = _writer.write
writes = _writer.writes
