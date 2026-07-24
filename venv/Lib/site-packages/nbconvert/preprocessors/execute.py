# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

"""Module containing a preprocessor that executes the code cells
and updates outputs"""

from __future__ import annotations

import typing as t
from warnings import warn

from jupyter_client.manager import KernelManager
from nbclient.client import NotebookClient
from nbclient.client import execute as _execute

# Backwards compatibility for imported name
from nbclient.exceptions import CellExecutionError  # noqa: F401
from nbformat import NotebookNode

from .base import Preprocessor


def executenb(*args, **kwargs):
    """DEPRECATED."""

    warn(
        "The 'nbconvert.preprocessors.execute.executenb' function was moved to nbclient.execute. "
        "We recommend importing that library directly.",
        FutureWarning,
        stacklevel=2,
    )
    return _execute(*args, **kwargs)


# We inherit from both classes to allow for traitlets to resolve as they did pre-6.0.
# This unfortunately makes for some ugliness around initialization as NotebookClient
# assumes it's a constructed class with a nb object that we have to hack around.
class ExecutePreprocessor(Preprocessor, NotebookClient):
    """
    Executes all the cells in a notebook
    """

    def __init__(self, **kw):
        """Initialize the preprocessor."""
        nb = kw.get("nb")
        if nb is None:
            nb = NotebookNode()
        Preprocessor.__init__(self, nb=nb, **kw)
        NotebookClient.__init__(self, nb, **kw)

    def _check_assign_resources(self, resources):
        if resources or not hasattr(self, "resources"):
            self.resources = resources

    def preprocess(
        self, nb: NotebookNode, resources: t.Any = None, km: KernelManager | None = None
    ) -> tuple[NotebookNode, dict[str, t.Any]]:
        """
        Preprocess notebook executing each code cell.

        The input argument *nb* is modified in-place.

        Note that this function recalls NotebookClient.__init__, which may look wrong.
        However since the preprocess call acts line an init on execution state it's expected.
        Therefore, we need to capture it here again to properly reset because traitlet
        assignments are not passed. There is a risk if traitlets apply any side effects for
        dual init.
        The risk should be manageable, and this approach minimizes side-effects relative
        to other alternatives.

        One alternative but rejected implementation would be to copy the client's init internals
        which has already gotten out of sync with nbclient 0.5 release before nbconvert 6.0 released.

        Parameters
        ----------
        nb : NotebookNode
            Notebook being executed.
        resources : dictionary (optional)
            Additional resources used in the conversion process. For example,
            passing ``{'metadata': {'path': run_path}}`` sets the
            execution path to ``run_path``.
        km: KernelManager (optional)
            Optional kernel manager. If none is provided, a kernel manager will
            be created.

        Returns
        -------
        nb : NotebookNode
            The executed notebook.
        resources : dictionary
            Additional resources used in the conversion process.
        """
        NotebookClient.__init__(self, nb, km)
        self.reset_execution_trackers()
        self._check_assign_resources(resources)

        with self.setup_kernel():
            assert self.kc
            info_msg = self.wait_for_reply(self.kc.kernel_info())
            assert info_msg
            self.nb.metadata["language_info"] = info_msg["content"]["language_info"]
            for index, cell in enumerate(self.nb.cells):
                self.preprocess_cell(cell, resources, index)
        self.set_widgets_metadata()

        return self.nb, self.resources

    def preprocess_cell(self, cell, resources, index):
        """
        Override if you want to apply some preprocessing to each cell.
        Must return modified cell and resource dictionary.

        Parameters
        ----------
        cell : NotebookNode cell
            Notebook cell being processed
        resources : dictionary
            Additional resources used in the conversion process.  Allows
            preprocessors to pass variables into the Jinja engine.
        index : int
            Index of the cell being processed
        """
        self._check_assign_resources(resources)
        cell = self.execute_cell(cell, index, store_history=True)
        return cell, self.resources
