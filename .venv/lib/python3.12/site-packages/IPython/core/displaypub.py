"""An interface for publishing rich data to frontends.

There are two components of the display system:

* Display formatters, which take a Python object and compute the
  representation of the object in various formats (text, HTML, SVG, etc.).
* The display publisher that is used to send the representation data to the
  various frontends.

This module defines the logic display publishing. The display publisher uses
the ``display_data`` message type that is defined in the IPython messaging
spec.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import sys

from traitlets.config.configurable import Configurable
from traitlets import List

# This used to be defined here - it is imported for backwards compatibility
from .display_functions import publish_display_data
from .history import HistoryOutput

import typing as t

# -----------------------------------------------------------------------------
# Main payload class
# -----------------------------------------------------------------------------

_sentinel = object()


class DisplayPublisher(Configurable):
    """A traited class that publishes display data to frontends.

    Instances of this class are created by the main IPython object and should
    be accessed there.
    """

    def __init__(self, shell=None, *args, **kwargs):
        self.shell = shell
        self._is_publishing = False
        self._in_post_execute = False
        if self.shell:
            self._setup_execution_tracking()
        super().__init__(*args, **kwargs)

    def _validate_data(self, data, metadata=None):
        """Validate the display data.

        Parameters
        ----------
        data : dict
            The formata data dictionary.
        metadata : dict
            Any metadata for the data.
        """

        if not isinstance(data, dict):
            raise TypeError("data must be a dict, got: %r" % data)
        if metadata is not None:
            if not isinstance(metadata, dict):
                raise TypeError("metadata must be a dict, got: %r" % data)

    def _setup_execution_tracking(self):
        """Set up hooks to track execution state"""
        self.shell.events.register("post_execute", self._on_post_execute)
        self.shell.events.register("pre_execute", self._on_pre_execute)

    def _on_post_execute(self):
        """Called at start of post_execute phase"""
        self._in_post_execute = True

    def _on_pre_execute(self):
        """Called at start of pre_execute phase"""
        self._in_post_execute = False

    # use * to indicate transient, update are keyword-only
    def publish(
        self,
        data,
        metadata=None,
        source=_sentinel,
        *,
        transient=None,
        update=False,
        **kwargs,
    ) -> None:
        """Publish data and metadata to all frontends.

        See the ``display_data`` message in the messaging documentation for
        more details about this message type.

        The following MIME types are currently implemented:

        * text/plain
        * text/html
        * text/markdown
        * text/latex
        * application/json
        * application/javascript
        * image/png
        * image/jpeg
        * image/svg+xml

        Parameters
        ----------
        data : dict
            A dictionary having keys that are valid MIME types (like
            'text/plain' or 'image/svg+xml') and values that are the data for
            that MIME type. The data itself must be a JSON'able data
            structure. Minimally all data should have the 'text/plain' data,
            which can be displayed by all frontends. If more than the plain
            text is given, it is up to the frontend to decide which
            representation to use.
        metadata : dict
            A dictionary for metadata related to the data. This can contain
            arbitrary key, value pairs that frontends can use to interpret
            the data.  Metadata specific to each mime-type can be specified
            in the metadata dict with the same mime-type keys as
            the data itself.
        source : str, deprecated
            Unused.
        transient : dict, keyword-only
            A dictionary for transient data.
            Data in this dictionary should not be persisted as part of saving this output.
            Examples include 'display_id'.
        update : bool, keyword-only, default: False
            If True, only update existing outputs with the same display_id,
            rather than creating a new output.
        """

        if source is not _sentinel:
            import warnings

            warnings.warn(
                "The 'source' parameter is deprecated since IPython 3.0 and will be ignored "
                "(this warning is present since 9.0). `source` parameter will be removed in the future.",
                DeprecationWarning,
                stacklevel=2,
            )

        handlers: t.Dict = {}
        if self.shell is not None:
            handlers = getattr(self.shell, "mime_renderers", {})

        outputs = self.shell.history_manager.outputs

        target_execution_count = self.shell.execution_count - 1
        if self._in_post_execute:
            # We're in post_execute, so this is likely a matplotlib flush
            # Use execution_count - 1 to associate with the cell that created the plot
            target_execution_count = self.shell.execution_count - 1

        outputs[target_execution_count].append(
            HistoryOutput(output_type="display_data", bundle=data)
        )

        for mime, handler in handlers.items():
            if mime in data:
                handler(data[mime], metadata.get(mime, None))
                return

        self._is_publishing = True
        if "text/plain" in data:
            print(data["text/plain"])
        self._is_publishing = False

    @property
    def is_publishing(self):
        return self._is_publishing

    def clear_output(self, wait=False):
        """Clear the output of the cell receiving output."""
        print("\033[2K\r", end="")
        sys.stdout.flush()
        print("\033[2K\r", end="")
        sys.stderr.flush()


class CapturingDisplayPublisher(DisplayPublisher):
    """A DisplayPublisher that stores"""

    outputs: List = List()

    def publish(
        self, data, metadata=None, source=None, *, transient=None, update=False
    ):
        self.outputs.append(
            {
                "data": data,
                "metadata": metadata,
                "transient": transient,
                "update": update,
            }
        )

    def clear_output(self, wait=False):
        super(CapturingDisplayPublisher, self).clear_output(wait)

        # empty the list, *do not* reassign a new list
        self.outputs.clear()
