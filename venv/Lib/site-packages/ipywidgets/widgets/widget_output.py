# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Output class.

Represents a widget that can be used to display output within the widget area.
"""

import sys
from functools import wraps

from .domwidget import DOMWidget
from .trait_types import TypedTuple
from .widget import register
from .._version import __jupyter_widgets_output_version__

from traitlets import Unicode, Dict
from IPython.core.interactiveshell import InteractiveShell
from IPython.display import clear_output
from IPython import get_ipython
import traceback

@register
class Output(DOMWidget):
    """Widget used as a context manager to display output.

    This widget can capture and display stdout, stderr, and rich output.  To use
    it, create an instance of it and display it.

    You can then use the widget as a context manager: any output produced while in the
    context will be captured and displayed in the widget instead of the standard output
    area.

    You can also use the .capture() method to decorate a function or a method. Any output
    produced by the function will then go to the output widget. This is useful for
    debugging widget callbacks, for example.

    Example::
        import ipywidgets as widgets
        from IPython.display import display
        out = widgets.Output()
        display(out)

        print('prints to output area')

        with out:
            print('prints to output widget')

        @out.capture()
        def func():
            print('prints to output widget')
    """
    _view_name = Unicode('OutputView').tag(sync=True)
    _model_name = Unicode('OutputModel').tag(sync=True)
    _view_module = Unicode('@jupyter-widgets/output').tag(sync=True)
    _model_module = Unicode('@jupyter-widgets/output').tag(sync=True)
    _view_module_version = Unicode(__jupyter_widgets_output_version__).tag(sync=True)
    _model_module_version = Unicode(__jupyter_widgets_output_version__).tag(sync=True)

    msg_id = Unicode('', help="Parent message id of messages to capture").tag(sync=True)
    outputs = TypedTuple(trait=Dict(), help="The output messages synced from the frontend.").tag(sync=True)

    __counter = 0

    def clear_output(self, *pargs, **kwargs):
        """
        Clear the content of the output widget.

        Parameters
        ----------

        wait: bool
            If True, wait to clear the output until new output is
            available to replace it. Default: False
        """
        with self:
            clear_output(*pargs, **kwargs)

    # PY3: Force passing clear_output and clear_kwargs as kwargs
    def capture(self, clear_output=False, *clear_args, **clear_kwargs):
        """
        Decorator to capture the stdout and stderr of a function.

        Parameters
        ----------

        clear_output: bool
            If True, clear the content of the output widget at every
            new function call. Default: False

        wait: bool
            If True, wait to clear the output until new output is
            available to replace it. This is only used if clear_output
            is also True.
            Default: False
        """
        def capture_decorator(func):
            @wraps(func)
            def inner(*args, **kwargs):
                if clear_output:
                    self.clear_output(*clear_args, **clear_kwargs)
                with self:
                    return func(*args, **kwargs)
            return inner
        return capture_decorator

    def __enter__(self):
        """Called upon entering output widget context manager."""
        self._flush()
        ip = get_ipython()
        kernel = None
        if ip and getattr(ip, "kernel", None) is not None:
            kernel = ip.kernel
        elif self.comm is not None and getattr(self.comm, 'kernel', None) is not None:
            kernel = self.comm.kernel

        if kernel:
            parent = None
            if hasattr(kernel, "get_parent"):
                parent = kernel.get_parent()
            elif hasattr(kernel, "_parent_header"):
                # ipykernel < 6: kernel._parent_header is the parent *request*
                parent = kernel._parent_header

            if parent and parent.get("header"):
                self.msg_id = parent["header"]["msg_id"]
                self.__counter += 1

    def __exit__(self, etype, evalue, tb):
        """Called upon exiting output widget context manager."""
        kernel = None
        if etype is not None:
            ip = get_ipython()
            if ip:
                kernel = ip
                ip.showtraceback((etype, evalue, tb), tb_offset=0)
            elif (self.comm is not None and
                    getattr(self.comm, "kernel", None) is not None and
                    # Check if it's ipykernel
                    getattr(self.comm.kernel, "send_response", None) is not None):
                kernel = self.comm.kernel
                kernel.send_response(kernel.iopub_socket,
                                     u'error',
                                     {
                    u'traceback': ["".join(traceback.format_exception(etype, evalue, tb))],
                    u'evalue': repr(evalue.args),
                    u'ename': etype.__name__
                    })
        self._flush()
        self.__counter -= 1
        if self.__counter == 0:
            self.msg_id = ''
        # suppress exceptions when in IPython, since they are shown above,
        # otherwise let someone else handle it
        return True if kernel else None

    def _flush(self):
        """Flush stdout and stderr buffers."""
        sys.stdout.flush()
        sys.stderr.flush()

    def _append_stream_output(self, text, stream_name):
        """Append a stream output."""
        self.outputs += (
            {'output_type': 'stream', 'name': stream_name, 'text': text},
        )

    def append_stdout(self, text):
        """Append text to the stdout stream."""
        self._append_stream_output(text, stream_name='stdout')

    def append_stderr(self, text):
        """Append text to the stderr stream."""
        self._append_stream_output(text, stream_name='stderr')

    def append_display_data(self, display_object):
        """Append a display object as an output.

        Parameters
        ----------
        display_object : IPython.core.display.DisplayObject
            The object to display (e.g., an instance of
            `IPython.display.Markdown` or `IPython.display.Image`).
        """
        fmt = InteractiveShell.instance().display_formatter.format
        data, metadata = fmt(display_object)
        self.outputs += (
            {
                'output_type': 'display_data',
                'data': data,
                'metadata': metadata
            },
        )
