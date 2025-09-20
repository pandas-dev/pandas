# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Interactive widgets for the Jupyter notebook.

Provide simple interactive controls in the notebook.
Each Widget corresponds to an object in Python and Javascript,
with controls on the page.

To put a Widget on the page, you can display it with Jupyter's display machinery::

    from ipywidgets import IntSlider
    slider = IntSlider(min=1, max=10)
    display(slider)

Moving the slider will change the value. Most Widgets have a current value,
accessible as a `value` attribute.
"""

# Must import __version__ first to avoid errors importing this file during the build process. See https://github.com/pypa/setuptools/issues/1724#issuecomment-627241822
from ._version import __version__, __protocol_version__, __jupyter_widgets_controls_version__, __jupyter_widgets_base_version__

import os
import sys

from traitlets import link, dlink
from IPython import get_ipython

from .widgets import *


def load_ipython_extension(ip):
    """Set up Jupyter to work with widgets"""
    if not hasattr(ip, 'kernel'):
        return
    register_comm_target()

def register_comm_target(kernel=None):
    """Register the jupyter.widget comm target"""
    from . import comm
    comm_manager = comm.get_comm_manager()
    if comm_manager is None:
        return
    comm_manager.register_target('jupyter.widget', Widget.handle_comm_opened)
    comm_manager.register_target('jupyter.widget.control', Widget.handle_control_comm_opened)

def _handle_ipython():
    """Register with the comm target at import if running in Jupyter"""
    ip = get_ipython()
    if ip is None:
        return
    register_comm_target()

_handle_ipython()
