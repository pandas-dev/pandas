# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Contains the Style class"""

from traitlets import Unicode
from .widget import Widget
from .._version import __jupyter_widgets_base_version__

class Style(Widget):
    """Style specification"""

    _model_name = Unicode('StyleModel').tag(sync=True)
    _view_name = Unicode('StyleView').tag(sync=True)
    _view_module = Unicode('@jupyter-widgets/base').tag(sync=True)
    _view_module_version = Unicode(__jupyter_widgets_base_version__).tag(sync=True)
