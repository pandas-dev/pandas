# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Color class.

Represents an HTML Color .
"""

from .widget_description import DescriptionWidget
from .valuewidget import ValueWidget
from .widget import register
from .widget_core import CoreWidget
from .trait_types import Color
from traitlets import Unicode, Bool


@register
class ColorPicker(DescriptionWidget, ValueWidget, CoreWidget):
    value = Color('black', help="The color value.").tag(sync=True)
    concise = Bool(help="Display short version with just a color selector.").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes.").tag(sync=True)

    _view_name = Unicode('ColorPickerView').tag(sync=True)
    _model_name = Unicode('ColorPickerModel').tag(sync=True)
