# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Contains the DOMWidget class"""

from traitlets import Bool, Unicode
from .widget import Widget, widget_serialization, register
from .trait_types import InstanceDict
from .widget_style import Style
from .widget_core import CoreWidget
from .domwidget import DOMWidget
from .utils import deprecation

import warnings

@register
class DescriptionStyle(Style, CoreWidget, Widget):
    """Description style widget."""
    _model_name = Unicode('DescriptionStyleModel').tag(sync=True)
    description_width = Unicode(help="Width of the description to the side of the control.").tag(sync=True)


class DescriptionWidget(DOMWidget, CoreWidget):
    """Widget that has a description label to the side."""
    _model_name = Unicode('DescriptionModel').tag(sync=True)
    description = Unicode('', help="Description of the control.").tag(sync=True)
    description_allow_html = Bool(False, help="Accept HTML in the description.").tag(sync=True)
    style = InstanceDict(DescriptionStyle, help="Styling customizations").tag(sync=True, **widget_serialization)

    def __init__(self, *args, **kwargs):
        if 'description_tooltip' in kwargs:
            deprecation("the description_tooltip argument is deprecated, use tooltip instead")
            kwargs.setdefault('tooltip', kwargs['description_tooltip'])
            del kwargs['description_tooltip']
        super().__init__(*args, **kwargs)

    def _repr_keys(self):
        for key in super()._repr_keys():
            # Exclude style if it had the default value
            if key == 'style':
                value = getattr(self, key)
                if repr(value) == '%s()' % value.__class__.__name__:
                    continue
            yield key

    @property
    def description_tooltip(self):
        """The tooltip information.
        .. deprecated :: 8.0.0
           Use tooltip attribute instead.
        """
        deprecation(".description_tooltip is deprecated, use .tooltip instead")
        return self.tooltip

    @description_tooltip.setter
    def description_tooltip(self, tooltip):
        deprecation(".description_tooltip is deprecated, use .tooltip instead")
        self.tooltip = tooltip
