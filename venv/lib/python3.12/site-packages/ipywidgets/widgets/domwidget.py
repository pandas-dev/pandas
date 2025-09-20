# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Contains the DOMWidget class"""

from traitlets import Bool, Unicode
from .widget import Widget, widget_serialization
from .trait_types import InstanceDict, TypedTuple
from .widget_layout import Layout
from .widget_style import Style


class DOMWidget(Widget):
    """Widget that can be inserted into the DOM

    Parameters
    ----------
    tooltip: str
       tooltip caption
    layout: InstanceDict(Layout)
       widget layout
    """

    _model_name = Unicode('DOMWidgetModel').tag(sync=True)
    _dom_classes = TypedTuple(trait=Unicode(), help="CSS classes applied to widget DOM element").tag(sync=True)
    tabbable = Bool(help="Is widget tabbable?", allow_none=True, default_value=None).tag(sync=True)
    tooltip = Unicode(None, allow_none=True, help="A tooltip caption.").tag(sync=True)
    layout = InstanceDict(Layout).tag(sync=True, **widget_serialization)

    def add_class(self, className):
        """
        Adds a class to the top level element of the widget.

        Doesn't add the class if it already exists.
        """
        if className not in self._dom_classes:
            self._dom_classes = list(self._dom_classes) + [className]
        return self

    def remove_class(self, className):
        """
        Removes a class from the top level element of the widget.

        Doesn't remove the class if it doesn't exist.
        """
        if className in self._dom_classes:
            self._dom_classes = [c for c in self._dom_classes if c != className]
        return self

    def focus(self):
        """
        Focus on the widget.
        """
        self.send({'do':'focus'})

    def blur(self):
        """
        Blur the widget.
        """
        self.send({'do':'blur'})

    def _repr_keys(self):
        for key in super()._repr_keys():
            # Exclude layout if it had the default value
            if key == 'layout':
                value = getattr(self, key)
                if repr(value) == '%s()' % value.__class__.__name__:
                    continue
            yield key
        # We also need to include _dom_classes in repr for reproducibility
        if self._dom_classes:
            yield '_dom_classes'
