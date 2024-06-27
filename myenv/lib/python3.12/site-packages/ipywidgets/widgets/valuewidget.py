# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Contains the ValueWidget class"""

from .widget import Widget
from traitlets import Any


class ValueWidget(Widget):
    """Widget that can be used for the input of an interactive function"""

    value = Any(help="The value of the widget.")

    def get_interact_value(self):
        """Return the value for this widget which should be passed to
        interactive functions. Custom widgets can change this method
        to process the raw value ``self.value``.
        """
        return self.value

    def _repr_keys(self):
        # Ensure value key comes first, and is always present
        yield 'value'
        for key in super()._repr_keys():
            if key != 'value':
                yield key
