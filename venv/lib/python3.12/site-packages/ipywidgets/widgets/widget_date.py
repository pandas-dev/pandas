# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Color class.

Represents an HTML Color .
"""

from .widget_description import DescriptionWidget
from .valuewidget import ValueWidget
from .widget import register
from .widget_core import CoreWidget
from .trait_types import Date, date_serialization
from traitlets import Unicode, Bool, Union, CInt, CaselessStrEnum, TraitError, validate


@register
class DatePicker(DescriptionWidget, ValueWidget, CoreWidget):
    """
    Display a widget for picking dates.

    Parameters
    ----------

    value: datetime.date
        The current value of the widget.

    disabled: bool
        Whether to disable user changes.

    Examples
    --------

    >>> import datetime
    >>> import ipywidgets as widgets
    >>> date_pick = widgets.DatePicker()
    >>> date_pick.value = datetime.date(2019, 7, 9)
    """

    _view_name = Unicode('DatePickerView').tag(sync=True)
    _model_name = Unicode('DatePickerModel').tag(sync=True)

    value = Date(None, allow_none=True).tag(sync=True, **date_serialization)
    disabled = Bool(False, help="Enable or disable user changes.").tag(sync=True)

    min = Date(None, allow_none=True).tag(sync=True, **date_serialization)
    max = Date(None, allow_none=True).tag(sync=True, **date_serialization)
    step = Union(
        (CInt(1), CaselessStrEnum(["any"])),
        help='The date step to use for the picker, in days, or "any".',
    ).tag(sync=True)

    @validate("value")
    def _validate_value(self, proposal):
        """Cap and floor value"""
        value = proposal["value"]
        if value is None:
            return value
        if self.min and self.min > value:
            value = max(value, self.min)
        if self.max and self.max < value:
            value = min(value, self.max)
        return value

    @validate("min")
    def _validate_min(self, proposal):
        """Enforce min <= value <= max"""
        min = proposal["value"]
        if min is None:
            return min
        if self.max and min > self.max:
            raise TraitError("Setting min > max")
        if self.value and min > self.value:
            self.value = min
        return min

    @validate("max")
    def _validate_max(self, proposal):
        """Enforce min <= value <= max"""
        max = proposal["value"]
        if max is None:
            return max
        if self.min and max < self.min:
            raise TraitError("setting max < min")
        if self.value and max < self.value:
            self.value = max
        return max
