# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""
Time and datetime picker widgets
"""

from traitlets import Unicode, Bool, validate, TraitError

from .trait_types import datetime_serialization, Datetime, naive_serialization
from .valuewidget import ValueWidget
from .widget import register
from .widget_core import CoreWidget
from .widget_description import DescriptionWidget


@register
class DatetimePicker(DescriptionWidget, ValueWidget, CoreWidget):
    """
    Display a widget for picking datetimes.

    Parameters
    ----------

    value: datetime.datetime
        The current value of the widget.

    disabled: bool
        Whether to disable user changes.

    min: datetime.datetime
        The lower allowed datetime bound

    max: datetime.datetime
        The upper allowed datetime bound

    Examples
    --------

    >>> import datetime
    >>> import ipydatetime
    >>> datetime_pick = ipydatetime.DatetimePicker()
    >>> datetime_pick.value = datetime.datetime(2018, 09, 5, 12, 34, 3)
    """

    _view_name = Unicode("DatetimeView").tag(sync=True)
    _model_name = Unicode("DatetimeModel").tag(sync=True)

    value = Datetime(None, allow_none=True).tag(sync=True, **datetime_serialization)
    disabled = Bool(False, help="Enable or disable user changes.").tag(sync=True)

    min = Datetime(None, allow_none=True).tag(sync=True, **datetime_serialization)
    max = Datetime(None, allow_none=True).tag(sync=True, **datetime_serialization)

    def _validate_tz(self, value):
        if value.tzinfo is None:
            raise TraitError('%s values needs to be timezone aware' % (self.__class__.__name__,))
        return value

    @validate("value")
    def _validate_value(self, proposal):
        """Cap and floor value"""
        value = proposal["value"]
        if value is None:
            return value
        value = self._validate_tz(value)
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
        min = self._validate_tz(min)
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
        max = self._validate_tz(max)
        if self.min and max < self.min:
            raise TraitError("setting max < min")
        if self.value and max < self.value:
            self.value = max
        return max


@register
class NaiveDatetimePicker(DatetimePicker):
    """
    Display a widget for picking naive datetimes (i.e. timezone unaware).

    Parameters
    ----------

    value: datetime.datetime
        The current value of the widget.

    disabled: bool
        Whether to disable user changes.

    min: datetime.datetime
        The lower allowed datetime bound

    max: datetime.datetime
        The upper allowed datetime bound

    Examples
    --------

    >>> import datetime
    >>> import ipydatetime
    >>> datetime_pick = ipydatetime.NaiveDatetimePicker()
    >>> datetime_pick.value = datetime.datetime(2018, 09, 5, 12, 34, 3)
    """

    # Replace the serializers and model names:

    _model_name = Unicode("NaiveDatetimeModel").tag(sync=True)

    value = Datetime(None, allow_none=True).tag(sync=True, **naive_serialization)

    min = Datetime(None, allow_none=True).tag(sync=True, **naive_serialization)
    max = Datetime(None, allow_none=True).tag(sync=True, **naive_serialization)

    def _validate_tz(self, value):
        if value.tzinfo is not None:
            raise TraitError('%s values needs to be timezone unaware' % (self.__class__.__name__,))
        return value
