"""Float class.

Represents an unbounded float using a widget.
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from traitlets import (
    Instance, Unicode, CFloat, Bool, CaselessStrEnum, Tuple, TraitError, validate, default
)
from .widget_description import DescriptionWidget
from .trait_types import InstanceDict, NumberFormat
from .valuewidget import ValueWidget
from .widget import register, widget_serialization
from .widget_core import CoreWidget
from .widget_int import ProgressStyle, SliderStyle


class _Float(DescriptionWidget, ValueWidget, CoreWidget):
    value = CFloat(0.0, help="Float value").tag(sync=True)

    def __init__(self, value=None, **kwargs):
        if value is not None:
            kwargs['value'] = value
        super().__init__(**kwargs)


class _BoundedFloat(_Float):
    max = CFloat(100.0, help="Max value").tag(sync=True)
    min = CFloat(0.0, help="Min value").tag(sync=True)

    @validate('value')
    def _validate_value(self, proposal):
        """Cap and floor value"""
        value = proposal['value']
        if self.min > value or self.max < value:
            value = min(max(value, self.min), self.max)
        return value

    @validate('min')
    def _validate_min(self, proposal):
        """Enforce min <= value <= max"""
        min = proposal['value']
        if min > self.max:
            raise TraitError('Setting min > max')
        if min > self.value:
            self.value = min
        return min

    @validate('max')
    def _validate_max(self, proposal):
        """Enforce min <= value <= max"""
        max = proposal['value']
        if max < self.min:
            raise TraitError('setting max < min')
        if max < self.value:
            self.value = max
        return max

class _BoundedLogFloat(_Float):
    max = CFloat(4.0, help="Max value for the exponent").tag(sync=True)
    min = CFloat(0.0, help="Min value for the exponent").tag(sync=True)
    base = CFloat(10.0, help="Base of value").tag(sync=True)
    value = CFloat(1.0, help="Float value").tag(sync=True)

    @validate('value')
    def _validate_value(self, proposal):
        """Cap and floor value"""
        value = proposal['value']
        if self.base ** self.min > value or self.base ** self.max < value:
            value = min(max(value, self.base **  self.min), self.base **  self.max)
        return value

    @validate('min')
    def _validate_min(self, proposal):
        """Enforce base ** min <= value <= base ** max"""
        min = proposal['value']
        if min > self.max:
            raise TraitError('Setting min > max')
        if self.base ** min > self.value:
            self.value = self.base ** min
        return min

    @validate('max')
    def _validate_max(self, proposal):
        """Enforce base ** min <= value <= base ** max"""
        max = proposal['value']
        if max < self.min:
            raise TraitError('setting max < min')
        if self.base ** max < self.value:
            self.value = self.base ** max
        return max


@register
class FloatText(_Float):
    """ Displays a float value within a textbox. For a textbox in
    which the value must be within a specific range, use BoundedFloatText.

    Parameters
    ----------
    value : float
        value displayed
    step : float
        step of the increment (if None, any step is allowed)
    description : str
        description displayed next to the text box
    """
    _view_name = Unicode('FloatTextView').tag(sync=True)
    _model_name = Unicode('FloatTextModel').tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    continuous_update = Bool(False, help="Update the value as the user types. If False, update on submission, e.g., pressing Enter or navigating away.").tag(sync=True)
    step = CFloat(None, allow_none=True, help="Minimum step to increment the value").tag(sync=True)


@register
class BoundedFloatText(_BoundedFloat):
    """ Displays a float value within a textbox. Value must be within the range specified.

    For a textbox in which the value doesn't need to be within a specific range, use FloatText.

    Parameters
    ----------
    value : float
        value displayed
    min : float
        minimal value of the range of possible values displayed
    max : float
        maximal value of the range of possible values displayed
    step : float
        step of the increment (if None, any step is allowed)
    description : str
        description displayed next to the textbox
    """
    _view_name = Unicode('FloatTextView').tag(sync=True)
    _model_name = Unicode('BoundedFloatTextModel').tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    continuous_update = Bool(False, help="Update the value as the user types. If False, update on submission, e.g., pressing Enter or navigating away.").tag(sync=True)
    step = CFloat(None, allow_none=True, help="Minimum step to increment the value").tag(sync=True)

@register
class FloatSlider(_BoundedFloat):
    """ Slider/trackbar of floating values with the specified range.

    Parameters
    ----------
    value : float
        position of the slider
    min : float
        minimal position of the slider
    max : float
        maximal position of the slider
    step : float
        step of the trackbar
    description : str
        name of the slider
    orientation : {'horizontal', 'vertical'}
        default is 'horizontal', orientation of the slider
    readout : {True, False}
        default is True, display the current value of the slider next to it
    behavior : str
        slider handle and connector dragging behavior. Default is 'drag-tap'.
    readout_format : str
        default is '.2f', specifier for the format function used to represent
        slider value for human consumption, modeled after Python 3's format
        specification mini-language (PEP 3101).
    """
    _view_name = Unicode('FloatSliderView').tag(sync=True)
    _model_name = Unicode('FloatSliderModel').tag(sync=True)
    step = CFloat(0.1, allow_none=True, help="Minimum step to increment the value").tag(sync=True)
    orientation = CaselessStrEnum(values=['horizontal', 'vertical'],
        default_value='horizontal', help="Vertical or horizontal.").tag(sync=True)
    readout = Bool(True, help="Display the current value of the slider next to it.").tag(sync=True)
    readout_format = NumberFormat(
        '.2f', help="Format for the readout").tag(sync=True)
    continuous_update = Bool(True, help="Update the value of the widget as the user is holding the slider.").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    style = InstanceDict(SliderStyle).tag(sync=True, **widget_serialization)
    behavior = CaselessStrEnum(values=['drag-tap', 'drag-snap', 'tap', 'drag', 'snap'],
        default_value='drag-tap', help="Slider dragging behavior.").tag(sync=True)


@register
class FloatLogSlider(_BoundedLogFloat):
    """ Slider/trackbar of logarithmic floating values with the specified range.

    Parameters
    ----------
    value : float
        position of the slider
    base : float
        base of the logarithmic scale. Default is 10
    min : float
        minimal position of the slider in log scale, i.e., actual minimum is base ** min
    max : float
        maximal position of the slider in log scale, i.e., actual maximum is base ** max
    step : float
        step of the trackbar, denotes steps for the exponent, not the actual value
    description : str
        name of the slider
    orientation : {'horizontal', 'vertical'}
        default is 'horizontal', orientation of the slider
    readout : {True, False}
        default is True, display the current value of the slider next to it
    behavior : str
        slider handle and connector dragging behavior. Default is 'drag-tap'.
    readout_format : str
        default is '.3g', specifier for the format function used to represent
        slider value for human consumption, modeled after Python 3's format
        specification mini-language (PEP 3101).
    """
    _view_name = Unicode('FloatLogSliderView').tag(sync=True)
    _model_name = Unicode('FloatLogSliderModel').tag(sync=True)
    step = CFloat(0.1, allow_none=True, help="Minimum step in the exponent to increment the value").tag(sync=True)
    orientation = CaselessStrEnum(values=['horizontal', 'vertical'],
        default_value='horizontal', help="Vertical or horizontal.").tag(sync=True)
    readout = Bool(True, help="Display the current value of the slider next to it.").tag(sync=True)
    readout_format = NumberFormat(
        '.3g', help="Format for the readout").tag(sync=True)
    continuous_update = Bool(True, help="Update the value of the widget as the user is holding the slider.").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    base = CFloat(10., help="Base for the logarithm").tag(sync=True)
    style = InstanceDict(SliderStyle).tag(sync=True, **widget_serialization)
    behavior = CaselessStrEnum(values=['drag-tap', 'drag-snap', 'tap', 'drag', 'snap'],
        default_value='drag-tap', help="Slider dragging behavior.").tag(sync=True)


@register
class FloatProgress(_BoundedFloat):
    """ Displays a progress bar.

    Parameters
    -----------
    value : float
        position within the range of the progress bar
    min : float
        minimal position of the slider
    max : float
        maximal position of the slider
    description : str
        name of the progress bar
    orientation : {'horizontal', 'vertical'}
        default is 'horizontal', orientation of the progress bar
    bar_style: {'success', 'info', 'warning', 'danger', ''}
        color of the progress bar, default is '' (blue)
        colors are: 'success'-green, 'info'-light blue, 'warning'-orange, 'danger'-red
    """
    _view_name = Unicode('ProgressView').tag(sync=True)
    _model_name = Unicode('FloatProgressModel').tag(sync=True)
    orientation = CaselessStrEnum(values=['horizontal', 'vertical'],
        default_value='horizontal', help="Vertical or horizontal.").tag(sync=True)

    bar_style = CaselessStrEnum(
        values=['success', 'info', 'warning', 'danger', ''],
        default_value='', allow_none=True,
        help="Use a predefined styling for the progress bar.").tag(sync=True)

    style = InstanceDict(ProgressStyle).tag(sync=True, **widget_serialization)


class _FloatRange(_Float):
    value = Tuple(CFloat(), CFloat(), default_value=(0.0, 1.0),
                  help="Tuple of (lower, upper) bounds").tag(sync=True)

    @property
    def lower(self):
        return self.value[0]

    @lower.setter
    def lower(self, lower):
        self.value = (lower, self.value[1])

    @property
    def upper(self):
        return self.value[1]

    @upper.setter
    def upper(self, upper):
        self.value = (self.value[0], upper)

    @validate('value')
    def _validate_value(self, proposal):
        lower, upper = proposal['value']
        if upper < lower:
            raise TraitError('setting lower > upper')
        return lower, upper


class _BoundedFloatRange(_FloatRange):
    step = CFloat(1.0, help="Minimum step that the value can take (ignored by some views)").tag(sync=True)
    max = CFloat(100.0, help="Max value").tag(sync=True)
    min = CFloat(0.0, help="Min value").tag(sync=True)

    def __init__(self, *args, **kwargs):
        min, max = kwargs.get('min', 0.0), kwargs.get('max', 100.0)
        if kwargs.get('value', None) is None:
            kwargs['value'] = (0.75 * min + 0.25 * max,
                               0.25 * min + 0.75 * max)
        elif not isinstance(kwargs['value'], tuple):
            try:
                kwargs['value'] = tuple(kwargs['value'])
            except:
                raise TypeError(
                    "A 'range' must be able to be cast to a tuple. The input of type"
                    " {} could not be cast to a tuple".format(type(kwargs['value']))
                )
        super().__init__(*args, **kwargs)

    @validate('min', 'max')
    def _validate_bounds(self, proposal):
        trait = proposal['trait']
        new = proposal['value']
        if trait.name == 'min' and new > self.max:
            raise TraitError('setting min > max')
        if trait.name == 'max' and new < self.min:
            raise TraitError('setting max < min')
        if trait.name == 'min':
            self.value = (max(new, self.value[0]), max(new, self.value[1]))
        if trait.name == 'max':
            self.value = (min(new, self.value[0]), min(new, self.value[1]))
        return new

    @validate('value')
    def _validate_value(self, proposal):
        lower, upper = super()._validate_value(proposal)
        lower, upper = min(lower, self.max), min(upper, self.max)
        lower, upper = max(lower, self.min), max(upper, self.min)
        return lower, upper


@register
class FloatRangeSlider(_BoundedFloatRange):
    """ Slider/trackbar that represents a pair of floats bounded by minimum and maximum value.

    Parameters
    ----------
    value : float tuple
        range of the slider displayed
    min : float
        minimal position of the slider
    max : float
        maximal position of the slider
    step : float
        step of the trackbar
    description : str
        name of the slider
    orientation : {'horizontal', 'vertical'}
        default is 'horizontal'
    readout : {True, False}
        default is True, display the current value of the slider next to it
    behavior : str
        slider handle and connector dragging behavior. Default is 'drag-tap'.
    readout_format : str
        default is '.2f', specifier for the format function used to represent
        slider value for human consumption, modeled after Python 3's format
        specification mini-language (PEP 3101).
    """
    _view_name = Unicode('FloatRangeSliderView').tag(sync=True)
    _model_name = Unicode('FloatRangeSliderModel').tag(sync=True)
    step = CFloat(0.1, allow_none=True, help="Minimum step to increment the value").tag(sync=True)
    orientation = CaselessStrEnum(values=['horizontal', 'vertical'],
        default_value='horizontal', help="Vertical or horizontal.").tag(sync=True)
    readout = Bool(True, help="Display the current value of the slider next to it.").tag(sync=True)
    readout_format = NumberFormat(
        '.2f', help="Format for the readout").tag(sync=True)
    continuous_update = Bool(True, help="Update the value of the widget as the user is sliding the slider.").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    style = InstanceDict(SliderStyle).tag(sync=True, **widget_serialization)
    behavior = CaselessStrEnum(values=['drag-tap', 'drag-snap', 'tap', 'drag', 'snap'],
        default_value='drag-tap', help="Slider dragging behavior.").tag(sync=True)
