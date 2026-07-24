# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Bool class.

Represents a boolean using a widget.
"""

from .widget_description import DescriptionStyle, DescriptionWidget
from .widget_core import CoreWidget
from .valuewidget import ValueWidget
from .widget import register, widget_serialization
from .trait_types import Color, InstanceDict
from traitlets import Unicode, Bool, CaselessStrEnum


@register
class CheckboxStyle(DescriptionStyle, CoreWidget):
    """Checkbox widget style."""
    _model_name = Unicode('CheckboxStyleModel').tag(sync=True)
    background = Unicode(None, allow_none=True, help="Background specifications.").tag(sync=True)


@register
class ToggleButtonStyle(DescriptionStyle, CoreWidget):
    """ToggleButton widget style."""
    _model_name = Unicode('ToggleButtonStyleModel').tag(sync=True)
    font_family = Unicode(None, allow_none=True, help="Toggle button text font family.").tag(sync=True)
    font_size = Unicode(None, allow_none=True, help="Toggle button text font size.").tag(sync=True)
    font_style = Unicode(None, allow_none=True, help="Toggle button text font style.").tag(sync=True)
    font_variant = Unicode(None, allow_none=True, help="Toggle button text font variant.").tag(sync=True)
    font_weight = Unicode(None, allow_none=True, help="Toggle button text font weight.").tag(sync=True)
    text_color = Color(None, allow_none=True, help="Toggle button text color").tag(sync=True)
    text_decoration = Unicode(None, allow_none=True, help="Toggle button text decoration.").tag(sync=True)


class _Bool(DescriptionWidget, ValueWidget, CoreWidget):
    """A base class for creating widgets that represent booleans."""
    value = Bool(False, help="Bool value").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes.").tag(sync=True)

    def __init__(self, value=None, **kwargs):
        if value is not None:
            kwargs['value'] = value
        super().__init__(**kwargs)

    _model_name = Unicode('BoolModel').tag(sync=True)


@register
class Checkbox(_Bool):
    """Displays a boolean `value` in the form of a checkbox.

    Parameters
    ----------
    value : {True,False}
        value of the checkbox: True-checked, False-unchecked
    description : str
        description displayed next to the checkbox
    indent : {True,False}
        indent the control to align with other controls with a description. The style.description_width attribute controls this width for consistence with other controls.
    """
    _view_name = Unicode('CheckboxView').tag(sync=True)
    _model_name = Unicode('CheckboxModel').tag(sync=True)
    indent = Bool(True, help="Indent the control to align with other controls with a description.").tag(sync=True)
    style = InstanceDict(CheckboxStyle, help="Styling customizations").tag(sync=True, **widget_serialization)



@register
class ToggleButton(_Bool):
    """Displays a boolean `value` in the form of a toggle button.

    Parameters
    ----------
    value : {True,False}
        value of the toggle button: True-pressed, False-unpressed
    description : str
        description displayed on the button
    icon: str
        font-awesome icon name
    style: instance of DescriptionStyle
        styling customizations
    button_style: enum
        button predefined styling
    """
    _view_name = Unicode('ToggleButtonView').tag(sync=True)
    _model_name = Unicode('ToggleButtonModel').tag(sync=True)

    icon = Unicode('', help= "Font-awesome icon.").tag(sync=True)

    button_style = CaselessStrEnum(
        values=['primary', 'success', 'info', 'warning', 'danger', ''], default_value='',
        help="""Use a predefined styling for the button.""").tag(sync=True)
    style = InstanceDict(ToggleButtonStyle, help="Styling customizations").tag(sync=True, **widget_serialization)


@register
class Valid(_Bool):
    """Displays a boolean `value` in the form of a green check (True / valid)
    or a red cross (False / invalid).

    Parameters
    ----------
    value: {True,False}
        value of the Valid widget
    """
    readout = Unicode('Invalid', help="Message displayed when the value is False").tag(sync=True)
    _view_name = Unicode('ValidView').tag(sync=True)
    _model_name = Unicode('ValidModel').tag(sync=True)
