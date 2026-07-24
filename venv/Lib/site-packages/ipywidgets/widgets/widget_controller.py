# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Controller class.

Represents a Gamepad or Joystick controller.
"""

from .valuewidget import ValueWidget
from .widget import register, widget_serialization
from .domwidget import DOMWidget
from .widget_core import CoreWidget
from .trait_types import TypedTuple
from traitlets import Bool, Int, Float, Unicode, Instance


@register
class Button(DOMWidget, ValueWidget, CoreWidget):
    """Represents a gamepad or joystick button."""
    value = Float(min=0.0, max=1.0, read_only=True, help="The value of the button.").tag(sync=True)
    pressed = Bool(read_only=True, help="Whether the button is pressed.").tag(sync=True)

    _view_name = Unicode('ControllerButtonView').tag(sync=True)
    _model_name = Unicode('ControllerButtonModel').tag(sync=True)


@register
class Axis(DOMWidget, ValueWidget, CoreWidget):
    """Represents a gamepad or joystick axis."""
    value = Float(min=-1.0, max=1.0, read_only=True, help="The value of the axis.").tag(sync=True)

    _view_name = Unicode('ControllerAxisView').tag(sync=True)
    _model_name = Unicode('ControllerAxisModel').tag(sync=True)


@register
class Controller(DOMWidget, CoreWidget):
    """Represents a game controller."""
    index = Int(help="The id number of the controller.").tag(sync=True)

    # General information about the gamepad, button and axes mapping, name.
    # These values are all read-only and set by the JavaScript side.
    name = Unicode(read_only=True, help="The name of the controller.").tag(sync=True)
    mapping = Unicode(read_only=True, help="The name of the control mapping.").tag(sync=True)
    connected = Bool(read_only=True, help="Whether the gamepad is connected.").tag(sync=True)
    timestamp = Float(read_only=True, help="The last time the data from this gamepad was updated.").tag(sync=True)

    # Buttons and axes - read-only
    buttons = TypedTuple(trait=Instance(Button), read_only=True, help="The buttons on the gamepad.").tag(sync=True, **widget_serialization)
    axes = TypedTuple(trait=Instance(Axis), read_only=True, help="The axes on the gamepad.").tag(sync=True, **widget_serialization)

    _view_name = Unicode('ControllerView').tag(sync=True)
    _model_name = Unicode('ControllerModel').tag(sync=True)
