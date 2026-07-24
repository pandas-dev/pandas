# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Button class.

Represents a button in the frontend using a widget.  Allows user to listen for
click events on the button and trigger backend code when the clicks are fired.
"""

from .utils import deprecation
from .domwidget import DOMWidget
from .widget import CallbackDispatcher, register, widget_serialization
from .widget_core import CoreWidget
from .widget_style import Style
from .trait_types import Color, InstanceDict

from traitlets import Unicode, Bool, CaselessStrEnum, Instance, validate, default


@register
class ButtonStyle(Style, CoreWidget):
    """Button style widget."""
    _model_name = Unicode('ButtonStyleModel').tag(sync=True)
    button_color = Color(None, allow_none=True, help="Color of the button").tag(sync=True)
    font_family = Unicode(None, allow_none=True, help="Button text font family.").tag(sync=True)
    font_size = Unicode(None, allow_none=True, help="Button text font size.").tag(sync=True)
    font_style = Unicode(None, allow_none=True, help="Button text font style.").tag(sync=True)
    font_variant = Unicode(None, allow_none=True, help="Button text font variant.").tag(sync=True)
    font_weight = Unicode(None, allow_none=True, help="Button text font weight.").tag(sync=True)
    text_color = Unicode(None, allow_none=True, help="Button text color.").tag(sync=True)
    text_decoration = Unicode(None, allow_none=True, help="Button text decoration.").tag(sync=True)


@register
class Button(DOMWidget, CoreWidget):
    """Button widget.

    This widget has an `on_click` method that allows you to listen for the
    user clicking on the button.  The click event itself is stateless.

    Parameters
    ----------
    description: str
       description displayed on the button
    icon: str
       font-awesome icon names, without the 'fa-' prefix
    disabled: bool
       whether user interaction is enabled
    """
    _view_name = Unicode('ButtonView').tag(sync=True)
    _model_name = Unicode('ButtonModel').tag(sync=True)

    description = Unicode(help="Button label.").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes.").tag(sync=True)
    icon = Unicode('', help="Font-awesome icon names, without the 'fa-' prefix.").tag(sync=True)

    button_style = CaselessStrEnum(
        values=['primary', 'success', 'info', 'warning', 'danger', ''], default_value='',
        help="""Use a predefined styling for the button.""").tag(sync=True)

    style = InstanceDict(ButtonStyle).tag(sync=True, **widget_serialization)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._click_handlers = CallbackDispatcher()
        self.on_msg(self._handle_button_msg)

    @validate('icon')
    def _validate_icon(self, proposal):
        """Strip 'fa-' if necessary'"""
        value = proposal['value']
        if 'fa-' in value:
            deprecation("icons names no longer need 'fa-', "
            "just use the class names themselves (for example, 'gear spin' instead of 'fa-gear fa-spin')",
            internal=['ipywidgets/widgets/', 'traitlets/traitlets.py', '/contextlib.py'])
            value = value.replace('fa-', '')
        return value

    def on_click(self, callback, remove=False):
        """Register a callback to execute when the button is clicked.

        The callback will be called with one argument, the clicked button
        widget instance.

        Parameters
        ----------
        remove: bool (optional)
            Set to true to remove the callback from the list of callbacks.
        """
        self._click_handlers.register_callback(callback, remove=remove)

    def click(self):
        """Programmatically trigger a click event.

        This will call the callbacks registered to the clicked button
        widget instance.
        """
        self._click_handlers(self)

    def _handle_button_msg(self, _, content, buffers):
        """Handle a msg from the front-end.

        Parameters
        ----------
        content: dict
            Content of the msg.
        """
        if content.get('event', '') == 'click':
            self.click()
