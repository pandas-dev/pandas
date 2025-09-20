# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""String class.

Represents a unicode string using a widget.
"""

from .widget_description import DescriptionStyle, DescriptionWidget
from .valuewidget import ValueWidget
from .widget import CallbackDispatcher, register, widget_serialization
from .widget_core import CoreWidget
from .trait_types import Color, InstanceDict, TypedTuple
from .utils import deprecation
from traitlets import Unicode, Bool, Int


class _StringStyle(DescriptionStyle, CoreWidget):
    """Text input style widget."""
    _model_name = Unicode('StringStyleModel').tag(sync=True)
    background = Unicode(None, allow_none=True, help="Background specifications.").tag(sync=True)
    font_size = Unicode(None, allow_none=True, help="Text font size.").tag(sync=True)
    text_color = Color(None, allow_none=True, help="Text color").tag(sync=True)


@register
class LabelStyle(_StringStyle):
    """Label style widget."""
    _model_name = Unicode('LabelStyleModel').tag(sync=True)
    font_family = Unicode(None, allow_none=True, help="Label text font family.").tag(sync=True)
    font_style = Unicode(None, allow_none=True, help="Label text font style.").tag(sync=True)
    font_variant = Unicode(None, allow_none=True, help="Label text font variant.").tag(sync=True)
    font_weight = Unicode(None, allow_none=True, help="Label text font weight.").tag(sync=True)
    text_decoration = Unicode(None, allow_none=True, help="Label text decoration.").tag(sync=True)


@register
class TextStyle(_StringStyle):
    """Text input style widget."""
    _model_name = Unicode('TextStyleModel').tag(sync=True)

@register
class HTMLStyle(_StringStyle):
    """HTML style widget."""
    _model_name = Unicode('HTMLStyleModel').tag(sync=True)

@register
class HTMLMathStyle(_StringStyle):
    """HTML with math style widget."""
    _model_name = Unicode('HTMLMathStyleModel').tag(sync=True)


class _String(DescriptionWidget, ValueWidget, CoreWidget):
    """Base class used to create widgets that represent a string."""

    value = Unicode(help="String value").tag(sync=True)

    # We set a zero-width space as a default placeholder to make sure the baseline matches
    # the text, not the bottom margin. See the last paragraph of
    # https://www.w3.org/TR/CSS2/visudet.html#leading
    placeholder = Unicode('\u200b', help="Placeholder text to display when nothing has been typed").tag(sync=True)
    style = InstanceDict(_StringStyle).tag(sync=True, **widget_serialization)

    def __init__(self, value=None, **kwargs):
        if value is not None:
            kwargs['value'] = value
        super().__init__(**kwargs)

    _model_name = Unicode('StringModel').tag(sync=True)

@register
class HTML(_String):
    """Renders the string `value` as HTML."""
    _view_name = Unicode('HTMLView').tag(sync=True)
    _model_name = Unicode('HTMLModel').tag(sync=True)
    style = InstanceDict(HTMLStyle).tag(sync=True, **widget_serialization)

@register
class HTMLMath(_String):
    """Renders the string `value` as HTML, and render mathematics."""
    _view_name = Unicode('HTMLMathView').tag(sync=True)
    _model_name = Unicode('HTMLMathModel').tag(sync=True)
    style = InstanceDict(HTMLMathStyle).tag(sync=True, **widget_serialization)


@register
class Label(_String):
    """Label widget.

    It also renders math inside the string `value` as Latex (requires $ $ or
    $$ $$ and similar latex tags).
    """
    _view_name = Unicode('LabelView').tag(sync=True)
    _model_name = Unicode('LabelModel').tag(sync=True)
    style = InstanceDict(LabelStyle).tag(sync=True, **widget_serialization)


@register
class Textarea(_String):
    """Multiline text area widget."""
    _view_name = Unicode('TextareaView').tag(sync=True)
    _model_name = Unicode('TextareaModel').tag(sync=True)
    rows = Int(None, allow_none=True, help="The number of rows to display.").tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    continuous_update = Bool(True, help="Update the value as the user types. If False, update on submission, e.g., pressing Enter or navigating away.").tag(sync=True)
    style = InstanceDict(TextStyle).tag(sync=True, **widget_serialization)

@register
class Text(_String):
    """Single line textbox widget."""
    _view_name = Unicode('TextView').tag(sync=True)
    _model_name = Unicode('TextModel').tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)
    continuous_update = Bool(True, help="Update the value as the user types. If False, update on submission, e.g., pressing Enter or navigating away.").tag(sync=True)
    style = InstanceDict(TextStyle).tag(sync=True, **widget_serialization)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._submission_callbacks = CallbackDispatcher()
        self.on_msg(self._handle_string_msg)

    def _handle_string_msg(self, _, content, buffers):
        """Handle a msg from the front-end.

        Parameters
        ----------
        content: dict
            Content of the msg.
        """
        if content.get('event', '') == 'submit':
            self._submission_callbacks(self)

    def on_submit(self, callback, remove=False):
        """(Un)Register a callback to handle text submission.

        Triggered when the user clicks enter.

        Parameters
        ----------
        callback: callable
            Will be called with exactly one argument: the Widget instance
        remove: bool (optional)
            Whether to unregister the callback
        """
        deprecation("on_submit is deprecated. Instead, set the .continuous_update attribute to False and observe the value changing with: mywidget.observe(callback, 'value').")
        self._submission_callbacks.register_callback(callback, remove=remove)


@register
class Password(Text):
    """Single line textbox widget."""
    _view_name = Unicode('PasswordView').tag(sync=True)
    _model_name = Unicode('PasswordModel').tag(sync=True)
    disabled = Bool(False, help="Enable or disable user changes").tag(sync=True)

    def _repr_keys(self):
        # Don't include password value in repr!
        super_keys = super()._repr_keys()
        for key in super_keys:
            if key != 'value':
                yield key


@register
class Combobox(Text):
    """Single line textbox widget with a dropdown and autocompletion.
    """
    _model_name = Unicode('ComboboxModel').tag(sync=True)
    _view_name = Unicode('ComboboxView').tag(sync=True)

    options = TypedTuple(
        trait=Unicode(),
        help="Dropdown options for the combobox"
    ).tag(sync=True)

    ensure_option = Bool(
        False,
        help='If set, ensure value is in options. Implies continuous_update=False.'
    ).tag(sync=True)
