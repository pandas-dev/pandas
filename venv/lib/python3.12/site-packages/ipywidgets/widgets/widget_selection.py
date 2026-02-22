# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Selection classes.

Represents an enumeration using a widget.
"""

from collections.abc import Iterable, Mapping
from itertools import chain

from .widget_description import DescriptionWidget, DescriptionStyle
from .valuewidget import ValueWidget
from .widget_core import CoreWidget
from .widget_style import Style
from .trait_types import InstanceDict, TypedTuple
from .widget import register, widget_serialization
from .widget_int import SliderStyle
from .docutils import doc_subst
from traitlets import (Unicode, Bool, Int, Any, Dict, TraitError, CaselessStrEnum,
                       Tuple, Union, observe, validate)

_doc_snippets = {}
_doc_snippets['selection_params'] = """
    options: list
        The options for the dropdown. This can either be a list of values, e.g.
        ``['Galileo', 'Brahe', 'Hubble']`` or ``[0, 1, 2]``, a list of
        (label, value) pairs, e.g.
        ``[('Galileo', 0), ('Brahe', 1), ('Hubble', 2)]``, or a Mapping between
        labels and values, e.g., ``{'Galileo': 0, 'Brahe': 1, 'Hubble': 2}``.

    index: int
        The index of the current selection.

    value: any
        The value of the current selection. When programmatically setting the
        value, a reverse lookup is performed among the options to check that
        the value is valid. The reverse lookup uses the equality operator by
        default, but another predicate may be provided via the ``equals``
        keyword argument. For example, when dealing with numpy arrays, one may
        set ``equals=np.array_equal``.

    label: str
        The label corresponding to the selected value.

    disabled: bool
        Whether to disable user changes.

    description: str
        Label for this input group. This should be a string
        describing the widget.
"""

_doc_snippets['multiple_selection_params'] = """
    options: dict or list
        The options for the dropdown. This can either be a list of values, e.g.
        ``['Galileo', 'Brahe', 'Hubble']`` or ``[0, 1, 2]``, or a list of
        (label, value) pairs, e.g.
        ``[('Galileo', 0), ('Brahe', 1), ('Hubble', 2)]``, or a Mapping between
        labels and values, e.g., ``{'Galileo': 0, 'Brahe': 1, 'Hubble': 2}``.
        The labels are the strings that will be displayed in the UI,
        representing the actual Python choices, and should be unique.

    index: iterable of int
        The indices of the options that are selected.

    value: iterable
        The values that are selected. When programmatically setting the
        value, a reverse lookup is performed among the options to check that
        the value is valid. The reverse lookup uses the equality operator by
        default, but another predicate may be provided via the ``equals``
        keyword argument. For example, when dealing with numpy arrays, one may
        set ``equals=np.array_equal``.

    label: iterable of str
        The labels corresponding to the selected value.

    disabled: bool
        Whether to disable user changes.

    description: str
        Label for this input group. This should be a string
        describing the widget.
"""

_doc_snippets['slider_params'] = """
    orientation: str
        Either ``'horizontal'`` or ``'vertical'``. Defaults to ``horizontal``.

    readout: bool
        Display the current label next to the slider. Defaults to ``True``.

    continuous_update: bool
        If ``True``, update the value of the widget continuously as the user
        holds the slider. Otherwise, the model is only updated after the
        user has released the slider. Defaults to ``True``.
"""


def _exhaust_iterable(x):
    """Exhaust any non-mapping iterable into a tuple"""
    if isinstance(x, Iterable) and not isinstance(x, Mapping):
        return tuple(x)
    return x


def _make_options(x):
    """Standardize the options tuple format.

    The returned tuple should be in the format (('label', value), ('label', value), ...).

    The input can be
    * an iterable of (label, value) pairs
    * an iterable of values, and labels will be generated
    * a Mapping between labels and values
    """
    if isinstance(x, Mapping):
        x = x.items()

    # only iterate once through the options.
    xlist = tuple(x)

    # Check if x is an iterable of (label, value) pairs
    if all((isinstance(i, (list, tuple)) and len(i) == 2) for i in xlist):
        return tuple((str(k), v) for k, v in xlist)

    # Otherwise, assume x is an iterable of values
    return tuple((str(i), i) for i in xlist)

def findvalue(array, value, compare = lambda x, y: x == y):
    "A function that uses the compare function to return a value from the list."
    try:
        return next(x for x in array if compare(x, value))
    except StopIteration:
        raise ValueError('%r not in array'%value)

class _Selection(DescriptionWidget, ValueWidget, CoreWidget):
    """Base class for Selection widgets

    ``options`` can be specified as a list of values or a list of (label, value)
    tuples. The labels are the strings that will be displayed in the UI,
    representing the actual Python choices, and should be unique.
    If labels are not specified, they are generated from the values.

    When programmatically setting the value, a reverse lookup is performed
    among the options to check that the value is valid. The reverse lookup uses
    the equality operator by default, but another predicate may be provided via
    the ``equals`` keyword argument. For example, when dealing with numpy arrays,
    one may set equals=np.array_equal.
    """

    value = Any(None, help="Selected value", allow_none=True)
    label = Unicode(None, help="Selected label", allow_none=True)
    index = Int(None, help="Selected index", allow_none=True).tag(sync=True)

    options = Any((),
    help="""Iterable of values, (label, value) pairs, or Mapping between labels and values that the user can select.

    The labels are the strings that will be displayed in the UI, representing the
    actual Python choices, and should be unique.
    """)

    _options_full = None

    # This being read-only means that it cannot be changed by the user.
    _options_labels = TypedTuple(trait=Unicode(), read_only=True, help="The labels for the options.").tag(sync=True)

    disabled = Bool(help="Enable or disable user changes").tag(sync=True)

    def __init__(self, *args, **kwargs):
        self.equals = kwargs.pop('equals', lambda x, y: x == y)
        # We have to make the basic options bookkeeping consistent
        # so we don't have errors the first time validators run
        self._initializing_traits_ = True
        kwargs['options'] = _exhaust_iterable(kwargs.get('options', ()))
        self._options_full = _make_options(kwargs['options'])
        self._propagate_options(None)

        # Select the first item by default, if we can
        if 'index' not in kwargs and 'value' not in kwargs and 'label' not in kwargs:
            options = self._options_full
            nonempty = (len(options) > 0)
            kwargs['index'] = 0 if nonempty else None
            kwargs['label'], kwargs['value'] = options[0] if nonempty else (None, None)

        super().__init__(*args, **kwargs)
        self._initializing_traits_ = False

    @validate('options')
    def _validate_options(self, proposal):
        # if an iterator is provided, exhaust it
        proposal.value = _exhaust_iterable(proposal.value)
        # throws an error if there is a problem converting to full form
        self._options_full = _make_options(proposal.value)
        return proposal.value

    @observe('options')
    def _propagate_options(self, change):
        "Set the values and labels, and select the first option if we aren't initializing"
        options = self._options_full
        self.set_trait('_options_labels', tuple(i[0] for i in options))
        self._options_values = tuple(i[1] for i in options)

        if self.index is None:
            # Do nothing, we don't want to force a selection if
            # the options list changed
            return

        if self._initializing_traits_ is not True:
            if len(options) > 0:
                if self.index == 0:
                    # Explicitly trigger the observers to pick up the new value and
                    # label. Just setting the value would not trigger the observers
                    # since traitlets thinks the value hasn't changed.
                    self._notify_trait('index', 0, 0)
                else:
                    self.index = 0
            else:
                self.index = None

    @validate('index')
    def _validate_index(self, proposal):
        if proposal.value is None or 0 <= proposal.value < len(self._options_labels):
            return proposal.value
        else:
            raise TraitError('Invalid selection: index out of bounds')

    @observe('index')
    def _propagate_index(self, change):
        "Propagate changes in index to the value and label properties"
        label = self._options_labels[change.new] if change.new is not None else None
        value = self._options_values[change.new] if change.new is not None else None
        if self.label is not label:
            self.label = label
        if self.value is not value:
            self.value = value

    @validate('value')
    def _validate_value(self, proposal):
        value = proposal.value
        try:
            return findvalue(self._options_values, value, self.equals) if value is not None else None
        except ValueError:
            raise TraitError('Invalid selection: value not found')

    @observe('value')
    def _propagate_value(self, change):
        if change.new is None:
            index = None
        elif self.index is not None and self.equals(self._options_values[self.index], change.new):
            index = self.index
        else:
            index = self._options_values.index(change.new)
        if self.index != index:
            self.index = index

    @validate('label')
    def _validate_label(self, proposal):
        if (proposal.value is not None) and (proposal.value not in self._options_labels):
            raise TraitError('Invalid selection: label not found')
        return proposal.value

    @observe('label')
    def _propagate_label(self, change):
        if change.new is None:
            index = None
        elif self.index is not None and self._options_labels[self.index] == change.new:
            index = self.index
        else:
            index = self._options_labels.index(change.new)
        if self.index != index:
            self.index = index

    def _repr_keys(self):
        keys = super()._repr_keys()
        # Include options manually, as it isn't marked as synced:
        for key in sorted(chain(keys, ('options',))):
            if key == 'index' and self.index == 0:
                # Index 0 is default when there are options
                continue
            yield key


class _MultipleSelection(DescriptionWidget, ValueWidget, CoreWidget):
    """Base class for multiple Selection widgets

    ``options`` can be specified as a list of values, list of (label, value)
    tuples, or a dict of {label: value}. The labels are the strings that will be
    displayed in the UI, representing the actual Python choices, and should be
    unique. If labels are not specified, they are generated from the values.

    When programmatically setting the value, a reverse lookup is performed
    among the options to check that the value is valid. The reverse lookup uses
    the equality operator by default, but another predicate may be provided via
    the ``equals`` keyword argument. For example, when dealing with numpy arrays,
    one may set equals=np.array_equal.
    """

    value = TypedTuple(trait=Any(), help="Selected values")
    label = TypedTuple(trait=Unicode(), help="Selected labels")
    index = TypedTuple(trait=Int(), help="Selected indices").tag(sync=True)

    options = Any((),
    help="""Iterable of values, (label, value) pairs, or Mapping between labels and values that the user can select.

    The labels are the strings that will be displayed in the UI, representing the
    actual Python choices, and should be unique.
    """)
    _options_full = None

    # This being read-only means that it cannot be changed from the frontend!
    _options_labels = TypedTuple(trait=Unicode(), read_only=True, help="The labels for the options.").tag(sync=True)

    disabled = Bool(help="Enable or disable user changes").tag(sync=True)

    def __init__(self, *args, **kwargs):
        self.equals = kwargs.pop('equals', lambda x, y: x == y)

        # We have to make the basic options bookkeeping consistent
        # so we don't have errors the first time validators run
        self._initializing_traits_ = True
        kwargs['options'] = _exhaust_iterable(kwargs.get('options', ()))
        self._options_full = _make_options(kwargs['options'])
        self._propagate_options(None)

        super().__init__(*args, **kwargs)
        self._initializing_traits_ = False

    @validate('options')
    def _validate_options(self, proposal):
        proposal.value = _exhaust_iterable(proposal.value)
        # throws an error if there is a problem converting to full form
        self._options_full = _make_options(proposal.value)
        return proposal.value

    @observe('options')
    def _propagate_options(self, change):
        "Unselect any option"
        options = self._options_full
        self.set_trait('_options_labels', tuple(i[0] for i in options))
        self._options_values = tuple(i[1] for i in options)
        if self._initializing_traits_ is not True:
            self.index = ()

    @validate('index')
    def _validate_index(self, proposal):
        "Check the range of each proposed index."
        if all(0 <= i < len(self._options_labels) for i in proposal.value):
            return proposal.value
        else:
            raise TraitError('Invalid selection: index out of bounds')

    @observe('index')
    def _propagate_index(self, change):
        "Propagate changes in index to the value and label properties"
        label = tuple(self._options_labels[i] for i in change.new)
        value = tuple(self._options_values[i] for i in change.new)
        # we check equality so we can avoid validation if possible
        if self.label != label:
            self.label = label
        if self.value != value:
            self.value = value

    @validate('value')
    def _validate_value(self, proposal):
        "Replace all values with the actual objects in the options list"
        try:
            return tuple(findvalue(self._options_values, i, self.equals) for i in proposal.value)
        except ValueError:
            raise TraitError('Invalid selection: value not found')

    @observe('value')
    def _propagate_value(self, change):
        index = tuple(self._options_values.index(i) for i in change.new)
        if self.index != index:
            self.index = index

    @validate('label')
    def _validate_label(self, proposal):
        if any(i not in self._options_labels for i in proposal.value):
            raise TraitError('Invalid selection: label not found')
        return proposal.value

    @observe('label')
    def _propagate_label(self, change):
        index = tuple(self._options_labels.index(i) for i in change.new)
        if self.index != index:
            self.index = index

    def _repr_keys(self):
        keys = super()._repr_keys()
        # Include options manually, as it isn't marked as synced:
        yield from sorted(chain(keys, ('options',)))


@register
class ToggleButtonsStyle(DescriptionStyle, CoreWidget):
    """Button style widget.

    Parameters
    ----------
    button_width: str
        The width of each button. This should be a valid CSS
        width, e.g. '10px' or '5em'.

    font_weight: str
        The text font weight of each button, This should be a valid CSS font
        weight unit, for example 'bold' or '600'
    """
    _model_name = Unicode('ToggleButtonsStyleModel').tag(sync=True)
    button_width = Unicode(help="The width of each button.").tag(sync=True)
    font_weight = Unicode(help="Text font weight of each button.").tag(sync=True)


@register
@doc_subst(_doc_snippets)
class ToggleButtons(_Selection):
    """Group of toggle buttons that represent an enumeration.

    Only one toggle button can be toggled at any point in time.

    Parameters
    ----------
    {selection_params}

    tooltips: list
        Tooltip for each button. If specified, must be the
        same length as `options`.

    icons: list
        Icons to show on the buttons. This must be the name
        of a font-awesome icon. See `http://fontawesome.io/icons/`
        for a list of icons.

    button_style: str
        One of 'primary', 'success', 'info', 'warning' or
        'danger'. Applies a predefined style to every button.

    style: ToggleButtonsStyle
        Style parameters for the buttons.
    """
    _view_name = Unicode('ToggleButtonsView').tag(sync=True)
    _model_name = Unicode('ToggleButtonsModel').tag(sync=True)

    tooltips = TypedTuple(Unicode(), help="Tooltips for each button.").tag(sync=True)
    icons = TypedTuple(Unicode(), help="Icons names for each button (FontAwesome names without the fa- prefix).").tag(sync=True)
    style = InstanceDict(ToggleButtonsStyle).tag(sync=True, **widget_serialization)

    button_style = CaselessStrEnum(
        values=['primary', 'success', 'info', 'warning', 'danger', ''],
        default_value='', allow_none=True, help="""Use a predefined styling for the buttons.""").tag(sync=True)


@register
@doc_subst(_doc_snippets)
class Dropdown(_Selection):
    """Allows you to select a single item from a dropdown.

    Parameters
    ----------
    {selection_params}
    """
    _view_name = Unicode('DropdownView').tag(sync=True)
    _model_name = Unicode('DropdownModel').tag(sync=True)


@register
@doc_subst(_doc_snippets)
class RadioButtons(_Selection):
    """Group of radio buttons that represent an enumeration.

    Only one radio button can be toggled at any point in time.

    Parameters
    ----------
    {selection_params}
    """
    _view_name = Unicode('RadioButtonsView').tag(sync=True)
    _model_name = Unicode('RadioButtonsModel').tag(sync=True)

    orientation = CaselessStrEnum(
        values=['horizontal', 'vertical'], default_value='vertical',
        help="Vertical or horizontal.").tag(sync=True)


@register
@doc_subst(_doc_snippets)
class Select(_Selection):
    """
    Listbox that only allows one item to be selected at any given time.

    Parameters
    ----------
    {selection_params}

    rows: int
        The number of rows to display in the widget.
    """
    _view_name = Unicode('SelectView').tag(sync=True)
    _model_name = Unicode('SelectModel').tag(sync=True)
    rows = Int(5, help="The number of rows to display.").tag(sync=True)

@register
@doc_subst(_doc_snippets)
class SelectMultiple(_MultipleSelection):
    """
    Listbox that allows many items to be selected at any given time.

    The ``value``, ``label`` and ``index`` attributes are all iterables.

    Parameters
    ----------
    {multiple_selection_params}

    rows: int
        The number of rows to display in the widget.
    """
    _view_name = Unicode('SelectMultipleView').tag(sync=True)
    _model_name = Unicode('SelectMultipleModel').tag(sync=True)
    rows = Int(5, help="The number of rows to display.").tag(sync=True)


class _SelectionNonempty(_Selection):
    """Selection that is guaranteed to have a value selected."""
    # don't allow None to be an option.
    value = Any(help="Selected value")
    label = Unicode(help="Selected label")
    index = Int(help="Selected index").tag(sync=True)

    def __init__(self, *args, **kwargs):
        if len(kwargs.get('options', ())) == 0:
            raise TraitError('options must be nonempty')
        super().__init__(*args, **kwargs)

    @validate('options')
    def _validate_options(self, proposal):
        proposal.value = _exhaust_iterable(proposal.value)
        self._options_full = _make_options(proposal.value)
        if len(self._options_full) == 0:
            raise TraitError("Option list must be nonempty")
        return proposal.value

    @validate('index')
    def _validate_index(self, proposal):
        if 0 <= proposal.value < len(self._options_labels):
            return proposal.value
        else:
            raise TraitError('Invalid selection: index out of bounds')

class _MultipleSelectionNonempty(_MultipleSelection):
    """Selection that is guaranteed to have an option available."""

    def __init__(self, *args, **kwargs):
        if len(kwargs.get('options', ())) == 0:
            raise TraitError('options must be nonempty')
        super().__init__(*args, **kwargs)

    @validate('options')
    def _validate_options(self, proposal):
        proposal.value = _exhaust_iterable(proposal.value)
        # throws an error if there is a problem converting to full form
        self._options_full = _make_options(proposal.value)
        if len(self._options_full) == 0:
            raise TraitError("Option list must be nonempty")
        return proposal.value

@register
@doc_subst(_doc_snippets)
class SelectionSlider(_SelectionNonempty):
    """
    Slider to select a single item from a list or dictionary.

    Parameters
    ----------
    {selection_params}

    {slider_params}
    """
    _view_name = Unicode('SelectionSliderView').tag(sync=True)
    _model_name = Unicode('SelectionSliderModel').tag(sync=True)

    orientation = CaselessStrEnum(
        values=['horizontal', 'vertical'], default_value='horizontal',
        help="Vertical or horizontal.").tag(sync=True)
    readout = Bool(True,
        help="Display the current selected label next to the slider").tag(sync=True)
    continuous_update = Bool(True,
        help="Update the value of the widget as the user is holding the slider.").tag(sync=True)
    behavior = CaselessStrEnum(values=['drag-tap', 'drag-snap', 'tap', 'drag', 'snap'],
        default_value='drag-tap', help="Slider dragging behavior.").tag(sync=True)

    style = InstanceDict(SliderStyle).tag(sync=True, **widget_serialization)


@register
@doc_subst(_doc_snippets)
class SelectionRangeSlider(_MultipleSelectionNonempty):
    """
    Slider to select multiple contiguous items from a list.

    The index, value, and label attributes contain the start and end of
    the selection range, not all items in the range.

    Parameters
    ----------
    {multiple_selection_params}

    {slider_params}
    """
    _view_name = Unicode('SelectionRangeSliderView').tag(sync=True)
    _model_name = Unicode('SelectionRangeSliderModel').tag(sync=True)

    value = Tuple(help="Min and max selected values")
    label = Tuple(help="Min and max selected labels")
    index = Tuple((0,0), help="Min and max selected indices").tag(sync=True)

    @observe('options')
    def _propagate_options(self, change):
        "Select the first range"
        options = self._options_full
        self.set_trait('_options_labels', tuple(i[0] for i in options))
        self._options_values = tuple(i[1] for i in options)
        if self._initializing_traits_ is not True:
            self.index = (0, 0)

    @validate('index')
    def _validate_index(self, proposal):
        "Make sure we have two indices and check the range of each proposed index."
        if len(proposal.value) != 2:
            raise TraitError('Invalid selection: index must have two values, but is %r'%(proposal.value,))
        if all(0 <= i < len(self._options_labels) for i in proposal.value):
            return proposal.value
        else:
            raise TraitError('Invalid selection: index out of bounds: %s'%(proposal.value,))

    orientation = CaselessStrEnum(
        values=['horizontal', 'vertical'], default_value='horizontal',
        help="Vertical or horizontal.").tag(sync=True)
    readout = Bool(True,
        help="Display the current selected label next to the slider").tag(sync=True)
    continuous_update = Bool(True,
        help="Update the value of the widget as the user is holding the slider.").tag(sync=True)

    style = InstanceDict(SliderStyle).tag(sync=True, **widget_serialization)
    behavior = CaselessStrEnum(values=['drag-tap', 'drag-snap', 'tap', 'drag', 'snap'],
        default_value='drag-tap', help="Slider dragging behavior.").tag(sync=True)
