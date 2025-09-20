# Copyright(c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""TagsInput class.

Represents a list of tags.
"""

from traitlets import (
    CaselessStrEnum, CInt, CFloat, Bool, Unicode, List, TraitError, validate
)

from .widget_description import DescriptionWidget
from .valuewidget import ValueWidget
from .widget_core import CoreWidget
from .widget import register
from .trait_types import Color, NumberFormat


class TagsInputBase(DescriptionWidget, ValueWidget, CoreWidget):
    _model_name = Unicode('TagsInputBaseModel').tag(sync=True)
    value = List().tag(sync=True)
    placeholder = Unicode('\u200b').tag(sync=True)
    allowed_tags = List().tag(sync=True)
    allow_duplicates = Bool(True).tag(sync=True)

    @validate('value')
    def _validate_value(self, proposal):
        if ('' in proposal['value']):
            raise TraitError('The value of a TagsInput widget cannot contain blank strings')

        if len(self.allowed_tags) == 0:
            return proposal['value']

        for tag_value in proposal['value']:
            if tag_value not in self.allowed_tags:
                raise TraitError('Tag value {} is not allowed, allowed tags are {}'.format(tag_value, self.allowed_tags))

        return proposal['value']


@register
class TagsInput(TagsInputBase):
    """
    List of string tags
    """
    _model_name = Unicode('TagsInputModel').tag(sync=True)
    _view_name = Unicode('TagsInputView').tag(sync=True)

    value = List(Unicode(), help='List of string tags').tag(sync=True)
    tag_style = CaselessStrEnum(
        values=['primary', 'success', 'info', 'warning', 'danger', ''], default_value='',
        help="""Use a predefined styling for the tags.""").tag(sync=True)


@register
class ColorsInput(TagsInputBase):
    """
    List of color tags
    """
    _model_name = Unicode('ColorsInputModel').tag(sync=True)
    _view_name = Unicode('ColorsInputView').tag(sync=True)

    value = List(Color(), help='List of string tags').tag(sync=True)


class NumbersInputBase(TagsInput):
    _model_name = Unicode('NumbersInputBaseModel').tag(sync=True)
    min = CFloat(default_value=None, allow_none=True).tag(sync=True)
    max = CFloat(default_value=None, allow_none=True).tag(sync=True)

    @validate('value')
    def _validate_numbers(self, proposal):
        for tag_value in proposal['value']:
            if self.min is not None and tag_value < self.min:
                raise TraitError('Tag value {} should be >= {}'.format(tag_value, self.min))
            if self.max is not None and tag_value > self.max:
                raise TraitError('Tag value {} should be <= {}'.format(tag_value, self.max))

        return proposal['value']


@register
class FloatsInput(NumbersInputBase):
    """
    List of float tags
    """
    _model_name = Unicode('FloatsInputModel').tag(sync=True)
    _view_name = Unicode('FloatsInputView').tag(sync=True)

    value = List(CFloat(), help='List of float tags').tag(sync=True)
    format = NumberFormat('.1f').tag(sync=True)


@register
class IntsInput(NumbersInputBase):
    """
    List of int tags
    """
    _model_name = Unicode('IntsInputModel').tag(sync=True)
    _view_name = Unicode('IntsInputView').tag(sync=True)

    value = List(CInt(), help='List of int tags').tag(sync=True)
    format = NumberFormat('d').tag(sync=True)
    min = CInt(default_value=None, allow_none=True).tag(sync=True)
    max = CInt(default_value=None, allow_none=True).tag(sync=True)
