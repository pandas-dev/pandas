# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Link and DirectionalLink classes.

Propagate changes between widgets on the javascript side.
"""

from .widget import Widget, register, widget_serialization
from .widget_core import CoreWidget

from traitlets import Unicode, Tuple, Instance, TraitError


class WidgetTraitTuple(Tuple):
    """Traitlet for validating a single (Widget, 'trait_name') pair"""

    info_text = "A (Widget, 'trait_name') pair"

    def __init__(self, **kwargs):
        super().__init__(Instance(Widget), Unicode(), **kwargs)
        if "default_value" not in kwargs and not kwargs.get("allow_none", False):
            # This is to keep consistent behavior for spec generation between traitlets 4 and 5
            # Having a default empty container is explicitly not allowed in traitlets 5 when
            # there are traits specified (as the default value will be invalid), but we do it
            # anyway as there is no empty "default" that makes sense.
            self.default_args = ()

    def validate_elements(self, obj, value):
        value = super().validate_elements(obj, value)
        widget, trait_name = value
        trait = widget.traits().get(trait_name)
        trait_repr = "{}.{}".format(widget.__class__.__name__, trait_name)
        # Can't raise TraitError because the parent will swallow the message
        # and throw it away in a new, less informative TraitError
        if trait is None:
            raise TypeError("No such trait: %s" % trait_repr)
        elif not trait.metadata.get('sync'):
            raise TypeError("%s cannot be synced" % trait_repr)
        return value


@register
class Link(CoreWidget):
    """Link Widget

    source: a (Widget, 'trait_name') tuple for the source trait
    target: a (Widget, 'trait_name') tuple that should be updated
    """

    _model_name = Unicode('LinkModel').tag(sync=True)
    target = WidgetTraitTuple(help="The target (widget, 'trait_name') pair").tag(sync=True, **widget_serialization)
    source = WidgetTraitTuple(help="The source (widget, 'trait_name') pair").tag(sync=True, **widget_serialization)

    def __init__(self, source, target, **kwargs):
        kwargs['source'] = source
        kwargs['target'] = target
        super().__init__(**kwargs)

    # for compatibility with traitlet links
    def unlink(self):
        self.close()


def jslink(attr1, attr2):
    """Link two widget attributes on the frontend so they remain in sync.

    The link is created in the front-end and does not rely on a roundtrip
    to the backend.

    Parameters
    ----------
    source : a (Widget, 'trait_name') tuple for the first trait
    target : a (Widget, 'trait_name') tuple for the second trait

    Examples
    --------

    >>> c = link((widget1, 'value'), (widget2, 'value'))
    """
    return Link(attr1, attr2)


@register
class DirectionalLink(Link):
    """A directional link

    source: a (Widget, 'trait_name') tuple for the source trait
    target: a (Widget, 'trait_name') tuple that should be updated
    when the source trait changes.
    """
    _model_name = Unicode('DirectionalLinkModel').tag(sync=True)


def jsdlink(source, target):
    """Link a source widget attribute with a target widget attribute.

    The link is created in the front-end and does not rely on a roundtrip
    to the backend.

    Parameters
    ----------
    source : a (Widget, 'trait_name') tuple for the source trait
    target : a (Widget, 'trait_name') tuple for the target trait

    Examples
    --------

    >>> c = dlink((src_widget, 'value'), (tgt_widget, 'value'))
    """
    return DirectionalLink(source, target)
