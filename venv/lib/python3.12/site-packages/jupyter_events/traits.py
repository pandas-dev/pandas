"""Trait types for events."""
from __future__ import annotations

import logging
import typing as t

from traitlets import TraitError, TraitType

baseclass = TraitType
if t.TYPE_CHECKING:
    baseclass = TraitType[t.Any, t.Any]  # type:ignore[misc]


class Handlers(baseclass):  # type:ignore[type-arg]
    """A trait that takes a list of logging handlers and converts
    it to a callable that returns that list (thus, making this
    trait pickleable).
    """

    info_text = "a list of logging handlers"

    def validate_elements(self, obj: t.Any, value: t.Any) -> None:
        """Validate the elements of an object."""
        if len(value) > 0:
            # Check that all elements are logging handlers.
            for el in value:
                if isinstance(el, logging.Handler) is False:
                    self.element_error(obj)

    def element_error(self, obj: t.Any) -> None:
        """Raise an error for bad elements."""
        msg = f"Elements in the '{self.name}' trait of an {obj.__class__.__name__} instance must be Python `logging` handler instances."
        raise TraitError(msg)

    def validate(self, obj: t.Any, value: t.Any) -> t.Any:
        """Validate an object."""
        # If given a callable, call it and set the
        # value of this trait to the returned list.
        # Verify that the callable returns a list
        # of logging handler instances.
        if callable(value):
            out = value()
            self.validate_elements(obj, out)
            return out
        # If a list, check it's elements to verify
        # that each element is a logging handler instance.
        if isinstance(value, list):
            self.validate_elements(obj, value)
            return value
        self.error(obj, value)
        return None  # type:ignore[unreachable]
