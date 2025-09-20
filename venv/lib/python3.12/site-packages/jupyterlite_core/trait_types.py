from pathlib import Path

from traitlets import Container, TraitType


class CPath(TraitType):
    """A trait for casting to a Path. It might not actually exist yet"""

    def validate(self, obj, value) -> Path:
        if isinstance(value, Path):
            return value.resolve()

        try:
            return Path(str(value)).resolve()
        except Exception:
            self.error(obj, value)


class TypedTuple(Container):
    """A trait for a tuple of any length with type-checked elements.

    From: https://github.com/jupyter-widgets/ipywidgets/blob/7.6.3/ipywidgets/widgets/trait_types.py#L213
    """

    klass = tuple
    _cast_types = (list,)
