from pathlib import Path

from traitlets import Container, TraitType


class CPath(TraitType):
    """A trait for casting to a Path. It might not actually exist yet"""

    def __init__(self, *args, resolve_relative=True, **kwargs):
        self.resolve_relative = resolve_relative
        super().__init__(*args, **kwargs)

    def validate(self, obj, value) -> Path:
        try:
            path = value if isinstance(value, Path) else Path(str(value))
            if path.is_absolute() or self.resolve_relative:
                return path.resolve()
            return path
        except Exception:
            self.error(obj, value)


class TypedTuple(Container):
    """A trait for a tuple of any length with type-checked elements.

    From: https://github.com/jupyter-widgets/ipywidgets/blob/7.6.3/ipywidgets/widgets/trait_types.py#L213
    """

    klass = tuple
    _cast_types = (list,)
