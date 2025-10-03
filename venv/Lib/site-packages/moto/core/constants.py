from typing import Final


class _Missing:
    """A singleton class to represent a sentinel for a missing value."""

    def __bool__(self) -> bool:
        return False

    def __copy__(self):  # type: ignore[no-untyped-def]
        return self

    def __deepcopy__(self, _):  # type: ignore[no-untyped-def]
        return self

    def __repr__(self) -> str:
        return "<moto.MISSING>"


MISSING: Final = _Missing()
