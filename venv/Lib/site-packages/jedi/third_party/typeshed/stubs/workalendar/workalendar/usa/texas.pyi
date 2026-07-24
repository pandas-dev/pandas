from typing import ClassVar

from .core import UnitedStates

class TexasBase(UnitedStates):
    texas_include_confederate_heroes: ClassVar[bool]
    texas_include_independance_day: ClassVar[bool]
    texas_san_jacinto_day: ClassVar[bool]
    texas_emancipation_day: ClassVar[bool]
    texas_lyndon_johnson_day: ClassVar[bool]

class Texas(TexasBase): ...
