from typing import Literal

from ._canonical_names import canonical_names as canonical_names, normalize_name as normalize_name

KEY_DOWN: Literal["down"]
KEY_UP: Literal["up"]

class KeyboardEvent:
    event_type: Literal["down", "up"] | None
    scan_code: int
    name: str | None
    time: float | None
    device: str | None
    modifiers: tuple[str, ...] | None
    is_keypad: bool | None

    def __init__(
        self,
        event_type: Literal["down", "up"] | None,
        scan_code: int,
        name: str | None = None,
        time: float | None = None,
        device: str | None = None,
        modifiers: tuple[str, ...] | None = None,
        is_keypad: bool | None = None,
    ) -> None: ...
    def to_json(self, ensure_ascii: bool = False) -> str: ...
    def __eq__(self, other: object) -> bool: ...
