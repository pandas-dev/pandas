from _typeshed import Unused
from typing import Any

DEFAULTS: dict[str, Any]  # Configuration values can be strings, booleans, callables, etc.
DEPRECATED_SETTINGS: list[str]

def is_callable(value: Any) -> bool: ...  # Accepts any value to test if it's callable

class Settings:
    # Setting values can be of any type, so getter and setter methods return/accept Any
    def __getattr__(self, name: str) -> Any: ...  # Returns setting values of various types
    def get_setting(self, setting: str) -> Any: ...  # Setting values vary by configuration option
    def change_setting(
        self, setting: str, value: Any, enter: bool, **kwargs: Unused
    ) -> None: ...  # Accepts any setting value type

settings: Settings
