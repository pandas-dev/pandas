from typing import Any, Final

# String constant used to indicate all model fields should be included
ALL_FIELDS: Final[str] = "__all__"
# Collection of values considered empty by Django filters - tuple type allows various empty containers
EMPTY_VALUES: Final[Any] = ...
