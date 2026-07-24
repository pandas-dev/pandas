from typing import Any


def clean_dict(obj: dict[str, Any]) -> dict[str, Any]:  # type: ignore[misc]
    return {
        key: value for key, value in obj.items() if value is not None and value != []
    }
