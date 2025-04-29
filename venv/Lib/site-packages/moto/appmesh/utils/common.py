from typing import Any, Dict


def clean_dict(obj: Dict[str, Any]) -> Dict[str, Any]:  # type: ignore[misc]
    return {
        key: value for key, value in obj.items() if value is not None and value != []
    }
