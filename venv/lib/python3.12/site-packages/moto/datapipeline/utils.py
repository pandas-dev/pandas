import collections.abc as collections_abc
from typing import Any

from moto.moto_api._internal import mock_random


def get_random_pipeline_id() -> str:
    return f"df-{mock_random.get_random_hex(length=19)}"


def remove_capitalization_of_dict_keys(obj: Any) -> Any:
    if isinstance(obj, collections_abc.Mapping):
        result = obj.__class__()
        for key, value in obj.items():
            normalized_key = key[:1].lower() + key[1:]
            result[normalized_key] = remove_capitalization_of_dict_keys(value)  # type: ignore[index]
        return result
    elif isinstance(obj, collections_abc.Iterable) and not isinstance(obj, str):
        result = obj.__class__()  # type: ignore[assignment]
        for item in obj:
            result += (remove_capitalization_of_dict_keys(item),)  # type: ignore[operator]
        return result
    else:
        return obj
