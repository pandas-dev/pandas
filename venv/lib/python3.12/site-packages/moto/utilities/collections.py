from typing import Any


def select_attributes(obj: Any, attributes: list[str]) -> Any:
    """Select a subset of attributes from the given dict (returns a copy)"""
    attributes = attributes if isinstance(attributes, (list, tuple)) else [attributes]  # type: ignore
    return {k: v for k, v in obj.items() if k in attributes}


def select_from_typed_dict(typed_dict: Any, obj: Any, filter: bool = False) -> Any:
    """
    Select a subset of attributes from a dictionary based on the keys of a given `TypedDict`.
    :param typed_dict: the `TypedDict` blueprint
    :param obj: the object to filter
    :param filter: if True, remove all keys with an empty (e.g., empty string or dictionary) or `None` value
    :return: the resulting dictionary (it returns a copy)
    """
    selection = select_attributes(
        obj, [*typed_dict.__required_keys__, *typed_dict.__optional_keys__]
    )
    if filter:
        selection = {k: v for k, v in selection.items() if v}
    return selection
