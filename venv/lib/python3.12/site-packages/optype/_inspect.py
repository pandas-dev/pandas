import typing

__all__ = ["_get_protocol_attrs"]


# same as from `typing_extensions._EXCLUDED_ATTRS`
_EXCLUDED_ATTRS: typing.Final = frozenset(getattr(typing, "EXCLUDED_ATTRIBUTES")) | {  # noqa: B009
    "__match_args__",
    "__protocol_attrs__",
    "__non_callable_proto_members__",
    "__final__",
}


# same as `typing_extensions._get_protocol_attrs`
def _get_protocol_attrs(cls: type) -> set[str]:
    attrs: set[str] = set()
    for base in cls.__mro__[:-1]:  # without object
        if base.__name__ in {"Protocol", "Generic"}:
            continue
        annotations = getattr(base, "__annotations__", {})
        for attr in (*base.__dict__, *annotations):
            if not attr.startswith("_abc_") and attr not in _EXCLUDED_ATTRS:
                attrs.add(attr)
    return attrs
