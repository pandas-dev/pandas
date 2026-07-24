rdates_property: property
exdates_property: property
rrules_property: property

def multi_language_text_property(main_prop: str, compatibility_prop: str, doc: str) -> property: ...
def single_int_property(prop: str, default: int, doc: str) -> property: ...
def single_utc_property(name: str, docs: str) -> property: ...
def single_string_property(name: str, docs: str, other_name: str | None = None) -> property: ...

color_property: property
sequence_property: property
categories_property: property
uid_property: property

__all__ = [
    "categories_property",
    "color_property",
    "exdates_property",
    "multi_language_text_property",
    "rdates_property",
    "rrules_property",
    "sequence_property",
    "single_int_property",
    "single_utc_property",
    "uid_property",
]
