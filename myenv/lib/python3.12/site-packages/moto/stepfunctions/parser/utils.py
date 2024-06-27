import re
from typing import Dict, List, Set, Type

TMP_THREADS = []


def is_list_or_tuple(obj) -> bool:
    return isinstance(obj, (list, tuple))


def select_attributes(obj: Dict, attributes: List[str]) -> Dict:
    """Select a subset of attributes from the given dict (returns a copy)"""
    attributes = attributes if is_list_or_tuple(attributes) else [attributes]
    return {k: v for k, v in obj.items() if k in attributes}


class SubtypesInstanceManager:
    """Simple instance manager base class that scans the subclasses of a base type for concrete named
    implementations, and lazily creates and returns (singleton) instances on demand."""

    _instances: Dict[str, "SubtypesInstanceManager"]

    @classmethod
    def get(cls, subtype_name: str):
        instances = cls.instances()
        base_type = cls.get_base_type()
        instance = instances.get(subtype_name)
        if instance is None:
            # lazily load subtype instance (required if new plugins are dynamically loaded at runtime)
            for clazz in get_all_subclasses(base_type):
                impl_name = clazz.impl_name()
                if impl_name not in instances:
                    instances[impl_name] = clazz()
            instance = instances.get(subtype_name)
        return instance

    @classmethod
    def instances(cls) -> Dict[str, "SubtypesInstanceManager"]:
        base_type = cls.get_base_type()
        if not hasattr(base_type, "_instances"):
            base_type._instances = {}
        return base_type._instances

    @staticmethod
    def impl_name() -> str:
        """Name of this concrete subtype - to be implemented by subclasses."""
        raise NotImplementedError

    @classmethod
    def get_base_type(cls) -> Type:
        """Get the base class for which instances are being managed - can be customized by subtypes."""
        return cls


def to_str(obj, encoding: str = "utf-8", errors="strict") -> str:
    """If ``obj`` is an instance of ``binary_type``, return
    ``obj.decode(encoding, errors)``, otherwise return ``obj``"""
    return obj.decode(encoding, errors) if isinstance(obj, bytes) else obj


_re_camel_to_snake_case = re.compile("((?<=[a-z0-9])[A-Z]|(?!^)[A-Z](?=[a-z]))")


def camel_to_snake_case(string: str) -> str:
    return _re_camel_to_snake_case.sub(r"_\1", string).replace("__", "_").lower()


def snake_to_camel_case(string: str, capitalize_first: bool = True) -> str:
    components = string.split("_")
    start_idx = 0 if capitalize_first else 1
    components = [x.title() for x in components[start_idx:]]
    return "".join(components)


def get_all_subclasses(clazz: Type) -> Set[Type]:
    """Recursively get all subclasses of the given class."""
    result = set()
    subs = clazz.__subclasses__()
    for sub in subs:
        result.add(sub)
        result.update(get_all_subclasses(sub))
    return result
