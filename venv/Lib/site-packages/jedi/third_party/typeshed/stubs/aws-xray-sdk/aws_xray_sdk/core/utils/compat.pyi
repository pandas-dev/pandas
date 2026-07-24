annotation_value_types: tuple[type, ...]

def is_classmethod(func: object) -> bool: ...  # argument func is passing to getattr() function
def is_instance_method(parent_class: type, func_name: str, func: object) -> bool: ...
