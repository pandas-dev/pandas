from _typeshed import Incomplete
from typing import Any

class Constants:
    @staticmethod
    def get_all_constant_values_from_class(klass: object) -> list[Any]: ...  # Any taken from klass.__dict__
    def get_all_enum_values(enum_class) -> list[Incomplete]: ...
