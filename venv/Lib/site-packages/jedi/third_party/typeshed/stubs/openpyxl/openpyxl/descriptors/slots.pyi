from _typeshed import Incomplete, Self

class AutoSlotProperties(type):
    def __new__(mcl: type[Self], classname: str, bases: tuple[type, ...], dictionary: dict[str, Incomplete]) -> Self: ...
