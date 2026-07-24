import abc
from abc import abstractmethod
from logging import Logger

class InvalidValidator(Exception): ...

class BaseValidator(metaclass=abc.ABCMeta):
    log: Logger
    def __init__(self) -> None: ...
    @abstractmethod
    def validate(self, input_string: str) -> bool: ...
