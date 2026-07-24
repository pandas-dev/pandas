from _typeshed import Incomplete
from logging import Logger
from typing import ClassVar

class LanguageValidator:
    VALID_KEYS: ClassVar[list[str]]
    logger: Logger | None
    @classmethod
    def get_logger(cls) -> Logger: ...
    @classmethod
    def validate_info(cls, language_id, info: dict[str, Incomplete]) -> bool: ...
