from typing import ClassVar, Final, Protocol, type_check_only
from typing_extensions import Self

from docutils.utils import Reporter

__docformat__: Final = "reStructuredText"

@type_check_only
class _LanguageModule(Protocol):
    labels: dict[str, str]
    author_separators: list[str]
    bibliographic_fields: list[str]

class LanguageImporter:
    packages: ClassVar[tuple[str, ...]]
    warn_msg: ClassVar[str]
    fallback: ClassVar[str]
    cache: dict[str, _LanguageModule]
    def __init__(self) -> None: ...
    def import_from_packages(self, name: str, reporter: Reporter | None = None) -> _LanguageModule: ...
    def check_content(self, module: _LanguageModule) -> None: ...
    def __call__(self, language_code: str, reporter: Reporter | None = None) -> _LanguageModule: ...
    def __class_getitem__(cls, name) -> type[Self]: ...

get_language: LanguageImporter
