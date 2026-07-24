from typing import ClassVar, Final, Protocol, type_check_only

from docutils.languages import LanguageImporter
from docutils.utils import Reporter

__docformat__: Final = "reStructuredText"

@type_check_only
class _RstLanguageModule(Protocol):
    directives: dict[str, str]
    roles: dict[str, str]

class RstLanguageImporter(LanguageImporter):
    cache: dict[str, _RstLanguageModule]  # type: ignore[assignment]
    fallback: ClassVar[None]  # type: ignore[assignment]
    def import_from_packages(self, name: str, reporter: Reporter | None = None) -> _RstLanguageModule: ...  # type: ignore[override]
    def check_content(self, module: _RstLanguageModule) -> None: ...  # type: ignore[override]
    def __call__(self, language_code: str, reporter: Reporter | None = None) -> _RstLanguageModule: ...  # type: ignore[override]

get_language: RstLanguageImporter
