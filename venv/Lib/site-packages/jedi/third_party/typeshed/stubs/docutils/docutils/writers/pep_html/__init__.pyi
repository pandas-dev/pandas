from typing import ClassVar, Final

from docutils.writers import html4css1

__docformat__: Final = "reStructuredText"

class Writer(html4css1.Writer):
    default_stylesheet: ClassVar[str]
    default_stylesheet_path: ClassVar[str]
    default_template_path: ClassVar[str]
    settings_default_overrides: ClassVar[dict[str, str]]
    relative_path_settings: ClassVar[tuple[str, ...]]
    config_section_dependencies: ClassVar[tuple[str, ...]]
    translator_class: type[HTMLTranslator]
    pepnum: str
    title: str
    def interpolation_dict(self) -> dict[str, str | int]: ...  # type: ignore[override]

class HTMLTranslator(html4css1.HTMLTranslator):
    def depart_field_list(self, node) -> None: ...
