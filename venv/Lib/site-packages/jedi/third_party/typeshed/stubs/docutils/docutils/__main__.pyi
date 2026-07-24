from typing import ClassVar, Final

import docutils

__docformat__: Final = "reStructuredText"

class CliSettingsSpec(docutils.SettingsSpec):
    config_section: ClassVar[str]
    config_section_dependencies: ClassVar[tuple[str, ...]]

def main() -> None: ...
