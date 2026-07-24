from consolemenu.console_menu import ConsoleMenu
from consolemenu.items import ExternalItem as ExternalItem

class CommandItem(ExternalItem):
    command: str
    arguments: list[str]
    exit_status: int | None
    def __init__(
        self,
        text: str,
        command: str,
        arguments: list[str] | None = None,
        menu: ConsoleMenu | None = None,
        should_exit: bool = False,
        menu_char: str | None = None,
    ) -> None: ...
    def action(self) -> None: ...
    def get_return(self) -> int: ...
